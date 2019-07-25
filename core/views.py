import django_filters
import numpy as np
import pandas as pd

from django.conf import settings
from django.db import connection
from django.views.decorators.cache import cache_page
from rest_framework import filters, generics, permissions, status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from core.models import Species, Assembly, Analysis, BackspliceJunction, Sample
from core.serializers import SpeciesDetailSerializer, SpeciesListSerializer, SampleListSerializer, SampleDetailsSerializer

caching_time = 60*5 if settings.DEBUG else 10*24*60*60


class SpeciesList(generics.ListAPIView):
    """
    View to list and search of species
    """

    filter_backends = (filters.SearchFilter,)
    ordering = ('scietific_name', )
    pagination_class = None
    queryset = Species.objects.all().order_by('scientific_name')
    serializer_class = SpeciesListSerializer
    search_fields = ['scientific_name', 'common_name']


class SpeciesDetail(generics.RetrieveAPIView):
    """
    View to get details of a species
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesDetailSerializer


class SampleFilter(django_filters.FilterSet):
    """
    Filter parameters for Samples
    """

    accession = django_filters.CharFilter(
        'accession', lookup_expr='icontains')
    source = django_filters.CharFilter('source', lookup_expr='icontains')
    description = django_filters.CharFilter(
        'description', lookup_expr='icontains')

    class Meta:
        model = Sample
        fields = ['accession', 'source', 'description']


class SampleList(generics.ListAPIView):
    """
    View to list and search of samples
    """

    serializer_class = SampleListSerializer
    pagination_class = None
    filter_backends = (django_filters.rest_framework.DjangoFilterBackend,
                       filters.SearchFilter)
    filter_class = SampleFilter

    def get_queryset(self):
        """
        This view should return a list of all the samples of
        the species as determined by the species portion of the URL.
        """
        species_id = self.kwargs['species_id']
        assembly_id = self.kwargs['assembly_id']
        sample_id_list = Analysis.objects.filter(
            assembly_id=assembly_id).values_list('sample_id', flat=True)
        return Sample.objects.filter(species_id=species_id, sample_id__in=sample_id_list)


class SampleDetail(generics.RetrieveAPIView):
    """
    View to get details of a sample
    """

    queryset = Sample.objects.all()
    serializer_class = SampleDetailsSerializer


@api_view(['GET'])
@cache_page(caching_time)
def species_view_stats(request, species_id, assembly_id):
    """
    View to get stats for the species view
    """

    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        return Response(data={'error': 'No species with the given id.'},
                        status=status.HTTP_404_NOT_FOUND)

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        return Response(data={'error': 'No assembly with the given id under the given species.'}, status=status.HTTP_404_NOT_FOUND)

    # List of all samples related to the given assembly
    sample_query = 'select source, sample_id from core_sample where sample_id in (select sample_id_id from core_analysis where assembly_id_id={});'.format(
        assembly.assembly_id)
    sample_df = pd.read_sql_query(sample_query, connection)

    # List of locus expression related to the given assembly
    locus_expression_query = 'select tpm, analysis_id_id from core_locusexpression where analysis_id_id in (select analysis_id from core_analysis where assembly_id_id={});'.format(
        assembly.assembly_id)
    locusexpression_df = pd.read_sql_query(locus_expression_query, connection)

    # List of all locus for a given assembly
    locus_query = 'select locus_id, nexons from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all analysis related to the assembly
    analysis_query = 'select sample_id_id, analysis_id from core_analysis where assembly_id_id={};'.format(
        assembly.assembly_id)
    analysis_df = pd.read_sql_query(analysis_query, connection)

    # List of all canonical junctions related to given assembly
    cj_query = 'select locus_id_id from core_canonicaljunction where locus_id_id in (select locus_id from core_locus where assembly_id_id={})'.format(
        assembly.assembly_id)
    cj_df = pd.read_sql_query(cj_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select locus_id_id from core_backsplicejunction where locus_id_id in (select locus_id from core_locus where assembly_id_id={})'.format(
        assembly.assembly_id)
    bj_df = pd.read_sql_query(bj_query, connection)

    # List of distinct tissues scanned
    distinct_tissues = sample_df['source'].unique().tolist()

    # Total number of samples for the given assembly
    count_total_samples = len(sample_df.index)

    # Count of the number of distinct tissues scanned for the assembly
    count_distinct_tissues = len(distinct_tissues)

    # Count of backsplice junctions related to the given assembly
    count_backsplice_junctions = len(bj_df)

    # Graph for circRNA for each locus
    circRNA_per_locus_df = bj_df[['locus_id_id']].groupby(
        ['locus_id_id']).size().reset_index(name='count')
    circRNA_per_locus = {
        "locus_id": circRNA_per_locus_df["locus_id_id"].apply(str).to_list(),
        "count": circRNA_per_locus_df["count"].to_list()
    }

    # Graph for circRNA vs linear transcripts for each locus
    cj_per_locus_df = cj_df[['locus_id_id']].groupby(
        ['locus_id_id']).size().reset_index(name='count')
    locus_cj_df = pd.merge(locus_df[['locus_id', 'nexons']], cj_per_locus_df[['locus_id_id', 'count']],
                           left_on='locus_id', right_on='locus_id_id', how='outer').fillna(0, downcast='infer')
    locus_cj_bj_df = pd.merge(locus_cj_df, circRNA_per_locus_df[['locus_id_id', 'count']],
                              left_on='locus_id', right_on='locus_id_id', how='outer', suffixes=('_cj', '_bj')).fillna(0, downcast='infer')
    circrna_vs_lt_per_locus = {
        "locus_id": locus_cj_bj_df["locus_id"].to_list(),
        "count_cj": locus_cj_bj_df["count_cj"].to_list(),
        "count_bj": locus_cj_bj_df["count_bj"].to_list(),
        "nexons": locus_cj_bj_df["nexons"].to_list(),
        "text": "Locus " + locus_cj_bj_df['locus_id'].map(str) + "<br>Exons: " + locus_cj_bj_df['nexons'].map(str)
    }

    # Graph for tissue TPM box plot
    sample_analysis_query = 'select analysis_id, sample_id_id from core_analysis where assembly_id_id={};'.format(
        assembly.assembly_id)
    sample_analysis_df = pd.read_sql_query(sample_analysis_query, connection)
    tpm_sample_analysis_df = pd.merge(
        locusexpression_df[['tpm', 'analysis_id_id']], sample_analysis_df, left_on='analysis_id_id', right_on='analysis_id')
    tpm_sample_merged_df = pd.merge(
        tpm_sample_analysis_df, sample_df, left_on='sample_id_id', right_on='sample_id')
    tpm_sample_merged_df.fillna(0)
    tpm_sample_merged_df['tpm'] = np.log10(tpm_sample_merged_df['tpm']+1)
    tpm_sample_grouped_df = tpm_sample_merged_df.groupby(['source'])[
        'tpm'].apply(list)
    tpm_grouped_list = tpm_sample_grouped_df.to_dict()
    source_list = tpm_sample_grouped_df.keys().to_list()
    tpm_tissue_boxplot = {
        'tpm_grouped': tpm_grouped_list,
        'tissue_list': source_list
    }

    data = {'species': species.scientific_name,
            'assembly': assembly.assembly_name,
            'count_total_samples': count_total_samples,
            'count_distinct_tissues': count_distinct_tissues,
            'count_backsplice_junctions': count_backsplice_junctions,
            'circRNA_per_locus': circRNA_per_locus,
            'circrna_vs_lt_per_locus': circrna_vs_lt_per_locus,
            'tpm_tissue_boxplot': tpm_tissue_boxplot
            }

    return Response(data=data, status=status.HTTP_200_OK)


@api_view(['GET'])
@cache_page(caching_time)
def sample_view_stats(request, species_id, assembly_id, sample_id):
    """
    View to get stats for the sample view
    """

    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        return Response(data={'error': 'No species with the given id.'},
                        status=status.HTTP_404_NOT_FOUND)

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        return Response(data={'error': 'No assembly with the given id under the given species.'}, status=status.HTTP_404_NOT_FOUND)

    try:
        sample = species.samples.get(sample_id=sample_id)
    except:
        return Response(data={'error': 'No sample with the given id under the given species.'}, status=status.HTTP_404_NOT_FOUND)

    library_size = str(int(sample.library_size or 0))
    mapped_reads = str(int(sample.mapped_reads or 0))
    total_spliced_reads = str(int(sample.total_spliced_reads or 0))
    canonical_spliced_reads = str(
        int(sample.total_spliced_reads or 0)-int(sample.backspliced_reads or 0))
    backspliced_reads = str(int(sample.backspliced_reads or 0))
    unmapped_reads = str(int(sample.library_size or 0) -
                         int(sample.mapped_reads or 0))

    sankey_labels = ['Library size ({})'.format(library_size),
                     'Mapped reads ({})'.format(mapped_reads),
                     'Total spliced reads ({})'.format(total_spliced_reads),
                     'Canonical spliced reads ({})'.format(
                         canonical_spliced_reads),
                     'Backspliced reads ({})'.format(backspliced_reads),
                     'Unmapped reads ({})'.format(unmapped_reads)
                     ]

    sankey_values = [int(sample.mapped_reads or 0), int(sample.total_spliced_reads or 0),
                     int(sample.total_spliced_reads or 0) -
                     int(sample.backspliced_reads or 0),
                     int(sample.backspliced_reads or 0),
                     int(sample.library_size or 0) -
                     int(sample.mapped_reads or 0)
                     ]
    sankey = {'sankey_labels': sankey_labels,
              'sankey_values': sankey_values}

    data = {'species': species.scientific_name,
            'assembly': assembly.assembly_name,
            'sankey': sankey
            }
    return Response(data=data, status=status.HTTP_200_OK)


@api_view(['GET'])
@cache_page(caching_time)
def location_view_stats(request, species_id, assembly_id):
    """
    View to get stats for the sample view
    """

    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        return Response(data={'error': 'No species with the given id.'},
                        status=status.HTTP_404_NOT_FOUND)

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        return Response(data={'error': 'No assembly with the given id under the given species.'}, status=status.HTTP_404_NOT_FOUND)

    return Response(data={'species': [species.scientific_name, species.description], 'genome': assembly.assembly_genoverse_genome}, status=status.HTTP_200_OK)
