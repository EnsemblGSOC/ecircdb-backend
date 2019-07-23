from django.db import connection
import pandas as pd
from rest_framework import filters, generics, permissions, status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from core.models import Species, Assembly, Analysis, BackspliceJunction
from core.serializers import SpeciesDetailSerializer, SpeciesListSerializer
import django_filters
from rest_framework import filters, generics, permissions

from core.models import Species, Sample
from core.serializers import SpeciesDetailSerializer, SpeciesListSerializer, SampleListSerializer, SampleDetailsSerializer


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

    # Get all of the analyses related to given assembly
    analyses = Analysis.objects.filter(assembly_id=assembly.assembly_id)

    # Get ids of all the samples for the given assembly
    sample_id_list = analyses.values_list('sample_id', flat=True)
    # All the samples particular to given assembly
    samples = Sample.objects.filter(
        species_id=species_id, sample_id__in=sample_id_list)
    count_total_samples = samples.count()

    # List of the distinct type of tissues scanned
    distinct_tissues = samples.values_list('source', flat=True).distinct()
    count_distinct_tissues = distinct_tissues.count()

    # Get ids of all analysis related to the given assembly
    analysis_id_list = analyses.values_list('analysis_id', flat=True)
    backsplice_junctions = BackspliceJunction.objects.filter(
        analysis_id__in=analysis_id_list).distinct()
    count_backsplice_junctions = backsplice_junctions.count()

    # Graph for circRNA for each locus
    circRNA_per_locus_query = 'select locus_id_id, count(*) as count from core_backsplicejunction where locus_id_id in (select locus_id from core_locus where assembly_id_id={}) group by locus_id_id;'.format(
        assembly.assembly_id)
    circRNA_per_locus_dataframe = pd.read_sql_query(
        circRNA_per_locus_query, connection)
    circRNA_per_locus = {
        "locus_id": circRNA_per_locus_dataframe["locus_id_id"].apply(str).to_list(),
        "count": circRNA_per_locus_dataframe["count"].to_list()
    }

    # Graph for circRNA vs linear transcripts for each locus
    locus_nexons_query = 'select locus_id, nexons from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    cj_query = 'select locus_id_id, count(*) as count from core_canonicaljunction where locus_id_id in (select locus_id from core_locus where assembly_id_id={}) group by locus_id_id;'.format(
        assembly.assembly_id)
    locus_nexons_df = pd.read_sql_query(locus_nexons_query, connection)
    cj_df = pd.read_sql_query(cj_query, connection)
    locus_cj_df = pd.merge(locus_nexons_df, cj_df,
                           left_on='locus_id', right_on='locus_id_id', how='outer').fillna(0, downcast='infer')
    locus_cj_bj_df = pd.merge(locus_cj_df, circRNA_per_locus_dataframe,
                              left_on='locus_id', right_on='locus_id_id', how='outer', suffixes=('_cj', '_bj')).fillna(0, downcast='infer')
    circrna_vs_lt_per_locus = {
        "locus_id": locus_cj_bj_df["locus_id"].to_list(),
        "count_cj": locus_cj_bj_df["count_cj"].to_list(),
        "count_bj": locus_cj_bj_df["count_bj"].to_list(),
        "nexons": locus_cj_bj_df["nexons"].to_list(),
        "text": "Locus " + locus_cj_bj_df['locus_id'].map(str) + "<br>Exons: " + locus_cj_bj_df['nexons'].map(str)
    }

    # Graph for tissue TPM box plot
    tpm_query = 'select tpm, analysis_id_id from core_locusexpression where analysis_id_id in (select analysis_id from core_analysis where assembly_id_id={});'.format(
        assembly.assembly_id)
    sample_query = 'select source, sample_id from core_sample where sample_id in (select sample_id_id from core_analysis where assembly_id_id={});'.format(
        assembly.assembly_id)
    sample_analysis_query = 'select analysis_id, sample_id_id from core_analysis where assembly_id_id={};'.format(
        assembly.assembly_id)
    tpm_df = pd.read_sql_query(tpm_query, connection)
    sample_df = pd.read_sql_query(sample_query, connection)
    sample_analysis_df = pd.read_sql_query(sample_analysis_query, connection)
    tpm_sample_analysis_df = pd.merge(
        tpm_df, sample_analysis_df, left_on='analysis_id_id', right_on='analysis_id')
    tpm_sample_merged_df = pd.merge(
        tpm_sample_analysis_df, sample_df, left_on='sample_id_id', right_on='sample_id')
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
