import datetime
import django_filters
import json
import numpy as np
import pandas as pd
import requests

from django.conf import settings
from django.db import connection
from django.http import HttpResponse, Http404
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

    # List of all analysis related to the assembly
    analysis_query = 'select sample_id_id, analysis_id from core_analysis where assembly_id_id={};'.format(
        assembly.assembly_id)
    analysis_df = pd.read_sql_query(analysis_query, connection)

    # List of all samples related to the given assembly
    sample_query = 'select * from core_sample where sample_id in ({});'.format(
        str(analysis_df['sample_id_id'].to_list())[1:-1])
    sample_df = pd.read_sql_query(sample_query, connection)

    # List of locus expression related to the given assembly
    locusexpression_query = 'select tpm, analysis_id_id from core_locusexpression where analysis_id_id in ({});'.format(
        str(analysis_df['analysis_id'].to_list())[1:-1])
    locusexpression_df = pd.read_sql_query(locusexpression_query, connection)

    # List of all locus for a given assembly
    locus_query = 'select * from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select * from core_backsplicejunction where locus_id_id in ({})'.format(
        str(locus_df["locus_id"].to_list())[1:-1])
    bj_df = pd.read_sql_query(bj_query, connection)

    # List of distinct tissues scanned
    distinct_tissues = sample_df['source'].unique().tolist()

    # Total number of samples for the given assembly
    count_total_samples = len(sample_df.index)

    # Count of the number of distinct tissues scanned for the assembly
    count_distinct_tissues = len(distinct_tissues)

    # Total number of the circRNAs
    total_circrnas = bj_df['coord_id'].nunique()

    # Sum of the library size
    sum_library_size = sample_df['library_size'].sum()

    # Total number of circRNA count / sum(library_size)
    circrna_per_library_size = str(
        round((total_circrnas/sum_library_size)*1000000, 2)) + ' X 10e-6'

    # circRNA per sample
    circrna_per_sample = total_circrnas/count_total_samples

    # Count of backsplice junctions related to the given assembly
    count_backsplice_junctions = len(bj_df)

    # Number of circRNA producing genes
    count_circrna_producing_genes = locus_df[['locus_id', 'gene_name']].merge(
        bj_df[['locus_id_id']], left_on='locus_id', right_on='locus_id_id')['gene_name'].nunique()

    # Histogram for number of genes per number of circRNAs
    locus_bj_merged_circRNA_count_df = locus_df[['locus_id']].merge(bj_df[['locus_id_id', 'coord_id']], left_on='locus_id', right_on='locus_id_id')[
        ['locus_id', 'coord_id']].groupby('locus_id').agg({'coord_id': pd.Series.nunique}).reset_index()
    gene_count_per_circrna_count_df = locus_bj_merged_circRNA_count_df.groupby(
        'coord_id').agg({'locus_id': pd.Series.nunique}).reset_index()
    gene_count_gte_30 = gene_count_per_circrna_count_df.loc[gene_count_per_circrna_count_df.coord_id > 30].locus_id.sum(
    )
    gene_count_per_circrna_count_df = gene_count_per_circrna_count_df.drop(
        gene_count_per_circrna_count_df[gene_count_per_circrna_count_df.coord_id > 30].index)
    gene_count_per_circrna_count_df = gene_count_per_circrna_count_df.append(
        {'locus_id': gene_count_gte_30, 'coord_id': '30+'}, ignore_index=True)
    gene_count_per_circrna_count = {
        'circrna_count': [str(x)+' circRNAs' for x in gene_count_per_circrna_count_df['coord_id'].tolist()],
        'labels': [str(x) for x in gene_count_per_circrna_count_df['coord_id'].tolist()],
        'genes_frequency': gene_count_per_circrna_count_df['locus_id'].tolist()
    }

    # Graph for circRNA for each locus
    circRNA_per_locus_df = bj_df[['locus_id_id', 'coord_id']].groupby(
        ['locus_id_id'])['coord_id'].nunique().reset_index(name='count')
    circRNA_per_locus = {
        "locus_id": circRNA_per_locus_df["locus_id_id"].apply(str).to_list(),
        "count": circRNA_per_locus_df["count"].to_list()
    }

    # Pie chart for circRNA classification
    values = bj_df[['classification', 'coord_id']].groupby(
        ['classification'])['coord_id'].nunique().tolist()
    labels = bj_df[['classification', 'coord_id']].groupby(
        ['classification'])['coord_id'].nunique().keys().tolist()
    circrna_classification = {
        'values': values,
        'labels': labels
    }

    # Pie chart/Treemap for circRNA n_methods
    n_methods_count_df = bj_df[['n_methods', 'coord_id']].groupby(
        ['n_methods'])['coord_id'].nunique().reset_index(name='count')
    count_n_methods_gt_3 = n_methods_count_df.loc[n_methods_count_df.n_methods > 3]['count'].sum(
    )
    n_methods_count_df = n_methods_count_df[n_methods_count_df.n_methods <= 3]
    n_methods_count_df = n_methods_count_df.append(
        {'n_methods': '>3', 'count': count_n_methods_gt_3}, ignore_index=True)
    values = n_methods_count_df['count'].tolist()
    labels = n_methods_count_df['n_methods'].tolist()
    circrna_n_methods = {
        'values': values,
        'labels': labels
    }

    # Histogram for circRNA n_exons count frequency
    circrna_n_exons_df = bj_df[['coord_id',
                                'predicted_exons']].drop_duplicates()
    circrna_n_exons_df.loc[:, 'count_predicted_exons'] = circrna_n_exons_df['predicted_exons'].apply(
        lambda x: len(x.split(',')))
    circrna_n_exons_series = circrna_n_exons_df[['count_predicted_exons']].groupby(
        ['count_predicted_exons']).size()
    circrna_n_exons_df = pd.DataFrame(
        {'nexons': circrna_n_exons_series.index, 'counts': circrna_n_exons_series.values})
    circrna_n_exons_gte_30 = circrna_n_exons_df.loc[circrna_n_exons_df.nexons > 30].counts.sum(
    )
    circrna_n_exons_df = circrna_n_exons_df.drop(
        circrna_n_exons_df[circrna_n_exons_df.nexons > 30].index)
    circrna_n_exons_df = circrna_n_exons_df.append(
        {'counts': circrna_n_exons_gte_30, 'nexons': '30+'}, ignore_index=True)
    circrna_n_exons = {
        'x': ['NExons-'+str(x) for x in circrna_n_exons_df['nexons'].tolist()],
        'labels': [str(x) for x in circrna_n_exons_df['nexons'].tolist()],
        'y': circrna_n_exons_df['counts'].tolist()
    }

    # Bar plot for circRNAs per chromosomes
    circrna_chromosomes = bj_df[['coord_id', 'seq_region_name']].groupby(['seq_region_name'])[
        'coord_id'].nunique().reset_index(name='circrnas').sort_values('circrnas', ascending=False).head(30)
    chromosomes_circrna_count = {
        'chromosomes': circrna_chromosomes['seq_region_name'].tolist(),
        'circrnas': circrna_chromosomes['circrnas'].tolist()
    }

    # DF, grouped by tissue for distinct circRNAs aggregating corresponding fields
    sample_analysis_query = 'select analysis_id, sample_id_id from core_analysis where assembly_id_id={};'.format(
        assembly.assembly_id)
    sample_analysis_df = pd.read_sql_query(sample_analysis_query, connection)
    bj_sample_analysis_df = pd.merge(
        bj_df[['tpm', 'jpm', 'abundance_ratio', 'genomic_size', 'spliced_size', 'analysis_id_id', 'coord_id']], sample_analysis_df, left_on='analysis_id_id', right_on='analysis_id')
    bj_sample_merged_df = pd.merge(
        bj_sample_analysis_df, sample_df[['sample_id', 'source']], left_on='sample_id_id', right_on='sample_id')
    bj_sample_merged_df = bj_sample_merged_df.fillna(0)
    bj_sample_merged_df = bj_sample_merged_df[['tpm', 'jpm', 'abundance_ratio', 'genomic_size', 'spliced_size', 'coord_id', 'source']].groupby(
        ['source', 'coord_id']).mean().reset_index('coord_id').reset_index('source')[['source', 'tpm', 'jpm', 'abundance_ratio', 'genomic_size', 'spliced_size']]
    source_list = bj_sample_merged_df[['source']].groupby(
        ['source']).apply(list).keys().to_list()

    # Graph for tissue circRNA size violin plot
    circrna_size_sample_merged_df = bj_sample_merged_df[[
        'spliced_size', 'genomic_size', 'source']]
    circrna_size_sample_merged_df.loc[:, 'spliced_size'] = np.log2(
        circrna_size_sample_merged_df[['spliced_size']] + 1)
    circrna_size_sample_merged_df.loc[:, 'genomic_size'] = np.log2(
        circrna_size_sample_merged_df[['genomic_size']] + 1)
    circrna_size_sample_merged_df = circrna_size_sample_merged_df.round(3)
    tissues_list = circrna_size_sample_merged_df['source'].tolist()
    spliced_sizes = circrna_size_sample_merged_df['spliced_size'].tolist()
    genomic_sizes = circrna_size_sample_merged_df['genomic_size'].tolist()
    size_tissue_boxplot = {
        'spliced_sizes': spliced_sizes,
        'genomic_sizes': genomic_sizes,
        'tissue_list': tissues_list
    }

    # Graph for tissue TPM box plot
    tpm_sample_merged_df = bj_sample_merged_df[['source', 'tpm']]
    tpm_sample_merged_df.loc[:, 'tpm'] = np.log2(
        tpm_sample_merged_df[['tpm']]+1)
    tpm_sample_merged_df = tpm_sample_merged_df.round(3)
    tpm_sample_grouped_df = tpm_sample_merged_df.groupby(['source'])[
        'tpm'].apply(list)
    tpm_grouped_list = tpm_sample_grouped_df.to_dict()
    tpm_tissue_boxplot = {
        'tpm_grouped': tpm_grouped_list,
        'tissue_list': source_list
    }

    # Graph for tissue JPM box plot
    jpm_sample_merged_df = bj_sample_merged_df[['source', 'jpm']]
    jpm_grouped = {}
    for tissue in source_list:
        data = jpm_sample_merged_df[jpm_sample_merged_df.source ==
                                    tissue]['jpm'].values
        data = np.log2(data)
        jpm_grouped[tissue] = [(x - np.mean(data)) /
                               np.std(data) for x in data]
    jpm_tissue_boxplot = {
        'jpm_grouped': jpm_grouped,
        'tissue_list': source_list
    }

    # Graph for tissue AR box plot
    ar_sample_merged_df = bj_sample_merged_df[['source', 'abundance_ratio']]
    ar_sample_merged_df = ar_sample_merged_df[ar_sample_merged_df['abundance_ratio'] != 0]
    ar_sample_merged_df.loc[:, 'abundance_ratio'] = np.log2(
        ar_sample_merged_df['abundance_ratio'])
    ar_sample_merged_df = ar_sample_merged_df.round(3)
    ar_sample_grouped_df = ar_sample_merged_df.groupby(['source'])[
        'abundance_ratio'].apply(list)
    ar_grouped_list = ar_sample_grouped_df.to_dict()
    ar_tissue_boxplot = {
        'ar_grouped': ar_grouped_list,
        'tissue_list': source_list
    }

    # List of distinct classifications among backsplice junctions
    distinct_classifications = bj_df['classification'].unique().tolist()
    distinct_classifications.sort(key=lambda x: (len(x), x))

    data = {'species': species.scientific_name,
            'assembly': assembly.assembly_name,
            'count_total_samples': count_total_samples,
            'count_distinct_tissues': count_distinct_tissues,
            'count_backsplice_junctions': count_backsplice_junctions,
            'circrna_per_library_size': circrna_per_library_size,
            'circrna_per_sample': circrna_per_sample,
            'count_circrna_producing_genes': count_circrna_producing_genes,
            'circRNA_per_locus': circRNA_per_locus,
            'gene_count_per_circrna_count': gene_count_per_circrna_count,
            'circrna_classification': circrna_classification,
            'circrna_n_methods': circrna_n_methods,
            'circrna_n_exons': circrna_n_exons,
            'tpm_tissue_boxplot': tpm_tissue_boxplot,
            'jpm_tissue_boxplot': jpm_tissue_boxplot,
            'ar_tissue_boxplot': ar_tissue_boxplot,
            'size_tissue_boxplot': size_tissue_boxplot,
            'chromosomes_circrna_count': chromosomes_circrna_count,
            'distinct_classifications': distinct_classifications,
            'distinct_tissues': distinct_tissues
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

    # List of all analysis related to the given sample
    analysis_query = 'select * from core_analysis where sample_id_id={}'.format(
        sample.sample_id)
    analysis_df = pd.read_sql_query(analysis_query, connection)

    # List of all the genes of the given assembly
    gene_query = 'select * from core_ensemblgene where assembly_id_id={}'.format(
        assembly.assembly_id
    )
    gene_df = pd.read_sql_query(gene_query, connection)

    # List of all locus expressions related to the given sample
    locusexpression_query = 'select * from core_locusexpression where analysis_id_id in ({})'.format(
        str(analysis_df['analysis_id'].to_list())[1:-1])
    locusexpression_df = pd.read_sql_query(locusexpression_query, connection)

    # List of all locus related to the given sample
    locus_query = 'select * from core_locus where locus_id in ({})'.format(
        str(locusexpression_df['locus_id_id'].to_list())[1:-1])
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all canonical junctions related to given assembly
    cj_query = 'select * from core_canonicaljunction where analysis_id_id in ({})'.format(
        str(analysis_df["analysis_id"].to_list())[1:-1])
    cj_df = pd.read_sql_query(cj_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select * from core_backsplicejunction where analysis_id_id in ({})'.format(
        str(analysis_df["analysis_id"].to_list())[1:-1])
    bj_df = pd.read_sql_query(bj_query, connection)

    lib_size = int(sample.library_size or 0)
    library_size = str(
        round((int(sample.library_size or 0)/lib_size*100), 2))+'%'
    mapped_reads = str(
        round((int(sample.mapped_reads or 0)/lib_size*100), 2))+'%'
    total_spliced_reads = str(
        round(((int(sample.total_spliced_reads or 0) + int(sample.chimeric_reads or 0))/lib_size*100), 2))+'%'
    canonical_spliced_reads = str(
        round(((int(sample.total_spliced_reads or 0)-int(sample.backspliced_reads or 0))/lib_size*100), 2))+'%'
    backspliced_reads = str(
        round((int(sample.backspliced_reads or bj_df['junction_id'].size)/lib_size*100), 2))+'%'
    unmapped_reads = str(round(((int(sample.library_size or 0) -
                                 int(sample.mapped_reads or 0))/lib_size*100), 2))+'%'
    chimeric_reads = str(
        round((int(sample.chimeric_reads or 0)/lib_size*100), 2))+'%'
    non_spliced_reads = str(round(((int(sample.mapped_reads or 0) -
                                    int(sample.total_spliced_reads or 0))/lib_size*100), 2))+'%'

    sankey_labels = ['Library size ({})'.format(library_size),
                     'Total spliced reads ({})'.format(total_spliced_reads),
                     'Canonical spliced reads ({})'.format(
        canonical_spliced_reads),
        'Backspliced reads ({})'.format(backspliced_reads),
        'Unmapped reads ({})'.format(unmapped_reads),
        'Chimeric reads ({})'.format(chimeric_reads),
        'Non spliced reads ({})'.format(non_spliced_reads)
    ]

    sankey_values = [int(sample.total_spliced_reads or 0),
                     int(sample.total_spliced_reads or 0) -
                     int(sample.backspliced_reads or 0),
                     int(sample.library_size or 0) -
                     int(sample.mapped_reads or 0),
                     int(sample.chimeric_reads or 0),
                     int(
        sample.backspliced_reads or bj_df['junction_id'].size),
        int(sample.mapped_reads or 0) -
        int(sample.total_spliced_reads or 0)
    ]
    sankey = {'sankey_labels': sankey_labels,
              'sankey_values': sankey_values}

    # Gene level AR dist plot
    gene_level_ar_sum = np.log2(locus_df[['gene_name', 'locus_id']].merge(bj_df[['locus_id_id', 'abundance_ratio', 'coord_id']], left_on='locus_id', right_on='locus_id_id')[
        ['gene_name', 'abundance_ratio', 'coord_id']].drop_duplicates()[['gene_name', 'abundance_ratio']].groupby('gene_name')['abundance_ratio'].sum().tolist())

    # JPM boxplot for circRNAs and canonical junction
    bj_jpm_list = np.log2(bj_df[['coord_id', 'jpm']].groupby(
        'coord_id').mean().reset_index()['jpm'].tolist())
    cj_jpm_list = np.log2(cj_df[['coord_id', 'jpm']].groupby(
        'coord_id').mean().reset_index()['jpm'].tolist())
    jpm_boxplot = {
        'bj_jpm_list': bj_jpm_list,
        'cj_jpm_list': cj_jpm_list
    }

    # circRNA vs canonical gene-level JPM scatterplot
    cj_unique_jpm = cj_df[['locus_id_id', 'coord_id', 'jpm']].drop_duplicates()
    bj_unique_jpm = bj_df[['locus_id_id', 'coord_id', 'jpm']].drop_duplicates()
    locus_cj_jpm_df = locus_df[['locus_id', 'gene_name']].merge(cj_unique_jpm[['locus_id_id', 'jpm']], left_on='locus_id', right_on='locus_id_id', suffixes=(
        'locustable_', 'cjtable_'))[['gene_name', 'jpm']].groupby('gene_name').sum().reset_index()
    locus_bj_jpm_df = locus_df[['locus_id', 'gene_name']].merge(bj_unique_jpm[['locus_id_id', 'jpm']], left_on='locus_id', right_on='locus_id_id', suffixes=(
        'locustable_', 'bjtable_'))[['gene_name', 'jpm']].groupby('gene_name').sum().reset_index()
    locus_cj_bj_jpm_df = locus_cj_jpm_df.merge(
        locus_bj_jpm_df, left_on='gene_name', right_on='gene_name', suffixes=('_cj', '_bj'), how='outer').dropna()
    gene_level_bj_cj_jpm = {
        "gene_name": locus_cj_bj_jpm_df["gene_name"].to_list(),
        "jpm_cj": np.log2(locus_cj_bj_jpm_df["jpm_cj"]).to_list(),
        "jpm_bj": np.log2(locus_cj_bj_jpm_df["jpm_bj"]).to_list(),
        "text": "Gene " + locus_cj_bj_jpm_df['gene_name'].map(str)
    }

    # Boxplot for TPM grouped by is_circrna_host
    locus_bj_df = locus_df[['locus_id']].merge(bj_df[['locus_id_id', 'coord_id']], left_on='locus_id', right_on='locus_id_id').groupby(
        'locus_id_id').agg({'coord_id': pd.Series.nunique}).reset_index()
    locus_bj_df = locus_bj_df.groupby('locus_id_id').sum().reset_index()
    locus_bj_df['circRNA_host'] = True
    locus_bj_tpm_df = locusexpression_df[['locus_id_id', 'tpm']].merge(
        locus_bj_df[['locus_id_id', 'circRNA_host']], left_on='locus_id_id', right_on='locus_id_id', how='outer').fillna(False)
    locus_bj_tpm_data = {
        'circRNA_host': np.log2(locus_bj_tpm_df[['tpm', 'circRNA_host']][locus_bj_tpm_df.circRNA_host]['tpm']+1).tolist(),
        'not_circRNA_host': np.log2(locus_bj_tpm_df[['tpm', 'circRNA_host']][locus_bj_tpm_df.circRNA_host == False]['tpm']+1).tolist()
    }

    # Top X structure level abundance
    top_x_structure_df = pd.merge(locus_df[['locus_id', 'stable_id', 'gene_name']], bj_df[[
        'locus_id_id', 'tpm', 'jpm', 'abundance_ratio', 'coord_id', 'gc_perc', 'raw_count', 'n_methods']], left_on='locus_id', right_on='locus_id_id')
    top_x_structure_df[['abundance_ratio']
                       ] = top_x_structure_df[['abundance_ratio']]*100
    top_x_structure_df = top_x_structure_df.round(3)
    top_x_structure_data = top_x_structure_df.to_dict(orient='records')

    # Top X gene level abundance
    locusexpression_max_tpm = locusexpression_df[['locus_id_id', 'tpm']].groupby(
        'locus_id_id').max().reset_index()
    locus_tpm = locus_df[['gene_name', 'locus_id', 'stable_id']].merge(
        locusexpression_max_tpm, left_on='locus_id', right_on='locus_id_id')
    bj_circrna_count_df = bj_df[['locus_id_id', 'coord_id', 'abundance_ratio']].groupby(
        'locus_id_id').agg({'coord_id': pd.Series.nunique, 'abundance_ratio': 'sum'}).reset_index()
    top_x_gene_level_abundance_df = locus_tpm.merge(bj_circrna_count_df, left_on='locus_id', right_on='locus_id_id', how='outer').fillna(
        0)[['gene_name', 'abundance_ratio', 'tpm', 'coord_id']].groupby('gene_name').agg({'abundance_ratio': 'sum', 'tpm': 'max', 'coord_id': 'sum'}).reset_index()
    top_x_gene_level_abundance_df[[
        'abundance_ratio']] = top_x_gene_level_abundance_df[['abundance_ratio']]*100
    top_x_gene_level_abundance_df = top_x_gene_level_abundance_df.round(3)
    top_x_gene_level_abundance_df = top_x_gene_level_abundance_df.merge(
        gene_df[['gene_name', 'stable_id', 'biotype', 'description']], left_on='gene_name', right_on='stable_id', suffixes=('_locus', ''))
    top_x_gene_level_abundance = top_x_gene_level_abundance_df.to_dict(
        orient='records')

    data = {'species': species.scientific_name,
            'assembly': assembly.assembly_name,
            'sankey': sankey,
            'gene_level_ar_sum': gene_level_ar_sum,
            'jpm_boxplot': jpm_boxplot,
            'gene_level_bj_cj_jpm': gene_level_bj_cj_jpm,
            'locus_bj_tpm_data': locus_bj_tpm_data,
            'top_x_structure_data': top_x_structure_data,
            'top_x_gene_level_abundance': top_x_gene_level_abundance
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

    # List of all locus for a given assembly
    locus_query = 'select locus_id, seq_region_name, seq_region_start, seq_region_end, gene_name from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select seq_region_name, abundance_ratio, coord_id, locus_id_id from core_backsplicejunction where locus_id_id in ({})'.format(
        str(locus_df["locus_id"].to_list())[1:-1])
    bj_df = pd.read_sql_query(bj_query, connection)

    chromosomes = bj_df[['seq_region_name', 'abundance_ratio']].groupby(
        ['seq_region_name']).mean().nlargest(10, ['abundance_ratio']).index.tolist()

    genes = locus_df[['seq_region_name', 'seq_region_start', 'seq_region_end', 'gene_name', 'locus_id']].merge(bj_df[['locus_id_id', 'coord_id']], left_on='locus_id', right_on='locus_id_id').groupby(
        ['gene_name', 'seq_region_name', 'seq_region_start', 'seq_region_end']).agg({'coord_id': 'count'}).reset_index().sort_values(by='coord_id', ascending=False)[['gene_name', 'seq_region_name', 'seq_region_start', 'seq_region_end']].to_dict('records')

    data = {
        'chromosomes': chromosomes,
        'genes': genes
    }

    return Response(data=data, status=status.HTTP_200_OK)


def export_species_view_list(request, species_id, assembly_id):
    """
    View to generate the export file in the species view
    """
    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        raise Http404('No species with the given id.')

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        raise Http404('No assembly with the given id under the given species.')

    # List of all locus for a given assembly
    locus_query = 'select * from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select * from core_backsplicejunction where locus_id_id in ({})'.format(
        str(locus_df["locus_id"].to_list())[1:-1])
    bj_df = pd.read_sql_query(bj_query, connection)

    # List of analysis combined with source field from sample table for given assembly
    sample_analysis_query = 'select sample_id, source, analysis_id from core_sample inner join core_analysis on core_sample.sample_id=core_analysis.sample_id_id where species_id_id={} and assembly_id_id={}'.format(
        species.taxon_id, assembly.assembly_id)
    sample_analysis_df = pd.read_sql_query(sample_analysis_query, connection)

    # Add tissues (source) column in the backsplice table
    bj_df = pd.merge(bj_df, sample_analysis_df,
                     left_on='analysis_id_id', right_on='analysis_id')

    # Get all the required get parameters
    chromosomes = request.GET.getlist('chromosome[]', [])
    classifications = request.GET.getlist('classification[]', [])
    tissues = request.GET.getlist('tissue[]', [])
    min_tpm = float(request.GET.get('tpm', 0))
    min_n_methods = float(request.GET.get('nMethods', 3))
    req_format = 'bed' if request.GET.get('format', 'csv') == 'bed' else 'csv'
    seperator = '\t' if req_format == 'bed' else ','

    # Modifying columns
    del bj_df['analysis_id']
    del bj_df['analysis_id_id']
    bj_df.rename(columns={'source': 'tissue'}, inplace=True)
    bj_df = bj_df[['seq_region_name', 'seq_region_start', 'seq_region_end', 'coord_id', 'raw_count',
                   'seq_region_strand', 'predicted_exons', 'sample_id', 'tissue', 'tpm', 'jpm', 'abundance_ratio', 'n_methods']]

    # Filtering according to the parameters
    if chromosomes:
        bj_df = bj_df[bj_df['seq_region_name'].isin(chromosomes)]
    if classifications:
        bj_df = bj_df[bj_df['classification'].isin(classifications)]
    if tissues:
        bj_df = bj_df[bj_df['tissue'].isin(tissues)]
    bj_df = bj_df[(bj_df['tpm'] >= min_tpm) & (
        bj_df['n_methods'] >= min_n_methods)]

    # Delete not required column according to the required format
    del bj_df['n_methods']
    if req_format == 'bed':
        bj_df = bj_df[['seq_region_name', 'seq_region_start', 'seq_region_end', 'coord_id', 'raw_count',
                       'seq_region_strand']]

    filename = 'ECIRCDB-'+str(species.scientific_name) + '-' + str(assembly.assembly_name) + \
        '-' + str(assembly.assembly_accession) + \
        '-' + str(datetime.datetime.now())

    response = HttpResponse(content_type='text/{}'.format(req_format))
    response['Content-Disposition'] = 'attachment; filename={}.{}'.format(
        filename, req_format)
    bj_df.to_csv(path_or_buf=response, sep=seperator, index=False)

    return response


def export_sample_view_list(request, species_id, assembly_id, sample_id):
    """
    View to generate the export file in the species view
    """
    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        raise Http404('No species with the given id.')

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        raise Http404('No assembly with the given id under the given species.')

    try:
        sample = species.samples.get(sample_id=sample_id)
    except:
        return Http404('No sample with the given id under the given species.')

    # List of the analysis related to the givem sample
    analysis_query = 'select analysis_id from core_analysis where sample_id_id={}'.format(
        sample.sample_id)
    analysis_df = pd.read_sql_query(analysis_query, connection)

    # List of all backsplice junctions related to the provided sample
    bj_query = 'select * from core_backsplicejunction where analysis_id_id in ({})'.format(
        str(analysis_df['analysis_id'].tolist())[1:-1])
    bj_df = pd.read_sql_query(bj_query, connection)

    # Add tissues (source) column in the backsplice table
    bj_df['tissue'] = sample.source
    bj_df['sample_id'] = sample.sample_id

    # Get the parameters
    req_format = 'bed' if request.GET.get('format', 'csv') == 'bed' else 'csv'
    seperator = '\t' if req_format == 'bed' else ','

    # Modifying columns
    del bj_df['analysis_id_id']
    if req_format == 'bed':
        bj_df = bj_df[['seq_region_name', 'seq_region_start', 'seq_region_end', 'coord_id', 'raw_count',
                       'seq_region_strand']]
    else:
        bj_df = bj_df[['seq_region_name', 'seq_region_start', 'seq_region_end', 'coord_id', 'raw_count',
                       'seq_region_strand', 'predicted_exons', 'sample_id', 'tissue', 'tpm', 'jpm', 'abundance_ratio']]

    filename = 'ECIRCDB-'+str(species.scientific_name) + '-' + str(assembly.assembly_name) + '-' + str(
        sample.accession) + '-' + str(assembly.assembly_accession) + '-' + str(datetime.datetime.now())

    response = HttpResponse(content_type='text/{}'.format(req_format))
    response['Content-Disposition'] = 'attachment; filename={}.{}'.format(
        filename, req_format)
    bj_df.to_csv(path_or_buf=response, sep=seperator, index=False)

    return response


@api_view(['GET'])
@cache_page(caching_time)
def get_genome(request, species_id, assembly_id):
    """
    View to get genome for the species
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

    genome = {}
    scientific_name = species.scientific_name
    scientific_name = scientific_name.lower().split(' ')
    scientific_name = '_'.join(scientific_name)

    url = 'https://rest.ensembl.org/info/assembly/{}?bands=1'.format(
        scientific_name)
    headers = {'Content-Type': 'application/json'}

    r = requests.get(url, headers=headers)
    for region in r.json()['top_level_region']:
        genome.setdefault(region['name'],
                          {'size': region['length'], 'bands': []})
        if region['coord_system'] == 'chromosome':
            if 'bands' in region:
                for band in region['bands']:
                    genome[band['seq_region_name']]['bands'].append({
                        'id': band['id'],
                        'start': band['start'],
                        'end': band['end'],
                        'type': band['stain']
                    })

    return Response(data=genome, status=status.HTTP_200_OK)


def gene_track_bed(request, species_id, assembly_id, position):
    """
    View to generate bed file to include gene tracks to the genoverse
    """

    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        raise Http404('No species with the given id.')

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        raise Http404('No assembly with the given id under the given species.')

    # Get location using url parameters
    chromosome = position.split(':')[0]
    start = position.split(':')[1].split('-')[0]
    end = position.split(':')[1].split('-')[1]

    # List of all locus for a given assembly
    locus_query = 'select browser_string from core_locus where assembly_id_id={} and seq_region_name={} and seq_region_start>={} and seq_region_end<={};'.format(
        assembly.assembly_id, chromosome, start, end)
    locus_df = pd.read_sql_query(locus_query, connection)

    def convert_to_josn(x):
        x = x.split('\t')
        data = {
            'chromosome': x[0],
            'start': int(x[1]),
            'end': int(x[2]),
            'name': 'ENSRNOT00000066509.4',
            'score': int(x[4]),
            'strand': x[5],
            'thickStart': int(x[1]) + 1,
            'thickEnd': int(x[2]) - 1,
            'itemRgb': '255,0,0',
            'blockCount': int(x[9]),
            'blockSizes': x[10] + ',',
            'blockStarts': x[11] + ',',
        }
        return data

    browser_string_array = locus_df['browser_string'].unique()
    browser_string_array = ['chr'+x for x in browser_string_array]
    browser_string_series = pd.Series(browser_string_array)
    browser_string_json_series = browser_string_series.apply(
        lambda x: convert_to_josn(x))
    browser_string_json = browser_string_json_series.tolist()
    browser_string_df = pd.DataFrame(browser_string_json)
    cols = ['chromosome', 'start', 'end', 'name', 'score', 'strand',
            'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    browser_string_df = browser_string_df[cols]

    response = HttpResponse(content_type='text/bed')
    browser_string_df.to_csv(
        path_or_buf=response, sep='\t', index=False, index_label=False, header=False)

    return response


def circrna_track_bed(request, species_id, assembly_id, position):
    """
    View to generate bed file to include circRNA tracks to the genoverse
    """

    try:
        species = Species.objects.get(taxon_id=species_id)
    except:
        raise Http404('No species with the given id.')

    try:
        assembly = species.assemblies.get(assembly_id=assembly_id)
    except:
        raise Http404('No assembly with the given id under the given species.')

    # Get location using url parameters
    chromosome = position.split(':')[0]
    start = position.split(':')[1].split('-')[0]
    end = position.split(':')[1].split('-')[1]

    # List of all locus for a given assembly
    locus_query = 'select * from core_locus where assembly_id_id={};'.format(
        assembly.assembly_id)
    locus_df = pd.read_sql_query(locus_query, connection)

    # List of all backsplice junctions related to given assembly
    bj_query = 'select * from core_backsplicejunction where locus_id_id in ({}) and seq_region_name={} and seq_region_start>={} and seq_region_end<={}'.format(
        str(locus_df["locus_id"].to_list())[1:-1], chromosome, start, end)
    bj_df = pd.read_sql_query(bj_query, connection)

    def convert_to_josn(x):
        x = x.split('\t')
        data = {
            'chromosome': x[0],
            'start': int(x[1]),
            'end': int(x[2]),
            'name': 'ENSRNOT00000066509.4',
            'score': int(x[4]),
            'strand': x[5],
            'thickStart': int(x[1]) + 1,
            'thickEnd': int(x[2]) - 1,
            'itemRgb': '52,0,255',
            'blockCount': int(x[9]),
            'blockSizes': x[10] + ',',
            'blockStarts': x[11] + ',',
        }
        return data

    browser_string_array = bj_df['browser_string'].unique()
    browser_string_array = ['chr'+x for x in browser_string_array]
    browser_string_series = pd.Series(browser_string_array)
    browser_string_json_series = browser_string_series.apply(
        lambda x: convert_to_josn(x))
    browser_string_json = browser_string_json_series.tolist()
    browser_string_df = pd.DataFrame(browser_string_json)
    cols = ['chromosome', 'start', 'end', 'name', 'score', 'strand',
            'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    browser_string_df = browser_string_df[cols]

    response = HttpResponse(content_type='text/bed')
    browser_string_df.to_csv(
        path_or_buf=response, sep='\t', index=False, index_label=False, header=False)

    return response
