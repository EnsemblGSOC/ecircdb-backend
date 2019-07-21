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

    data = {'species': species.scientific_name,
            'assembly': assembly.assembly_name,
            'count_total_samples': count_total_samples,
            'count_distinct_tissues': count_distinct_tissues,
            'count_backsplice_junctions': count_backsplice_junctions}

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

    return Response(data={'species': [species.scientific_name, species.description], 'assembly': assembly.assembly_name}, status=status.HTTP_200_OK)


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
