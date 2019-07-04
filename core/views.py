from rest_framework import filters, generics, permissions, status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from core.models import Species, Assembly
from core.serializers import SpeciesDetailSerializer, SpeciesListSerializer


class SpeciesList(generics.ListAPIView):
    """
    View to list and search of species
    """

    filter_backends = (filters.SearchFilter,)
    ordering = ('name', )
    pagination_class = None
    queryset = Species.objects.all().order_by('name')
    serializer_class = SpeciesListSerializer
    search_fields = ['name']


class SpeciesDetail(generics.RetrieveAPIView):
    """
    View to get details of a species
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesDetailSerializer


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

    return Response(data={'species': [species.name, species.description], 'assembly': assembly.assembly_name}, status=status.HTTP_200_OK)


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

    return Response(data={'species': [species.name, species.description], 'assembly': assembly.assembly_name}, status=status.HTTP_200_OK)
