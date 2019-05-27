from rest_framework import filters, generics, permissions

from core.models import Species
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
