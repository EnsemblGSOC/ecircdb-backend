from rest_framework import filters, generics, permissions

from core.models import Species
from core.serializers import SpeciesSerializer


class SpeciesList(generics.ListAPIView):
    """
    View to list and search of species
    """

    filter_backends = (filters.SearchFilter,)
    pagination_class = None
    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer
    search_fields = ['name']


class SpeciesDetail(generics.RetrieveAPIView):
    """
    View to get details of a species
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer
