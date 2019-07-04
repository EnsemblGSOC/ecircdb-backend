from rest_framework import serializers

from core.models import Assembly, Species, Sample


class SpeciesListSerializer(serializers.ModelSerializer):
    """
    Serializers for List view of Species object
    """

    class Meta:
        model = Species
        fields = "__all__"


class AssemblySerializer(serializers.ModelSerializer):
    """
    Serializers for Assembly Object
    """

    class Meta:
        model = Assembly
        exclude = ('species_id',)


class SampleListSerializer(serializers.ModelSerializer):
    """
    Serializers for Sample Object
    """

    class Meta:
        model = Sample
        fields = ('sample_id', 'accession', 'source', 'description',)


class SampleDetailsSerializer(serializers.ModelSerializer):
    """
    Serializers for Sample Object
    """

    class Meta:
        model = Sample
        fields = "__all__"


class SpeciesDetailSerializer(serializers.ModelSerializer):
    """
    Serializers for Detail view of Species object having list of assemblies
    """

    assemblies = AssemblySerializer(many=True, read_only=True)

    class Meta:
        model = Species
        fields = "__all__"
