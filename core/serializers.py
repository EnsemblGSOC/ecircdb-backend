from rest_framework import serializers

from core.models import Species, Sample


class SpeciesListSerializer(serializers.ModelSerializer):
    """
    Serializers for List view of Species object
    """

    class Meta:
        model = Species
        fields = "__all__"


class SampleSerializer(serializers.ModelSerializer):
    """
    Serializers for Sample Object
    """

    class Meta:
        model = Sample
        fields = "__all__"


class SpeciesDetailSerializer(serializers.ModelSerializer):
    """
    Serializers for Detail view of Species object
    """

    samples = SampleSerializer(many=True, read_only=True)

    class Meta:
        model = Species
        fields = "__all__"
