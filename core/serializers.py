from rest_framework import serializers

from core.models import Species


class SpeciesSerializer(serializers.ModelSerializer):
    """
    Serializers for Species object
    """

    class Meta:
        model = Species
        fields = ("id", "name", "thumbnail", "description")
