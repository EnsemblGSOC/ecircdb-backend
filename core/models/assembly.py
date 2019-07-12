from django.db import models

from core.models import Species


class Assembly(models.Model):
    """
    This model holds the information about an assembly
    """

    assembly_id = models.AutoField(
        primary_key=True,
    )
    assembly_accession = models.CharField(
        max_length=255
    )
    assembly_name = models.CharField(
        max_length=255
    )
    assembly_genoverse_genome = models.CharField(
        max_length=255
    )
    species_id = models.ForeignKey(
        Species,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='assemblies'
    )

    class Meta():
        verbose_name_plural = 'Assemblies'

    def __str__(self):
        return self.assembly_name
