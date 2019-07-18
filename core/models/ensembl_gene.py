from django.db import models

from core.models import Assembly


class EnsemblGene(models.Model):
    """
    This model holds the information about ensembl gene
    """

    gene_id = models.AutoField(
        primary_key=True
    )
    coord_id = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    seq_region_end = models.IntegerField(
        blank=True,
        null=True
    )
    seq_region_name = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    seq_region_start = models.IntegerField(
        blank=True,
        null=True
    )
    seq_region_strand = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    gene_name = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    description = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    is_circrna_host = models.BooleanField(
        blank=True,
        null=True
    )
    biotype = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    stable_id = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    assembly_id = models.ForeignKey(
        Assembly,
        on_delete=models.CASCADE,
        related_name='assemblies'
    )

    def __str__(self):
        return self.seq_region_name
