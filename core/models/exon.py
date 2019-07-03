from django.db import models

from core.models import Locus


class Exon(models.Model):
    """
    This model holds the information about exon
    """

    coord_id = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    genomic_size = models.IntegerField(
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
    rank = models.IntegerField(
        blank=True,
        null=True
    )
    stable_id = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    locus_id = models.ForeignKey(
        Locus,
        on_delete=models.CASCADE,
        related_name='exons'
    )

    def __str__(self):
        return self.locus_id
