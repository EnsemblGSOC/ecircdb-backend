from django.db import models

from core.models import Analysis, Locus


class CanonicalJunction(models.Model):
    """
    This model holds the information about a Canonical Junction
    """

    junction_id = models.AutoField(
        primary_key=True
    )
    browser_string = models.TextField(
        blank=True,
        null=True
    )
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
    jpm = models.IntegerField(
        blank=True,
        null=True
    )
    predicted_exons = models.CharField(
        max_length=255,
        blank=True,
        null=True,
    )
    raw_count = models.IntegerField(
        blank=True,
        null=True
    )
    splice_signal = models.IntegerField(
        blank=True,
        null=True
    )
    analysis_id = models.ForeignKey(
        Analysis,
        on_delete=models.CASCADE,
        related_name='canonical_junctions'
    )
    locus_id = models.ForeignKey(
        Locus,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='canonical_junctions'
    )
    splice_type = models.CharField(
        max_length=255,
        blank=True,
        null=True,
    )

    def __str__(self):
        return self.junction_id
