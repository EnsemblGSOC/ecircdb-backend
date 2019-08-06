from django.db import models

from core.models import Analysis, Locus


class BackspliceJunction(models.Model):
    """
    This model holds the information about a Backsplice expression
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
        related_name='backsplice_junctions'
    )
    locus_id = models.ForeignKey(
        Locus,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='backsplice_junctions'
    )
    classification = models.CharField(
        max_length=255,
        blank=True,
        null=True,
    )
    in_platelets = models.BooleanField(
        blank=True,
        null=True
    )
    is_published = models.BooleanField(
        blank=True,
        null=True
    )
    junction_repeat_coverage = models.IntegerField(
        blank=True,
        null=True
    )
    n_methods = models.IntegerField(
        blank=True,
        null=True
    )
    spliced_size = models.IntegerField(
        blank=True,
        null=True
    )
    splice_type = models.CharField(
        max_length=255,
        blank=True,
        null=True,
    )
    splice_3_support = models.IntegerField(
        blank=True,
        null=True
    )
    splice_5_support = models.IntegerField(
        blank=True,
        null=True
    )
    splice_site_count = models.IntegerField(
        blank=True,
        null=True
    )
    tpm = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    transcript_count = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    gc_perc = models.IntegerField(
        blank=True,
        null=True
    )
    abundance_ratio = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )

    def __str__(self):
        return str(self.junction_id)
