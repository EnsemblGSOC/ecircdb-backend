from django.db import models

from core.models import Assembly


class Locus(models.Model):
    """
    This model holds the information about the Locus
    """

    locus_id = models.AutoField(
        primary_key=True
    )
    browser_string = models.TextField(
        blank=True
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
    gene_name = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    is_circrna_host = models.BooleanField(
        blank=True,
        null=True
    )
    nexons = models.IntegerField(
        blank=True,
        null=True
    )
    source = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    spliced_size = models.IntegerField(
        blank=True,
        null=True
    )
    stable_id = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    total_backsplice_reads = models.IntegerField(
        blank=True,
        null=True
    )
    total_splice_reads = models.IntegerField(
        blank=True,
        null=True
    )
    assembly_id = models.ForeignKey(
        Assembly,
        on_delete=models.CASCADE,
        related_name='loci'
    )

    class Meta():
        verbose_name_plural = 'Loci'

    def __str__(self):
        return str(self.locus_id)
