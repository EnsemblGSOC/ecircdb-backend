from django.db import models

from core.models import Species


class Sample(models.Model):
    """
    This model holds the information about a sample
    """

    sample_id = models.AutoField(
        primary_key=True
    )
    accession = models.CharField(
        max_length=55
    )
    backspliced_reads = models.IntegerField(
        blank=True,
        null=True
    )
    chimeric_reads = models.IntegerField(
        blank=True,
        null=True
    )
    circrna_count = models.IntegerField(
        blank=True,
        null=True
    )
    description = models.TextField(
        blank=True,
        null=True
    )
    fastqc_path = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    bigwig_path = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    ftp_path = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    instrument = models.CharField(
        max_length=103,
        blank=True,
        null=True
    )
    library_size = models.IntegerField(
        blank=True,
        null=True
    )
    library_strategy = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    project = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    gtag_cjuncs_reads = models.IntegerField(
        blank=True,
        null=True
    )
    mapped_reads = models.IntegerField(
        blank=True,
        null=True
    )
    read_length = models.IntegerField(
        blank=True,
        null=True
    )
    fraction = models.CharField(
        max_length=103,
        blank=True,
        null=True
    )
    submitter = models.CharField(
        max_length=103,
        blank=True,
        null=True
    )
    source = models.CharField(
        max_length=55,
        blank=True,
        null=True
    )
    total_spliced_reads = models.IntegerField(
        blank=True,
        null=True
    )
    TIN = models.IntegerField(
        blank=True,
        null=True
    )
    species_id = models.ForeignKey(
        Species,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='samples'
    )

    def __str__(self):
        return self.accession
