from django.db import models


class Species(models.Model):
    """
    This model holds general information about a spieces
    """

    taxon_id = models.IntegerField(
        primary_key=True
    )
    name = models.CharField(
        max_length=127
    )
    thumbnail = models.FileField(
        upload_to='species_thumbnails/',
        null=True,
        blank=True
    )
    description = models.TextField(
        blank=True
    )

    class Meta:
        verbose_name_plural = "species"
        ordering = ('name', )

    def __str__(self):
        return self.name


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
    backspliced_reads = models.IntegerField()
    chimeric_reads = models.IntegerField()
    circrna_count = models.IntegerField()
    description = models.TextField()
    fastq_path = models.CharField(
        max_length=255
    )
    ftp_path = models.CharField(
        max_length=255
    )
    instrument = models.CharField(
        max_length=103
    )
    library_size = models.IntegerField()
    library_stratergy = models.CharField(
        max_length=255
    )
    gtag_cjunks_reads = models.IntegerField()
    mapped_reads = models.IntegerField()
    read_length = models.IntegerField()
    fraction = models.CharField(
        max_length=103
    )
    submitter = models.CharField(
        max_length=103
    )
    source = models.CharField(
        max_length=55
    )
    total_splice_reads = models.IntegerField()
    source_id = models.ForeignKey(
        Species,
        on_delete=models.CASCADE,
        related_name='samples'
    )

    def __str__(self):
        return self.source+'-'+self.submitter
