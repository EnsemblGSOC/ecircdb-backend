from django.db import models


class Species(models.Model):
    """
    This model holds general information about a species
    """

    taxon_id = models.AutoField(
        primary_key=True
    )
    scientific_name = models.CharField(
        max_length=127
    )
    common_name = models.CharField(
        max_length=127,
        null=True,
        blank=True
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
        ordering = ('scientific_name', )

    def __str__(self):
        return self.scientific_name
