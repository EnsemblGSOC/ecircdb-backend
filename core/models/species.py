from django.db import models


class Species(models.Model):
    """
    This model holds general information about a species
    """

    taxon_id = models.AutoField(
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
