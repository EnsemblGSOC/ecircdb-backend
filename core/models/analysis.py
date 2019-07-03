from django.db import models

from core.models import Assembly, Sample


class Analysis(models.Model):
    """
    This model holds the information about analysis and relate assembly and 
    sample
    """

    analysis_id = models.AutoField(
        primary_key=True,
    )
    run_date = models.CharField(
        max_length=255
    )
    logic_name = models.CharField(
        max_length=255
    )
    parameters = models.CharField(
        max_length=103,
        blank=True,
        null=True
    )
    assembly_id = models.ForeignKey(
        Assembly,
        on_delete=models.CASCADE,
        related_name='analyses'
    )
    sample_id = models.ForeignKey(
        Sample,
        on_delete=models.CASCADE,
        related_name='analyses'
    )

    class Meta():
        verbose_name_plural = 'Analyses'

    def __str__(self):
        return self.logic_name
