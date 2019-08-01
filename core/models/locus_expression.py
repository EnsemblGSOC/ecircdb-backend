from django.db import models

from core.models import Analysis, Locus


class LocusExpression(models.Model):
    """
    This model holds the information about the Locus expression
    """

    expression_id = models.AutoField(
        primary_key=True
    )
    rpkm = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    rpkm_external = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    rpkm_internal = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    rpkm_ratio = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    raw_count = models.IntegerField(
        blank=True,
        null=True
    )
    tpm = models.DecimalField(
        max_digits=19,
        decimal_places=6,
        blank=True,
        null=True
    )
    analysis_id = models.ForeignKey(
        Analysis,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='locus_expressions'
    )
    locus_id = models.ForeignKey(
        Locus,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name='locus_expressions'
    )

    def __str__(self):
        return self.expression_id
