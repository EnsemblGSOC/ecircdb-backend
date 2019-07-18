from django.db import models


class HibernateSequences(models.Model):
    """
    This model holds the information about hibernate sequences
    """

    sequence_name = models.CharField(
        max_length=255,
        blank=True,
        null=True
    )
    sequence_next_hi_value = models.IntegerField(
        blank=True,
        null=True
    )

    def __str__(self):
        return self.sequence_name
