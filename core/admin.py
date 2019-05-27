from django.contrib import admin

from core.models import Species, Sample

admin.site.register((Species, Sample))
