from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns

from core import views

urlpatterns = [
    path('species/', views.SpeciesList.as_view()),
    path('species/<int:pk>/', views.SpeciesDetail.as_view()),
    path('species/samples/<int:species_id>/<int:assembly_id>/',
         views.SampleList.as_view()),
    path('samples/<int:pk>/', views.SampleDetail.as_view()),
    path('species_view_stats/<int:species_id>/<int:assembly_id>/',
         views.species_view_stats),
    path('sample_view_stats/<int:species_id>/<int:assembly_id>/<int:sample_id>/',
         views.sample_view_stats),
    path('location_view_stats/<int:species_id>/<int:assembly_id>/',
         views.location_view_stats),
    path('export_species_view_list/<int:species_id>/<int:assembly_id>/',
         views.export_species_view_list),
    path('export_sample_view_list/<int:species_id>/<int:assembly_id>/<int:sample_id>/',
         views.export_sample_view_list),
    path('get_genome/<int:species_id>/<int:assembly_id>/', views.get_genome),
    path('gene_track_bed/<int:species_id>/<int:assembly_id>/<str:position>/',
         views.gene_track_bed),
    path('circrna_track_bed/<int:species_id>/<int:assembly_id>/<str:position>/',
         views.circrna_track_bed)
]

# urlpatterns = format_suffix_patterns(urlpatterns)
