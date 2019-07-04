from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns

from core import views

urlpatterns = [
    path('species/', views.SpeciesList.as_view()),
    path('species/<int:pk>/', views.SpeciesDetail.as_view()),
    path('species/samples/<int:species>/', views.SampleList.as_view()),
    path('samples/<int:pk>/', views.SampleDetail.as_view()),
    path('species_view_stats/<int:species_id>/<int:assembly_id>/',
         views.species_view_stats),
    path('sample_view_stats/<int:species_id>/<int:assembly_id>/<int:sample_id>/',
         views.sample_view_stats)
]

# urlpatterns = format_suffix_patterns(urlpatterns)
