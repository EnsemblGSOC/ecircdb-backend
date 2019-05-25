from django.urls import path
from rest_framework.urlpatterns import format_suffix_patterns

from core import views

urlpatterns = [
    path('species/', views.SpeciesList.as_view()),
    path('species/<int:pk>', views.SpeciesDetail.as_view())
]

# urlpatterns = format_suffix_patterns(urlpatterns)
