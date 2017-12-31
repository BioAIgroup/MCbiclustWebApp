from django.urls import path

from . import views

app_name = 'mcbiclustweb'
urlpatterns = [
    path('', views.index, name='index'),
    path('results/', views.results, name="results"),
    path('run/', views.run, name="run"),
]
