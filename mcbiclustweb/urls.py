from django.urls import path

from django.contrib.auth import views as auth_views
from . import views

app_name = 'mcbiclustweb'
urlpatterns = [
    path('', views.IndexView.as_view(), name='index'),
    path('analysis/<int:analysis_id>', views.analysis, name="analysis"),
    path('register/', views.RegisterFormView.as_view(), name='register'),
    path('login/', auth_views.LoginView.as_view(template_name='mcbiclustweb/login.html', redirect_authenticated_user=True), name='login'),
    path('logout/', auth_views.LogoutView.as_view(template_name='mcbiclustweb/logout.html'), name='logout'),
    path('analysis/<int:analysis_id>/start', views.start, name="start"),
    path('analysis/<int:analysis_id>/delete', views.delete, name="delete"),
    path('run/', views.run, name="run"),
]