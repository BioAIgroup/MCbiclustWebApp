from django.urls import path

from django.contrib.auth import views as auth_views
from . import views

app_name = 'mcbiclustweb'
urlpatterns = [
    path('', views.index, name='index'),
    path('register/', views.UserFormView.as_view(), name='register'),
    path('login/', auth_views.LoginView.as_view(template_name='mcbiclustweb/login.html', redirect_authenticated_user=True), name='login'),
    path('logout/', auth_views.LogoutView.as_view(template_name='mcbiclustweb/logout.html'), name='logout'),
    path('results/', views.results, name="results"),
    path('run/', views.run, name="run"),
]
