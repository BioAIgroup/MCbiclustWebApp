from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver

class Profile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    user_type = models.CharField(max_length=100, default="normal")

    @receiver(post_save, sender=User)
    def create_user_profile(sender, instance, created, **kwargs):
        if created:
            Profile.objects.create(user=instance)

# Create your models here.
class Analysis(models.Model):
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=300)
    gene_expr_mat = models.CharField(max_length=100) # File Link
    user = models.ForeignKey(Profile, on_delete=models.CASCADE)
    status = models.CharField(max_length=20)

class Biclusters(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE)
    biclusters = models.CharField(max_length=100)