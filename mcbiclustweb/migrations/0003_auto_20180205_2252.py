# Generated by Django 2.0 on 2018-02-05 22:52

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('mcbiclustweb', '0002_auto_20180205_2218'),
    ]

    operations = [
        migrations.RenameField(
            model_name='analysis',
            old_name='gene_expr_mat',
            new_name='gem',
        ),
    ]