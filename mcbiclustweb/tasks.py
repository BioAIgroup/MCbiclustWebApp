from celery import shared_task
from django.conf import settings
import rpy2.robjects as robjects
import os

from mcbiclustweb.models import Profile, Analysis, Biclusters

@shared_task
def runFindSeed(analysis_id):
    r = robjects.r
    script_dir = os.path.join(settings.STATIC_ROOT, 'mcbiclustweb/scripts/myscript.R')
    
    r.library('MCbiclust')
    r.library('ggplot2')
    r.library('gplots')
    r.library('dplyr')
    r.library('gProfileR')
    r.library('MASS')
    r.library('devtools')

    r.data('CCLE_small')
    r.data('Mitochondrial_genes')
    r.source(script_dir)
    seed = r['CCLE.seed']

    a = Analysis.objects.get(id=analysis_id)
    
    biclusters = ','.join(map(str, list(seed))) 
    b = Biclusters(analysis=a, biclusters=biclusters)
    b.save()

    a.status = "completed"
    a.save()
    

    