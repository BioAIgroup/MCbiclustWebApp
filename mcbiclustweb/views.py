from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect, HttpResponse
from django.urls import reverse
from django.templatetags.static import static
import rpy2.robjects as robjects

def index(request):
    return render(request, 'mcbiclustweb/index.html')

def results(request):
    gem = request.POST['gem']
    seed = gem
    return render(request, 'mcbiclustweb/results.html', {'gem': gem, 'seed': seed})

def run(request):
    r = robjects.r
    script_dir = static('mcbiclustweb/scripts/myscript.R')
    r.source('/var/www/static/mcbiclustweb/scripts/myscript.R')
    seed = r['CCLE.seed']
    return HttpResponse(seed)

