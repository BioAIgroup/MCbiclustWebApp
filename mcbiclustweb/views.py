from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect, HttpResponse
from django.urls import reverse
from django.templatetags.static import static
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.views.generic import View
import os, shutil

from mcbiclustweb.models import Profile, Analysis, Biclusters

from .forms import RegisterForm, CreateAnalysisForm

from .tasks import *

from django.contrib import messages 

class RegisterFormView(View):
    form_class = RegisterForm
    template_name = 'mcbiclustweb/register.html'

    # new user signing up
    def get(self, request):
        if request.user.is_authenticated:
            return redirect ('mcbiclustweb:index')

        form = self.form_class(None)
        return render(request, self.template_name, {'form': form})

    # process signup form data
    def post(self, request):
        form = self.form_class(request.POST)

        if form.is_valid():
            user = form.save(commit=False)
            
            # cleaned data
            username = form.cleaned_data['username']
            password = form.cleaned_data['password1']
            user.set_password(password)
            user.save()

            # return User objects
            user = authenticate(username=username, password=password)

            if user is not None:
                if user.is_active:
                    login(request, user)
                    return redirect('mcbiclustweb:index')

        return render(request, self.template_name, {'form': form})

class IndexView(View):
    form_class = CreateAnalysisForm
    template_name = 'mcbiclustweb/index.html'

    def get(self, request):
        if request.user.is_authenticated:
            # form = self.form_class(None)
            profile = Profile.objects.get(user=self.request.user)
            a = Analysis.objects.order_by('id').filter(user=profile)
            return render(request, self.template_name, {'analyses': a})
        else:
            return redirect('mcbiclustweb:login')

    def post(self, request):
        if request.user.is_authenticated:
            form = self.form_class(request.POST)
            print(form)
            if form.is_valid():
                analysis = form.save(commit=False)
                analysis.user = Profile.objects.get(user=self.request.user)
                analysis.status = "1. Preprocessing"
                analysis.save()
                analysis.gem = request.FILES['gem']
                analysis.save()
                # try:
                #     shutil.rmtree(os.path.join(settings.MEDIA_ROOT, 'analyses/user_{0}/None'.format(analysis.user.id)))
                # except:
                #     pass
                preprocess.delay(analysis.id)

                return redirect('mcbiclustweb:index')
            else:
                return render(request, 'mcbiclustweb/index.html', {'form':form})
        else:
            return redirect('mcbiclustweb:login')

def analysis(request, analysis_id):
    a = Analysis.objects.get(id=analysis_id)
    status = int(a.status.split('.')[0])
    fig_dir_url = os.path.join(settings.MEDIA_URL, *a.gem.name.split("/")[:-1])
    fig_dir_root = os.path.join(settings.MEDIA_ROOT, *a.gem.name.split("/")[:-1])
    
    nbiclusters = 0
    chars = []
    if status >= 11:
        nbiclusters = sum(os.path.isdir(os.path.join(fig_dir_root, i)) for i in os.listdir(fig_dir_root))

        fork_dir = os.path.join(fig_dir_root, "1")
        for i in os.listdir(fork_dir):
            chars.append(i[9:-4])
    
    return render(request, "mcbiclustweb/analysis.html", {'analysis': a, 'status': status, 'fig_dir': fig_dir_url, 'nbiclusters': range(1, nbiclusters + 1), 'chars': chars})

def start(request, analysis_id):
    a = Analysis.objects.get(id=analysis_id)
    a.status = "3. Started analysis"
    a.save()

    seed_size = int(request.POST['seedSize'])
    init_seed = request.POST['initSeed']
    iterations = int(request.POST['iterations'])
    num_runs = int(request.POST['numRuns'])
    geneset_file = request.FILES['geneset']
    geneset = geneset_file.read().decode("utf-8").split(',')

    runFindSeed.delay(analysis_id, seed_size, init_seed, geneset, iterations, num_runs)

    return redirect('mcbiclustweb:analysis', analysis_id=analysis_id)
    
def delete(request, analysis_id):
    a = Analysis.objects.get(id=analysis_id)
    user_id = a.user.id
    a.delete()

    try:
        shutil.rmtree(os.path.join(settings.MEDIA_ROOT, 'analyses/user_{0}/{1}'.format(user_id, analysis_id)))
    except:
        pass

    return redirect('mcbiclustweb:index')

def run(request):
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
    return HttpResponse(seed)

