from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect, HttpResponse
from django.urls import reverse
from django.templatetags.static import static
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.views.generic import View
import os

from .forms import UserForm

import rpy2.robjects as robjects

class UserFormView(View):
    form_class = UserForm
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
            password = form.cleaned_data['password']
            user.set_password(password)
            user.save()

            # return User objects
            user = authenticate(username=username, password=password)

            if user is not None:
                if user.is_active:
                    login(request, user)
                    return redirect('mcbiclustweb:index')

        return render(request, self.template_name, {'form': form})
        


def index(request):
    if request.user.is_authenticated:
        return render(request, 'mcbiclustweb/index.html')
    else:
        return redirect('mcbiclustweb:login')

def results(request):
    if len(request.FILES) != 0:
        analysis_dir = os.path.join(settings.BASE_DIR, 'analysis/', str(request.user.id))
        try:
            os.makedirs(analysis_dir)
        except:
            pass
            
        gem_file = request.FILES['gem']
        dest_file = os.path.join(analysis_dir, gem_file.name)

        with open(dest_file, 'wb+') as destination:
            for chunk in gem_file.chunks():
                destination.write(chunk)

        gem = "abc"
        seed = "abc"
        return render(request, 'mcbiclustweb/results.html', {'gem': gem, 'seed': seed})
    else:
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

