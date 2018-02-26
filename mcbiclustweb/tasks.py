from celery import shared_task
from django.conf import settings
import os
import random
import numpy as np

from mcbiclustweb.models import Profile, Analysis, Biclusters

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.rinterface as ri

@shared_task
def runFindSeed(analysis_id, seed_size, init_seed, geneset, iterations, num_runs):
    a = Analysis.objects.get(id=analysis_id)
    # Get GEM directory
    gem_dir = os.path.join(settings.MEDIA_ROOT, a.gem.name)
    # Directory to store figures
    fig_dir = os.path.join(settings.MEDIA_ROOT, 'analyses/user_{0}/{1}'.format(a.user.id, a.id))

    ri.initr()

    # Imports libraries
    mcbiclust = importr('MCbiclust')
    geoquery = importr('GEOquery')
    gplots = importr('gplots')
    ggplot2 = importr('ggplot2')
    grdevices = importr('grDevices')

    # Get GEO series matrix and extract GEM
    gsm = geoquery.getGEO(filename=gem_dir, getGPL=False)
    gem = gsm.slots['assayData']['exprs']
    # Remove any genes with NAs
    # is_na = ro.r['is.na']
    row_keep = ro.IntVector(np.argwhere(np.array(ro.r.rowSums(ro.r['is.na'](gem)))==0)+1)
    gem = gem.rx(row_keep, True)
    # Remove anay genes with all 0s
    row_keep = ro.IntVector(np.nonzero(np.array(ro.r.rowSums(gem)))[0]+1)
    gem = gem.rx(row_keep, True)
    # Keep the genes in the geneset of interest
    for x in geneset:
        if x not in np.array(ro.r.rownames(gem)):
            print(x)
            a.status = "failed: geneset of interest contains genes that are either not in the series matrix or have NA or 0 as value"
            a.save()
            return "geneset of interest contains genes that are either not in the series matrix or have NA or 0 as value"

    gem_sub = gem.rx(ro.StrVector(geneset), True)
    # Choose at most 500 genes randomly to perform FindSeed on
    # gem_num_genes = min(500, gem.nrow)
    # gem_genes = ri.IntSexpVector(random.sample(list(range(1, gem.nrow + 1)), gem_num_genes))
    # gem_samples = ri.IntSexpVector(range(1, gem.ncol + 1))
    # gem_sub = gem.rx(gem_genes, gem_samples)

    # Checks if the seed size user inputed is bigger than the number of samples in GEM
    if gem.ncol < seed_size:
        a.status = "failed: seed size bigger than sample number"
        a.save()
        return "seed size bigger than number of samples"
    # Checks if one of the initial seeds specified is bigger than the number of samples in GEM
    init_seed_specified = False
    if init_seed != "":
        init_seed_specified = True
        temp = np.array(ro.StrVector(init_seed.split(',')))
        print(temp)
        for s in temp:
            print(s)
            if s not in np.array(ro.r.colnames(gem)):
                a.status = "failed: one initial seed is not in gene expression matrix"
                a.save()
                return "one initial seed is not in gene expression matrix"
        init_seed = []
        for s in temp:
            print(np.where(np.array(ro.r.colnames(gem))==s))
            init_seed.append(np.where(np.array(ro.r.colnames(gem))==s)[0][0]+1)
    
    # Find seeds
    multi_seed = ro.ListVector({})
    for i in range(num_runs):
        # Create parameters for FindSeed
        if init_seed_specified == True:
            init_seed = ro.IntVector(init_seed)
            params = {'gem': gem_sub, 'seed.size': seed_size, 'initial.seed': init_seed, 'iterations': iterations}
        else:
            random.seed(i)
            init_seed = ri.IntSexpVector(random.sample(list(range(1, gem.ncol + 1)), seed_size))
            params = {'gem': gem_sub, 'seed.size': seed_size, 'initial.seed': init_seed, 'iterations': iterations}
        # Perform FindSeed
        seed = mcbiclust.FindSeed(**params)
        print(seed)
        multi_seed.rx2[i + 1] = seed
    print("Found seeds.")

    # Calculate correlation vector for each run
    multi_corvec = ro.ListVector({})
    for i in range(num_runs):
        params = {'gem.part': gem_sub, 'gem.all': gem, 'seed': multi_seed.rx2(i + 1), 'splits': 10}
        corvec = mcbiclust.CVEval(**params)
        multi_corvec.rx2[i + 1] = corvec
    print("Calculated correlation vector.")

    # Turn correlation vectors to correlation matrix
    len_a = ro.r.length(multi_corvec.rx2(1))[0]
    len_b = ro.r.length(multi_corvec)[0]
    multi_cormat=ro.r.matrix(0.0, len_a, len_b)
    for i in range(num_runs):
        multi_cormat.rx[True, i + 1] = multi_corvec.rx2(i + 1)
    multi_corvec = None
    print("Created correlation matrix.")

    # Plot correlation heatmap
    cormat1 = ro.r.abs(ro.r.cor(multi_cormat, use='complete.obs'))
    cordist = ro.r("""
        function(c){
            as.dist(1 - abs(c))
        }
    """)
    grdevices.png(file=os.path.join(fig_dir, "cor_heatmap.png"), width=1400, height=875, pointsize=25)
    cormat_heat = gplots.heatmap_2(cormat1, trace="none", distfun=cordist)
    grdevices.dev_off()
    print("Plotted heatmap")

    # Find clusters and plot them
    grdevices.png(file=os.path.join(fig_dir, "sil_clust%02d.png"), width=1400, height=875, pointsize=25)
    params = {'cor.vec.mat': multi_cormat, 'max.clusters': 20, 'plots': True, 'rand.vec': False}
    multi_clust_group = mcbiclust.SilhouetteClustGroups(**params)
    grdevices.dev_off()
    print("Plotted silhouette")

    # CVPlot
    gene_names = ro.r.rownames(gem)
    average_corvec = ro.ListVector({})
    for i in range(ro.r.length(multi_clust_group)[0]):
        x = ro.IntVector(np.argwhere(np.array(multi_clust_group.rx2(i + 1)) == 1) + 1)
        if ro.r.length(x)[0] == 1:
            average_corvec.rx2[i + 1] = multi_cormat.rx(True,x)
        else:
            average_corvec.rx2[i + 1] = ro.r.rowMeans(multi_cormat.rx(True,x))
    geneset_loc = []
    for gene in geneset:
        geneset_loc.append(np.where(np.array(ro.r.rownames(gem))==gene)[0][0])
    params = {'cv.df': ro.r['as.data.frame'](average_corvec), 'geneset.loc': ro.IntVector(geneset_loc), 'geneset.name': 'Interest', 'alpha1': 0.1}
    cvplot = mcbiclust.CVPlot(**params)
    ggplot2.ggsave("cvplot.png", plot=cvplot, device='png', path=fig_dir, width=4.7, height=2.9)

    # biclusters = ','.join(map(str, list(seed))) 
    # b = Biclusters(analysis=a, biclusters=biclusters)
    # b.save()

    # a.status = "completed"
    # a.save()
    return "success"

    