from celery import shared_task
from django.conf import settings
import os
import random
import numpy as np
from math import ceil

from mcbiclustweb.models import Profile, Analysis

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.rinterface as ri

@shared_task
def preprocess(analysis_id):
    a = Analysis.objects.get(id=analysis_id)
    # Get GEM directory
    gem_dir = os.path.join(settings.MEDIA_ROOT, a.gem.name)
    # Directory to store processed data
    store_dir = os.path.join(settings.MEDIA_ROOT, 'analyses/user_{0}/{1}'.format(a.user.id, a.id))

    ri.initr()

    # Import libraries
    base = importr('base')
    geoquery = importr('GEOquery')

    # Check if file starts with !Series_title, otherwise getGEO never stops
    with open(gem_dir) as f:
        first_line = f.readline()
    if not first_line[:14].startswith("!Series_title\t"):
        a.status = "-2. Preprocessing failed: invalid gene expression matrix format"
        a.save()
        return "invalid gene expression matrix format"

    # Get GEO series matrix and extract GEM
    try:
        gsm = geoquery.getGEO(filename=gem_dir, getGPL=False)
    except:
        a.status = "-2. Preprocessing failed: invalid gene expression matrix format"
        a.save()
        return "invalid gene expression matrix format"

    try:
        gem = gsm.slots['assayData']['exprs']
        # Remove any genes with NAs
        row_keep = ro.IntVector(np.argwhere(np.array(ro.r.rowSums(ro.r['is.na'](gem))) == 0) + 1)
        gem = gem.rx(row_keep, True)
        # Remove any genes with all 0s
        row_keep = ro.IntVector(np.nonzero(np.array(ro.r.rowSums(gem)))[0] + 1)
        gem = gem.rx(row_keep, True)
        # Write to CSV
        ro.r['write.table'](gem, file=os.path.join(store_dir, 'gem.csv'))
    except:
        a.status = "-2. Preprocessing failed: invalid gene expression matrix format"
        a.save()
        return "invalid gene expression matrix format"

    try:
        # Get pheno data
        pheno_data = gsm.slots['phenoData']
        pheno_data = pheno_data.slots['data']

        # Extract explicitly defined characteristics
        char_index = ro.r['!'](ro.r.grepl('characteristics|date', ro.r.names(pheno_data)))
        char = pheno_data.rx(True, char_index)
        gene_name = ro.r.rownames(char)
        char = ro.r['data.frame'](ro.r.lapply(char, ro.r['as.character']), stringsAsFactors=False)
        char = ro.r.cbind(char, **{'gene.name': gene_name}, stringsAsFactors=False)

        # Data cleaning
        for i in range(char.nrow):
            for j in range(char.ncol - 1):
                char_val = char.rx(i + 1, j + 1)[0].strip()
                if char_val == "":
                    char.rx[i + 1, j + 1] = "unknown"
                elif char_val == "None" or char_val == "NONE" or char_val == "none":
                    char.rx[i + 1, j + 1] = "none"
                else:
                    char.rx[i + 1, j + 1] = char_val

        # Remove columns where more than 80% of unique values have less than 5 occurrences
        col_keep = np.repeat([True], char.ncol)
        for i in range(char.ncol - 1):
            unique_count = ro.r.table(char.rx(True, i + 1))
            unique_length = ro.r.length(unique_count)[0]
            if unique_length == 1:
                col_keep[i] = False
                continue
            count = 0
            for j in range(unique_length):
                if ro.r.names(unique_count).rx(j + 1)[0] == "unknown":
                    unique_length = unique_length - 1
                elif unique_count.rx(j + 1)[0] < 5:
                    count = count + 1
            if count / unique_length > 0.8:
                col_keep[i] = False

        char = char.rx(True,ro.BoolVector(col_keep))

        # Change all values with less than 5 occurences to "Other"
        for i in range(char.ncol - 1):
            rare_char = np.where(np.array(ro.r.table(char.rx2(i + 1))) < 5)[0] + 1
            rare_char = ro.r.names(ro.r.table(char.rx2(i + 1))).rx(ro.IntVector(rare_char))
            if ro.r.length(rare_char)[0] == 1:
                rare_char = ro.IntVector([])
            if ro.r.length(rare_char)[0] != 0:
                for j in range(char.nrow):
                    if str(char.rx(j + 1, i + 1)[0]) in np.array(rare_char):
                        print(np.array(rare_char))
                        char.rx[j + 1, i + 1] = "Other"

        # Write to CSV
        ro.r['write.table'](char, file=os.path.join(store_dir, 'characteristics.csv'))
    except:
        a.char_ok = False
        a.save()

    a.status = "2. Ready for analysis"
    a.save()

    return "success"


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
    dplyr = importr('dplyr', on_conflict="warn")
    gtools = importr("gtools")
    stringr = importr("stringr")

    gem = ro.r['read.csv'](os.path.join(fig_dir, 'gem.csv'), header=True, sep=" ")

    # Check that gene expression matrix contains the genes in gene set of interest 
    for x in geneset:
        if x not in np.array(ro.r.rownames(gem)):
            print(x)
            a.status = "-1. Failed: geneset of interest contains genes that are either not in the series matrix or have NA or 0 as value"
            a.save()
            return "geneset of interest contains genes that are either not in the series matrix or have NA or 0 as value"

    gem_sub = gem.rx(ro.StrVector(geneset), True)

    # Checks if the seed size user inputed is bigger than the number of samples in GEM
    if gem.ncol < seed_size:
        a.status = "-1. Failed: sample seed size bigger than sample number"
        a.save()
        return "sample seed size bigger than number of samples"
    # Checks if one of the initial seeds specified is bigger than the number of samples in GEM
    init_seed_specified = False
    if init_seed != "":
        init_seed_specified = True
        temp = np.array(ro.StrVector(init_seed.split(',')))
        print(temp)
        for s in temp:
            print(s)
            if s not in np.array(ro.r.colnames(gem)):
                a.status = "-1. Failed: one initial seed sample is not in gene expression matrix"
                a.save()
                return "one initial seed sample is not in gene expression matrix"
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
    a.status = "4. Started analysis: found seeds"
    a.save()

    # Calculate correlation vector for each run
    multi_corvec = ro.ListVector({})
    for i in range(num_runs):
        params = {'gem.part': gem_sub, 'gem.all': gem, 'seed': multi_seed.rx2(i + 1), 'splits': 10}
        corvec = mcbiclust.CVEval(**params)
        multi_corvec.rx2[i + 1] = corvec
    print("Calculated correlation vector.")
    a.status = "5. Started analysis: calculated correlation vector"
    a.save()

    # Turn correlation vectors to correlation matrix
    len_a = ro.r.length(multi_corvec.rx2(1))[0]
    len_b = ro.r.length(multi_corvec)[0]
    multi_cormat=ro.r.matrix(0.0, len_a, len_b)
    for i in range(num_runs):
        multi_cormat.rx[True, i + 1] = multi_corvec.rx2(i + 1)
    multi_corvec = None
    print("Created correlation matrix.")
    a.status = "6. Started analysis: created correlation matrix"
    a.save()

    # Plot correlation heatmap
    try:
        cormat1 = ro.r.abs(ro.r.cor(multi_cormat, use='complete.obs'))
        cordist = ro.r("""
            function(c){
                as.dist(1 - abs(c))
            }
        """)
        grdevices.png(file=os.path.join(fig_dir, "cor_heatmap.png"), width=1400, height=875, pointsize=25)
        cormat_heat = gplots.heatmap_2(cormat1, trace="none", distfun=cordist)
        grdevices.dev_off()
    except:
        a.status = "-3. Failed: error plotting correlation heatmap, your gene expression matrix may not be suitable for analysis"
        a.save()
        return "error plotting correlation heatmap, your gene expression matrix may not be suitable for analysis"
    print("Plotted heatmap")
    a.status = "7. Started analysis: plotted heatmap"
    a.save()

    # Find clusters and plot them
    try:
        grdevices.png(file=os.path.join(fig_dir, "sil_clust%02d.png"), width=1400, height=875, pointsize=25)
        params = {'cor.vec.mat': multi_cormat, 'max.clusters': 20, 'plots': True, 'rand.vec': False}
        multi_clust_group = mcbiclust.SilhouetteClustGroups(**params)
        grdevices.dev_off()
    except:
        a.status = "-4. Failed: error finding distinct biclusters, your gene expression matrix may not be suitable for analysis"
        a.save()
        return "error finding distinct biclusters, your gene expression matrix may not be suitable for analysis"
    print("Plotted silhouette")
    a.status = "8. Started analysis: plotted silhouette"
    a.save()

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
    print("Plotted CVPlot")
    a.status = "9. Started analysis: plotted CVPlot"
    a.save()

    # Gene set enrichment
    # corvec_gsea = ro.ListVector({})
    # for i in range(ro.r.length(average_corvec)[0]):
    #     params = {'gene.names': gene_names, 'gene.values': average_corvec.rx2(i + 1), 'sig.rate': 0.05}
    #     corvec_gsea.rx2[i + 1] = mcbiclust.GOEnrichmentAnalysis(**params)
    # for i in range(ro.r.length(corvec_gsea)[0]):
    #     filename = "gsea" + "%02d" % (i + 1,) + ".csv"
    #     ro.r['write.table'](corvec_gsea.rx2(i + 1),file=os.path.join(fig_dir, filename))
    # print("Calculated gene set enrichment")
    # a.status = "10. Started analysis: calculated gene set enrichment"
    # a.save()

    # Sample sorting
    multi_samp_sort = ro.ListVector({})
    params = {'gem': gem, 'av.corvec': average_corvec, 'top.genes.num': 750, 'groups': multi_clust_group, 'initial.seeds': multi_seed}
    multi_prep = mcbiclust.MultiSampleSortPrep(**params)
    for i in range(ro.r.length(multi_prep.rx2(1))[0]):
        print(i)
        params = {'gem': multi_prep.rx2(1).rx2(i + 1), 'seed': multi_prep.rx2(2).rx2(i + 1)}
        multi_samp_sort.rx2[i + 1] = mcbiclust.SampleSort(**params)
    print("Sample sorting finished")
    a.status = "10. Started analysis: sample sorting finished"
    a.save()

    # Calculate PC1 values and threshold new biclusters
    multi_pc1_vec = ro.ListVector({})
    multi_bic = ro.ListVector({})
    try:
        for i in range(ro.r.length(multi_prep.rx2(1))[0]):
            # Calculate PC1
            params = {'top.gem': multi_prep.rx2(1).rx2(i + 1), 'seed.sort': multi_samp_sort.rx2(i + 1), 'n': min(ceil(ro.r.ncol(gem)[0]/20), 10)}
            multi_pc1_vec.rx2[i + 1] = mcbiclust.PC1VecFun(**params)
            # Threshold biclusters
            params = {'cor.vec': average_corvec.rx2(i + 1), 'sort.order': multi_samp_sort.rx2(i + 1), 'pc1': multi_pc1_vec.rx2(i + 1), 'samp.sig': 0.05}
            multi_bic.rx2[i + 1] = mcbiclust.ThresholdBic(**params)
            # Align PC1
            params = {'gem': gem, 'pc1': multi_pc1_vec.rx2(i + 1), 'sort.order': multi_samp_sort.rx2(i + 1), 'cor.vec': average_corvec.rx2(i + 1), 'bic': multi_bic.rx2(i + 1)}
            multi_pc1_vec.rx2[i + 1] = mcbiclust.PC1Align(**params)
    except:
        a.status = "-5. Failed: error extending distinct biclusters, please try to restart the analysis. If this error is shown repeatedly, your gene expression matrix may not be suitable for analysis"
        a.save()
        return "error extending distinct biclusters, please try to restart the analysis. If this error is shown repeatedly, your gene expression matrix may not be suitable for analysis"

    multi_df_args = ro.ListVector({})
    multi_df_args.rx2['gene.name'] = ro.r.colnames(gem)
    for i in range(ro.r.length(multi_samp_sort)[0]):
        multi_df_args.rx2["Bic" + str(i + 1) + ".order"] = ro.r.order(multi_samp_sort.rx2(i + 1))
        multi_df_args.rx2["Bic" + str(i + 1) + ".PC1"] = multi_pc1_vec.rx2(i + 1).rx(ro.r.order(multi_samp_sort.rx2(i + 1)))
    multi_df = ro.r['data.frame'](multi_df_args)

    if a.char_ok == False:
        a.status = "12. Analysis completed"
        a.save()
        return "success"

    char = ro.r['read.csv'](os.path.join(fig_dir,'characteristics.csv'), header=True, sep=" ", stringsAsFactors=True)
    multi_df_char = dplyr.inner_join(multi_df,char,by="gene.name")

    start = 1 + ro.r.length(multi_samp_sort)[0] * 2
    end = multi_df_char.ncol
    for b in range(ro.r.length(multi_samp_sort)[0]):
        for i in range(start, end):
            temp_name = ro.r.colnames(multi_df_char).rx(i + 1)[0]
            
            legend = gtools.mixedsort(ro.r.names(ro.r.table(multi_df_char.rx(True, i + 1))))
            legend_title = gtools.mixedsort(ro.r.names(ro.r.table(multi_df_char.rx(True, i + 1))))
            for j in range(ro.r.length(legend_title)[0]):
                legend_title.rx[j + 1] = stringr.str_wrap(legend_title.rx(j + 1),20)
            
            ro.globalenv['multi_df_char'] = multi_df_char
            forkplot = ro.r('ggplot(multi_df_char, aes(Bic%i.order,Bic%i.PC1)) + geom_point(aes(colour=%s)) + ylab("Bic%i PC1") + labs(colour="%s") + scale_color_discrete(breaks=%s,labels=%s)'%(b + 1, b + 1, temp_name, b + 1, temp_name,legend.r_repr(), legend_title.r_repr()))
            del(ro.globalenv['multi_df_char'])
            
            try:
                os.makedirs(os.path.join(fig_dir,str(b + 1)))
            except:
                pass
            
            ggplot2.ggsave("forkplot_%s.png"%(temp_name.replace(".", "_")), plot=forkplot, device='png', path=os.path.join(fig_dir,str(b + 1)), scale=1.5, width=12, height=7.4, units="cm")
    
    print("Plotted forks")
    a.status = "11. Started analysis: plotted forks"
    a.save()

    a.status = "12. Analysis completed"
    a.save()
    return "success"

