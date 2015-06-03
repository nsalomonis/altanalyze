### clustering.py
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#import matplotlib
#matplotlib.use('GTKAgg')

import sys, os, string
command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>0 and '--' in command_args: commandLine=True
else: commandLine=False

import traceback
try:
    import math
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        try: matplotlib.use('TkAgg')
        except Exception: pass
        if commandLine==False:
            try: matplotlib.rcParams['backend'] = 'TkAgg'
            except Exception: pass
        try:
            import matplotlib.pyplot as pylab
            import matplotlib.colors as mc
            import matplotlib.mlab as mlab
            from matplotlib import mpl
            import matplotlib.ticker as tic
            from matplotlib.patches import Circle
            from mpl_toolkits.mplot3d import Axes3D
            mpl.rcParams['axes.linewidth'] = 0.5
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.rcParams['font.family'] = 'sans-serif'
            mpl.rcParams['font.sans-serif'] = 'Arial'
        except Exception:
            print traceback.format_exc()
            print 'Matplotlib support not enabled'
        import scipy
        from scipy.sparse.csgraph import _validation
        from scipy.linalg import svd
        import scipy.cluster.hierarchy as sch
        import scipy.spatial.distance as dist
        try: import numpy; np = numpy
        except Exception:
            print 'Numpy import error...'
            print traceback.format_exc()
        try:
            import igraph.vendor.texttable
        except ImportError: pass
        try:
            from sklearn.decomposition import PCA, FastICA
        except Exception: pass
        #pylab.ion() # closes Tk window after show - could be nice to include
except Exception:
    print traceback.format_exc()
    pass

import time
import unique
import statistics
import os
import export
import webbrowser
import warnings
import UI

try:
    warnings.simplefilter("ignore", numpy.ComplexWarning)
    warnings.simplefilter("ignore", DeprecationWarning) ### Annoying depreciation warnings (occurs in sch somewhere)
    #This shouldn't be needed in python 2.7 which suppresses DeprecationWarning - Larsson
except Exception: None

import WikiPathways_webservice

try:
    import fastcluster as fc
    #print 'Using fastcluster instead of scipy hierarchical cluster'
    #fc = sch
except Exception:
    #print 'Using scipy insteady of fastcluster (not installed)'
    try: fc = sch ### fastcluster uses the same convention names for linkage as sch
    except Exception: print 'Scipy support not present...'
    
def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin
    
def heatmap(x, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient,
            dataset_name, display=False, contrast=None, allowAxisCompression=True,Normalize=True):
    print "Performing hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)
    show_color_bars = True ### Currently, the color bars don't exactly reflect the dendrogram colors
    try: ExportCorreleationMatrix = exportCorreleationMatrix
    except Exception: ExportCorreleationMatrix = False

    try: os.mkdir(root_dir) ### May need to create this directory
    except Exception: None
    if display == False:
        pylab.figure() ### Add this to avoid a Tkinter bug after running MarkerFinder (not sure why it is needed) - creates a second empty window when display == True
                
    if row_method == 'hopach' or column_method == 'hopach':
        try:
            from pyper import R
            r = R(use_numpy=True)
            print_out = r('library("hopach")')
            if "Error" in print_out:
                print 'Installing the R package "hopach" in Config/R'
                print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("hopach")')
                if "Error" in print_out: print 'unable to download the package "hopach"'; forceError
        except Exception,e:
            print 'Failed to install hopach'
            row_method = 'average'; column_method = 'average'
        if len(column_header)==2: column_method = 'average'
        if len(row_header)==2: row_method = 'average'
    """
    Prototype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre
    Scaling the color gradient so that zero is white:
    http://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set
    Other cluster methods:
    http://stackoverflow.com/questions/9362304/how-to-get-centroids-from-scipys-hierarchical-agglomerative-clustering
    
    x is a m by n ndarray, m observations, n genes
    """
    ### Perform the associated clustering by HOPACH via PYPE or Rpy to R
    print row_method, column_method
    if row_method == 'hopach' or column_method == 'hopach':
        try:
    
            """ HOPACH is a clustering method implemented in R that builds a hierarchical tree of clusters by recursively
            partitioning a data set, while ordering and possibly collapsing clusters at each level:
            http://www.bioconductor.org/packages/release/bioc/html/hopach.html
            """
            
            import R_interface
    
            if row_method == 'hopach' and column_method == 'hopach': cluster_method = 'both'
            elif row_method == 'hopach': cluster_method = 'gene'
            else: cluster_method = 'array'
            
            if row_metric == 'cosine': metric_gene = "euclid"
            elif row_metric == 'euclidean': metric_gene = "cosangle"
            elif row_metric == 'correlation': metric_gene = "cor"
            else: metric_gene = "cosangle"
            
            if column_metric == 'cosine': metric_array = "euclid"
            elif column_metric == 'euclidean': metric_array = "cosangle"
            elif column_metric == 'correlation': metric_array = "cor"
            else: metric_array = "euclid"
            
            ### Returned are the row_order and column_order in the Scipy clustering output format
            newFilename, Z1, Z2 = R_interface.remoteHopach(inputFilename,cluster_method,metric_gene,metric_array)
            if newFilename != inputFilename:
                ### If there were duplicates, re-import the matrix data for the cleaned up filename
                try:
                    matrix, column_header, row_header, dataset_name, group_db = importData(newFilename,Normalize=normalize,reverseOrder=False)
                except Exception:
                    matrix, column_header, row_header, dataset_name, group_db = importData(newFilename)
                x = numpy.array(matrix)
        except Exception:
            row_method = 'average'; column_method = 'average'
            print traceback.format_exc()
            print 'hopach failed... continue with an alternative method'
            
    n = len(x[0]); m = len(x)
    if color_gradient == 'red_white_blue':
        cmap=pylab.cm.bwr
    if color_gradient == 'red_black_sky':
        cmap=RedBlackSkyBlue()
    if color_gradient == 'red_black_blue':
        cmap=RedBlackBlue()
    if color_gradient == 'red_black_green':
        cmap=RedBlackGreen()
    if color_gradient == 'yellow_black_blue':
        cmap=YellowBlackBlue()
    if color_gradient == 'black_yellow_blue':
        cmap=BlackYellowBlue()
    if color_gradient == 'seismic':
        cmap=pylab.cm.seismic
    if color_gradient == 'green_white_purple':
        cmap=pylab.cm.PiYG_r
    if color_gradient == 'coolwarm':
        cmap=pylab.cm.coolwarm
    
    vmin=x.min()
    vmax=x.max()
    vmax = max([vmax,abs(vmin)])
    if Normalize != False:
        vmin = vmax*-1
    elif 'Clustering-Zscores-' in dataset_name:
        vmin = vmax*-1
    default_window_hight = 8.5
    default_window_width = 12
    if len(column_header)>80:
        default_window_width = 14
    if len(column_header)>100:
        default_window_width = 16
    if contrast == None:
        scaling_factor = 2.5 #2.5
    else:
        try: scaling_factor = float(contrast)
        except Exception: scaling_factor = 2.5
    
    #print vmin/scaling_factor
    norm = mpl.colors.Normalize(vmin/scaling_factor, vmax/scaling_factor) ### adjust the max and min to scale these colors by 2.5 (1 scales to the highest change)
    fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    pylab.rcParams['font.size'] = 7.5

    if show_color_bars == False:
        color_bar_w = 0.000001 ### Invisible but not gone (otherwise an error persists)
    else:
        color_bar_w = 0.0125 ### Sufficient size to show

    bigSampleDendrogram = True
    if bigSampleDendrogram == True and row_method==None and column_method != None and allowAxisCompression == True:
        dg2 = 0.30
        dg1 = 0.43
    else: dg2 = 0.1; dg1 = 0.63
    
    try:
        if EliteGeneSets != [''] and EliteGeneSets !=[]:
            matrix_horiz_pos = 0.22
        else:
            matrix_horiz_pos = 0.14
    except Exception:
        matrix_horiz_pos = 0.14

    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.235,matrix_horiz_pos,dg1]   ### The last controls matrix hight, second value controls the position of the matrix relative to the bottom of the view [0.05,0.22,0.2,0.6] 
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
    # axr, placement of row side colorbar    
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w-0.002,0.6] ### second to last controls the width of the side color bar - 0.015 when showing [0.31,0.1,color_bar_w,0.6]
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004
    
    # axc, placement of column side colorbar (3rd value controls the width of the matrix!)
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.6,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing [0.4,0.63,0.5,color_bar_w]
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004
    
    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5] #[0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w
        
    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,dg2] ### last one controls hight of the dendrogram [0.3,0.72,0.6,0.135]
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2 
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.02,0.93,0.17,0.025] ### Last one controls the hight [0.07,0.88,0.18,0.076]
    
    # axcc - placement of the colum colormap legend colormap (distinct map)
    [axcc_x, axcc_y, axcc_w, axcc_h] = [0.02,0.12,0.17,0.025] ### Last one controls the hight [0.07,0.88,0.18,0.076]
    
    # Compute and plot top dendrogram
    if column_method == 'hopach':
        ind2 = numpy.array(Z2['level']) ### from R_interface - hopach root cluster level
    elif column_method != None:
        start_time = time.time()
        #print x;sys.exit()
        d2 = dist.pdist(x.T)
        #print d2
        #import mdistance2
        #d2 = mdistance2.mpdist(x.T)
        #print d2;sys.exit()
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
        if ExportCorreleationMatrix:
            new_matrix=[]
            for i in D2:
                #string.join(map(inverseDist,i),'\t')
                log2_data = map(inverseDist,i)
                avg = statistics.avg(log2_data)
                log2_norm = map(lambda x: x-avg,log2_data)
                new_matrix.append(log2_norm)
            x = numpy.array(new_matrix)
            row_header = column_header
            #sys.exit()
        Y2 = fc.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        #Y2 = sch.fcluster(Y2, 10, criterion = "maxclust")
        try: Z2 = sch.dendrogram(Y2)
        except Exception:
            if column_method == 'average':
                column_metric = 'euclidean'
            else: column_method = 'average'
            Y2 = fc.linkage(D2, method=column_method, metric=column_metric)
            Z2 = sch.dendrogram(Y2)
        #ind2 = sch.fcluster(Y2,0.6*D2.max(), 'distance') ### get the correlations
        #ind2 = sch.fcluster(Y2,0.2*D2.max(), 'maxclust') ### alternative method biased based on number of clusters to obtain (like K-means)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
    # Compute and plot left dendrogram
    if row_method == 'hopach':
        ind1 = numpy.array(Z1['level']) ### from R_interface - hopach root cluster level
    elif row_method != None:
        start_time = time.time()
        d1 = dist.pdist(x)
        D1 = dist.squareform(d1)  # full matrix
        # postion = [left(x), bottom(y), width, height]
        #print D1;sys.exit()
        Y1 = fc.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        no_plot=False ### Indicates that we want to show the dendrogram
        try:
            if runGOElite: no_plot = True
            else:
                ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False - this window conflicts with GO-Elite labels
        except Exception:
            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
        try: Z1 = sch.dendrogram(Y1, orientation='right',no_plot=no_plot) ### This is where plotting occurs
        except Exception:
            row_method = 'average'
            try:
                Y1 = fc.linkage(D1, method=row_method, metric=row_metric)
                Z1 = sch.dendrogram(Y1, orientation='right',no_plot=no_plot)
            except Exception:
                row_method = 'ward'
                Y1 = fc.linkage(D1, method=row_method, metric=row_metric)
                Z1 = sch.dendrogram(Y1, orientation='right',no_plot=no_plot)
        #ind1 = sch.fcluster(Y1,0.6*D1.max(),'distance') ### get the correlations
        #ind1 = sch.fcluster(Y1,0.2*D1.max(),'maxclust')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        if ExportCorreleationMatrix:
            Z1 = sch.dendrogram(Y2, orientation='right')
            Y1 = Y2
            d1 = d2
            D1 = D2
            ind1 = ind2
        try: ax1.set_xticks([]); ax1.set_yticks([]) ### Hides ticks
        except Exception: pass
        
        
        time_diff = str(round(time.time()-start_time,1))
        print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
        

    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = x
    if column_method != None:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        #ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
        ind2 = [ind2[i] for i in idx2] ### replaces the above due to numpy specific windows version issue
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        prior_xt = xt
        xt = xt[idx1,:]   # xt is transformed x
        #ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
        ind1 = [ind1[i] for i in idx1] ### replaces the above due to numpy specific windows version issue
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])
    #axm.set_axis_off() ### Hide border
    #fix_verts(ax1,1)
    #fix_verts(ax2,0)

    ### Adjust the size of the fonts for genes and arrays based on size and character length
    row_fontsize = 5
    column_fontsize = 5
    column_text_max_len = max(map(lambda x: len(x), column_header)) ### Get the maximum length of a column annotation
    if len(row_header)<75:
        row_fontsize = 6.5
        if len(row_header)<50:
            row_fontsize = 8
            if len(row_header)<25:
                row_fontsize = 11
    if len(column_header)<75:
        column_fontsize = 6.5
        if len(column_header)<50:
            column_fontsize = 8
            if len(column_header)<25:
                column_fontsize = 11
                if column_text_max_len < 15:
                    column_fontsize = 15
                elif column_text_max_len > 30:
                    column_fontsize = 6.5
                else:
                    column_fontsize = 10
    
    try:
        if len(justShowTheseIDs)>50:
            column_fontsize = 7
        elif len(justShowTheseIDs)>0:
            column_fontsize = 10
        if len(justShowTheseIDs)>0:
            additional_symbols=[]
            import gene_associations, OBO_import
            try:
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                #symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
            except Exception: gene_to_symbol={}; symbol_to_gene={}
    except Exception: pass

    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            new_row_header.append(row_header[idx1[i]])
        else:
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if column_method != None:
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            new_column_header.append(column_header[i])
            
    dataset_name = string.replace(dataset_name,'Clustering-','')### clean up the name if already a clustered file
    if '-hierarchical' in dataset_name:
        dataset_name = string.split(dataset_name,'-hierarchical')[0]
    filename = 'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name,column_metric,row_metric)

    elite_dir, cdt_file = exportFlatClusterData(root_dir + filename, root_dir, dataset_name, new_row_header,new_column_header,xt,ind1,ind2,display)

    def ViewPNG(png_file_dir):
        if os.name == 'nt':
            try: os.startfile('"'+png_file_dir+'"')
            except Exception:  os.system('open "'+png_file_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')
        
    try:
        if 'monocle' in justShowTheseIDs:
            import R_interface
            R_interface.performMonocleAnalysisFromHeatmap(species,cdt_file[:-3]+'txt',cdt_file[:-3]+'txt')
            png_file_dir = root_dir+'/Monocle/monoclePseudotime.png'
            print png_file_dir
            ViewPNG(png_file_dir)
    except Exception:
        print traceback.format_exc()
    
    cluster_elite_terms={}; ge_fontsize=12; top_genes=[]; proceed=True
    try:
        try:
            if 'driver' in justShowTheseIDs: proceed = False
        except Exception: pass
        if proceed:
            cluster_elite_terms,top_genes = remoteGOElite(elite_dir)
            if cluster_elite_terms['label-size']>40: ge_fontsize = 10
    except Exception: pass  #print traceback.format_exc()
    #exportGOEliteInputs(root_dir,dataset_name,new_row_header,xt,ind1)

    try:
        if len(justShowTheseIDs)<1 and len(top_genes) > 0 and column_fontsize < 9:
            column_fontsize = 10
        if len(justShowTheseIDs)<1:
            additional_symbols=[]
            import gene_associations, OBO_import
            try:
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                #symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
            except Exception: gene_to_symbol={}; symbol_to_gene={}
    except Exception: pass
    
    # Add text
    new_row_header=[]
    new_column_header=[]
    ci=0 ### index of entries in the cluster
    last_cluster=1
    interval = int(float(string.split(str(len(row_header)/40.0),'.')[0]))+1 ### for enrichment term labels with over 100 genes
    increment=interval-2
    if len(row_header)<100: increment = interval-1
    label_pos=-0.03*len(column_header)-.5
    #print label_pos
    try:
        if 'top' in justShowTheseIDs: justShowTheseIDs.remove('top')
        if 'positive' in justShowTheseIDs: justShowTheseIDs.remove('positive')
        if 'amplify' in justShowTheseIDs: justShowTheseIDs.remove('amplify')
        if 'IntraCorrelatedOnly' in justShowTheseIDs: justShowTheseIDs.remove('IntraCorrelatedOnly')
    except Exception:
        pass

    for i in range(x.shape[0]):
        if len(row_header)<40:
            radj = len(row_header)*0.009 ### row offset value to center the vertical position of the row label
        elif len(row_header)<70:
            radj = len(row_header)*0.007 ### row offset value to center the vertical position of the row label
        else:
            radj = len(row_header)*0.005
        cluster = str(ind1[i])
        if cluster != last_cluster:
            ci=0
            increment=0
        last_cluster = cluster
        #print cluster,i,row_header[idx1[i]]
        color = 'black'
        if row_method != None:
            try:
                if row_header[idx1[i]] in justShowTheseIDs:
                    if len(row_header)>len(justShowTheseIDs):
                        color = 'red'
                else: color = 'black'
            except Exception: pass
            if len(row_header)<106: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[idx1[i]],fontsize=row_fontsize, color=color, picker=True)
            new_row_header.append(row_header[idx1[i]])
            new_index = idx1[i]
        else:
            try:
                if row_header[i] in justShowTheseIDs: color = 'red'
                else: color = 'black'
            except Exception: pass
            if len(row_header)<106: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[i],fontsize=row_fontsize, color=color, picker=True) ### When not clustering rows
            new_row_header.append(row_header[i])
            new_index = i ### This is different when clustering rows versus not
        if len(row_header)<106:
            """
            if cluster in cluster_elite_terms:
                try:
                    term = cluster_elite_terms[cluster][ci][1]
                    axm.text(-1.5, i-radj, term,horizontalalignment='right',fontsize=row_fontsize)
                except Exception: pass
            ci+=1
            """
            pass
        else:
            feature_id = row_header[new_index]
            if ':' in feature_id:
                if 'ENS' != feature_id[:3]:
                    feature_id = string.split(feature_id,':')[1]
                else:
                    feature_id = string.split(feature_id,':')[0]
                    try: feature_id = gene_to_symbol[feature_id][0]
                    except Exception: pass
            if ' ' in feature_id and 'ENS' in feature_id:
                feature_id = string.split(feature_id,' ')[1]
            try:
                if feature_id in justShowTheseIDs: color = 'red'
                else: color = 'black'
            except Exception: pass
            try:
                if feature_id in justShowTheseIDs or (len(justShowTheseIDs)<5 and feature_id in top_genes):
                    axm.text(x.shape[1]-0.5, i-radj, '  '+feature_id,fontsize=column_fontsize, color=color,picker=True) ### When not clustering rows
                elif ' ' in row_header[new_index]:
                    symbol = string.split(row_header[new_index], ' ')[-1]
                    if symbol in justShowTheseIDs:
                        axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[new_index],fontsize=column_fontsize, color=color,picker=True)
            except Exception: pass
        
        if cluster in cluster_elite_terms:
                try:
                    increment+=1
                    #print [increment,interval,cluster],cluster_elite_terms[cluster][ci][1];sys.exit()
                    if increment == interval:
                        original_term = cluster_elite_terms[cluster][ci][1]
                        term = original_term
                        if 'GO:' in term:
                            term = string.split(term, '(')[0]
                        if ':WP' in term:
                            term = string.split(term, ':WP')[0]
                        term += ' (c'+str(cluster)+')'
                        try: cluster_elite_terms[term] = cluster_elite_terms[cluster,original_term]  ### store the new term name with the associated genes
                        except Exception: pass
                        axm.text(label_pos, i-radj, term,horizontalalignment='right',fontsize=ge_fontsize, picker=True, color = 'blue', zorder=11)
                        increment=0
                        ci+=1
                except Exception,e: increment=0
    
    def onpick1(event):
        text = event.artist
        print('onpick1 text:', text.get_text())
        if '(c' not in text.get_text():
            webbrowser.open('http://www.genecards.org/cgi-bin/carddisp.pl?gene='+text.get_text())
        elif 'TreeView' in text.get_text():
            try: openTreeView(cdt_file)
            except Exception: print 'Failed to open TreeView'
        else:
            #"""
            import TableViewer
            header = ['Associated Genes']
            tuple_list = []
            
            for gene in cluster_elite_terms[text.get_text()]:
                tuple_list.append([(gene)])
            TableViewer.viewTable(text.get_text(),header,tuple_list) #"""
            
            cluster_prefix = 'c'+string.split(text.get_text(),'(c')[1][:-1]+'-'
            for geneSet in EliteGeneSets:
                if geneSet == 'GeneOntology':
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+'GO'+'.png'
                if geneSet == 'WikiPathways':
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+'local'+'.png'
                elif len(geneSet)>1:
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+geneSet+'.png'
            #try: UI.GUI(root_dir,'ViewPNG',[],png_file_dir)
            #except Exception: print traceback.format_exc()
            
            if os.name == 'nt':
                try: os.startfile('"'+png_file_dir+'"')
                except Exception:  os.system('open "'+png_file_dir+'"')
            elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
            elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')   
            
            #print cluster_elite_terms[text.get_text()]
            
    fig.canvas.mpl_connect('pick_event', onpick1)
            
    for i in range(x.shape[1]):
        adji = i
        ### Controls the vertical position of the column (array) labels
        if len(row_header)<3:
            cadj = len(row_header)*-0.26 ### column offset value
        elif len(row_header)<4:
            cadj = len(row_header)*-0.23 ### column offset value  
        elif len(row_header)<6:
            cadj = len(row_header)*-0.18 ### column offset value     
        elif len(row_header)<10:
            cadj = len(row_header)*-0.08 ### column offset value     
        elif len(row_header)<15:
            cadj = len(row_header)*-0.04 ### column offset value
        elif len(row_header)<20:
            cadj = len(row_header)*-0.05 ### column offset value
        elif len(row_header)<22:
            cadj = len(row_header)*-0.06 ### column offset value
        elif len(row_header)<23:
            cadj = len(row_header)*-0.08 ### column offset value
        elif len(row_header)>200:
            cadj = -2
        else:
            cadj = -0.9
        #cadj = -1
        if len(column_header)>15:
            adji = i-0.1 ### adjust the relative position of the column label horizontally
        if len(column_header)>20:
            adji = i-0.2 ### adjust the relative position of the column label horizontally
        if len(column_header)>25:
            adji = i-0.2 ### adjust the relative position of the column label horizontally
        if len(column_header)>30:
            adji = i-0.25 ### adjust the relative position of the column label horizontally
        if len(column_header)>35:
            adji = i-0.3 ### adjust the relative position of the column label horizontally
        if column_method != None:
            axm.text(adji, cadj, ''+column_header[idx2[i]], rotation=270, verticalalignment="top",fontsize=column_fontsize) # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            axm.text(adji, cadj, ''+column_header[i], rotation=270, verticalalignment="top",fontsize=column_fontsize)
            new_column_header.append(column_header[i])


    # Plot colside colors
    # axc --> axes for column side colorbar

    group_name_list=[]
    ind1_clust,ind2_clust = ind1,ind2
    ind1,ind2,group_name_list,cb_status = updateColorBarData(ind1,ind2,new_column_header,new_row_header,row_method)
    if (column_method != None or 'column' in cb_status) and show_color_bars == True:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
        #cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
        if len(unique.unique(ind2))==2: ### cmap_c is too few colors
            cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF'])
        elif len(unique.unique(ind2))==4: ### cmap_c is too few colors
            cmap_c = mpl.colors.ListedColormap(['#88BF47', '#3D3181', '#FEBC18', '#EE2C3C'])
        elif len(unique.unique(ind2))==5: ### cmap_c is too few colors
            cmap_c = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
        elif len(unique.unique(ind2))==6: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
        elif len(unique.unique(ind2))==7: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
        elif len(unique.unique(ind2))>0: ### cmap_c is too few colors
            cmap_c = pylab.cm.gist_rainbow
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2)) 
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        if 'hopach' == column_method and len(group_name_list)>0:
            axc.set_yticklabels(['','Groups'],fontsize=10)
        else:
            axc.set_yticks([])
        #axc.set_frame_on(False) ### Hide border

        if len(group_name_list)>0: ### Add a group color legend key
            if 'hopach' == column_method: ### allows us to add the second color bar       
                axcd = fig.add_axes([ax2_x, ax2_y, ax2_w, color_bar_w])  # dendrogram coordinates with color_bar_w substituted - can use because dendrogram is not used
                cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
                #cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
                if len(unique.unique(ind2_clust))==2: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#00FF00', '#1E90FF'])
                elif len(unique.unique(ind2_clust))==4: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#3D3181', '#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2_clust))==5: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2_clust))==6: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2_clust))==7: ### cmap_c is too few colors
                    cmap_c = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2_clust))>0: ### cmap_c is too few colors
                    cmap_c = pylab.cm.gist_rainbow
                dc = numpy.array(ind2_clust, dtype=int)
                dc.shape = (1,len(ind2_clust)) 
                im_cd = axcd.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
                #axcd.text(-1,-1,'clusters')
                axcd.set_yticklabels(['','Clusters'],fontsize=10)
                #pylab.yticks(range(1),['HOPACH clusters'])
                axcd.set_xticks([]) ### Hides ticks
                #axcd.set_yticks([])
   
            axd = fig.add_axes([axcc_x, axcc_y, axcc_w, axcc_h])
            group_name_list.sort()
            group_colors = map(lambda x: x[0],group_name_list)
            group_names = map(lambda x: x[1],group_name_list)
            cmap_d = mpl.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
            #cmap_d = mpl.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
            if len(unique.unique(ind2))==2: ### cmap_c is too few colors
                cmap_d = mpl.colors.ListedColormap(['#00FF00', '#1E90FF'])
            elif len(unique.unique(ind2))==4: ### cmap_c is too few colors
                cmap_d = mpl.colors.ListedColormap(['#88BF47', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))==5: ### cmap_c is too few colors
                cmap_d = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))==6: ### cmap_c is too few colors
                cmap_d = mpl.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))==7: ### cmap_c is too few colors
                cmap_d = mpl.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))>0: ### cmap_c is too few colors
                cmap_d = pylab.cm.gist_rainbow
            dc = numpy.array(group_colors, dtype=int)
            dc.shape = (1,len(group_colors)) 
            im_c = axd.matshow(dc, aspect='auto', origin='lower', cmap=cmap_d)
            axd.set_yticks([])
            #axd.set_xticklabels(group_names, rotation=45, ha='left')
            pylab.xticks(range(len(group_names)),group_names,rotation=45,ha='left')
            #cmap_c = mpl.colors.ListedColormap(map(lambda x: GroupDB[x][-1], new_column_header))

    if show_color_bars == False:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        axc.set_frame_on(False)
    
    # Plot rowside colors
    # axr --> axes for row side colorbar
    if (row_method != None or 'row' in cb_status) and show_color_bars == True:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        try: dr = numpy.array(ind1, dtype=int)
        except Exception:
            print ind1;kill
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['#00FF00', '#1E90FF', '#FFFF00', '#FF1493'])
        if len(unique.unique(ind1))>4: ### cmap_r is too few colors
            cmap_r = pylab.cm.gist_rainbow
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])
        #axr.set_frame_on(False) ### Hide border
        
    if show_color_bars == False:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        axr.set_frame_on(False)

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
    #axcb.set_title("colorkey",fontsize=14)
    
    if 'LineageCorrelations' in dataset_name:
        cb.set_label("Lineage Correlation Z Scores",fontsize=11)
    elif 'Heatmap' in root_dir:
        cb.set_label("GO-Elite Z Scores",fontsize=11)
    else:
        cb.set_label("Differential Expression (log2)",fontsize=10)

    ### Add filename label to the heatmap
    if len(dataset_name)>30:fontsize = 10
    else: fontsize = 12.5
    fig.text(0.015, 0.970, dataset_name, fontsize = fontsize)
    
    ### Render and save the graphic
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename, dpi=100) #,dpi=200
    
    fig.text(0.020, 0.070, 'Open heatmap in TreeView (click here)', fontsize = 11.5, picker=True,color = 'red', backgroundcolor='white')
    if 'Outlier' in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Outlier Genes Genes',root_dir+filename])
    elif 'Relative' in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Significant Genes (Relative comparisons)',root_dir+filename])
    elif 'LineageCorrelations' in filename:
        graphic_link.append(['Hierarchical Clustering - Lineage Correlations',root_dir+filename])
    elif 'MarkerGenes' in filename:
        graphic_link.append(['Hierarchical Clustering - MarkerFinder',root_dir+filename])
    elif 'AltExonConfirmed' in filename:
        graphic_link.append(['Hierarchical Clustering - AltExonConfirmed',root_dir+filename])
    elif 'AltExon' in filename:
        graphic_link.append(['Hierarchical Clustering - AltExon',root_dir+filename])
    else:
        graphic_link.append(['Hierarchical Clustering - Significant Genes',root_dir+filename])
    if display:
        proceed=True
        try:
            if 'driver' in justShowTheseIDs:
                proceed = False
        except Exception: pass
        if proceed:
            print 'Exporting:',filename
            try: pylab.show()
            except Exception: None ### when run in headless mode
    fig.clf()
    #fig.close() causes segfault
    #pylab.close() causes segfault

def openTreeView(filename):
    import subprocess
    fn = filepath("AltDatabase/TreeView/TreeView.jar")
    retcode = subprocess.call(['java', "-Xmx500m", '-jar', fn, "-r", filename])

def remoteGOElite(elite_dir):
    mod = 'Ensembl'
    pathway_permutations = 'FisherExactTest'
    filter_method = 'z-score'
    z_threshold = 1.96
    p_val_threshold = 0.05
    change_threshold = 2
    if runGOElite:
        resources_to_analyze = EliteGeneSets
        if 'all' in resources_to_analyze:
            resources_to_analyze = 'all'
        returnPathways = 'no'
        root = None
        import GO_Elite
        
        input_files = dir_list = unique.read_directory(elite_dir) ### Are there any files to analyze?
        if len(input_files)>0 and resources_to_analyze !=['']:
            print '\nBeginning to run GO-Elite analysis on all results' 
            file_dirs = elite_dir,None,elite_dir
            variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,root
            try: GO_Elite.remoteAnalysis(variables,'non-UI Heatmap')
            except Exception: 'GO-Elite failed for:',elite_dir
            try: UI.openDirectory(elite_dir+'/GO-Elite_results')
            except Exception: None
            cluster_elite_terms,top_genes = importGOEliteResults(elite_dir)
            return cluster_elite_terms,top_genes
        else:
            return {},[]
    else:
        return {},[]
    
def importGOEliteResults(elite_dir):
    pruned_results = elite_dir+'/GO-Elite_results/pruned-results_z-score_elite.txt'
    firstLine=True
    cluster_elite_terms={}
    all_term_length=[0]
    for line in open(pruned_results,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        if firstLine:
            firstLine=False
            try: symbol_index = values.index('gene symbols')
            except Exception: symbol_index = None
        else:
            try: symbol_index = values.index('gene symbols')
            except Exception: pass
            try:
                cluster = str(int(float(string.split(values[0][1:],'-')[0])))
                term = values[2]
                all_term_length.append(len(term))
                pval = float(values[9])
                try: cluster_elite_terms[cluster].append([pval,term])
                except Exception: cluster_elite_terms[cluster] = [[pval,term]]
                if symbol_index!=None:
                    symbols = string.split(values[symbol_index],'|')
                    cluster_elite_terms[cluster,term] = symbols
            except Exception,e: pass
    for cluster in cluster_elite_terms:
        cluster_elite_terms[cluster].sort()
    cluster_elite_terms['label-size'] = max(all_term_length)
    
    top_genes = []; count=0
    ranked_genes = elite_dir+'/GO-Elite_results/CompleteResults/ORA_pruned/gene_associations/pruned-gene-ranking.txt'
    for line in open(ranked_genes,'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data,'\t')
        count+=1
        if len(values)>2:
            if values[2]!='Symbol':
                try: top_genes.append((int(values[4]),values[2]))
                except Exception: pass
    top_genes.sort(); top_genes.reverse()
    top_genes = map(lambda x: x[1],top_genes[:21])
    return cluster_elite_terms,top_genes
    
def mergeRotateAroundPointPage(page, page2, rotation, tx, ty):
    from pyPdf import PdfFileWriter, PdfFileReader
    translation = [[1, 0, 0],
                   [0, 1, 0],
                   [-tx,-ty,1]]
    rotation = math.radians(rotation)
    rotating = [[math.cos(rotation), math.sin(rotation),0],
                [-math.sin(rotation),math.cos(rotation), 0],
                [0,                  0,                  1]]
    rtranslation = [[1, 0, 0],
                   [0, 1, 0],
                   [tx,ty,1]]
    ctm = numpy.dot(translation, rotating)
    ctm = numpy.dot(ctm, rtranslation)

    return page.mergeTransformedPage(page2, [ctm[0][0], ctm[0][1],
                                             ctm[1][0], ctm[1][1],
                                             ctm[2][0], ctm[2][1]])

def mergePDFs2(pdf1,pdf2,outPdf):
    from pyPdf import PdfFileWriter, PdfFileReader
    input1 = PdfFileReader(file(pdf1, "rb"))
    page1 = input1.getPage(0)

    input2 = PdfFileReader(file(pdf2, "rb"))
    page2 = input2.getPage(0)
    
    page3 = mergeRotateAroundPointPage(page1, page2, 
                    page1.get('/Rotate') or 0, 
                    page2.mediaBox.getWidth()/2, page2.mediaBox.getWidth()/2)
    output = PdfFileWriter()
    output.addPage(page3)
    outputStream = file(outPdf, "wb")
    output.write(outputStream)
    outputStream.close()    
    
def mergePDFs(pdf1,pdf2,outPdf):
    # http://stackoverflow.com/questions/6041244/how-to-merge-two-landscape-pdf-pages-using-pypdf
    from pyPdf import PdfFileWriter, PdfFileReader

    input1 = PdfFileReader(file(pdf1, "rb"))
    page1 = input1.getPage(0)
    page1.mediaBox.upperRight = (page1.mediaBox.getUpperRight_x(), page1.mediaBox.getUpperRight_y())
    
    input2 = PdfFileReader(file(pdf2, "rb"))
    page2 = input2.getPage(0)
    page2.mediaBox.getLowerLeft_x = (page2.mediaBox.getLowerLeft_x(), page2.mediaBox.getLowerLeft_y())
    # Merge
    page2.mergePage(page1)

    # Output
    output = PdfFileWriter()
    output.addPage(page1)
    outputStream = file(outPdf, "wb")
    output.write(outputStream)
    outputStream.close()    

def merge_horizontal(out_filename, left_filename, right_filename):
    """ Merge the first page of two PDFs side-to-side """
    import pyPdf
    # open the PDF files to be merged
    with open(left_filename) as left_file, open(right_filename) as right_file, open(out_filename, 'w') as output_file:
        left_pdf = pyPdf.PdfFileReader(left_file)
        right_pdf = pyPdf.PdfFileReader(right_file)
        output = pyPdf.PdfFileWriter()

        # get the first page from each pdf
        left_page = left_pdf.pages[0]
        right_page = right_pdf.pages[0]

        # start a new blank page with a size that can fit the merged pages side by side
        page = output.addBlankPage(
            width=left_page.mediaBox.getWidth() + right_page.mediaBox.getWidth(),
            height=max(left_page.mediaBox.getHeight(), right_page.mediaBox.getHeight()),
        )

        # draw the pages on that new page
        page.mergeTranslatedPage(left_page, 0, 0)
        page.mergeTranslatedPage(right_page, left_page.mediaBox.getWidth(), 0)

        # write to file
        output.write(output_file)
        
def inverseDist(value):
    if value == 0: value = 1
    return math.log(value,2)

def getGOEliteExportDir(root_dir,dataset_name):
    if 'AltResults' in root_dir:
        root_dir = string.split(root_dir,'AltResults')[0]
    if 'ExpressionInput' in root_dir:
        root_dir = string.split(root_dir,'ExpressionInput')[0]
    if 'ExpressionOutput' in root_dir:
        root_dir = string.split(root_dir,'ExpressionOutput')[0]
    if 'DataPlots' in root_dir:
        root_dir = string.replace(root_dir,'DataPlots','GO-Elite')
        elite_dir = root_dir
    else:
        elite_dir = root_dir+'/GO-Elite'
    try: os.mkdir(elite_dir)
    except Exception: pass
    return elite_dir+'/clustering/'+dataset_name

def systemCodeCheck(IDs):
    import gene_associations
    id_type_db={}
    for id in IDs:
        id_type = gene_associations.predictIDSourceSimple(id)
        try: id_type_db[id_type]+=1
        except Exception: id_type_db[id_type]=1
    id_type_count=[]
    for i in id_type_db:
        id_type_count.append((id_type_db[i],i))
    id_type_count.sort()
    id_type = id_type_count[-1][-1]
    return id_type

def exportFlatClusterData(filename, root_dir, dataset_name, new_row_header,new_column_header,xt,ind1,ind2,display):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
    filename = string.replace(filename,'.pdf','.txt')
    export_text = export.ExportFile(filename)
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','']+ map(str, ind2),'\t')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)

    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    #new_row_header = new_row_header[::-1]
    #xt = xt[::-1]
    
    try: elite_dir = getGOEliteExportDir(root_dir,dataset_name)
    except Exception: elite_dir = None
    
    elite_columns = string.join(['InputID','SystemCode'])
    try: sy = systemCodeCheck(new_row_header)
    except Exception: sy = None
    
    ### Export each row in the clustered data matrix xt
    i=0
    cluster_db={}
    export_lines = []
    for row in xt:
        id = new_row_header[i]
        if sy == '$En:Sy':
            cluster = 'cluster-'+string.split(id,':')[0]
        elif sy == 'S' and ':' in id:
            cluster = 'cluster-'+string.split(id,':')[0]
        elif sy == 'Sy' and ':' in id:
            cluster = 'cluster-'+string.split(id,':')[0]
        else:
            cluster = 'c'+str(ind1[i])
        try: cluster_db[cluster].append(new_row_header[i])
        except Exception: cluster_db[cluster] = [new_row_header[i]]
        export_lines.append(string.join([new_row_header[i],str(ind1[i])]+map(str, row),'\t')+'\n')
        i+=1
        
    ### Reverse the order of the file
    export_lines.reverse()
    for line in export_lines:
        export_text.write(line)
    export_text.close()


    ### Export GO-Elite input files
    allGenes={}
    for cluster in cluster_db:
        export_elite = export.ExportFile(elite_dir+'/'+cluster+'.txt')
        if sy==None:
            export_elite.write('ID\n')
        else:
            export_elite.write('ID\tSystemCode\n')
        for id in cluster_db[cluster]:
            if sy == '$En:Sy':
                id = string.split(id,':')[1]
                ids = string.split(id,' ')
                if 'ENS' in ids[0]: id = ids[0]
                else: id = ids[-1]
                sc = 'En'
            elif sy == 'Sy' and ':' in id:
                id = string.split(id,':')[1]
                ids = string.split(id,' ')
                sc = 'Sy'
            elif sy == 'En:Sy':
                id = string.split(id,' ')[0]
                sc = 'En'
            elif sy == 'Ae':
                l = string.split(id,':')
                if len(l)==2:
                    id = string.split(id,':')[0] ### Use the Ensembl
                if len(l) == 3:
                    id = string.split(id,':')[1] ### Use the Ensembl
                sc = 'En'
                if ' ' in id:
                    ids = string.split(id,' ')
                    if 'ENS' in ids[-1]: id = ids[-1]
                    else: id = ids[0]
            else:
                sc = sy
            if sy == 'S':
                if ':' in id:
                    id = string.split(id,':')[-1]
                    sc = 'Ae'
            try: export_elite.write(id+'\t'+sc+'\n')
            except Exception: export_elite.write(id+'\n') ### if no System Code known
            allGenes[id]=[]
        export_elite.close()
    try:
        if storeGeneSetName != None:
            if len(storeGeneSetName)>0 and 'driver' not in justShowTheseIDs:
                exportCustomGeneSet(storeGeneSetName,species,allGenes)
                print 'Exported geneset to "StoredGeneSets"'
    except Exception: pass
    
    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    if display:
        try: exportJTV(filename, new_column_header, new_row_header)
        except Exception: pass
    export_cdt = export.ExportFile(filename)
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)
    
    ### Export each row in the clustered data matrix xt
    i=0; cdt_lines=[]
    for row in xt:
        cdt_lines.append(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
        i+=1
        
    ### Reverse the order of the file
    cdt_lines.reverse()
    for line in cdt_lines:
        export_cdt.write(line)
    
    export_cdt.close()
    return elite_dir, filename

def exportJTV(cdt_dir, column_header, row_header):
    ### This is a config file for TreeView
    filename = string.replace(cdt_dir,'.cdt','.jtv')
    export_jtv = export.ExportFile(filename)
    cscale = '3'
    if len(column_header)>100:
        cscale = '1.5'
    if len(column_header)>200:
        cscale = '1.1'
    if len(column_header)>300:
        cscale = '0.6'
    if len(column_header)>400:
        cscale = '0.3'
        
    hscale = '5'
    if len(row_header)< 50:
        hscale = '10'
    if len(row_header)>100:
        hscale = '3'
    if len(row_header)>500:
        hscale = '1'
    if len(row_header)>1000:
        hscale = '0.5'
    config = '<DocumentConfig><UrlExtractor/><ArrayUrlExtractor/><MainView><ColorExtractor>'
    config+= '<ColorSet down="#00FFFF"/></ColorExtractor><ArrayDrawer/><GlobalXMap>'
    config+= '<FixedMap type="Fixed" scale="'+cscale+'"/><FillMap type="Fill"/><NullMap type="Null"/>'
    config+= '</GlobalXMap><GlobalYMap><FixedMap type="Fixed" scale="'+hscale+'"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></GlobalYMap><ZoomXMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></ZoomXMap><ZoomYMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></ZoomYMap><TextView><TextView><GeneSummary/></TextView><TextView>'
    config+= '<GeneSummary/></TextView><TextView><GeneSummary/></TextView></TextView><ArrayNameView>'
    config+= '<ArraySummary included="0"/></ArrayNameView><AtrSummary/><GtrSummary/></MainView></DocumentConfig>'
    export_jtv.write(config)
    
### How to create custom colors - http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

def updateColorBarData(ind1,ind2,column_header,row_header,row_method):
    """ Replace the top-level cluster information with group assignments for color bar coloring (if group data present)"""
    cb_status = 'original'
    group_number_list=[]
    group_name_list=[]
    try: ### Error if GroupDB not recognized as global
        if column_header[0] in GroupDB: ### Thus group assignments exist for column headers
            cb_status = 'column'
            for header in column_header:
                group,color,color_num = GroupDB[header]
                group_number_list.append(color_num) ### will replace ind2
                if (color_num,group) not in group_name_list:
                    group_name_list.append((color_num,group))
            ind2 = group_number_list
        if row_header[0] in GroupDB and row_method == None: ### Thus group assignments exist for row headers
            group_number_list=[]
            if cb_status == 'column': cb_status = 'column-row'
            else: cb_status = 'row'
            for header in row_header:
                group,color,color_num = GroupDB[header]
                group_number_list.append(color_num) ### will replace ind2
            #group_number_list.reverse()
            ind1 = group_number_list
    except Exception: None
    return ind1,ind2,group_name_list,cb_status

def RedBlackSkyBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.9),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackGreen():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    
    my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def YellowBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.8),
                       (0.5, 0.1, 0.0),
                       (1.0, 1.0, 1.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def BlackYellowBlue():
    cdict = {'red':   ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green': ((0.0, 0.0, 0.8),
                       (0.5, 0.1, 0.0),
                       (1.0, 1.0, 1.0)),
    
             'blue':  ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mc.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def importData(filename,Normalize=False,reverseOrder=True,geneFilter=None):
    start_time = time.time()
    fn = filepath(filename)
    matrix=[]
    original_matrix=[]
    row_header=[]
    x=0; inputMax=0; inputMin=100
    filename = string.replace(filename,'\\','/')
    dataset_name = string.split(filename,'/')[-1][:-4]
    if '.cdt' in filename: start = 3
    else: start = 1
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            if '.cdt' in filename: t = [t[0]]+t[3:]
            ### color samples by annotated groups if an expression file
            if 'exp.' in filename and ':' not in data:
                filename = string.replace(filename,'-steady-state.txt','.txt')
                try:
                    import ExpressionBuilder
                    sample_group_db = ExpressionBuilder.simplerGroupImport(filename)
                    new_headers = []
                    for v in t:
                        if v in sample_group_db:
                            v = sample_group_db[v]+':'+v
                        new_headers.append(v)
                    t = new_headers
                except Exception:
                    #print traceback.format_exc()
                    pass
            group_db, column_header = assignGroupColors(t[1:])
            x=1
        elif 'column_clusters-flat' in t:
            start = 2
        elif 'EWEIGHT' in t: pass
        else:
                nullsPresent = False
                #if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                try: s = map(float,t[start:])
                except Exception:
                    nullsPresent=True
                    s=[]
                    for value in t[start:]:
                        try: s.append(float(value))
                        except Exception: s.append(0.000101)
                    #s = numpy.ma.masked_values(s, 0.000101)
                original_matrix.append(s)
                if max(s)>inputMax: inputMax = max(s)
                if min(s)<inputMin: inputMin = min(s)
                #if (abs(max(s)-min(s)))>2:
                if Normalize!=False:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
                        if Normalize=='row mean':
                            #avg = min(s)
                            avg = numpy.mean(s)
                        else: avg = avg = numpy.median(s)
                    if nullsPresent:
                        s=[] ### Needs to be done to zero out the values
                        for value in t[start:]:
                            try: s.append(float(value)-avg)
                            except Exception: s.append(0.000101)
                        #s = numpy.ma.masked_values(s, 0.000101)
                    else:
                        s = map(lambda x: x-avg,s) ### normalize to the mean
                
                if geneFilter==None:
                    matrix.append(s)
                    row_header.append(t[0])
                else:
                    if t[0] in geneFilter:
                        matrix.append(s)
                        row_header.append(t[0])
                x+=1
            
    if inputMax>100: ### Thus, not log values
        print 'Converting values to log2...'
        matrix=[]
        k=0
        if inputMin==0: increment = 1#0.01
        else: increment = 1
        for s in original_matrix:
            if 'counts.' in filename:
                s = map(lambda x: math.log(x+1,2),s)
            else:
                try: s = map(lambda x: math.log(x+increment,2),s)
                except Exception:
                    print filename
                    print Normalize
                    print row_header[k], min(s),max(s); kill
            if Normalize!=False:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
                    if Normalize=='row mean':
                        avg = numpy.average(s)
                    else: avg = avg = numpy.median(s)
                s = map(lambda x: x-avg,s) ### normalize to the mean
            matrix.append(s)
            k+=1
        del original_matrix
        
    #reverseOrder = True ### Cluster order is background (this is a temporary workaround)
    if reverseOrder == True:
        matrix.reverse(); row_header.reverse()
        
    time_diff = str(round(time.time()-start_time,1))
    try:
        print '%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error

    ### Add groups for column pre-clustered samples if there

    group_db2, row_header2 = assignGroupColors(list(row_header)) ### row_header gets sorted in this function and will get permenantly screwed up if not mutated

    #if '.cdt' in filename: matrix.reverse(); row_header.reverse()
    for i in group_db2:
        if i not in group_db: group_db[i] = group_db2[i]

    return matrix, column_header, row_header, dataset_name, group_db

def importSIF(filename):
    fn = filepath(filename)
    edges=[]
    x=0
    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        parent,type,child = string.split(data,'\t')
        if 'AltAnalyze' in dataset_name:
            ### This is the order for proper directed interactions in the AltAnalyze-interaction viewer
            edges.append([parent,child,type])
        else:
            if '(' in parent: ### for TF-target annotations
                parent = string.split(parent,'(')[0]
            if ':' in child:
                child = string.split(child,':')[1]
    
            if 'TF' in dataset_name or 'UserSuppliedAssociations' in dataset_name or 'WGRV' in dataset_name:
                edges.append([parent,child,type]) ### Do this to indicate that the TF is regulating the target
            else:
                edges.append([child,parent,type])
    edges = unique.unique(edges)
    return edges

def assignGroupColors(t):
    """ Assign a unique color to each group. Optionally used for cluster display. """
    column_header=[]; group_number_db={}
    for i in t:
        repls = {'.2txt' : '', '.2bed' : '', '.2tab' : ''}
        i=reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        column_header.append(i)
        if ':' in i:
            group,j = string.split(i,':')[:2]
            group_number_db[group]=[]
            
    #import random
    k = 0
    group_db={}; color_db={}
    color_list = ['r', 'b', 'y', 'g', 'w', 'k', 'm']

    if len(group_number_db)>3:
        color_list = []
        cm = pylab.cm.get_cmap('gist_rainbow') #gist_ncar # binary
        for i in range(len(group_number_db)):
            color_list.append(cm(1.*i/len(group_number_db)))  # color will now be an RGBA tuple
    #color_list=[]
    #color_template = [1,1,1,0,0,0,0.5,0.5,0.5,0.25,0.25,0.25,0.75,0.75,0.75]
    t.sort() ### Ensure that all clusters have the same order of groups
    for i in t:
        repls = {'.2txt' : '', '.2bed' : '', '.2tab' : ''}
        i=reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        if ':' in i:
            group,j = string.split(i,':')[:2]
            try: color,ko = color_db[group]
            except Exception:
                try: color_db[group] = color_list[k],k
                except Exception:
                    ### If not listed in the standard color set add a new random color
                    rgb = tuple(scipy.rand(3)) ### random color
                    #rgb = tuple(random.sample(color_template,3)) ### custom alternative method
                    color_list.append(rgb)
                    color_db[group] = color_list[k], k
                color,ko = color_db[group]
                k+=1
            group_db[i] = group, color, ko
        #column_header.append(i)
    return group_db, column_header

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

def AppendOrWrite(export_path):
    export_path = filepath(export_path)
    status = verifyFile(export_path)
    if status == 'not found':
        export_data = export.ExportFile(export_path) ### Write this new file
    else:
        export_data = open(export_path,'a') ### Appends to existing file
    return export_path, export_data, status

def exportCustomGeneSet(geneSetName,species,allGenes):
    for gene in allGenes:break
    if 'ENS' not in gene:
        try:
            import gene_associations; import OBO_import
            gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        except Exception: symbol_to_gene={}
    
    if species != None:
        export_path, export_data, status = AppendOrWrite('AltDatabase/goelite/'+species+'/gene-mapp/Ensembl-StoredGeneSets.txt')

        stored_lines=[]
        for line in open(export_path,'rU').xreadlines(): stored_lines.append(line)
    
        if status == 'not found':
            export_data.write('GeneID\tEmpty\tGeneSetName\n')
        for gene in allGenes:
            if ' ' in gene:
                a,b=string.split(gene,' ')
                if 'ENS' in a: gene = a
                else: gene = b
            if 'ENS' not in gene and gene in symbol_to_gene:
                gene = symbol_to_gene[gene][0]
            line = gene+'\t\t'+geneSetName+'\n'
            if line not in stored_lines:
                export_data.write(line)
        export_data.close()
    else:
        print 'Could not store since no species name provided.'

    

def PrincipalComponentAnalysis(matrix, column_header, row_header, dataset_name, group_db, display=False, showLabels=True, algorithm='SVD',geneSetName=None, species=None, pcA=1,pcB=2):
    print "Performing Principal Component Analysis..."
    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank

    pcA-=1
    pcB-=1
    
    """ Based in part on code from:
    http://glowingpython.blogspot.com/2011/07/principal-component-analysis-with-numpy.html

    Performs performs principal components analysis 
    (PCA) on the n-by-p data matrix A
    Rows of A correspond to observations, columns to variables. 

    Returns :  
      coeff :
        is a p-by-p matrix, each column containing coefficients 
        for one principal component.
      score : 
        the principal component scores; that is, the representation 
        of A in the principal component space. Rows of SCORE 
        correspond to observations, columns to components.
    
      latent : 
        a vector containing the eigenvalues 
        of the covariance matrix of A.
    """
    # computing eigenvalues and eigenvectors of covariance matrix
    
    if algorithm == 'SVD': use_svd = True
    else: use_svd = False
    M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    Mdif = matrix/matrix.std()
    Mdif = Mdif.T
    u, s, vt = svd(Mdif, 0)
    fracs = s**2/np.sum(s**2)
    entropy = -sum(fracs*np.log(fracs))/np.log(np.min(vt.shape))
    
    label1 = 'PC%i (%2.1f%%)' %(pcA+1, fracs[0]*100)
    label2 = 'PC%i (%2.1f%%)' %(pcB+1, fracs[1]*100)

    ####  FROM LARSSON ########
    #100 most correlated Genes with PC1
    idx = numpy.argsort(u[0])
    idx2 = numpy.argsort(u[1])
    idx3 = numpy.argsort(u[2])
    idx4 = numpy.argsort(u[3])

    correlated_genes = map(lambda i: row_header[i],idx[:200])
    anticorrelated_genes = map(lambda i: row_header[i],idx[-200:])

    correlated_genes2 = map(lambda i: row_header[i],idx2[:200])
    anticorrelated_genes2 = map(lambda i: row_header[i],idx2[-200:])

    correlated_genes3 = map(lambda i: row_header[i],idx3[:200])
    anticorrelated_genes3 = map(lambda i: row_header[i],idx3[-200:])
    
    correlated_genes4 = map(lambda i: row_header[i],idx4[:200])
    anticorrelated_genes4 = map(lambda i: row_header[i],idx4[-200:])
    
    print 'exporting PCA driver genes to:',root_dir+'/PCA/correlated.txt'
    exportData = export.ExportFile(root_dir+'/PCA/correlated.txt')
    allGenes={}
    for gene in correlated_genes: exportData.write(gene+'\tcorrelated-PC1\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes: exportData.write(gene+'\tanticorrelated-PC1\n'); allGenes[gene]=[]
    for gene in correlated_genes2: exportData.write(gene+'\tcorrelated-PC2\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes2: exportData.write(gene+'\tanticorrelated-PC2\n'); allGenes[gene]=[]
    for gene in correlated_genes3: exportData.write(gene+'\tcorrelated-PC3\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes3: exportData.write(gene+'\tanticorrelated-PC3\n'); allGenes[gene]=[]
    for gene in correlated_genes4: exportData.write(gene+'\tcorrelated-PC4\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes4: exportData.write(gene+'\tanticorrelated-PC4\n'); allGenes[gene]=[]
    exportData.close()
    if geneSetName != None:
        if len(geneSetName)>0:
            exportCustomGeneSet(geneSetName,species,allGenes)
            print 'Exported geneset to "StoredGeneSets"'
        
    ###########################
    
    #if len(row_header)>20000:
    #print '....Using eigenvectors of the real symmetric square matrix for efficiency...'
    #[latent,coeff] = scipy.sparse.linalg.eigsh(cov(M))
    #scores=mlab.PCA(scores)
    
    if use_svd == False:
        [latent,coeff] = linalg.eig(cov(M))
        scores = dot(coeff.T,M) # projection of the data in the new space
    else:
        ### transform u into the same structure as the original scores from linalg.eig coeff
        scores = vt
    
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    pylab.xlabel(label1)
    pylab.ylabel(label2)
    pylab.title('Principal Component Analysis - '+dataset_name)
            
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)

    marker_size = 15
    if len(column_header)>20:
        marker_size = 12
    if len(column_header)>40:
        marker_size = 10
 
    group_names={}
    i=0
    for sample_name in column_header: #scores[0]
        ### Add the text labels for each
        try:
            ### Get group name and color information
            group_name,color,k = group_db[sample_name]
            if group_name not in group_names:
                label = group_name ### Only add once for each group
            else: label = None
            group_names[group_name] = color
        except Exception:
            color = 'r'; label=None
        try: ax.plot(scores[pcA][i],scores[1][i],color=color,marker='o',markersize=marker_size,label=label)
        except Exception: print i, len(scores[pcB]);sys.exit()
        if showLabels:
            try: sample_name = '   '+string.split(sample_name,':')[1]
            except Exception: pass
            ax.text(scores[pcA][i],scores[pcB][i],sample_name,fontsize=11)
        i+=1

    group_count = []
    for i in group_db:
        if group_db[i][0] not in group_count:
            group_count.append(group_db[i][0])
    
    #print len(group_count)
    Lfontsize = 8
    if len(group_count)>20:
        Lfontsize = 10
    if len(group_count)>30:
        Lfontsize = 8
    if len(group_count)>40:
        Lfontsize = 6
    if len(group_count)>50:
        Lfontsize = 5
    i=0
    
    box = ax.get_position()
    if len(group_count) > 4:
        # Shink current axis by 20%
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        try: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = Lfontsize) ### move the legend over to the right of the plot
        except Exception: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        ax.set_position([box.x0, box.y0, box.width, box.height])
        pylab.legend(loc="upper left", prop={'size': 10})
        
    filename = 'Clustering-%s-PCA.pdf' % dataset_name
    try: pylab.savefig(root_dir + filename)
    except Exception: None ### Rare error
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    try: pylab.savefig(root_dir + filename) #dpi=200
    except Exception: None ### Rare error
    graphic_link.append(['Principal Component Analysis',root_dir+filename])
    if display:
        print 'Exporting:',filename
        try:
            pylab.show()
        except Exception:
            pass### when run in headless mode
    fig.clf()

def ica(filename):
    X, column_header, row_header, dataset_name, group_db = importData(filename)
    X = map(numpy.array, zip(*X)) ### coverts these to tuples
    column_header, row_header = row_header, column_header

    ica = FastICA()
    S_ica_ = ica.fit(X).transform(X)  # Estimate the sources
    
    S_ica_ /= S_ica_.std(axis=0)
    
    pylab.plot()
    #plot_samples(S_ica_ / numpy.std(S_ica_))
    
    group_names={}
    i=0
    for sample_name in column_header: #scores[0]
        ### Add the text labels for each
        try:
            ### Get group name and color information
            group_name,color,k = group_db[sample_name]
            if group_name not in group_names:
                label = group_name ### Only add once for each group
            else: label = None
            group_names[group_name] = color
        except Exception:
            color = 'r'; label=None
        try: ax.plot(S_ica_[0][i],S_ica_[1][i],color=color,marker='o',markersize=marker_size,label=label)
        except Exception:
            print len(S_ica_)
            print i, len(S_ica_[0]);sys.exit()
        if showLabels:
            ax.text(scores[0][i],S_ica_[1][i],sample_name,fontsize=8)
        i+=1
        
    pylab.title('ICA recovered signals')
    
    pylab.show()

def plot_samples(S, axis_list=None):
    pylab.scatter(S[:, 0], S[:, 1], s=20, marker='o', linewidths=0, zorder=10,
                color='red', alpha=0.5)
    if axis_list is not None:
        colors = ['orange', 'red']
        for color, axis in zip(colors, axis_list):
            axis /= axis.std()
            x_axis, y_axis = axis
            # Trick to get legend to work
            pylab.plot(0.1 * x_axis, 0.1 * y_axis, linewidth=2, color=color)
            pylab.quiver(0, 0, x_axis, y_axis, zorder=11, width=2, scale=6,
                       color=color)

    pylab.xlabel('x')
    pylab.ylabel('y')
    

def PCA3D(matrix, column_header, row_header, dataset_name, group_db, display=False, showLabels=True, algorithm='SVD',geneSetName=None, species=None):
    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    start = time.time()
    M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    
    if algorithm == 'SVD': use_svd = True
    else: use_svd = False
    Mdif = matrix/matrix.std()
    Mdif = Mdif.T
    u, s, vt = svd(Mdif, 0)

    fracs = s**2/np.sum(s**2)
    entropy = -sum(fracs*numpy.oldnumeric.log(fracs))/np.log(np.min(vt.shape))
    
    label1 = 'PC%i (%2.1f%%)' %(0+1, fracs[0]*100)
    label2 = 'PC%i (%2.1f%%)' %(1+1, fracs[1]*100)
    label3 = 'PC%i (%2.1f%%)' %(2+1, fracs[2]*100)

    ####  FROM LARSSON ########
    #100 most correlated Genes with PC1
    idx = numpy.argsort(u[0])
    idx2 = numpy.argsort(u[1])
    idx3 = numpy.argsort(u[2])
    idx4 = numpy.argsort(u[3])
    correlated_genes = map(lambda i: row_header[i],idx[:200])
    anticorrelated_genes = map(lambda i: row_header[i],idx[-200:])

    correlated_genes2 = map(lambda i: row_header[i],idx2[:200])
    anticorrelated_genes2 = map(lambda i: row_header[i],idx2[-200:])

    correlated_genes3 = map(lambda i: row_header[i],idx3[:200])
    anticorrelated_genes3 = map(lambda i: row_header[i],idx3[-200:])
    
    correlated_genes4 = map(lambda i: row_header[i],idx4[:200])
    anticorrelated_genes4 = map(lambda i: row_header[i],idx4[-200:])
    
    print 'exporting PCA driver genes to:',root_dir+'/PCA/correlated.txt'
    exportData = export.ExportFile(root_dir+'/PCA/correlated.txt')
    
    allGenes={}
    for gene in correlated_genes: exportData.write(gene+'\tcorrelated-PC1\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes: exportData.write(gene+'\tanticorrelated-PC1\n'); allGenes[gene]=[]
    for gene in correlated_genes2: exportData.write(gene+'\tcorrelated-PC2\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes2: exportData.write(gene+'\tanticorrelated-PC2\n'); allGenes[gene]=[]
    for gene in correlated_genes3: exportData.write(gene+'\tcorrelated-PC3\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes3: exportData.write(gene+'\tanticorrelated-PC3\n'); allGenes[gene]=[]
    for gene in correlated_genes4: exportData.write(gene+'\tcorrelated-PC4\n'); allGenes[gene]=[]
    for gene in anticorrelated_genes4: exportData.write(gene+'\tanticorrelated-PC4\n'); allGenes[gene]=[]
    exportData.close()
    if geneSetName != None:
        if len(geneSetName)>0:
            exportCustomGeneSet(geneSetName,species,allGenes)
            print 'Exported geneset to "StoredGeneSets"'

    #numpy.Mdiff.toFile(root_dir+'/PCA/correlated.txt','\t')
    if use_svd == False:
        [latent,coeff] = linalg.eig(cov(M))
        scores = dot(coeff.T,M) # projection of the data in the new space
    else:
        ### transform u into the same structure as the original scores from linalg.eig coeff
        scores = vt

    end = time.time()
    print 'PCA completed in', end-start, 'seconds.'
    ### Hide the axis number labels
    #ax.w_xaxis.set_ticklabels([])
    #ax.w_yaxis.set_ticklabels([])
    #ax.w_zaxis.set_ticklabels([])

    #"""
    #ax.set_xticks([]) ### Hides ticks
    #ax.set_yticks([])
    #ax.set_zticks([])    
    
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    ax.set_zlabel(label3)
    #"""
    #pylab.title('Principal Component Analysis\n'+dataset_name)
    """
    pylab.figure()
    pylab.xlabel('Principal Component 1')
    pylab.ylabel('Principal Component 2')

    """
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)
    
    Lfontsize = 8
    group_count = []
    for i in group_db:
        if group_db[i][0] not in group_count:
            group_count.append(group_db[i][0])
    
    #print len(group_count)
    if len(group_count)>20:
        Lfontsize = 10
    if len(group_count)>30:
        Lfontsize = 8
    if len(group_count)>40:
        Lfontsize = 6
    if len(group_count)>50:
        Lfontsize = 5
    i=0
    group_names={}
    for x in scores[0]:
        ### Add the text labels for each
        sample_name = column_header[i]
        try:
            ### Get group name and color information
            group_name,color, k = group_db[sample_name]
            if group_name not in group_names:
                label = group_name ### Only add once for each group
            else: label = None
            group_names[group_name] = color, k
        except Exception:
            color = 'r'; label=None

        ax.plot([scores[0][i]],[scores[1][i]],[scores[2][i]],color=color,marker='o',markersize=7,label=label,markeredgewidth=0.1) #markeredgecolor=color
        if showLabels:
            #try: sample_name = '   '+string.split(sample_name,':')[1]
            #except Exception: pass
            ax.text(scores[0][i],scores[1][i],scores[2][i], '   '+sample_name,fontsize=9)
        i+=1

    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    #pylab.legend(loc="upper left", prop={'size': 10})
    try: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = Lfontsize) ### move the legend over to the right of the plot
    except Exception: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    filename = 'Clustering-%s-3D-PCA.pdf' % dataset_name
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename) #dpi=200
    graphic_link.append(['Principal Component Analysis',root_dir+filename])
    if display:
        print 'Exporting:',filename
        try: pylab.show()
        except Exception: None ### when run in headless mode
    fig.clf()
    
def getAxes(scores):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    try:
        x_range = max(scores[0])-min(scores[0])
        y_range = max(scores[1])-min(scores[1])
        x_axis_min = min(scores[0])-(x_range/1.5)
        x_axis_max = max(scores[0])+(x_range/1.5)
        y_axis_min = min(scores[1])-(y_range/5)
        y_axis_max = max(scores[1])+(y_range/5) 
    except KeyError:
        None
    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]

    
def Kmeans(features, column_header, row_header):
    #http://www.janeriksolem.net/2009/04/clustering-using-scipys-k-means.html
    #class1 = numpy.array(numpy.random.standard_normal((100,2))) + numpy.array([5,5]) 
    #class2 = 1.5 * numpy.array(numpy.random.standard_normal((100,2)))
    features = numpy.vstack((class1,class2))
    centroids,variance = scipy.cluster.vq.kmeans(features,2)
    code,distance = scipy.cluster.vq.vq(features,centroids)
    """
    This generates two normally distributed classes in two dimensions. To try and cluster the points, run k-means with k=2 like this.
    The variance is returned but we don't really need it since the SciPy implementation computes several runs (default is 20) and selects the one with smallest variance for us. Now you can check where each data point is assigned using the vector quantization function in the SciPy package.
    By checking the value of code we can see if there are any incorrect assignments. To visualize, we can plot the points and the final centroids.
    """
    import pylab
    pylab.plot([p[0] for p in class1],[p[1] for p in class1],'*')
    pylab.plot([p[0] for p in class2],[p[1] for p in class2],'r*') 
    pylab.plot([p[0] for p in centroids],[p[1] for p in centroids],'go') 
    pylab.show()

"""
def displaySimpleNetworkX():
    import networkx as nx
    print 'Graphing output with NetworkX'
    gr = nx.Graph(rotate=90,bgcolor='white') ### commands for neworkx and pygraphviz are the same or similiar

    edges = importSIF('Config/TissueFateMap.sif')

    ### Add nodes and edges
    for (node1,node2,type) in edges:
        gr.add_edge(node1,node2)
        draw_networkx_edges
   
    #gr['Myometrium']['color']='red'
    
    # Draw as PNG
    nx.draw_shell(gr) #wopi, gvcolor, wc, ccomps, tred, sccmap, fdp, circo, neato, acyclic, nop, gvpr, dot, sfdp. - fdp
    pylab.savefig('LineageNetwork.png')



def displaySimpleNetwork(sif_filename,fold_db,pathway_name):
    import pygraphviz as pgv
    #print 'Graphing output with PygraphViz'

    gr = pgv.AGraph(bgcolor='white',directed=True) ### Graph creation and setting of attributes - directed indicates arrows should be added
    #gr = pgv.AGraph(rotate='90',bgcolor='lightgray')

    ### Set graph attributes
    gr.node_attr['style']='filled'
    gr.graph_attr['label']='%s Network' % pathway_name

    edges = importSIF(sif_filename)
    if len(edges) > 700:
        print sif_filename, 'too large to visualize...'
    else:
        ### Add nodes and edges
        for (node1,node2,type) in edges:
            nodes = (node1,node2)
            gr.add_edge(nodes)
            child, parent = nodes
            edge = gr.get_edge(nodes[0],nodes[1])
            if 'TF' in pathway_name or 'WGRV' in pathway_name:
                node = child ### This is the regulating TF
            else:
                node = parent ### This is the pathway
            n=gr.get_node(node)
            ### http://www.graphviz.org/doc/info/attrs.html
            n.attr['penwidth'] = 4
            n.attr['fillcolor']= '#FFFF00' ### yellow
            n.attr['shape']='rectangle'
            #n.attr['weight']='yellow'
            #edge.attr['arrowhead'] = 'diamond' ### set the arrow type
        
        id_color_db = WikiPathways_webservice.getHexadecimalColorRanges(fold_db,'Genes')
        for gene_symbol in id_color_db:
            color_code = id_color_db[gene_symbol]
            try:
                n=gr.get_node(gene_symbol)
                n.attr['fillcolor']= '#'+string.upper(color_code) #'#FF0000'
                #n.attr['rotate']=90 
            except Exception: None
                   
        
        # Draw as PNG
        #gr.layout(prog='dot') #fdp (spring embedded), sfdp (OK layout), neato (compressed), circo (lots of empty space), dot (hierarchical - linear)
        gr.layout(prog='neato')
        output_filename = '%s.png' % sif_filename[:-4]
        #print output_filename
        gr.draw(output_filename)
"""
    
def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1  ### get just the parent directory
    return filename[:x]
    
def findFilename(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[x:]

def runHierarchicalClustering(matrix, row_header, column_header, dataset_name,
                              row_method, row_metric, column_method, column_metric,
                              color_gradient, display=False, contrast=None,
                              allowAxisCompression=True,Normalize=True):

    """ Running with cosine or other distance metrics can often produce negative Z scores
        during clustering, so adjustments to the clustering may be required.
        
    === Options Include ===
    row_method = 'average'
    column_method = 'single'
    row_metric = 'cosine'
    column_metric = 'euclidean'
    
    color_gradient = 'red_white_blue'
    color_gradient = 'red_black_sky'
    color_gradient = 'red_black_blue'
    color_gradient = 'red_black_green'
    color_gradient = 'yellow_black_blue'
    color_gradient == 'coolwarm'
    color_gradient = 'seismic'
    color_gradient = 'green_white_purple'
    """
    try:
        if allowLargeClusters: maxSize = 20000
        else: maxSize = 7000
    except Exception: maxSize = 7000
    
    run = False
    print 'max allowed cluster size:',maxSize
    if len(matrix)>0 and (len(matrix)<maxSize or row_method == None):
        #if len(matrix)>5000: row_metric = 'euclidean'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
            try:
                ### Default for display is False, when set to True, Pylab will render the image
                heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method,
                        row_metric, column_metric, color_gradient, dataset_name, display=display,
                        contrast=contrast,allowAxisCompression=allowAxisCompression,Normalize=Normalize)
                run = True
            except Exception:
                #print traceback.format_exc()
                try:
                    pylab.clf()
                    pylab.close() ### May result in TK associated errors later on
                    import gc
                    gc.collect()
                except Exception: None
                if len(matrix)<10000:
                    print 'Error using %s ... trying euclidean instead' % row_metric
                    row_metric = 'cosine'; row_method = 'average' ### cityblock
                else:
                    print 'Error with hierarchical clustering... only clustering arrays'
                    row_method = None ### Skip gene clustering
                try:
                    heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method,
                            row_metric, column_metric, color_gradient, dataset_name, display=display,
                            contrast=contrast,allowAxisCompression=allowAxisCompression,Normalize=Normalize)
                    run = True
                except Exception:
                    print traceback.format_exc()
                    print 'Unable to generate cluster due to dataset incompatibilty.'
    elif len(matrix)==0:
        print_out = 'SKIPPING HIERARCHICAL CLUSTERING!!! - Your dataset file has no associated rows.'
        print print_out
    else:
        print_out = 'SKIPPING HIERARCHICAL CLUSTERING!!! - Your dataset file is over the recommended size limit for clustering ('+str(maxSize)+' rows). Please cluster later using "Additional Analyses"'
        print print_out

    try:
        pylab.clf()
        pylab.close() ### May result in TK associated errors later on
        import gc
        gc.collect()
    except Exception: None
    return run
    
def debugTKBug():
    return None

def runHCexplicit(filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient,
                  extra_params, display=True, contrast=None, Normalize=False, JustShowTheseIDs=[],compressAxis=True):
    """ Explicit method for hiearchical clustering with defaults defined by the user (see below function) """
    
    global root_dir
    global inputFilename
    global originalFilename
    global graphic_link
    global allowLargeClusters
    global GroupDB
    global justShowTheseIDs
    global normalize
    global rho_cutoff
    global species
    global runGOElite
    global EliteGeneSets
    global storeGeneSetName
    EliteGeneSets=[]
    targetGene=[]
    filterByPathways=False
    
    runGOElite = False
    justShowTheseIDs = JustShowTheseIDs
    allowLargeClusters = True
    if compressAxis:
        allowAxisCompression = True
    else:
        allowAxisCompression = False
        
    graphic_link=graphics ### Store all locations of pngs
    inputFilename = filename ### Used when calling R
    filterIDs = False
    normalize = Normalize

    try:
        ### Specific additional optional parameters for filtering
        transpose = extra_params.Transpose()
        try:
            rho_cutoff = extra_params.RhoCutoff()
            print 'Setting correlation cutoff to a rho of',rho_cutoff
        except Exception: rho_cutoff = 0.5
        PathwayFilter = extra_params.PathwaySelect()
        GeneSet = extra_params.GeneSet()
        OntologyID = extra_params.OntologyID()
        Normalize = extra_params.Normalize()
        normalize = Normalize
        filterIDs = True
        species = extra_params.Species()
        platform = extra_params.Platform()
        vendor = extra_params.Vendor()
        newInput = findParentDir(inputFilename)+'/GeneSetClustering/'+findFilename(inputFilename)
        targetGene = extra_params.GeneSelection() ### Select a gene or ID to get the top correlating genes
        getGeneCorrelations = extra_params.GetGeneCorrelations() ### Select a gene or ID to get the top correlating genes
        filterByPathways = extra_params.FilterByPathways()
        PathwayFilter, filterByPathways = verifyPathwayName(PathwayFilter,GeneSet,OntologyID,filterByPathways)
        justShowTheseIDs_var = extra_params.JustShowTheseIDs()
        if len(justShowTheseIDs_var)>0:
            justShowTheseIDs = justShowTheseIDs_var
        elif len(targetGene)>0:
            targetGene = string.replace(targetGene,'\n',' ')
            targetGene = string.replace(targetGene,'\r',' ')
            justShowTheseIDs = string.split(targetGene,' ')
        try:
            EliteGeneSets = extra_params.ClusterGOElite()
            if EliteGeneSets != ['']: runGOElite = True
            #print EliteGeneSets
        except Exception: pass
        try:
            storeGeneSetName = extra_params.StoreGeneSetName()
        except Exception:
            storeGeneSetName = ''
    except Exception,e:
        #print traceback.format_exc();sys.exit()
        transpose = extra_params

    root_dir = findParentDir(filename)
    if 'ExpressionOutput/Clustering' in root_dir:
        root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
    elif 'ExpressionOutput' in root_dir:
        root_dir = string.replace(root_dir,'ExpressionOutput','DataPlots') ### Applies to clustering of LineageProfiler results
        root_dir = string.replace(root_dir,'/Clustering','') ### Applies to clustering of MarkerFinder results
    else:
        root_dir += '/DataPlots/'
        try: os.mkdir(root_dir) ### May need to create this directory
        except Exception: None
    if row_method == 'hopach': reverseOrder = False
    else: reverseOrder = True
    #"""
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,Normalize=Normalize,reverseOrder=reverseOrder)
    GroupDB = group_db
    inputFilename = string.replace(inputFilename,'.cdt','.txt')
    originalFilename = inputFilename
    #print len(matrix),;print len(column_header),;print len(row_header)
    if filterIDs:
        transpose_update = True ### Since you can filterByPathways and getGeneCorrelations, only transpose once
        if filterByPathways: ### Restrict analyses to only a single pathway/gene-set/ontology term
            if isinstance(PathwayFilter, tuple) or isinstance(PathwayFilter, list):
                 FileName = string.join(list(PathwayFilter),' ')
                 FileName = string.replace(FileName,':','-')
            else: FileName = PathwayFilter
            if len(FileName)>40:
                FileName = FileName[:40]
            try: inputFilename = string.replace(newInput,'.txt','_'+FileName+'.txt') ### update the pathway reference for HOPACH
            except Exception: inputFilename = string.replace(newInput,'.txt','_GeneSets.txt')
            vars = filterByPathway(matrix,row_header,column_header,species,platform,vendor,GeneSet,PathwayFilter,OntologyID,transpose)
            try: dataset_name += '-'+FileName
            except Exception: dataset_name += '-GeneSets'
            transpose_update = False
            if 'amplify' in targetGene:
                targetGene = string.join(vars[1],' ')+' amplify '+targetGene ### amplify the gene sets, but need the original matrix and headers (not the filtered)
            else: matrix,row_header,column_header = vars
            
        try:
            alt_targetGene = string.replace(targetGene,'amplify','')
            alt_targetGene = string.replace(alt_targetGene,'amplify','')
            alt_targetGene = string.replace(alt_targetGene,'driver','')
            alt_targetGene = string.replace(alt_targetGene,'top','')
            alt_targetGene = string.replace(alt_targetGene,'positive','')
            alt_targetGene = string.replace(alt_targetGene,' ','')  
        except Exception:
            alt_targetGene = ''
        if getGeneCorrelations and targetGene != 'driver' and targetGene !='excludeCellCycle' and targetGene !='top' and targetGene != 'monocle' and targetGene !='positive' and len(alt_targetGene)>0: ###Restrict analyses to only genes that correlate with the target gene of interest
            allowAxisCompression = False
            if transpose and transpose_update == False: transpose_update = False ### If filterByPathways selected
            elif transpose and transpose_update: transpose_update = True ### If filterByPathways not selected
            else: transpose_update = False ### If transpose == False
            if '\r' in targetGene or '\n' in targetGene:
                targetGene = string.replace(targetGene, '\r',' ')
                targetGene = string.replace(targetGene, '\n',' ')
            if len(targetGene)>15:
                inputFilename = string.replace(newInput,'.txt','_'+targetGene[:50]+'.txt') ### update the pathway reference for HOPACH
                dataset_name += '-'+targetGene[:50]
            else:
                inputFilename = string.replace(newInput,'.txt','_'+targetGene+'.txt') ### update the pathway reference for HOPACH
                dataset_name += '-'+targetGene
            inputFilename = root_dir+'/'+string.replace(findFilename(inputFilename),'|',' ')
            inputFilename = root_dir+'/'+string.replace(findFilename(inputFilename),':',' ') ### need to be careful of C://
            dataset_name = string.replace(dataset_name,'|',' ')
            dataset_name = string.replace(dataset_name,':',' ')
            try:
                matrix,row_header,column_header,row_method = getAllCorrelatedGenes(matrix,row_header,column_header,species,platform,vendor,targetGene,row_method,transpose_update)
            except Exception:
                print traceback.format_exc()
                print targetGene, 'not found in input expression file. Exiting. \n\n'
                badExit
            exportTargetGeneList(targetGene)      
    else:
        if transpose: ### Transpose the data matrix
            print 'Transposing the data matrix'
            matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
            column_header, row_header = row_header, column_header
    #print len(matrix),;print len(column_header),;print len(row_header)
    
    if len(column_header)>1000 or len(row_header)>1000:
        print 'Performing hierarchical clustering (please be patient)...'
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric,
                              column_method, column_metric, color_gradient, display=display,contrast=contrast,
                              allowAxisCompression=allowAxisCompression, Normalize=Normalize)
    #"""
    #graphic_link = [root_dir+'Clustering-exp.myeloid-steady-state-amplify positive Mki67 Clec4a2 Gria3 Ifitm6 Gfi1b -hierarchical_cosine_cosine.txt']

    if 'driver' in targetGene:
        import RNASeq
        input_file = graphic_link[-1][-1][:-4]+'.txt'
        if 'excludeCellCycle' in targetGene: excludeCellCycle = True
        else: excludeCellCycle = False
        print 'excludeCellCycle',excludeCellCycle
        targetGene = RNASeq.remoteGetDriverGenes(species,platform,input_file,excludeCellCycle=excludeCellCycle)
        extra_params.setGeneSelection(targetGene) ### force correlation to these
        extra_params.setGeneSet('None Selected') ### silence this
        graphic_link= runHCexplicit(filename, graphic_link, row_method, row_metric, column_method, column_metric, color_gradient,
                extra_params, display=display, contrast=contrast, Normalize=Normalize, JustShowTheseIDs=JustShowTheseIDs,compressAxis=compressAxis)
    return graphic_link

def exportTargetGeneList(targetGene):
    exclude=['positive','top','driver', 'amplify']
    exportFile = originalFilename[:-4]+'-targetGenes.txt'
    eo = export.ExportFile(root_dir+findFilename(exportFile))
    targetGenes = string.split(targetGene,' ')
    for gene in targetGenes:
        if gene not in exclude:
            eo.write(gene+'\n')
    eo.close()
    
def debugPylab():
    pylab.figure()
    pylab.close()
    pylab.figure()

def verifyPathwayName(PathwayFilter,GeneSet,OntologyID,filterByPathways):
    import gene_associations
    ### If the user supplied an Ontology ID rather than a Ontology term name, lookup the term name and return this as the PathwayFilter
    if len(OntologyID)>0:
        PathwayFilter = gene_associations.lookupOntologyID(GeneSet,OntologyID,type='ID')
        filterByPathways = True
    return PathwayFilter, filterByPathways

def filterByPathway(matrix,row_header,column_header,species,platform,vendor,GeneSet,PathwayFilter,OntologyID,transpose):
    ### Filter all the matrix and header entries for IDs in the selected pathway
    import gene_associations
    import OBO_import

    exportData = export.ExportFile(inputFilename)
    
    matrix2=[]; row_header2=[]
    if 'Ontology' in GeneSet: directory = 'nested'
    else: directory = 'gene-mapp'
    
    print "GeneSet(s) to analyze:",PathwayFilter
    if isinstance(PathwayFilter, tuple) or isinstance(PathwayFilter, list): ### see if it is one or more pathways
        associated_IDs={}
        for p in PathwayFilter:
            associated = gene_associations.simpleGenePathwayImport(species,GeneSet,p,OntologyID,directory)
            for i in associated:associated_IDs[i]=[]        
    else:
        associated_IDs = gene_associations.simpleGenePathwayImport(species,GeneSet,PathwayFilter,OntologyID,directory)
    gene_annotations = gene_associations.importGeneData(species,'Ensembl')
    vendor = string.replace(vendor,'other:','') ### For other IDs
    try: array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,associated_IDs)
    except Exception: array_to_ens={}
    
    if platform == "3'array":
        ### IDs thus won't be Ensembl - need to translate
        try:
            #ens_to_array = gene_associations.getGeneToUidNoExon(species,'Ensembl-'+vendor); print vendor, 'IDs imported...'
            array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,associated_IDs)
        except Exception:
            pass
            #print platform, vendor, 'not found!!! Exiting method'; badExit
        #array_to_ens = gene_associations.swapKeyValues(ens_to_array)
    try:
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    except Exception:
        pass
        
    i=0
    original_rows={} ### Don't add the same original ID twice if it associates with different Ensembl IDs
    for row_id in row_header:
        original_id = row_id; symbol = row_id
        if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename or 'AltConfirmed' in inputFilename or 'MarkerGenes' in inputFilename or 'blah' not in inputFilename:
            try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
            except Exception:
                try: symbol = gene_to_symbol[row_id][0]
                except Exception: None
            if len(symbol)==0: symbol = row_id
            if ':' in row_id:
                cluster,row_id = string.split(row_id,':')
                updated_row_id = cluster+':'+symbol
            else:
                updated_row_id = symbol
            original_id = updated_row_id
        if platform == "3'array":
            try:
                try: row_ids = array_to_ens[row_id]
                except Exception: row_ids = symbol_to_gene[symbol]
            except Exception:
                row_ids = [row_id]
        else:
            try:
                try: row_ids = array_to_ens[row_id]
                except Exception: row_ids = symbol_to_gene[symbol]
            except Exception:
                row_ids = [row_id]
        for row_id in row_ids:
            if row_id in associated_IDs:
                if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename:
                    if original_id != symbol:
                        row_id = original_id+' '+symbol
                    else: row_id = symbol
                else:
                    try: row_id = gene_annotations[row_id].Symbol()
                    except Exception: None ### If non-Ensembl data
                if original_id not in original_rows: ### Don't add the same ID twice if associated with mult. Ensembls
                    matrix2.append(matrix[i])
                    #row_header2.append(row_id)
                    row_header2.append(original_id)
                    original_rows[original_id]=None
        i+=1
        
    if transpose:
        matrix2 = map(numpy.array, zip(*matrix2)) ### coverts these to tuples
        column_header, row_header2 = row_header2, column_header

    exportData.write(string.join(['UID']+column_header,'\t')+'\n') ### title row export
    i=0
    for row_id in row_header2:
        exportData.write(string.join([row_id]+map(str,matrix2[i]),'\t')+'\n') ### export values
        i+=1
        
    print len(row_header2), 'filtered IDs'
    exportData.close()
    return matrix2,row_header2,column_header

def getAllCorrelatedGenes(matrix,row_header,column_header,species,platform,vendor,targetGene,row_method,transpose):
    ### Filter all the matrix and header entries for IDs in the selected targetGene
    resort_by_ID_name=False
    if resort_by_ID_name:
        index=0; new_row_header=[]; new_matrix=[]; temp_row_header = []
        for name in row_header: temp_row_header.append((name,index)); index+=1
        temp_row_header.sort()
        for (name,index) in temp_row_header:
            new_row_header.append(name)
            new_matrix.append(matrix[index])
        matrix = new_matrix
        row_header = new_row_header
        
    exportData = export.ExportFile(inputFilename)
    
    try:
        import gene_associations
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        #import OBO_import; symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    except Exception:
        pass
    
    if platform == "3'array":
        ### IDs thus won't be Ensembl - need to translate
        try:
            if ':' in vendor:
                vendor = string.split(vendor,':')[1]
            #ens_to_array = gene_associations.getGeneToUidNoExon(species,'Ensembl-'+vendor); print vendor, 'IDs imported...'
            array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,{})
        except Exception,e:
            array_to_ens={}

        for uid in array_to_ens:
            for gid in array_to_ens[uid]:
                if gid in gene_to_symbol:
                    symbol = gene_to_symbol[gid][0]
                    try: gene_to_symbol[uid].append(symbol)
                    except Exception: gene_to_symbol[uid] = [symbol]
        
    matrix2=[]
    row_header2=[]
    matrix_db={} ### Used to optionally sort according to the original order
    multipleGenes = False
    i=0
    ### If multiple genes entered, just display these
    if ' ' in targetGene or ',' in targetGene or '|' in targetGene or '\n' in targetGene or '\r' in targetGene:
        multipleGenes = True
        if ' ' in targetGene: delim = ' '
        if ',' in targetGene: delim = ','
        if '|' in targetGene and 'alt_junction' not in originalFilename: delim = '|'
        if '\n' in targetGene: delim = '\n'
        if '\r' in targetGene: delim = '\r'

        targetGenes = string.split(targetGene,delim)
        if row_method != None: targetGenes.sort()
        for row_id in row_header:
            original_rowid = row_id
            if ':' in row_id:
                a,b = string.split(row_id,':')[:2]
                if 'ENS' in a:
                    try:
                        row_id = a
                        symbol = gene_to_symbol[row_id][0]
                    except Exception: symbol =''
                elif 'ENS' not in b:
                    row_id = b
            try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
            except Exception:
                try: symbol = gene_to_symbol[row_id][0]
                except Exception:
                    row_id, symbol = row_id, row_id
            if 'ENS' not in original_rowid:
                if original_rowid != symbol:
                    symbol = original_rowid+' '+symbol

            for gene in targetGenes:
                if string.lower(gene) == string.lower(row_id) or string.lower(gene) == string.lower(symbol) or string.lower(original_rowid)==string.lower(gene):
                    matrix2.append(matrix[i]) ### Values for the row
                    row_header2.append(symbol)
                    matrix_db[symbol]=matrix[i]
            i+=1
        i=0
        #for gene in targetGenes:
        #    if gene not in matrix_db: print gene
    else:
        i=0
        original_rows={} ### Don't add the same original ID twice if it associates with different Ensembl IDs
        for row_id in row_header:
            original_id = row_id
            symbol = 'NA'
            if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename or 'blah' not in inputFilename:
                try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
                except Exception:
                    try: symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        row_id, symbol = row_id, row_id
                original_id = row_id
            
            if row_id == targetGene or symbol == targetGene:
                targetGeneValues = matrix[i] ### Values for the row
                break
            i+=1
        i=0
    if multipleGenes==False: limit = 50
    else: limit = 140 # lower limit is 132
    print 'limit:',limit
    
    if multipleGenes==False or 'amplify' in targetGene or 'correlated' in targetGene:
        row_header3=[] ### Convert to symbol if possible
        if multipleGenes==False:
            targetGeneValue_array = [targetGeneValues]
        else:
            targetGeneValue_array = matrix2
            if len(row_header2)>4:
                print 'Performing all pairwise corelations...',
                corr_matrix = numpyCorrelationMatrixGene(matrix,row_header,row_header2,gene_to_symbol)
                print 'complete'
            matrix2=[]; original_headers=row_header2; row_header2 = []
            matrix2_alt=[]; row_header2_alt=[]
        ### If one gene entered, display the most positive and negative correlated
        import markerFinder; k=0
        for targetGeneValues in targetGeneValue_array:
            try: targetGeneID = original_headers[k]
            except Exception: targetGeneID=''
            try:
                rho_results = list(corr_matrix[targetGeneID])
            except Exception:
                #print traceback.format_exc()
                rho_results = markerFinder.simpleScipyPearson(matrix,targetGeneValues)
            correlated_symbols={}
            #print targetGeneID, rho_results[:130][-1];sys.exit()
            for (rho,ind) in rho_results[:limit]: ### Get the top-50 correlated plus the gene of interest
                proceed = True
                if 'top' in targetGene:
                    if rho_results[4][0]<rho_cutoff: proceed = False
                if rho>rho_cutoff and proceed: #and rho_results[3][0]>rho_cutoff:# ensures only clustered genes considered
                    if len(row_header2)<100 or multipleGenes:
                        rh = row_header[ind]
                        #print rh, rho # Ly6c1, S100a8
                        if matrix[ind] not in matrix2:
                            if 'correlated' in targetGene:
                                if rho!=1:
                                    matrix2.append(matrix[ind])
                                    row_header2.append(rh)
                                    if targetGeneValues not in matrix2: ### gene ID systems can be different between source and query
                                        matrix2.append(targetGeneValues)
                                        row_header2.append(targetGeneID)
                                        try:correlated_symbols[gene_to_symbol[rh][0]]=ind
                                        except Exception: correlated_symbols[rh]=ind
                                        #print targetGeneValues, targetGene;sys.exit()
                            else:
                                matrix2.append(matrix[ind])
                                row_header2.append(rh)
                                try: correlated_symbols[gene_to_symbol[rh][0]]=ind
                                except Exception: correlated_symbols[rh]=ind
                            #if rho!=1: print gene_to_symbol[rh][0],'pos',targetGeneID
            #sys.exit()
            rho_results.reverse()
            for (rho,ind) in rho_results[:limit]: ### Get the top-50 anti-correlated plus the gene of interest
                if rho<-1*rho_cutoff and 'positive' not in targetGene:
                    if len(row_header2)<100 or multipleGenes:
                        rh = row_header[ind]
                        if matrix[ind] not in matrix2:
                            if 'correlated' in targetGene:
                                if rho!=1:
                                    matrix2.append(matrix[ind])
                                    row_header2.append(rh)
                                    if targetGeneValues not in matrix2:
                                        matrix2.append(targetGeneValues)
                                        row_header2.append(targetGeneID)
                                        try: correlated_symbols[gene_to_symbol[rh][0]]=ind
                                        except Exception: correlated_symbols[rh]=ind
                                        #print targetGeneValues, targetGene;sys.exit()
                            else:
                                matrix2.append(matrix[ind])
                                row_header2.append(rh)
                                try: correlated_symbols[gene_to_symbol[rh][0]]=ind
                                except Exception: correlated_symbols[rh]=ind
                            #if rho!=1: print gene_to_symbol[rh][0],'neg',targetGeneID
            try:
                ### print overlapping input genes that are correlated
                if len(correlated_symbols)>0:
                    potentially_redundant=[]
                    for i in targetGenes:
                        if i in correlated_symbols:
                            if i != targetGeneID: potentially_redundant.append((i,correlated_symbols[i]))
                    if len(potentially_redundant)>0:
                        ### These are intra-correlated genes based on the original filtered query
                        #print targetGeneID, potentially_redundant
                        for (rh,ind) in potentially_redundant:
                            matrix2_alt.append(matrix[ind])
                            row_header2_alt.append(rh)
                    rho_results.reverse()
                    #print targetGeneID, correlated_symbols, rho_results[:5]            
            except Exception:
                pass
            k+=1
            
        if 'IntraCorrelatedOnly' in targetGene:
            matrix2 = matrix2_alt
            row_header2 = row_header2_alt
            
        for r in row_header2:
            try:
                row_header3.append(gene_to_symbol[r][0])
            except Exception: row_header3.append(r) 
        row_header2 = row_header3
        #print len(row_header2),len(row_header3),len(matrix2);sys.exit()

        matrix2.reverse() ### Display from top-to-bottom rather than bottom-to-top (this is how the clusters are currently ordered in the heatmap)
        row_header2.reverse()
        if 'amplify' not in targetGene:
            row_method = None ### don't cluster the rows (row_method)
    else:
        ### reorder according to orignal
        matrix_temp=[]
        header_temp=[]
        for symbol in targetGenes:
            if symbol in matrix_db:
                matrix_temp.append(matrix_db[symbol]); header_temp.append(symbol)
        #print len(header_temp), len(matrix_db)
        if len(header_temp) >= len(matrix_db): ### Hence it worked and all IDs are the same type
            matrix2 = matrix_temp
            row_header2 = header_temp
        
    if transpose:
        matrix2 = map(numpy.array, zip(*matrix2)) ### coverts these to tuples
        column_header, row_header2 = row_header2, column_header

    exportData.write(string.join(['UID']+column_header,'\t')+'\n') ### title row export
    i=0
    for row_id in row_header2:
        if ':' in row_id:
            a,b = string.split(row_id,':')[:2]
            if 'ENS' in a:
                try: row_id=string.replace(row_id,a,gene_to_symbol[a][0])
                except Exception,e: pass
                row_header2[i] = row_id
        elif 'ENS' in row_id and ' ' in row_id:
            row_id = string.split(row_id, ' ')[1]
            row_header2[i] = row_id
        elif ' ' in row_id:
            try: a,b = string.split(row_id, ' ')
            except Exception: a = 1; b=2
            if a==b:
                row_id = a
        exportData.write(string.join([row_id]+map(str,matrix2[i]),'\t')+'\n') ### export values
        i+=1

    print len(row_header2), 'top-correlated IDs'
    exportData.close()

    return matrix2,row_header2,column_header,row_method

def numpyCorrelationMatrixGeneStore(x,rows,genes,gene_to_symbol):
    ### Decided not to use since it would require writing out the whole correlation matrix which is huge (1+GB) and time-intensive to import
    start = time.time()
    output_file = string.replace(originalFilename,'.txt','.corrmatrix')
    status = verifyFile(output_file)
    gene_correlations={}
    if status == 'found':
        try: symbol = gene_to_symbol[rows[i]][0]
        except Exception: symbol = '$'
        def splitInt(x):
            rho,ind = string.split(x,'|')
            return (float(rho),int(float(ind)))
        
        for line in open(output_file,'rU').xreadlines():         
            data = line.rstrip()
            t = string.split(data,'\t')
            scores = map(lambda x: splitInt(x), t[1:])
            gene_correlations[t[0]] = scores

    else:
        eo=export.ExportFile(output_file)
        D1 = numpy.ma.corrcoef(x)
        i=0
        for score_ls in D1:
            scores = []
            try: symbol = gene_to_symbol[rows[i]][0]
            except Exception: symbol = '$'
            if rows[i] in genes or symbol in genes:
                k=0
                for v in score_ls:
                    if str(v)!='nan':
                        scores.append((v,k))
                    k+=1    
                scores.sort()
                scores.reverse()
                if len(symbol)==1: symbol = rows[i]
                gene_correlations[symbol] = scores
                export_values = [symbol]
                for (v,k) in scores: ### re-import next time to save time
                    export_values.append(str(v)[:5]+'|'+str(k))
                eo.write(string.join(export_values,'\t')+'\n')
            i+=1
        eo.close()
    print len(gene_correlations)
    print time.time() - start, 'seconds';sys.exit()
    return gene_correlations

def numpyCorrelationMatrixGene(x,rows,genes,gene_to_symbol):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
        D1 = numpy.ma.corrcoef(x)
    i=0
    gene_correlations={}
    for score_ls in D1:
        scores = []
        try: symbol = gene_to_symbol[rows[i]][0]
        except Exception: symbol = '$'
        if rows[i] in genes or symbol in genes:
            k=0
            for v in score_ls:
                if str(v)!='nan':
                    scores.append((v,k))
                k+=1    
            scores.sort()
            scores.reverse()
            if len(symbol)==1: symbol = rows[i]
            gene_correlations[symbol] = scores 
        i+=1
    return gene_correlations

def runHCOnly(filename,graphics,Normalize=False):
    """ Simple method for hiearchical clustering with defaults defined by the function rather than the user (see above function) """
    
    global root_dir
    global graphic_link
    global inputFilename
    global GroupDB
    global allowLargeClusters
    allowLargeClusters = False
    
    graphic_link=graphics ### Store all locations of pngs
    inputFilename = filename ### Used when calling R
    
    root_dir = findParentDir(filename)
    if 'ExpressionOutput/Clustering' in root_dir:
        root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
    elif 'ExpressionOutput' in root_dir:
        root_dir = string.replace(root_dir,'ExpressionOutput','DataPlots') ### Applies to clustering of LineageProfiler results
    else: 
        root_dir += '/DataPlots/'
        try: os.mkdir(root_dir) ### May need to create this directory
        except Exception: None
        
    row_method = 'average'
    column_method = 'weighted'
    row_metric = 'cosine'
    column_metric = 'cosine'
    if 'Lineage' in filename or 'Elite' in filename:
        color_gradient = 'red_white_blue'
    else:
        color_gradient = 'yellow_black_blue'
        color_gradient = 'red_black_sky'
    
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,Normalize=Normalize)
    GroupDB = group_db
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name,
                row_method, row_metric, column_method, column_metric, color_gradient, display=False, Normalize=Normalize)
    return graphic_link

def runPCAonly(filename,graphics,transpose,showLabels=True,plotType='3D',display=True,algorithm='SVD',geneSetName=None, species=None):
    global root_dir
    global graphic_link
    graphic_link=graphics ### Store all locations of pngs
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
    root_dir = string.replace(root_dir,'ExpressionInput','DataPlots')
    if 'DataPlots' not in root_dir:
        root_dir += '/DataPlots/'
    try: os.mkdir(root_dir) ### May need to create this directory
    except Exception: None
        
    ### Transpose matrix and build PCA
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    
    if transpose == False: ### We normally transpose the data, so if True, we don't transpose (I know, it's confusing)
        matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
        column_header, row_header = row_header, column_header
    if len(column_header)>1000 or len(row_header)>1000:
        print 'Performing Principal Component Analysis (please be patient)...'
    #PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        if plotType == '3D':
            PCA3D(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName, species=species)
        else:
            PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName, species=species)
    return graphic_link

def outputClusters(filenames,graphics,Normalize=False):
    """ Peforms PCA and Hiearchical clustering on exported log-folds from AltAnalyze """
    
    global root_dir
    global graphic_link
    global inputFilename
    global GroupDB
    global allowLargeClusters
    global EliteGeneSets
    EliteGeneSets=[]
    global runGOElite
    runGOElite = False
    
    allowLargeClusters=False

    graphic_link=graphics ### Store all locations of pngs
    filename = filenames[0] ### This is the file to cluster with "significant" gene changes
    inputFilename = filename ### Used when calling R
    
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')

    ### Transpose matrix and build PCA
    original = importData(filename,Normalize=Normalize)
    matrix, column_header, row_header, dataset_name, group_db = original
    matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
    column_header, row_header = row_header, column_header
    if len(row_header)<700000 and len(column_header)<700000:
        PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db)
    else:
        print 'SKIPPING PCA!!! - Your dataset file is over the recommended size limit for clustering (>7000 rows). Please cluster later using "Additional Analyses".'

    row_method = 'average'
    column_method = 'average'
    row_metric = 'cosine'
    column_metric = 'cosine'
    color_gradient = 'red_white_blue'
    color_gradient = 'red_black_sky'
    
    ### Generate Significant Gene HeatMap
    matrix, column_header, row_header, dataset_name, group_db = original
    GroupDB = group_db
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, Normalize=Normalize)
    
    ### Generate Outlier and other Significant Gene HeatMap
    for filename in filenames[1:]:
        inputFilename = filename
        matrix, column_header, row_header, dataset_name, group_db = importData(filename,Normalize=Normalize)
        GroupDB = group_db
        try:
            runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, Normalize=Normalize)
        except Exception: print 'Could not cluster',inputFilename,', file not found'
    return graphic_link

def importEliteGeneAssociations(gene_filename):
    fn = filepath(gene_filename)
    x=0; fold_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': x=0
        elif x==0: x=1
        else:
            geneid=t[0];symbol=t[1]
            fold = 0
            try:
                if '|' in t[6]:
                    fold = float(string.split(t[6])[0]) ### Sometimes there are multiple folds for a gene (multiple probesets)
            except Exception:
                None
            try: fold=float(t[6])
            except Exception: None
            fold_db[symbol] = fold
    return fold_db

def importPathwayLevelFolds(filename):
    
    fn = filepath(filename)
    x=0
    folds_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if len(data)==0: x=0
        elif x==0:
            z_score_indexes = []; i=0
            z_headers = []
            for header in t:
                if 'z_score.' in header:
                    z_score_indexes.append(i)
                    header = string.split(header,'z_score.')[1] ### Get rid of z_score.
                    if 'AS.' in header:
                        header = string.split(header,'.p')[0] ### Remove statistics details
                        header = 'AS.'+string.join(string.split(header,'_')[2:],'_') ### species and array type notation
                    else:
                        header = string.join(string.split(header,'-')[:-2],'-')
                        if '-fold' in header:
                            header = string.join(string.split(header,'-')[:-1],'-')
                    z_headers.append(header)
                i+=1
            headers = string.join(['Gene-Set Name']+z_headers,'\t')+'\n'
            x=1
        else:
            term_name=t[1];geneset_type=t[2]
            zscores = map(lambda x: t[x], z_score_indexes)
            max_z = max(map(float, zscores)) ### If there are a lot of terms, only show the top 70
            line = string.join([term_name]+zscores,'\t')+'\n'
            try: zscore_db[geneset_type].append((max_z,line))
            except Exception: zscore_db[geneset_type] = [(max_z,line)]
    exported_files = []
    for geneset_type in zscore_db:
        ### Create an input file for hierarchical clustering in a child directory (Heatmaps)
        clusterinput_filename = findParentDir(filename)+'/Heatmaps/Clustering-Zscores-'+geneset_type+'.txt'
        exported_files.append(clusterinput_filename)
        export_text = export.ExportFile(clusterinput_filename)
        export_text.write(headers) ### Header is the same for each file
        zscore_db[geneset_type].sort()
        zscore_db[geneset_type].reverse()
        i=0 ### count the entries written
        for (max_z,line) in zscore_db[geneset_type]:
            if i<60:
                export_text.write(line) ### Write z-score values and row names
            i+=1
        export_text.close()
    return exported_files

def importOverlappingEliteScores(filename):
    fn = filepath(filename)
    x=0
    zscore_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if len(data)==0: x=0
        elif x==0:
            z_score_indexes = []; i=0
            z_headers = []
            for header in t:
                if 'z_score.' in header:
                    z_score_indexes.append(i)
                    header = string.split(header,'z_score.')[1] ### Get rid of z_score.
                    if 'AS.' in header:
                        header = string.split(header,'.p')[0] ### Remove statistics details
                        header = 'AS.'+string.join(string.split(header,'_')[2:],'_') ### species and array type notation
                    else:
                        header = string.join(string.split(header,'-')[:-2],'-')
                        if '-fold' in header:
                            header = string.join(string.split(header,'-')[:-1],'-')
                    z_headers.append(header)
                i+=1
            headers = string.join(['Gene-Set Name']+z_headers,'\t')+'\n'
            x=1
        else:
            term_name=t[1];geneset_type=t[2]
            zscores = map(lambda x: t[x], z_score_indexes)
            max_z = max(map(float, zscores)) ### If there are a lot of terms, only show the top 70
            line = string.join([term_name]+zscores,'\t')+'\n'
            try: zscore_db[geneset_type].append((max_z,line))
            except Exception: zscore_db[geneset_type] = [(max_z,line)]
    exported_files = []
    for geneset_type in zscore_db:
        ### Create an input file for hierarchical clustering in a child directory (Heatmaps)
        clusterinput_filename = findParentDir(filename)+'/Heatmaps/Clustering-Zscores-'+geneset_type+'.txt'
        exported_files.append(clusterinput_filename)
        export_text = export.ExportFile(clusterinput_filename)
        export_text.write(headers) ### Header is the same for each file
        zscore_db[geneset_type].sort()
        zscore_db[geneset_type].reverse()
        i=0 ### count the entries written
        for (max_z,line) in zscore_db[geneset_type]:
            if i<60:
                export_text.write(line) ### Write z-score values and row names
            i+=1
        export_text.close()
    return exported_files

def buildGraphFromSIF(mod,species,sif_filename,ora_input_dir):
    """ Imports a SIF and corresponding gene-association file to get fold changes for standardized gene-symbols """
    global SpeciesCode; SpeciesCode = species
    mod = 'Ensembl'
    if sif_filename == None:
        ### Used for testing only
        sif_filename = '/Users/nsalomonis/Desktop/dataAnalysis/collaborations/WholeGenomeRVista/Alex-Figure/GO-Elite_results/CompleteResults/ORA_pruned/up-2f_p05-WGRV.sif'
        ora_input_dir = '/Users/nsalomonis/Desktop/dataAnalysis/collaborations/WholeGenomeRVista/Alex-Figure/up-stringent/up-2f_p05.txt'
        #sif_filename = 'C:/Users/Nathan Salomonis/Desktop/Endothelial_Kidney/GO-Elite/GO-Elite_results/CompleteResults/ORA_pruned/GE.b_vs_a-fold2.0_rawp0.05-local.sif'
        #ora_input_dir = 'C:/Users/Nathan Salomonis/Desktop/Endothelial_Kidney/GO-Elite/input/GE.b_vs_a-fold2.0_rawp0.05.txt'
        
    gene_filename = string.replace(sif_filename,'.sif','_%s-gene-associations.txt') % mod
    gene_filename = string.replace(gene_filename,'ORA_pruned','ORA_pruned/gene_associations')
    pathway_name = string.split(sif_filename,'/')[-1][:-4]
    output_filename = None
    try: fold_db = importEliteGeneAssociations(gene_filename)
    except Exception: fold_db={}
    if ora_input_dir != None:
        ### This is an optional accessory function that adds fold changes from genes that are NOT in the GO-Elite pruned results (TFs regulating these genes)
        try: fold_db = importDataSimple(ora_input_dir,species,fold_db,mod)
        except Exception: None
    try:
        ### Alternative Approaches dependening on the availability of GraphViz
        #displaySimpleNetXGraph(sif_filename,fold_db,pathway_name)
        output_filename = iGraphSimple(sif_filename,fold_db,pathway_name)
    except Exception:
        print traceback.format_exc()
        try: displaySimpleNetwork(sif_filename,fold_db,pathway_name)
        except Exception: None ### GraphViz problem
    return output_filename

def iGraphSimple(sif_filename,fold_db,pathway_name):
    """ Build a network export using iGraph and Cairo """
    edges = importSIF(sif_filename)
    id_color_db = WikiPathways_webservice.getHexadecimalColorRanges(fold_db,'Genes')
    output_filename = iGraphDraw(edges,pathway_name,filePath=sif_filename,display=True,graph_layout='spring',colorDB=id_color_db)
    return output_filename

def iGraphDraw(edges, pathway_name, labels=None, graph_layout='shell', display=False,
               node_size=700, node_color='yellow', node_alpha=0.5, node_text_size=7,
               edge_color='black', edge_alpha=0.5, edge_thickness=2, edges_pos=.3,
               text_font='sans-serif',filePath='test',colorDB=None):
    ### Here node = vertex
    output_filename=None
    if len(edges) > 700 and 'AltAnalyze' not in pathway_name:
        print findFilename(filePath), 'too large to visualize...'
    else:
        arrow_scaler = 1 ### To scale the arrow
        if edges>40: arrow_scaler = .9
        
        vars = formatiGraphEdges(edges,pathway_name,colorDB,arrow_scaler)
        vertices,iGraph_edges,vertice_db,label_list,shape_list,vertex_size, color_list, vertex_label_colors, arrow_width, edge_colors = vars
        if vertices>0:
            import igraph
            gr = igraph.Graph(vertices, directed=True)
        
            canvas_scaler = 0.8 ### To scale the canvas size (bounding box)
            if vertices<15: canvas_scaler = 0.5
            elif vertices<25: canvas_scaler = .70
            elif vertices>35:
                canvas_scaler += len(iGraph_edges)/400.00
                
            filePath,canvas_scaler = correctedFilePath(filePath,canvas_scaler) ### adjust for GO-Elite
            #print vertices, len(iGraph_edges), pathway_name, canvas_scaler
            
            canvas_size = (600*canvas_scaler,600*canvas_scaler)
            gr.add_edges(iGraph_edges)
            gr.vs["label"] = label_list
            gr.vs["shape"] = shape_list
            gr.vs["size"] = vertex_size
            gr.vs["label_dist"] = [1.3]*vertices
            gr.vs["label_size"] = [12]*vertices
            gr.vs["color"]=color_list
            gr.vs["label_color"]=vertex_label_colors
            gr.es["color"] = edge_colors
            gr.es["arrow_size"]=arrow_width
    
            output_filename = '%s.pdf' % filePath[:-4]
            output_filename = output_filename.encode('ascii','ignore') ### removes the damned unicode u proceeding the filename
            layout = "kk"
            visual_style = {}
            #visual_style["layout"] = layout #The default is auto, which selects a layout algorithm automatically based on the size and connectedness of the graph
            visual_style["margin"] = 50 ### white-space around the network (see vertex size)
            visual_style["bbox"] = canvas_size
            igraph.plot(gr,output_filename, **visual_style)
            
            output_filename = '%s.png' % filePath[:-4]
            output_filename = output_filename.encode('ascii','ignore') ### removes the damned unicode u proceeding the filename
            if vertices <15: gr,visual_style = increasePlotSize(gr,visual_style)
            igraph.plot(gr,output_filename, **visual_style)
            #surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    return output_filename

def correctedFilePath(filePath,canvas_scaler):
    """ Move this file to it's own network directory for GO-Elite """
    if 'ORA_pruned' in filePath:
        filePath = string.replace(filePath,'CompleteResults/ORA_pruned','networks')
        try: os.mkdir(findParentDir(filePath))
        except Exception: pass
        canvas_scaler = canvas_scaler*1.3 ### These graphs tend to be more dense and difficult to read
    return filePath,canvas_scaler
    
def increasePlotSize(gr,visual_style):
    ### To display the plot better, need to manually increase the size of everything
    factor = 2
    object_list = ["size","label_size"]
    for i in object_list:
        new=[]
        for k in gr.vs[i]:
            new.append(k*factor)
        gr.vs[i] = new

    new=[]
    for i in gr.es["arrow_size"]:
        new.append(i*factor)
        
    new=[]
    for i in visual_style["bbox"]:
        new.append(i*factor)
    visual_style["bbox"] = new
    visual_style["margin"]=visual_style["margin"]*factor
    return gr,visual_style
    
def getHMDBDataSimple():
    ### Determine which IDs are metabolites
    program_type,database_dir = unique.whatProgramIsThis()
    filename = database_dir+'/'+SpeciesCode+'/gene/HMDB.txt'
    symbol_hmdb_db={}
    x=0
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            hmdb_id = t[0]; symbol = t[1]; ProteinNames = t[-1]
            symbol_hmdb_db[symbol]=hmdb_id
    return symbol_hmdb_db

def formatiGraphEdges(edges,pathway_name,colorDB,arrow_scaler):
    ### iGraph appears to require defined vertice number and edges as numbers corresponding to these vertices
    edge_db={}
    edges2=[]
    vertice_db={}
    shape_list=[] ### node shape in order
    label_list=[] ### Names of each vertix aka node
    vertex_size=[]
    color_list=[]
    vertex_label_colors=[]
    arrow_width=[] ### Indicates the presence or absence of an arrow
    edge_colors=[]
    k=0
    
    try: symbol_hmdb_db = getHMDBDataSimple()
    except Exception: symbol_hmdb_db={}

    for (node1,node2,type) in edges:
        edge_color = 'grey'
        ### Assign nodes to a numeric vertix ID
        if 'TF' in pathway_name or 'WGRV' in pathway_name:
            pathway = node1 ### This is the regulating TF
        else:
            pathway = node2 ### This is the pathway
        
        if 'drugInteraction' == type: edge_color = "purple"
        elif 'TBar' == type: edge_color = 'blue'
        elif 'microRNAInteraction' == type: edge_color = '#53A26D'
        elif 'transcription' in type: edge_color = '#FF7D7D'
        if 'AltAnalyze' in pathway_name: default_node_color = 'grey'
        else: default_node_color = "yellow"
        if node1 in vertice_db: v1=vertice_db[node1]
        else: #### Left hand node
            ### Only time the vertex is added to the below attribute lists
            v1=k; label_list.append(node1)
            rs = 1 ### relative size
        
            if 'TF' in pathway_name or 'WGRV' in pathway_name and 'AltAnalyze' not in pathway_name:
                shape_list.append('rectangle')
                vertex_size.append(15)
                vertex_label_colors.append('blue')    
            else:
                if 'drugInteraction' == type:
                    rs = 0.75
                    shape_list.append('rectangle')
                    vertex_label_colors.append('purple')
                    default_node_color = "purple"
                elif 'Metabolic' == type and node1 in symbol_hmdb_db:
                    shape_list.append('triangle-up')
                    vertex_label_colors.append('blue') #dark green
                    default_node_color = 'grey' #'#008000'  
                elif 'microRNAInteraction' == type:
                    rs = 0.75
                    shape_list.append('triangle-up')
                    vertex_label_colors.append('#008000') #dark green
                    default_node_color = 'grey' #'#008000'
                else:
                    shape_list.append('circle')
                    vertex_label_colors.append('black')
                vertex_size.append(10*rs)
            vertice_db[node1]=v1; k+=1
            try:
                color = '#'+string.upper(colorDB[node1])
                color_list.append(color) ### Hex color
            except Exception:
                color_list.append(default_node_color)
        if node2 in vertice_db: v2=vertice_db[node2]
        else: #### Right hand node
            ### Only time the vertex is added to the below attribute lists
            v2=k; label_list.append(node2)
            if 'TF' in pathway_name or 'WGRV' in pathway_name:
                shape_list.append('circle')
                vertex_size.append(10)
                vertex_label_colors.append('black')
                default_node_color = "grey"
            elif 'AltAnalyze' not in pathway_name:
                shape_list.append('rectangle')
                vertex_size.append(15)
                vertex_label_colors.append('blue')
                default_node_color = "grey"
            elif 'Metabolic' == type and node2 in symbol_hmdb_db:
                shape_list.append('triangle-up')
                vertex_label_colors.append('blue') #dark green
                default_node_color = 'grey' #'#008000'   
            else:
                shape_list.append('circle')
                vertex_size.append(10)
                vertex_label_colors.append('black')
                default_node_color = "grey"
            vertice_db[node2]=v2; k+=1
            try:
                color = '#'+string.upper(colorDB[node2])
                color_list.append(color) ### Hex color
            except Exception: color_list.append(default_node_color)
        edges2.append((v1,v2))

        if type == 'physical': arrow_width.append(0)
        else: arrow_width.append(arrow_scaler)
        try: edge_db[v1].append(v2)
        except Exception: edge_db[v1]=[v2]
        try: edge_db[v2].append(v1)
        except Exception: edge_db[v2]=[v1]
        edge_colors.append(edge_color)
    vertices = len(edge_db) ### This is the number of nodes
    edge_db = eliminate_redundant_dict_values(edge_db)

    vertice_db2={} ### Invert
    for node in vertice_db:
        vertice_db2[vertice_db[node]] = node
    #print len(edges2), len(edge_colors)
    print vertices, 'and', len(edges2),'edges in the iGraph network.'
    return vertices,edges2,vertice_db2, label_list, shape_list, vertex_size, color_list, vertex_label_colors, arrow_width, edge_colors

def eliminate_redundant_dict_values(database):
    db1={}
    for key in database: list = unique.unique(database[key]); list.sort(); db1[key] = list
    return db1

def importDataSimple(filename,species,fold_db,mod):
    """ Imports an input ID file and converts those IDs to gene symbols for analysis with folds """
    import GO_Elite
    import OBO_import
    import gene_associations
    fn = filepath(filename)
    x=0
    metabolite_codes = ['Ck','Ca','Ce','Ch','Cp']
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': x=0
        
        elif x==0: x=1
        else:
            if x == 1:
                system_code = t[1]
                if system_code in metabolite_codes:
                    mod = 'HMDB'
                system_codes,source_types,mod_types = GO_Elite.getSourceData()
                try: source_data = system_codes[system_code]
                except Exception:
                    source_data = None
                    if 'ENS' in t[0]: source_data = system_codes['En']
                    else: ### Assume the file is composed of gene symbols
                        source_data = system_codes['Sy']
                if source_data == mod:
                    source_is_mod = True
                elif source_data==None:
                    None ### Skip this
                else:
                    source_is_mod = False
                    mod_source = mod+'-'+source_data+'.txt'
                    gene_to_source_id = gene_associations.getGeneToUid(species,('hide',mod_source))
                    source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)

                try: gene_to_symbol = gene_associations.getGeneToUid(species,('hide',mod+'-Symbol'))
                except Exception: gene_to_symbol={}
                try: met_to_symbol = gene_associations.importGeneData(species,'HMDB',simpleImport=True)
                except Exception: met_to_symbol={}
                for i in met_to_symbol: gene_to_symbol[i] = met_to_symbol[i] ### Add metabolite names

            x+=1
            if source_is_mod == True:
                if t[0] in gene_to_symbol:
                    symbol = gene_to_symbol[t[0]][0]
                    try: fold_db[symbol] = float(t[2])
                    except Exception: fold_db[symbol] = 0
                else:
                    fold_db[t[0]] = 0 ### If not found (wrong ID with the wrong system) still try to color the ID in the network as yellow
            elif t[0] in source_to_gene:
                mod_ids = source_to_gene[t[0]]
                try: mod_ids+=source_to_gene[t[2]] ###If the file is a SIF
                except Exception:
                    try: mod_ids+=source_to_gene[t[1]] ###If the file is a SIF
                    except Exception: None
                for mod_id in mod_ids:
                    if mod_id in gene_to_symbol:
                        symbol = gene_to_symbol[mod_id][0]
                        try: fold_db[symbol] = float(t[2]) ### If multiple Ensembl IDs in dataset, only record the last associated fold change
                        except Exception: fold_db[symbol] = 0
            else: fold_db[t[0]] = 0
    return fold_db

def clusterPathwayZscores(filename):
    """ Imports a overlapping-results file and exports an input file for hierarchical clustering and clusters """
    ### This method is not fully written or in use yet - not sure if needed
    if filename == None:
        ### Only used for testing
        filename = '/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/GO-Elite/TF-enrichment2/GO-Elite_results/overlapping-results_z-score_elite.txt'
    exported_files = importOverlappingEliteScores(filename)
    graphic_links=[]
    for file in exported_files:
        try: graphic_links = runHCOnly(file,graphic_links)
        except Exception,e:
            #print e
            print 'Unable to generate cluster due to dataset incompatibilty.'
    print 'Clustering of overlapping-results_z-score complete (see "GO-Elite_results/Heatmaps" directory)'
    
def clusterPathwayMeanFolds():
    """ Imports the pruned-results file and exports an input file for hierarchical clustering and clusters """
    
    filename = '/Users/nsalomonis/Desktop/User Diagnostics/Mm_spinal_cord_injury/GO-Elite/GO-Elite_results/pruned-results_z-score_elite.txt'
    exported_files = importPathwayLevelFolds(filename)

def VennDiagram():
    f = pylab.figure()
    ax = f.gca()
    rad = 1.4
    c1 = Circle((-1,0),rad, alpha=.2, fc ='red',label='red')
    c2 = Circle((1,0),rad, alpha=.2, fc ='blue',label='blue')
    c3 = Circle((0,1),rad, alpha=.2, fc ='green',label='g')
    #pylab.plot(c1,color='green',marker='o',markersize=7,label='blue')
    #ax.add_patch(c1)
    ax.add_patch(c2)
    ax.add_patch(c3)
    ax.set_xlim(-3,3)
    ax.set_ylim(-3,3)
    pylab.show()
    
def plotHistogram(filename):
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    transpose=True
    if transpose: ### Transpose the data matrix
        print 'Transposing the data matrix'
        matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
        column_header, row_header = row_header, column_header
            
    pylab.figure()
    for i in matrix:
        pylab.hist(i, 200, normed=0, histtype='step', cumulative=-1)
    #pylab.hist(matrix, 50, cumulative=-1)
    pylab.show()

def multipleSubPlots(filename,uids,SubPlotType='column'):
    uids = [uids[-1]]+uids[:-1]
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,geneFilter=uids)
    fig = pylab.figure()
    def ReplaceZeros(val,min_val):
        if val == 0: 
            return min_val
        else: return val
        
    ### Order the graphs based on the original gene order
    new_row_header=[]
    matrix2 = []
    for uid in uids:
        if uid in row_header:
            ind = row_header.index(uid)
            new_row_header.append(uid)
            try: update_exp_vals = map(lambda x: ReplaceZeros(x,0.0001),matrix[ind])
            except Exception: print uid, len(matrix[ind]);sys.exit()
            matrix2.append(update_exp_vals)
    matrix = numpy.array(matrix2)
    row_header = new_row_header
    color_list = ['r', 'b', 'y', 'g', 'w', 'k', 'm']
    fontsize=10
    if len(row_header)>3:
        color_list = []
        cm = pylab.cm.get_cmap('gist_rainbow') #gist_ncar
        for i in range(len(row_header)):
            color_list.append(cm(1.*i/len(row_header)))  # color will now be an RGBA tuple
            
    for i in range(len(matrix)):
        temp = 510 + i
        try: ax = pylab.subplot(temp)
        except Exception: break
        OY = matrix[i]
        pylab.xlim(0,len(OY))
        pylab.subplots_adjust(right=0.85)
        ind = np.arange(len(OY))
        if SubPlotType=='column':
            pylab.bar(ind, OY,edgecolor='black',linewidth=0,color=color_list[i])
            width = .35
            print i ,row_header[i]
        if SubPlotType=='plot':
            pylab.plot(x,y)

        ax.text(matrix.shape[1]-0.5, i, '  '+row_header[i],fontsize=16)
        
        fig.autofmt_xdate()
        pylab.subplots_adjust(hspace = .001)
        temp = tic.MaxNLocator(3)
        ax.yaxis.set_major_locator(temp)
        ax.set_xticks([])
        #ax.title.set_visible(False)
        #pylab.xticks(ind + width / 2, column_header)
        #ax.set_xticklabels(column_header)
        #ax.xaxis.set_ticks([-1]+range(len(OY)+1))
        #xtickNames = pylab.setp(pylab.gca(), xticklabels=['']+column_header)
        #pylab.setp(xtickNames, rotation=90, fontsize=10)
        
    pylab.show()

def simpleTranspose(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,' ')
        matrix.append(t)

    matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
    filename = filename[:-4]+'-transposed.txt'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i,'\t')+'\n')
    ea.close()
 
def CorrdinateToBed(filename):
    fn = filepath(filename)
    matrix = []
    translation={}
    multiExon={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        t = string.split(data,'\t')
        if '.gtf' in filename:
            if 'chr' not in t[0]: chr = 'chr'+t[0]
            else: chr = t[0]
            start = t[3]; end = t[4]; strand = t[6]; annotation = t[8]
            annotation = string.replace(annotation,'gene_id','')
            annotation = string.replace(annotation,'transcript_id','')
            annotation = string.replace(annotation,'gene_name','')
            geneIDs = string.split(annotation,';')
            geneID = geneIDs[0]; symbol = geneIDs[3]
        else:
            chr = t[4]; strand = t[5]; start = t[6]; end = t[7]
        #if 'ENS' not in annotation:
        t = [chr,start,end,geneID,'0',strand]
        #matrix.append(t)
        translation[geneID] = symbol
        try: multiExon[geneID]+=1
        except Exception: multiExon[geneID]=1
    filename = filename[:-4]+'-new.bed'
    ea = export.ExportFile(filename)
    for i in translation:
        #ea.write(string.join(i,'\t')+'\n')
        ea.write(i+'\t'+translation[i]+'\t'+str(multiExon[i])+'\n')
    ea.close()
    
def SimpleCorrdinateToBed(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        t = string.split(data,'\t')
        if '.bed' in filename:
            print t;sys.exit()
        chr = t[4]; strand = t[5]; start = t[6]; end = t[7]
        if 'ENS' in t[0]:
            t = [chr,start,end,t[0],'0',strand]
            matrix.append(t)

    filename = filename[:-4]+'-new.bed'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i,'\t')+'\n')
    ea.close()
    
def simpleIntegrityCheck(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        t = string.split(data,'\t')
        matrix.append(t)

    filename = filename[:-4]+'-new.bed'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i,'\t')+'\n')
    ea.close()

def simpleFilter(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        t = string.split(data,'\t')
        if t[1] != '-':
            matrix.append(t)

    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i,'\t')+'\n')
    ea.close()
    
def test(filename):
    fn = filepath(filename)
    for line in open(fn,'rU').xreadlines():
        print [line];sys.exit()
        data = cleanUpLine(line)
        data = string.replace(data,' ','')
        t = string.split(data,'\t')
        if t[1] != '-':
            matrix.append(t)

    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i,'\t')+'\n')
    ea.close()
    
    
if __name__ == '__main__':
    import UI
    folder = 'clustering.py'
    test(folder);sys.exit()
    files = UI.read_directory(folder)
    for file in files:
        SimpleCorrdinateToBed(folder+'/'+file)
    filename = '/Users/saljh8/Desktop/bed/RREs0.5_exons_unique.txt'
    #simpleIntegrityCheck(filename);sys.exit()
    gene_list = ['Cx3cr1','Csf1r','Flt3','Dnmt3a','Gfi1']
    filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/genes.tpm_tracking-output.txt'
    multipleSubPlots(filename,gene_list,SubPlotType='column');sys.exit()

    plotHistogram(filename);sys.exit()
    filename = '/Users/saljh8/Desktop/Grimes/Expression_final_files/ExpressionInput/amplify-wt/DataPlots/Clustering-exp.myeloid-steady-state-PCA-all_wt_myeloid_SingleCell-Klhl7 Dusp7 Slc25a33 H6pd Bcorl1 Sdpr Ypel3 251000-hierarchical_cosine_cosine.cdt'
    openTreeView(filename);sys.exit()
    pdf1 = "/Users/saljh8/Desktop/Grimes/1.pdf"
    pdf2 = "/Users/saljh8/Desktop/Grimes/2.pdf"
    outPdf = "/Users/saljh8/Desktop/Grimes/3.pdf"
    merge_horizontal(outPdf, pdf1, pdf2);sys.exit()
    mergePDFs(pdf1,pdf2,outPdf);sys.exit()
    filename = '/Volumes/SEQ-DATA/CardiacRNASeq/BedFiles/ExpressionOutput/Clustering/SampleLogFolds-CardiacRNASeq.txt'
    ica(filename);sys.exit()
    features = 5
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    Kmeans(features, column_header, row_header); sys.exit()
    #graphViz();sys.exit()
    filename = '/Users/saljh8/Desktop/delete.txt'
    runHCOnly(filename,[]); sys.exit()
    filenames = [filename]
    outputClusters(filenames,[]); sys.exit()
    #runPCAonly(filename,[],False);sys.exit()
    #VennDiagram(); sys.exit()
    #buildGraphFromSIF('Ensembl','Mm',None,None); sys.exit()
    #clusterPathwayZscores(None); sys.exit()
    pruned_folder = '/Users/nsalomonis/Desktop/CBD/LogTransformed/GO-Elite/GO-Elite_results/CompleteResults/ORA_pruned/'
    input_ora_folder = '/Users/nsalomonis/Desktop/CBD/LogTransformed/GO-Elite/input/'

    files = UI.read_directory(pruned_folder)
    for file in files:
        if '.sif' in file:
            input_file = string.join(string.split(file,'-')[:-1],'-')+'.txt'
            sif_file = pruned_folder+file
            input_file = input_ora_folder+input_file
            buildGraphFromSIF('Ensembl','Hs',sif_file,input_file)
    sys.exit()
    filenames = [filename]
    outputClusters(filenames,[])