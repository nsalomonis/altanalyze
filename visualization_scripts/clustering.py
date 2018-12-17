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

import sys,string,os,copy
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>0 and '--' in command_args: commandLine=True
else: commandLine=False

display_label_names = True

import traceback
try:
    import math
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        if commandLine and 'linux' in sys.platform:
            ### TkAgg doesn't work when AltAnalyze is run remotely (ssh or sh script)
            try: matplotlib.use('Agg')
            except Exception: pass
            try: matplotlib.rcParams['backend'] = 'Agg'
            except Exception: pass
        else:
            try: matplotlib.use('TkAgg')
            except Exception: pass
            try: matplotlib.rcParams['backend'] = 'TkAgg'
            except Exception: pass
        try:
            import matplotlib.pyplot as pylab
            import matplotlib.colors as mc
            import matplotlib.mlab as mlab
            import matplotlib.ticker as tic
            from matplotlib.patches import Circle
            from mpl_toolkits.mplot3d import Axes3D
            try: from matplotlib.cbook import _string_to_bool
            except: pass
            matplotlib.rcParams['axes.linewidth'] = 0.5
            matplotlib.rcParams['pdf.fonttype'] = 42
            matplotlib.rcParams['font.family'] = 'sans-serif'
            matplotlib.rcParams['font.sans-serif'] = 'Arial'
            matplotlib.rcParams['figure.facecolor'] = 'white' ### Added in 2.1.2
        except Exception:
            print traceback.format_exc()
            print 'Matplotlib support not enabled'
        import scipy
        try: from scipy.sparse.csgraph import _validation
        except Exception: pass
        from scipy.linalg import svd
        import scipy.cluster.hierarchy as sch
        import scipy.spatial.distance as dist
        try: import numpy; np = numpy
        except Exception:
            print 'Numpy import error...'
            print traceback.format_exc()
        try: import umap
        except: pass
        try:
            from cairo import ImageSurface
        except: pass
        try:
            import igraph.vendor.texttable
        except ImportError: pass
        try:
            from sklearn.decomposition import PCA, FastICA
        except Exception: pass
        
        from sklearn.neighbors import quad_tree
        from sklearn.neighbors import *; from sklearn.manifold.t_sne import *
        from sklearn.tree import *; from sklearn.tree import _utils
        from sklearn.manifold.t_sne import _utils
        from sklearn.manifold import TSNE
        try:
            import numba
            import llvmlite; from llvmlite import binding; from llvmlite.binding import *
            from llvmlite.binding import ffi; from llvmlite.binding import dylib
        except:
            pass
        #pylab.ion() # closes Tk window after show - could be nice to include
except Exception:
    print traceback.format_exc()
    pass

import time
import unique
from stats_scripts import statistics
import os
import export
import webbrowser
import warnings
import UI
use_default_colors = False

try:
    warnings.simplefilter("ignore", numpy.ComplexWarning)
    warnings.simplefilter("ignore", DeprecationWarning) ### Annoying depreciation warnings (occurs in sch somewhere)
    #This shouldn't be needed in python 2.7 which suppresses DeprecationWarning - Larsson
except Exception: None

from visualization_scripts import WikiPathways_webservice

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
            dataset_name, display=False, contrast=None, allowAxisCompression=True,Normalize=True,
            PriorColumnClusters=None, PriorRowClusters=None):
    print "Performing hieararchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)
    show_color_bars = True ### Currently, the color bars don't exactly reflect the dendrogram colors
    try: ExportCorreleationMatrix = exportCorreleationMatrix
    except Exception: ExportCorreleationMatrix = False
    
    try: os.mkdir(root_dir) ### May need to create this directory
    except Exception: None
    if display == False:
        pylab.figure() ### Add this to avoid a Tkinter bug after running MarkerFinder (not sure why it is needed) - creates a second empty window when display == True
        
    if row_method == 'hopach' or column_method == 'hopach':
        ### Test R and hopach
        """
        try:
            import R_test
        except Exception,e:
            #print traceback.format_exc()
            print 'Failed to install hopach or R not installed (install R before using hopach)'
            row_method = 'average'; column_method = 'average'
        if len(column_header)==2: column_method = 'average'
        if len(row_header)==2: row_method = 'average'
        """
        pass
    
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

    if row_method == 'hopach' or column_method == 'hopach':
        try:
    
            """ HOPACH is a clustering method implemented in R that builds a hierarchical tree of clusters by recursively
            partitioning a data set, while ordering and possibly collapsing clusters at each level:
            http://www.bioconductor.org/packages/release/bioc/html/hopach.html
            """
            
            import R_interface
            #reload(R_interface)
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
       
    skipClustering = False
    try:
        if len(PriorColumnClusters)>0 and len(PriorRowClusters)>0 and row_method==None and column_method == None:
            print 'Prior generated clusters being used rather re-clustering'
            """
            try:
                if len(targetGeneIDs)>0:
                    PriorColumnClusters=[] ### If orderded genes input, we want to retain this order rather than change
            except Exception: pass
            """
            if len(PriorColumnClusters)>0: ### this corresponds to the above line
                Z1={}; Z2={}
                Z1['level'] = PriorRowClusters; Z1['level'].reverse()
                Z2['level'] = PriorColumnClusters; #Z2['level'].reverse()
                Z1['leaves'] = range(0,len(row_header)); #Z1['leaves'].reverse()
                Z2['leaves'] = range(0,len(column_header)); #Z2['leaves'].reverse()
                skipClustering = True
                ### When clusters are imported, you need something other than None, otherwise, you need None (need to fix here)
                row_method = None
                column_method = None
                row_method = 'hopach'
                column_method = 'hopach'
    except Exception,e:
        #print traceback.format_exc()
        pass
    
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
    if color_gradient == 'Greys':
        cmap=pylab.cm.Greys
    if color_gradient == 'yellow_orange_red':
        cmap=pylab.cm.YlOrRd
        
    vmin=x.min()
    vmax=x.max()
    vmax = max([vmax,abs(vmin)])
    if Normalize != False:
        vmin = vmax*-1
    elif 'Clustering-Zscores-' in dataset_name:
        vmin = vmax*-1
    elif vmin<0 and vmax>0 and Normalize==False:
        vmin = vmax*-1
    #vmin = vmax*-1
    #print vmax, vmin
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
    norm = matplotlib.colors.Normalize(vmin/scaling_factor, vmax/scaling_factor) ### adjust the max and min to scale these colors by 2.5 (1 scales to the highest change)
    #fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    fig = pylab.figure() ### could use m,n to scale here - figsize=(12,10)
    fig.set_figwidth(12)
    fig.set_figheight(7)
    fig.patch.set_facecolor('white')
    pylab.rcParams['font.size'] = 7.5
    #pylab.rcParams['axes.facecolor'] = 'white' ### Added in 2.1.2

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
            matrix_horiz_pos = 0.27
        elif skipClustering:
            if len(row_header)<100:
                matrix_horiz_pos = 0.20
            else:
                matrix_horiz_pos = 0.27
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
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.02,0.938,0.17,0.025] ### Last one controls the hight [0.07,0.88,0.18,0.076]
    
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
            elif len(PriorColumnClusters)>0 and len(PriorRowClusters)>0 and row_method==None and column_method == None:
                no_plot = True ### If trying to instantly view prior results, no dendrogram will be display, but prior GO-Elite can
            else:
                ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False - this window conflicts with GO-Elite labels
        except Exception:
            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
        try: Z1 = sch.dendrogram(Y1, orientation='left',no_plot=no_plot) ### This is where plotting occurs - orientation 'right' in old matplotlib
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
        """ Error can occur here if hopach was selected in a prior run but now running NONE """
        try: ind2 = [ind2[i] for i in idx2] ### replaces the above due to numpy specific windows version issue
        except Exception:
            column_method=None
            xt = x
            ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
            ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
            
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        prior_xt = xt
        xt = xt[idx1,:]   # xt is transformed x
        #ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
        try: ind1 = [ind1[i] for i in idx1] ### replaces the above due to numpy specific windows version issue
        except Exception:
            if 'MarkerGenes' in dataset_name:
                ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
                row_method = None
                
        
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
            import gene_associations
            try:
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            except Exception: gene_to_symbol={}; symbol_to_gene={}
        JustShowTheseIDs = copy.deepcopy(justShowTheseIDs)
    except Exception:
        JustShowTheseIDs=[]

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
    if 'MarkerGenes' in dataset_name:
        time_stamp = timestamp() ### Don't overwrite the previous version
        filename = string.replace(filename,'hierarchical',time_stamp)

    elite_dir, cdt_file, SystemCode = exportFlatClusterData(root_dir + filename, root_dir, dataset_name, new_row_header,new_column_header,xt,ind1,ind2,vmax,display)

    def ViewPNG(png_file_dir):
        if os.name == 'nt':
            try: os.startfile('"'+png_file_dir+'"')
            except Exception:  os.system('open "'+png_file_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')
        
    try:
        try:
            temp1=len(justShowTheseIDs)
            if 'monocle' in justShowTheseIDs and ('guide' not in justShowTheseIDs):
                import R_interface
                print 'Running Monocle through R (be patient, this can take 20 minutes+)'
                R_interface.performMonocleAnalysisFromHeatmap(species,cdt_file[:-3]+'txt',cdt_file[:-3]+'txt')
                png_file_dir = root_dir+'/Monocle/monoclePseudotime.png'
                #print png_file_dir
                ViewPNG(png_file_dir)
        except Exception: pass # no justShowTheseIDs
    except Exception:
        print '...Monocle error:'
        print traceback.format_exc()
        pass
    
    cluster_elite_terms={}; ge_fontsize=11.5; top_genes=[]; proceed=True
    try:
        try:
            if 'guide' in justShowTheseIDs: proceed = False
        except Exception: pass
        if proceed:
            try:
                cluster_elite_terms,top_genes = remoteGOElite(elite_dir,SystemCode=SystemCode)
                if cluster_elite_terms['label-size']>40: ge_fontsize = 9.5
            except Exception:
                pass

    except Exception: pass  #print traceback.format_exc()
    if len(cluster_elite_terms)<1:
        try:
            elite_dirs = string.split(elite_dir,'GO-Elite')
            old_elite_dir = elite_dirs[0]+'GO-Elite'+elite_dirs[-1] ### There are actually GO-Elite/GO-Elite directories for the already clustered
            old_elite_dir = string.replace(old_elite_dir,'ICGS/','')
            if len(PriorColumnClusters)>0 and len(PriorRowClusters)>0 and skipClustering:
                cluster_elite_terms,top_genes = importGOEliteResults(old_elite_dir)
        except Exception,e:
            #print traceback.format_exc()
            pass
    try:
        if len(justShowTheseIDs)<1 and len(top_genes) > 0 and column_fontsize < 9:
            column_fontsize = 10
        if len(justShowTheseIDs)<1:
            additional_symbols=[]
            import gene_associations; from import_scripts import OBO_import
            try:
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                #symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
            except Exception: gene_to_symbol={}; symbol_to_gene={}
    except Exception: pass
    
    def formatpval(p):
        if '-' in p: p1=p[:1]+p[-4:]
        else:
            p1 = '{number:.{digits}f}'.format(number=float(p), digits=3)
            p1=str(p1)
        #print traceback.format_exc();sys.exit()
        return p1
    
    # Add text
    new_row_header=[]
    new_column_header=[]
    ci=0 ### index of entries in the cluster
    last_cluster=1
    interval = int(float(string.split(str(len(row_header)/35.0),'.')[0]))+1 ### for enrichment term labels with over 100 genes
    increment=interval-2
    if len(row_header)<100: increment = interval-1
    label_pos=-0.03*len(column_header)-.5
    #print label_pos
    try:
        if 'top' in justShowTheseIDs: justShowTheseIDs.remove('top')
        if 'positive' in justShowTheseIDs: justShowTheseIDs.remove('positive')
        if 'amplify' in justShowTheseIDs: justShowTheseIDs.remove('amplify')
        if 'IntraCorrelatedOnly' in justShowTheseIDs: justShowTheseIDs.remove('IntraCorrelatedOnly')
        if 'GuideOnlyCorrelation' in justShowTheseIDs: justShowTheseIDs.remove('GuideOnlyCorrelation')
    except Exception:
        pass

    for i in range(x.shape[0]):
        if len(row_header)<40:
            radj = len(row_header)*0.009 ### row offset value to center the vertical position of the row label
        elif len(row_header)<70:
            radj = len(row_header)*0.007 ### row offset value to center the vertical position of the row label
        else:
            radj = len(row_header)*0.005
        try: cluster = str(ind1[i])
        except Exception: cluster = 'NA'
        if cluster == 'NA':
            new_index = i
            try: cluster = 'cluster-'+string.split(row_header[new_index],':')[0]
            except Exception: pass

        if cluster != last_cluster:
            ci=0
            increment=0
        #print cluster,i,row_header[idx1[i]]
        color = 'black'
        
        if row_method != None:
            try:
                if row_header[idx1[i]] in JustShowTheseIDs:
                    if len(row_header)>len(justShowTheseIDs):
                        color = 'red'
                else: color = 'black'
            except Exception: pass
            if len(row_header)<106: ### Don't visualize gene associations when more than 100 rows
                if display_label_names == False or 'ticks' in JustShowTheseIDs:
                    if color=='red':
                        axm.text(x.shape[1]-0.5, i-radj, '  '+'-',fontsize=row_fontsize, color=color, picker=True)
                    else:
                        axm.text(x.shape[1]-0.5, i-radj, '  '+'',fontsize=row_fontsize, color=color, picker=True)
                else:
                    axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[idx1[i]],fontsize=row_fontsize, color=color, picker=True)
            new_row_header.append(row_header[idx1[i]])
            new_index = idx1[i]
        else:
            try:
                if row_header[i] in JustShowTheseIDs: color = 'red'
                else: color = 'black'
            except Exception: pass
            if len(row_header)<106: ### Don't visualize gene associations when more than 100 rows
                if display_label_names == False or 'ticks' in JustShowTheseIDs:
                    if color=='red':
                        axm.text(x.shape[1]-0.5, i-radj, '  '+'-',fontsize=row_fontsize, color=color, picker=True)
                    else:
                        axm.text(x.shape[1]-0.5, i-radj, '  '+'',fontsize=row_fontsize, color=color, picker=True)
                else:
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
            original_feature_id = feature_id
            if ':' in feature_id:
                if 'ENS' != feature_id[:3] or 'G0000' in feature_id:
                    feature_id = string.split(feature_id,':')[1]
                else:
                    feature_id = string.split(feature_id,':')[0]
                    try: feature_id = gene_to_symbol[feature_id][0]
                    except Exception: pass
            if (' ' in feature_id and ('ENS' in feature_id or 'G0000' in feature_id)):
                feature_id = string.split(feature_id,' ')[1]
            try:
                if feature_id in JustShowTheseIDs or original_feature_id in JustShowTheseIDs: color = 'red'
                else: color = 'black'
            except Exception: pass
            try:
                if feature_id in justShowTheseIDs or (len(justShowTheseIDs)<1 and feature_id in top_genes) or original_feature_id in justShowTheseIDs:
                    if original_feature_id in justShowTheseIDs:
                        feature_id = original_feature_id
                    if display_label_names and 'ticks' not in justShowTheseIDs:
                        axm.text(x.shape[1]-0.5, i-radj, '  '+feature_id,fontsize=column_fontsize, color=color,picker=True) ### When not clustering rows
                    else:
                        axm.text(x.shape[1]-0.5, i-radj, '  '+"-",fontsize=column_fontsize, color=color,picker=True) ### When not clustering rows
                elif ' ' in row_header[new_index]:
                    symbol = string.split(row_header[new_index], ' ')[-1]
                    if symbol in justShowTheseIDs:
                        if display_label_names and 'ticks' not in justShowTheseIDs:
                            axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[new_index],fontsize=column_fontsize, color=color,picker=True)
                        else:
                            axm.text(x.shape[1]-0.5, i-radj, '  '+"-",fontsize=column_fontsize, color=color,picker=True)
            except Exception: pass
                    
        if cluster in cluster_elite_terms or 'cluster-'+cluster in cluster_elite_terms:
                if 'cluster-'+cluster in cluster_elite_terms:
                    new_cluster_id = 'cluster-'+cluster
                else:
                    new_cluster_id = cluster
                if cluster != last_cluster:
                    cluster_intialized = False
                try:
                    increment+=1
                    #print [increment,interval,cluster],cluster_elite_terms[cluster][ci][1];sys.exit()
                    #if increment == interval or (
                    #print increment,interval,len(row_header),cluster_intialized
                    if (increment == interval) or (len(row_header)>200 and increment == (interval-9) and cluster_intialized==False): ### second argument brings the label down
                        cluster_intialized=True
                        atypical_cluster = False
                        if ind1[i+9] == 'NA': ### This occurs for custom cluster, such MarkerFinder (not cluster numbers)
                            atypical_cluster = True
                            cluster9 = 'cluster-'+string.split(row_header[new_index+9],':')[0]
                            if (len(row_header)>200 and str(cluster9)!=cluster): continue
                        elif (len(row_header)>200 and str(ind1[i+9])!=cluster): continue ### prevents the last label in a cluster from overlapping with the first in the next cluster
                        pvalue,original_term = cluster_elite_terms[new_cluster_id][ci]
                        term = original_term
                        if 'GO:' in term:
                            term = string.split(term, '(')[0]
                        if ':WP' in term:
                            term = string.split(term, ':WP')[0]
                        pvalue = formatpval(str(pvalue))
                        term += ' p='+pvalue
                        if atypical_cluster == False:
                            term += ' (c'+str(cluster)+')'
                        try: cluster_elite_terms[term] = cluster_elite_terms[cluster,original_term]  ### store the new term name with the associated genes
                        except Exception: pass
                        axm.text(label_pos, i-radj, term,horizontalalignment='right',fontsize=ge_fontsize, picker=True, color = 'blue')
                        increment=0
                        ci+=1
                except Exception,e:
                    #print traceback.format_exc();sys.exit()
                    increment=0
        last_cluster = cluster
    
    def onpick1(event):
        text = event.artist
        print('onpick1 text:', text.get_text())
        if 'TreeView' in text.get_text():
            try: openTreeView(cdt_file)
            except Exception: print 'Failed to open TreeView'
        elif 'p=' not in text.get_text():
            webbrowser.open('http://www.genecards.org/cgi-bin/carddisp.pl?gene='+string.replace(text.get_text(),' ',''))
        else:
            #"""
            from visualization_scripts import TableViewer
            header = ['Associated Genes']
            tuple_list = []
            
            for gene in cluster_elite_terms[text.get_text()]:
                tuple_list.append([(gene)])
            TableViewer.viewTable(text.get_text(),header,tuple_list) #"""
            
            cluster_prefix = 'c'+string.split(text.get_text(),'(c')[1][:-1]+'-'

            for geneSet in EliteGeneSets:
                if geneSet == 'GeneOntology':
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+'GO'+'.png'
                elif geneSet == 'WikiPathways':
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+'local'+'.png'
                elif len(geneSet)>1:
                    png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+geneSet+'.png'
            #try: UI.GUI(root_dir,'ViewPNG',[],png_file_dir)
            #except Exception: print traceback.format_exc()

            try:
                alt_png_file_dir = elite_dir+'/GO-Elite_results/networks/'+cluster_prefix+eliteGeneSet+'.png'
                png_file_dirs = string.split(alt_png_file_dir,'GO-Elite/')
                alt_png_file_dir = png_file_dirs[0]+'GO-Elite/'+png_file_dirs[-1]
            except Exception: pass  
            if os.name == 'nt':
                try: os.startfile('"'+png_file_dir+'"')
                except Exception:
                    try: os.system('open "'+png_file_dir+'"')
                    except Exception: os.startfile('"'+alt_png_file_dir+'"')
            elif 'darwin' in sys.platform:
                try: os.system('open "'+png_file_dir+'"')
                except Exception: os.system('open "'+alt_png_file_dir+'"')
            elif 'linux' in sys.platform:
                try: os.system('xdg-open "'+png_file_dir+'"')
                except Exception: os.system('xdg-open "'+alt_png_file_dir+'"')

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
        if len(column_header)>200:
            column_fontsize = 2
        if column_method != None:
            if len(column_header)<300: ### Don't show the headers when too many values exist
                axm.text(adji, cadj, ''+column_header[idx2[i]], rotation=270, verticalalignment="top",fontsize=column_fontsize) # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            if len(column_header)<300: ### Don't show the headers when too many values exist
                axm.text(adji, cadj, ''+column_header[i], rotation=270, verticalalignment="top",fontsize=column_fontsize)
            new_column_header.append(column_header[i])

    # Plot colside colors
    # axc --> axes for column side colorbar

    group_name_list=[]
    ind1_clust,ind2_clust = ind1,ind2
    ind1,ind2,group_name_list,cb_status = updateColorBarData(ind1,ind2,new_column_header,new_row_header,row_method)

    if (column_method != None or 'column' in cb_status) and show_color_bars == True:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
        if use_default_colors:
            cmap_c = pylab.cm.gist_rainbow
        else:
            #cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
            if len(unique.unique(ind2))==2: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                #cmap_c = matplotlib.colors.ListedColormap(['#7CFC00','k'])
                cmap_c = matplotlib.colors.ListedColormap(['w', 'k'])
            elif len(unique.unique(ind2))==3: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                #cmap_c = matplotlib.colors.ListedColormap(['r', 'y', 'b'])
            elif len(unique.unique(ind2))==4: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C','#FEBC18'])
                #cmap_c = matplotlib.colors.ListedColormap(['k', 'w', 'w', 'w'])
            elif len(unique.unique(ind2))==5: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))==6: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                #cmap_c = matplotlib.colors.ListedColormap(['black', '#1DA532', '#88BF47','b', 'grey','r'])
                #cmap_c = matplotlib.colors.ListedColormap(['w', '#0B9B48', 'w', '#5D82C1','#4CB1E4','#71C065'])
                #cmap_c = matplotlib.colors.ListedColormap(['w', 'w', 'k', 'w','w','w'])
            elif len(unique.unique(ind2))==7: ### cmap_c is too few colors
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                #cmap_c = matplotlib.colors.ListedColormap(['w', 'w', 'w', 'k', 'w','w','w'])
                #cmap_c = matplotlib.colors.ListedColormap(['w','w', '#0B9B48', 'w', '#5D82C1','#4CB1E4','#71C065'])
            #elif len(unique.unique(ind2))==9:  cmap_c = matplotlib.colors.ListedColormap(['k', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'])
            #elif len(unique.unique(ind2))==11: 
            #cmap_c = matplotlib.colors.ListedColormap(['w', '#DC2342', '#0B9B48', '#FDDF5E', '#E0B724', 'k', '#5D82C1', '#F79020', '#4CB1E4', '#983894', '#71C065'])
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
                cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
                #cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
                if use_default_colors:
                    cmap_c = pylab.cm.gist_rainbow
                else:
                    if len(unique.unique(ind2_clust))==2: ### cmap_c is too few colors
                        #cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                        cmap_c = matplotlib.colors.ListedColormap(['w', 'k'])
                    elif len(unique.unique(ind2_clust))==3: ### cmap_c is too few colors
                        cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                        #cmap_c = matplotlib.colors.ListedColormap(['r', 'y', 'b'])
                    elif len(unique.unique(ind2_clust))==4: ### cmap_c is too few colors
                        cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
                        #cmap_c = matplotlib.colors.ListedColormap(['black', '#1DA532', 'b','r'])
                    elif len(unique.unique(ind2_clust))==5: ### cmap_c is too few colors
                        cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
                    elif len(unique.unique(ind2_clust))==6: ### cmap_c is too few colors
                        cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                    elif len(unique.unique(ind2_clust))==7: ### cmap_c is too few colors
                        cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
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
            cmap_d = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0','#000066','#FFFF00', '#FF1493'])
            #cmap_d = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF','#FFFF00', '#FF1493'])
            if len(unique.unique(ind2))==2: ### cmap_c is too few colors
                #cmap_d = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                cmap_d = matplotlib.colors.ListedColormap(['w', 'k'])
            elif len(unique.unique(ind2))==3: ### cmap_c is too few colors
                cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                #cmap_d = matplotlib.colors.ListedColormap(['r', 'y', 'b'])
            elif len(unique.unique(ind2))==4: ### cmap_c is too few colors
                cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif len(unique.unique(ind2))==5: ### cmap_c is too few colors
                cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2))==6: ### cmap_c is too few colors
                cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                #cmap_d = matplotlib.colors.ListedColormap(['black', '#1DA532', '#88BF47','b', 'grey','r'])
                #cmap_d = matplotlib.colors.ListedColormap(['w', '#0B9B48', 'w', '#5D82C1','#4CB1E4','#71C065'])
                #cmap_d = matplotlib.colors.ListedColormap(['w', 'w', 'k', 'w', 'w','w','w'])
            elif len(unique.unique(ind2))==7: ### cmap_c is too few colors
                cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
                #cmap_d = matplotlib.colors.ListedColormap(['w', 'w', 'w', 'k', 'w','w','w'])
                #cmap_d = matplotlib.colors.ListedColormap(['w','w', '#0B9B48', 'w', '#5D82C1','#4CB1E4','#71C065'])
            #elif len(unique.unique(ind2))==10: cmap_d = matplotlib.colors.ListedColormap(['w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'k'])
            #elif len(unique.unique(ind2))==11:
            #Eryth Gfi1 Gran HSCP-1 HSCP-2 IG2 MDP Meg Mono Multi-Lin Myelo
            #cmap_d = matplotlib.colors.ListedColormap(['#DC2342', 'k', '#0B9B48', '#FDDF5E', '#E0B724', 'w', '#5D82C1', '#F79020', '#4CB1E4', '#983894', '#71C065'])
            elif len(unique.unique(ind2))>0: ### cmap_c is too few colors
                cmap_d = pylab.cm.gist_rainbow
            dc = numpy.array(group_colors, dtype=int)
            dc.shape = (1,len(group_colors)) 
            im_c = axd.matshow(dc, aspect='auto', origin='lower', cmap=cmap_d)
            axd.set_yticks([])
            #axd.set_xticklabels(group_names, rotation=45, ha='left')
            #if len(group_names)<200:
            pylab.xticks(range(len(group_names)),group_names,rotation=45,ha='left')
            #cmap_c = matplotlib.colors.ListedColormap(map(lambda x: GroupDB[x][-1], new_column_header))

    if show_color_bars == False:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        axc.set_frame_on(False)
    
    # Plot rowside colors
    # axr --> axes for row side colorbar
    if (row_method != None or 'row' in cb_status) and show_color_bars == True:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        try:
            dr = numpy.array(ind1, dtype=int)
            dr.shape = (len(ind1),1)
            #print ind1, len(ind1)
            cmap_r = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#FFFF00', '#FF1493'])
            if len(unique.unique(ind1))>4: ### cmap_r is too few colors
                cmap_r = pylab.cm.gist_rainbow_r
            if len(unique.unique(ind1))==2:
                cmap_r = matplotlib.colors.ListedColormap(['w', 'k'])
            im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
            axr.set_xticks([]) ### Hides ticks
            axr.set_yticks([])
            #axr.set_frame_on(False) ### Hide border
        except Exception:
            row_method = None
            pass ### likely occurs for some reason when row_method should be None
        
    if show_color_bars == False:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        axr.set_frame_on(False)

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = matplotlib.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
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
    
    includeBackground=False
    try:
        if 'TkAgg' != matplotlib.rcParams['backend']:
            includeBackground = False
    except Exception: pass
    
    if includeBackground:
        fig.text(0.020, 0.070, 'Open heatmap in TreeView (click here)', fontsize = 11.5, picker=True,color = 'red', backgroundcolor='white')
    else:
        fig.text(0.020, 0.070, 'Open heatmap in TreeView (click here)', fontsize = 11.5, picker=True,color = 'red')
    if 'Outlier' in dataset_name and 'Removed' not in dataset_name:
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
    elif 'alt_junction' in filename:
        graphic_link.append(['Hierarchical Clustering - Variable Splice-Events',root_dir+filename])
    else:
        graphic_link.append(['Hierarchical Clustering - Significant Genes',root_dir+filename])
    
    if display:
        proceed=True
        try:
            if 'guide' in justShowTheseIDs:
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
    retcode = subprocess.Popen(['java', "-Xmx500m", '-jar', fn, "-r", filename])

def remoteGOElite(elite_dir,SystemCode = None):
    mod = 'Ensembl'

    if SystemCode == 'Ae':    
        mod = 'AltExon'

    pathway_permutations = 'FisherExactTest'
    filter_method = 'z-score'
    z_threshold = 1.96
    p_val_threshold = 0.05
    change_threshold = 0
    if runGOElite:
        resources_to_analyze = EliteGeneSets
        if 'all' in resources_to_analyze:
            resources_to_analyze = 'all'
        returnPathways = 'no'
        root = None
        import GO_Elite
        reload(GO_Elite)
        
        input_files = dir_list = unique.read_directory(elite_dir) ### Are there any files to analyze?
        if len(input_files)>0 and resources_to_analyze !=['']:
            print '\nBeginning to run GO-Elite analysis on all results' 
            file_dirs = elite_dir,None,elite_dir
            enrichmentAnalysisType = 'ORA'
            #enrichmentAnalysisType = 'URA'
            variables = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,enrichmentAnalysisType,root
            try: GO_Elite.remoteAnalysis(variables,'non-UI Heatmap')
            except Exception: 'GO-Elite failed for:',elite_dir
            if commandLine==False:
                try: UI.openDirectory(elite_dir+'/GO-Elite_results') 
                except Exception: None
            cluster_elite_terms,top_genes = importGOEliteResults(elite_dir)
            return cluster_elite_terms,top_genes
        else:
            return {},[]
    else:
        return {},[]
    
def importGOEliteResults(elite_dir):
    global eliteGeneSet
    pruned_results = elite_dir+'/GO-Elite_results/CompleteResults/ORA_pruned/pruned-results_z-score_elite.txt' ### This is the exception (not moved)
    if os.path.isfile(pruned_results) == False:
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
                eliteGeneSet = string.split(values[0][1:],'-')[1][:-4]
                try: cluster = str(int(float(string.split(values[0][1:],'-')[0])))
                except Exception:
                    cluster = string.join(string.split(values[0],'-')[:-1],'-')
                term = values[2]
                num_genes_changed = int(values[3])
                all_term_length.append(len(term))
                pval = float(values[9]) ### adjusted is values[10]
                #pval = float(values[10]) ### adjusted is values[10]
                if num_genes_changed>2:
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

"""
def merge_horizontal(out_filename, left_filename, right_filename):
    #Merge the first page of two PDFs side-to-side
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
"""     
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

def exportFlatClusterData(filename, root_dir, dataset_name, new_row_header,new_column_header,xt,ind1,ind2,vmax,display):
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
        try:
            id = new_row_header[i]
            original_id = str(id)
            if sy == '$En:Sy':
                cluster = 'cluster-'+string.split(id,':')[0]
            elif sy == 'S' and ':' in id:
                cluster = 'cluster-'+string.split(id,':')[0]
            elif sy == 'Sy' and ':' in id:
                cluster = 'cluster-'+string.split(id,':')[0]
            else:
                cluster = 'c'+str(ind1[i])
        except Exception:
            pass
    
        try:
            if 'MarkerGenes' in originalFilename:
                cluster = 'cluster-'+string.split(id,':')[0]
                id = string.split(id,':')[1]
                if ' ' in id:
                    id = string.split(id,' ')[0]
                if 'G000' in id: sy = 'En'
                else: sy = 'Sy'
        except Exception: pass
        try: cluster_db[cluster].append(id)
        except Exception: cluster_db[cluster] = [id]
        try: export_lines.append(string.join([original_id,str(ind1[i])]+map(str, row),'\t')+'\n')
        except Exception:
            export_lines.append(string.join([original_id,'NA']+map(str, row),'\t')+'\n')
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
            try:
                i1,i2 = string.split(id,' ')
                if i1==i2: id = i1
            except Exception: pass  
            if sy == '$En:Sy':
                id = string.split(id,':')[1]
                ids = string.split(id,' ')
                if 'ENS' in ids[0] or 'G0000' in ids[0]: id = ids[0]
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
                    if 'ENS' in ids[-1] or 'G0000' in ids[-1]: id = ids[-1]
                    else: id = ids[0]
            elif sy == 'En' and '&' in id:
                for i in string.split(id,'&'):
                    if 'G0000' in i: id = i; sc = 'En'; break
            elif sy == 'Sy' and 'EFN' in id:
                sc = 'En'
            else:
                sc = sy
            if sy == 'S':
                if ':' in id:
                    id = string.split(id,':')[-1]
                    sc = 'Ae'
            if '&' in id:
                sc = 'Ae'
            if (len(id)==9 and "SRS" in id) or (len(id)==15 and "TCGA-" in id):
                sc = 'En'
            try: export_elite.write(id+'\t'+sc+'\n')
            except Exception: export_elite.write(id+'\n') ### if no System Code known
            allGenes[id]=[]
        export_elite.close()
    
    try:
        if storeGeneSetName != None:
            if len(storeGeneSetName)>0 and ('guide' not in justShowTheseIDs):
                exportCustomGeneSet(storeGeneSetName,species,allGenes)
                print 'Exported geneset to "StoredGeneSets"'
    except Exception: pass
    
    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    if display:
        try: exportJTV(filename, new_column_header, new_row_header,vmax=vmax)
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
    return elite_dir, filename, sc

def exportJTV(cdt_dir, column_header, row_header,vmax=None):
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
    contrast = str(float(vmax)/4)[:4] ### base the contrast on the heatmap vmax variable
    """
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
    """
    
    config = '<DocumentConfig><UrlExtractor/><ArrayUrlExtractor/><MainView><ColorExtractor>'
    config+= '<ColorSet down="#00FFFF"/></ColorExtractor><ArrayDrawer/><GlobalXMap>'
    config+= '<FixedMap type="Fixed" scale="'+cscale+'"/><FillMap type="Fill"/><NullMap type="Null"/>'
    config+= '</GlobalXMap><GlobalYMap><FixedMap type="Fixed" scale="'+hscale+'"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></GlobalYMap><ZoomXMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></ZoomXMap><ZoomYMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></ZoomYMap><TextView><TextView><GeneSummary/></TextView><TextView>'
    config+= '<GeneSummary/></TextView><TextView><GeneSummary/></TextView></TextView><ArrayNameView>'
    config+= '<ArraySummary included="0"/></ArrayNameView><AtrSummary/><GtrSummary/></MainView><Views>'
    config+= '<View type="Dendrogram" dock="1"><ColorExtractor contrast="'+contrast+'"><ColorSet up="#FFFF00" down="#00CCFF"/>'
    config+= '</ColorExtractor><ArrayDrawer/><GlobalXMap current="Fill"><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></GlobalXMap><GlobalYMap current="Fill"><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config+= '<NullMap type="Null"/></GlobalYMap><ZoomXMap><FixedMap type="Fixed"/><FillMap type="Fill"/><NullMap type="Null"/>'
    config+= '</ZoomXMap><ZoomYMap current="Fixed"><FixedMap type="Fixed"/><FillMap type="Fill"/><NullMap type="Null"/></ZoomYMap>'
    config+= '<TextView><TextView><GeneSummary/></TextView><TextView><GeneSummary/></TextView><TextView><GeneSummary/></TextView>'
    config+= '</TextView><ArrayNameView><ArraySummary included="0"/></ArrayNameView><AtrSummary/><GtrSummary/></View></Views></DocumentConfig>'
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


def ConvertFromHex(color1,color2,color3):
    c1tuple = tuple(ord(c) for c in color1.lsstrip('0x').decode('hex'))
    c2tuple = tuple(ord(c) for c in color2.lsstrip('0x').decode('hex'))
    c3tuple = tuple(ord(c) for c in color3.lsstrip('0x').decode('hex'))

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

def remoteImportData(filename,geneFilter=None,reverseOrder=True):
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,geneFilter=geneFilter,reverseOrder=reverseOrder)
    try:
        return matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters
    except:
        return matrix, column_header, row_header, dataset_name, group_db, [], []

def importData(filename,Normalize=False,reverseOrder=True,geneFilter=None,zscore=False):
    
    global priorColumnClusters 
    global priorRowClusters
    try: 
        if len(priorColumnClusters)>0:
            priorColumnClusters = None
            priorRowClusters = None
    except Exception: pass
    getRowClusters=False
    start_time = time.time()
    fn = filepath(filename)
    matrix=[]
    original_matrix=[]
    row_header=[]
    overwriteGroupNotations=True
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
            if t[1] == 'row_clusters-flat':
                t = t = [t[0]]+t[2:]
            ### color samples by annotated groups if an expression file
            new_headers=[]
            temp_groups={}
            original_headers=t[1:]
            if ('exp.' in filename or 'filteredExp.' in filename or 'MarkerGene' in filename):# and ':' not in data:
                if overwriteGroupNotations:
                    ### Use groups file annotations over any header sample separation with a ":"
                    for i in t:
                        if ':' in i:
                            group,i = string.split(i,':')
                            new_headers.append(i)
                            temp_groups[i] = group
                        else: new_headers.append(i)
                filename = string.replace(filename,'-steady-state.txt','.txt')
                try:
                    import ExpressionBuilder
                    try: sample_group_db = ExpressionBuilder.simplerGroupImport(filename)
                    except Exception: sample_group_db={}
                    if len(temp_groups)>0 and len(sample_group_db)==0:
                        sample_group_db = temp_groups

                    if len(new_headers)>0:
                        t = new_headers
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
            try:
                if 'NA' in t:
                    kill
                try: prior = map(lambda x: int(float(x)),t[2:])
                except Exception:
                    ### Replace the cluster string with number
                    index=0
                    c=1; prior=[]; clusters={}
                    for i in t[2:]:
                        original_headers[index] = i+':'+original_headers[index]
                        if i in clusters:
                            c1 = clusters[i]
                        else:
                            c1 = c; clusters[i]=c1
                            c+=1
                        prior.append(c1)
                        index+=1
                    #prior=[]
                    if len(temp_groups)==0: ### Hence, no secondary group label combined with the sample name
                        group_db, column_header = assignGroupColors(original_headers)
                #priorColumnClusters = dict(zip(column_header, prior))
                priorColumnClusters = prior
            except Exception:
                #print traceback.format_exc() 
                pass
            start = 2
            getRowClusters = True
            priorRowClusters=[]
        elif 'EWEIGHT' in t: pass
        else:
            gene = t[0]
            if geneFilter==None:
                proceed = True
            elif gene in geneFilter:
                proceed = True
            else:
                proceed = False
            if proceed:
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
                try:
                    if max(s)>inputMax: inputMax = max(s)
                except:
                    continue ### empty row
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
                
                if ' ' in gene:
                    try:
                        g1,g2 = string.split(gene,' ')
                        if g1 == g2: gene = g1
                    except Exception: pass
                    
                if getRowClusters:
                    try:
                        #priorRowClusters[gene]=int(float(t[1]))
                        priorRowClusters.append(int(float(t[1])))
                    except Exception: pass
                
                if zscore:
                    ### convert to z-scores for normalization prior to PCA
                    avg = numpy.mean(s)
                    std = numpy.std(s)
                    if std ==0:
                        std = 0.1
                    try: s = map(lambda x: (x-avg)/std,s)
                    except Exception: pass
                    
                if geneFilter==None:
                    matrix.append(s)
                    row_header.append(gene)
                else:
                    if gene in geneFilter:
                        matrix.append(s)
                        row_header.append(gene)
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
            if zscore: ### The above z-score does not impact the original_matrix which is analyzed
                ### convert to z-scores for normalization prior to PCA
                avg = numpy.mean(s)
                std = numpy.std(s)
                if std ==0:
                    std = 0.1
                try: s = map(lambda x: (x-avg)/std,s)
                except Exception: pass
            matrix.append(s)
            k+=1
        del original_matrix

    if zscore: print 'Converting values to normalized z-scores...'
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
    groupNamesPresent=False # Some samples may have missing group names which will result in a clustering error
    for i in t:
        if ':' in i: groupNamesPresent = True

    for i in t:
        repls = {'.2txt' : '', '.2bed' : '', '.2tab' : ''}
        i=reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        if ':' in i:
            group,j = string.split(i,':')[:2]
            group_number_db[group]=[]
        elif groupNamesPresent:
            group_number_db['UNK']=[]
            i = 'UNK:'+i
        column_header.append(i)
            
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
    for i in column_header:
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
            import gene_associations; from import_scripts import OBO_import
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

def writetSNEScores(scores,outputdir):
    export_obj = export.ExportFile(outputdir)
    for matrix_row in scores:
        matrix_row = map(str,matrix_row)
        export_obj.write(string.join(matrix_row,'\t')+'\n')
    export_obj.close()

def importtSNEScores(inputdir):
    #print inputdir
    scores=[]
    ### Imports tSNE scores to allow for different visualizations of the same scatter plot
    for line in open(inputdir,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        t=map(float,t)
        scores.append(t)
    return scores

def runUMAP(matrix, column_header,dataset_name,group_db,display=False,showLabels=False,
         row_header=None,colorByGene=None,species=None,reimportModelScores=True,method="UMAP",rootDir='',finalOutputDir=''):
    global root_dir
    global graphic_link
    graphic_link=[]
    root_dir = rootDir
    
    tSNE(matrix, column_header,dataset_name,group_db,display=False,showLabels=False,
         row_header=None,colorByGene=None,species=None,reimportModelScores=True,method="UMAP")
    import shutil
    filename = 'Clustering-'+dataset_name+'-'+method+'.pdf'
    filename = string.replace(filename,'Clustering-Clustering','Clustering')
    new_file=finalOutputDir + filename
    new_file=string.replace(new_file,'Clustering-','')
    new_file=string.replace(new_file,'exp.','')
    old_file=root_dir+filename
    shutil.move(old_file,new_file)
    
    filename = filename[:-3]+'png'
    new_file=finalOutputDir + filename
    new_file=string.replace(new_file,'Clustering-','')
    new_file=string.replace(new_file,'exp.','')
    old_file=root_dir+filename
    shutil.move(old_file,new_file)
    old_file=root_dir+dataset_name+'-'+method+'_scores.txt'
    new_file=finalOutputDir+dataset_name+'-'+method+'_coordinates.txt'
    new_file=string.replace(new_file,'exp.','')
    shutil.move(old_file,new_file)

def tSNE(matrix, column_header,dataset_name,group_db,display=True,showLabels=False,
         row_header=None,colorByGene=None,species=None,reimportModelScores=True,method="tSNE"):

    try: prior_clusters = priorColumnClusters
    except Exception: prior_clusters=[]
    try:
        if priorColumnClusters==None: prior_clusters=[]
    except:
        pass
    try:
        if len(prior_clusters)>0 and len(group_db)==0:
            newColumnHeader=[]
            i=0
            for sample_name in column_header:
                newColumnHeader.append(str(prior_clusters[i])+':'+sample_name)
                i+=1
            group_db, column_header = assignGroupColors(newColumnHeader)    
    except Exception,e:
        print traceback.format_exc()
        #print e
        group_db={}

    if reimportModelScores:
        print 'Re-importing',method,'model scores rather than calculating from scratch',
        print root_dir+dataset_name+'-'+method+'_scores.txt'
        try: scores = importtSNEScores(root_dir+dataset_name+'-'+method+'_scores.txt'); print '...import finished'
        except Exception:
            reimportModelScores=False; print '...import failed'
        
    if reimportModelScores==False:
        X=matrix.T
        """
        from tsne import bh_sne
        X = np.asarray(X).astype('float64')
        X = X.reshape((X.shape[0], -1))
        x_data = x_data.reshape((x_data.shape[0], -1))
        scores = bh_sne(X)"""
        
        #model = TSNE(n_components=2, random_state=0,init='pca',early_exaggeration=4.0,perplexity=20)
        print "Performing",method
        if method=="tSNE" or method=="t-SNE":
            from sklearn.manifold import TSNE
            model = TSNE(n_components=2)
        if method=="UMAP":
            try:
                import umap
            except: 
                from visualization_scripts.umap_learn import umap ### Bypasses issues with Py2app importing umap (and secondarily numba/llvmlite)
            model=umap.UMAP(n_neighbors=50,min_dist=0.75,metric='correlation')
            print 'UMAP run'
        #model = TSNE(n_components=2,init='pca', random_state=0, verbose=1, perplexity=40, n_iter=300)
        #model = TSNE(n_components=2,verbose=1, perplexity=40, n_iter=300)
        #model = TSNE(n_components=2, random_state=0, n_iter=10000, early_exaggeration=10)
        scores=model.fit_transform(X)
        
        ### Export the results for optional re-import later
        writetSNEScores(scores,root_dir+dataset_name+'-'+method+'_scores.txt')
        #pylab.scatter(scores[:,0], scores[:,1], 20, labels);
    
    ### Exclude samples with high TSNE deviations
    scoresT = zip(*scores)
    exclude={}
    
    try:
        for vector in scoresT:
            lower1th,median_val,upper99th,int_qrt_range = statistics.iqr(list(vector),k1=99.9,k2=0.1)
            index=0
            for i in vector:
                if (i > upper99th+1) or (i<lower1th-1):
                    exclude[index]=None
                index+=1    
    except Exception:
        pass

    print 'Not showing',len(exclude),'outlier samples.'
            
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    pylab.xlabel(method.upper()+'-X')
    pylab.ylabel(method.upper()+'-Y')

    axes = getAxesTransposed(scores,exclude=exclude) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)

    marker_size = 15
    if len(column_header)>20:
        marker_size = 12
    if len(column_header)>40:
        marker_size = 10
    if len(column_header)>150:
        marker_size = 7
    if len(column_header)>500:
        marker_size = 5
    if len(column_header)>1000:
        marker_size = 4
    if len(column_header)>2000:
        marker_size = 3
        
    ### Color By Gene
    if colorByGene != None and len(matrix)==0:
        print 'Gene %s not found in the imported dataset... Coloring by groups.' % colorByGene
    if colorByGene != None and len(matrix)>0:
        gene_translation_db={}
        matrix = numpy.array(matrix)
        min_val = matrix.min()  ### min val
        if ' ' in colorByGene:
            genes = string.split(colorByGene,' ')
        else:
            genes = [colorByGene]
        genePresent=False
        numberGenesPresent=[]
        
        for gene in genes:
            if gene in row_header:
                numberGenesPresent.append(gene)
                genePresent = True
        ### Translate symbol to Ensembl
        if len(numberGenesPresent)==0:
            try:
                import gene_associations; from import_scripts import OBO_import
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)

                for symbol in genes:
                    if symbol in symbol_to_gene:
                        gene = symbol_to_gene[symbol][0]
                        if gene in row_header:
                            numberGenesPresent.append(gene)
                            genePresent = True
                            gene_translation_db[symbol]=gene    
            except Exception: pass
            
        numberGenesPresent = len(numberGenesPresent)
        if numberGenesPresent==1:
            cm = pylab.cm.get_cmap('Reds')
        else:
            if numberGenesPresent==2:
                cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                #cm = matplotlib.colors.ListedColormap(['w', 'k']) ### If you want to hide one of the groups
            elif numberGenesPresent==3: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif numberGenesPresent==4:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif numberGenesPresent==5:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==6: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==7:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            else:
                cm = pylab.cm.get_cmap('gist_rainbow')
                
        def get_cmap(n, name='hsv'):
            '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
            RGB color; the keyword argument name must be a standard mpl colormap name.'''
            return pylab.cm.get_cmap(name, n)
        
        
        if genePresent:
            dataset_name+='-'+colorByGene
            group_db={}
            bestGeneAssociated={}
            k=0
            for gene in genes:
                try:
                    try: i = row_header.index(gene)
                    except Exception: i = row_header.index(gene_translation_db[gene])
                    values = map(float,matrix[i])
                    min_val = min(values)
                    bin_size = (max(values)-min_val)/8
                    max_val = max(values)
                    
                    ranges = []
                    iz=min_val
                    while iz < (max(values)-bin_size/100):
                        r = iz,iz+bin_size
                        if len(ranges)==7:
                            r = iz,max_val
                        ranges.append(r)
                        iz+=bin_size
                    color_db = {}
                    colors = get_cmap(len(genes))
                    for i in range(len(ranges)):
                        if i==0:
                            color = '#C0C0C0'
                        else:
                            if numberGenesPresent==1:
                                ### use a single color gradient
                                color = cm(1.*i/len(ranges))
                                #color = cm(1.*(i+1)/len(ranges))
                            else:
                                if i>2:
                                    if len(genes)<8:
                                        color = cm(k)
                                    else:
                                        color = colors(k)
                                else:
                                    color = '#C0C0C0'
                        color_db[ranges[i]] = color
                    i=0
                    for val in values:
                        sample = column_header[i]
                        for (l,u) in color_db:
                            range_index = ranges.index((l,u)) ### what is the ranking of this range
                            if val>=l and val<=u:
                                color = color_db[(l,u)]
                                color_label = [gene+'-range: '+str(l)[:4]+'-'+str(u)[:4],color,'']
                                group_db[sample] = color_label
                                try: bestGeneAssociated[sample].append([range_index,val,color_label])
                                except Exception: bestGeneAssociated[sample] = [[range_index,val,color_label]]
                        i+=1
                    #print min(values),min_val,bin_size,max_val
                    if len(genes)>1:
                        ### Collapse and rank multiple gene results
                        for sample in bestGeneAssociated:
                            bestGeneAssociated[sample].sort()
                            color_label = bestGeneAssociated[sample][-1][-1]
                            if numberGenesPresent>1:
                                index = bestGeneAssociated[sample][-1][0]
                                if index > 2:
                                    gene = string.split(color_label[0],'-')[0]
                                else:
                                    gene = 'Null'
                                color_label[0] = gene
                            group_db[sample] = color_label
                except Exception:
                    print [gene], 'not found in rows...'
                    #print traceback.format_exc()
                k+=1
        else:
            print [colorByGene], 'not found in rows...'
            
    pylab.title(method+' - '+dataset_name)
    
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
        if i not in exclude:
            ax.plot(scores[i][0],scores[i][1],color=color,marker='o',markersize=marker_size,label=label,markeredgewidth=0,picker=True)
            #except Exception: print i, len(scores[pcB]);kill
            if showLabels:
                try: sample_name = '   '+string.split(sample_name,':')[1]
                except Exception: pass
                ax.text(scores[i][0],scores[i][1],sample_name,fontsize=11)
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
    if len(group_count) > 0: ### Make number larger to get the legend in the plot -- BUT, the axis buffer above has been disabled
        # Shink current axis by 20%
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        try: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = Lfontsize) ### move the legend over to the right of the plot
        except Exception: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        ax.set_position([box.x0, box.y0, box.width, box.height])
        pylab.legend(loc="upper left", prop={'size': 10})
  
    filename = 'Clustering-'+dataset_name+'-'+method+'.pdf'
    filename = string.replace(filename,'Clustering-Clustering','Clustering')
    try: pylab.savefig(root_dir + filename)
    except Exception: None ### Rare error
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'

    try: pylab.savefig(root_dir + filename) #dpi=200, transparent=True
    except Exception: None ### Rare error
    graphic_link.append(['Principal Component Analysis',root_dir+filename])

    if display:
        print 'Exporting:',filename
        try:
            pylab.show()
        except Exception:
            #print traceback.format_exc()
            pass### when run in headless mode

def excludeHighlyCorrelatedHits(x,row_header):
    ### For methylation data or other data with redundant signatures, remove these and only report the first one
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
        D1 = numpy.corrcoef(x)
    i=0
    exclude={}
    gene_correlations={}
    include = []
    for score_ls in D1:
        k=0
        for v in score_ls:
            if str(v)!='nan':
                if v>1.00 and k!=i:
                    #print row_header[i], row_header[k], v
                    if row_header[i] not in exclude:
                        exclude[row_header[k]]=[]
                #if k not in exclude: include.append(row_header[k])
            k+=1    
        #if i not in exclude: include.append(row_header[i])
        i+=1
    #print len(exclude),len(row_header);sys.exit()
    return exclude

def PrincipalComponentAnalysis(matrix, column_header, row_header, dataset_name,
            group_db, display=False, showLabels=True, algorithm='SVD', geneSetName=None,
            species=None, pcA=1,pcB=2, colorByGene=None):
    print "Performing Principal Component Analysis..."
    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank

    try: prior_clusters = priorColumnClusters
    except Exception: prior_clusters=[]
    if prior_clusters == None: prior_clusters=[]
    try:
        if len(prior_clusters)>0 and len(group_db)==0:
            newColumnHeader=[]
            i=0
            for sample_name in column_header:
                newColumnHeader.append(str(prior_clusters[i])+':'+sample_name)
                i+=1
            group_db, column_header = assignGroupColors(newColumnHeader)    
    except Exception,e:
        print traceback.format_exc()
        group_db={}
        
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
    #Mdif = matrix-matrix.mean(axis=0)# subtract the mean (along columns)
    #M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    Mdif = matrix/matrix.std()
    Mdif = Mdif.T
    u, s, vt = svd(Mdif, 0)
    fracs = s**2/np.sum(s**2)
    entropy = -sum(fracs*np.log(fracs))/np.log(np.min(vt.shape))
    
    label1 = 'PC%i (%2.1f%%)' %(pcA+1, fracs[0]*100)
    label2 = 'PC%i (%2.1f%%)' %(pcB+1, fracs[1]*100)

    #http://docs.scipy.org/doc/scipy/reference/sparse.html
    #scipy.sparse.linalg.svds - sparse svd
    #idx = numpy.argsort(vt[0,:])
    #print idx;sys.exit() # Use this as your cell order or use a density analysis to get groups
    
    ####  FROM LARSSON ########
    #100 most correlated Genes with PC1
    #print vt
    PCsToInclude = 4
    correlated_db={}
    allGenes={}
    new_matrix = []
    new_headers = []
    added_indexes=[]
    x = 0
    #100 most correlated Genes with PC1
    print 'exporting PCA loading genes to:',root_dir+'/PCA/correlated.txt'
    exportData = export.ExportFile(root_dir+'/PCA/correlated.txt')
    
    matrix = zip(*matrix) ### transpose this back to normal
    try:
        while x<PCsToInclude:
            idx = numpy.argsort(u[:,x])         
            correlated = map(lambda i: row_header[i],idx[:300])
            anticorrelated = map(lambda i: row_header[i],idx[-300:])
            correlated_db[x] = correlated,anticorrelated
            ### Create a new filtered matrix of loading gene indexes
            fidx = list(idx[:300])+list(idx[-300:])
            for i in fidx:
                if i not in added_indexes:
                    added_indexes.append(i)
                    new_headers.append(row_header[i])
                    new_matrix.append(matrix[i])
            x+=1
                    
        #redundant_genes = excludeHighlyCorrelatedHits(numpy.array(new_matrix),new_headers)
        redundant_genes = []
        
        for x in correlated_db:
            correlated,anticorrelated = correlated_db[x]
            count=0
            for gene in correlated:
                if gene not in redundant_genes and count<100:
                    exportData.write(gene+'\tcorrelated-PC'+str(x+1)+'\n'); allGenes[gene]=[]
                    count+=1
            count=0
            for gene in anticorrelated:
                if gene not in redundant_genes and count<100:
                    exportData.write(gene+'\tanticorrelated-PC'+str(x+1)+'\n'); allGenes[gene]=[]
                    count+=1
        exportData.close()
        
        if geneSetName != None:
            if len(geneSetName)>0:
                exportCustomGeneSet(geneSetName,species,allGenes)
                print 'Exported geneset to "StoredGeneSets"'
    except Exception:
        pass       

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
            
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)

    marker_size = 15
    if len(column_header)>20:
        marker_size = 12
    if len(column_header)>40:
        marker_size = 10
    if len(column_header)>150:
        marker_size = 7
    if len(column_header)>500:
        marker_size = 5
    if len(column_header)>1000:
        marker_size = 4
    if len(column_header)>2000:
        marker_size = 3
    #marker_size = 9
        
    #samples = list(column_header)
    ### Color By Gene
    if colorByGene != None:
        gene_translation_db={}
        matrix = numpy.array(matrix)
        min_val = matrix.min()  ### min val
        if ' ' in colorByGene:
            genes = string.split(colorByGene,' ')
        else:
            genes = [colorByGene]
        genePresent=False
        numberGenesPresent=[]
        
        for gene in genes:
            if gene in row_header:
                numberGenesPresent.append(gene)
                genePresent = True
        ### Translate symbol to Ensembl
        if len(numberGenesPresent)==0:
            try:
                import gene_associations; from import_scripts import OBO_import
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
                for symbol in genes:
                    if symbol in symbol_to_gene:
                        gene = symbol_to_gene[symbol][0]
                        if gene in row_header:
                            numberGenesPresent.append(gene)
                            genePresent = True
                            gene_translation_db[symbol]=gene
            except Exception: pass
            
        numberGenesPresent = len(numberGenesPresent)
        if numberGenesPresent==1:
            cm = pylab.cm.get_cmap('Reds')
        else:
            if numberGenesPresent==2:
                cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                #cm = matplotlib.colors.ListedColormap(['w', 'k'])
            elif numberGenesPresent==3: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif numberGenesPresent==4:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif numberGenesPresent==5:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==6: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==7:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            else:
                cmap_c = pylab.cm.gist_rainbow
        if genePresent:
            dataset_name+='-'+colorByGene
            group_db={}
            bestGeneAssociated={}
            k=0
            for gene in genes:
                try:
                    try: i = row_header.index(gene)
                    except Exception: i = row_header.index(gene_translation_db[gene])
                    values = map(float,matrix[i])
                    min_val = min(values)
                    bin_size = (max(values)-min_val)/8
                    max_val = max(values)
                    
                    ranges = []
                    iz=min_val
                    while iz < (max(values)-bin_size/100):
                        r = iz,iz+bin_size
                        if len(ranges)==7:
                            r = iz,max_val
                        ranges.append(r)
                        iz+=bin_size
                    color_db = {}
                    for i in range(len(ranges)):
                        if i==0:
                            color = '#C0C0C0'
                        else:
                            if numberGenesPresent==1:
                                ### use a single color gradient
                                color = cm(1.*i/len(ranges))
                                #color = cm(1.*(i+1)/len(ranges))
                            else:
                                if i>2:
                                    color = cm(k)
                                else:
                                    color = '#C0C0C0'
                        color_db[ranges[i]] = color
                    i=0
                    for val in values:
                        sample = column_header[i]
                        for (l,u) in color_db:
                            range_index = ranges.index((l,u)) ### what is the ranking of this range
                            if val>=l and val<=u:
                                color = color_db[(l,u)]
                                color_label = [gene+'-range: '+str(l)[:4]+'-'+str(u)[:4],color,'']
                                group_db[sample] = color_label
                                try: bestGeneAssociated[sample].append([range_index,val,color_label])
                                except Exception: bestGeneAssociated[sample] = [[range_index,val,color_label]]
                        i+=1
                    #print min(values),min_val,bin_size,max_val
                    if len(genes)>1:
                        ### Collapse and rank multiple gene results
                        for sample in bestGeneAssociated:
                            bestGeneAssociated[sample].sort()
                            color_label = bestGeneAssociated[sample][-1][-1]
                            if numberGenesPresent>1:
                                index = bestGeneAssociated[sample][-1][0]
                                if index > 2:
                                    gene = string.split(color_label[0],'-')[0]
                                else:
                                    gene = 'Null'
                                color_label[0] = gene
                            group_db[sample] = color_label
                except Exception:
                    print [gene], 'not found in rows...'
                    #print traceback.format_exc()
                k+=1
        else:
            print [colorByGene], 'not found in rows...'
            
    pylab.title('Principal Component Analysis - '+dataset_name)
    
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
        try: ax.plot(scores[pcA][i],scores[1][i],color=color,marker='o',markersize=marker_size,label=label,markeredgewidth=0,picker=True)
        except Exception, e: print e; print i, len(scores[pcB]);kill
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
    #group_count = group_count*10 ### force the legend box out of the PCA core plot
    box = ax.get_position()
    if len(group_count) > 0:
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
  
def ViolinPlot():
    def readData(filename):
      all_data = {}
      headers={}
      groups=[]
      firstRow=True
      for line in open(filename,'rU').xreadlines():         
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if firstRow:
              firstRow=False
              i=0
              for x in t[1:]:
                try: g,h = string.split(x,':')
                except Exception: g=x; h=x
                headers[i] = g      
                if g not in groups: groups.append(g)
                i+=1
            else:
              #all_data.append(map(lambda x: math.log(math.pow(2,float(x))-1+0.001,2), t[1:]))
              t = map(lambda x: float(x), t[1:])
              i = 0
              for x in t:
                try: g = headers[i]
                except Exception: print i;sys.exit()
                try: all_data[g].append(x)
                except Exception: all_data[g] = [x]
                i+=1
      all_data2=[]
      print groups
      for group in groups:
        all_data2.append(all_data[group])
      return all_data2
    
    def violin_plot(ax, data, pos, bp=False):
        '''      
        create violin plots on an axis   
        '''
        from scipy.stats import gaussian_kde
        from numpy import arange
        dist = max(pos)-min(pos)
        w = min(0.15*max(dist,1.0),0.5)
        for d,p in zip(data,pos):
            k = gaussian_kde(d) #calculates the kernel density   
            m = k.dataset.min() #lower bound of violin           
            M = k.dataset.max() #upper bound of violin           
            x = arange(m,M,(M-m)/100.) # support for violin      
            v = k.evaluate(x) #violin profile (density curve)    
            v = v/v.max()*w #scaling the violin to the available space       
            ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
            ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
        if bp:
            ax.boxplot(data,notch=1,positions=pos,vert=1)
        
    def draw_all(data, output):
      pos = [1,2,3]
      fig = pylab.figure()
      ax = fig.add_subplot(111)
      violin_plot(ax, data, pos)
      pylab.show()
      pylab.savefig(output+'.pdf')
  
    all_data = []
    all_data = readData('/Users/saljh8/Downloads/TPM_cobound.txt')
    import numpy
    #all_data = map(numpy.array, zip(*all_data))
    #failed_data = map(numpy.array, zip(*failed_data))
    draw_all(all_data, 'alldata')   


def simpleScatter(fn):
    import matplotlib.patches as mpatches
    values=[]
    legends={}
    colors={}
    skip=True
    scale = 100.0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if skip:
            x_header, y_header, color_header,label_header, shape_header = string.split(data,'\t')
            skip=False
        else:
            x, y, color,label,shape = string.split(data,'\t')
            if color in colors:
                xval,yval,label,shape = colors[color]
                xval.append(float(x)); yval.append(float(y))
            else:
                xval = [float(x)]; yval = [float(y)]
                colors[color] = xval,yval,label,shape
    
    for color in colors:
        xval,yval,label,shape = colors[color]
        pylab.scatter(xval, yval, s=scale, c=color, alpha=0.75, label=label, marker=shape,edgecolor="none")

    pylab.legend(loc='upper left')
    
    pylab.title(fn)
    pylab.xlabel(x_header, fontsize=15)
    pylab.ylabel(y_header, fontsize=15)
    marker_size = 7
    
    #pylab.grid(True)
    pylab.show()
    
def ica(filename):
    showLabels=True
    X, column_header, row_header, dataset_name, group_db = importData(filename)
    X = map(numpy.array, zip(*X)) ### coverts these to tuples
    column_header, row_header = row_header, column_header

    ica = FastICA()
    scores = ica.fit(X).transform(X)  # Estimate the sources
    
    scores /= scores.std(axis=0)
        
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    pylab.xlabel('ICA-X')
    pylab.ylabel('ICA-Y')
    pylab.title('ICA - '+dataset_name)
            
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)
    
    marker_size = 15
    if len(column_header)>20:
        marker_size = 12
    if len(column_header)>40:
        marker_size = 10
    if len(column_header)>150:
        marker_size = 7
    if len(column_header)>500:
        marker_size = 5
    if len(column_header)>1000:
        marker_size = 4
    if len(column_header)>2000:
        marker_size = 3
        
    group_names={}
    i=0
    for sample_name in row_header: #scores[0]
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
        ax.plot(scores[0][i],scores[1][i],color=color,marker='o',markersize=marker_size,label=label)
        if showLabels:
            ax.text(scores[0][i],scores[1][i],sample_name,fontsize=8)
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
    
def PCA3D(matrix, column_header, row_header, dataset_name, group_db,
          display=False, showLabels=True, algorithm='SVD',geneSetName=None,
          species=None,colorByGene=None):

    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    start = time.time()
    #M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    
    try: prior_clusters = priorColumnClusters
    except Exception: prior_clusters=[]
    if prior_clusters == None: prior_clusters=[]
    try:
        if len(prior_clusters)>0 and len(group_db)==0:
            newColumnHeader=[]
            i=0
            for sample_name in column_header:
                newColumnHeader.append(str(prior_clusters[i])+':'+sample_name)
                i+=1
            group_db, column_header = assignGroupColors(newColumnHeader)    
    except Exception,e:
        #print e
        group_db={}
        
    if algorithm == 'SVD': use_svd = True
    else: use_svd = False
    Mdif = matrix/matrix.std()
    Mdif = Mdif.T
    u, s, vt = svd(Mdif, 0)

    fracs = s**2/np.sum(s**2)
    entropy = -sum(fracs*np.log(fracs))/np.log(np.min(vt.shape))
    
    label1 = 'PC%i (%2.1f%%)' %(0+1, fracs[0]*100)
    label2 = 'PC%i (%2.1f%%)' %(1+1, fracs[1]*100)
    label3 = 'PC%i (%2.1f%%)' %(2+1, fracs[2]*100)

    PCsToInclude = 4
    correlated_db={}
    allGenes={}
    new_matrix = []
    new_headers = []
    added_indexes=[]
    x = 0
    #100 most correlated Genes with PC1
    print 'exporting PCA loading genes to:',root_dir+'/PCA/correlated.txt'
    exportData = export.ExportFile(root_dir+'/PCA/correlated.txt')
    
    matrix = zip(*matrix) ### transpose this back to normal
    try:
        while x<PCsToInclude:
            idx = numpy.argsort(u[:,x])         
            correlated = map(lambda i: row_header[i],idx[:300])
            anticorrelated = map(lambda i: row_header[i],idx[-300:])
            correlated_db[x] = correlated,anticorrelated
            ### Create a new filtered matrix of loading gene indexes
            fidx = list(idx[:300])+list(idx[-300:])
            for i in fidx:
                if i not in added_indexes:
                    added_indexes.append(i)
                    new_headers.append(row_header[i])
                    new_matrix.append(matrix[i])
            x+=1
                    
        #redundant_genes = excludeHighlyCorrelatedHits(numpy.array(new_matrix),new_headers)
        redundant_genes = []
        
        for x in correlated_db:
            correlated,anticorrelated = correlated_db[x]
            count=0
            for gene in correlated:
                if gene not in redundant_genes and count<100:
                    exportData.write(gene+'\tcorrelated-PC'+str(x+1)+'\n'); allGenes[gene]=[]
                    count+=1
            count=0
            for gene in anticorrelated:
                if gene not in redundant_genes and count<100:
                    exportData.write(gene+'\tanticorrelated-PC'+str(x+1)+'\n'); allGenes[gene]=[]
                    count+=1
        exportData.close()
        
        if geneSetName != None:
            if len(geneSetName)>0:
                exportCustomGeneSet(geneSetName,species,allGenes)
                print 'Exported geneset to "StoredGeneSets"'
    except ZeroDivisionError:
        pass       

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
    axes = getAxes(scores,PlotType='3D') ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)
    
    Lfontsize = 8
    group_count = []
    for i in group_db:
        if group_db[i][0] not in group_count:
            group_count.append(group_db[i][0])
    
    ### Color By Gene
    if colorByGene != None:
        gene_translation_db={}
        matrix = numpy.array(matrix)
        min_val = matrix.min()  ### min val
        if ' ' in colorByGene:
            genes = string.split(colorByGene,' ')
        else:
            genes = [colorByGene]
        genePresent=False
        numberGenesPresent=[]
        
        for gene in genes:
            if gene in row_header:
                numberGenesPresent.append(gene)
                genePresent = True
        ### Translate symbol to Ensembl
        if len(numberGenesPresent)==0:
            try:
                import gene_associations; from import_scripts import OBO_import
                gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
                for symbol in genes:
                    if symbol in symbol_to_gene:
                        gene = symbol_to_gene[symbol][0]
                        if gene in row_header:
                            numberGenesPresent.append(gene)
                            genePresent = True
                            gene_translation_db[symbol]=gene
            except Exception: pass
            
        numberGenesPresent = len(numberGenesPresent)
        if numberGenesPresent==1:
            cm = pylab.cm.get_cmap('Reds')
        else:
            if numberGenesPresent==2:
                cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                cm = matplotlib.colors.ListedColormap(['w', 'k'])
            elif numberGenesPresent==3: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif numberGenesPresent==4:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif numberGenesPresent==5:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==6: 
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            elif numberGenesPresent==7:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976','#FEBC18', '#EE2C3C'])
            else:
                cm = pylab.cm.get_cmap('gist_rainbow')
        if genePresent:
            dataset_name+='-'+colorByGene
            group_db={}
            bestGeneAssociated={}
            k=0
            for gene in genes:
                try:
                    try: i = row_header.index(gene)
                    except Exception: i = row_header.index(gene_translation_db[gene])
                    values = map(float,matrix[i])
                    min_val = min(values)
                    bin_size = (max(values)-min_val)/8
                    max_val = max(values)
                    
                    ranges = []
                    iz=min_val
                    while iz < (max(values)-bin_size/100):
                        r = iz,iz+bin_size
                        if len(ranges)==7:
                            r = iz,max_val
                        ranges.append(r)
                        iz+=bin_size
                    color_db = {}
                    for i in range(len(ranges)):
                        if i==0:
                            color = '#C0C0C0'
                        else:
                            if numberGenesPresent==1:
                                ### use a single color gradient
                                color = cm(1.*i/len(ranges))
                                #color = cm(1.*(i+1)/len(ranges))
                            else:
                                if i>2:
                                    color = cm(k)
                                else:
                                    color = '#C0C0C0'
                        color_db[ranges[i]] = color
                    i=0
                    for val in values:
                        sample = column_header[i]
                        for (l,u) in color_db:
                            range_index = ranges.index((l,u)) ### what is the ranking of this range
                            if val>=l and val<=u:
                                color = color_db[(l,u)]
                                color_label = [gene+'-range: '+str(l)[:4]+'-'+str(u)[:4],color,'']
                                group_db[sample] = color_label
                                try: bestGeneAssociated[sample].append([range_index,val,color_label])
                                except Exception: bestGeneAssociated[sample] = [[range_index,val,color_label]]
                        i+=1
                    #print min(values),min_val,bin_size,max_val
                    if len(genes)>1:
                        ### Collapse and rank multiple gene results
                        for sample in bestGeneAssociated:
                            bestGeneAssociated[sample].sort()
                            color_label = bestGeneAssociated[sample][-1][-1]
                            if numberGenesPresent>1:
                                index = bestGeneAssociated[sample][-1][0]
                                if index > 2:
                                    gene = string.split(color_label[0],'-')[0]
                                else:
                                    gene = 'Null'
                                color_label[0] = gene
                            group_db[sample] = color_label
                except Exception:
                    print [gene], 'not found in rows...'
                    #print traceback.format_exc()
                k+=1
        else:
            print [colorByGene], 'not found in rows...'
            
    #print len(group_count)
    if len(group_count)>20:
        Lfontsize = 10
    if len(group_count)>30:
        Lfontsize = 8
    if len(group_count)>40:
        Lfontsize = 6
    if len(group_count)>50:
        Lfontsize = 5
    
    if len(scores[0])>150:
        markersize = 7
    else:
        markersize = 10
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

        ax.plot([scores[0][i]],[scores[1][i]],[scores[2][i]],color=color,marker='o',markersize=markersize,label=label,markeredgewidth=0,picker=True) #markeredgecolor=color
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
    
def getAxes1(scores,PlotType=None):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    try:
        x_range = max(scores[0])-min(scores[0])
        y_range = max(scores[1])-min(scores[1])
        if PlotType == '3D':
            x_axis_min = min(scores[0])-(x_range/10)
            x_axis_max = max(scores[0])+(x_range/10)
            y_axis_min = min(scores[1])-(y_range/10)
            y_axis_max = max(scores[1])+(y_range/10) 
        else:
            x_axis_min = min(scores[0])-(x_range/10)
            x_axis_max = max(scores[0])+(x_range/10)
            y_axis_min = min(scores[1])-(y_range/10)
            y_axis_max = max(scores[1])+(y_range/10) 
    except KeyError:
        None
    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]


def getAxes(scores,PlotType=None):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    try:
        x_range = max(scores[0])-min(scores[0])
        y_range = max(scores[1])-min(scores[1])
        if PlotType == '3D':
            x_axis_min = min(scores[0])-(x_range/1.5)
            x_axis_max = max(scores[0])+(x_range/1.5)
            y_axis_min = min(scores[1])-(y_range/5)
            y_axis_max = max(scores[1])+(y_range/5)
        else:
            x_axis_min = min(scores[0])-(x_range/10)
            x_axis_max = max(scores[0])+(x_range/10)
            y_axis_min = min(scores[1])-(y_range/10)
            y_axis_max = max(scores[1])+(y_range/10) 
    except KeyError:
        None
    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]

def getAxesTransposed(scores,exclude={}):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    scores_filtered=[]

    for i in range(len(scores)):
        if i not in exclude:
            scores_filtered.append(scores[i])
    scores = scores_filtered
    scores = map(numpy.array, zip(*scores))

    try:
        x_range = max(scores[0])-min(scores[0])
        y_range = max(scores[1])-min(scores[1])
        x_axis_min = min(scores[0])-int((float(x_range)/7))
        x_axis_max = max(scores[0])+int((float(x_range)/7))
        y_axis_min = min(scores[1])-int(float(y_range/7))
        y_axis_max = max(scores[1])+int(float(y_range/7))
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
        if allowLargeClusters: maxSize = 50000
        else: maxSize = 7000
    except Exception: maxSize = 7000
    
    try:
        PriorColumnClusters=priorColumnClusters
        PriorRowClusters=priorRowClusters
    except Exception:
        PriorColumnClusters=None
        PriorRowClusters=None
        
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
                        contrast=contrast,allowAxisCompression=allowAxisCompression,Normalize=Normalize,
                        PriorColumnClusters=PriorColumnClusters,PriorRowClusters=PriorRowClusters)
                run = True
            except Exception:
                print traceback.format_exc()
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
                            contrast=contrast,allowAxisCompression=allowAxisCompression,Normalize=Normalize,
                            PriorColumnClusters=PriorColumnClusters,PriorRowClusters=PriorRowClusters)
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
    """ Explicit method for hieararchical clustering with defaults defined by the user (see below function) """
    #print [filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient, contrast, Normalize]
    
    global root_dir
    global inputFilename
    global originalFilename
    global graphic_link
    global allowLargeClusters
    global GroupDB
    global justShowTheseIDs
    global targetGeneIDs
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
        except Exception:
            rho_cutoff = 0.5 ### Always done if no rho, but only used if getGeneCorrelations == True
            #print 'Setting correlation cutoff to a rho of',rho_cutoff
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
        except Exception:
            #print traceback.format_exc()
            pass
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
    if len(justShowTheseIDs)==0:
        try:
            if len(priorColumnClusters)>0 and priorRowClusters>0 and row_method==None and column_method == None:
                try: justShowTheseIDs = importPriorDrivers(inputFilename)
                except Exception: pass #justShowTheseIDs=[]
        except Exception:
            #print traceback.format_exc()
            pass

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
            alt_targetGene = string.replace(alt_targetGene,'guide','')
            alt_targetGene = string.replace(alt_targetGene,'top','')
            alt_targetGene = string.replace(alt_targetGene,'positive','')
            alt_targetGene = string.replace(alt_targetGene,'excludeCellCycle','')
            alt_targetGene = string.replace(alt_targetGene,'monocle','')
            alt_targetGene = string.replace(alt_targetGene,'GuideOnlyCorrelation','')
            alt_targetGene = string.replace(alt_targetGene,' ','')  
        except Exception:
            alt_targetGene = ''
        if getGeneCorrelations and targetGene != 'driver' and targetGene != 'GuideOnlyCorrelation' and \
                    targetGene != 'guide' and targetGene !='excludeCellCycle' and \
                    targetGene !='top' and targetGene != ' monocle' and \
                    targetGene !='positive' and len(alt_targetGene)>0: ###Restrict analyses to only genes that correlate with the target gene of interest
            allowAxisCompression = False
            if transpose and transpose_update == False: transpose_update = False ### If filterByPathways selected
            elif transpose and transpose_update: transpose_update = True ### If filterByPathways not selected
            else: transpose_update = False ### If transpose == False
            if '\r' in targetGene or '\n' in targetGene:
                targetGene = string.replace(targetGene, '\r',' ')
                targetGene = string.replace(targetGene, '\n',' ')
            if len(targetGene)>15:
                inputFilename = string.replace(newInput,'.txt','-'+targetGene[:50]+'.txt') ### update the pathway reference for HOPACH
                dataset_name += '-'+targetGene[:50]
            else:
                inputFilename = string.replace(newInput,'.txt','-'+targetGene+'.txt') ### update the pathway reference for HOPACH
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
            targetGeneIDs = targetGene
            exportTargetGeneList(targetGene,inputFilename)      
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

    if 'guide' in targetGene:
        import RNASeq
        input_file = graphic_link[-1][-1][:-4]+'.txt'
        if 'excludeCellCycle' in targetGene: excludeCellCycle = True
        else: excludeCellCycle = False
        print 'excludeCellCycle',excludeCellCycle
        targetGene = RNASeq.remoteGetDriverGenes(species,platform,input_file,excludeCellCycle=excludeCellCycle,ColumnMethod=column_method)
        extra_params.setGeneSelection(targetGene) ### force correlation to these
        extra_params.setGeneSet('None Selected') ### silence this
        graphic_link= runHCexplicit(filename, graphic_link, row_method, row_metric, column_method, column_metric, color_gradient,
                extra_params, display=display, contrast=contrast, Normalize=Normalize, JustShowTheseIDs=JustShowTheseIDs,compressAxis=compressAxis)
    return graphic_link

def importPriorDrivers(inputFilename):
    filename = string.replace(inputFilename,'Clustering-','')
    filename = string.split(filename,'-hierarchical')[0]+'-targetGenes.txt'
    genes = open(filename, "rU")
    genes = map(lambda x: cleanUpLine(x),genes)
    return genes

def exportTargetGeneList(targetGene,inputFilename):
    exclude=['positive','top','driver', 'guide', 'amplify','GuideOnlyCorrelation']
    exportFile = inputFilename[:-4]+'-targetGenes.txt'
    eo = export.ExportFile(root_dir+findFilename(exportFile))
    targetGenes = string.split(targetGene,' ')
    for gene in targetGenes:
        if gene not in exclude:
            try: eo.write(gene+'\n')
            except Exception: print 'Error export out gene (bad ascii):', [gene]
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
    from import_scripts import OBO_import

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
                try:
                    cluster,row_id = string.split(row_id,':')
                    updated_row_id = cluster+':'+symbol
                except Exception:
                    pass
            else:
                updated_row_id = symbol
            try: original_id = updated_row_id
            except Exception: pass
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
        #from import_scripts import OBO_import; symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    except Exception:
        print 'No Ensembl-Symbol database available for',species
    
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
    intersecting_ids=[]
    i=0
        
    ### If multiple genes entered, just display these
    targetGenes=[targetGene]
    if ' ' in targetGene or ',' in targetGene or '|' in targetGene or '\n' in targetGene or '\r' in targetGene:
        multipleGenes = True
        if '  ' in targetGene: targetGene = string.replace(targetGene,'  ', ' ')
        if ',' in targetGene: targetGene = string.replace(targetGene,',', ' ')
        #if '|' in targetGene and 'alt_junction' not in originalFilename: targetGene = string.replace(targetGene,'|', ' ')
        if '\n' in targetGene: targetGene = string.replace(targetGene,'\n', ' ')
        if '\r' in targetGene: targetGene = string.replace(targetGene,'\r', ' ')

        targetGenes = string.split(targetGene,' ')

        if row_method != None: targetGenes.sort()
        
        intersecting_ids = [val for val in targetGenes if val in row_header]
        
        for row_id in row_header:
            original_rowid = row_id
            symbol=row_id
            new_symbol = symbol
            rigorous_search = True
            if ':' in row_id and '|' in row_id:
                rigorous_search = False
            elif ':' in row_id and '|' not in row_id:
                a,b = string.split(row_id,':')[:2]
                if 'ENS' in a or len(a)==17:
                    try:
                        row_id = a
                        symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        symbol =''
                elif 'ENS' not in b and len(a)!=17:
                    row_id = b
                elif 'ENS' in b:
                    symbol = original_rowid
                    row_id = a
            if rigorous_search:
                try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
                except Exception:
                    try: symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        if 'ENS' not in original_rowid:
                            row_id, symbol = row_id, row_id
                    new_symbol = symbol
                if 'ENS' not in original_rowid and len(original_rowid)!=17:
                    if original_rowid != symbol:
                        symbol = original_rowid+' '+symbol
                for gene in targetGenes:
                    if string.lower(gene) == string.lower(row_id) or string.lower(gene) == string.lower(symbol) or string.lower(original_rowid)==string.lower(gene) or string.lower(gene) == string.lower(new_symbol):
                        matrix2.append(matrix[i]) ### Values for the row
                        row_header2.append(symbol)
                        matrix_db[symbol]=matrix[i]
            else:
                if row_id in targetGenes:
                    matrix2.append(matrix[i])
                    row_header2.append(row_id)
                    matrix_db[row_id]=matrix[i]
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

    #print len(intersecting_ids),len(targetGenes), multipleGenes
    if multipleGenes==False or 'amplify' in targetGene or 'correlated' in targetGene:
        row_header3=[] ### Convert to symbol if possible
        if multipleGenes==False:
            targetGeneValue_array = [targetGeneValues]
        else:
            targetGeneValue_array = matrix2
            if (len(row_header2)>4) or len(row_header)>15000: # and len(row_header)<20000
                print 'Performing all iterative pairwise corelations...',
                corr_matrix = numpyCorrelationMatrixGene(matrix,row_header,row_header2,gene_to_symbol)
                print 'complete'
            matrix2=[]; original_headers=row_header2; row_header2 = []
            matrix2_alt=[]; row_header2_alt=[]
        ### If one gene entered, display the most positive and negative correlated
        import markerFinder; k=0
        for targetGeneValues in targetGeneValue_array:
            correlated=[]
            anticorrelated=[]
            try: targetGeneID = original_headers[k]
            except Exception: targetGeneID=''
            try:
                rho_results = list(corr_matrix[targetGeneID])
            except Exception:
                #print traceback.format_exc()
                rho_results = markerFinder.simpleScipyPearson(matrix,targetGeneValues)
            correlated_symbols={}
            #print targetGeneID, rho_results[:100]
            #print targetGeneID, rho_results[-100:];sys.exit()
            for (rho,ind) in rho_results[:limit]: ### Get the top-50 correlated plus the gene of interest
                proceed = True
                try:
                    if len(rho)==2: rho = rho[0]
                except: pass
                if 'top' in targetGene:
                    if rho_results[4][0]<rho_cutoff: proceed = False
                if rho>rho_cutoff and proceed: #and rho_results[3][0]>rho_cutoff:# ensures only clustered genes considered
                    rh = row_header[ind]
                    #if gene_to_symbol[rh][0] in targetGenes:correlated.append(gene_to_symbol[rh][0])
                    #correlated.append(gene_to_symbol[rh][0])
                    
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
                try:
                    if len(rho)==2: rho = rho[0]
                except: pass
                if rho<-1*rho_cutoff and 'positive' not in targetGene:
                    rh = row_header[ind]
                    #if gene_to_symbol[rh][0] in targetGenes:anticorrelated.append(gene_to_symbol[rh][0])
                    #anticorrelated.append(gene_to_symbol[rh][0])
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
            #print targetGeneID+'\t'+str(len(correlated))+'\t'+str(len(anticorrelated))
        #sys.exit()
        
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
    if 'amplify' not in targetGene and 'correlated' not in targetGene:
        ### reorder according to orignal
        matrix_temp=[]
        header_temp=[]
        #print targetGenes
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

    exclude=[]
    #exclude = excludeHighlyCorrelatedHits(numpy.array(matrix2),row_header2)
    exportData.write(string.join(['UID']+column_header,'\t')+'\n') ### title row export
    i=0

    for row_id in row_header2:
        if ':' in row_id and '|' not in row_id:
            a,b = string.split(row_id,':')[:2]
            if 'ENS' in a:
                try: row_id=string.replace(row_id,a,gene_to_symbol[a][0])
                except Exception,e: pass
                row_header2[i] = row_id
        elif 'ENS' in row_id and ' ' in row_id and '|' not in row_id:
            row_id = string.split(row_id, ' ')[1]
            row_header2[i] = row_id
        elif ' ' in row_id:
            try: a,b = string.split(row_id, ' ')
            except Exception: a = 1; b=2
            if a==b:
                row_id = a
        if row_id not in exclude:
            exportData.write(string.join([row_id]+map(str,matrix2[i]),'\t')+'\n') ### export values
        i+=1
    if 'amplify' not in targetGene and 'correlated' not in targetGene:
        print len(row_header2), 'input gene IDs found'
    else:
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
        #D1 = numpy.ma.corrcoef(x)
        D1 = numpy.corrcoef(x)
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
        #D1 = numpy.ma.corrcoef(x)
        D1 = numpy.corrcoef(x)
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
    """ Simple method for hieararchical clustering with defaults defined by the function rather than the user (see above function) """
    
    global root_dir
    global graphic_link
    global inputFilename
    global GroupDB
    global allowLargeClusters
    global runGOElite
    global EliteGeneSets
    runGOElite = False
    EliteGeneSets=[]
    allowLargeClusters = False
    
    ###############
    global inputFilename
    global originalFilename
    global GroupDB
    global justShowTheseIDs
    global targetGeneIDs
    global normalize
    global species
    global storeGeneSetName
    targetGene=[]
    filterByPathways=False
    justShowTheseIDs=[]
    ###############
    
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

def timestamp():
    import datetime
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
    time_stamp = string.replace(time.ctime(),':','')
    time_stamp = string.replace(time_stamp,'  ',' ')
    time_stamp = string.split(time_stamp,' ') ###Use a time-stamp as the output dir (minus the day)
    time_stamp = today+'-'+time_stamp[3]
    return time_stamp

def runPCAonly(filename,graphics,transpose,showLabels=True,plotType='3D',display=True,
               algorithm='SVD',geneSetName=None, species=None, zscore=True, colorByGene=None,
               reimportModelScores=True, separateGenePlots=False):

    global root_dir
    global graphic_link
    graphic_link=graphics ### Store all locations of pngs
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'/DataPlots','')
    root_dir = string.replace(root_dir,'/amplify','')
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
    root_dir = string.replace(root_dir,'ExpressionInput','DataPlots')
    root_dir = string.replace(root_dir,'ICGS-NMF','DataPlots')
    if 'DataPlots' not in root_dir:
        root_dir += '/DataPlots/'
    try: os.mkdir(root_dir) ### May need to create this directory
    except Exception: None
        
    ### Transpose matrix and build PCA
    geneFilter=None
    if (algorithm == 't-SNE' or algorithm == 'UMAP') and reimportModelScores:
        dataset_name = string.split(filename,'/')[-1][:-4]
        try:
            ### if the scores are present, we only need to import the genes of interest (save time importing large matrices)
            if algorithm == 't-SNE':
                importtSNEScores(root_dir+dataset_name+'-t-SNE_scores.txt')
            if algorithm == 'UMAP':
                importtSNEScores(root_dir+dataset_name+'-UMAP_scores.txt')
            if len(colorByGene)==None:
                geneFilter = [''] ### It won't import the matrix, basically
            elif ' ' in colorByGene or ',' in colorByGene:
                colorByGene = string.replace(colorByGene,',',' ')
                geneFilter = string.split(colorByGene,' ')
            else:
                geneFilter = [colorByGene]
        except Exception:
            #print traceback.format_exc();sys.exit()
            geneFilter = None ### It won't import the matrix, basically

    matrix, column_header, row_header, dataset_name, group_db = importData(filename,zscore=zscore,geneFilter=geneFilter)
    if transpose == False: ### We normally transpose the data, so if True, we don't transpose (I know, it's confusing)
        matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
        column_header, row_header = row_header, column_header
    if (len(column_header)>1000 or len(row_header)>1000) and algorithm != 't-SNE' and algorithm != 'UMAP':
        print 'Performing Principal Component Analysis (please be patient)...'
    #PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        if algorithm == 't-SNE' or algorithm == 'UMAP':
            matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
            column_header, row_header = row_header, column_header
            if separateGenePlots and (len(colorByGene)>0 or colorByGene==None):
                for gene in geneFilter:
                    tSNE(numpy.array(matrix),column_header,dataset_name,group_db,display=False,
                         showLabels=showLabels,row_header=row_header,colorByGene=gene,species=species,
                         reimportModelScores=reimportModelScores,method=algorithm)
                if display:
                    ### Show the last one
                    tSNE(numpy.array(matrix),column_header,dataset_name,group_db,display=True,
                        showLabels=showLabels,row_header=row_header,colorByGene=gene,species=species,
                        reimportModelScores=reimportModelScores,method=algorithm)    
            else:
                tSNE(numpy.array(matrix),column_header,dataset_name,group_db,display=display,
                     showLabels=showLabels,row_header=row_header,colorByGene=colorByGene,species=species,
                     reimportModelScores=reimportModelScores,method=algorithm)
        elif plotType == '3D':
            try: PCA3D(numpy.array(matrix), row_header, column_header, dataset_name, group_db,
                       display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName,
                       species=species, colorByGene=colorByGene)
            except Exception:
                print traceback.format_exc()
                PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header,
                        dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm,
                        geneSetName=geneSetName, species=species, colorByGene=colorByGene)
        else:
            PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name,
                        group_db, display=display, showLabels=showLabels, algorithm=algorithm,
                        geneSetName=geneSetName, species=species, colorByGene=colorByGene)

    return graphic_link

def outputClusters(filenames,graphics,Normalize=False,Species=None,platform=None,vendor=None):
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
    if len(row_header)<700000 and len(column_header)<700000 and len(column_header)>2:
        PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db)
    else:
        print 'SKIPPING PCA!!! - Your dataset file is over or under the recommended size limit for clustering (>7000 rows). Please cluster later using "Additional Analyses".'

    row_method = 'average'
    column_method = 'average'
    row_metric = 'cosine'
    column_metric = 'cosine'
    color_gradient = 'red_white_blue'
    color_gradient = 'red_black_sky'
    
    global species
    species = Species
    if 'LineageCorrelations' not in filename and 'Zscores' not in filename:
        EliteGeneSets=['GeneOntology']
        runGOElite = True

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
    elif len(edges) > 3000:
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
    from import_scripts import OBO_import
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

def stackedbarchart(filename,display=False,output=False):
    header=[]
    conditions = []
    data_matrix=[]
    for line in open(filename,'rU').xreadlines():         
        cd = cleanUpLine(line)
        t = string.split(cd,'\t')
        if len(header)==0:
            header = t[4:]
            exc_indexes = [0,2,4,6,8,10,12]
            inc_indexes = [1,3,5,7,9,11,13]
            inlc_header = map(lambda i: string.split(header[i],'_')[0],inc_indexes)
            header = inlc_header            
        else:
            condition = t[0]
            data = t[4:]
            conditions.append(condition+'-inclusion ')
            data_matrix.append(map(lambda i: float(data[i]),inc_indexes))
            conditions.append(condition+'-exclusion ')
            data_matrix.append(map(lambda i: float(data[i]),exc_indexes))
            
    data_matrix = map(numpy.array, zip(*data_matrix))
    
    #https://www.w3resource.com/graphics/matplotlib/barchart/matplotlib-barchart-exercise-16.php
    # multi-dimensional data_matrix 
    y_pos = np.arange(len(conditions))
    
    fig, ax = pylab.subplots()
    #fig = pylab.figure(figsize=(10,8))
    #ax = fig.add_subplot(111)
    #pos1 = ax.get_position() # get the original position 
    #pos2 = [pos1.x0 + 0.2, pos1.y0 - 0.2,  pos1.width / 1.2, pos1.height / 1.2 ] 
    #ax.set_position(pos2) # set a new position
    
    colors =['royalblue','salmon','grey','gold','cornflowerblue','mediumseagreen','navy']
    patch_handles = []
    # left alignment of data_matrix starts at zero
    left = np.zeros(len(conditions))
    index=0
    for i, d in enumerate(data_matrix):
        patch_handles.append(ax.barh(y_pos, d, 0.3,
            color=colors[index], align='center', 
            left=left,label = header[index]))
        left += d
        index+=1
    
    # search all of the bar segments and annotate
    """
    for j in range(len(patch_handles)):
        for i, patch in enumerate(patch_handles[j].get_children()):
            bl = patch.get_xy()
            x = 0.5*patch.get_width() + bl[0]
            y = 0.5*patch.get_height() + bl[1]
            #ax.text(x,y, "%d%%" % (percentages[i,j]), ha='center')
    """
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(conditions)
    ax.set_xlabel('Events')
    ax.legend(loc="best", bbox_to_anchor=(1.0, 1.0))

    box = ax.get_position()

    # Shink current axis by 20%
    ax.set_position([box.x0+0.2, box.y0, box.width * 0.6, box.height])
    try: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 10) ### move the legend over to the right of the plot
    except Exception: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title('MultiPath-PSI Splicing Event Types')
    
    #pylab.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    if output==False:
        pylab.savefig(filename[:-4]+'.pdf')
        pylab.savefig(filename[:-4]+'.png')
    else:
        pylab.savefig(output[:-4]+'.pdf')
        pylab.savefig(output[:-4]+'.png')
    
    if display:
        print 'Exporting:',filename
        try: pylab.show()
        except Exception: None ### when run in headless mode    
    
    
def barchart(filename,index1,index2,x_axis,y_axis,title,display=False,color1='SkyBlue',color2='IndianRed',output=False):

    header=[]
    reference_data=[]
    query_data=[]
    groups=[]
    for line in open(filename,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if len(header)==0:
            header = t
            header1=header[index1]
            header2=header[index2]
        else:
            reference_data.append(int(t[index1]))
            query_data.append(int(t[index2]))
            name = t[0]
            if '_vs_' in name and 'event_summary' not in filename:
                name = string.split(name,'_vs_')[0]
                if '_' in name:
                    name = string.split(name,'_')[:-1]
                    name = string.join(name,'_')
            groups.append(name)

    fig, ax = pylab.subplots()
    
    pos1 = ax.get_position() # get the original position 
    pos2 = [pos1.x0 + 0.1, pos1.y0 + 0.1,  pos1.width / 1.2, pos1.height / 1.2 ] 
    ax.set_position(pos2) # set a new position
    
    ind = np.arange(len(groups))  # the x locations for the groups
    width = 0.35  # the width of the bars
    query_data.reverse()
    reference_data.reverse()
    groups.reverse()

    ax.barh(ind - width/2, query_data, width, color=color2, label=header2)
    ax.barh(ind + width/2, reference_data, width,color=color1, label=header1)

    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    ax.set_yticks(ind+0.175)
    ax.set_yticklabels(groups)
    ax.set_title(title)
    ax.legend()

    if output==False:
        pylab.savefig(filename[:-4]+'.pdf')
        pylab.savefig(filename[:-4]+'.png')
    else:
        pylab.savefig(output[:-4]+'.pdf')
        pylab.savefig(output[:-4]+'.png')
        
    if display:
        print 'Exporting:',filename
        try: pylab.show()
        except Exception: None ### when run in headless mode    
    
def multipleSubPlots(filename,uids,SubPlotType='column',n=20):
    #uids = [uids[-1]]+uids[:-1]
    str_uids = string.join(uids,'_')
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,geneFilter=uids)
    for uid in uids:
        if uid not in row_header:
            print uid,"is missing from the expression file."
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
    #print row_header
    color_list = ['r', 'b', 'y', 'g', 'w', 'k', 'm']
    
    groups=[]
    for sample in column_header:
        try: group = group_db[sample][0]
        except: group = '1'
        if group not in groups:
            groups.append(group)
            
    fontsize=10
    if len(groups)>0:
        color_list = []
        if len(groups)==9:
            cm = matplotlib.colors.ListedColormap(['#80C241', '#118943', '#6FC8BB', '#ED1D30', '#F26E21','#8051A0', '#4684C5', '#FBD019','#3A52A4'])
        elif len(groups)==3:
            cm = matplotlib.colors.ListedColormap(['#4684C4','#FAD01C','#7D7D7F'])
        #elif len(groups)==5: cm = matplotlib.colors.ListedColormap(['#41449B','#6182C1','#9DDAEA','#42AED0','#7F7F7F'])
        else:
            cm = pylab.cm.get_cmap('gist_rainbow') #gist_ncar
        for i in range(len(groups)):
            color_list.append(cm(1.*i/len(groups)))  # color will now be an RGBA tuple
            
    for i in range(len(matrix)):
        ax = pylab.subplot(n,1,1+i)
        OY = matrix[i]
        pylab.xlim(0,len(OY))
        pylab.subplots_adjust(right=0.85)
        ind = np.arange(len(OY))
        index_list = []
        v_list = []
        colors_list = []
        if SubPlotType=='column':
            index=-1
            for v in OY:
                index+=1
                try: group = group_db[column_header[index]][0]
                except: group = '1'
                index_list.append(index)
                v_list.append(v)
                colors_list.append(color_list[groups.index(group)])
                #pylab.bar(index, v,edgecolor='black',linewidth=0,color=color_list[groups.index(group)])
                width = .35
            #print i ,row_header[i]
            print 1
            barlist = pylab.bar(index_list, v_list,edgecolor='black',linewidth=0)
            ci = 0

            for cs in barlist:
                barlist[ci].set_color(colors_list[ci])
                ci+=1
        if SubPlotType=='plot':
            pylab.plot(x,y)

        ax.text(matrix.shape[1]-0.5, i, '  '+row_header[i],fontsize=8)
        
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
        
    #pylab.show()
    if len(str_uids)>50:
        str_uids = str_uids[:50]
    pylab.savefig(filename[:-4]+'-1'+str_uids+'.pdf')

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

def BedFileCheck(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.bed'
    ea = export.ExportFile(filename)
    
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
        else:
            #if len(t) != 12: print len(t);sys.exit()
            ea.write(string.join(t,'\t')+'\n')
    ea.close()
    
def simpleFilter(filename):
    fn = filepath(filename)
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)

    matrix = []
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,',')

        uid = t[0]
        #if '=chr' in t[0]:
        if 1==2:
            a,b = string.split(t[0],'=')
            b = string.replace(b,'_',':')
            uid = a+ '='+b
            matrix.append(t)
        ea.write(string.join([uid]+t[1:],'\t')+'\n')
    ea.close()
    
def test(filename):
    symbols2={}
    firstLine=True
    fn = filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine=False
            header = t
            i=0; start=None; alt_start=None
            value_indexes=[]
            groups = {}
            group = 0
            for h in header:
                if h == 'WikiPathways': start=i
                if h == 'Select Protein Classes': alt_start=i
                i+=1
            if start == None: start = alt_start
            
            for h in header:
                if h>i:
                    group[i]
                i+=1
            if start == None: start = alt_start     
        
            
        else:
            
            uniprot = t[0]
            symbols = string.replace(t[-1],';;',';')
            symbols = string.split(symbols,';')
            for s in symbols:
                if len(s)>0:
                    symbols2[string.upper(s),uniprot]=[]
    for (s,u) in symbols2:
        ea.write(string.join([s,u],'\t')+'\n')    
    ea.close()
        
def coincentIncedenceTest(exp_file,TFs):
    fn = filepath(TFs)
    tfs={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        tfs[data]=[]
    
    comparisons={}
    for tf1 in tfs:
        for tf2 in tfs:
            if tf1!=tf2:
                temp = [tf1,tf2]
                temp.sort()
                comparisons[tuple(temp)]=[]

    gene_data={}
    firstLine=True
    fn = filepath(exp_file)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine:
            firstLine=False
            header = string.split(data,'\t')[1:]
        else:
            t = string.split(data,'\t')
            gene = t[0]
            values = map(float,t[1:])
            gene_data[gene] = values
                
    filename = TFs[:-4]+'-all-coincident-5z.txt'
    ea = export.ExportFile(filename)    
    comparison_db={}
    for comparison in comparisons:
        vals1 = gene_data[comparison[0]]
        vals2 = gene_data[comparison[1]]
        i=0
        coincident=[]
        for v1 in vals1:
            v2 = vals2[i]
            #print v1,v2
            if v1>1 and v2>1:
                coincident.append(i)
            i+=1
        i=0
        population_db={}; coincident_db={}
        for h in header:
            population=string.split(h,':')[0]
            if i in coincident:
                try: coincident_db[population]+=1
                except Exception: coincident_db[population]=1
            try: population_db[population]+=1
            except Exception: population_db[population]=1
            i+=1
            
        import mappfinder
        
        final_population_percent=[]
        for population in population_db:
            d = population_db[population]
            try: c = coincident_db[population]
            except Exception: c = 0
            
            N = float(len(header))  ### num all samples examined
            R = float(len(coincident))  ### num all coincedent samples for the TFs
            n = float(d) ### num all samples in cluster
            r = float(c) ### num all coincident samples in cluster
            try: z = mappfinder.Zscore(r,n,N,R)
            except Exception: z=0
            #if 'Gfi1b' in comparison and 'Gata1' in comparison: print N, R, n, r, z
            final_population_percent.append([population,str(c),str(d),str(float(c)/float(d)),str(z)])

        comparison_db[comparison]=final_population_percent

    filtered_comparison_db={}
    top_scoring_population={}
    for comparison in comparison_db:
        max_group=[]
        for population_stat in comparison_db[comparison]:
            z = float(population_stat[-1])
            c = float(population_stat[1])
            population = population_stat[0]
            max_group.append([z,population])
        max_group.sort()
        z = max_group[-1][0]
        pop = max_group[-1][1]
        if z>(1.96)*2 and c>3:
            filtered_comparison_db[comparison]=comparison_db[comparison]
            top_scoring_population[comparison] = pop,z
        
    firstLine = True
    for comparison in filtered_comparison_db:
        comparison_alt = string.join(list(comparison),'|')
        all_percents=[]
        for line in filtered_comparison_db[comparison]:
            all_percents.append(line[3])
        if firstLine:
            all_headers=[]
            for line in filtered_comparison_db[comparison]:
                all_headers.append(line[0])
            ea.write(string.join(['gene-pair']+all_headers+['Top Population','Top Z'],'\t')+'\n')    
            firstLine=False
        pop,z = top_scoring_population[comparison]
        ea.write(string.join([comparison_alt]+all_percents+[pop,str(z)],'\t')+'\n')    
    ea.close()

def getlastexon(filename):
    filename2 = filename[:-4]+'-last-exon.txt'
    ea = export.ExportFile(filename2)    
    firstLine=True
    fn = filepath(filename)
    last_gene = 'null'; last_exon=''
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine=False
        else:
            gene = t[2]
            if gene != last_gene:
                if ':E' in last_exon:
                    gene,exon = last_exon = string.split(':E')
                    block,region = string.split(exon,'.')
                    try: ea.write(last_exon+'\n')
                    except: pass
            last_gene = gene
            last_exon = t[0]
    ea.close()

def replaceWithBinary(filename):
    filename2 = filename[:-4]+'-binary.txt'
    ea = export.ExportFile(filename2)    
    firstLine=True
    fn = filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            ea.write(line)
            firstLine=False
        else:
            try: values = map(float,t[1:])
            except Exception: print t[1:];sys.exit()
            values2=[]
            for v in values:
                if v == 0: values2.append('0')
                else: values2.append('1')
            ea.write(string.join([t[0]]+values2,'\t')+'\n')
    ea.close()

def geneMethylationOutput(filename):
    filename2 = filename[:-4]+'-binary.txt'
    ea = export.ExportFile(filename2)    
    firstLine=True
    fn = filepath(filename)
    db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        values = (t[20],t[3]+'-methylation')
        db[values]=[]
    for value in db:
        ea.write(string.join(list(value),'\t')+'\n')
    ea.close()
    
def coincidentIncedence(filename,genes):
    exportPairs=True
    gene_data=[]
    firstLine=True
    fn = filepath(filename)
    if exportPairs:
        filename = filename[:-4]+'_'+genes[0]+'-'+genes[1]+'2.txt'
        ea = export.ExportFile(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine:
            firstLine=False
            header = string.split(data,'\t')[1:]
        else:
            t = string.split(data,'\t')
            gene = t[0]
            if gene in genes:
                values = map(float,t[1:])
                gene_data.append(values)

    vals1 = gene_data[0]
    vals2 = gene_data[1]
    i=0
    coincident=[]
    for v1 in vals1:
        v2 = vals2[i]
        #print v1,v2
        if v1>1 and v2>1:
            coincident.append(i)
        i+=1
    i=0
    population_db={}; coincident_db={}
    for h in header:
        population=string.split(h,':')[0]
        if i in coincident:
            try: coincident_db[population]+=1
            except Exception: coincident_db[population]=1
        try: population_db[population]+=1
        except Exception: population_db[population]=1
        i+=1
        
    import mappfinder
    
    final_population_percent=[]
    for population in population_db:
        d = population_db[population]
        try: c = coincident_db[population]
        except Exception: c = 0
        
        N = float(len(header))  ### num all samples examined
        R = float(len(coincident))  ### num all coincedent samples for the TFs
        n = d ### num all samples in cluster
        r = c ### num all coincident samples in cluster
        try: z = mappfinder.zscore(r,n,N,R)
        except Exception: z = 0
        
        final_population_percent.append([population,str(c),str(d),str(float(c)/float(d)),str(z)])
        
    if exportPairs:
        for line in final_population_percent:
            ea.write(string.join(line,'\t')+'\n')    
        ea.close()
    else:
        return final_population_percent

def extractFeatures(countinp,IGH_gene_file):
    import export
    ExonsPresent=False
    igh_genes=[]
    firstLine = True
    for line in open(IGH_gene_file,'rU').xreadlines():
        if firstLine: firstLine=False
        else:
            data = cleanUpLine(line)
            gene = string.split(data,'\t')[0]
            igh_genes.append(gene)
            
    if 'counts.' in countinp:
        feature_file = string.replace(countinp,'counts.','IGH.')
        fe = export.ExportFile(feature_file)
        firstLine = True
        for line in open(countinp,'rU').xreadlines():
            if firstLine:
                fe.write(line)
                firstLine=False
            else:
                feature_info = string.split(line,'\t')[0]
                gene = string.split(feature_info,':')[0]
                if gene in igh_genes:
                    fe.write(line)
        fe.close()
                       
def filterForJunctions(countinp):
    import export
    ExonsPresent=False
    igh_genes=[]
    firstLine = True
    count = 0
    if 'counts.' in countinp:
        feature_file = countinp[:-4]+'-output.txt'
        fe = export.ExportFile(feature_file)
        firstLine = True
        for line in open(countinp,'rU').xreadlines():
            if firstLine:
                fe.write(line)
                firstLine=False
            else:
                feature_info = string.split(line,'\t')[0]
                junction = string.split(feature_info,'=')[0]
                if '-' in junction:
                    fe.write(line)
                    count+=1
        fe.close()
        
    print count
        
def countIntronsExons(filename):
    import export
    exon_db={}
    intron_db={}
    firstLine = True
    last_transcript=None
    for line in open(filename,'rU').xreadlines():
        if firstLine:
            firstLine=False
        else:
            line = line.rstrip()
            t = string.split(line,'\t')
            transcript = t[-1]
            chr = t[1]
            strand = t[2]
            start = t[3]
            end = t[4]
            exon_db[chr,start,end]=[]
            if transcript==last_transcript:
                if strand == '1':
                    intron_db[chr,last_end,start]=[]
                else:
                    intron_db[chr,last_start,end]=[]
            last_end = end
            last_start = start
            last_transcript = transcript
            
    print len(exon_db)+1, len(intron_db)+1
             
def importGeneList(gene_list_file,n=20):
    genesets=[]
    genes=[]
    for line in open(gene_list_file,'rU').xreadlines():
        gene = line.rstrip()
        gene = string.split(gene,'\t')[0]
        genes.append(gene)
        if len(genes)==n:
            genesets.append(genes)
            genes=[]
    if len(genes)>0 and len(genes)<(n+1):
        genes+=(n-len(genes))*[gene]
        genesets.append(genes)
    return genesets
            
def customClean(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
            #print len(t)
            ea.write(string.join(['UID']+t,'\t')+'\n')
        else:
            if ';' in t[0]:
                uid = string.split(t[0],';')[0]
            else:
                uid = t[0]
            values = map(lambda x: float(x),t[1:])
            values.sort()
            if values[3]>=1:
                ea.write(string.join([uid]+t[1:],'\t')+'\n')
    ea.close()

def MakeJunctionFasta(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'.fasta'
    ea = export.ExportFile(filename)
    
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        probeset, seq = string.split(data,'\t')[:2]
        ea.write(">"+probeset+'\n')
        ea.write(string.upper(seq)+'\n')
    ea.close()
    
def ToppGeneFilter(filename):
    import gene_associations; from import_scripts import OBO_import
    gene_to_symbol = gene_associations.getGeneToUid('Hs',('hide','Ensembl-Symbol'))
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
            #print len(t)
            ea.write(string.join(['Ensembl\t\tCategory'],'\t')+'\n')
        else:
            symbol = t[1]; category = t[3]
            symbol = symbol[0]+string.lower(symbol[1:]) ### Mouse
            category = category[:100]
            if symbol in symbol_to_gene:
                ensembl = symbol_to_gene[symbol][0]
                ea.write(string.join([ensembl,symbol,category],'\t')+'\n')
    ea.close()

def CountKallistoAlignedJunctions(filename):
    fn = filepath(filename)
    firstRow=True
    #filename = filename[:-4]+'.fasta'
    ea = export.ExportFile(filename)
    
    found = False
    counts=0
    unique={}
    ea = export.ExportFile(filename[:-4]+'-Mpo.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if 'ENS' in line and 'JUNC1201' in line:
            ea.write(line)
            unique[t[0]]=[]
            counts+=1
    print counts, len(unique)
    ea.close()
    
def filterRandomFile(filename,col1,col2):
    fn = filepath(filename)
    firstRow=True
    
    counts=0
    ea = export.ExportFile(filename[:-4]+'-columns.txt')
    for line in open(fn,'rU').xreadlines():
        if line[0]!='#':
            data = line.rstrip()
            t = string.split(data,',')
            #print t[col1-1]+'\t'+t[col2-1];sys.exit()
            if ' ' in t[col2-1]:
                t[col2-1] = string.split(t[col2-1],' ')[2]
            ea.write(t[col1-1]+'\t'+t[col2-1]+'\n')
            counts+=1
    #print counts, len(unique)
    ea.close()
    
def getBlockExonPositions():
    fn = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart65/ensembl/Mm/Mm_Ensembl_exon.txt'
    firstRow=True
    filename = fn[:-4]+'.block.txt'
    ea = export.ExportFile(filename)
    
    found = False
    lines=0
    exon_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene,exonid,chromosome,strand,start,stop, a, b, c, d = string.split(data,'\t')
        exonid = string.split(exonid,'.')[0]
        uid = gene+':'+exonid
        if lines>0:
            try:
                exon_db[uid,strand].append(int(start))
                exon_db[uid,strand].append(int(stop))
            except Exception:
                exon_db[uid,strand] = [int(start)]
                exon_db[uid,strand].append(int(stop))
        lines+=1
    print len(exon_db)
    for (uid,strand) in exon_db:
        exon_db[uid,strand].sort()
        if strand == '-':
            exon_db[uid,strand].reverse()
        start = str(exon_db[uid,strand][0])
        stop = str(exon_db[uid,strand][1])
        coord = [start,stop]; coord.sort()
        
        ea.write(uid+'\t'+strand+'\t'+coord[0]+'\t'+coord[1]+'\n')
    
    ea.close()
    
def combineVariants(fn):
    firstRow=True
    filename = fn[:-4]+'.gene-level.txt'
    ea = export.ExportFile(filename)
    
    found = False
    lines=0
    gene_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        gene = t[9]
        if lines == 0:
            header = ['UID']+t[16:]
            header = string.join(header,'\t')+'\n'
            ea.write(header)
            lines+=1
        else:
            var_calls = map(float,t[16:])
            if gene in gene_db:
                count_sum_array = gene_db[gene]
                count_sum_array = [sum(value) for value in zip(*[count_sum_array,var_calls])]
                gene_db[gene] = count_sum_array
            else:
                gene_db[gene] = var_calls
    
    for gene in gene_db:
        var_calls = gene_db[gene]
        var_calls2=[]
        for i in var_calls:
            if i==0: var_calls2.append('0')
            else: var_calls2.append('1')
        ea.write(gene+'\t'+string.join(var_calls2,'\t')+'\n')
    ea.close()

def compareFusions(fn):
    firstRow=True
    filename = fn[:-4]+'.matrix.txt'
    ea = export.ExportFile(filename)
    
    found = False
    lines=0
    fusion_db={}
    sample_list=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if 'Gene_Fusion_Pair' in line:
            headers = string.split(data,'\t')[1:]
        try:
            sample, fusion = string.split(data,'\t')
            try: fusion_db[fusion].append(sample)
            except Exception: fusion_db[fusion] = [sample]
            if sample not in sample_list: sample_list.append(sample)    
        except Exception:
            t = string.split(data,'\t')
            fusion = t[0]
            index=0
            for i in t[1:]:
                if i=='1':
                    sample = headers[index]
                    try: fusion_db[fusion].append(sample)
                    except Exception: fusion_db[fusion] = [sample]
                    if sample not in sample_list: sample_list.append(sample)    
                index+=1 
    fusion_db2=[]
    for fusion in fusion_db:
        samples = fusion_db[fusion]
        samples2=[]
        for s in sample_list:
            if s in samples: samples2.append('1')
            else: samples2.append('0')
        fusion_db[fusion] = samples2
        
    ea.write(string.join(['Fusion']+sample_list,'\t')+'\n')
    for fusion in fusion_db:
        print [fusion]
        ea.write(fusion+'\t'+string.join(fusion_db[fusion],'\t')+'\n')
    ea.close()

def customCleanSupplemental(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        line = string.split(data,', ')
        gene_data=[]
        for gene in line:
            gene = string.replace(gene,' ','')
            if '/' in gene:
                genes = string.split(gene,'/')
                gene_data.append(genes[0])
                for i in genes[1:]:
                    gene_data.append(genes[0][:len(genes[1])*-1]+i)
            elif '(' in gene:
                genes = string.split(gene[:-1],'(')
                gene_data+=genes
            else:
                gene_data.append(gene)
            
        ea.write(string.join(gene_data,' ')+'\n')
    ea.close()

def customCleanBinomial(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    from stats_scripts import statistics
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            headers = t
            firstRow = False
            ea.write(string.join(['uid']+headers,'\t')+'\n')
        else:
            gene = t[0]
            values = map(float,t[1:])
            min_val = abs(min(values))
            values = map(lambda x: x+min_val,values)
            values = map(str,values)
            ea.write(string.join([gene]+values,'\t')+'\n')
    ea.close()

class MarkerFinderInfo:
    def __init__(self,gene,rho,tissue):
        self.gene = gene
        self.rho = rho
        self.tissue = tissue
    def Gene(self): return self.gene
    def Rho(self): return self.rho
    def Tissue(self): return self.tissue
    
def ReceptorLigandCellInteractions(species,lig_receptor_dir,cell_type_gene_dir):
    ligand_db={}
    receptor_db={}
    fn = filepath(lig_receptor_dir)    
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        ligand,receptor = string.split(data,'\t')
        if species=='Mm':
            ligand = ligand[0]+string.lower(ligand[1:])
            receptor = receptor[0]+string.lower(receptor[1:])
        try: ligand_db[ligand].apepnd(receptor)
        except Exception: ligand_db[ligand] = [receptor]
        try: receptor_db[receptor].append(ligand)
        except Exception: receptor_db[receptor] = [ligand]
    
    firstRow=True
    filename = cell_type_gene_dir[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    
    found = False
    cell_specific_ligands={}
    cell_specific_receptor={}
    fn = filepath(cell_type_gene_dir) 
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        gene, rho, tissue, notes, order = string.split(data,'\t')
        mf = MarkerFinderInfo(gene, rho, tissue)
        if gene in ligand_db:
            cell_specific_ligands[gene]=mf
        if gene in receptor_db:
            cell_specific_receptor[gene]=mf
    
    ligand_receptor_pairs=[]
    for gene in cell_specific_ligands:
        receptors = ligand_db[gene]
        for receptor in receptors:
            if receptor in cell_specific_receptor:
                rmf = cell_specific_receptor[receptor]
                lmf = cell_specific_ligands[gene]
                gene_data = [gene,lmf.Tissue(),lmf.Rho(),receptor,rmf.Tissue(),rmf.Rho()]
                pair = gene,receptor
                if pair not in ligand_receptor_pairs:
                    ea.write(string.join(gene_data,'\t')+'\n')
                    ligand_receptor_pairs.append(pair)
    for receptor in cell_specific_receptor:
        ligands = receptor_db[receptor]
        for gene in ligands:
            if gene in cell_specific_ligands:
                rmf = cell_specific_receptor[receptor]
                lmf = cell_specific_ligands[gene]
                gene_data = [gene,lmf.Tissue(),lmf.Rho(),receptor,rmf.Tissue(),rmf.Rho()]
                pair = gene,receptor
                if pair not in ligand_receptor_pairs:
                    ea.write(string.join(gene_data,'\t')+'\n')
                    ligand_receptor_pairs.append(pair)
    ea.close()
    
def findReciprocal(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-filtered.txt'
    ea = export.ExportFile(filename)
    
    found = False
    gene_ko={}; gene_oe={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
            headers={}
            TFs={}
            i=0
            for v in t[1:]:
                TF,direction = string.split(v,'-')
                headers[i]=TF,direction,v
                i+=1
                if v not in TFs:
                    f = filename[:-4]+'-'+v+'-up.txt'
                    tea = export.ExportFile(f)
                    TFs[v+'-up']=tea
                    tea.write('GeneID\tEn\n')
                    f = filename[:-4]+'-'+v+'-down.txt'
                    tea = export.ExportFile(f)
                    TFs[v+'-down']=tea
                    tea.write('GeneID\tEn\n')
        else:
            values = map(float,t[1:])
            gene = t[0]
            i=0
            for v in values:
                TF,direction,name = headers[i]
                if 'KO' in direction:
                    if v > 1:
                        gene_ko[gene,TF,1]=[]
                        tea = TFs[name+'-up']
                        tea.write(gene+'\tEn\n')
                    else:
                        gene_ko[gene,TF,-1]=[]
                        tea = TFs[name+'-down']
                        tea.write(gene+'\tEn\n')
                if 'OE' in direction:
                    if v > 1:
                        gene_oe[gene,TF,1]=[]
                        tea = TFs[name+'-up']
                        tea.write(gene+'\tEn\n')
                    else:
                        gene_oe[gene,TF,-1]=[]
                        tea = TFs[name+'-down']
                        tea.write(gene+'\tEn\n')
                i+=1

    print len(gene_oe)
    for (gene,TF,direction) in gene_oe:
        alt_dir=direction*-1
        if (gene,TF,alt_dir) in gene_ko:
            ea.write(string.join([TF,gene,str(direction)],'\t')+'\n')
    
    ea.close()
    for TF in TFs:
        TFs[TF].close()
    
def effectsPrioritization(filename):
    fn = filepath(filename)
    firstRow=True
    filename = filename[:-4]+'-new.txt'
    ea = export.ExportFile(filename)
    from stats_scripts import statistics
    found = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            headers = t[1:]
            firstRow = False
        else:
            gene = t[0]
            values = map(float,t[1:])
            max_val = abs(max(values))
            max_header = headers[values.index(max_val)]
            ea.write(gene+'\t'+max_header+'\t'+str(max_val)+'\n')
    ea.close()
    
def simpleCombine(folder):
    filename = folder+'/combined/combined.txt'
    ea = export.ExportFile(filename)
    headers=['UID']
    data_db={}
    files = UI.read_directory(folder)
    for file in files: #:70895507-70895600
        if '.txt' in file:
            fn = filepath(folder+'/'+file)
            print fn
            firstRow=True
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if firstRow:
                    for i in t[1:]:
                        headers.append(i+'.'+file[:-4])
                    firstRow = False
                else:
                    gene = t[0]
                    try: data_db[gene]+=t[1:]
                    except Exception: data_db[gene] = t[1:]

    len_db={}
    ea.write(string.join(headers,'\t')+'\n')
    for gene in data_db:
        if len(data_db[gene])==(len(headers)-1):
            values = map(float,data_db[gene])
            count=0
            for i in values:
                if i>0.9: count+=1
            if count>7:
                ea.write(string.join([gene]+data_db[gene],'\t')+'\n')
        len_db[len(data_db[gene])]=[]
    print len(len_db)
    for i in len_db:
        print i
    
    ea.close()
    
def simpleCombineFiles(folder):
    filename = folder+'/combined/combined.txt'
    ea = export.ExportFile(filename)
    files = UI.read_directory(folder)
    firstRow=True
    for file in files:
        if '.txt' in file:
            fn = filepath(folder+'/'+file)
            firstRow=True
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if firstRow:
                    t.append('Comparison')
                    ea.write(string.join(t,'\t')+'\n')
                    firstRow = False
                else:
                    t.append(file[:-4])
                    ea.write(string.join(t,'\t')+'\n')
    ea.close()
                    
def evaluateMultiLinRegulatoryStructure(all_genes_TPM,MarkerFinder,SignatureGenes,state,query=None):
    """Predict multi-lineage cells and their associated coincident lineage-defining TFs"""
    
    ICGS_State_as_Row = True
    ### Import all genes with TPM values for all cells
    matrix, column_header, row_header, dataset_name, group_db = importData(all_genes_TPM)
    group_index={}
    all_indexes=[]
    for sampleName in group_db:
        ICGS_state = group_db[sampleName][0]
        try: group_index[ICGS_state].append(column_header.index(sampleName))
        except Exception: group_index[ICGS_state] = [column_header.index(sampleName)]
        all_indexes.append(column_header.index(sampleName))
    for ICGS_state in group_index:
        group_index[ICGS_state].sort()
    all_indexes.sort()
        
    def importGeneLists(fn):
        genes={}
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            gene,cluster = string.split(data,'\t')[0:2]
            genes[gene]=cluster
        return genes
    
    def importMarkerFinderHits(fn):
        genes={}
        skip=True
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            if skip: skip=False
            else:
                gene,symbol,rho,ICGS_State = string.split(data,'\t')
                #if ICGS_State!=state and float(rho)>0.0:
                if float(rho)>0.0:
                    genes[gene]=float(rho),ICGS_State ### Retain all population specific genes (lax)
                    genes[symbol]=float(rho),ICGS_State
        return genes
    
    def importQueryDataset(fn):
        matrix, column_header, row_header, dataset_name, group_db = importData(fn)
        return matrix, column_header, row_header, dataset_name, group_db
    
    signatureGenes = importGeneLists(SignatureGenes)
    markerFinderGenes = importMarkerFinderHits(MarkerFinder)
    #print len(signatureGenes),len(markerFinderGenes)

    ### Determine for each gene, its population frequency per cell state
    index=0
    expressedGenesPerState={}
    
    def freqCutoff(x,cutoff):
        if x>cutoff: return 1 ### minimum expression cutoff
        else: return 0
        
    for row in matrix:
        ICGS_state_gene_frq={}
        gene = row_header[index]
        for ICGS_state in group_index:
            state_values = map(lambda i: row[i],group_index[ICGS_state])    
            def freqCheck(x):
                if x>1: return 1 ### minimum expression cutoff
                else: return 0
                
            expStateCells = sum(map(lambda x: freqCheck(x),state_values))
            statePercentage = (float(expStateCells)/len(group_index[ICGS_state]))
            ICGS_state_gene_frq[ICGS_state] = statePercentage
        
        multilin_frq = ICGS_state_gene_frq[state]
 
        datasets_values = map(lambda i: row[i],all_indexes)    
        all_cells_frq = sum(map(lambda x: freqCheck(x),datasets_values))/(len(datasets_values)*1.0)
        all_states_frq = map(lambda x: ICGS_state_gene_frq[x],ICGS_state_gene_frq)
        all_states_frq.sort() ### frequencies of all non-multilin states
        rank = all_states_frq.index(multilin_frq)
        states_expressed = sum(map(lambda x: freqCutoff(x,0.5),all_states_frq))/(len(all_states_frq)*1.0)

        if multilin_frq > 0.25 and rank>0: #and states_expressed<0.75 #and all_cells_frq>0.75
            if 'Rik' not in gene and 'Gm' not in gene:
                if gene in signatureGenes:# and gene in markerFinderGenes:
                    if ICGS_State_as_Row:
                        ICGS_State = signatureGenes[gene]
                    if gene in markerFinderGenes:
                        if ICGS_State_as_Row == False:
                            rho, ICGS_State = markerFinderGenes[gene]
                        else:
                            rho, ICGS_Cell_State = markerFinderGenes[gene]
                        score = int(rho*100*multilin_frq)*(float(rank)/len(all_states_frq))
                        try: expressedGenesPerState[ICGS_State].append((score,gene))
                        except Exception: expressedGenesPerState[ICGS_State]=[(score,gene)] #(rank*multilin_frq)
        index+=1

    if query!=None:
        matrix, column_header, row_header, dataset_name, group_db = importQueryDataset(query)
        
    createPseudoCell=True
    ### The expressedGenesPerState defines genes and modules co-expressed in the multi-Lin
    ### Next, find the cells that are most frequent in mulitple states
    representativeMarkers={}
    for ICGS_State in expressedGenesPerState:
        expressedGenesPerState[ICGS_State].sort()
        expressedGenesPerState[ICGS_State].reverse()
        if '1Multi' not in ICGS_State:
            markers = expressedGenesPerState[ICGS_State][:5]
            print ICGS_State,":",string.join(map(lambda x: x[1],list(markers)),', ')
            if createPseudoCell:
                for gene in markers:
                    def getBinary(x):
                        if x>1: return 1
                        else: return 0
                    if gene[1] in row_header: ### Only for query datasets
                        row_index = row_header.index(gene[1])
                        binaryValues = map(lambda x: getBinary(x), matrix[row_index])
                        #if gene[1]=='S100a8': print binaryValues;sys.exit()
                        try: representativeMarkers[ICGS_State].append(binaryValues)
                        except Exception: representativeMarkers[ICGS_State] = [binaryValues]    
            else:
                representativeMarkers[ICGS_State]=markers[0][-1]
        #int(len(markers)*.25)>5:
        #print ICGS_State, markers
    #sys.exit()
    
    for ICGS_State in representativeMarkers:
        if createPseudoCell:
            signature_values = representativeMarkers[ICGS_State]
            signature_values = [int(numpy.median(value)) for value in zip(*signature_values)]
            representativeMarkers[ICGS_State] = signature_values
        else:
            gene = representativeMarkers[ICGS_State]
            row_index = row_header.index(gene)
            gene_values = matrix[row_index]
            representativeMarkers[ICGS_State] = gene_values
        
    ### Determine for each gene, its population frequency per cell state
    expressedStatesPerCell={}
    for ICGS_State in representativeMarkers:
        gene_values = representativeMarkers[ICGS_State]
        index=0
        for cell in column_header:
            log2_tpm = gene_values[index]
            if log2_tpm>=1:
                try: expressedStatesPerCell[cell].append(ICGS_State)
                except Exception: expressedStatesPerCell[cell] = [ICGS_State]
            index+=1
    
    cell_mutlilin_ranking=[]   
    for cell in expressedStatesPerCell:
        lineageCount = expressedStatesPerCell[cell]
        cell_mutlilin_ranking.append((len(lineageCount),cell))
    cell_mutlilin_ranking.sort()
    cell_mutlilin_ranking.reverse()
    for cell in cell_mutlilin_ranking:
        print cell[0], cell[1], string.join(expressedStatesPerCell[cell[1]],'|')
    
def compareGenomicLocationAndICGSClusters():
    species = 'Mm'
    array_type = 'RNASeq'
    from build_scripts import EnsemblImport
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,array_type,'key_by_array')

    markerfinder = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/ExpressionOutput/MarkerFinder/AllCorrelationsAnnotated-ProteinCodingOnly.txt'
    eo = export.ExportFile(markerfinder[:-4]+'-bidirectional_promoters.txt')
    firstRow=True
    chr_cellTypeSpecific={}
    for line in open(markerfinder,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        symbol = t[1]
        ensembl = t[0]
        try: rho = float(t[6])
        except Exception: pass
        cellType = t[7]
        if firstRow:
            firstRow = False
        else:
            if ensembl in gene_location_db and rho>0.2:
                chr,strand,start,end = gene_location_db[ensembl]
                start = int(start)
                end = int(end)
                #region = start[:-5]
                try:
                    db = chr_cellTypeSpecific[chr,cellType]
                    try: db[strand].append([start,end,symbol,ensembl])
                    except Exception: db[strand] = [[start,end,symbol,ensembl]]
                except Exception:
                    db={}
                    db[strand] = [[start,end,symbol,ensembl]]
                    chr_cellTypeSpecific[chr,cellType] = db

    bidirectional={}
    eo.write(string.join(['CellType','Chr','Ensembl1','Symbol1','Start1','End1','Strand1','Ensembl2','Symbol2','Start2','End2','Strand2'],'\t')+'\n')
    for (chr,cellType) in chr_cellTypeSpecific:
        db = chr_cellTypeSpecific[chr,cellType]
        if len(db)>1: ### hence two strands
            for (start,end,symbol,ens) in db['+']:
                for (start2,end2,symbol2,ens2) in db['-']:
                    if abs(start-end2)<100000 and start>end2:
                        eo.write(string.join([cellType,chr,ens,symbol,str(start),str(end),'+',ens2,symbol2,str(end2),str(start2),'-'],'\t')+'\n')
                        try: bidirectional[chr,cellType].append([start,end,symbol,ens,start2,end2,symbol2,ens2])
                        except Exception: bidirectional[chr,cellType] = [[start,end,symbol,ens,start2,end2,symbol2,ens2]]
    eo.close()
    
def filterCountsFile(filename):
    fn = filepath(filename)
    firstRow=True
    
    def countif(value,cutoff=9):
        if float(value)>cutoff: return 1
        else: return 0
        
    header = True
    unique_genes = {}
    ea = export.ExportFile(filename[:-4]+'-filtered.txt')

    for line in open(fn,'rU').xreadlines():
        if header:
            header = False
            ea.write(line)
        else:
            data = line.rstrip()
            t = string.split(data,'\t')
            gene = string.split(t[0],':')[0]
            unique_genes[gene]=[]
            expressedSamples = map(countif,t[1:])
            if sum(expressedSamples)>2: 
                ea.write(line)
    ea.close()
    print len(unique_genes),'unique genes.'
    
def filterPSIValues(filename):
    fn = filepath(filename)
    firstRow=True
          
    header = True
    rows=0
    filtered=0
    new_file = filename[:-4]+'-75p.txt'
    new_file_clust = new_file[:-4]+'-clustID.txt'
    ea = export.ExportFile(new_file)
    eac = export.ExportFile(new_file_clust)
    added=[]
    for line in open(fn,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
            t = [t[1]]+t[8:]
            header_length = len(t)-1
            minimum_values_present = int(0.75*int(header_length))
            not_detected = header_length-minimum_values_present
            new_line = string.join(t,'\t')+'\n'
            ea.write(new_line)
        else:
            cID = t[5]
            t = [t[1]]+t[8:]
            missing_values_at_the_end = (header_length+1)-len(t)
            missing = missing_values_at_the_end+t.count('')
            if missing<not_detected:
                #if cID not in added:
                added.append(cID)
                new_line = string.join(t,'\t')+'\n'
                ea.write(new_line)
                eac.write(t[0]+'\t'+cID+'\n')
                filtered+=1
        rows+=1
    print rows, filtered
    ea.close()
    eac.close()
    #removeRedundantCluster(new_file,new_file_clust)

def removeRedundantCluster(filename,clusterID_file):
    from scipy import stats
    import ExpressionBuilder
    sort_col=0
    export_count=0
    ### Sort the filtered PSI model by gene name
    ExpressionBuilder.exportSorted(filename, sort_col, excludeHeader=True)
    new_file = filename[:-4]+'-unique.txt'
    ea = export.ExportFile(new_file)
    
    event_clusterID_db={}
    for line in open(clusterID_file,'rU').xreadlines():
        data = line.rstrip()
        eventID,clusterID = string.split(data,'\t')
        event_clusterID_db[eventID]=clusterID

        def compareEvents(events_to_compare,export_count):
                ### This is where we compare the events and write out the unique entries
                if len(events_to_compare)==1:
                    ea.write(events_to_compare[0][-1])
                    export_count+=1
                else:
                    exclude={}
                    compared={}
                    for event1 in events_to_compare:
                        if event1[0] not in exclude:
                            ea.write(event1[-1])
                            exclude[event1[0]]=[]
                            export_count+=1
                        for event2 in events_to_compare:
                            if event2[0] not in exclude:
                                if event1[0] != event2[0] and (event1[0],event2[0]) not in compared:
                                    uid1,values1,line1 = event1
                                    uid2,values2,line2 = event2
                                    coefr=numpy.ma.corrcoef(values1,values2)
                                    #rho,p = stats.pearsonr(values1,values2)
                                    rho = coefr[0][1]
                                    if rho>0.6 or rho<-0.6:
                                        exclude[event2[0]]=[]
                                    compared[event1[0],event2[0]]=[]
                                    compared[event2[0],event1[0]]=[]
                                    
                    for event in events_to_compare:
                        if event[0] not in exclude:
                            ea.write(event[-1]) ### write out the line
                            exclude.append(event[0])
                            export_count+=1
                return export_count

    header = True
    rows=0
    filtered=0
    prior_cID = 0
    events_to_compare=[]
    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            ea.write(line)
            header_row = t
            header=False
        else:
            uid = t[0]
            cID = event_clusterID_db[uid]
            empty_offset = len(header_row)-len(t)
            t+=['']*empty_offset
            values = ['0.000101' if x=='' else x for x in t[1:]]
            values = map(float,values)
            values = numpy.ma.masked_values(values,0.000101)
            if prior_cID==0: prior_cID = cID ### Occurs for the first entry
            if cID == prior_cID:
                ### Replace empty values with 0
                events_to_compare.append((uid,values,line))
            else:
                export_count = compareEvents(events_to_compare,export_count)
                events_to_compare=[(uid,values,line)]
            prior_cID = cID

    if len(events_to_compare)>0: ### If the laster cluster set not written out yet
        export_count = compareEvents(events_to_compare,export_count)

    ea.close()
    print export_count,'Non-redundant splice-events exported'
    
def convertToGOElite(folder):
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file:
            gene_count=0; up_count=0; down_count=0
            new_filename = string.split(file[3:],"_")[0]+'.txt'
            ea = export.ExportFile(folder+'/GO-Elite/'+new_filename)
            fn = folder+'/'+file
            ea.write('GeneID\tSystemCode\n')
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                if firstLine:
                    firstLine= False
                    continue
                data = line.rstrip()
                t = string.split(data,'\t')
                if ':' in t[0]:
                    ea.write(string.split(t[0],':')[0]+'\tSy\n')
                else:
                    gene_count+=1
                    if '-' in t[2]: down_count+=1
                    else: up_count+=1
            ea.close()
            print file,'\t',gene_count,'\t',up_count,'\t',down_count
            
def geneExpressionSummary(folder):
    import collections
    event_db = collections.OrderedDict()
    groups_list=['']
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file and 'GE.' in file:
            ls=[]
            event_db[file[:-4]]=ls
            groups_list.append(file[:-4])
            fn = folder+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    fold_index = t.index('LogFold')
                    firstLine= False
                    continue
                uid = t[0]
                if float(t[fold_index])>0:
                    fold_dir = 1
                else:
                    fold_dir = -1
                ls.append((uid,fold_dir))
    for file in event_db:
        print file,'\t',len(event_db[file])
                    
def compareEventLists(folder):
    import collections
    event_db = collections.OrderedDict()
    groups_list=['']
    files = UI.read_directory(folder)
    file_headers = {}
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            ls={}
            event_db[file[:-4]]=ls
            groups_list.append(file[:-4])
            fn = folder+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    file_headers[file[:-4]] = t ### Store the headers
                    try: event_index = t.index('Event-Direction')
                    except:
                        try: event_index = t.index('Inclusion-Junction') ### legacy
                        except: print file, 'Event-Direction error';sys.exit()
                    firstLine= False
                    continue
                uid = t[0]
                uid = string.split(uid,'|')[0]
                if 'U2AF1-l' in file or 'U2AF1-E' in file:
                    if t[2] == "inclusion":
                        ls[(uid,t[event_index])]=t ### Keep the event data for output
                else:
                    ls[(uid,t[event_index])]=t ### Keep the event data for output

    def convertEvents(events):
        opposite_events=[]
        for (event,direction) in events:
            if direction == 'exclusion':
                direction = 'inclusion'
            else:
                direction = 'exclusion'
            opposite_events.append((event,direction))
        return opposite_events
        
    ea1 = export.ExportFile(folder+'/overlaps-same-direction.txt')
    ea2 = export.ExportFile(folder+'/overlaps-opposite-direction.txt')
    ea3 = export.ExportFile(folder+'/concordance.txt')
    #ea4 = export.ExportFile(folder+'/overlap-same-direction-events.txt')
    ea1.write(string.join(groups_list,'\t')+'\n')
    ea2.write(string.join(groups_list,'\t')+'\n')
    ea3.write(string.join(groups_list,'\t')+'\n')
    
    comparison_db={}
    best_hits={}
    for comparison1 in event_db:
        events1 = event_db[comparison1]
        hits1=[comparison1]
        hits2=[comparison1]
        hits3=[comparison1]
        best_hits[comparison1]=[]
        for comparison2 in event_db:
            events2 = event_db[comparison2]
            events3 = convertEvents(events2)
            overlapping_events = list(set(events1).intersection(events2))                  
            overlap = len(overlapping_events)
            inverse_overlap = len(set(events1).intersection(events3)) ### Get opposite events
            ### Calculate ratios based on the size of the smaller set
            min_events1 = min([len(events1),len(events2)]) 
            min_events2 = min([len(events1),len(events3)])
            denom = overlap+inverse_overlap
            if denom == 0: denom = 0.00001
            #comparison_db[comparison1,comparison2]=overlap
            if min_events1 == 0: min_events1 = 1
            if (overlap+inverse_overlap)<20:
                hits1.append('0.5')
                hits2.append('0.5')
                hits3.append('0.5|0.5')
            else:
                hits1.append(str((1.00*overlap)/min_events1))
                hits2.append(str((1.00*inverse_overlap)/min_events1))
                hits3.append(str(1.00*overlap/denom)+'|'+str(1.00*inverse_overlap/denom)+':'+str(overlap+inverse_overlap))
                if 'Leu' not in comparison2:
                    comp_name = string.split(comparison2,'_vs')[0]
                    best_hits[comparison1].append([abs(1.00*overlap/denom),'cor',comp_name])
                    best_hits[comparison1].append([abs(1.00*inverse_overlap/denom),'anti',comp_name])
            if comparison1 != comparison2:
                if len(overlapping_events)>0:
                    #ea4.write(string.join(['UID',comparison1]+file_headers[comparison1]+[comparison2]+file_headers[comparison2],'\t')+'\n')
                    pass
                overlapping_events.sort()
                for event in overlapping_events:
                    vals = string.join([event[0],comparison1]+event_db[comparison1][event]+[comparison2]+event_db[comparison2][event],'\t')
                    #ea4.write(vals+'\n')
                    pass
        ea1.write(string.join(hits1,'\t')+'\n')
        ea2.write(string.join(hits2,'\t')+'\n')
        ea3.write(string.join(hits3,'\t')+'\n')
    ea1.close()
    ea2.close()
    ea3.close()
    #ea4.close()
    for comparison in best_hits:
        best_hits[comparison].sort()
        best_hits[comparison].reverse()
        hits = best_hits[comparison][:10]
        hits2=[]
        for (score,dir,comp) in hits:
            h = str(score)[:4]+'|'+dir+'|'+comp
            hits2.append(h)
        print comparison,'\t',string.join(hits2,', ')
     
def convertGroupsToBinaryMatrix(groups_file,sample_order,cellHarmony=False):
    eo = export.ExportFile(groups_file[:-4]+'-matrix.txt')
    print groups_file[:-4]+'-matrix.txt'
    firstRow=True
    samples = []
    ### Import a file with the sample names in the groups file in the correct order
    for line in open(sample_order,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if 'row_clusters-flat' in t:
            samples=[]
            samples1 = t[2:]
            for name in samples1:
                if ':' in name:
                    group,name = string.split(name,':')
                samples.append(name)
            if cellHarmony==False:
                break
        elif 'column_clusters-flat' in t and cellHarmony:
            clusters = t[2:]
        elif groups_file == sample_order:
            samples.append(t[0])
        elif firstRow:
            samples = t[1:]
            firstRow=False
    
    ### Import a groups file
    import collections
    sample_groups = collections.OrderedDict()
    for line in open(groups_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        sample, groupNum, groupName = t[:3]
        if cellHarmony == False: ### JUST USE THE PROVIDED GROUPS FOR SAMPLES FOUND IN BOTH FILES
            if sample in samples:
                si=samples.index(sample) ### Index of the sample
                try: sample_groups[groupName][si] = '1' ### set that sample to 1
                except Exception:
                    sample_groups[groupName] = ['0']*len(samples)
                    sample_groups[groupName][si] = '1' ### set that sample to 1
        else: ### JUST GRAB THE GROUP NAMES FOR THE SAMPLE GROUPS NOT THE SAMPES
            sample_groups[groupNum]=groupName
            
    if cellHarmony:
        i=0
        for sample in samples1:
            cluster = clusters[i]
            group_name = sample_groups[cluster]
            eo.write(sample+'\t'+cluster+'\t'+group_name+'\n')
            i+=1
        eo.close()
        
    else:
        eo.write(string.join(['GroupName']+samples,'\t')+'\n')
        for group in sample_groups:
            eo.write(string.join([group]+sample_groups[group],'\t')+'\n')
        eo.close()
        
def returnIntronJunctionRatio(counts_file,species = 'Mm'):
    eo = export.ExportFile(counts_file[:-4]+'-intron-ratios.txt')
    ### Import a groups file
    header=True
    prior_gene=[]
    exon_junction_values=[]
    intron_junction_values=[]
    eoi = export.ExportFile(counts_file[:-4]+'-intron-ratios-gene.txt')
    rows=0
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    def logratio(list):
        try: return list[0]/list[1]
        except Exception: return 0
    for line in open(counts_file,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        junctionID = t[0]
        if header:
            eoi.write(line)
            samples = t[1:]
            #zero_ref =[0]*len(samples)
            global_intron_ratios={}
            i=0
            for val in samples:
                global_intron_ratios[i]=[]
                i+=1
            header = False
            continue
        else:
            uid,coords = string.split(junctionID,'=')
            genes = string.split(uid,':') ### can indicate trans-splicing
            if len(genes)>2: trans_splicing = True
            else: trans_splicing = False
            coords = string.split(coords,':')[1]
            coords = string.split(coords,'-')
            coords = map(int,coords)
            coord_diff = abs(coords[1]-coords[0])
        #ENSMUSG00000027770:I23.1-E24.1=chr3:62470748-62470747
        gene = string.split(junctionID,':')[0]
        rows+=1
        if rows == 1:
            prior_gene = gene
        if gene != prior_gene:
            #print gene
            ### merge all of the gene level counts for all samples
            if len(intron_junction_values)==0:
                #global_intron_ratios = [sum(value) for value in zip(*[global_intron_ratios,zero_ref])]
                pass
            else:
                intron_junction_values_original = list(intron_junction_values)
                exon_junction_values_original = list(exon_junction_values)
                intron_junction_values = [sum(i) for i in zip(*intron_junction_values)]
                exon_junction_values = [sum(i) for i in zip(*exon_junction_values)]
                intron_ratios = [logratio(value) for value in zip(*[intron_junction_values,exon_junction_values])]
                #if sum(intron_ratios)>3:

                intron_ratios2=[]
                if prior_gene in gene_to_symbol:
                    symbol = gene_to_symbol[prior_gene][0]
                else:
                    symbol = prior_gene
                i=0
                #"""
                if symbol == 'Pi4ka':
                    print samples[482:487]
                    for x in exon_junction_values_original:
                        print x[482:487]
                    print exon_junction_values[482:487]
                    print intron_ratios[482:487]
                #"""
                for val in intron_ratios:
                    if exon_junction_values[i]>9:
                        if val>0:
                            ### stringent requirement - make sure it's not just a few reads
                            if intron_junction_values[i]>9: 
                                intron_ratios2.append(val)
                            else:
                                intron_ratios2.append(0)
                        else:
                            intron_ratios2.append(0)
                    else:
                        """
                        if val>0:
                            print val
                            print intron_junction_values
                            print exon_junction_values;sys.exit()"""
                        intron_ratios2.append('')
                    i+=1

                eoi.write(string.join([symbol]+map(str,intron_ratios2),'\t')+'\n')
                i = 0
                for val in intron_ratios:
                    if exon_junction_values[i]!=0: ### Only consider values with a non-zero denominator
                        global_intron_ratios[i].append(intron_ratios[i])
                    i+=1

            exon_junction_values = []
            intron_junction_values = []
            prior_gene = gene
        values = map(float,t[1:])
        if 'I' in junctionID and '_' not in junctionID and coord_diff==1 and trans_splicing == False:
            intron_junction_values.append(values)
            exon_junction_values.append(values)
        elif trans_splicing == False:
            exon_junction_values.append(values)
    print rows, 'processed'

    import numpy
    i=0; global_intron_ratios_values=[]
    for val in samples:
        global_intron_ratios_values.append(100*numpy.mean(global_intron_ratios[i])) ### list of lists
        i+=1
        
    eo.write(string.join(['UID']+samples,'\t')+'\n')
    eo.write(string.join(['Global-Intron-Retention-Ratio']+map(str,global_intron_ratios_values),'\t')+'\n')
    eo.close()
    eoi.close()

def convertSymbolLog(input_file,ensembl_symbol):
    
    gene_symbol_db={}
    for line in open(ensembl_symbol,'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl,symbol = string.split(data,'\t')
        gene_symbol_db[ensembl]=symbol
    convert = False
    eo = export.ExportFile(input_file[:-4]+'-log2.txt')
    header=0
    added_symbols=[]
    not_found=[]
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data,'\t')
        gene = values[0]
        if header == 0:
            #eo.write(line)
            data = cleanUpLine(line)
            headers = []
            values = string.split(data,'\t')
            for v in values:
                if "exp." in v:
                    headers.append(string.split(v,'.exp.')[0])
                else:
                    headers.append(v)
            eo.write(string.join(headers,'\t')+'\n')
        header +=1
        if gene in gene_symbol_db:
            symbol = gene_symbol_db[gene]
            if symbol not in added_symbols: 
                added_symbols.append(symbol)
                values = map(lambda x: math.log(float(x)+1,2),values[1:])
                if max(values)> 0.5:
                    values = map(lambda x: str(x)[:5],values)
                    eo.write(string.join([symbol]+values,'\t')+'\n')
        elif convert==False and header>1:
            values = map(lambda x: math.log(float(x)+1,2),values[1:])
            if max(values)> 0.5:
                values = map(lambda x: str(x)[:5],values)
                eo.write(string.join([gene]+values,'\t')+'\n')
        else:
            not_found.append(gene)
    print len(not_found),not_found[:10]

    eo.close()

def outputForGOElite(folds_dir):
    
    matrix, column_header, row_header, dataset_name, group_db = importData(folds_dir,Normalize=False)
    matrix = zip(*matrix) ### transpose
    
    ci=0
    root_dir = findParentDir(folds_dir)
    for group_data in matrix:
        group_name = column_header[ci]
        eo = export.ExportFile(root_dir+'/folds/'+group_name+'.txt')
        gi=0
        eo.write('geneID'+'\tSy\t'+'log2-fold'+'\n')
        for fold in group_data:
            gene = row_header[gi]
            if fold>0:
                eo.write(gene+'\tSy\t'+str(fold)+'\n')
            gi+=1
        eo.close()
        ci+=1
    
def transposeMatrix(input_file):
    arrays=[]
    eo = export.ExportFile(input_file[:-4]+'-transposed.txt')
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data,'\t')
        arrays.append(values)
    t_arrays = zip(*arrays)
    for t in t_arrays:
        eo.write(string.join(t,'\t')+'\n')
    eo.close()

def simpleStatsSummary(input_file):
    cluster_counts={}
    header=True
    for line in open(input_file,'rU').xreadlines():
        data = cleanUpLine(line)
        if header:
            header = False
        else:
            sample,cluster,counts = string.split(data,'\t')
            try: cluster_counts[cluster].append(float(counts))
            except Exception: cluster_counts[cluster]=[float(counts)]
    
    for cluster in cluster_counts:
        avg = statistics.avg(cluster_counts[cluster])
        stdev = statistics.stdev(cluster_counts[cluster])
        print cluster+'\t'+str(avg)+'\t'+str(stdev)
        
def latteralMerge(file1, file2):
    import collections
    cluster_db = collections.OrderedDict()
    eo = export.ExportFile(file2[:-4]+'combined.txt')
    for line in open(file1,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        cluster_db[t[0]]=t
    for line in open(file2,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if t[0] in cluster_db:
            t1=cluster_db[t[0]]
            eo.write(string.join(t1+t[2:],'\t')+'\n')
    eo.close()
    
def removeMarkerFinderDoublets(heatmap_file,diff=1):
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = remoteImportData(heatmap_file)
        
    priorRowClusters.reverse()
    if len(priorColumnClusters)==0:
        for c in column_header:
            cluster = string.split(c,':')[0]
            priorColumnClusters.append(cluster)
        for r in row_header:
            cluster = string.split(r,':')[0]
            priorRowClusters.append(cluster)

    import collections
    cluster_db = collections.OrderedDict()
    i=0
    for cluster in priorRowClusters:
        try: cluster_db[cluster].append(matrix[i])
        except: cluster_db[cluster] = [matrix[i]]
        i+=1
    
    transposed_data_matrix=[]
    clusters=[]
    for cluster in cluster_db:
        cluster_cell_means = numpy.mean(cluster_db[cluster],axis=0)
        cluster_db[cluster] = cluster_cell_means
        transposed_data_matrix.append(cluster_cell_means)
        if cluster not in clusters:
            clusters.append(cluster)
    transposed_data_matrix = zip(*transposed_data_matrix)
    
    i=0
    cell_max_scores=[]
    cell_max_score_db = collections.OrderedDict()

    for cell_scores in transposed_data_matrix:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ci = clusters.index(cluster)
        #print ci, cell, cluster, cell_scores;sys.exit()
        cell_state_score = cell_scores[ci] ### This is the score for that cell for it's assigned MarkerFinder cluster
        alternate_state_scores=[]
        for score in cell_scores:
            if score != cell_state_score:
                alternate_state_scores.append(score)
        alt_max_score = max(alternate_state_scores)
        alt_sum_score = sum(alternate_state_scores)
        cell_max_scores.append([cell_state_score,alt_max_score,alt_sum_score]) ### max and secondary max score - max for the cell-state should be greater than secondary max
        try: cell_max_score_db[cluster].append(([cell_state_score,alt_max_score,alt_sum_score]))
        except: cell_max_score_db[cluster] = [[cell_state_score,alt_max_score,alt_sum_score]]
        i+=1
    
    
    for cluster in cell_max_score_db:
        cluster_cell_means = numpy.median(cell_max_score_db[cluster],axis=0)
        cell_max_score_db[cluster] = cluster_cell_means ### This is the cell-state mean score for all cells in that cluster and the alternative max mean score (difference gives you the threshold for detecting double)
    i=0
    print len(cell_max_scores)
    keep=['row_clusters-flat']
    keep_alt=['row_clusters-flat']
    remove = ['row_clusters-flat']
    remove_alt = ['row_clusters-flat']
    for (cell_score,alt_score,alt_sum) in cell_max_scores:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ref_max, ref_alt, ref_sum = cell_max_score_db[cluster]
        ci = clusters.index(cluster)
        ref_diff= math.pow(2,(ref_max-ref_alt))*diff #1.1
        ref_alt = math.pow(2,(ref_alt))
        cell_diff = math.pow(2,(cell_score-alt_score))
        cell_score = math.pow(2,cell_score)
        if cell_diff>ref_diff and cell_diff>diff: #cell_score cutoff removes some, but cell_diff is more crucial
                #if alt_sum<cell_score:
                assignment=0 #1.2
                keep.append(cell)
                try: keep_alt.append(string.split(cell,':')[1]) ### if prefix added
                except Exception:
                    keep_alt.append(cell)
        else:
            remove.append(cell)
            try: remove_alt.append(string.split(cell,':')[1])
            except Exception: remove_alt.append(cell)
            assignment=1

        #print assignment
        i+=1
    print len(keep), len(remove)
    from import_scripts import sampleIndexSelection
    input_file=heatmap_file
    output_file = heatmap_file[:-4]+'-Singlets.txt'
    try: sampleIndexSelection.filterFile(input_file,output_file,keep)
    except: sampleIndexSelection.filterFile(input_file,output_file,keep_alt)

    output_file = heatmap_file[:-4]+'-Multiplets.txt'
    try: sampleIndexSelection.filterFile(input_file,output_file,remove)
    except: sampleIndexSelection.filterFile(input_file,output_file,remove_alt)
    
if __name__ == '__main__':
    filename=''
    index1=2;index2=3; x_axis='Number of Differentially Expressed Genes'; y_axis = 'Comparisons'; title='Hippocampus - Number of Differentially Expressed Genes'
    #OutputFile = export.findParentDir(filename)
    #OutputFile = export.findParentDir(OutputFile[:-1])+'/test.pdf'
    
    #stackedbarchart(filename,display=True,output=OutputFile);sys.exit()
    barchart(filename,index1,index2,x_axis,y_axis,title,display=True,color1='IndianRed',color2='SkyBlue');sys.exit()
    diff=0.7
    print 'diff:',diff
    #latteralMerge(file1, file2);sys.exit()
    #removeMarkerFinderDoublets('/Volumes/salomonis2/Erica-data/Demuxlet8Human/Seurat/ICGS2/CellHarmonyReference/DataPlots/Clustering-exp.ICGS-cellHarmony-reference-filtered-hierarchical_cosine_cosine2.txt',diff=diff);sys.exit()
    #outputForGOElite('/Users/saljh8/Desktop/R412X/completed/centroids.WT.R412X.median.txt');sys.exit()
    #simpleStatsSummary('/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/HCA/Mean-Comparisons/ExpressionInput/MergedFiles.Counts.UMI.txt');sys.exit()
    a = '/Volumes/salomonis2/Grimes/Andre-10X/PROJ-00504/Project_s1115g01001_6lib_11lane_BCL/Combined_Captures/Day0/cellHarmony/exp.MarkerFinder-cellHarmony-reference__Day0-CPTT-log2-ReOrdered-Query.txt'
    b = '/Volumes/salomonis2/Immune-10x-data-Human-Atlas/Bone-Marrow/Stuart/Browser/ExpressionInput/HS-compatible_symbols.txt'
    #b = '/data/salomonis2/GSE107727_RAW-10X-Mm/filtered-counts/ExpressionInput/Mm_compatible_symbols.txt'
    #a = '/Volumes/salomonis2/Immune-10x-data-Human-Atlas/Bone-Marrow/Stuart/Browser/head.txt'
    transposeMatrix(a);sys.exit()
    convertSymbolLog(a,b);sys.exit()
    #returnIntronJunctionRatio('/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Fluidigm_scRNA-Seq/12.09.2107/counts.WT-R412X.txt');sys.exit()
    #geneExpressionSummary('/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/ExpressionInput/DEGs-LogFold_1.0_rawp');sys.exit()
    b = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/HCA/Plasma-cell/ExpressionInput/DataPlots/a.txt'
    a = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/HCA/Plasma-cell/ICGS-NMF1/groups.FinalMarkerHeatmap.txt'
    convertGroupsToBinaryMatrix(a,b,cellHarmony=False);sys.exit()
    a = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/tests/events.txt'
    b = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/tests/clusters.txt'
    #simpleCombineFiles('/Users/saljh8/Desktop/dataAnalysis/Collaborative/Jose/NewTranscriptome/CombinedDataset/ExpressionInput/Events-LogFold_0.58_rawp')
    #removeRedundantCluster(a,b);sys.exit()
    a = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/SpliceICGS.R1.Depleted.12.27.17/all-depleted-and-KD'
    #a = '/Users/saljh8/Desktop/circadian_splicing_concordance'
    compareEventLists(a);sys.exit()
    #filterPSIValues('/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/CORNEL-AML/PSI/exp.Cornell-Bulk.txt');sys.exit()
    #compareGenomicLocationAndICGSClusters();sys.exit()
    #ViolinPlot();sys.exit()
    #simpleScatter('/Users/saljh8/Downloads/CMdiff_paper/calcium_data-KO4.txt');sys.exit()
    query_dataset = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/exp.GSE81682_HTSeq-cellHarmony-filtered.txt'
    all_tpm = '/Users/saljh8/Desktop/demo/BoneMarrow/ExpressionInput/exp.BoneMarrow-scRNASeq.txt'
    markerfinder = '/Users/saljh8/Desktop/demo/BoneMarrow/ExpressionOutput1/MarkerFinder/AllGenes_correlations-ReplicateBasedOriginal.txt'
    signature_genes = '/Users/saljh8/Desktop/Grimes/KashishNormalization/test/Panorama.txt'
    state = 'Multi-Lin'
    #evaluateMultiLinRegulatoryStructure(all_tpm,markerfinder,signature_genes,state);sys.exit()
    query_dataset = None
    all_tpm = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/MultiLin/Gottgens_HarmonizeReference.txt'
    signature_genes = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/MultiLin/Gottgens_HarmonizeReference.txt'
    markerfinder = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    state = 'Eryth_Multi-Lin'
    #evaluateMultiLinRegulatoryStructure(all_tpm,markerfinder,signature_genes,state,query = query_dataset);sys.exit()
    #simpleCombine("/Volumes/My Passport/Ari-10X/input");sys.exit()
    #effectsPrioritization('/Users/saljh8/Documents/1-dataAnalysis/RBM20-collaboration/RBM20-BAG3_splicing/Missing Values-Splicing/Effects.txt');sys.exit()
    #customCleanBinomial('/Volumes/salomonis2-1/Lab backup/Theresa-Microbiome-DropSeq/NegBinomial/ExpressionInput/exp.Instesinal_microbiome2.txt');sys.exit()
    #findReciprocal('/Volumes/HomeBackup/CCHMC/Jared-KO/BatchCorrectedFiltered/exp.CM-KO-steady-state.txt');sys.exit()
    #ReceptorLigandCellInteractions('Mm','/Users/saljh8/Downloads/ncomms8866-s3.txt','/Users/saljh8/Downloads/Round3-MarkerFinder_All-Genes.txt');sys.exit()
    #compareFusions('/Volumes/salomonis2-2/CPMC_Melanoma-GBM/Third-batch-files/Complete_analysis/temp/Combined_Fusion_GBM.txt');sys.exit()
    #combineVariants('/Volumes/salomonis2/CPMC_Melanoma-GBM/Third-batch-files/Complete_analysis/Variant_results/GBM/Variants_HighModerate-GBM_selected.txt');sys.exit()
    #customCleanSupplemental('/Users/saljh8/Desktop/dataAnalysis/CPMC/TCGA_MM/MM_genes_published.txt');sys.exit()
    #customClean('/Users/saljh8/Desktop/dataAnalysis/Driscoll/R3/2000_run1708A_normalized.txt');sys.exit()
    #simpleFilter('/Volumes/SEQ-DATA 1/all_10.5_mapped_norm_GC.csv');sys.exit()
    #filterRandomFile('/Users/saljh8/Downloads/HuGene-1_1-st-v1.na36.hg19.transcript2.csv',1,8);sys.exit()
    filename = '/Users/saljh8/Desktop/Grimes/GEC14078/MergedFiles.txt'
    #CountKallistoAlignedJunctions(filename);sys.exit()
    filename = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart72/Mm/junction1/junction_critical-junction-seq.txt'
    #MakeJunctionFasta(filename);sys.exit()
    filename = '/Users/saljh8/Downloads/CoexpressionAtlas.txt'
    #ToppGeneFilter(filename); sys.exit()
    #countIntronsExons(filename);sys.exit()
    #filterForJunctions(filename);sys.exit()
    #filename = '/Users/saljh8/Desktop/Grimes/GEC14074/ExpressionOutput/LineageCorrelations-test-protein_coding-zscores.txt'
    #runHCOnly(filename,[]); sys.exit()
    folder = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart72/ensembl/Hs'
    try:
        files = UI.read_directory(folder)
        for file in files: #:70895507-70895600
            if '.bed' in file:
                #BedFileCheck(folder+'/'+file)
                pass
    except Exception: pass
    #sys.exit()
    #runPCAonly(filename,[],False,showLabels=False,plotType='2D');sys.exit()
    
    countinp = '/Volumes/salomonis2/SinghLab/20150715_single_GCBCell/bams/ExpressionInput/counts.Bcells.txt'
    IGH_gene_file = '/Volumes/salomonis2/SinghLab/20150715_single_GCBCell/bams/ExpressionInput/IGH_genes.txt'
    #extractFeatures(countinp,IGH_gene_file);sys.exit()
    
    import UI
    #geneMethylationOutput(filename);sys.exit()
    #ica(filename);sys.exit()
    #replaceWithBinary('/Users/saljh8/Downloads/Neg_Bi_wholegenome.txt');sys.exit()
    #simpleFilter('/Volumes/SEQ-DATA/AML-TCGA/ExpressionInput/counts.LAML1.txt');sys.exit()
    filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/genes.tpm_tracking-ordered.txt'
    
    #filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/6-5-2015/ExpressionInput/amplify/exp.All-wt-output.txt'
    #getlastexon(filename);sys.exit()
    TFs = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/TF-by-gene_matrix/all-TFs2.txt'
    folder = '/Users/saljh8/Downloads/BLASTX2_Gecko.tab'
    genes = ['Gfi1', 'Irf8'] #'Cebpe', 'Mecom', 'Vwf', 'Itga2b', 'Meis1', 'Gata2','Ctsg','Elane', 'Klf4','Gata1']
    #genes = ['Gata1','Gfi1b']
    #coincentIncedenceTest(filename,TFs);sys.exit()
    #coincidentIncedence(filename,genes);sys.exit()
    #test(folder);sys.exit()
    #files = UI.read_directory(folder)
    #for file in files: SimpleCorrdinateToBed(folder+'/'+file)
    #filename = '/Users/saljh8/Desktop/bed/RREs0.5_exons_unique.txt'
    #simpleIntegrityCheck(filename);sys.exit()
    
    gene_list = ['S100a8','Chd7','Ets1','Chd7','S100a8']
    gene_list_file = '/Users/saljh8/Desktop/demo/Amit/ExpressionInput/genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/Grimes/Comb-plots/AML_genes-interest.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Grimes/Mm_Sara-single-cell-AML/alt/AdditionalHOPACH/ExpressionInput/AML_combplots.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Grimes/MDS-array/Comb-plot genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Grimes/All-Fluidigm/ExpressionInput/comb_plot3.txt'
    gene_list_file = '/Users/saljh8/Desktop/Grimes/MultiLin-Code/MultiLin-TFs.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/ExpressionInput/genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/10X-DropSeq-comparison/Final-Classifications/genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/TFs/Myelo_TFs2.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/R412X/customGenes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/ExpressionInput/genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/R412X/genes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/HCA/BM1-8_CD34+/ExpressionInput/MixedLinPrimingGenes.txt'
    gene_list_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/Grimes/DEGs.txt'
    genesets = importGeneList(gene_list_file,n=37)
    filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/comb-plots/exp.IG2_GG1-extended-output.txt'
    filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/comb-plots/genes.tpm_tracking-ordered.txt'
    filename = '/Users/saljh8/Desktop/demo/Amit/ExpressedCells/GO-Elite_results/3k_selected_LineageGenes-CombPlotInput2.txt'
    filename = '/Users/saljh8/Desktop/Grimes/Comb-plots/exp.AML_single-cell-output.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Grimes/Mm_Sara-single-cell-AML/alt/AdditionalHOPACH/ExpressionInput/exp.AML.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Grimes/MDS-array/comb-plot/input.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Grimes/All-Fluidigm/ExpressionInput/exp.Lsk_panorama.txt'
    filename = '/Users/saljh8/Desktop/demo/BoneMarrow/ExpressionInput/exp.BoneMarrow-scRNASeq.txt'
    filename = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/exp.GSE81682_HTSeq-cellHarmony-filtered.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Harinder/scRNASeq_Mm-Plasma/PCA-loading/ExpressionInput/exp.PCA-Symbol.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/10X-DropSeq-comparison/Final-Classifications/cellHarmony/MF-analysis/ExpressionInput/exp.Fluidigm-log2-NearestNeighbor-800.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/10X-DropSeq-comparison/Final-Classifications/cellHarmony/MF-analysis/ExpressionInput/exp.10X-log2-NearestNeighbor.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/10X-DropSeq-comparison/DropSeq/MultiLinDetect/ExpressionInput/DataPlots/exp.DropSeq-2k-log2.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/R412X/exp.allcells-v2.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/HCA/BM1-8_CD34+/ExpressionInput/exp.CD34+.v5-log2.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/cellHarmony-evaluation/Grimes/exp.AML-relative.txt'
    #filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-10x/CITE-Seq-MF-indexed/ExpressionInput/exp.cellHarmony.v3.txt'
    #filename = '/Volumes/salomonis2/Theodosia-Kalfa/Combined-10X-CPTT/ExpressionInput/exp.MergedFiles-ICGS.txt'
    #filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/Ly6g/combined-ICGS-Final/R412X/exp.cellHarmony-WT-R412X-relative.txt'
    #filename = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/Ly6g/CodingOnly/Guide3-Kallisto-Coding-NatureAugmented/SubClustering/Nov-27-Final-version/ExpressionInput/exp.wt-panorama.txt'
    #filename = '/Volumes/salomonis2/Harinder-singh/Run2421-10X/10X_IRF4_Lo/outs/filtered_gene_bc_matrices/ExpressionInput/exp.10X_IRF4_Lo_matrix_CPTT-ICGS.txt'
    #filename = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/Ly6g/CodingOnly/Guide3-Kallisto-Coding-NatureAugmented/SubClustering/Nov-27-Final-version/R412X/exp.R412X-RSEM-order.txt'

    print genesets
    for gene_list in genesets:
        multipleSubPlots(filename,gene_list,SubPlotType='column',n=37)
    sys.exit()

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