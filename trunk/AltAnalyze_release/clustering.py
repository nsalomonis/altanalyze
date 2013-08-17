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

import traceback
try:
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        matplotlib.rcParams['backend'] = 'TkAgg'
        import matplotlib.pyplot as pylab
        import matplotlib.colors as mc
        import matplotlib.mlab as mlab
        from matplotlib import mpl
        from matplotlib.patches import Circle
        from mpl_toolkits.mplot3d import Axes3D
        import scipy
        import scipy.cluster.hierarchy as sch
        import scipy.spatial.distance as dist
        try: import numpy
        except Exception:
            print traceback.format_exc()
        try: import igraph.vendor.texttable
        except Exception: None
        #pylab.ion() # closes Tk window after show - could be nice to include
except Exception:
    print traceback.format_exc()
    None ### Not needed for buildGraphFromSIF

import string
import time

import unique
import sys, os
import warnings
try:
    warnings.simplefilter("ignore", numpy.ComplexWarning)
    warnings.simplefilter("ignore", DeprecationWarning) ### Annoying depreciation warnings (occurs in sch somewhere)
except Exception: None

import WikiPathways_webservice

try:
    import fastcluster as fc
    print 'Using fastcluster instead of scipy hierarchical cluster'
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
    
def heatmap(x, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name,display=False,contrast=None):
    print "Performing hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)
    show_color_bars = True ### Currently, the color bars don't exactly reflect the dendrogram colors

    try: os.mkdir(root_dir) ### May need to create this directory
    except Exception: None
    if display == False:
        pylab.figure() ### Add this to avoid a Tkinter bug after running MarkerFinder (not sure why it is needed) - creates a second empty window when display == True
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
            matrix, column_header, row_header, dataset_name, group_db = importData(newFilename)
            x = numpy.array(matrix)
            
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
    if color_gradient == 'seismic':
        cmap=pylab.cm.seismic
    if color_gradient == 'green_white_purple':
        cmap=pylab.cm.PiYG_r
    if color_gradient == 'coolwarm':
        cmap=pylab.cm.coolwarm

    vmin=x.min()
    vmax=x.max()
    vmax = max([vmax,abs(vmin)])
    vmin = vmax*-1
    default_window_hight = 8.5
    default_window_width = 12
    
    if contrast == None:
        scaling_factor = 2.5 #2.5
    else:
        try: scaling_factor = float(contrast)
        except Exception: scaling_factor = 2.5
    
    norm = mpl.colors.Normalize(vmin/scaling_factor, vmax/scaling_factor) ### adjust the max and min to scale these colors by 2.5 (1 scales to the highest change)
    fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    pylab.rcParams['font.size'] = 7.5

    if show_color_bars == False:
        color_bar_w = 0.000001 ### Invisible but not gone (otherwise an error persists)
    else:
        color_bar_w = 0.015 ### Sufficient size to show
        
    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.135] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.076] ### Last one controls the hight

    # Compute and plot top dendrogram
    if column_method == 'hopach':
        ind2 = numpy.array(Z2['level']) ### from R_interface - hopach root cluster level
    elif column_method != None:
        start_time = time.time()
        d2 = dist.pdist(x.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
        Y2 = fc.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        #Y2 = sch.fcluster(Y2, 10, criterion = "maxclust")
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
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
        Y1 = fc.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='right')
        #ind1 = sch.fcluster(Y1,0.6*D1.max(),'distance') ### get the correlations
        #ind1 = sch.fcluster(Y1,0.2*D1.max(),'maxclust') 
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
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
        ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])
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
    
    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if len(row_header)<40:
            radj = len(row_header)*0.009 ### row offset value to center the vertical position of the row label
        elif len(row_header)<70:
            radj = len(row_header)*0.007 ### row offset value to center the vertical position of the row label
        else:
            radj = len(row_header)*0.005
        if row_method != None:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[idx1[i]],fontsize=row_fontsize)
            new_row_header.append(row_header[idx1[i]])
        else:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i-radj, '  '+row_header[i],fontsize=row_fontsize) ### When not clustering rows
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if len(row_header)<6:
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
        elif len(row_header)<25:
            cadj = len(row_header)*-0.11 ### column offset value
        elif len(row_header)>200:
            cadj = -2
        else:
            cadj = -0.9
        if column_method != None:
            axm.text(i, cadj, ''+column_header[idx2[i]], rotation=270, verticalalignment="top",fontsize=column_fontsize) # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            axm.text(i, cadj, ''+column_header[i], rotation=270, verticalalignment="top",fontsize=column_fontsize)
            new_column_header.append(column_header[i])

    # Plot colside colors
    # axc --> axes for column side colorbar

    ind1,ind2,cb_status = updateColorBarData(ind1,ind2,new_column_header,new_row_header)
    if (column_method != None or cb_status=='column') and show_color_bars == True:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        if len(unique.unique(ind2))>7: ### cmap_c is too few colors
            cmap_c = pylab.cm.gist_ncar
        dc = numpy.array(ind2, dtype=int)
        #print ind2, len(ind2)
        dc.shape = (1,len(ind2)) 
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
        axc.set_xticks([]) ### Hides ticks
        axc.set_yticks([])
    
    #cmap_c = mpl.colors.ListedColormap(map(lambda x: GroupDB[x][-1], new_column_header))

    if show_color_bars == False:
        axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
        axc.set_frame_on(False)
    
    # Plot rowside colors
    # axr --> axes for row side colorbar
    if (row_method != None or cb_status=='row') and show_color_bars == True:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        try: dr = numpy.array(ind1, dtype=int)
        except Exception:
            print ind1;kill
        dr.shape = (len(ind1),1)
        #print ind1, len(ind1)
        cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
        if len(unique.unique(ind1))>7: ### cmap_r is too few colors
            cmap_r = pylab.cm.gist_ncar
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
        axr.set_xticks([]) ### Hides ticks
        axr.set_yticks([])
        
    if show_color_bars == False:
        axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
        axr.set_frame_on(False)

    # Plot color legend
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
    axcb.set_title("colorkey",fontsize=14)
    
    filename = 'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name,column_metric,row_metric)
    
    if 'LineageCorrelations' in dataset_name:
        cb.set_label("Lineage Correlation Z Scores",fontsize=11)
    elif 'Heatmap' in root_dir:
        cb.set_label("GO-Elite Z Scores",fontsize=11)
    else:
        cb.set_label("Differential Expression (log2 fold)",fontsize=11)
        exportFlatClusterData(root_dir + filename, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render and save the graphic
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename, dpi=100) #,dpi=200
    if 'Outlier' in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Outlier Genes Genes',root_dir+filename])
    elif 'Relative' in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Significant Genes (Relative comparisons)',root_dir+filename])
    elif 'LineageCorrelations' in filename:
        graphic_link.append(['Hierarchical Clustering - Lineage Correlations',root_dir+filename])
    else:
        graphic_link.append(['Hierarchical Clustering - Significant Genes',root_dir+filename])
    if display:
        print 'Exporting:',filename
        try: pylab.show()
        except Exception: None ### when run in headless mode
    #pylab.close()
    #sys.exit()
    
def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
    import export
    filename = string.replace(filename,'.pdf','.txt')
    export_text = export.ExportFile(filename)
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','']+ map(str, ind2),'\t')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)
    
    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i],str(ind1[i])]+map(str, row),'\t')+'\n')
        i+=1
    export_text.close()
    
    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    export_cdt = export.ExportFile(filename)
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_cdt.write(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
        i+=1
    export_cdt.close()
    
### How to create custom colors - http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

def updateColorBarData(ind1,ind2,column_header,row_header):
    """ Replace the top-level cluster information with group assignments for color bar coloring (if group data present)"""
    cb_status = 'original'
    group_number_db={}
    group_number_list=[]
    try: ### Error if GroupDB not recognized as global
        if column_header[0] in GroupDB: ### Thus group assignments exist for column headers
            cb_status = 'column'
            k=0
            for header in column_header:
                group,color = GroupDB[header]
                try: value = group_number_db[group]
                except KeyError:
                    group_number_db[group] = k
                    value = k
                    k+=1
                group_number_list.append(value) ### will replace ind2
            ind2 = group_number_list
        if row_header[0] in GroupDB: ### Thus group assignments exist for column headers
            cb_status = 'row'
            k=0
            for header in row_header:
                group,color = GroupDB[header]
                try: value = group_number_db[group]
                except KeyError:
                    group_number_db[group] = k
                    value = k
                    k+=1
                group_number_list.append(k) ### will replace ind2
            ind1 = group_number_list
    except Exception: None
    
    return ind1,ind2,cb_status

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

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def importData(filename,Normalize=False):
    start_time = time.time()
    fn = filepath(filename)
    matrix=[]
    row_header=[]
    x=0
    filename = string.replace(filename,'\\','/')
    dataset_name = string.split(filename,'/')[-1][:-4]
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            group_db, column_header = assignGroupColors(t[1:])
            x=1
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                s = map(float,t[1:])
                #if (abs(max(s)-min(s)))>2:
                if Normalize!=False:
                    if Normalize=='row mean':
                        #avg = min(s)
                        avg = numpy.average(s)
                    else: avg = avg = numpy.median(s)
                    s = map(lambda x: x-avg,s) ### normalize to the mean
                matrix.append(s)
                row_header.append(t[0])
            x+=1
            
    time_diff = str(round(time.time()-start_time,1))
    try:
        print '%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error
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
    
    k = 0
    #import random
    column_header=[]; group_db={}; color_db={}
    color_list = ['r', 'b', 'y', 'g', 'w', 'k', 'm']
    #color_list=[]
    #color_template = [1,1,1,0,0,0,0.5,0.5,0.5,0.25,0.25,0.25,0.75,0.75,0.75]
    
    for i in t:
        repls = {'.txt' : '', '.bed' : '', '.tab' : ''}
        i=reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        if ':' in i:
            group,j = string.split(i,':')
            try: color = color_db[group]
            except Exception:
                try: color_db[group] = color_list[k]
                except Exception:
                    ### If not listed in the standard color set add a new random color
                    rgb = tuple(scipy.rand(3)) ### random color
                    #rgb = tuple(random.sample(color_template,3)) ### custom alternative method
                    color_list.append(rgb)
                    color_db[group] = color_list[k]
                color = color_db[group]
                k+=1
            group_db[i] = group, color
        column_header.append(i)
    return group_db, column_header

def PrincipalComponentAnalysis(matrix, column_header, row_header, dataset_name, group_db, display=False):
    print "Performing Principal Component Analysis..."
    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
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
    M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    #if len(row_header)>20000:
    #print '....Using eigenvectors of the real symmetric square matrix for efficiency...'
    #[latent,coeff] = scipy.sparse.linalg.eigsh(cov(M))
    [latent,coeff] = linalg.eig(cov(M))
    scores = dot(coeff.T,M) # projection of the data in the new space
    #scores=mlab.PCA(scores)
    pylab.figure()
    pylab.xlabel('Principal Component 1')
    pylab.ylabel('Principal Component 2')
    pylab.title('Principal Component Analysis - '+dataset_name)
    
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)

    i=0
    group_names={}
    for x in scores[0]:
        ### Add the text labels for each
        sample_name = column_header[i]
        try:
            ### Get group name and color information
            group_name,color = group_db[sample_name]
            if group_name not in group_names:
                label = group_name ### Only add once for each group
            else: label = None
            group_names[group_name] = color
        except Exception:
            color = 'r'; label=None
        pylab.plot(scores[0][i],scores[1][i],color=color,marker='o',markersize=7,label=label)
        pylab.text(scores[0][i],scores[1][i],sample_name,fontsize=8)
        i+=1

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
        try: pylab.show()
        except Exception: None ### when run in headless mode
    #pylab.close()

def PCA3D(matrix, column_header, row_header, dataset_name, group_db, display=False, showLabels=True):
    from numpy import mean,cov,double,cumsum,dot,linalg,array,rank
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    start = time.time()
    M = (matrix-mean(matrix.T,axis=1)).T # subtract the mean (along columns)
    [latent,coeff] = linalg.eig(cov(M))
    #[latent,coeff] = scipy.sparse.linalg.eigsh(cov(M))
    scores = dot(coeff.T,M) # projection of the data in the new space
    #scores = mlab.PCA(scores)
    end = time.time()
    print 'PCA completed in', end-start, 'seconds.'
    ### Hide the axis number labels
    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])

    """
    ax.set_xticks([]) ### Hides ticks
    ax.set_yticks([])
    ax.set_zticks([])    
    
    ax.set_xlabel('Component 1')
    ax.set_ylabel('Component 2')
    ax.set_zlabel('Component 3')
    """
    #pylab.title('Principal Component Analysis\n'+dataset_name)
    """
    pylab.figure()
    pylab.xlabel('Principal Component 1')
    pylab.ylabel('Principal Component 2')

    """
    axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    pylab.axis(axes)
    
    i=0
    group_names={}
    for x in scores[0]:
        ### Add the text labels for each
        sample_name = column_header[i]
        try:
            ### Get group name and color information
            group_name,color = group_db[sample_name]
            if group_name not in group_names:
                label = group_name ### Only add once for each group
            else: label = None
            group_names[group_name] = color
        except Exception:
            color = 'r'; label=None

        ax.plot([scores[0][i]],[scores[1][i]],[scores[2][i]],color=color,marker='o',markersize=9,label=label,markeredgewidth=0.1) #markeredgecolor=color
        if showLabels:
            ax.text(scores[0][i],scores[1][i],scores[2][i], '   '+sample_name,fontsize=8)
        i+=1

    #pylab.legend(loc="upper left", prop={'size': 10})
    ax.legend()
    filename = 'Clustering-%s-PCA.pdf' % dataset_name
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename) #dpi=200
    graphic_link.append(['Principal Component Analysis',root_dir+filename])
    if display:
        print 'Exporting:',filename
        try: pylab.show()
        except Exception: None ### when run in headless mode
        pylab.close()
        
def getAxes(scores):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    x_range = max(scores[0])-min(scores[0])
    y_range = max(scores[1])-min(scores[1])
    x_axis_min = min(scores[0])-(x_range/1.5)
    x_axis_max = max(scores[0])+(x_range/1.5)
    y_axis_min = min(scores[1])-(y_range/5)
    y_axis_max = max(scores[1])+(y_range/5)
    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]
    
def Kmeans(features, column_header, row_header):
    #http://www.janeriksolem.net/2009/04/clustering-using-scipys-k-means.html
    class1 = array(random.standard_normal((100,2))) + array([5,5])
    class2 = 1.5 * array(random.standard_normal((100,2)))
    #features = vstack((class1,class2))
    #This generates two normally distributed classes in two dimensions. To try and cluster the points, run k-means with k=2 like this.
    centroids,variance = kmeans(features,2)
    #The variance is returned but we don't really need it since the SciPy implementation computes several runs (default is 20) and selects the one with smallest variance for us. Now you can check where each data point is assigned using the vector quantization function in the SciPy package.
    code,distance = vq(features,centroids)
    
    #By checking the value of code we can see if there are any incorrect assignments. To visualize, we can plot the points and the final centroids.
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
                              color_gradient, display=False, contrast=None):
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
        if allowLargeClusters: maxSize = 200000
        else: maxSize = 7000
    except Exception: maxSize = 7000
    
    run = False
    print 'max allowed cluster size:',maxSize
    if len(matrix)>0 and len(matrix)<maxSize:
        if len(matrix)>5000:
            row_metric = 'euclidean'
        try:
            ### Default for display is False, when set to True, Pylab will render the image
            heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name, display=display,contrast=contrast)
            run = True
        except Exception:
            try:
                pylab.clf()
                pylab.close() ### May result in TK associated errors later on
                import gc
                gc.collect()
            except Exception: None
            if len(matrix)<5000:
                print 'Error using %s ... trying euclidean instead' % row_metric
                row_metric = 'euclidean' ### cityblock
            else:
                print 'Error with hierarchical clustering... only clustering arrays'
                row_method = None ### Skip gene clustering
            try:
                heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name, display=display,contrast=contrast)
                run = True
            except Exception:
                run = traceback.format_exc()
                print run
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

def runHCexplicit(filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient, extra_params, display=True, contrast=None):
    """ Explicit method for hiearchical clustering with defaults defined by the user (see below function) """
    
    global root_dir
    global inputFilename
    global graphic_link
    global allowLargeClusters
    global GroupDB
    allowLargeClusters = True
    Normalize = False ### if dealing with expression values not folds
    
    graphic_link=graphics ### Store all locations of pngs
    inputFilename = filename ### Used when calling R
    filterIDs = False
    
    try:
        ### Specific additional optional parameters for filtering
        transpose = extra_params.Transpose()
        PathwayFilter = extra_params.PathwaySelect()
        GeneSet = extra_params.GeneSet()
        OntologyID = extra_params.OntologyID()
        Normalize = extra_params.Normalize()
        filterIDs = True
        species = extra_params.Species()
        platform = extra_params.Platform()
        vendor = extra_params.Vendor()
        newInput = findParentDir(inputFilename)+'/GeneSetClustering/'+findFilename(inputFilename)
        PathwayFilter = verifyPathwayName(PathwayFilter,GeneSet,OntologyID)
        targetGene = extra_params.GeneSelection() ### Select a gene or ID to get the top correlating genes
        getGeneCorrelations = extra_params.GetGeneCorrelations() ### Select a gene or ID to get the top correlating genes
        filterByPathways = extra_params.FilterByPathways()
    except Exception,e:
        transpose = extra_params
        Normalize=False
        
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
    
    matrix, column_header, row_header, dataset_name, group_db = importData(filename,Normalize=Normalize)
    GroupDB = group_db
    #print len(matrix),;print len(column_header),;print len(row_header)
    
    if filterIDs:
        transpose_update = True ### Since you can filterByPathways and getGeneCorrelations, only transpose once
        if filterByPathways: ### Restrict analyses to only a single pathway/gene-set/ontology term
            inputFilename = string.replace(newInput,'.txt','_'+PathwayFilter+'.txt') ### update the pathway reference for HOPACH
            matrix,row_header,column_header = filterByPathway(matrix,row_header,column_header,species,platform,vendor,GeneSet,PathwayFilter,OntologyID,transpose)
            dataset_name += '-'+PathwayFilter
            transpose_update = False
        if getGeneCorrelations: ###Restrict analyses to only genes that correlate with the target gene of interest
            if transpose and transpose_status == False: transpose_update = False ### If filterByPathways selected
            elif transpose and transpose_update: transpose_update = True ### If filterByPathways not selected
            else: transpose_update = False ### If transpose == False
            inputFilename = string.replace(newInput,'.txt','_'+targetGene+'.txt') ### update the pathway reference for HOPACH
            try: matrix,row_header,column_header,row_method = getAllCorrelatedGenes(matrix,row_header,column_header,species,platform,vendor,targetGene,transpose_update)
            except Exception:
                print targetGene, 'not found in input expression file. Exiting. \n\n'
                badExit
            dataset_name += '-'+targetGene
    else:
        if transpose: ### Transpose the data matrix
            print 'Transposing the data matrix'
            matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
            column_header, row_header = row_header, column_header
    #print len(matrix),;print len(column_header),;print len(row_header)
    
    if len(column_header)>1000 or len(row_header)>1000:
        print 'Performing hierarchical clustering (please be patient)...'
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, display=display,contrast=contrast)
    #debugPylab()
    return graphic_link

def debugPylab():
    pylab.figure()
    pylab.close()
    pylab.figure()

def verifyPathwayName(PathwayFilter,GeneSet,OntologyID):
    ### If the user supplied an Ontology ID rather than a Ontology term name, lookup the term name and return this as the PathwayFilter
    if len(OntologyID)>0:
        PathwayFilter = lookupOntologyID(GeneSet,OntologyID,'ID')
    return PathwayFilter

def filterByPathway(matrix,row_header,column_header,species,platform,vendor,GeneSet,PathwayFilter,OntologyID,transpose):
    ### Filter all the matrix and header entries for IDs in the selected pathway
    import gene_associations
    import export
    exportData = export.ExportFile(inputFilename)
    
    matrix2=[]; row_header2=[]
    if 'Ontology' in GeneSet: directory = 'nested'
    else: directory = 'gene-mapp'
    associated_IDs = gene_associations.simpleGenePathwayImport(species,GeneSet,PathwayFilter,OntologyID,directory)
    gene_annotations = gene_associations.importGeneData(species,'Ensembl')

    if platform == "3'array":
        ### IDs thus won't be Ensembl - need to translate
        try:
            #ens_to_array = gene_associations.getGeneToUidNoExon(species,'Ensembl-'+vendor); print vendor, 'IDs imported...'
            array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,associated_IDs)
        except Exception: print platform, vendor, 'not found!!! Exiting method'; badExit
        #array_to_ens = gene_associations.swapKeyValues(ens_to_array)
        
    i=0
    original_rows={} ### Don't add the same original ID twice if it associates with different Ensembl IDs
    for row_id in row_header:
        original_id = row_id
        if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename:
            try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
            except Exception:
                if i==0:
                    import gene_associations
                    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                try: symbol = gene_to_symbol[row_id][0]
                except Exception: None
            original_id = row_id
        if platform == "3'array":
            try:
                try: row_ids = array_to_ens[row_id]
                except Exception: row_ids = [row_id]
            except Exception:
                print len(array_to_ens[row_id]);kill
            #print row_id, row_ids
        else:
            row_ids = [row_id]
        for row_id in row_ids:
            if row_id in associated_IDs:
                if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename:
                    row_id = original_id+' '+symbol
                else:
                    try: row_id = gene_annotations[row_id].Symbol()
                    except Exception: None ### If non-Ensembl data
                if original_id not in original_rows: ### Don't add the same ID twice if associated with mult. Ensembls
                    matrix2.append(matrix[i])
                    row_header2.append(row_id)
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

def getAllCorrelatedGenes(matrix,row_header,column_header,species,platform,vendor,targetGene,transpose):
    ### Filter all the matrix and header entries for IDs in the selected targetGene
    import export
    exportData = export.ExportFile(inputFilename)
    
    i=0
    original_rows={} ### Don't add the same original ID twice if it associates with different Ensembl IDs
    for row_id in row_header:
        original_id = row_id
        symbol = 'NA'
        if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename:
            try: row_id,symbol = string.split(row_id,' ')[:2] ### standard ID convention is ID space symbol
            except Exception: row_id, symbol = row_id, row_id
            original_id = row_id
        if row_id == targetGene or symbol == targetGene:
            targetGeneValues = matrix[i] ### Values for the row
            break
        i+=1
    i=0

    import markerFinder
    rho_results = markerFinder.simpleScipyPearson(matrix,targetGeneValues)
    
    matrix2=[]
    row_header2=[]
    for (rho,ind) in rho_results[:551]: ### Get the top-50 correlated plus the gene of interest
        if rho>0.6:
            matrix2.append(matrix[ind])
            row_header2.append(row_header[ind])

    matrix2.reverse() ### Display from top-to-bottom rather than bottom-to-top (this is how the clusters are currently ordered in the heatmap)
    row_header2.reverse()

    if transpose:
        matrix2 = map(numpy.array, zip(*matrix2)) ### coverts these to tuples
        column_header, row_header2 = row_header2, column_header

    exportData.write(string.join(['UID']+column_header,'\t')+'\n') ### title row export
    i=0
    for row_id in row_header2:
        exportData.write(string.join([row_id]+map(str,matrix2[i]),'\t')+'\n') ### export values
        i+=1
        
    print len(row_header2), 'top-correlated IDs'
    exportData.close()
    row_method = None ### don't cluster the rows (row_method)
    return matrix2,row_header2,column_header,row_method

def runHCOnly(filename,graphics):
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
        
    row_method = 'weighted'
    column_method = 'average'
    row_metric = 'cosine'
    column_metric = 'cosine'
    if 'Lineage' in filename or 'Elite' in filename:
        color_gradient = 'red_white_blue'
    else:
        color_gradient = 'yellow_black_blue'
        color_gradient = 'red_black_green'
    
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    GroupDB = group_db
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, display=False)
    return graphic_link

def runPCAonly(filename,graphics,transpose,showLabels=True,plotType='3D'):
    global root_dir
    global graphic_link
    graphic_link=graphics ### Store all locations of pngs
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')
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
    if plotType == '3D':
        PCA3D(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=True, showLabels=showLabels)
    else:
        PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=True)
        
def outputClusters(filenames,graphics):
    """ Peforms PCA and Hiearchical clustering on exported log-folds from AltAnalyze """
    
    global root_dir
    global graphic_link
    global inputFilename
    global GroupDB
    global allowLargeClusters
    allowLargeClusters=False

    graphic_link=graphics ### Store all locations of pngs
    filename = filenames[0] ### This is the file to clustering with "significant" gene changes
    inputFilename = filename ### Used when calling R
    
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionOutput/Clustering','DataPlots')

    ### Transpose matrix and build PCA
    original = importData(filename)
    matrix, column_header, row_header, dataset_name, group_db = original
    matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
    column_header, row_header = row_header, column_header
    if len(row_header)<7000 and len(column_header)<7000:
        PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db)
    else:
        print 'SKIPPING PCA!!! - Your dataset file is over the recommended size limit for clustering (>7000 rows). Please cluster later using "Additional Analyses".'

    row_method = 'weighted'
    column_method = 'average'
    row_metric = 'cosine'
    column_metric = 'cosine'
    color_gradient = 'red_white_blue'
    color_gradient = 'red_black_green'
    
    ### Generate Significant Gene HeatMap
    matrix, column_header, row_header, dataset_name, group_db = original
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient)
    
    ### Generate Outlier and other Significant Gene HeatMap
    for filename in filenames[1:]:
        inputFilename = filename
        matrix, column_header, row_header, dataset_name, group_db = importData(filename)
        GroupDB = group_db
        try:
            runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient)
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
    import export
    
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
    import export
    
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

            layout = "kk"
            visual_style = {}
            #visual_style["layout"] = layout #The default is auto, which selects a layout algorithm automatically based on the size and connectedness of the graph
            visual_style["margin"] = 50 ### white-space around the network (see vertex size)
            visual_style["bbox"] = canvas_size
            igraph.plot(gr,output_filename, **visual_style)
            
            output_filename = '%s.png' % filePath[:-4]
            if vertices <15: gr,visual_style = increasePlotSize(gr,visual_style)
            igraph.plot(gr,output_filename, **visual_style)
            #surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    return output_filename

def correctedFilePath(filePath,canvas_scaler):
    """ Move this file to it's own network directory for GO-Elite """
    if 'ORA_pruned' in filePath:
        filePath = string.replace(filePath,'CompleteResults/ORA_pruned','networks')
        try: os.mkdir(findParentDir(filePath))
        except Exception: None
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
            print e
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
    
if __name__ == '__main__':
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
    import UI
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