### hierarchical_clustering.py
#Copyright 2005-2012 J. David Gladstone Institutes, San Francisco California
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

#################
### Imports an tab-delimited expression matrix and produces and hierarchically clustered heatmap
#################

import matplotlib.pyplot as pylab
import matplotlib as mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy
import string
import time
import sys, os
import getopt
import warnings
import fastcluster as fc

################# Perform the hierarchical clustering #################

def heatmap(x, row_header, column_header, row_method,
            column_method, row_metric, column_metric,
            color_gradient, filename, graphics=True,
            contrast=2.5):
    
    print "\nPerforming hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)
        
    ### Define the color gradient to use based on the provided name
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
    
    ### Scale the max and min colors so that 0 is white/black
    vmin=x.min()
    vmax=x.max()
    vmax = max([vmax,abs(vmin)])
    vmin = vmax*-1
    scaling_factor = float(contrast)
    norm = mpl.colors.Normalize(vmin/scaling_factor, vmax/scaling_factor)

    ### Scale the Matplotlib window size
    default_window_hight = 8.5
    default_window_width = 12
    fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show
        
    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 = 
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
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.02,0.938,0.17,0.025] ### Last one controls the hight [0.07,0.88,0.18,0.076]

    # Compute and plot top dendrogram

    if column_method != None:
        start_time = time.time()
        d2 = dist.pdist(x.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
        Y2 = fc.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z2 = sch.dendrogram(Y2)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
    # Compute and plot left dendrogram.
    if row_method != None:
        start_time = time.time()
        d1 = dist.pdist(x)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
        Y1 = fc.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='left') ### This is where plotting occurs - orientation 'right' in old matplotlib
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

    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, '  '+row_header[idx1[i]])
            new_row_header.append(row_header[idx1[i]])
        else:
            if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
                axm.text(x.shape[1]-0.5, i, '  '+row_header[i]) ### When not clustering rows
            new_row_header.append(row_header[i])
    for i in range(x.shape[1]):
        if column_method != None:
            axm.text(i, -0.9, ' '+column_header[idx2[i]], rotation=270, verticalalignment="top") # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            axm.text(i, -0.9, ' '+column_header[i], rotation=270, verticalalignment="top")
            new_column_header.append(column_header[i])

    exportFlatClusterData(filename[:-4]+'-clustered.txt', new_row_header,new_column_header,xt,ind1,ind2)
    
    if graphics:
        # Plot colside colors
        # axc --> axes for column side colorbar
        if column_method != None:
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
            cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            dc = numpy.array(ind2, dtype=int)
            dc.shape = (1,len(ind2)) 
            im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
            axc.set_xticks([]) ### Hides ticks
            axc.set_yticks([])
        
        # Plot rowside colors
        # axr --> axes for row side colorbar
        if row_method != None:
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
            dr = numpy.array(ind1, dtype=int)
            dr.shape = (len(ind1),1)
            #print ind1, len(ind1)
            cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
            axr.set_xticks([]) ### Hides ticks
            axr.set_yticks([])
    
        # Plot color legend
        axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
        axcb.set_title("Normalized Expression")
    
        if '/' in filename:
            dataset_name = string.split(filename,'/')[-1][:-4]
            root_dir = string.join(string.split(filename,'/')[:-1],'/')+'/'
        else:
            dataset_name = string.split(filename,'\\')[-1][:-4]
            root_dir = string.join(string.split(filename,'\\')[:-1],'\\')+'\\'
        filename = root_dir+'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name,column_metric,row_metric)
        cb.set_label("Differential Expression (log2 fold)")

        ### Render the graphic
        if len(row_header)>50 or len(column_header)>50:
            pylab.rcParams['font.size'] = 5
        else:
            pylab.rcParams['font.size'] = 8
    
        pylab.savefig(filename)
        print 'Exporting:',filename
        filename = filename[:-3]+'png'
        pylab.savefig(filename, dpi=100) #,dpi=200
        pylab.show()

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
    
################# Export the flat cluster data #################

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    filename = string.replace(filename,'.pdf','.txt')
    export_text = open(filename,'w')
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
    export_cdt = open(filename,'w')
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

################# Create Custom Color Gradients #################
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

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

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
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

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
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
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

################# General data import methods #################

def importData(filename,normalize=False):
    start_time = time.time()
    matrix=[]
    row_header=[]
    first_row=True

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
        
    for line in open(filename,'rU').xreadlines():         
        t = string.split(line[:-1],'\t') ### remove end-of-line character - file is tab-delimited
        if first_row:
            column_header = t[1:]
            first_row=False
        else:
            try:
                s = map(float,t[1:])
            except: ### If missing values
                s=[]
                for v in t[1:]:
                    try: s.append(float(v))
                    except: s.append(0)
            
            if normalize!=False:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
                    avg = numpy.median(s)
                s = map(lambda x: x-avg,s) ### normalize to the mean
                matrix.append(s)
                row_header.append(t[0])
                
    time_diff = str(round(time.time()-start_time,1))
    try:
        print '\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error
    return numpy.array(matrix), column_header, row_header
  
if __name__ == '__main__':
    
    ################  Default Methods ################
    row_method = 'ward'
    column_method = 'ward'
    row_metric = 'euclidean'
    column_metric = 'euclidean'
    color_gradient = 'red_white_blue'
    graphics = True
    normalize = True
    contrast = 2.5
    
    """ Running with cosine or other distance metrics can often produce negative Z scores
        during clustering, so adjustments to the clustering may be required.
        
    see: http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    see: http://docs.scipy.org/doc/scipy/reference/spatial.distance.htm  
    color_gradient = red_white_blue|red_black_sky|red_black_blue|red_black_green|yellow_black_blue|green_white_purple'
    """
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print "Example:"
        print "python hierarchical_clustering_basic.py --i FinalMarkerHeatmap.txt --row_method ward --column_method ward --row_metric cosine --column_metric cosine --normalize True --color_gradient 'yellow_black_blue' --contrast 5 --graphics True"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','row_method=','column_method=',
                                                    'row_metric=','column_metric=','color_gradient=',
                                                    'graphics=','normalize=','contrast='])
        for opt, arg in options:
            if opt == '--i': filename=arg
            elif opt == '--row_method':
                row_method=arg
                if row_method=='None': row_method = None
            elif opt == '--column_method':
                column_method=arg
                if column_method=='None': column_method = None
            elif opt == '--row_metric': row_metric=arg
            elif opt == '--column_metric': column_metric=arg
            elif opt == '--color_gradient': color_gradient=arg
            elif opt == '--contrast': contrast=float(contrast)
            elif opt == '--graphics':
                graphics = arg
                if graphics == 'no' or graphics == 'false' or graphics == 'False' or grapics == 'FALSE':
                    graphics = False
            elif opt == '--normalize':
                normalize = arg
                if normalize == 'no' or normalize == 'false' or normalize == 'False' or normalize == 'FALSE':
                    normalize = False

            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
    matrix, column_header, row_header = importData(filename,normalize=normalize)

    if len(matrix)>0:
        try:
            heatmap(matrix, row_header, column_header, row_method, column_method, row_metric,
                    column_metric, color_gradient, filename, graphics=graphics, contrast=contrast)
        except Exception:
            print 'Error using %s ... trying euclidean instead' % row_metric
            row_metric = 'euclidean'
            try:
                heatmap(matrix, row_header, column_header, row_method, column_method, row_metric,
                        column_metric, color_gradient, filename, graphics=graphics, contrast=contrast)
            except IOError:
                print 'Error with clustering encountered'
