import sys, string, os, copy
sys.path.insert(1, os.path.join(sys.path[0], '..'))
command_args = string.join(sys.argv, ' ')
if len(sys.argv[1:]) > 0 and '--' in command_args:
    commandLine = True
else:
    commandLine = False
display_label_names = True
import traceback
try:
    import math, warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        import matplotlib
        if commandLine and 'linux' in sys.platform:
            try:
                matplotlib.use('Agg')
            except Exception:
                pass
            else:
                try:
                    matplotlib.rcParams['backend'] = 'Agg'
                except Exception:
                    pass

        else:
            try:
                import matplotlib.backends.backend_tkagg
                matplotlib.use('TkAgg')
            except Exception:
                pass
            else:
                try:
                    matplotlib.rcParams['backend'] = 'TkAgg'
                except Exception:
                    pass
                else:
                    try:
                        import matplotlib.pyplot as pylab, matplotlib.colors as mc, matplotlib.mlab as mlab, matplotlib.ticker as tic
                        from matplotlib.patches import Circle
                        import mpl_toolkits
                        from mpl_toolkits import mplot3d
                        from mpl_toolkits.mplot3d import Axes3D
                        try:
                            from matplotlib.cbook import _string_to_bool
                        except:
                            pass

                        matplotlib.rcParams['axes.linewidth'] = 0.5
                        matplotlib.rcParams['pdf.fonttype'] = 42
                        matplotlib.rcParams['font.family'] = 'sans-serif'
                        matplotlib.rcParams['font.sans-serif'] = 'Arial'
                        matplotlib.rcParams['figure.facecolor'] = 'white'
                    except Exception:
                        print traceback.format_exc()
                        print 'Matplotlib support not enabled'
                    else:
                        import scipy
                        try:
                            from scipy.sparse.csgraph import _validation
                        except Exception:
                            pass
                        else:
                            try:
                                from scipy.linalg import svd
                                import scipy.special._ufuncs_cxx
                                from scipy.spatial import _voronoi
                                from scipy.spatial import _spherical_voronoi
                                from scipy.spatial import qhull
                                import scipy._lib.messagestream
                            except Exception:
                                pass
                            else:
                                import scipy.cluster.hierarchy as sch, scipy.spatial.distance as dist
                                try:
                                    import numpy
                                    np = numpy
                                except:
                                    print 'Numpy import error...'
                                    print traceback.format_exc()
                                else:
                                    if 'darwin' in sys.platform:
                                        try:
                                            import umap
                                        except:
                                            pass

                                    try:
                                        from cairo import ImageSurface
                                    except:
                                        pass

                                    try:
                                        import igraph.vendor.texttable
                                    except ImportError:
                                        pass

                                try:
                                    from sklearn.decomposition import PCA, FastICA
                                except Exception:
                                    pass

                            try:
                                from sklearn.neighbors import quad_tree
                            except:
                                pass

                        try:
                            import sklearn.utils.sparsetools._graph_validation
                        except:
                            pass

                    try:
                        import sklearn.utils.weight_vector
                    except:
                        pass

                from sklearn.neighbors import *
                from sklearn.manifold.t_sne import *
                from sklearn.tree import *
                from sklearn.tree import _utils
                from sklearn.manifold.t_sne import _utils
                from sklearn.manifold import TSNE
                from sklearn.neighbors import NearestNeighbors
                import sklearn.linear_model.sgd_fast, sklearn.utils.lgamma
                try:
                    import scipy.special.cython_special
                except:
                    pass

            import sklearn.neighbors.typedefs, sklearn.neighbors.ball_tree
            try:
                import numba, numba.config, llvmlite
                from llvmlite import binding
                from llvmlite.binding import *
                from llvmlite.binding import ffi
                from llvmlite.binding import dylib
            except:
                pass

except Exception:
    print traceback.format_exc()
else:
    try:
        import numpy
    except:
        pass
    else:
        import time, unique
        from stats_scripts import statistics
        import os, export, webbrowser, warnings, UI
        use_default_colors = False
        try:
            warnings.simplefilter('ignore', numpy.ComplexWarning)
            warnings.simplefilter('ignore', DeprecationWarning)
        except Exception:
            None

    from visualization_scripts import WikiPathways_webservice
    try:
        import fastcluster as fc
    except Exception:
        try:
            fc = sch
        except Exception:
            print 'Scipy support not present...'

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax = x.max()
    vmin = x.min()
    if vmax < 0 and vmin < 0:
        direction = 'negative'
    elif vmax > 0 and vmin > 0:
        direction = 'positive'
    else:
        direction = 'both'
    if direction == 'both':
        vmax = max([vmax, abs(vmin)])
        vmin = -1 * vmax
        return (
         vmax, vmin)
    else:
        return (
         vmax, vmin)


def heatmap(x, row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name, display=False, contrast=None, allowAxisCompression=True, Normalize=True, PriorColumnClusters=None, PriorRowClusters=None):
    print 'Performing hieararchical clustering using %s for columns and %s for rows' % (column_metric, row_metric)
    show_color_bars = True
    try:
        ExportCorreleationMatrix = exportCorreleationMatrix
    except Exception:
        ExportCorreleationMatrix = False
    else:
        try:
            os.mkdir(root_dir)
        except Exception:
            None

        if display == False:
            pylab.figure()
        if row_method == 'hopach' or column_method == 'hopach':
            pass
        if row_method == 'hopach' or column_method == 'hopach':
            try:
                import R_interface
                if row_method == 'hopach' and column_method == 'hopach':
                    cluster_method = 'both'
                elif row_method == 'hopach':
                    cluster_method = 'gene'
                else:
                    cluster_method = 'array'
                if row_metric == 'cosine':
                    metric_gene = 'euclid'
                elif row_metric == 'euclidean':
                    metric_gene = 'cosangle'
                elif row_metric == 'correlation':
                    metric_gene = 'cor'
                else:
                    metric_gene = 'cosangle'
                if column_metric == 'cosine':
                    metric_array = 'euclid'
                elif column_metric == 'euclidean':
                    metric_array = 'cosangle'
                elif column_metric == 'correlation':
                    metric_array = 'cor'
                else:
                    metric_array = 'euclid'
                newFilename, Z1, Z2 = R_interface.remoteHopach(inputFilename, cluster_method, metric_gene, metric_array)
                if newFilename != inputFilename:
                    try:
                        matrix, column_header, row_header, dataset_name, group_db = importData(newFilename, Normalize=normalize, reverseOrder=False)
                    except Exception:
                        matrix, column_header, row_header, dataset_name, group_db = importData(newFilename)

                    x = numpy.array(matrix)
            except Exception:
                row_method = 'average'
                column_method = 'average'
                print traceback.format_exc()
                print 'hopach failed... continue with an alternative method'

        skipClustering = False
        try:
            if len(PriorColumnClusters) > 0 and len(PriorRowClusters) > 0 and row_method == None and column_method == None:
                print 'Prior generated clusters being used rather re-clustering'
                if len(PriorColumnClusters) > 0:
                    Z1 = {}
                    Z2 = {}
                    Z1['level'] = PriorRowClusters
                    Z1['level'].reverse()
                    Z2['level'] = PriorColumnClusters
                    Z1['leaves'] = range(0, len(row_header))
                    Z2['leaves'] = range(0, len(column_header))
                    skipClustering = True
                    row_method = None
                    column_method = None
                    row_method = 'hopach'
                    column_method = 'hopach'
        except Exception as e:
            pass

    n = len(x[0])
    m = len(x)
    if color_gradient == 'red_white_blue':
        cmap = pylab.cm.bwr
    if color_gradient == 'red_black_sky':
        cmap = RedBlackSkyBlue()
    if color_gradient == 'red_black_blue':
        cmap = RedBlackBlue()
    if color_gradient == 'red_black_green':
        cmap = RedBlackGreen()
    if color_gradient == 'yellow_black_blue':
        cmap = YellowBlackBlue()
    if color_gradient == 'black_yellow_blue':
        cmap = BlackYellowBlue()
    if color_gradient == 'seismic':
        cmap = pylab.cm.seismic
    if color_gradient == 'green_white_purple':
        cmap = pylab.cm.PiYG_r
    if color_gradient == 'coolwarm':
        cmap = pylab.cm.coolwarm
    if color_gradient == 'Greys':
        cmap = pylab.cm.Greys
    if color_gradient == 'yellow_orange_red':
        cmap = pylab.cm.YlOrRd
    if color_gradient == 'Spectral':
        cmap = pylab.cm.Spectral_r
    vmin = x.min()
    vmax = x.max()
    vmax = max([vmax, abs(vmin)])
    if Normalize != False:
        vmin = vmax * -1
    else:
        if 'Clustering-Zscores-' in dataset_name:
            vmin = vmax * -1
        elif vmin < 0 and vmax > 0 and Normalize == False:
            vmin = vmax * -1
        else:
            default_window_hight = 8.5
            default_window_width = 12
            if len(column_header) > 80:
                default_window_width = 14
            if len(column_header) > 100:
                default_window_width = 16
            if contrast == None:
                scaling_factor = 2.5
            else:
                try:
                    scaling_factor = float(contrast)
                except Exception:
                    scaling_factor = 2.5

        norm = matplotlib.colors.Normalize(vmin / scaling_factor, vmax / scaling_factor)
        fig = pylab.figure()
        fig.set_figwidth(12)
        fig.set_figheight(7)
        fig.patch.set_facecolor('white')
        pylab.rcParams['font.size'] = 7.5
        if show_color_bars == False:
            color_bar_w = 1e-06
        else:
            color_bar_w = 0.0125
        bigSampleDendrogram = True
        if bigSampleDendrogram == True and row_method == None and column_method != None and allowAxisCompression == True:
            dg2 = 0.3
            dg1 = 0.43
        else:
            dg2 = 0.1
            dg1 = 0.63
        try:
            if EliteGeneSets != [''] and EliteGeneSets != []:
                matrix_horiz_pos = 0.27
            elif skipClustering:
                if len(row_header) < 100:
                    matrix_horiz_pos = 0.2
                else:
                    matrix_horiz_pos = 0.27
            else:
                matrix_horiz_pos = 0.14
        except Exception:
            matrix_horiz_pos = 0.14

    if len(column_header) < 50:
        matrix_horiz_pos += 0.1
    ax1_x, ax1_y, ax1_w, ax1_h = [
     0.05, 0.235, matrix_horiz_pos, dg1]
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004
    axr_x, axr_y, axr_w, axr_h = [
     0.31, 0.1, color_bar_w - 0.002, 0.6]
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y
    axr_h = ax1_h
    width_between_axr_axm = 0.004
    axc_x, axc_y, axc_w, axc_h = [
     0.5, 0.63, 0.5, color_bar_w]
    if len(column_header) < 50:
        axc_w = 0.3
        if len(column_header) < 20:
            axc_w = 0.2
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004
    axm_x, axm_y, axm_w, axm_h = [
     0.4, 0.9, 2.5, 0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y
    axm_h = ax1_h
    axm_w = axc_w
    ax2_x, ax2_y, ax2_w, ax2_h = [
     0.3, 0.72, 0.6, dg2]
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w
    axcb_x, axcb_y, axcb_w, axcb_h = [
     0.02, 0.938, 0.17, 0.025]
    axcc_x, axcc_y, axcc_w, axcc_h = [
     0.02, 0.12, 0.17, 0.025]
    if column_method == 'hopach':
        ind2 = numpy.array(Z2['level'])
    elif column_method != None:
        start_time = time.time()
        d2 = dist.pdist(x.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
        if ExportCorreleationMatrix:
            new_matrix = []
            for i in D2:
                log2_data = map(inverseDist, i)
                avg = statistics.avg(log2_data)
                log2_norm = map(lambda x: x - avg, log2_data)
                new_matrix.append(log2_norm)

            x = numpy.array(new_matrix)
            row_header = column_header
        Y2 = fc.linkage(D2, method=column_method, metric=column_metric)
        try:
            Z2 = sch.dendrogram(Y2)
        except Exception:
            if column_method == 'average':
                column_metric = 'euclidean'
            else:
                column_method = 'average'
            Y2 = fc.linkage(D2, method=column_method, metric=column_metric)
            Z2 = sch.dendrogram(Y2)

        ind2 = sch.fcluster(Y2, 0.7 * max(Y2[:, 2]), 'distance')
        ax2.set_xticks([])
        ax2.set_yticks([])
        time_diff = str(round(time.time() - start_time, 1))
        print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = [
         'NA'] * len(column_header)
    if row_method == 'hopach':
        ind1 = numpy.array(Z1['level'])
    elif row_method != None:
        start_time = time.time()
        d1 = dist.pdist(x)
        D1 = dist.squareform(d1)
        Y1 = fc.linkage(D1, method=row_method, metric=row_metric)
        no_plot = False
        try:
            if runGOElite:
                no_plot = True
            elif len(PriorColumnClusters) > 0 and len(PriorRowClusters) > 0 and row_method == None and column_method == None:
                no_plot = True
            else:
                ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False)
        except Exception:
            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False)
        else:
            try:
                Z1 = sch.dendrogram(Y1, orientation='left', no_plot=no_plot)
            except Exception:
                row_method = 'average'
                try:
                    Y1 = fc.linkage(D1, method=row_method, metric=row_metric)
                    Z1 = sch.dendrogram(Y1, orientation='right', no_plot=no_plot)
                except Exception:
                    row_method = 'ward'
                    Y1 = fc.linkage(D1, method=row_method, metric=row_metric)
                    Z1 = sch.dendrogram(Y1, orientation='right', no_plot=no_plot)

            ind1 = sch.fcluster(Y1, 0.7 * max(Y1[:, 2]), 'distance')
            if ExportCorreleationMatrix:
                Z1 = sch.dendrogram(Y2, orientation='right')
                Y1 = Y2
                d1 = d2
                D1 = D2
                ind1 = ind2
            try:
                ax1.set_xticks([])
                ax1.set_yticks([])
            except Exception:
                pass

        time_diff = str(round(time.time() - start_time, 1))
        print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = [
         'NA'] * len(row_header)
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])
    xt = x
    if column_method != None:
        idx2 = Z2['leaves']
        xt = xt[:, idx2]
        try:
            ind2 = [ ind2[i] for i in idx2 ]
        except Exception:
            column_method = None
            xt = x
            ind2 = ['NA'] * len(column_header)
            ind1 = ['NA'] * len(row_header)

    if row_method != None:
        idx1 = Z1['leaves']
        prior_xt = xt
        xt = xt[idx1, :]
        try:
            ind1 = [ ind1[i] for i in idx1 ]
        except Exception:
            if 'MarkerGenes' in dataset_name:
                ind1 = [
                 'NA'] * len(row_header)
                row_method = None

    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm)
    axm.set_xticks([])
    axm.set_yticks([])
    row_fontsize = 5
    column_fontsize = 5
    column_text_max_len = max(map(lambda x: len(x), column_header))
    if len(row_header) < 75:
        row_fontsize = 6.5
        if len(row_header) < 50:
            row_fontsize = 8
            if len(row_header) < 25:
                row_fontsize = 11
    if len(column_header) < 75:
        column_fontsize = 6.5
        if len(column_header) < 50:
            column_fontsize = 8
            if len(column_header) < 25:
                column_fontsize = 11
                if column_text_max_len < 15:
                    column_fontsize = 15
                elif column_text_max_len > 30:
                    column_fontsize = 6.5
                else:
                    column_fontsize = 10
    try:
        if len(justShowTheseIDs) > 50:
            column_fontsize = 7
        elif len(justShowTheseIDs) > 0:
            column_fontsize = 10
        if len(justShowTheseIDs) > 0:
            additional_symbols = []
            import gene_associations
            try:
                gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', 'Ensembl-Symbol'))
            except Exception:
                gene_to_symbol = {}
                symbol_to_gene = {}

        JustShowTheseIDs = copy.deepcopy(justShowTheseIDs)
    except Exception:
        JustShowTheseIDs = []
    else:
        new_row_header = []
        new_column_header = []
        for i in range(x.shape[0]):
            if row_method != None:
                new_row_header.append(row_header[idx1[i]])
            else:
                new_row_header.append(row_header[i])

        for i in range(x.shape[1]):
            if column_method != None:
                new_column_header.append(column_header[idx2[i]])
            else:
                new_column_header.append(column_header[i])

        dataset_name = string.replace(dataset_name, 'Clustering-', '')
        if '-hierarchical' in dataset_name:
            dataset_name = string.split(dataset_name, '-hierarchical')[0]
        filename = 'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name, column_metric, row_metric)
        if 'MarkerGenes' in dataset_name:
            time_stamp = timestamp()
            filename = string.replace(filename, 'hierarchical', time_stamp)
        elite_dir, cdt_file, SystemCode = exportFlatClusterData(root_dir + filename, root_dir, dataset_name, new_row_header, new_column_header, xt, ind1, ind2, vmax, display)

        def ViewPNG(png_file_dir):
            if os.name == 'nt':
                try:
                    os.startfile('"' + png_file_dir + '"')
                except Exception:
                    os.system('open "' + png_file_dir + '"')

            elif 'darwin' in sys.platform:
                os.system('open "' + png_file_dir + '"')
            elif 'linux' in sys.platform:
                os.system('xdg-open "' + png_file_dir + '"')

        try:
            try:
                temp1 = len(justShowTheseIDs)
                if 'monocle' in justShowTheseIDs and 'guide' not in justShowTheseIDs:
                    import R_interface
                    print 'Running Monocle through R (be patient, this can take 20 minutes+)'
                    R_interface.performMonocleAnalysisFromHeatmap(species, cdt_file[:-3] + 'txt', cdt_file[:-3] + 'txt')
                    png_file_dir = root_dir + '/Monocle/monoclePseudotime.png'
                    ViewPNG(png_file_dir)
            except Exception:
                pass

        except Exception:
            print '...Monocle error:'
            print traceback.format_exc()
        else:
            cluster_elite_terms = {}
            ge_fontsize = 11.5
            top_genes = []
            proceed = True
            try:
                try:
                    if 'guide' in justShowTheseIDs:
                        proceed = False
                except Exception:
                    pass

                if proceed:
                    try:
                        cluster_elite_terms, top_genes = remoteGOElite(elite_dir, SystemCode=SystemCode)
                        if cluster_elite_terms['label-size'] > 40:
                            ge_fontsize = 9.5
                    except Exception:
                        pass

            except Exception:
                pass

        if len(cluster_elite_terms) < 1:
            try:
                elite_dirs = string.split(elite_dir, 'GO-Elite')
                old_elite_dir = elite_dirs[0] + 'GO-Elite' + elite_dirs[(-1)]
                old_elite_dir = string.replace(old_elite_dir, 'ICGS/', '')
                if len(PriorColumnClusters) > 0 and len(PriorRowClusters) > 0 and skipClustering:
                    cluster_elite_terms, top_genes = importGOEliteResults(old_elite_dir)
            except Exception as e:
                pass

        try:
            if len(justShowTheseIDs) < 1 and len(top_genes) > 0 and column_fontsize < 9:
                column_fontsize = 10
            if len(justShowTheseIDs) < 1:
                additional_symbols = []
                import gene_associations
                from import_scripts import OBO_import
                try:
                    gene_to_symbol = gene_associations.getGeneToUid(species, ('hide',
                                                                              'Ensembl-Symbol'))
                except Exception:
                    gene_to_symbol = {}
                    symbol_to_gene = {}

        except Exception:
            pass

    def formatpval(p):
        if '-' in p:
            p1 = p[:1] + p[-4:]
        else:
            p1 = ('{number:.{digits}f}').format(number=float(p), digits=3)
            p1 = str(p1)
        return p1

    new_row_header = []
    new_column_header = []
    ci = 0
    last_cluster = 1
    interval = int(float(string.split(str(len(row_header) / 35.0), '.')[0])) + 1
    increment = interval - 2
    if len(row_header) < 100:
        increment = interval - 1
    label_pos = -0.03 * len(column_header) - 0.8
    alternate = 1
    try:
        if 'top' in justShowTheseIDs:
            justShowTheseIDs.remove('top')
        if 'positive' in justShowTheseIDs:
            justShowTheseIDs.remove('positive')
        if 'amplify' in justShowTheseIDs:
            justShowTheseIDs.remove('amplify')
        if 'IntraCorrelatedOnly' in justShowTheseIDs:
            justShowTheseIDs.remove('IntraCorrelatedOnly')
        if 'GuideOnlyCorrelation' in justShowTheseIDs:
            justShowTheseIDs.remove('GuideOnlyCorrelation')
    except Exception:
        pass
    else:
        for i in range(x.shape[0]):
            if len(row_header) < 40:
                radj = len(row_header) * 0.009
            elif len(row_header) < 70:
                radj = len(row_header) * 0.007
            else:
                radj = len(row_header) * 0.005
            try:
                cluster = str(ind1[i])
            except Exception:
                cluster = 'NA'
            else:
                if cluster == 'NA':
                    new_index = i
                    try:
                        cluster = 'cluster-' + string.split(row_header[new_index], ':')[0]
                    except Exception:
                        pass

                if cluster != last_cluster:
                    ci = 0
                    increment = 0
                color = 'black'
                if row_method != None:
                    try:
                        if row_header[idx1[i]] in JustShowTheseIDs:
                            if len(row_header) > len(justShowTheseIDs):
                                color = 'red'
                        else:
                            color = 'black'
                    except Exception:
                        pass

                    if len(row_header) < 106:
                        if display_label_names == False or 'ticks' in JustShowTheseIDs:
                            if color == 'red':
                                axm.text(x.shape[1] - 0.5, i - radj, '  -', fontsize=row_fontsize, color=color, picker=True)
                            else:
                                axm.text(x.shape[1] - 0.5, i - radj, '  ', fontsize=row_fontsize, color=color, picker=True)
                        else:
                            axm.text(x.shape[1] - 0.5, i - radj, '  ' + row_header[idx1[i]], fontsize=row_fontsize, color=color, picker=True)
                    new_row_header.append(row_header[idx1[i]])
                    new_index = idx1[i]
                else:
                    try:
                        feature_id = row_header[i]
                        if ':' in feature_id:
                            feature_id = string.split(feature_id, ':')[1]
                            if feature_id[(-1)] == ' ':
                                feature_id = feature_id[:-1]
                        if feature_id in JustShowTheseIDs:
                            color = 'red'
                        else:
                            color = 'black'
                    except Exception:
                        pass

                    if len(row_header) < 106:
                        if display_label_names == False or 'ticks' in JustShowTheseIDs:
                            if color == 'red':
                                axm.text(x.shape[1] - 0.5, i - radj, '  -', fontsize=row_fontsize, color=color, picker=True)
                            else:
                                axm.text(x.shape[1] - 0.5, i - radj, '  ', fontsize=row_fontsize, color=color, picker=True)
                        else:
                            axm.text(x.shape[1] - 0.5, i - radj, '  ' + row_header[i], fontsize=row_fontsize, color=color, picker=True)
                    new_row_header.append(row_header[i])
                    new_index = i
                if len(row_header) < 106:
                    pass
                else:
                    feature_id = row_header[new_index]
                    original_feature_id = feature_id
                    if ':' in feature_id:
                        if 'ENS' != feature_id[:3] or 'G0000' in feature_id:
                            feature_id = string.split(feature_id, ':')[1]
                            if feature_id[(-1)] == ' ':
                                feature_id = feature_id[:-1]
                        else:
                            feature_id = string.split(feature_id, ':')[0]
                            try:
                                feature_id = gene_to_symbol[feature_id][0]
                            except Exception:
                                pass

                    if ' ' in feature_id and ('ENS' in feature_id or 'G0000' in feature_id):
                        feature_id = string.split(feature_id, ' ')[1]
                    try:
                        if feature_id in JustShowTheseIDs or original_feature_id in JustShowTheseIDs:
                            color = 'red'
                        else:
                            color = 'black'
                    except Exception:
                        pass

                try:
                    if feature_id in justShowTheseIDs or len(justShowTheseIDs) < 1 and feature_id in top_genes or original_feature_id in justShowTheseIDs:
                        if original_feature_id in justShowTheseIDs:
                            feature_id = original_feature_id
                        if display_label_names and 'ticks' not in justShowTheseIDs:
                            if alternate == 1:
                                buffer = 1.2
                                alternate = 2
                            elif alternate == 2:
                                buffer = 2.4
                                alternate = 3
                            elif alternate == 3:
                                buffer = 3.6
                                alternate = 4
                            elif alternate == 4:
                                buffer = 0
                                alternate = 1
                            axm.text(x.shape[1] - 0.4 + buffer, i - radj, feature_id, fontsize=column_fontsize, color=color, picker=True)
                        else:
                            axm.text(x.shape[1] - 0.5, i - radj, '  -', fontsize=column_fontsize, color=color, picker=True)
                    elif ' ' in row_header[new_index]:
                        symbol = string.split(row_header[new_index], ' ')[(-1)]
                        if len(symbol) > 0:
                            if symbol in justShowTheseIDs:
                                if display_label_names and 'ticks' not in justShowTheseIDs:
                                    axm.text(x.shape[1] - 0.5, i - radj, '  ' + row_header[new_index], fontsize=column_fontsize, color=color, picker=True)
                                else:
                                    axm.text(x.shape[1] - 0.5, i - radj, '  -', fontsize=column_fontsize, color=color, picker=True)
                except Exception:
                    pass
                else:
                    if cluster in cluster_elite_terms or 'cluster-' + cluster in cluster_elite_terms:
                        if 'cluster-' + cluster in cluster_elite_terms:
                            new_cluster_id = 'cluster-' + cluster
                        else:
                            new_cluster_id = cluster
                        if cluster != last_cluster:
                            cluster_intialized = False
                        try:
                            increment += 1
                            if increment == interval or len(row_header) > 200 and increment == interval - 9 and cluster_intialized == False:
                                cluster_intialized = True
                                atypical_cluster = False
                                if ind1[(i + 9)] == 'NA':
                                    atypical_cluster = True
                                    cluster9 = 'cluster-' + string.split(row_header[(new_index + 9)], ':')[0]
                                    if len(row_header) > 200 and str(cluster9) != cluster:
                                        continue
                                elif len(row_header) > 200 and str(ind1[(i + 9)]) != cluster:
                                    continue
                                else:
                                    pvalue, original_term = cluster_elite_terms[new_cluster_id][ci]
                                    term = original_term
                                    if 'GO:' in term:
                                        term = string.split(term, '(')[0]
                                    if ':WP' in term:
                                        term = string.split(term, ':WP')[0]
                                    pvalue = formatpval(str(pvalue))
                                    term += ' p=' + pvalue
                                    if atypical_cluster == False:
                                        term += ' (c' + str(cluster) + ')'
                                    try:
                                        cluster_elite_terms[term] = cluster_elite_terms[(cluster, original_term)]
                                    except Exception:
                                        pass

                                axm.text(label_pos, i - radj, term, horizontalalignment='right', fontsize=ge_fontsize, picker=True, color='blue')
                                increment = 0
                                ci += 1
                        except Exception as e:
                            increment = 0

                    last_cluster = cluster

        def onpick1(event):
            text = event.artist
            print ('onpick1 text:', text.get_text())
            if 'TreeView' in text.get_text():
                try:
                    openTreeView(cdt_file)
                except Exception:
                    print 'Failed to open TreeView'

            else:
                if 'p=' not in text.get_text():
                    webbrowser.open('http://www.genecards.org/cgi-bin/carddisp.pl?gene=' + string.replace(text.get_text(), ' ', ''))
                else:
                    from visualization_scripts import TableViewer
                    header = [
                     'Associated Genes']
                    tuple_list = []
                    for gene in cluster_elite_terms[text.get_text()]:
                        tuple_list.append([gene])

                TableViewer.viewTable(text.get_text(), header, tuple_list)
                cluster_prefix = 'c' + string.split(text.get_text(), '(c')[1][:-1] + '-'
                for geneSet in EliteGeneSets:
                    if geneSet == 'GeneOntology':
                        png_file_dir = elite_dir + '/GO-Elite_results/networks/' + cluster_prefix + 'GO' + '.png'
                    elif geneSet == 'WikiPathways':
                        png_file_dir = elite_dir + '/GO-Elite_results/networks/' + cluster_prefix + 'local' + '.png'
                    elif len(geneSet) > 1:
                        png_file_dir = elite_dir + '/GO-Elite_results/networks/' + cluster_prefix + geneSet + '.png'

                try:
                    alt_png_file_dir = elite_dir + '/GO-Elite_results/networks/' + cluster_prefix + eliteGeneSet + '.png'
                    png_file_dirs = string.split(alt_png_file_dir, 'GO-Elite/')
                    alt_png_file_dir = png_file_dirs[0] + 'GO-Elite/' + png_file_dirs[(-1)]
                except Exception:
                    pass

            if os.name == 'nt':
                try:
                    os.startfile('"' + png_file_dir + '"')
                except Exception:
                    try:
                        os.system('open "' + png_file_dir + '"')
                    except Exception:
                        os.startfile('"' + alt_png_file_dir + '"')

            elif 'darwin' in sys.platform:
                try:
                    os.system('open "' + png_file_dir + '"')
                except Exception:
                    os.system('open "' + alt_png_file_dir + '"')

            elif 'linux' in sys.platform:
                try:
                    os.system('xdg-open "' + png_file_dir + '"')
                except Exception:
                    os.system('xdg-open "' + alt_png_file_dir + '"')

        fig.canvas.mpl_connect('pick_event', onpick1)
        for i in range(x.shape[1]):
            adji = i
            if len(row_header) < 3:
                cadj = len(row_header) * -0.26
            elif len(row_header) < 4:
                cadj = len(row_header) * -0.23
            elif len(row_header) < 6:
                cadj = len(row_header) * -0.18
            elif len(row_header) < 10:
                cadj = len(row_header) * -0.08
            elif len(row_header) < 15:
                cadj = len(row_header) * -0.04
            elif len(row_header) < 20:
                cadj = len(row_header) * -0.05
            elif len(row_header) < 22:
                cadj = len(row_header) * -0.06
            elif len(row_header) < 23:
                cadj = len(row_header) * -0.08
            elif len(row_header) > 200:
                cadj = -2
            else:
                cadj = -0.9
            if len(column_header) > 15:
                adji = i - 0.1
            if len(column_header) > 20:
                adji = i - 0.2
            if len(column_header) > 25:
                adji = i - 0.2
            if len(column_header) > 30:
                adji = i - 0.25
            if len(column_header) > 35:
                adji = i - 0.3
            if len(column_header) > 200:
                column_fontsize = 2
            if column_method != None:
                if len(column_header) < 300:
                    axm.text(adji, cadj, '' + column_header[idx2[i]], rotation=270, verticalalignment='top', fontsize=column_fontsize)
                new_column_header.append(column_header[idx2[i]])
            elif len(column_header) < 300:
                axm.text(adji, cadj, '' + column_header[i], rotation=270, verticalalignment='top', fontsize=column_fontsize)
            else:
                new_column_header.append(column_header[i])

        group_name_list = []
        ind1_clust, ind2_clust = ind1, ind2
        ind1, ind2, group_name_list, cb_status = updateColorBarData(ind1, ind2, new_column_header, new_row_header, row_method)
        if (column_method != None or 'column' in cb_status) and show_color_bars == True:
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])
            cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0', '#000066', '#FFFF00', '#FF1493'])
            if use_default_colors:
                cmap_c = pylab.cm.nipy_spectral
            elif len(unique.unique(ind2)) == 2:
                cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                cmap_c = matplotlib.colors.ListedColormap(['w', 'k'])
            elif len(unique.unique(ind2)) == 3:
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif len(unique.unique(ind2)) == 4:
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif len(unique.unique(ind2)) == 5:
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2)) == 6:
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2)) == 7:
                cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            elif len(unique.unique(ind2)) > 0:
                cmap_c = pylab.cm.nipy_spectral
            dc = numpy.array(ind2, dtype=int)
            dc.shape = (1, len(ind2))
            im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
            axc.set_xticks([])
            if 'hopach' == column_method and len(group_name_list) > 0:
                axc.set_yticklabels(['', 'Groups'], fontsize=10)
            else:
                axc.set_yticks([])
            if len(group_name_list) > 0:
                if 'hopach' == column_method:
                    axcd = fig.add_axes([ax2_x, ax2_y, ax2_w, color_bar_w])
                    cmap_c = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0', '#000066', '#FFFF00', '#FF1493'])
                    if use_default_colors:
                        cmap_c = pylab.cm.nipy_spectral
                    else:
                        if len(unique.unique(ind2_clust)) == 2:
                            cmap_c = matplotlib.colors.ListedColormap(['w', 'k'])
                        elif len(unique.unique(ind2_clust)) == 3:
                            cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                        elif len(unique.unique(ind2_clust)) == 4:
                            cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
                        elif len(unique.unique(ind2_clust)) == 5:
                            cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
                        elif len(unique.unique(ind2_clust)) == 6:
                            cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                        elif len(unique.unique(ind2_clust)) == 7:
                            cmap_c = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                        elif len(unique.unique(ind2_clust)) > 0:
                            cmap_c = pylab.cm.nipy_spectral
                        dc = numpy.array(ind2_clust, dtype=int)
                        dc.shape = (1, len(ind2_clust))
                        im_cd = axcd.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
                        axcd.set_yticklabels(['', 'Clusters'], fontsize=10)
                        axcd.set_xticks([])
                axd = fig.add_axes([axcc_x, axcc_y, axcc_w, axcc_h])
                group_name_list.sort()
                group_colors = map(lambda x: x[0], group_name_list)
                group_names = map(lambda x: x[1], group_name_list)
                cmap_d = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#CCCCE0', '#000066', '#FFFF00', '#FF1493'])
                if len(unique.unique(ind2)) == 2:
                    cmap_d = matplotlib.colors.ListedColormap(['w', 'k'])
                elif len(unique.unique(ind2)) == 3:
                    cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                elif len(unique.unique(ind2)) == 4:
                    cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
                elif len(unique.unique(ind2)) == 5:
                    cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2)) == 6:
                    cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2)) == 7:
                    cmap_d = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                elif len(unique.unique(ind2)) > 0:
                    cmap_d = pylab.cm.nipy_spectral
                dc = numpy.array(group_colors, dtype=int)
                dc.shape = (1, len(group_colors))
                im_c = axd.matshow(dc, aspect='auto', origin='lower', cmap=cmap_d)
                axd.set_yticks([])
                pylab.xticks(range(len(group_names)), group_names, rotation=45, ha='left')
        if show_color_bars == False:
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])
            axc.set_frame_on(False)
        if (row_method != None or 'row' in cb_status) and show_color_bars == True:
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])
            try:
                dr = numpy.array(ind1, dtype=int)
                dr.shape = (len(ind1), 1)
                cmap_r = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF', '#FFFF00', '#FF1493'])
                if len(unique.unique(ind1)) > 4:
                    cmap_r = pylab.cm.nipy_spectral_r
                if len(unique.unique(ind1)) == 2:
                    cmap_r = matplotlib.colors.ListedColormap(['w', 'k'])
                im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
                axr.set_xticks([])
                axr.set_yticks([])
            except Exception:
                row_method = None

        if show_color_bars == False:
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])
            axr.set_frame_on(False)
        axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)
        cb = matplotlib.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
        if 'LineageCorrelations' in dataset_name:
            cb.set_label('Lineage Correlation Z Scores', fontsize=11)
        elif 'Heatmap' in root_dir:
            cb.set_label('GO-Elite Z Scores', fontsize=11)
        else:
            cb.set_label('Differential Expression (log2)', fontsize=10)
        if len(dataset_name) > 30:
            fontsize = 10
        else:
            fontsize = 12.5
        fig.text(0.015, 0.97, dataset_name, fontsize=fontsize)
        pylab.savefig(root_dir + filename, dpi=1000)
        filename = filename[:-3] + 'png'
        pylab.savefig(root_dir + filename, dpi=100)
        includeBackground = False
        try:
            if 'TkAgg' != matplotlib.rcParams['backend']:
                includeBackground = False
        except Exception:
            pass

    if includeBackground:
        fig.text(0.02, 0.07, 'Open heatmap in TreeView (click here)', fontsize=11.5, picker=True, color='red', backgroundcolor='white')
    else:
        fig.text(0.02, 0.07, 'Open heatmap in TreeView (click here)', fontsize=11.5, picker=True, color='red')
    if 'Outlier' in dataset_name and 'Removed' not in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Outlier Genes Genes', root_dir + filename])
    elif 'Relative' in dataset_name:
        graphic_link.append(['Hierarchical Clustering - Significant Genes (Relative comparisons)', root_dir + filename])
    elif 'LineageCorrelations' in filename:
        graphic_link.append(['Hierarchical Clustering - Lineage Correlations', root_dir + filename])
    elif 'MarkerGenes' in filename:
        graphic_link.append(['Hierarchical Clustering - MarkerFinder', root_dir + filename])
    elif 'AltExonConfirmed' in filename:
        graphic_link.append(['Hierarchical Clustering - AltExonConfirmed', root_dir + filename])
    elif 'AltExon' in filename:
        graphic_link.append(['Hierarchical Clustering - AltExon', root_dir + filename])
    elif 'alt_junction' in filename:
        graphic_link.append(['Hierarchical Clustering - Variable Splice-Events', root_dir + filename])
    else:
        graphic_link.append(['Hierarchical Clustering - Significant Genes', root_dir + filename])
    if display:
        proceed = True
        try:
            if 'guide' in justShowTheseIDs:
                proceed = False
        except Exception:
            pass

        if proceed:
            print 'Exporting:', filename
            try:
                pylab.show()
            except Exception:
                None

    fig.clf()
    return


def openTreeView(filename):
    import subprocess
    fn = filepath('AltDatabase/TreeView/TreeView.jar')
    retcode = subprocess.Popen(['java', '-Xmx500m', '-jar', fn, '-r', filename])


def remoteGOElite(elite_dir, SystemCode=None):
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
        input_files = dir_list = unique.read_directory(elite_dir)
        if len(input_files) > 0 and resources_to_analyze != ['']:
            print '\nBeginning to run GO-Elite analysis on all results'
            file_dirs = (elite_dir, None, elite_dir)
            enrichmentAnalysisType = 'ORA'
            variables = (
             species, mod, pathway_permutations, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, returnPathways, file_dirs, enrichmentAnalysisType, root)
            try:
                GO_Elite.remoteAnalysis(variables, 'non-UI Heatmap')
            except Exception:
                print 'GO-Elite failed for:', elite_dir
                print traceback.format_exc()

            if commandLine == False:
                try:
                    UI.openDirectory(elite_dir + '/GO-Elite_results')
                except Exception:
                    None

            cluster_elite_terms, top_genes = importGOEliteResults(elite_dir)
            return (
             cluster_elite_terms, top_genes)
        return ({}, [])
    else:
        return ({}, [])
    return


def importGOEliteResults(elite_dir):
    global eliteGeneSet
    pruned_results = elite_dir + '/GO-Elite_results/CompleteResults/ORA_pruned/pruned-results_z-score_elite.txt'
    if os.path.isfile(pruned_results) == False:
        pruned_results = elite_dir + '/GO-Elite_results/pruned-results_z-score_elite.txt'
    firstLine = True
    cluster_elite_terms = {}
    all_term_length = [0]
    for line in open(pruned_results, 'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data, '\t')
        if firstLine:
            firstLine = False
            try:
                symbol_index = values.index('gene symbols')
            except Exception:
                symbol_index = None

        else:
            try:
                symbol_index = values.index('gene symbols')
            except Exception:
                pass
            else:
                try:
                    eliteGeneSet = string.split(values[0][1:], '-')[1][:-4]
                    try:
                        cluster = str(int(float(string.split(values[0][1:], '-')[0])))
                    except Exception:
                        cluster = string.join(string.split(values[0], '-')[:-1], '-')

                    term = values[2]
                    num_genes_changed = int(values[3])
                    all_term_length.append(len(term))
                    pval = float(values[9])
                    if num_genes_changed > 2:
                        try:
                            cluster_elite_terms[cluster].append([pval, term])
                        except Exception:
                            cluster_elite_terms[cluster] = [[pval, term]]

                        if symbol_index != None:
                            symbols = string.split(values[symbol_index], '|')
                            cluster_elite_terms[(cluster, term)] = symbols
                except Exception as e:
                    pass

    for cluster in cluster_elite_terms:
        cluster_elite_terms[cluster].sort()

    cluster_elite_terms['label-size'] = max(all_term_length)
    top_genes = []
    count = 0
    ranked_genes = elite_dir + '/GO-Elite_results/CompleteResults/ORA_pruned/gene_associations/pruned-gene-ranking.txt'
    for line in open(ranked_genes, 'rU').xreadlines():
        data = line.rstrip()
        values = string.split(data, '\t')
        count += 1
        if len(values) > 2:
            if values[2] != 'Symbol':
                try:
                    top_genes.append((int(values[4]), values[2]))
                except Exception:
                    pass

    top_genes.sort()
    top_genes.reverse()
    top_genes = map(lambda x: x[1], top_genes[:21])
    return (
     cluster_elite_terms, top_genes)


def mergeRotateAroundPointPage(page, page2, rotation, tx, ty):
    from pyPdf import PdfFileWriter, PdfFileReader
    translation = [
     [
      1, 0, 0],
     [
      0, 1, 0],
     [
      -tx, -ty, 1]]
    rotation = math.radians(rotation)
    rotating = [[math.cos(rotation), math.sin(rotation), 0],
     [
      -math.sin(rotation), math.cos(rotation), 0],
     [
      0, 0, 1]]
    rtranslation = [[1, 0, 0],
     [
      0, 1, 0],
     [
      tx, ty, 1]]
    ctm = numpy.dot(translation, rotating)
    ctm = numpy.dot(ctm, rtranslation)
    return page.mergeTransformedPage(page2, [ctm[0][0], ctm[0][1],
     ctm[1][0], ctm[1][1],
     ctm[2][0], ctm[2][1]])


def mergePDFs2(pdf1, pdf2, outPdf):
    from pyPdf import PdfFileWriter, PdfFileReader
    input1 = PdfFileReader(file(pdf1, 'rb'))
    page1 = input1.getPage(0)
    input2 = PdfFileReader(file(pdf2, 'rb'))
    page2 = input2.getPage(0)
    page3 = mergeRotateAroundPointPage(page1, page2, page1.get('/Rotate') or 0, page2.mediaBox.getWidth() / 2, page2.mediaBox.getWidth() / 2)
    output = PdfFileWriter()
    output.addPage(page3)
    outputStream = file(outPdf, 'wb')
    output.write(outputStream)
    outputStream.close()


def mergePDFs(pdf1, pdf2, outPdf):
    from pyPdf import PdfFileWriter, PdfFileReader
    input1 = PdfFileReader(file(pdf1, 'rb'))
    page1 = input1.getPage(0)
    page1.mediaBox.upperRight = (page1.mediaBox.getUpperRight_x(), page1.mediaBox.getUpperRight_y())
    input2 = PdfFileReader(file(pdf2, 'rb'))
    page2 = input2.getPage(0)
    page2.mediaBox.getLowerLeft_x = (page2.mediaBox.getLowerLeft_x(), page2.mediaBox.getLowerLeft_y())
    page2.mergePage(page1)
    output = PdfFileWriter()
    output.addPage(page1)
    outputStream = file(outPdf, 'wb')
    output.write(outputStream)
    outputStream.close()


def inverseDist(value):
    if value == 0:
        value = 1
    return math.log(value, 2)


def getGOEliteExportDir(root_dir, dataset_name):
    if 'AltResults' in root_dir:
        root_dir = string.split(root_dir, 'AltResults')[0]
    if 'ExpressionInput' in root_dir:
        root_dir = string.split(root_dir, 'ExpressionInput')[0]
    if 'ExpressionOutput' in root_dir:
        root_dir = string.split(root_dir, 'ExpressionOutput')[0]
    if 'DataPlots' in root_dir:
        root_dir = string.replace(root_dir, 'DataPlots', 'GO-Elite')
        elite_dir = root_dir
    else:
        elite_dir = root_dir + '/GO-Elite'
    try:
        os.mkdir(elite_dir)
    except Exception:
        pass

    return elite_dir + '/clustering/' + dataset_name


def systemCodeCheck(IDs):
    import gene_associations
    id_type_db = {}
    for id in IDs:
        id_type = gene_associations.predictIDSourceSimple(id)
        try:
            id_type_db[id_type] += 1
        except Exception:
            id_type_db[id_type] = 1

    id_type_count = []
    for i in id_type_db:
        id_type_count.append((id_type_db[i], i))

    id_type_count.sort()
    id_type = id_type_count[(-1)][(-1)]
    return id_type


def exportFlatClusterData(filename, root_dir, dataset_name, new_row_header, new_column_header, xt, ind1, ind2, vmax, display):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    filename = string.replace(filename, '.pdf', '.txt')
    export_text = export.ExportFile(filename)
    column_header = string.join(['UID', 'row_clusters-flat'] + new_column_header, '\t') + '\n'
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat', ''] + map(str, ind2), '\t') + '\n'
    export_text.write(column_clusters)
    try:
        elite_dir = getGOEliteExportDir(root_dir, dataset_name)
    except Exception:
        elite_dir = None

    elite_columns = string.join(['InputID', 'SystemCode'])
    try:
        sy = systemCodeCheck(new_row_header)
    except Exception:
        sy = None
    else:
        i = 0
        cluster_db = {}
        export_lines = []
        for row in xt:
            try:
                id = new_row_header[i]
                original_id = str(id)
                if sy == 'Ae' and '--' in id:
                    cluster = 'cluster-' + string.split(id, ':')[0]
                elif sy == '$En:Sy':
                    cluster = 'cluster-' + string.split(id, ':')[0]
                elif sy == 'S' and ':' in id:
                    cluster = 'cluster-' + string.split(id, ':')[0]
                elif sy == 'Sy' and ':' in id:
                    cluster = 'cluster-' + string.split(id, ':')[0]
                else:
                    cluster = 'c' + str(ind1[i])
            except Exception:
                pass
            else:
                try:
                    if 'MarkerGenes' in originalFilename:
                        cluster = 'cluster-' + string.split(id, ':')[0]
                        id = string.split(id, ':')[1]
                        if ' ' in id:
                            id = string.split(id, ' ')[0]
                        if 'G000' in id:
                            sy = 'En'
                        else:
                            sy = 'Sy'
                except Exception:
                    pass
                else:
                    try:
                        cluster_db[cluster].append(id)
                    except Exception:
                        cluster_db[cluster] = [id]
                    else:
                        try:
                            export_lines.append(string.join([original_id, str(ind1[i])] + map(str, row), '\t') + '\n')
                        except Exception:
                            export_lines.append(string.join([original_id, 'NA'] + map(str, row), '\t') + '\n')
                        else:
                            i += 1

        export_lines.reverse()
        for line in export_lines:
            export_text.write(line)

        export_text.close()
        allGenes = {}
        for cluster in cluster_db:
            export_elite = export.ExportFile(elite_dir + '/' + cluster + '.txt')
            if sy == None:
                export_elite.write('ID\n')
            else:
                export_elite.write('ID\tSystemCode\n')
            for id in cluster_db[cluster]:
                try:
                    i1, i2 = string.split(id, ' ')
                    if i1 == i2:
                        id = i1
                except Exception:
                    pass
                else:
                    if sy == '$En:Sy':
                        id = string.split(id, ':')[1]
                        ids = string.split(id, ' ')
                        if 'ENS' in ids[0] or 'G0000' in ids[0]:
                            id = ids[0]
                        else:
                            id = ids[(-1)]
                        sc = 'En'
                    elif sy == 'Sy' and ':' in id:
                        id = string.split(id, ':')[1]
                        ids = string.split(id, ' ')
                        sc = 'Sy'
                    elif sy == 'En:Sy':
                        id = string.split(id, ' ')[0]
                        sc = 'En'
                    elif sy == 'Ae':
                        if '--' in id:
                            sc = 'En'
                            id = string.split(id, ':')[(-1)]
                        else:
                            l = string.split(id, ':')
                            if len(l) == 2:
                                id = string.split(id, ':')[0]
                            if len(l) == 3:
                                id = string.split(id, ':')[1]
                            sc = 'En'
                            if ' ' in id:
                                ids = string.split(id, ' ')
                                if 'ENS' in ids[(-1)] or 'G0000' in ids[(-1)]:
                                    id = ids[(-1)]
                                else:
                                    id = ids[0]
                    elif sy == 'En' and '&' in id:
                        for i in string.split(id, '&'):
                            if 'G0000' in i:
                                id = i
                                sc = 'En'
                                break

                    elif sy == 'Sy' and 'EFN' in id:
                        sc = 'En'
                    else:
                        sc = sy
                    if sy == 'S':
                        if ':' in id:
                            id = string.split(id, ':')[(-1)]
                            sc = 'Ae'
                    if '&' in id:
                        sc = 'Ae'
                    if len(id) == 9 and 'SRS' in id or len(id) == 15 and 'TCGA-' in id:
                        sc = 'En'
                    try:
                        export_elite.write(id + '\t' + sc + '\n')
                    except Exception:
                        export_elite.write(id + '\n')
                    else:
                        allGenes[id] = []

            export_elite.close()

        try:
            if storeGeneSetName != None:
                if len(storeGeneSetName) > 0 and 'guide' not in justShowTheseIDs:
                    exportCustomGeneSet(storeGeneSetName, species, allGenes)
                    print 'Exported geneset to "StoredGeneSets"'
        except Exception:
            pass

    filename = string.replace(filename, '.txt', '.cdt')
    if display:
        try:
            exportJTV(filename, new_column_header, new_row_header, vmax=vmax)
        except Exception:
            pass

    export_cdt = export.ExportFile(filename)
    column_header = string.join(['UNIQID', 'NAME', 'GWEIGHT'] + new_column_header, '\t') + '\n'
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT', '', ''] + ['1'] * len(new_column_header), '\t') + '\n'
    export_cdt.write(eweight)
    i = 0
    cdt_lines = []
    for row in xt:
        cdt_lines.append(string.join([new_row_header[i]] * 2 + ['1'] + map(str, row), '\t') + '\n')
        i += 1

    cdt_lines.reverse()
    for line in cdt_lines:
        export_cdt.write(line)

    export_cdt.close()
    return (
     elite_dir, filename, sc)


def exportJTV(cdt_dir, column_header, row_header, vmax=None):
    filename = string.replace(cdt_dir, '.cdt', '.jtv')
    export_jtv = export.ExportFile(filename)
    cscale = '3'
    if len(column_header) > 100:
        cscale = '1.5'
    if len(column_header) > 200:
        cscale = '1.1'
    if len(column_header) > 300:
        cscale = '0.6'
    if len(column_header) > 400:
        cscale = '0.3'
    hscale = '5'
    if len(row_header) < 50:
        hscale = '10'
    if len(row_header) > 100:
        hscale = '3'
    if len(row_header) > 500:
        hscale = '1'
    if len(row_header) > 1000:
        hscale = '0.5'
    contrast = str(float(vmax) / 4)[:4]
    config = '<DocumentConfig><UrlExtractor/><ArrayUrlExtractor/><MainView><ColorExtractor>'
    config += '<ColorSet down="#00FFFF"/></ColorExtractor><ArrayDrawer/><GlobalXMap>'
    config += '<FixedMap type="Fixed" scale="' + cscale + '"/><FillMap type="Fill"/><NullMap type="Null"/>'
    config += '</GlobalXMap><GlobalYMap><FixedMap type="Fixed" scale="' + hscale + '"/><FillMap type="Fill"/>'
    config += '<NullMap type="Null"/></GlobalYMap><ZoomXMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config += '<NullMap type="Null"/></ZoomXMap><ZoomYMap><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config += '<NullMap type="Null"/></ZoomYMap><TextView><TextView><GeneSummary/></TextView><TextView>'
    config += '<GeneSummary/></TextView><TextView><GeneSummary/></TextView></TextView><ArrayNameView>'
    config += '<ArraySummary included="0"/></ArrayNameView><AtrSummary/><GtrSummary/></MainView><Views>'
    config += '<View type="Dendrogram" dock="1"><ColorExtractor contrast="' + contrast + '"><ColorSet up="#FFFF00" down="#00CCFF"/>'
    config += '</ColorExtractor><ArrayDrawer/><GlobalXMap current="Fill"><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config += '<NullMap type="Null"/></GlobalXMap><GlobalYMap current="Fill"><FixedMap type="Fixed"/><FillMap type="Fill"/>'
    config += '<NullMap type="Null"/></GlobalYMap><ZoomXMap><FixedMap type="Fixed"/><FillMap type="Fill"/><NullMap type="Null"/>'
    config += '</ZoomXMap><ZoomYMap current="Fixed"><FixedMap type="Fixed"/><FillMap type="Fill"/><NullMap type="Null"/></ZoomYMap>'
    config += '<TextView><TextView><GeneSummary/></TextView><TextView><GeneSummary/></TextView><TextView><GeneSummary/></TextView>'
    config += '</TextView><ArrayNameView><ArraySummary included="0"/></ArrayNameView><AtrSummary/><GtrSummary/></View></Views></DocumentConfig>'
    export_jtv.write(config)


def updateColorBarData(ind1, ind2, column_header, row_header, row_method):
    """ Replace the top-level cluster information with group assignments for color bar coloring (if group data present)"""
    cb_status = 'original'
    group_number_list = []
    group_name_list = []
    try:
        if column_header[0] in GroupDB:
            cb_status = 'column'
            for header in column_header:
                group, color, color_num = GroupDB[header]
                group_number_list.append(color_num)
                if (color_num, group) not in group_name_list:
                    group_name_list.append((color_num, group))

            ind2 = group_number_list
        if row_header[0] in GroupDB and row_method == None:
            group_number_list = []
            if cb_status == 'column':
                cb_status = 'column-row'
            else:
                cb_status = 'row'
            for header in row_header:
                group, color, color_num = GroupDB[header]
                group_number_list.append(color_num)

            ind1 = group_number_list
    except Exception:
        None

    return (
     ind1, ind2, group_name_list, cb_status)


def ConvertFromHex(color1, color2, color3):
    c1tuple = tuple(ord(c) for c in color1.lsstrip('0x').decode('hex'))
    c2tuple = tuple(ord(c) for c in color2.lsstrip('0x').decode('hex'))
    c3tuple = tuple(ord(c) for c in color3.lsstrip('0x').decode('hex'))


def RedBlackSkyBlue():
    cdict = {}
    my_cmap = mc.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


def RedBlackBlue():
    cdict = {}
    my_cmap = mc.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


def RedBlackGreen():
    cdict = {}
    my_cmap = mc.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


def YellowBlackBlue():
    cdict = {}
    my_cmap = mc.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


def BlackYellowBlue():
    cdict = {}
    my_cmap = mc.LinearSegmentedColormap('my_colormap', cdict, 256)
    return my_cmap


def cleanUpLine(line):
    line = string.replace(line, '\n', '')
    line = string.replace(line, '\\c', '')
    data = string.replace(line, '\r', '')
    data = string.replace(data, '"', '')
    return data


def filepath(filename):
    fn = unique.filepath(filename)
    return fn


def remoteImportData(filename, geneFilter=None, reverseOrder=True):
    matrix, column_header, row_header, dataset_name, group_db = importData(filename, geneFilter=geneFilter, reverseOrder=reverseOrder)
    try:
        return (matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters)
    except:
        return (
         matrix, column_header, row_header, dataset_name, group_db, [], [])


def importData(filename, Normalize=False, reverseOrder=True, geneFilter=None, zscore=False, forceClusters=False):
    global priorColumnClusters
    global priorRowClusters
    try:
        if len(priorColumnClusters) > 0:
            priorColumnClusters = None
            priorRowClusters = None
    except Exception:
        pass
    else:
        getRowClusters = False
        start_time = time.time()
        fn = filepath(filename)
        matrix = []
        original_matrix = []
        row_header = []
        overwriteGroupNotations = True
        x = 0
        inputMax = 0
        inputMin = 100
        filename = string.replace(filename, '\\', '/')
        dataset_name = string.split(filename, '/')[(-1)][:-4]
        if '.cdt' in filename:
            start = 3
        else:
            start = 1
        for line in open(fn, 'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data, '\t')
            if x == 0:
                if '.cdt' in filename:
                    t = [t[0]] + t[3:]
                if t[1] == 'row_clusters-flat':
                    t = t = [
                     t[0]] + t[2:]
                new_headers = []
                temp_groups = {}
                original_headers = t[1:]
                if ('exp.' in filename or 'filteredExp.' in filename or 'MarkerGene' in filename) and forceClusters == False:
                    if overwriteGroupNotations:
                        for i in t:
                            if ':' in i:
                                group, i = string.split(i, ':')
                                new_headers.append(i)
                                temp_groups[i] = group
                            else:
                                new_headers.append(i)

                    filename = string.replace(filename, '-steady-state.txt', '.txt')
                    try:
                        import ExpressionBuilder
                        try:
                            sample_group_db = ExpressionBuilder.simplerGroupImport(filename)
                        except Exception:
                            sample_group_db = {}

                        if len(temp_groups) > 0 and len(sample_group_db) == 0:
                            sample_group_db = temp_groups
                        if len(new_headers) > 0:
                            t = new_headers
                        new_headers = []
                        for v in t:
                            if v in sample_group_db:
                                v = sample_group_db[v] + ':' + v
                            new_headers.append(v)

                        t = new_headers
                    except Exception:
                        pass

                group_db, column_header = assignGroupColors(t[1:])
                x = 1
            elif 'column_clusters-flat' in t:
                try:
                    if 'NA' in t:
                        kill
                    try:
                        if forceClusters == False:
                            prior = map(lambda x: int(float(x)), t[2:])
                        else:
                            prior = map(lambda x: x, t[2:])
                    except Exception:
                        index = 0
                        c = 1
                        prior = []
                        clusters = {}
                        for i in t[2:]:
                            original_headers[index] = i + ':' + original_headers[index]
                            if i in clusters:
                                c1 = clusters[i]
                            else:
                                c1 = c
                                clusters[i] = c1
                                c += 1
                            prior.append(c1)
                            index += 1

                        if len(temp_groups) == 0:
                            if '-ReOrdered.txt' not in filename:
                                group_db, column_header = assignGroupColors(original_headers)

                    priorColumnClusters = prior
                except Exception:
                    pass

                start = 2
                getRowClusters = True
                priorRowClusters = []
            elif 'EWEIGHT' in t:
                pass
            else:
                gene = t[0]
                if geneFilter == None:
                    proceed = True
                elif gene in geneFilter:
                    proceed = True
                else:
                    proceed = False
                if proceed:
                    nullsPresent = False
                    try:
                        s = map(float, t[start:])
                    except Exception:
                        nullsPresent = True
                        s = []
                        for value in t[start:]:
                            try:
                                s.append(float(value))
                            except Exception:
                                s.append(0.000101)

                    else:
                        original_matrix.append(s)
                        try:
                            if max(s) > inputMax:
                                inputMax = max(s)
                        except:
                            continue

                    if min(s) < inputMin:
                        inputMin = min(s)
                    if Normalize != False:
                        with warnings.catch_warnings():
                            warnings.filterwarnings('ignore', category=UserWarning)
                            if Normalize == 'row mean':
                                avg = numpy.mean(s)
                            else:
                                avg = numpy.median(s)
                        if nullsPresent:
                            s = []
                            for value in t[start:]:
                                try:
                                    s.append(float(value) - avg)
                                except Exception:
                                    s.append(0.000101)

                        else:
                            s = map(lambda x: x - avg, s)
                    if ' ' in gene:
                        try:
                            g1, g2 = string.split(gene, ' ')
                            if g1 == g2:
                                gene = g1
                        except Exception:
                            pass

                    if getRowClusters:
                        try:
                            priorRowClusters.append(int(float(t[1])))
                        except Exception:
                            pass

                    if zscore:
                        avg = numpy.mean(s)
                        std = numpy.std(s)
                        if std == 0:
                            std = 0.1
                        try:
                            s = map(lambda x: (x - avg) / std, s)
                        except Exception:
                            pass

                    if geneFilter == None:
                        matrix.append(s)
                        row_header.append(gene)
                    elif gene in geneFilter:
                        matrix.append(s)
                        row_header.append(gene)
                    x += 1

        if inputMax > 100:
            print 'Converting values to log2...'
            matrix = []
            k = 0
            if inputMin == 0:
                increment = 1
            else:
                increment = 1
            for s in original_matrix:
                if 'counts.' in filename:
                    s = map(lambda x: math.log(x + 1, 2), s)
                else:
                    try:
                        s = map(lambda x: math.log(x + increment, 2), s)
                    except Exception:
                        print filename
                        print Normalize
                        print row_header[k], min(s), max(s)
                        kill

                if Normalize != False:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=UserWarning)
                        if Normalize == 'row mean':
                            avg = numpy.average(s)
                        else:
                            avg = avg = numpy.median(s)
                    s = map(lambda x: x - avg, s)
                if zscore:
                    avg = numpy.mean(s)
                    std = numpy.std(s)
                    if std == 0:
                        std = 0.1
                    try:
                        s = map(lambda x: (x - avg) / std, s)
                    except Exception:
                        pass

                matrix.append(s)
                k += 1

            del original_matrix
        if zscore:
            print 'Converting values to normalized z-scores...'
        if reverseOrder == True:
            matrix.reverse()
            row_header.reverse()
        time_diff = str(round(time.time() - start_time, 1))
        try:
            print '%d rows and %d columns imported for %s in %s seconds...' % (len(matrix), len(column_header), dataset_name, time_diff)
        except Exception:
            print 'No data in input file.'
            force_error

        group_db2, row_header2 = assignGroupColors(list(row_header))
        for i in group_db2:
            if i not in group_db:
                group_db[i] = group_db2[i]

    return (
     matrix, column_header, row_header, dataset_name, group_db)


def importSIF(filename):
    fn = filepath(filename)
    edges = []
    x = 0
    if '/' in filename:
        dataset_name = string.split(filename, '/')[(-1)][:-4]
    else:
        dataset_name = string.split(filename, '\\')[(-1)][:-4]
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        parent, type, child = string.split(data, '\t')
        if 'AltAnalyze' in dataset_name:
            edges.append([parent, child, type])
        else:
            if '(' in parent:
                parent = string.split(parent, '(')[0]
            if ':' in child:
                child = string.split(child, ':')[1]
            if 'TF' in dataset_name or 'UserSuppliedAssociations' in dataset_name or 'WGRV' in dataset_name:
                edges.append([parent, child, type])
            else:
                edges.append([child, parent, type])

    edges = unique.unique(edges)
    return edges


def assignGroupColors(t):
    """ Assign a unique color to each group. Optionally used for cluster display. """
    column_header = []
    group_number_db = {}
    groupNamesPresent = False
    for i in t:
        if ':' in i:
            groupNamesPresent = True

    for i in t:
        repls = {}
        i = reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        if ':' in i:
            group, j = string.split(i, ':')[:2]
            group_number_db[group] = []
        elif groupNamesPresent:
            group_number_db['UNK'] = []
            i = 'UNK:' + i
        column_header.append(i)

    k = 0
    group_db = {}
    color_db = {}
    color_list = [
     'r', 'b', 'y', 'g', 'w', 'k', 'm']
    if len(group_number_db) > 3:
        color_list = []
        cm = pylab.cm.get_cmap('nipy_spectral')
        for i in range(len(group_number_db)):
            color_list.append(cm(1.0 * i / len(group_number_db)))

    t.sort()
    for i in column_header:
        repls = {}
        i = reduce(lambda a, kv: a.replace(*kv), repls.iteritems(), i)
        if ':' in i:
            group, j = string.split(i, ':')[:2]
            try:
                color, ko = color_db[group]
            except Exception:
                try:
                    color_db[group] = (color_list[k], k)
                except Exception:
                    rgb = tuple(scipy.rand(3))
                    color_list.append(rgb)
                    color_db[group] = (color_list[k], k)

                color, ko = color_db[group]
                k += 1

            group_db[i] = (
             group, color, ko)

    return (
     group_db, column_header)


def verifyFile(filename):
    status = 'not found'
    try:
        fn = filepath(filename)
        for line in open(fn, 'rU').xreadlines():
            status = 'found'
            break

    except Exception:
        status = 'not found'

    return status


def AppendOrWrite(export_path):
    export_path = filepath(export_path)
    status = verifyFile(export_path)
    if status == 'not found':
        export_data = export.ExportFile(export_path)
    else:
        export_data = open(export_path, 'a')
    return (export_path, export_data, status)


def exportCustomGeneSet(geneSetName, species, allGenes):
    for gene in allGenes:
        break

    if 'ENS' not in gene:
        try:
            import gene_associations
            from import_scripts import OBO_import
            gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', 'Ensembl-Symbol'))
            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        except Exception:
            symbol_to_gene = {}

    if species != None:
        export_path, export_data, status = AppendOrWrite('AltDatabase/goelite/' + species + '/gene-mapp/Ensembl-StoredGeneSets.txt')
        stored_lines = []
        for line in open(export_path, 'rU').xreadlines():
            stored_lines.append(line)

        if status == 'not found':
            export_data.write('GeneID\tEmpty\tGeneSetName\n')
        for gene in allGenes:
            if ' ' in gene:
                a, b = string.split(gene, ' ')
                if 'ENS' in a:
                    gene = a
                else:
                    gene = b
            if 'ENS' not in gene and gene in symbol_to_gene:
                gene = symbol_to_gene[gene][0]
            line = gene + '\t\t' + geneSetName + '\n'
            if line not in stored_lines:
                export_data.write(line)

        export_data.close()
    else:
        print 'Could not store since no species name provided.'
    return


def writetSNEScores(scores, outputdir):
    export_obj = export.ExportFile(outputdir)
    for matrix_row in scores:
        matrix_row = map(str, matrix_row)
        export_obj.write(string.join(matrix_row, '\t') + '\n')

    export_obj.close()


def importtSNEScores(inputdir):
    scores = []
    for line in open(inputdir, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        t = map(float, t)
        scores.append(t)

    return scores


def runUMAP(matrix, column_header, dataset_name, group_db, display=False, showLabels=False, row_header=None, colorByGene=None, species=None, reimportModelScores=True, method='UMAP', rootDir='', finalOutputDir=''):
    global graphic_link
    global root_dir
    graphic_link = []
    root_dir = rootDir
    tSNE(matrix, column_header, dataset_name, group_db, display=False, showLabels=False, row_header=None, colorByGene=None, species=None, reimportModelScores=True, method='UMAP')
    import shutil
    filename = 'Clustering-' + dataset_name + '-' + method + '.pdf'
    filename = string.replace(filename, 'Clustering-Clustering', 'Clustering')
    new_file = finalOutputDir + filename
    new_file = string.replace(new_file, 'Clustering-', '')
    new_file = string.replace(new_file, 'exp.', '')
    old_file = root_dir + filename
    shutil.move(old_file, new_file)
    filename = filename[:-3] + 'png'
    new_file = finalOutputDir + filename
    new_file = string.replace(new_file, 'Clustering-', '')
    new_file = string.replace(new_file, 'exp.', '')
    old_file = root_dir + filename
    shutil.move(old_file, new_file)
    old_file = root_dir + dataset_name + '-' + method + '_scores.txt'
    new_file = finalOutputDir + dataset_name + '-' + method + '_coordinates.txt'
    new_file = string.replace(new_file, 'exp.', '')
    shutil.move(old_file, new_file)
    return


def tSNE(matrix, column_header, dataset_name, group_db, display=True, showLabels=False, row_header=None, colorByGene=None, species=None, reimportModelScores=True, method='tSNE', maskGroups=None):
    try:
        prior_clusters = priorColumnClusters
    except Exception:
        prior_clusters = []
    else:
        try:
            if priorColumnClusters == None:
                prior_clusters = []
        except:
            pass
        else:
            try:
                if len(prior_clusters) > 0 and len(group_db) == 0:
                    newColumnHeader = []
                    i = 0
                    for sample_name in column_header:
                        newColumnHeader.append(str(prior_clusters[i]) + ':' + sample_name)
                        i += 1

                    group_db, column_header = assignGroupColors(newColumnHeader)
            except Exception as e:
                print traceback.format_exc()
                group_db = {}
            else:
                if reimportModelScores:
                    print 'Re-importing', method, 'model scores rather than calculating from scratch',
                    print root_dir + dataset_name + '-' + method + '_scores.txt'
                    try:
                        scores = importtSNEScores(root_dir + dataset_name + '-' + method + '_scores.txt')
                        print '...import finished'
                    except Exception:
                        reimportModelScores = False
                        print '...no existing score file found'

                if reimportModelScores == False:
                    X = matrix.T
                    print 'Performing', method
                    if method == 'tSNE' or method == 't-SNE':
                        from sklearn.manifold import TSNE
                        model = TSNE(n_components=2)
                    if method == 'UMAP':
                        try:
                            import umap
                            model = umap.UMAP(n_neighbors=50, min_dist=0.75, metric='correlation')
                        except:
                            try:
                                from visualization_scripts.umap_learn import umap
                                model = umap.UMAP(n_neighbors=50, min_dist=0.75, metric='correlation')
                            except:
                                from visualization_scripts.umap_learn_single import umap
                                model = umap.UMAP(n_neighbors=50, min_dist=0.75, metric='correlation')

                        print 'UMAP run'
                    scores = model.fit_transform(X)
                    writetSNEScores(scores, root_dir + dataset_name + '-' + method + '_scores.txt')
                if maskGroups != None:
                    group_name, restricted_samples = maskGroups
                    dataset_name += '-' + group_name
                scoresT = zip(*scores)
                exclude = {}
                try:
                    for vector in scoresT:
                        lower1th, median_val, upper99th, int_qrt_range = statistics.iqr(list(vector), k1=99.9, k2=0.1)
                        index = 0
                        for i in vector:
                            if i > upper99th + 1 or i < lower1th - 1:
                                exclude[index] = None
                            index += 1

                except Exception:
                    pass

                print 'Not showing', len(exclude), 'outlier samples.'
                fig = pylab.figure()
                ax = fig.add_subplot(111)
                pylab.xlabel(method.upper() + '-X')
                pylab.ylabel(method.upper() + '-Y')
                axes = getAxesTransposed(scores, exclude=exclude)
                pylab.axis(axes)
                marker_size = 15
                if len(column_header) > 20:
                    marker_size = 12
                if len(column_header) > 40:
                    marker_size = 10
                if len(column_header) > 150:
                    marker_size = 7
                if len(column_header) > 500:
                    marker_size = 5
                if len(column_header) > 1000:
                    marker_size = 4
                if len(column_header) > 2000:
                    marker_size = 3
                if len(column_header) > 4000:
                    marker_size = 2
                if len(column_header) > 6000:
                    marker_size = 1
                if colorByGene != None and len(matrix) == 0:
                    print 'Gene %s not found in the imported dataset... Coloring by groups.' % colorByGene
                if colorByGene != None and len(matrix) > 0:
                    gene_translation_db = {}
                    matrix = numpy.array(matrix)
                    min_val = matrix.min()
                    if ' ' in colorByGene:
                        genes = string.split(colorByGene, ' ')
                    else:
                        genes = [
                         colorByGene]
                    genePresent = False
                    numberGenesPresent = []
                    for gene in genes:
                        if gene in row_header:
                            numberGenesPresent.append(gene)
                            genePresent = True

                    if len(numberGenesPresent) == 0:
                        try:
                            import gene_associations
                            from import_scripts import OBO_import
                            gene_to_symbol = gene_associations.getGeneToUid(species, ('hide',
                                                                                      'Ensembl-Symbol'))
                            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
                            for symbol in genes:
                                if symbol in symbol_to_gene:
                                    gene = symbol_to_gene[symbol][0]
                                    if gene in row_header:
                                        numberGenesPresent.append(gene)
                                        genePresent = True
                                        gene_translation_db[symbol] = gene

                        except Exception:
                            pass

                    numberGenesPresent = len(numberGenesPresent)
                    if numberGenesPresent == 1:
                        cm = pylab.cm.get_cmap('Reds')
                    elif numberGenesPresent == 2:
                        cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                    elif numberGenesPresent == 3:
                        cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
                    elif numberGenesPresent == 4:
                        cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
                    elif numberGenesPresent == 5:
                        cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
                    elif numberGenesPresent == 6:
                        cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                    elif numberGenesPresent == 7:
                        cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
                    else:
                        cm = pylab.cm.get_cmap('gist_rainbow')

                    def get_cmap(n, name='hsv'):
                        """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
            RGB color; the keyword argument name must be a standard mpl colormap name."""
                        return pylab.cm.get_cmap(name, n)

                    if genePresent:
                        dataset_name += '-' + colorByGene
                        group_db = {}
                        bestGeneAssociated = {}
                        k = 0
                        for gene in genes:
                            try:
                                try:
                                    i = row_header.index(gene)
                                except Exception:
                                    i = row_header.index(gene_translation_db[gene])
                                else:
                                    values = map(float, matrix[i])
                                    min_val = min(values)
                                    bin_size = (max(values) - min_val) / 8
                                    max_val = max(values)
                                    ranges = []
                                    iz = min_val
                                    while iz < max(values) - bin_size / 100:
                                        r = (
                                         iz, iz + bin_size)
                                        if len(ranges) == 7:
                                            r = (
                                             iz, max_val)
                                        ranges.append(r)
                                        iz += bin_size

                                color_db = {}
                                colors = get_cmap(len(genes))
                                for i in range(len(ranges)):
                                    if i == 0:
                                        color = '#C0C0C0'
                                    elif numberGenesPresent == 1:
                                        color = cm(1.0 * i / len(ranges))
                                    elif i > 2:
                                        if len(genes) < 8:
                                            color = cm(k)
                                        else:
                                            color = colors(k)
                                    else:
                                        color = '#C0C0C0'
                                    color_db[ranges[i]] = color

                                i = 0
                                for val in values:
                                    sample = column_header[i]
                                    for l, u in color_db:
                                        range_index = ranges.index((l, u))
                                        if val >= l and val <= u:
                                            color = color_db[(l, u)]
                                            color_label = [gene + '-range: ' + str(l)[:4] + '-' + str(u)[:4], color, '']
                                            group_db[sample] = color_label
                                            try:
                                                bestGeneAssociated[sample].append([range_index, val, color_label])
                                            except Exception:
                                                bestGeneAssociated[sample] = [[range_index, val, color_label]]

                                    i += 1

                                if len(genes) > 1:
                                    for sample in bestGeneAssociated:
                                        bestGeneAssociated[sample].sort()
                                        color_label = bestGeneAssociated[sample][(-1)][(-1)]
                                        if numberGenesPresent > 1:
                                            index = bestGeneAssociated[sample][(-1)][0]
                                            if index > 2:
                                                gene = string.split(color_label[0], '-')[0]
                                            else:
                                                gene = 'Null'
                                            color_label[0] = gene
                                        group_db[sample] = color_label

                            except Exception:
                                print [
                                 gene], 'not found in rows...'
                            else:
                                k += 1

                    else:
                        print [
                         colorByGene], 'not found in rows...'
                pylab.title(method + ' - ' + dataset_name)
                group_names = {}
                i = 0
                for sample_name in column_header:
                    if maskGroups != None:
                        base_name = sample_name
                        if ':' in sample_name:
                            base_name = string.split(base_name, ':')[1]
                        if base_name not in restricted_samples:
                            exclude[i] = None
                    if i not in exclude:
                        try:
                            group_name, color, k = group_db[sample_name]
                            if group_name not in group_names:
                                label = group_name
                            else:
                                label = None
                            group_names[group_name] = color
                        except Exception:
                            color = 'r'
                            label = None

                        ax.plot(scores[i][0], scores[i][1], color=color, marker='o', markersize=marker_size, label=label, markeredgewidth=0, picker=True)
                        if showLabels:
                            try:
                                sample_name = '   ' + string.split(sample_name, ':')[1]
                            except Exception:
                                pass

                            ax.text(scores[i][0], scores[i][1], sample_name, fontsize=11)
                    i += 1

            group_count = []
            for i in group_db:
                if group_db[i][0] not in group_count:
                    group_count.append(group_db[i][0])

            Lfontsize = 8
            if len(group_count) > 20:
                Lfontsize = 10
            if len(group_count) > 30:
                Lfontsize = 8
            if len(group_count) > 40:
                Lfontsize = 6
            if len(group_count) > 50:
                Lfontsize = 5
            i = 0
            box = ax.get_position()
            if len(group_count) > 0:
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                try:
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=Lfontsize)
                except Exception:
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            else:
                ax.set_position([box.x0, box.y0, box.width, box.height])
                pylab.legend(loc='upper left', prop={})
            filename = 'Clustering-' + dataset_name + '-' + method + '.pdf'
            filename = string.replace(filename, 'Clustering-Clustering', 'Clustering')
            try:
                pylab.savefig(root_dir + filename)
            except Exception:
                None

        filename = filename[:-3] + 'png'
        try:
            pylab.savefig(root_dir + filename)
        except Exception:
            None

    graphic_link.append(['Principal Component Analysis', root_dir + filename])
    if display:
        print 'Exporting:', filename
        try:
            pylab.show()
        except Exception:
            pass

    return


def excludeHighlyCorrelatedHits(x, row_header):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        D1 = numpy.corrcoef(x)
    i = 0
    exclude = {}
    gene_correlations = {}
    include = []
    for score_ls in D1:
        k = 0
        for v in score_ls:
            if str(v) != 'nan':
                if v > 1.0 and k != i:
                    if row_header[i] not in exclude:
                        exclude[row_header[k]] = []
            k += 1

        i += 1

    return exclude


def PrincipalComponentAnalysis(matrix, column_header, row_header, dataset_name, group_db, display=False, showLabels=True, algorithm='SVD', geneSetName=None, species=None, pcA=1, pcB=2, colorByGene=None, reimportModelScores=True):
    print 'Performing Principal Component Analysis...'
    from numpy import mean, cov, double, cumsum, dot, linalg, array, rank
    try:
        prior_clusters = priorColumnClusters
    except Exception:
        prior_clusters = []

    if prior_clusters == None:
        prior_clusters = []
    try:
        if len(prior_clusters) > 0 and len(group_db) == 0:
            newColumnHeader = []
            i = 0
            for sample_name in column_header:
                newColumnHeader.append(str(prior_clusters[i]) + ':' + sample_name)
                i += 1

            group_db, column_header = assignGroupColors(newColumnHeader)
    except Exception as e:
        print traceback.format_exc()
        group_db = {}
    else:
        pcA -= 1
        pcB -= 1
        label1 = ''
        label2 = ''
        if algorithm == 'SVD':
            use_svd = True
        else:
            use_svd = False
        if reimportModelScores:
            print 'Re-importing PCA model scores rather than calculating from scratch',
            print root_dir + dataset_name + '-PCA_scores.txt'
            try:
                scores = importtSNEScores(root_dir + dataset_name + '-PCA_scores.txt')
                print '...import finished'
                matrix = zip(*matrix)
            except Exception:
                reimportModelScores = False
                print '...no existing score file found'

        if reimportModelScores == False:
            Mdif = matrix / matrix.std()
            Mdif = Mdif.T
            u, s, vt = svd(Mdif, 0)
            fracs = s ** 2 / np.sum(s ** 2)
            entropy = -sum(fracs * np.log(fracs)) / np.log(np.min(vt.shape))
            label1 = 'PC%i (%2.1f%%)' % (pcA + 1, fracs[0] * 100)
            label2 = 'PC%i (%2.1f%%)' % (pcB + 1, fracs[1] * 100)
            PCsToInclude = 4
            correlated_db = {}
            allGenes = {}
            new_matrix = []
            new_headers = []
            added_indexes = []
            x = 0
            print 'exporting PCA loading genes to:', root_dir + '/PCA/correlated.txt'
            exportData = export.ExportFile(root_dir + '/PCA/correlated.txt')
            matrix = zip(*matrix)
            try:
                while x < PCsToInclude:
                    idx = numpy.argsort(u[:, x])
                    correlated = map(lambda i: row_header[i], idx[:300])
                    anticorrelated = map(lambda i: row_header[i], idx[-300:])
                    correlated_db[x] = (correlated, anticorrelated)
                    fidx = list(idx[:300]) + list(idx[-300:])
                    for i in fidx:
                        if i not in added_indexes:
                            added_indexes.append(i)
                            new_headers.append(row_header[i])
                            new_matrix.append(matrix[i])

                    x += 1

                redundant_genes = []
                for x in correlated_db:
                    correlated, anticorrelated = correlated_db[x]
                    count = 0
                    for gene in correlated:
                        if gene not in redundant_genes and count < 100:
                            exportData.write(gene + '\tcorrelated-PC' + str(x + 1) + '\n')
                            allGenes[gene] = []
                            count += 1

                    count = 0
                    for gene in anticorrelated:
                        if gene not in redundant_genes and count < 100:
                            exportData.write(gene + '\tanticorrelated-PC' + str(x + 1) + '\n')
                            allGenes[gene] = []
                            count += 1

                exportData.close()
                if geneSetName != None:
                    if len(geneSetName) > 0:
                        exportCustomGeneSet(geneSetName, species, allGenes)
                        print 'Exported geneset to "StoredGeneSets"'
            except Exception:
                pass

            if use_svd == False:
                latent, coeff = linalg.eig(cov(M))
                scores = dot(coeff.T, M)
            else:
                scores = vt
            writetSNEScores(scores, root_dir + dataset_name + '-PCA_scores.txt')
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        pylab.xlabel(label1)
        pylab.ylabel(label2)
        axes = getAxes(scores)
        pylab.axis(axes)
        marker_size = 15
        if len(column_header) > 20:
            marker_size = 12
        if len(column_header) > 40:
            marker_size = 10
        if len(column_header) > 150:
            marker_size = 7
        if len(column_header) > 500:
            marker_size = 5
        if len(column_header) > 1000:
            marker_size = 4
        if len(column_header) > 2000:
            marker_size = 3
        if colorByGene != None:
            print 'Coloring based on feature expression.'
            gene_translation_db = {}
            matrix = numpy.array(matrix)
            min_val = matrix.min()
            if ' ' in colorByGene:
                genes = string.split(colorByGene, ' ')
            else:
                genes = [
                 colorByGene]
            genePresent = False
            numberGenesPresent = []
            for gene in genes:
                if gene in row_header:
                    numberGenesPresent.append(gene)
                    genePresent = True

            if len(numberGenesPresent) == 0:
                try:
                    import gene_associations
                    from import_scripts import OBO_import
                    gene_to_symbol = gene_associations.getGeneToUid(species, ('hide',
                                                                              'Ensembl-Symbol'))
                    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
                    for symbol in genes:
                        if symbol in symbol_to_gene:
                            gene = symbol_to_gene[symbol][0]
                            if gene in row_header:
                                numberGenesPresent.append(gene)
                                genePresent = True
                                gene_translation_db[symbol] = gene

                except Exception:
                    pass

            numberGenesPresent = len(numberGenesPresent)
            if numberGenesPresent == 1:
                cm = pylab.cm.get_cmap('Reds')
            elif numberGenesPresent == 2:
                cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
            elif numberGenesPresent == 3:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif numberGenesPresent == 4:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif numberGenesPresent == 5:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent == 6:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent == 7:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            else:
                cm = pylab.cm.get_cmap('gist_rainbow')

            def get_cmap(n, name='hsv'):
                """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
            RGB color; the keyword argument name must be a standard mpl colormap name."""
                return pylab.cm.get_cmap(name, n)

            if genePresent:
                dataset_name += '-' + colorByGene
                group_db = {}
                bestGeneAssociated = {}
                k = 0
                for gene in genes:
                    try:
                        try:
                            i = row_header.index(gene)
                        except Exception:
                            i = row_header.index(gene_translation_db[gene])
                        else:
                            values = map(float, matrix[i])
                            min_val = min(values)
                            bin_size = (max(values) - min_val) / 8
                            max_val = max(values)
                            ranges = []
                            iz = min_val
                            while iz < max(values) - bin_size / 100:
                                r = (
                                 iz, iz + bin_size)
                                if len(ranges) == 7:
                                    r = (
                                     iz, max_val)
                                ranges.append(r)
                                iz += bin_size

                        color_db = {}
                        colors = get_cmap(len(genes))
                        for i in range(len(ranges)):
                            if i == 0:
                                color = '#C0C0C0'
                            elif numberGenesPresent == 1:
                                color = cm(1.0 * i / len(ranges))
                            elif i > 2:
                                if len(genes) < 8:
                                    color = cm(k)
                                else:
                                    color = colors(k)
                            else:
                                color = '#C0C0C0'
                            color_db[ranges[i]] = color

                        i = 0
                        for val in values:
                            sample = column_header[i]
                            for l, u in color_db:
                                range_index = ranges.index((l, u))
                                if val >= l and val <= u:
                                    color = color_db[(l, u)]
                                    color_label = [gene + '-range: ' + str(l)[:4] + '-' + str(u)[:4], color, '']
                                    group_db[sample] = color_label
                                    try:
                                        bestGeneAssociated[sample].append([range_index, val, color_label])
                                    except Exception:
                                        bestGeneAssociated[sample] = [[range_index, val, color_label]]

                            i += 1

                        if len(genes) > 1:
                            for sample in bestGeneAssociated:
                                bestGeneAssociated[sample].sort()
                                color_label = bestGeneAssociated[sample][(-1)][(-1)]
                                if numberGenesPresent > 1:
                                    index = bestGeneAssociated[sample][(-1)][0]
                                    if index > 2:
                                        gene = string.split(color_label[0], '-')[0]
                                    else:
                                        gene = 'Null'
                                    color_label[0] = gene
                                group_db[sample] = color_label

                    except Exception:
                        print [
                         gene], 'not found in rows...'
                    else:
                        k += 1

            else:
                print [
                 colorByGene], 'not found in rows...'
        pylab.title('Principal Component Analysis - ' + dataset_name)
        group_names = {}
        i = 0
        for sample_name in column_header:
            try:
                group_name, color, k = group_db[sample_name]
                if group_name not in group_names:
                    label = group_name
                else:
                    label = None
                group_names[group_name] = color
            except Exception:
                color = 'r'
                label = None
            else:
                try:
                    ax.plot(scores[pcA][i], scores[1][i], color=color, marker='o', markersize=marker_size, label=label, markeredgewidth=0, picker=True)
                except Exception as e:
                    print e
                    print i, len(scores[pcB])
                    kill
                else:
                    if showLabels:
                        try:
                            sample_name = '   ' + string.split(sample_name, ':')[1]
                        except Exception:
                            pass

                        ax.text(scores[pcA][i], scores[pcB][i], sample_name, fontsize=11)
                    i += 1

        group_count = []
        for i in group_db:
            if group_db[i][0] not in group_count:
                group_count.append(group_db[i][0])

        Lfontsize = 8
        if len(group_count) > 20:
            Lfontsize = 10
        if len(group_count) > 30:
            Lfontsize = 8
        if len(group_count) > 40:
            Lfontsize = 6
        if len(group_count) > 50:
            Lfontsize = 5
        i = 0
        box = ax.get_position()
        if len(group_count) > 0:
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            try:
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=Lfontsize)
            except Exception:
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        else:
            ax.set_position([box.x0, box.y0, box.width, box.height])
            pylab.legend(loc='upper left', prop={})
        filename = 'Clustering-%s-PCA.pdf' % dataset_name
        try:
            pylab.savefig(root_dir + filename)
        except Exception:
            None

        filename = filename[:-3] + 'png'
        try:
            pylab.savefig(root_dir + filename)
        except Exception:
            None

    graphic_link.append(['Principal Component Analysis', root_dir + filename])
    if display:
        print 'Exporting:', filename
        try:
            pylab.show()
        except Exception:
            pass

    fig.clf()
    return


def ViolinPlot():

    def readData(filename):
        all_data = {}
        headers = {}
        groups = []
        firstRow = True
        for line in open(filename, 'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data, '\t')
            if firstRow:
                firstRow = False
                i = 0
                for x in t[1:]:
                    try:
                        g, h = string.split(x, ':')
                    except Exception:
                        g = x
                        h = x
                    else:
                        headers[i] = g
                        if g not in groups:
                            groups.append(g)
                        i += 1

            else:
                t = map(lambda x: float(x), t[1:])
                i = 0
                for x in t:
                    try:
                        g = headers[i]
                    except Exception:
                        print i
                        sys.exit()
                    else:
                        try:
                            all_data[g].append(x)
                        except Exception:
                            all_data[g] = [x]
                        else:
                            i += 1

        all_data2 = []
        print groups
        for group in groups:
            all_data2.append(all_data[group])

        return all_data2

    def violin_plot(ax, data, pos, bp=False):
        """      
        create violin plots on an axis   
        """
        from scipy.stats import gaussian_kde
        from numpy import arange
        dist = max(pos) - min(pos)
        w = min(0.15 * max(dist, 1.0), 0.5)
        for d, p in zip(data, pos):
            k = gaussian_kde(d)
            m = k.dataset.min()
            M = k.dataset.max()
            x = arange(m, M, (M - m) / 100.0)
            v = k.evaluate(x)
            v = v / v.max() * w
            ax.fill_betweenx(x, p, v + p, facecolor='y', alpha=0.3)
            ax.fill_betweenx(x, p, -v + p, facecolor='y', alpha=0.3)

        if bp:
            ax.boxplot(data, notch=1, positions=pos, vert=1)

    def draw_all(data, output):
        pos = [1, 2, 3]
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        violin_plot(ax, data, pos)
        pylab.show()
        pylab.savefig(output + '.pdf')

    all_data = []
    all_data = readData('/Users/saljh8/Downloads/temp3.txt')
    import numpy
    draw_all(all_data, 'alldata')


def simpleScatter(fn):
    import matplotlib.patches as mpatches
    values = []
    legends = {}
    colors = {}
    skip = True
    scale = 100.0
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        if skip:
            x_header, y_header, color_header, label_header, shape_header = string.split(data, '\t')
            skip = False
        else:
            x, y, color, label, shape = string.split(data, '\t')
            if color in colors:
                xval, yval, label, shape = colors[color]
                xval.append(float(x))
                yval.append(float(y))
            else:
                xval = [
                 float(x)]
                yval = [float(y)]
                colors[color] = (
                 xval, yval, label, shape)

    for color in colors:
        xval, yval, label, shape = colors[color]
        pylab.scatter(xval, yval, s=scale, c=color, alpha=0.75, label=label, marker=shape, edgecolor='none')

    pylab.legend(loc='upper left')
    pylab.title(fn)
    pylab.xlabel(x_header, fontsize=15)
    pylab.ylabel(y_header, fontsize=15)
    marker_size = 7
    pylab.show()


def ica(filename):
    showLabels = True
    X, column_header, row_header, dataset_name, group_db = importData(filename)
    X = map(numpy.array, zip(*X))
    column_header, row_header = row_header, column_header
    ica = FastICA()
    scores = ica.fit(X).transform(X)
    scores /= scores.std(axis=0)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    pylab.xlabel('ICA-X')
    pylab.ylabel('ICA-Y')
    pylab.title('ICA - ' + dataset_name)
    axes = getAxes(scores)
    pylab.axis(axes)
    marker_size = 15
    if len(column_header) > 20:
        marker_size = 12
    if len(column_header) > 40:
        marker_size = 10
    if len(column_header) > 150:
        marker_size = 7
    if len(column_header) > 500:
        marker_size = 5
    if len(column_header) > 1000:
        marker_size = 4
    if len(column_header) > 2000:
        marker_size = 3
    group_names = {}
    i = 0
    for sample_name in row_header:
        try:
            group_name, color, k = group_db[sample_name]
            if group_name not in group_names:
                label = group_name
            else:
                label = None
            group_names[group_name] = color
        except Exception:
            color = 'r'
            label = None
        else:
            ax.plot(scores[0][i], scores[1][i], color=color, marker='o', markersize=marker_size, label=label)
            if showLabels:
                ax.text(scores[0][i], scores[1][i], sample_name, fontsize=8)
            i += 1

    pylab.title('ICA recovered signals')
    pylab.show()
    return


def plot_samples(S, axis_list=None):
    pylab.scatter(S[:, 0], S[:, 1], s=20, marker='o', linewidths=0, zorder=10, color='red', alpha=0.5)
    if axis_list is not None:
        colors = [
         'orange', 'red']
        for color, axis in zip(colors, axis_list):
            axis /= axis.std()
            x_axis, y_axis = axis
            pylab.plot(0.1 * x_axis, 0.1 * y_axis, linewidth=2, color=color)
            pylab.quiver(0, 0, x_axis, y_axis, zorder=11, width=2, scale=6, color=color)

    pylab.xlabel('x')
    pylab.ylabel('y')
    return


def PCA3D(matrix, column_header, row_header, dataset_name, group_db, display=False, showLabels=True, algorithm='SVD', geneSetName=None, species=None, colorByGene=None):
    from numpy import mean, cov, double, cumsum, dot, linalg, array, rank
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    start = time.time()
    try:
        prior_clusters = priorColumnClusters
    except Exception:
        prior_clusters = []
    else:
        if prior_clusters == None:
            prior_clusters = []
        try:
            if len(prior_clusters) > 0 and len(group_db) == 0:
                newColumnHeader = []
                i = 0
                for sample_name in column_header:
                    newColumnHeader.append(str(prior_clusters[i]) + ':' + sample_name)
                    i += 1

                group_db, column_header = assignGroupColors(newColumnHeader)
        except Exception as e:
            group_db = {}

    if algorithm == 'SVD':
        use_svd = True
    else:
        use_svd = False
    Mdif = matrix / matrix.std()
    Mdif = Mdif.T
    u, s, vt = svd(Mdif, 0)
    fracs = s ** 2 / np.sum(s ** 2)
    entropy = -sum(fracs * np.log(fracs)) / np.log(np.min(vt.shape))
    label1 = 'PC%i (%2.1f%%)' % (1, fracs[0] * 100)
    label2 = 'PC%i (%2.1f%%)' % (2, fracs[1] * 100)
    label3 = 'PC%i (%2.1f%%)' % (3, fracs[2] * 100)
    PCsToInclude = 4
    correlated_db = {}
    allGenes = {}
    new_matrix = []
    new_headers = []
    added_indexes = []
    x = 0
    print 'exporting PCA loading genes to:', root_dir + '/PCA/correlated.txt'
    exportData = export.ExportFile(root_dir + '/PCA/correlated.txt')
    matrix = zip(*matrix)
    try:
        while x < PCsToInclude:
            idx = numpy.argsort(u[:, x])
            correlated = map(lambda i: row_header[i], idx[:300])
            anticorrelated = map(lambda i: row_header[i], idx[-300:])
            correlated_db[x] = (correlated, anticorrelated)
            fidx = list(idx[:300]) + list(idx[-300:])
            for i in fidx:
                if i not in added_indexes:
                    added_indexes.append(i)
                    new_headers.append(row_header[i])
                    new_matrix.append(matrix[i])

            x += 1

        redundant_genes = []
        for x in correlated_db:
            correlated, anticorrelated = correlated_db[x]
            count = 0
            for gene in correlated:
                if gene not in redundant_genes and count < 100:
                    exportData.write(gene + '\tcorrelated-PC' + str(x + 1) + '\n')
                    allGenes[gene] = []
                    count += 1

            count = 0
            for gene in anticorrelated:
                if gene not in redundant_genes and count < 100:
                    exportData.write(gene + '\tanticorrelated-PC' + str(x + 1) + '\n')
                    allGenes[gene] = []
                    count += 1

        exportData.close()
        if geneSetName != None:
            if len(geneSetName) > 0:
                exportCustomGeneSet(geneSetName, species, allGenes)
                print 'Exported geneset to "StoredGeneSets"'
    except ZeroDivisionError:
        pass
    else:
        if use_svd == False:
            latent, coeff = linalg.eig(cov(M))
            scores = dot(coeff.T, M)
        else:
            scores = vt
        end = time.time()
        print 'PCA completed in', end - start, 'seconds.'
        ax.set_xlabel(label1)
        ax.set_ylabel(label2)
        ax.set_zlabel(label3)
        axes = getAxes(scores, PlotType='3D')
        pylab.axis(axes)
        Lfontsize = 8
        group_count = []
        for i in group_db:
            if group_db[i][0] not in group_count:
                group_count.append(group_db[i][0])

        if colorByGene != None:
            gene_translation_db = {}
            matrix = numpy.array(matrix)
            min_val = matrix.min()
            if ' ' in colorByGene:
                genes = string.split(colorByGene, ' ')
            else:
                genes = [
                 colorByGene]
            genePresent = False
            numberGenesPresent = []
            for gene in genes:
                if gene in row_header:
                    numberGenesPresent.append(gene)
                    genePresent = True

            if len(numberGenesPresent) == 0:
                try:
                    import gene_associations
                    from import_scripts import OBO_import
                    gene_to_symbol = gene_associations.getGeneToUid(species, ('hide',
                                                                              'Ensembl-Symbol'))
                    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
                    for symbol in genes:
                        if symbol in symbol_to_gene:
                            gene = symbol_to_gene[symbol][0]
                            if gene in row_header:
                                numberGenesPresent.append(gene)
                                genePresent = True
                                gene_translation_db[symbol] = gene

                except Exception:
                    pass

            numberGenesPresent = len(numberGenesPresent)
            if numberGenesPresent == 1:
                cm = pylab.cm.get_cmap('Reds')
            elif numberGenesPresent == 2:
                cm = matplotlib.colors.ListedColormap(['#00FF00', '#1E90FF'])
                cm = matplotlib.colors.ListedColormap(['w', 'k'])
            elif numberGenesPresent == 3:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C'])
            elif numberGenesPresent == 4:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#3D3181', '#EE2C3C', '#FEBC18'])
            elif numberGenesPresent == 5:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#3D3181', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent == 6:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            elif numberGenesPresent == 7:
                cm = matplotlib.colors.ListedColormap(['#88BF47', '#63C6BB', '#29C3EC', '#3D3181', '#7B4976', '#FEBC18', '#EE2C3C'])
            else:
                cm = pylab.cm.get_cmap('gist_rainbow')
            if genePresent:
                dataset_name += '-' + colorByGene
                group_db = {}
                bestGeneAssociated = {}
                k = 0
                for gene in genes:
                    try:
                        try:
                            i = row_header.index(gene)
                        except Exception:
                            i = row_header.index(gene_translation_db[gene])
                        else:
                            values = map(float, matrix[i])
                            min_val = min(values)
                            bin_size = (max(values) - min_val) / 8
                            max_val = max(values)
                            ranges = []
                            iz = min_val
                            while iz < max(values) - bin_size / 100:
                                r = (
                                 iz, iz + bin_size)
                                if len(ranges) == 7:
                                    r = (
                                     iz, max_val)
                                ranges.append(r)
                                iz += bin_size

                        color_db = {}
                        for i in range(len(ranges)):
                            if i == 0:
                                color = '#C0C0C0'
                            elif numberGenesPresent == 1:
                                color = cm(1.0 * i / len(ranges))
                            elif i > 2:
                                color = cm(k)
                            else:
                                color = '#C0C0C0'
                            color_db[ranges[i]] = color

                        i = 0
                        for val in values:
                            sample = column_header[i]
                            for l, u in color_db:
                                range_index = ranges.index((l, u))
                                if val >= l and val <= u:
                                    color = color_db[(l, u)]
                                    color_label = [gene + '-range: ' + str(l)[:4] + '-' + str(u)[:4], color, '']
                                    group_db[sample] = color_label
                                    try:
                                        bestGeneAssociated[sample].append([range_index, val, color_label])
                                    except Exception:
                                        bestGeneAssociated[sample] = [[range_index, val, color_label]]

                            i += 1

                        if len(genes) > 1:
                            for sample in bestGeneAssociated:
                                bestGeneAssociated[sample].sort()
                                color_label = bestGeneAssociated[sample][(-1)][(-1)]
                                if numberGenesPresent > 1:
                                    index = bestGeneAssociated[sample][(-1)][0]
                                    if index > 2:
                                        gene = string.split(color_label[0], '-')[0]
                                    else:
                                        gene = 'Null'
                                    color_label[0] = gene
                                group_db[sample] = color_label

                    except Exception:
                        print [
                         gene], 'not found in rows...'
                    else:
                        k += 1

            else:
                print [
                 colorByGene], 'not found in rows...'
        if len(group_count) > 20:
            Lfontsize = 10
        if len(group_count) > 30:
            Lfontsize = 8
        if len(group_count) > 40:
            Lfontsize = 6
        if len(group_count) > 50:
            Lfontsize = 5
        if len(scores[0]) > 150:
            markersize = 7
        else:
            markersize = 10
        i = 0
        group_names = {}
        for x in scores[0]:
            sample_name = column_header[i]
            try:
                group_name, color, k = group_db[sample_name]
                if group_name not in group_names:
                    label = group_name
                else:
                    label = None
                group_names[group_name] = (
                 color, k)
            except Exception:
                color = 'r'
                label = None
            else:
                ax.plot([scores[0][i]], [scores[1][i]], [scores[2][i]], color=color, marker='o', markersize=markersize, label=label, markeredgewidth=0, picker=True)
                if showLabels:
                    ax.text(scores[0][i], scores[1][i], scores[2][i], '   ' + sample_name, fontsize=9)
                i += 1

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        try:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=Lfontsize)
        except Exception:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    filename = 'Clustering-%s-3D-PCA.pdf' % dataset_name
    pylab.savefig(root_dir + filename)
    filename = filename[:-3] + 'png'
    pylab.savefig(root_dir + filename)
    graphic_link.append(['Principal Component Analysis', root_dir + filename])
    if display:
        print 'Exporting:', filename
        try:
            pylab.show()
        except Exception:
            None

    fig.clf()
    return


def getAxes1(scores, PlotType=None):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    try:
        x_range = max(scores[0]) - min(scores[0])
        y_range = max(scores[1]) - min(scores[1])
        if PlotType == '3D':
            x_axis_min = min(scores[0]) - x_range / 10
            x_axis_max = max(scores[0]) + x_range / 10
            y_axis_min = min(scores[1]) - y_range / 10
            y_axis_max = max(scores[1]) + y_range / 10
        else:
            x_axis_min = min(scores[0]) - x_range / 10
            x_axis_max = max(scores[0]) + x_range / 10
            y_axis_min = min(scores[1]) - y_range / 10
            y_axis_max = max(scores[1]) + y_range / 10
    except KeyError:
        None

    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]


def getAxes(scores, PlotType=None):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    try:
        x_range = max(scores[0]) - min(scores[0])
        y_range = max(scores[1]) - min(scores[1])
        if PlotType == '3D':
            x_axis_min = min(scores[0]) - x_range / 1.5
            x_axis_max = max(scores[0]) + x_range / 1.5
            y_axis_min = min(scores[1]) - y_range / 5
            y_axis_max = max(scores[1]) + y_range / 5
        else:
            x_axis_min = min(scores[0]) - x_range / 10
            x_axis_max = max(scores[0]) + x_range / 10
            y_axis_min = min(scores[1]) - y_range / 10
            y_axis_max = max(scores[1]) + y_range / 10
    except KeyError:
        None

    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]


def getAxesTransposed(scores, exclude={}):
    """ Adjust these axes to account for (A) legend size (left hand upper corner)
    and (B) long sample name extending to the right
    """
    scores_filtered = []
    for i in range(len(scores)):
        if i not in exclude:
            scores_filtered.append(scores[i])

    scores = scores_filtered
    scores = map(numpy.array, zip(*scores))
    try:
        x_range = max(scores[0]) - min(scores[0])
        y_range = max(scores[1]) - min(scores[1])
        x_axis_min = min(scores[0]) - int(float(x_range) / 7)
        x_axis_max = max(scores[0]) + int(float(x_range) / 7)
        y_axis_min = min(scores[1]) - int(float(y_range / 7))
        y_axis_max = max(scores[1]) + int(float(y_range / 7))
    except KeyError:
        None

    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]


def Kmeans(features, column_header, row_header):
    features = numpy.vstack((class1, class2))
    centroids, variance = scipy.cluster.vq.kmeans(features, 2)
    code, distance = scipy.cluster.vq.vq(features, centroids)
    pylab.plot([ p[0] for p in class1 ], [ p[1] for p in class1 ], '*')
    pylab.plot([ p[0] for p in class2 ], [ p[1] for p in class2 ], 'r*')
    pylab.plot([ p[0] for p in centroids ], [ p[1] for p in centroids ], 'go')
    pylab.show()


def findParentDir(filename):
    filename = string.replace(filename, '//', '/')
    filename = string.replace(filename, '\\', '/')
    x = string.find(filename[::-1], '/') * -1
    return filename[:x]


def findFilename(filename):
    filename = string.replace(filename, '//', '/')
    filename = string.replace(filename, '\\', '/')
    x = string.find(filename[::-1], '/') * -1
    return filename[x:]


def runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, display=False, contrast=None, allowAxisCompression=True, Normalize=True):
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
        if allowLargeClusters:
            maxSize = 50000
        else:
            maxSize = 7000
    except Exception:
        maxSize = 7000
    else:
        try:
            PriorColumnClusters = priorColumnClusters
            PriorRowClusters = priorRowClusters
        except Exception:
            PriorColumnClusters = None
            PriorRowClusters = None

        run = False
        print 'max allowed cluster size:', maxSize
        if len(matrix) > 0 and (len(matrix) < maxSize or row_method == None):
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=UserWarning)
                try:
                    heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name, display=display, contrast=contrast, allowAxisCompression=allowAxisCompression, Normalize=Normalize, PriorColumnClusters=PriorColumnClusters, PriorRowClusters=PriorRowClusters)
                    run = True
                except Exception:
                    print traceback.format_exc()
                    try:
                        pylab.clf()
                        pylab.close()
                        import gc
                        gc.collect()
                    except Exception:
                        None
                    else:
                        if len(matrix) < 10000:
                            print 'Error using %s ... trying euclidean instead' % row_metric
                            row_metric = 'cosine'
                            row_method = 'average'
                        else:
                            print 'Error with hierarchical clustering... only clustering arrays'
                            row_method = None
                        try:
                            heatmap(numpy.array(matrix), row_header, column_header, row_method, column_method, row_metric, column_metric, color_gradient, dataset_name, display=display, contrast=contrast, allowAxisCompression=allowAxisCompression, Normalize=Normalize, PriorColumnClusters=PriorColumnClusters, PriorRowClusters=PriorRowClusters)
                            run = True
                        except Exception:
                            print traceback.format_exc()
                            print 'Unable to generate cluster due to dataset incompatibilty.'

        elif len(matrix) == 0:
            print_out = 'SKIPPING HIERARCHICAL CLUSTERING!!! - Your dataset file has no associated rows.'
            print print_out
        else:
            print_out = 'SKIPPING HIERARCHICAL CLUSTERING!!! - Your dataset file is over the recommended size limit for clustering (' + str(maxSize) + ' rows). Please cluster later using "Additional Analyses"'
            print print_out
        try:
            pylab.clf()
            pylab.close()
            import gc
            gc.collect()
        except Exception:
            None

    return run


def debugTKBug():
    return


def runHCexplicit(filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient, extra_params, display=True, contrast=None, Normalize=False, JustShowTheseIDs=[], compressAxis=True):
    """ Explicit method for hieararchical clustering with defaults defined by the user (see below function) """
    global EliteGeneSets
    global GroupDB
    global allowLargeClusters
    global graphic_link
    global inputFilename
    global justShowTheseIDs
    global normalize
    global originalFilename
    global rho_cutoff
    global root_dir
    global runGOElite
    global species
    global storeGeneSetName
    global targetGeneIDs
    EliteGeneSets = []
    targetGene = []
    filterByPathways = False
    runGOElite = False
    justShowTheseIDs = JustShowTheseIDs
    allowLargeClusters = True
    if compressAxis:
        allowAxisCompression = True
    else:
        allowAxisCompression = False
    graphic_link = graphics
    inputFilename = filename
    filterIDs = False
    normalize = Normalize
    try:
        transpose = extra_params.Transpose()
        try:
            rho_cutoff = extra_params.RhoCutoff()
            print 'Setting correlation cutoff to a rho of', rho_cutoff
        except Exception:
            rho_cutoff = 0.5

        PathwayFilter = extra_params.PathwaySelect()
        GeneSet = extra_params.GeneSet()
        OntologyID = extra_params.OntologyID()
        Normalize = extra_params.Normalize()
        normalize = Normalize
        filterIDs = True
        species = extra_params.Species()
        platform = extra_params.Platform()
        vendor = extra_params.Vendor()
        newInput = findParentDir(inputFilename) + '/GeneSetClustering/' + findFilename(inputFilename)
        targetGene = extra_params.GeneSelection()
        getGeneCorrelations = extra_params.GetGeneCorrelations()
        filterByPathways = extra_params.FilterByPathways()
        PathwayFilter, filterByPathways = verifyPathwayName(PathwayFilter, GeneSet, OntologyID, filterByPathways)
        justShowTheseIDs_var = extra_params.JustShowTheseIDs()
        if len(justShowTheseIDs_var) > 0:
            justShowTheseIDs = justShowTheseIDs_var
        elif len(targetGene) > 0:
            targetGene = string.replace(targetGene, '\n', ' ')
            targetGene = string.replace(targetGene, '\r', ' ')
            justShowTheseIDs = string.split(targetGene, ' ')
        else:
            try:
                EliteGeneSets = extra_params.ClusterGOElite()
                if EliteGeneSets != ['']:
                    runGOElite = True
            except Exception:
                pass

            try:
                storeGeneSetName = extra_params.StoreGeneSetName()
            except Exception:
                storeGeneSetName = ''

    except Exception as e:
        transpose = extra_params
    else:
        root_dir = findParentDir(filename)
        if 'ExpressionOutput/Clustering' in root_dir:
            root_dir = string.replace(root_dir, 'ExpressionOutput/Clustering', 'DataPlots')
        elif 'ExpressionOutput' in root_dir:
            root_dir = string.replace(root_dir, 'ExpressionOutput', 'DataPlots')
            root_dir = string.replace(root_dir, '/Clustering', '')
        else:
            root_dir += '/DataPlots/'
            try:
                os.mkdir(root_dir)
            except Exception:
                None

    if row_method == 'hopach':
        reverseOrder = False
    else:
        reverseOrder = True
    matrix, column_header, row_header, dataset_name, group_db = importData(filename, Normalize=Normalize, reverseOrder=reverseOrder)
    GroupDB = group_db
    inputFilename = string.replace(inputFilename, '.cdt', '.txt')
    originalFilename = inputFilename
    if len(justShowTheseIDs) == 0:
        try:
            if len(priorColumnClusters) > 0 and priorRowClusters > 0 and row_method == None and column_method == None:
                try:
                    justShowTheseIDs = importPriorDrivers(inputFilename)
                except Exception:
                    pass

        except Exception:
            pass

    if filterIDs:
        transpose_update = True
        if filterByPathways:
            if isinstance(PathwayFilter, tuple) or isinstance(PathwayFilter, list):
                FileName = string.join(list(PathwayFilter), ' ')
                FileName = string.replace(FileName, ':', '-')
            else:
                FileName = PathwayFilter
            if len(FileName) > 40:
                FileName = FileName[:40]
            try:
                inputFilename = string.replace(newInput, '.txt', '_' + FileName + '.txt')
            except Exception:
                inputFilename = string.replace(newInput, '.txt', '_GeneSets.txt')
            else:
                vars = filterByPathway(matrix, row_header, column_header, species, platform, vendor, GeneSet, PathwayFilter, OntologyID, transpose)
                try:
                    dataset_name += '-' + FileName
                except Exception:
                    dataset_name += '-GeneSets'

            transpose_update = False
            if 'amplify' in targetGene:
                targetGene = string.join(vars[1], ' ') + ' amplify ' + targetGene
            else:
                matrix, row_header, column_header = vars
        try:
            alt_targetGene = string.replace(targetGene, 'amplify', '')
            alt_targetGene = string.replace(alt_targetGene, 'amplify', '')
            alt_targetGene = string.replace(alt_targetGene, 'driver', '')
            alt_targetGene = string.replace(alt_targetGene, 'guide', '')
            alt_targetGene = string.replace(alt_targetGene, 'top', '')
            alt_targetGene = string.replace(alt_targetGene, 'positive', '')
            alt_targetGene = string.replace(alt_targetGene, 'excludeCellCycle', '')
            alt_targetGene = string.replace(alt_targetGene, 'monocle', '')
            alt_targetGene = string.replace(alt_targetGene, 'GuideOnlyCorrelation', '')
            alt_targetGene = string.replace(alt_targetGene, ' ', '')
        except Exception:
            alt_targetGene = ''

        if getGeneCorrelations and targetGene != 'driver' and targetGene != 'GuideOnlyCorrelation' and targetGene != 'guide' and targetGene != 'excludeCellCycle' and targetGene != 'top' and targetGene != ' monocle' and targetGene != 'positive' and len(alt_targetGene) > 0:
            allowAxisCompression = False
            if transpose and transpose_update == False:
                transpose_update = False
            elif transpose and transpose_update:
                transpose_update = True
            else:
                transpose_update = False
            if '\r' in targetGene or '\n' in targetGene:
                targetGene = string.replace(targetGene, '\r', ' ')
                targetGene = string.replace(targetGene, '\n', ' ')
            if len(targetGene) > 15:
                inputFilename = string.replace(newInput, '.txt', '-' + targetGene[:50] + '.txt')
                dataset_name += '-' + targetGene[:50]
            else:
                inputFilename = string.replace(newInput, '.txt', '-' + targetGene + '.txt')
                dataset_name += '-' + targetGene
            inputFilename = root_dir + '/' + string.replace(findFilename(inputFilename), '|', ' ')
            inputFilename = root_dir + '/' + string.replace(findFilename(inputFilename), ':', ' ')
            dataset_name = string.replace(dataset_name, '|', ' ')
            dataset_name = string.replace(dataset_name, ':', ' ')
            try:
                matrix, row_header, column_header, row_method = getAllCorrelatedGenes(matrix, row_header, column_header, species, platform, vendor, targetGene, row_method, transpose_update)
            except Exception:
                print traceback.format_exc()
                print targetGene, 'not found in input expression file. Exiting. \n\n'
                badExit

            targetGeneIDs = targetGene
            exportTargetGeneList(targetGene, inputFilename)
    elif transpose:
        print 'Transposing the data matrix'
        matrix = map(numpy.array, zip(*matrix))
        column_header, row_header = row_header, column_header
    if len(column_header) > 1000 or len(row_header) > 1000:
        print 'Performing hierarchical clustering (please be patient)...'
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, display=display, contrast=contrast, allowAxisCompression=allowAxisCompression, Normalize=Normalize)
    if 'guide' in targetGene:
        import RNASeq
        input_file = graphic_link[(-1)][(-1)][:-4] + '.txt'
        if 'excludeCellCycle' in targetGene:
            excludeCellCycle = True
        else:
            excludeCellCycle = False
        print 'excludeCellCycle', excludeCellCycle
        targetGene = RNASeq.remoteGetDriverGenes(species, platform, input_file, excludeCellCycle=excludeCellCycle, ColumnMethod=column_method)
        extra_params.setGeneSelection(targetGene)
        extra_params.setGeneSet('None Selected')
        graphic_link = runHCexplicit(filename, graphic_link, row_method, row_metric, column_method, column_metric, color_gradient, extra_params, display=display, contrast=contrast, Normalize=Normalize, JustShowTheseIDs=JustShowTheseIDs, compressAxis=compressAxis)
    return graphic_link


def importPriorDrivers(inputFilename):
    filename = string.replace(inputFilename, 'Clustering-', '')
    filename = string.split(filename, '-hierarchical')[0] + '-targetGenes.txt'
    genes = open(filename, 'rU')
    genes = map(lambda x: cleanUpLine(x), genes)
    return genes


def exportTargetGeneList(targetGene, inputFilename):
    exclude = [
     'positive', 'top', 'driver', 'guide', 'amplify', 'GuideOnlyCorrelation']
    exportFile = inputFilename[:-4] + '-targetGenes.txt'
    eo = export.ExportFile(root_dir + findFilename(exportFile))
    targetGenes = string.split(targetGene, ' ')
    for gene in targetGenes:
        if gene not in exclude:
            try:
                eo.write(gene + '\n')
            except Exception:
                print 'Error export out gene (bad ascii):', [gene]

    eo.close()


def debugPylab():
    pylab.figure()
    pylab.close()
    pylab.figure()


def verifyPathwayName(PathwayFilter, GeneSet, OntologyID, filterByPathways):
    import gene_associations
    if len(OntologyID) > 0:
        PathwayFilter = gene_associations.lookupOntologyID(GeneSet, OntologyID, type='ID')
        filterByPathways = True
    return (PathwayFilter, filterByPathways)


def filterByPathway(matrix, row_header, column_header, species, platform, vendor, GeneSet, PathwayFilter, OntologyID, transpose):
    import gene_associations
    from import_scripts import OBO_import
    exportData = export.ExportFile(inputFilename)
    matrix2 = []
    row_header2 = []
    if 'Ontology' in GeneSet:
        directory = 'nested'
    else:
        directory = 'gene-mapp'
    print 'GeneSet(s) to analyze:', PathwayFilter
    if isinstance(PathwayFilter, tuple) or isinstance(PathwayFilter, list):
        associated_IDs = {}
        for p in PathwayFilter:
            associated = gene_associations.simpleGenePathwayImport(species, GeneSet, p, OntologyID, directory)
            for i in associated:
                associated_IDs[i] = []

    else:
        associated_IDs = gene_associations.simpleGenePathwayImport(species, GeneSet, PathwayFilter, OntologyID, directory)
    gene_annotations = gene_associations.importGeneData(species, 'Ensembl')
    vendor = string.replace(vendor, 'other:', '')
    try:
        array_to_ens = gene_associations.filterGeneToUID(species, 'Ensembl', vendor, associated_IDs)
    except Exception:
        array_to_ens = {}
    else:
        if platform == "3'array":
            try:
                array_to_ens = gene_associations.filterGeneToUID(species, 'Ensembl', vendor, associated_IDs)
            except Exception:
                pass

        try:
            gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', 'Ensembl-Symbol'))
            symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
        except Exception:
            pass

    i = 0
    original_rows = {}
    for row_id in row_header:
        original_id = row_id
        symbol = row_id
        if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename or 'AltConfirmed' in inputFilename or 'MarkerGenes' in inputFilename or 'blah' not in inputFilename:
            try:
                row_id, symbol = string.split(row_id, ' ')[:2]
            except Exception:
                try:
                    symbol = gene_to_symbol[row_id][0]
                except Exception:
                    None

            else:
                if len(symbol) == 0:
                    symbol = row_id
                if ':' in row_id:
                    try:
                        cluster, row_id = string.split(row_id, ':')
                        updated_row_id = cluster + ':' + symbol
                    except Exception:
                        pass

                else:
                    updated_row_id = symbol
                try:
                    original_id = updated_row_id
                except Exception:
                    pass

        if platform == "3'array":
            try:
                try:
                    row_ids = array_to_ens[row_id]
                except Exception:
                    row_ids = symbol_to_gene[symbol]

            except Exception:
                row_ids = [
                 row_id]

        else:
            try:
                try:
                    row_ids = array_to_ens[row_id]
                except Exception:
                    row_ids = symbol_to_gene[symbol]

            except Exception:
                row_ids = [
                 row_id]

        for row_id in row_ids:
            if row_id in associated_IDs:
                if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename:
                    if original_id != symbol:
                        row_id = original_id + ' ' + symbol
                    else:
                        row_id = symbol
                else:
                    try:
                        row_id = gene_annotations[row_id].Symbol()
                    except Exception:
                        None

                if original_id not in original_rows:
                    matrix2.append(matrix[i])
                    row_header2.append(original_id)
                    original_rows[original_id] = None

        i += 1

    if transpose:
        matrix2 = map(numpy.array, zip(*matrix2))
        column_header, row_header2 = row_header2, column_header
    exportData.write(string.join(['UID'] + column_header, '\t') + '\n')
    i = 0
    for row_id in row_header2:
        exportData.write(string.join([row_id] + map(str, matrix2[i]), '\t') + '\n')
        i += 1

    print len(row_header2), 'filtered IDs'
    exportData.close()
    return (
     matrix2, row_header2, column_header)


def getAllCorrelatedGenes(matrix, row_header, column_header, species, platform, vendor, targetGene, row_method, transpose):
    resort_by_ID_name = False
    if resort_by_ID_name:
        index = 0
        new_row_header = []
        new_matrix = []
        temp_row_header = []
        for name in row_header:
            temp_row_header.append((name, index))
            index += 1

        temp_row_header.sort()
        for name, index in temp_row_header:
            new_row_header.append(name)
            new_matrix.append(matrix[index])

        matrix = new_matrix
        row_header = new_row_header
    exportData = export.ExportFile(inputFilename)
    try:
        import gene_associations
        gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', 'Ensembl-Symbol'))
    except Exception:
        print 'No Ensembl-Symbol database available for', species

    if platform == "3'array":
        try:
            if ':' in vendor:
                vendor = string.split(vendor, ':')[1]
            array_to_ens = gene_associations.filterGeneToUID(species, 'Ensembl', vendor, {})
        except Exception as e:
            array_to_ens = {}

        for uid in array_to_ens:
            for gid in array_to_ens[uid]:
                if gid in gene_to_symbol:
                    symbol = gene_to_symbol[gid][0]
                    try:
                        gene_to_symbol[uid].append(symbol)
                    except Exception:
                        gene_to_symbol[uid] = [symbol]

    matrix2 = []
    row_header2 = []
    matrix_db = {}
    multipleGenes = False
    intersecting_ids = []
    i = 0
    targetGenes = [
     targetGene]
    if ' ' in targetGene or ',' in targetGene or '|' in targetGene or '\n' in targetGene or '\r' in targetGene:
        multipleGenes = True
        if '  ' in targetGene:
            targetGene = string.replace(targetGene, '  ', ' ')
        if ',' in targetGene:
            targetGene = string.replace(targetGene, ',', ' ')
        if '\n' in targetGene:
            targetGene = string.replace(targetGene, '\n', ' ')
        if '\r' in targetGene:
            targetGene = string.replace(targetGene, '\r', ' ')
        targetGenes = string.split(targetGene, ' ')
        if row_method != None:
            targetGenes.sort()
        intersecting_ids = [ val for val in targetGenes if val in row_header ]
        for row_id in row_header:
            original_rowid = row_id
            symbol = row_id
            new_symbol = symbol
            rigorous_search = True
            if ':' in row_id and '|' in row_id:
                rigorous_search = False
            elif ':' in row_id and '|' not in row_id:
                a, b = string.split(row_id, ':')[:2]
                if 'ENS' in a or len(a) == 17:
                    try:
                        row_id = a
                        symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        symbol = ''

                elif 'ENS' not in b and len(a) != 17:
                    row_id = b
                elif 'ENS' in b:
                    symbol = original_rowid
                    row_id = a
            if rigorous_search:
                try:
                    row_id, symbol = string.split(row_id, ' ')[:2]
                except Exception:
                    try:
                        symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        if 'ENS' not in original_rowid:
                            row_id, symbol = row_id, row_id

                    new_symbol = symbol

                if 'ENS' not in original_rowid and len(original_rowid) != 17:
                    if original_rowid != symbol:
                        symbol = original_rowid + ' ' + symbol
                for gene in targetGenes:
                    if string.lower(gene) == string.lower(row_id) or string.lower(gene) == string.lower(symbol) or string.lower(original_rowid) == string.lower(gene) or string.lower(gene) == string.lower(new_symbol):
                        matrix2.append(matrix[i])
                        row_header2.append(symbol)
                        matrix_db[symbol] = matrix[i]

            elif row_id in targetGenes:
                matrix2.append(matrix[i])
                row_header2.append(row_id)
                matrix_db[row_id] = matrix[i]
            i += 1

        i = 0
    else:
        i = 0
        original_rows = {}
        for row_id in row_header:
            original_id = row_id
            symbol = 'NA'
            if 'SampleLogFolds' in inputFilename or 'RelativeLogFolds' in inputFilename or 'blah' not in inputFilename:
                try:
                    row_id, symbol = string.split(row_id, ' ')[:2]
                except Exception:
                    try:
                        symbol = gene_to_symbol[row_id][0]
                    except Exception:
                        row_id, symbol = row_id, row_id

                original_id = row_id
            if row_id == targetGene or symbol == targetGene:
                targetGeneValues = matrix[i]
                break
            i += 1

        i = 0
    if multipleGenes == False:
        limit = 50
    else:
        limit = 140
    print 'limit:', limit
    if multipleGenes == False or 'amplify' in targetGene or 'correlated' in targetGene:
        row_header3 = []
        if multipleGenes == False:
            targetGeneValue_array = [
             targetGeneValues]
        else:
            targetGeneValue_array = matrix2
            if len(row_header2) > 4 or len(row_header) > 15000:
                print 'Performing all iterative pairwise corelations...',
                corr_matrix = numpyCorrelationMatrixGene(matrix, row_header, row_header2, gene_to_symbol)
                print 'complete'
            matrix2 = []
            original_headers = row_header2
            row_header2 = []
            matrix2_alt = []
            row_header2_alt = []
        import markerFinder
        k = 0
        for targetGeneValues in targetGeneValue_array:
            correlated = []
            anticorrelated = []
            try:
                targetGeneID = original_headers[k]
            except Exception:
                targetGeneID = ''

            try:
                rho_results = list(corr_matrix[targetGeneID])
            except Exception:
                rho_results = markerFinder.simpleScipyPearson(matrix, targetGeneValues)
            else:
                correlated_symbols = {}
                for rho, ind in rho_results[:limit]:
                    proceed = True
                    try:
                        if len(rho) == 2:
                            rho = rho[0]
                    except:
                        pass
                    else:
                        if 'top' in targetGene:
                            if rho_results[4][0] < rho_cutoff:
                                proceed = False
                        if rho > rho_cutoff and proceed:
                            rh = row_header[ind]
                            if len(row_header2) < 100 or multipleGenes:
                                rh = row_header[ind]
                                if matrix[ind] not in matrix2:
                                    if 'correlated' in targetGene:
                                        if rho != 1:
                                            matrix2.append(matrix[ind])
                                            row_header2.append(rh)
                                            if targetGeneValues not in matrix2:
                                                matrix2.append(targetGeneValues)
                                                row_header2.append(targetGeneID)
                                                try:
                                                    correlated_symbols[gene_to_symbol[rh][0]] = ind
                                                except Exception:
                                                    correlated_symbols[rh] = ind

                                    else:
                                        matrix2.append(matrix[ind])
                                        row_header2.append(rh)
                                        try:
                                            correlated_symbols[gene_to_symbol[rh][0]] = ind
                                        except Exception:
                                            correlated_symbols[rh] = ind

                rho_results.reverse()
                for rho, ind in rho_results[:limit]:
                    try:
                        if len(rho) == 2:
                            rho = rho[0]
                    except:
                        pass
                    else:
                        if rho < -1 * rho_cutoff and 'positive' not in targetGene:
                            rh = row_header[ind]
                            if len(row_header2) < 100 or multipleGenes:
                                rh = row_header[ind]
                                if matrix[ind] not in matrix2:
                                    if 'correlated' in targetGene:
                                        if rho != 1:
                                            matrix2.append(matrix[ind])
                                            row_header2.append(rh)
                                            if targetGeneValues not in matrix2:
                                                matrix2.append(targetGeneValues)
                                                row_header2.append(targetGeneID)
                                                try:
                                                    correlated_symbols[gene_to_symbol[rh][0]] = ind
                                                except Exception:
                                                    correlated_symbols[rh] = ind

                                    else:
                                        matrix2.append(matrix[ind])
                                        row_header2.append(rh)
                                        try:
                                            correlated_symbols[gene_to_symbol[rh][0]] = ind
                                        except Exception:
                                            correlated_symbols[rh] = ind

                try:
                    if len(correlated_symbols) > 0:
                        potentially_redundant = []
                        for i in targetGenes:
                            if i in correlated_symbols:
                                if i != targetGeneID:
                                    potentially_redundant.append((i, correlated_symbols[i]))

                        if len(potentially_redundant) > 0:
                            for rh, ind in potentially_redundant:
                                matrix2_alt.append(matrix[ind])
                                row_header2_alt.append(rh)

                        rho_results.reverse()
                except Exception:
                    pass
                else:
                    k += 1

        if 'IntraCorrelatedOnly' in targetGene:
            matrix2 = matrix2_alt
            row_header2 = row_header2_alt
        for r in row_header2:
            try:
                row_header3.append(gene_to_symbol[r][0])
            except Exception:
                row_header3.append(r)

        row_header2 = row_header3
        matrix2.reverse()
        row_header2.reverse()
        if 'amplify' not in targetGene:
            row_method = None
    if 'amplify' not in targetGene and 'correlated' not in targetGene:
        matrix_temp = []
        header_temp = []
        for symbol in targetGenes:
            if symbol in matrix_db:
                matrix_temp.append(matrix_db[symbol])
                header_temp.append(symbol)

        if len(header_temp) >= len(matrix_db):
            matrix2 = matrix_temp
            row_header2 = header_temp
    if transpose:
        matrix2 = map(numpy.array, zip(*matrix2))
        column_header, row_header2 = row_header2, column_header
    exclude = []
    exportData.write(string.join(['UID'] + column_header, '\t') + '\n')
    i = 0
    for row_id in row_header2:
        if ':' in row_id and '|' not in row_id:
            a, b = string.split(row_id, ':')[:2]
            if 'ENS' in a:
                try:
                    row_id = string.replace(row_id, a, gene_to_symbol[a][0])
                except Exception as e:
                    pass

                row_header2[i] = row_id
        elif 'ENS' in row_id and ' ' in row_id and '|' not in row_id:
            row_id = string.split(row_id, ' ')[1]
            row_header2[i] = row_id
        elif ' ' in row_id:
            try:
                a, b = string.split(row_id, ' ')
            except Exception:
                a = 1
                b = 2

            if a == b:
                row_id = a
        if row_id not in exclude:
            exportData.write(string.join([row_id] + map(str, matrix2[i]), '\t') + '\n')
        i += 1

    if 'amplify' not in targetGene and 'correlated' not in targetGene:
        print len(row_header2), 'input gene IDs found'
    else:
        print len(row_header2), 'top-correlated IDs'
    exportData.close()
    return (
     matrix2, row_header2, column_header, row_method)


def numpyCorrelationMatrixGeneStore(x, rows, genes, gene_to_symbol):
    start = time.time()
    output_file = string.replace(originalFilename, '.txt', '.corrmatrix')
    status = verifyFile(output_file)
    gene_correlations = {}
    if status == 'found':
        try:
            symbol = gene_to_symbol[rows[i]][0]
        except Exception:
            symbol = '$'

        def splitInt(x):
            rho, ind = string.split(x, '|')
            return (
             float(rho), int(float(ind)))

        for line in open(output_file, 'rU').xreadlines():
            data = line.rstrip()
            t = string.split(data, '\t')
            scores = map(lambda x: splitInt(x), t[1:])
            gene_correlations[t[0]] = scores

    else:
        eo = export.ExportFile(output_file)
        D1 = numpy.corrcoef(x)
        i = 0
        for score_ls in D1:
            scores = []
            try:
                symbol = gene_to_symbol[rows[i]][0]
            except Exception:
                symbol = '$'
            else:
                if rows[i] in genes or symbol in genes:
                    k = 0
                    for v in score_ls:
                        if str(v) != 'nan':
                            scores.append((v, k))
                        k += 1

                    scores.sort()
                    scores.reverse()
                    if len(symbol) == 1:
                        symbol = rows[i]
                    gene_correlations[symbol] = scores
                    export_values = [symbol]
                    for v, k in scores:
                        export_values.append(str(v)[:5] + '|' + str(k))

                    eo.write(string.join(export_values, '\t') + '\n')
                i += 1

        eo.close()
    print len(gene_correlations)
    print time.time() - start, 'seconds'
    sys.exit()
    return gene_correlations


def numpyCorrelationMatrixGene(x, rows, genes, gene_to_symbol):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        D1 = numpy.corrcoef(x)
    i = 0
    gene_correlations = {}
    for score_ls in D1:
        scores = []
        try:
            symbol = gene_to_symbol[rows[i]][0]
        except Exception:
            symbol = '$'
        else:
            if rows[i] in genes or symbol in genes:
                k = 0
                for v in score_ls:
                    if str(v) != 'nan':
                        scores.append((v, k))
                    k += 1

                scores.sort()
                scores.reverse()
                if len(symbol) == 1:
                    symbol = rows[i]
                gene_correlations[symbol] = scores
            i += 1

    return gene_correlations


def runHCOnly(filename, graphics, Normalize=False):
    """ Simple method for hieararchical clustering with defaults defined by the function rather than the user (see above function) """
    global EliteGeneSets
    global GroupDB
    global allowLargeClusters
    global graphic_link
    global inputFilename
    global justShowTheseIDs
    global root_dir
    global runGOElite
    runGOElite = False
    EliteGeneSets = []
    allowLargeClusters = False
    targetGene = []
    filterByPathways = False
    justShowTheseIDs = []
    graphic_link = graphics
    inputFilename = filename
    root_dir = findParentDir(filename)
    if 'ExpressionOutput/Clustering' in root_dir:
        root_dir = string.replace(root_dir, 'ExpressionOutput/Clustering', 'DataPlots')
    elif 'ExpressionOutput' in root_dir:
        root_dir = string.replace(root_dir, 'ExpressionOutput', 'DataPlots')
    else:
        root_dir += '/DataPlots/'
        try:
            os.mkdir(root_dir)
        except Exception:
            None

    row_method = 'average'
    column_method = 'weighted'
    row_metric = 'cosine'
    column_metric = 'cosine'
    if 'Lineage' in filename or 'Elite' in filename:
        color_gradient = 'red_white_blue'
    else:
        color_gradient = 'yellow_black_blue'
        color_gradient = 'red_black_sky'
    matrix, column_header, row_header, dataset_name, group_db = importData(filename, Normalize=Normalize)
    GroupDB = group_db
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, display=False, Normalize=Normalize)
    return graphic_link


def timestamp():
    import datetime
    today = str(datetime.date.today())
    today = string.split(today, '-')
    today = today[0] + '' + today[1] + '' + today[2]
    time_stamp = string.replace(time.ctime(), ':', '')
    time_stamp = string.replace(time_stamp, '  ', ' ')
    time_stamp = string.split(time_stamp, ' ')
    time_stamp = today + '-' + time_stamp[3]
    return time_stamp


def runPCAonly(filename, graphics, transpose, showLabels=True, plotType='3D', display=True, algorithm='SVD', geneSetName=None, species=None, zscore=True, colorByGene=None, reimportModelScores=True, separateGenePlots=False, forceClusters=False, maskGroups=None):
    global graphic_link
    global root_dir
    graphic_link = graphics
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir, '/DataPlots', '')
    root_dir = string.replace(root_dir, '/amplify', '')
    root_dir = string.replace(root_dir, 'ExpressionOutput/Clustering', 'DataPlots')
    root_dir = string.replace(root_dir, 'ExpressionInput', 'DataPlots')
    root_dir = string.replace(root_dir, 'ICGS-NMF', 'DataPlots')
    if 'DataPlots' not in root_dir:
        root_dir += '/DataPlots/'
    try:
        os.mkdir(root_dir)
    except Exception:
        None

    geneFilter = None
    if (algorithm == 't-SNE' or algorithm == 'UMAP') and reimportModelScores:
        dataset_name = string.split(filename, '/')[(-1)][:-4]
        try:
            if algorithm == 't-SNE':
                importtSNEScores(root_dir + dataset_name + '-t-SNE_scores.txt')
            if algorithm == 'UMAP':
                importtSNEScores(root_dir + dataset_name + '-UMAP_scores.txt')
            if len(colorByGene) == None:
                geneFilter = [
                 '']
            elif ' ' in colorByGene or ',' in colorByGene:
                colorByGene = string.replace(colorByGene, ',', ' ')
                geneFilter = string.split(colorByGene, ' ')
            else:
                geneFilter = [
                 colorByGene]
        except Exception:
            geneFilter = None

    matrix, column_header, row_header, dataset_name, group_db = importData(filename, zscore=zscore, geneFilter=geneFilter, forceClusters=forceClusters)
    if transpose == False:
        matrix = map(numpy.array, zip(*matrix))
        column_header, row_header = row_header, column_header
    if (len(column_header) > 1000 or len(row_header) > 1000) and algorithm != 't-SNE' and algorithm != 'UMAP':
        print 'Performing Principal Component Analysis (please be patient)...'
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        if algorithm == 't-SNE' or algorithm == 'UMAP':
            matrix = map(numpy.array, zip(*matrix))
            column_header, row_header = row_header, column_header
            if separateGenePlots and (len(colorByGene) > 0 or colorByGene == None):
                for gene in geneFilter:
                    tSNE(numpy.array(matrix), column_header, dataset_name, group_db, display=False, showLabels=showLabels, row_header=row_header, colorByGene=gene, species=species, reimportModelScores=reimportModelScores, method=algorithm)

                if display:
                    tSNE(numpy.array(matrix), column_header, dataset_name, group_db, display=True, showLabels=showLabels, row_header=row_header, colorByGene=gene, species=species, reimportModelScores=reimportModelScores, method=algorithm)
            elif maskGroups != None:
                import ExpressionBuilder
                sample_group_db = ExpressionBuilder.simplerGroupImport(maskGroups)
                group_sample_db = {}
                for sample in sample_group_db:
                    try:
                        group_sample_db[sample_group_db[sample]].append(sample)
                    except:
                        group_sample_db[sample_group_db[sample]] = [sample]

                for group in group_sample_db:
                    restricted_samples = group_sample_db[group]
                    tSNE(numpy.array(matrix), column_header, dataset_name, group_db, display=display, showLabels=showLabels, row_header=row_header, colorByGene=colorByGene, species=species, reimportModelScores=reimportModelScores, method=algorithm, maskGroups=(group, restricted_samples))

            else:
                tSNE(numpy.array(matrix), column_header, dataset_name, group_db, display=display, showLabels=showLabels, row_header=row_header, colorByGene=colorByGene, species=species, reimportModelScores=reimportModelScores, method=algorithm)
        elif plotType == '3D':
            try:
                PCA3D(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName, species=species, colorByGene=colorByGene)
            except Exception:
                print traceback.format_exc()
                PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName, species=species, colorByGene=colorByGene, reimportModelScores=reimportModelScores)

        else:
            PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db, display=display, showLabels=showLabels, algorithm=algorithm, geneSetName=geneSetName, species=species, colorByGene=colorByGene, reimportModelScores=reimportModelScores)
    return graphic_link


def outputClusters(filenames, graphics, Normalize=False, Species=None, platform=None, vendor=None):
    """ Peforms PCA and Hiearchical clustering on exported log-folds from AltAnalyze """
    global EliteGeneSets
    global GroupDB
    global allowLargeClusters
    global graphic_link
    global inputFilename
    global root_dir
    global runGOElite
    global species
    EliteGeneSets = []
    runGOElite = False
    allowLargeClusters = False
    graphic_link = graphics
    filename = filenames[0]
    inputFilename = filename
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir, 'ExpressionOutput/Clustering', 'DataPlots')
    original = importData(filename, Normalize=Normalize)
    matrix, column_header, row_header, dataset_name, group_db = original
    matrix = map(numpy.array, zip(*matrix))
    column_header, row_header = row_header, column_header
    if len(row_header) < 700000 and len(column_header) < 700000 and len(column_header) > 2:
        PrincipalComponentAnalysis(numpy.array(matrix), row_header, column_header, dataset_name, group_db)
    else:
        print 'SKIPPING PCA!!! - Your dataset file is over or under the recommended size limit for clustering (>7000 rows). Please cluster later using "Additional Analyses".'
    row_method = 'average'
    column_method = 'average'
    row_metric = 'cosine'
    column_metric = 'cosine'
    color_gradient = 'red_white_blue'
    color_gradient = 'red_black_sky'
    species = Species
    if 'LineageCorrelations' not in filename and 'Zscores' not in filename:
        EliteGeneSets = [
         'GeneOntology']
        runGOElite = True
    matrix, column_header, row_header, dataset_name, group_db = original
    GroupDB = group_db
    runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, Normalize=Normalize)
    for filename in filenames[1:]:
        inputFilename = filename
        matrix, column_header, row_header, dataset_name, group_db = importData(filename, Normalize=Normalize)
        GroupDB = group_db
        try:
            runHierarchicalClustering(matrix, row_header, column_header, dataset_name, row_method, row_metric, column_method, column_metric, color_gradient, Normalize=Normalize)
        except Exception:
            print 'Could not cluster', inputFilename, ', file not found'

    return graphic_link


def importEliteGeneAssociations(gene_filename):
    fn = filepath(gene_filename)
    x = 0
    fold_db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if data[0] == '#':
            x = 0
        elif x == 0:
            x = 1
        else:
            geneid = t[0]
            symbol = t[1]
            fold = 0
            try:
                if '|' in t[6]:
                    fold = float(string.split(t[6])[0])
            except Exception:
                None
            else:
                try:
                    fold = float(t[6])
                except Exception:
                    None
                else:
                    fold_db[symbol] = fold

    return fold_db


def importPathwayLevelFolds(filename):
    fn = filepath(filename)
    x = 0
    folds_db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if len(data) == 0:
            x = 0
        elif x == 0:
            z_score_indexes = []
            i = 0
            z_headers = []
            for header in t:
                if 'z_score.' in header:
                    z_score_indexes.append(i)
                    header = string.split(header, 'z_score.')[1]
                    if 'AS.' in header:
                        header = string.split(header, '.p')[0]
                        header = 'AS.' + string.join(string.split(header, '_')[2:], '_')
                    else:
                        header = string.join(string.split(header, '-')[:-2], '-')
                        if '-fold' in header:
                            header = string.join(string.split(header, '-')[:-1], '-')
                    z_headers.append(header)
                i += 1

            headers = string.join(['Gene-Set Name'] + z_headers, '\t') + '\n'
            x = 1
        else:
            term_name = t[1]
            geneset_type = t[2]
            zscores = map(lambda x: t[x], z_score_indexes)
            max_z = max(map(float, zscores))
            line = string.join([term_name] + zscores, '\t') + '\n'
            try:
                zscore_db[geneset_type].append((max_z, line))
            except Exception:
                zscore_db[geneset_type] = [(max_z, line)]

    exported_files = []
    for geneset_type in zscore_db:
        clusterinput_filename = findParentDir(filename) + '/Heatmaps/Clustering-Zscores-' + geneset_type + '.txt'
        exported_files.append(clusterinput_filename)
        export_text = export.ExportFile(clusterinput_filename)
        export_text.write(headers)
        zscore_db[geneset_type].sort()
        zscore_db[geneset_type].reverse()
        i = 0
        for max_z, line in zscore_db[geneset_type]:
            if i < 60:
                export_text.write(line)
            i += 1

        export_text.close()

    return exported_files


def importOverlappingEliteScores(filename):
    fn = filepath(filename)
    x = 0
    zscore_db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if len(data) == 0:
            x = 0
        elif x == 0:
            z_score_indexes = []
            i = 0
            z_headers = []
            for header in t:
                if 'z_score.' in header:
                    z_score_indexes.append(i)
                    header = string.split(header, 'z_score.')[1]
                    if 'AS.' in header:
                        header = string.split(header, '.p')[0]
                        header = 'AS.' + string.join(string.split(header, '_')[2:], '_')
                    else:
                        header = string.join(string.split(header, '-')[:-2], '-')
                        if '-fold' in header:
                            header = string.join(string.split(header, '-')[:-1], '-')
                    z_headers.append(header)
                i += 1

            headers = string.join(['Gene-Set Name'] + z_headers, '\t') + '\n'
            x = 1
        else:
            term_name = t[1]
            geneset_type = t[2]
            zscores = map(lambda x: t[x], z_score_indexes)
            max_z = max(map(float, zscores))
            line = string.join([term_name] + zscores, '\t') + '\n'
            try:
                zscore_db[geneset_type].append((max_z, line))
            except Exception:
                zscore_db[geneset_type] = [(max_z, line)]

    exported_files = []
    for geneset_type in zscore_db:
        clusterinput_filename = findParentDir(filename) + '/Heatmaps/Clustering-Zscores-' + geneset_type + '.txt'
        exported_files.append(clusterinput_filename)
        export_text = export.ExportFile(clusterinput_filename)
        export_text.write(headers)
        zscore_db[geneset_type].sort()
        zscore_db[geneset_type].reverse()
        i = 0
        for max_z, line in zscore_db[geneset_type]:
            if i < 60:
                export_text.write(line)
            i += 1

        export_text.close()

    return exported_files


def buildGraphFromSIF(mod, species, sif_filename, ora_input_dir):
    """ Imports a SIF and corresponding gene-association file to get fold changes for standardized gene-symbols """
    global SpeciesCode
    SpeciesCode = species
    mod = 'Ensembl'
    if sif_filename == None:
        sif_filename = '/Users/nsalomonis/Desktop/dataAnalysis/collaborations/WholeGenomeRVista/Alex-Figure/GO-Elite_results/CompleteResults/ORA_pruned/up-2f_p05-WGRV.sif'
        ora_input_dir = '/Users/nsalomonis/Desktop/dataAnalysis/collaborations/WholeGenomeRVista/Alex-Figure/up-stringent/up-2f_p05.txt'
    gene_filename = string.replace(sif_filename, '.sif', '_%s-gene-associations.txt') % mod
    gene_filename = string.replace(gene_filename, 'ORA_pruned', 'ORA_pruned/gene_associations')
    pathway_name = string.split(sif_filename, '/')[(-1)][:-4]
    output_filename = None
    try:
        fold_db = importEliteGeneAssociations(gene_filename)
    except Exception:
        fold_db = {}
    else:
        if ora_input_dir != None:
            try:
                fold_db = importDataSimple(ora_input_dir, species, fold_db, mod)
            except Exception:
                None

        try:
            output_filename = iGraphSimple(sif_filename, fold_db, pathway_name)
        except Exception:
            print 'igraph export failed due to an unknown error (not installed)'
            print traceback.format_exc()
            try:
                displaySimpleNetwork(sif_filename, fold_db, pathway_name)
            except Exception:
                pass

    return output_filename


def iGraphSimple(sif_filename, fold_db, pathway_name):
    """ Build a network export using iGraph and Cairo """
    edges = importSIF(sif_filename)
    id_color_db = WikiPathways_webservice.getHexadecimalColorRanges(fold_db, 'Genes')
    output_filename = iGraphDraw(edges, pathway_name, filePath=sif_filename, display=True, graph_layout='spring', colorDB=id_color_db)
    return output_filename


def iGraphDraw(edges, pathway_name, labels=None, graph_layout='shell', display=False, node_size=700, node_color='yellow', node_alpha=0.5, node_text_size=7, edge_color='black', edge_alpha=0.5, edge_thickness=2, edges_pos=0.3, text_font='sans-serif', filePath='test', colorDB=None):
    output_filename = None
    if len(edges) > 700 and 'AltAnalyze' not in pathway_name:
        print findFilename(filePath), 'too large to visualize...'
    elif len(edges) > 3000:
        print findFilename(filePath), 'too large to visualize...'
    else:
        arrow_scaler = 1
        if edges > 40:
            arrow_scaler = 0.9
        vars = formatiGraphEdges(edges, pathway_name, colorDB, arrow_scaler)
        vertices, iGraph_edges, vertice_db, label_list, shape_list, vertex_size, color_list, vertex_label_colors, arrow_width, edge_colors = vars
        if vertices > 0:
            import igraph
            gr = igraph.Graph(vertices, directed=True)
            canvas_scaler = 0.8
            if vertices < 15:
                canvas_scaler = 0.5
            elif vertices < 25:
                canvas_scaler = 0.7
            elif vertices > 35:
                canvas_scaler += len(iGraph_edges) / 400.0
            filePath, canvas_scaler = correctedFilePath(filePath, canvas_scaler)
            canvas_size = (
             600 * canvas_scaler, 600 * canvas_scaler)
            gr.add_edges(iGraph_edges)
            gr.vs['label'] = label_list
            gr.vs['shape'] = shape_list
            gr.vs['size'] = vertex_size
            gr.vs['label_dist'] = [1.3] * vertices
            gr.vs['label_size'] = [12] * vertices
            gr.vs['color'] = color_list
            gr.vs['label_color'] = vertex_label_colors
            gr.es['color'] = edge_colors
            gr.es['arrow_size'] = arrow_width
            output_filename = '%s.pdf' % filePath[:-4]
            output_filename = output_filename.encode('ascii', 'ignore')
            layout = 'kk'
            visual_style = {}
            visual_style['margin'] = 50
            visual_style['bbox'] = canvas_size
            igraph.plot(gr, output_filename, **visual_style)
            output_filename = '%s.png' % filePath[:-4]
            output_filename = output_filename.encode('ascii', 'ignore')
            if vertices < 15:
                gr, visual_style = increasePlotSize(gr, visual_style)
            igraph.plot(gr, output_filename, **visual_style)
    return output_filename


def correctedFilePath(filePath, canvas_scaler):
    """ Move this file to it's own network directory for GO-Elite """
    if 'ORA_pruned' in filePath:
        filePath = string.replace(filePath, 'CompleteResults/ORA_pruned', 'networks')
        try:
            os.mkdir(findParentDir(filePath))
        except Exception:
            pass

        canvas_scaler = canvas_scaler * 1.3
    return (filePath, canvas_scaler)


def increasePlotSize(gr, visual_style):
    factor = 2
    object_list = ['size', 'label_size']
    for i in object_list:
        new = []
        for k in gr.vs[i]:
            new.append(k * factor)

        gr.vs[i] = new

    new = []
    for i in gr.es['arrow_size']:
        new.append(i * factor)

    new = []
    for i in visual_style['bbox']:
        new.append(i * factor)

    visual_style['bbox'] = new
    visual_style['margin'] = visual_style['margin'] * factor
    return (
     gr, visual_style)


def getHMDBDataSimple():
    program_type, database_dir = unique.whatProgramIsThis()
    filename = database_dir + '/' + SpeciesCode + '/gene/HMDB.txt'
    symbol_hmdb_db = {}
    x = 0
    fn = filepath(filename)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0:
            x = 1
        else:
            t = string.split(data, '\t')
            hmdb_id = t[0]
            symbol = t[1]
            ProteinNames = t[(-1)]
            symbol_hmdb_db[symbol] = hmdb_id

    return symbol_hmdb_db


def formatiGraphEdges(edges, pathway_name, colorDB, arrow_scaler):
    edge_db = {}
    edges2 = []
    vertice_db = {}
    shape_list = []
    label_list = []
    vertex_size = []
    color_list = []
    vertex_label_colors = []
    arrow_width = []
    edge_colors = []
    k = 0
    try:
        symbol_hmdb_db = getHMDBDataSimple()
    except Exception:
        symbol_hmdb_db = {}

    for node1, node2, type in edges:
        edge_color = 'grey'
        if 'TF' in pathway_name or 'WGRV' in pathway_name:
            pathway = node1
        else:
            pathway = node2
        if 'drugInteraction' == type:
            edge_color = 'purple'
        else:
            if 'TBar' == type:
                edge_color = 'blue'
            else:
                if 'microRNAInteraction' == type:
                    edge_color = '#53A26D'
                elif 'transcription' in type:
                    edge_color = '#FF7D7D'
                if 'AltAnalyze' in pathway_name:
                    default_node_color = 'grey'
                else:
                    default_node_color = 'yellow'
                if node1 in vertice_db:
                    v1 = vertice_db[node1]
                else:
                    v1 = k
                    label_list.append(node1)
                    rs = 1
                    if 'TF' in pathway_name or 'WGRV' in pathway_name and 'AltAnalyze' not in pathway_name:
                        shape_list.append('rectangle')
                        vertex_size.append(15)
                        vertex_label_colors.append('blue')
                    else:
                        if 'drugInteraction' == type:
                            rs = 0.75
                            shape_list.append('rectangle')
                            vertex_label_colors.append('purple')
                            default_node_color = 'purple'
                        elif 'Metabolic' == type and node1 in symbol_hmdb_db:
                            shape_list.append('triangle-up')
                            vertex_label_colors.append('blue')
                            default_node_color = 'grey'
                        elif 'microRNAInteraction' == type:
                            rs = 0.75
                            shape_list.append('triangle-up')
                            vertex_label_colors.append('#008000')
                            default_node_color = 'grey'
                        else:
                            shape_list.append('circle')
                            vertex_label_colors.append('black')
                        vertex_size.append(10 * rs)
                    vertice_db[node1] = v1
                    k += 1
                    try:
                        color = '#' + string.upper(colorDB[node1])
                        color_list.append(color)
                    except Exception:
                        color_list.append(default_node_color)

            if node2 in vertice_db:
                v2 = vertice_db[node2]
            else:
                v2 = k
                label_list.append(node2)
                if 'TF' in pathway_name or 'WGRV' in pathway_name:
                    shape_list.append('circle')
                    vertex_size.append(10)
                    vertex_label_colors.append('black')
                    default_node_color = 'grey'
                elif 'AltAnalyze' not in pathway_name:
                    shape_list.append('rectangle')
                    vertex_size.append(15)
                    vertex_label_colors.append('blue')
                    default_node_color = 'grey'
                elif 'Metabolic' == type and node2 in symbol_hmdb_db:
                    shape_list.append('triangle-up')
                    vertex_label_colors.append('blue')
                    default_node_color = 'grey'
                else:
                    shape_list.append('circle')
                    vertex_size.append(10)
                    vertex_label_colors.append('black')
                    default_node_color = 'grey'
                vertice_db[node2] = v2
                k += 1
                try:
                    color = '#' + string.upper(colorDB[node2])
                    color_list.append(color)
                except Exception:
                    color_list.append(default_node_color)

        edges2.append((v1, v2))
        if type == 'physical':
            arrow_width.append(0)
        else:
            arrow_width.append(arrow_scaler)
        try:
            edge_db[v1].append(v2)
        except Exception:
            edge_db[v1] = [v2]
        else:
            try:
                edge_db[v2].append(v1)
            except Exception:
                edge_db[v2] = [v1]
            else:
                edge_colors.append(edge_color)

    vertices = len(edge_db)
    edge_db = eliminate_redundant_dict_values(edge_db)
    vertice_db2 = {}
    for node in vertice_db:
        vertice_db2[vertice_db[node]] = node

    print vertices, 'and', len(edges2), 'edges in the iGraph network.'
    return (
     vertices, edges2, vertice_db2, label_list, shape_list, vertex_size, color_list, vertex_label_colors, arrow_width, edge_colors)


def eliminate_redundant_dict_values(database):
    db1 = {}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list

    return db1


def importDataSimple(filename, species, fold_db, mod):
    """ Imports an input ID file and converts those IDs to gene symbols for analysis with folds """
    import GO_Elite
    from import_scripts import OBO_import
    import gene_associations
    fn = filepath(filename)
    x = 0
    metabolite_codes = ['Ck', 'Ca', 'Ce', 'Ch', 'Cp']
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if data[0] == '#':
            x = 0
        elif x == 0:
            si = None
            symbol_present = False
            try:
                si = t.index('Symbol')
                symbol_present = True
            except:
                pass

            x = 1
        else:
            if x == 1:
                system_code = t[1]
                if system_code in metabolite_codes:
                    mod = 'HMDB'
                system_codes, source_types, mod_types = GO_Elite.getSourceData()
                try:
                    source_data = system_codes[system_code]
                except Exception:
                    source_data = None
                    if 'ENS' in t[0]:
                        source_data = system_codes['En']
                    else:
                        source_data = system_codes['Sy']
                else:
                    if source_data == mod:
                        source_is_mod = True
                    elif source_data == None:
                        None
                    else:
                        source_is_mod = False
                        mod_source = mod + '-' + source_data + '.txt'
                        gene_to_source_id = gene_associations.getGeneToUid(species, ('hide', mod_source))
                        source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
                    try:
                        gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', mod + '-Symbol'))
                    except Exception:
                        gene_to_symbol = {}

                    try:
                        met_to_symbol = gene_associations.importGeneData(species, 'HMDB', simpleImport=True)
                    except Exception:
                        met_to_symbol = {}

                for i in met_to_symbol:
                    gene_to_symbol[i] = met_to_symbol[i]

            x += 1
            if source_is_mod == True:
                if t[0] in gene_to_symbol:
                    symbol = gene_to_symbol[t[0]][0]
                    try:
                        fold_db[symbol] = float(t[2])
                    except Exception:
                        fold_db[symbol] = 0

                else:
                    fold_db[t[0]] = 0
            elif symbol_present:
                fold_db[t[si]] = 0
                try:
                    fold_db[t[si]] = float(t[2])
                except Exception:
                    try:
                        fold_db[t[si]] = 0
                    except:
                        fold_db[t[0]] = 0

            elif t[0] in source_to_gene:
                mod_ids = source_to_gene[t[0]]
                try:
                    mod_ids += source_to_gene[t[2]]
                except Exception:
                    try:
                        mod_ids += source_to_gene[t[1]]
                    except Exception:
                        None

                for mod_id in mod_ids:
                    if mod_id in gene_to_symbol:
                        symbol = gene_to_symbol[mod_id][0]
                        try:
                            fold_db[symbol] = float(t[2])
                        except Exception:
                            fold_db[symbol] = 0

            else:
                fold_db[t[0]] = 0

    return fold_db


def clusterPathwayZscores(filename):
    """ Imports a overlapping-results file and exports an input file for hierarchical clustering and clusters """
    if filename == None:
        filename = '/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/GO-Elite/TF-enrichment2/GO-Elite_results/overlapping-results_z-score_elite.txt'
    exported_files = importOverlappingEliteScores(filename)
    graphic_links = []
    for file in exported_files:
        try:
            graphic_links = runHCOnly(file, graphic_links)
        except Exception as e:
            print 'Unable to generate cluster due to dataset incompatibilty.'

    print 'Clustering of overlapping-results_z-score complete (see "GO-Elite_results/Heatmaps" directory)'
    return


def clusterPathwayMeanFolds():
    """ Imports the pruned-results file and exports an input file for hierarchical clustering and clusters """
    filename = '/Users/nsalomonis/Desktop/User Diagnostics/Mm_spinal_cord_injury/GO-Elite/GO-Elite_results/pruned-results_z-score_elite.txt'
    exported_files = importPathwayLevelFolds(filename)


def VennDiagram():
    f = pylab.figure()
    ax = f.gca()
    rad = 1.4
    c1 = Circle((-1, 0), rad, alpha=0.2, fc='red', label='red')
    c2 = Circle((1, 0), rad, alpha=0.2, fc='blue', label='blue')
    c3 = Circle((0, 1), rad, alpha=0.2, fc='green', label='g')
    ax.add_patch(c2)
    ax.add_patch(c3)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    pylab.show()


def plotHistogram(filename):
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    transpose = True
    if transpose:
        print 'Transposing the data matrix'
        matrix = map(numpy.array, zip(*matrix))
        column_header, row_header = row_header, column_header
    pylab.figure()
    for i in matrix:
        pylab.hist(i, 200, normed=0, histtype='step', cumulative=-1)

    pylab.show()


def stackedbarchart(filename, display=False, output=False):
    header = []
    conditions = []
    data_matrix = []
    for line in open(filename, 'rU').xreadlines():
        cd = cleanUpLine(line)
        t = string.split(cd, '\t')
        if len(header) == 0:
            header = t[4:]
            exc_indexes = [0, 2, 4, 6, 8, 10, 12]
            inc_indexes = [1, 3, 5, 7, 9, 11, 13]
            inlc_header = map(lambda i: string.split(header[i], '_')[0], inc_indexes)
            header = inlc_header
        else:
            condition = t[0]
            data = t[4:]
            conditions.append(condition + '-inclusion ')
            data_matrix.append(map(lambda i: float(data[i]), inc_indexes))
            conditions.append(condition + '-exclusion ')
            data_matrix.append(map(lambda i: float(data[i]), exc_indexes))

    data_matrix = map(numpy.array, zip(*data_matrix))
    y_pos = np.arange(len(conditions))
    fig, ax = pylab.subplots()
    colors = [
     'royalblue', 'salmon', 'grey', 'gold', 'cornflowerblue', 'mediumseagreen', 'navy']
    patch_handles = []
    left = np.zeros(len(conditions))
    index = 0
    for i, d in enumerate(data_matrix):
        patch_handles.append(ax.barh(y_pos, d, 0.3, color=colors[index], align='center', left=left, label=header[index]))
        left += d
        index += 1

    ax.set_yticks(y_pos)
    ax.set_yticklabels(conditions)
    ax.set_xlabel('Events')
    ax.legend(loc='best', bbox_to_anchor=(1.0, 1.0))
    box = ax.get_position()
    ax.set_position([box.x0 + 0.2, box.y0, box.width * 0.6, box.height])
    try:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    except Exception:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_title('MultiPath-PSI Splicing Event Types')
    if output == False:
        pylab.savefig(filename[:-4] + '.pdf')
        pylab.savefig(filename[:-4] + '.png')
    else:
        pylab.savefig(output[:-4] + '.pdf')
        pylab.savefig(output[:-4] + '.png')
    if display:
        print 'Exporting:', filename
        try:
            pylab.show()
        except Exception:
            None

    return


def barchart(filename, index1, index2, x_axis, y_axis, title, display=False, color1='gold', color2='darkviolet', output=False):
    header = []
    reference_data = []
    query_data = []
    groups = []
    for line in open(filename, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if len(header) == 0:
            header = t
            header1 = header[index1]
            header2 = header[index2]
        else:
            reference_data.append(float(t[index1]))
            q_value = float(t[index2])
            if 'frequen' not in filename:
                q_value = q_value * -1
            query_data.append(q_value)
            name = t[0]
            if '_vs_' in name and 'event_summary' not in filename:
                name = string.split(name, '_vs_')[0]
                suffix = None
                if '__' in name:
                    suffix = string.split(name, '__')[(-1)]
                if '_' in name:
                    name = string.split(name, '_')[:-1]
                    name = string.join(name, '_')
                    if len(name) > 20:
                        name = string.split(name, '_')[0]
                if suffix != None:
                    name += '_' + suffix
            groups.append(name)

    fig, ax = pylab.subplots()
    pos1 = ax.get_position()
    pos2 = [pos1.x0 + 0.2, pos1.y0 + 0.1, pos1.width / 1.2, pos1.height / 1.2]
    ax.set_position(pos2)
    ind = np.arange(len(groups))
    width = 0.35
    query_data.reverse()
    reference_data.reverse()
    groups.reverse()
    ax.barh(ind - width / 2, query_data, width, color=color2, label=header2)
    ax.barh(ind + width / 2, reference_data, width, color=color1, label=header1)
    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    ax.set_yticks(ind + 0.175)
    ax.set_yticklabels(groups)
    ax.set_title(title)
    ax.legend()
    if output == False:
        pylab.savefig(filename[:-4] + '.pdf')
    else:
        pylab.savefig(output[:-4] + '.pdf')
    if display:
        print 'Exporting:', filename
        try:
            pylab.show()
        except Exception:
            None

    return


def multipleSubPlots(filename, uids, SubPlotType='column', n=20):
    str_uids = string.join(uids, '_')
    matrix, column_header, row_header, dataset_name, group_db = importData(filename, geneFilter=uids)
    for uid in uids:
        if uid not in row_header:
            print uid, 'is missing from the expression file.'

    fig = pylab.figure()

    def ReplaceZeros(val, min_val):
        if val == 0:
            return min_val
        else:
            return val

    new_row_header = []
    matrix2 = []
    for uid in uids:
        if uid in row_header:
            ind = row_header.index(uid)
            new_row_header.append(uid)
            try:
                update_exp_vals = map(lambda x: ReplaceZeros(x, 0.0001), matrix[ind])
            except Exception:
                print uid, len(matrix[ind])
                sys.exit()

            matrix2.append(update_exp_vals)

    matrix = numpy.array(matrix2)
    row_header = new_row_header
    color_list = [
     'r', 'b', 'y', 'g', 'w', 'k', 'm']
    groups = []
    for sample in column_header:
        try:
            group = group_db[sample][0]
        except:
            group = '1'
        else:
            if group not in groups:
                groups.append(group)

    fontsize = 10
    if len(groups) > 0:
        color_list = []
        if len(groups) == 9:
            cm = matplotlib.colors.ListedColormap(['#80C241', '#118943', '#6FC8BB', '#ED1D30', '#F26E21', '#8051A0', '#4684C5', '#FBD019', '#3A52A4'])
        elif len(groups) == 3:
            cm = matplotlib.colors.ListedColormap(['#4684C4', '#FAD01C', '#7D7D7F'])
        else:
            cm = pylab.cm.get_cmap('gist_rainbow')
        for i in range(len(groups)):
            color_list.append(cm(1.0 * i / len(groups)))

    for i in range(len(matrix)):
        ax = pylab.subplot(n, 1, 1 + i)
        OY = matrix[i]
        pylab.xlim(0, len(OY))
        pylab.subplots_adjust(right=0.85)
        ind = np.arange(len(OY))
        index_list = []
        v_list = []
        colors_list = []
        if SubPlotType == 'column':
            index = -1
            for v in OY:
                index += 1
                try:
                    group = group_db[column_header[index]][0]
                except:
                    group = '1'
                else:
                    index_list.append(index)
                    v_list.append(v)
                    colors_list.append(color_list[groups.index(group)])
                    width = 0.35

            print 1
            barlist = pylab.bar(index_list, v_list, edgecolor='black', linewidth=0)
            ci = 0
            for cs in barlist:
                barlist[ci].set_color(colors_list[ci])
                ci += 1

        if SubPlotType == 'plot':
            pylab.plot(x, y)
        ax.text(matrix.shape[1] - 0.5, i, '  ' + row_header[i], fontsize=8)
        fig.autofmt_xdate()
        pylab.subplots_adjust(hspace=0.001)
        temp = tic.MaxNLocator(3)
        ax.yaxis.set_major_locator(temp)
        ax.set_xticks([])

    if len(str_uids) > 50:
        str_uids = str_uids[:50]
    pylab.savefig(filename[:-4] + '-1' + str_uids + '.pdf')


def simpleTranspose(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, ' ')
        matrix.append(t)

    matrix = map(numpy.array, zip(*matrix))
    filename = filename[:-4] + '-transposed.txt'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i, '\t') + '\n')

    ea.close()


def CorrdinateToBed(filename):
    fn = filepath(filename)
    matrix = []
    translation = {}
    multiExon = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data, ' ', '')
        t = string.split(data, '\t')
        if '.gtf' in filename:
            if 'chr' not in t[0]:
                chr = 'chr' + t[0]
            else:
                chr = t[0]
            start = t[3]
            end = t[4]
            strand = t[6]
            annotation = t[8]
            annotation = string.replace(annotation, 'gene_id', '')
            annotation = string.replace(annotation, 'transcript_id', '')
            annotation = string.replace(annotation, 'gene_name', '')
            geneIDs = string.split(annotation, ';')
            geneID = geneIDs[0]
            symbol = geneIDs[3]
        else:
            chr = t[4]
            strand = t[5]
            start = t[6]
            end = t[7]
        t = [
         chr, start, end, geneID, '0', strand]
        translation[geneID] = symbol
        try:
            multiExon[geneID] += 1
        except Exception:
            multiExon[geneID] = 1

    filename = filename[:-4] + '-new.bed'
    ea = export.ExportFile(filename)
    for i in translation:
        ea.write(i + '\t' + translation[i] + '\t' + str(multiExon[i]) + '\n')

    ea.close()


def SimpleCorrdinateToBed(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data, ' ', '')
        t = string.split(data, '\t')
        if '.bed' in filename:
            print t
            sys.exit()
        chr = t[4]
        strand = t[5]
        start = t[6]
        end = t[7]
        if 'ENS' in t[0]:
            t = [
             chr, start, end, t[0], '0', strand]
            matrix.append(t)

    filename = filename[:-4] + '-new.bed'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i, '\t') + '\n')

    ea.close()


def simpleIntegrityCheck(filename):
    fn = filepath(filename)
    matrix = []
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        data = string.replace(data, ' ', '')
        t = string.split(data, '\t')
        matrix.append(t)

    filename = filename[:-4] + '-new.bed'
    ea = export.ExportFile(filename)
    for i in matrix:
        ea.write(string.join(i, '\t') + '\n')

    ea.close()


def BedFileCheck(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.bed'
    ea = export.ExportFile(filename)
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            firstRow = False
        else:
            ea.write(string.join(t, '\t') + '\n')

    ea.close()


def simpleFilter(filename):
    fn = filepath(filename)
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    matrix = []
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, ',')
        uid = t[0]
        if 1 == 2:
            a, b = string.split(t[0], '=')
            b = string.replace(b, '_', ':')
            uid = a + '=' + b
            matrix.append(t)
        ea.write(string.join([uid] + t[1:], '\t') + '\n')

    ea.close()


def test(filename):
    symbols2 = {}
    firstLine = True
    fn = filepath(filename)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstLine:
            firstLine = False
            header = t
            i = 0
            start = None
            alt_start = None
            value_indexes = []
            groups = {}
            group = 0
            for h in header:
                if h == 'WikiPathways':
                    start = i
                if h == 'Select Protein Classes':
                    alt_start = i
                i += 1

            if start == None:
                start = alt_start
            for h in header:
                if h > i:
                    group[i]
                i += 1

            if start == None:
                start = alt_start
        else:
            uniprot = t[0]
            symbols = string.replace(t[(-1)], ';;', ';')
            symbols = string.split(symbols, ';')
            for s in symbols:
                if len(s) > 0:
                    symbols2[(string.upper(s), uniprot)] = []

    for s, u in symbols2:
        ea.write(string.join([s, u], '\t') + '\n')

    ea.close()
    return


def coincentIncedenceTest(exp_file, TFs):
    fn = filepath(TFs)
    tfs = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        tfs[data] = []

    comparisons = {}
    for tf1 in tfs:
        for tf2 in tfs:
            if tf1 != tf2:
                temp = [
                 tf1, tf2]
                temp.sort()
                comparisons[tuple(temp)] = []

    gene_data = {}
    firstLine = True
    fn = filepath(exp_file)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine:
            firstLine = False
            header = string.split(data, '\t')[1:]
        else:
            t = string.split(data, '\t')
            gene = t[0]
            values = map(float, t[1:])
            gene_data[gene] = values

    filename = TFs[:-4] + '-all-coincident-5z.txt'
    ea = export.ExportFile(filename)
    comparison_db = {}
    for comparison in comparisons:
        vals1 = gene_data[comparison[0]]
        vals2 = gene_data[comparison[1]]
        i = 0
        coincident = []
        for v1 in vals1:
            v2 = vals2[i]
            if v1 > 1 and v2 > 1:
                coincident.append(i)
            i += 1

        i = 0
        population_db = {}
        coincident_db = {}
        for h in header:
            population = string.split(h, ':')[0]
            if i in coincident:
                try:
                    coincident_db[population] += 1
                except Exception:
                    coincident_db[population] = 1

            try:
                population_db[population] += 1
            except Exception:
                population_db[population] = 1
            else:
                i += 1

        import mappfinder
        final_population_percent = []
        for population in population_db:
            d = population_db[population]
            try:
                c = coincident_db[population]
            except Exception:
                c = 0
            else:
                N = float(len(header))
                R = float(len(coincident))
                n = float(d)
                r = float(c)
                try:
                    z = mappfinder.Zscore(r, n, N, R)
                except Exception:
                    z = 0
                else:
                    final_population_percent.append([population, str(c), str(d), str(float(c) / float(d)), str(z)])

        comparison_db[comparison] = final_population_percent

    filtered_comparison_db = {}
    top_scoring_population = {}
    for comparison in comparison_db:
        max_group = []
        for population_stat in comparison_db[comparison]:
            z = float(population_stat[(-1)])
            c = float(population_stat[1])
            population = population_stat[0]
            max_group.append([z, population])

        max_group.sort()
        z = max_group[(-1)][0]
        pop = max_group[(-1)][1]
        if z > 3.92 and c > 3:
            filtered_comparison_db[comparison] = comparison_db[comparison]
            top_scoring_population[comparison] = (pop, z)

    firstLine = True
    for comparison in filtered_comparison_db:
        comparison_alt = string.join(list(comparison), '|')
        all_percents = []
        for line in filtered_comparison_db[comparison]:
            all_percents.append(line[3])

        if firstLine:
            all_headers = []
            for line in filtered_comparison_db[comparison]:
                all_headers.append(line[0])

            ea.write(string.join(['gene-pair'] + all_headers + ['Top Population', 'Top Z'], '\t') + '\n')
            firstLine = False
        pop, z = top_scoring_population[comparison]
        ea.write(string.join([comparison_alt] + all_percents + [pop, str(z)], '\t') + '\n')

    ea.close()


def getlastexon(filename):
    filename2 = filename[:-4] + '-last-exon.txt'
    ea = export.ExportFile(filename2)
    firstLine = True
    fn = filepath(filename)
    last_gene = 'null'
    last_exon = ''
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstLine:
            firstLine = False
        else:
            gene = t[2]
            if gene != last_gene:
                if ':E' in last_exon:
                    gene, exon = last_exon = string.split(':E')
                    block, region = string.split(exon, '.')
                    try:
                        ea.write(last_exon + '\n')
                    except:
                        pass

            last_gene = gene
            last_exon = t[0]

    ea.close()


def replaceWithBinary(filename):
    filename2 = filename[:-4] + '-binary.txt'
    ea = export.ExportFile(filename2)
    firstLine = True
    fn = filepath(filename)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstLine:
            ea.write(line)
            firstLine = False
        else:
            try:
                values = map(float, t[1:])
            except Exception:
                print t[1:]
                sys.exit()
            else:
                values2 = []
                for v in values:
                    if v == 0:
                        values2.append('0')
                    else:
                        values2.append('1')

                ea.write(string.join([t[0]] + values2, '\t') + '\n')

    ea.close()


def geneMethylationOutput(filename):
    filename2 = filename[:-4] + '-binary.txt'
    ea = export.ExportFile(filename2)
    firstLine = True
    fn = filepath(filename)
    db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        values = (t[20], t[3] + '-methylation')
        db[values] = []

    for value in db:
        ea.write(string.join(list(value), '\t') + '\n')

    ea.close()


def coincidentIncedence(filename, genes):
    exportPairs = True
    gene_data = []
    firstLine = True
    fn = filepath(filename)
    if exportPairs:
        filename = filename[:-4] + '_' + genes[0] + '-' + genes[1] + '2.txt'
        ea = export.ExportFile(filename)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine:
            firstLine = False
            header = string.split(data, '\t')[1:]
        else:
            t = string.split(data, '\t')
            gene = t[0]
            if gene in genes:
                values = map(float, t[1:])
                gene_data.append(values)

    vals1 = gene_data[0]
    vals2 = gene_data[1]
    i = 0
    coincident = []
    for v1 in vals1:
        v2 = vals2[i]
        if v1 > 1 and v2 > 1:
            coincident.append(i)
        i += 1

    i = 0
    population_db = {}
    coincident_db = {}
    for h in header:
        population = string.split(h, ':')[0]
        if i in coincident:
            try:
                coincident_db[population] += 1
            except Exception:
                coincident_db[population] = 1

        try:
            population_db[population] += 1
        except Exception:
            population_db[population] = 1
        else:
            i += 1

    import mappfinder
    final_population_percent = []
    for population in population_db:
        d = population_db[population]
        try:
            c = coincident_db[population]
        except Exception:
            c = 0
        else:
            N = float(len(header))
            R = float(len(coincident))
            n = d
            r = c
            try:
                z = mappfinder.zscore(r, n, N, R)
            except Exception:
                z = 0
            else:
                final_population_percent.append([population, str(c), str(d), str(float(c) / float(d)), str(z)])

    if exportPairs:
        for line in final_population_percent:
            ea.write(string.join(line, '\t') + '\n')

        ea.close()
    else:
        return final_population_percent


def extractFeatures(countinp, IGH_gene_file):
    import export
    ExonsPresent = False
    igh_genes = []
    firstLine = True
    for line in open(IGH_gene_file, 'rU').xreadlines():
        if firstLine:
            firstLine = False
        else:
            data = cleanUpLine(line)
            gene = string.split(data, '\t')[0]
            igh_genes.append(gene)

    if 'counts.' in countinp:
        feature_file = string.replace(countinp, 'counts.', 'IGH.')
        fe = export.ExportFile(feature_file)
        firstLine = True
        for line in open(countinp, 'rU').xreadlines():
            if firstLine:
                fe.write(line)
                firstLine = False
            else:
                feature_info = string.split(line, '\t')[0]
                gene = string.split(feature_info, ':')[0]
                if gene in igh_genes:
                    fe.write(line)

        fe.close()


def filterForJunctions(countinp):
    import export
    ExonsPresent = False
    igh_genes = []
    firstLine = True
    count = 0
    if 'counts.' in countinp:
        feature_file = countinp[:-4] + '-output.txt'
        fe = export.ExportFile(feature_file)
        firstLine = True
        for line in open(countinp, 'rU').xreadlines():
            if firstLine:
                fe.write(line)
                firstLine = False
            else:
                feature_info = string.split(line, '\t')[0]
                junction = string.split(feature_info, '=')[0]
                if '-' in junction:
                    fe.write(line)
                    count += 1

        fe.close()
    print count


def countIntronsExons(filename):
    import export
    exon_db = {}
    intron_db = {}
    firstLine = True
    last_transcript = None
    for line in open(filename, 'rU').xreadlines():
        if firstLine:
            firstLine = False
        else:
            line = line.rstrip()
            t = string.split(line, '\t')
            transcript = t[(-1)]
            chr = t[1]
            strand = t[2]
            start = t[3]
            end = t[4]
            exon_db[(chr, start, end)] = []
            if transcript == last_transcript:
                if strand == '1':
                    intron_db[(chr, last_end, start)] = []
                else:
                    intron_db[(chr, last_start, end)] = []
            last_end = end
            last_start = start
            last_transcript = transcript

    print len(exon_db) + 1, len(intron_db) + 1
    return


def importGeneList(gene_list_file, n=20):
    genesets = []
    genes = []
    for line in open(gene_list_file, 'rU').xreadlines():
        gene = line.rstrip()
        gene = string.split(gene, '\t')[0]
        genes.append(gene)
        if len(genes) == n:
            genesets.append(genes)
            genes = []

    if len(genes) > 0 and len(genes) < n + 1:
        genes += (n - len(genes)) * [gene]
        genesets.append(genes)
    return genesets


def simpleListImport(filename):
    genesets = []
    genes = []
    for line in open(filename, 'rU').xreadlines():
        gene = line.rstrip()
        gene = string.split(gene, '\t')[0]
        genes.append(gene)

    return genes


def customClean(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            firstRow = False
            ea.write(string.join(['UID'] + t, '\t') + '\n')
        else:
            if ';' in t[0]:
                uid = string.split(t[0], ';')[0]
            else:
                uid = t[0]
            values = map(lambda x: float(x), t[1:])
            values.sort()
            if values[3] >= 1:
                ea.write(string.join([uid] + t[1:], '\t') + '\n')

    ea.close()


def MakeJunctionFasta(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '.fasta'
    ea = export.ExportFile(filename)
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        probeset, seq = string.split(data, '\t')[:2]
        ea.write('>' + probeset + '\n')
        ea.write(string.upper(seq) + '\n')

    ea.close()


def ToppGeneFilter(filename):
    import gene_associations
    from import_scripts import OBO_import
    gene_to_symbol = gene_associations.getGeneToUid('Hs', ('hide', 'Ensembl-Symbol'))
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            firstRow = False
            ea.write(string.join(['Ensembl\t\tCategory'], '\t') + '\n')
        else:
            symbol = t[1]
            category = t[3]
            symbol = symbol[0] + string.lower(symbol[1:])
            category = category[:100]
            if symbol in symbol_to_gene:
                ensembl = symbol_to_gene[symbol][0]
                ea.write(string.join([ensembl, symbol, category], '\t') + '\n')

    ea.close()


def CountKallistoAlignedJunctions(filename):
    fn = filepath(filename)
    firstRow = True
    ea = export.ExportFile(filename)
    found = False
    counts = 0
    unique = {}
    ea = export.ExportFile(filename[:-4] + '-Mpo.txt')
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if 'ENS' in line and 'JUNC1201' in line:
            ea.write(line)
            unique[t[0]] = []
            counts += 1

    print counts, len(unique)
    ea.close()


def filterRandomFile(filename, col1, col2):
    fn = filepath(filename)
    firstRow = True
    counts = 0
    ea = export.ExportFile(filename[:-4] + '-columns.txt')
    for line in open(fn, 'rU').xreadlines():
        if line[0] != '#':
            data = line.rstrip()
            t = string.split(data, ',')
            if ' ' in t[(col2 - 1)]:
                t[col2 - 1] = string.split(t[(col2 - 1)], ' ')[2]
            ea.write(t[(col1 - 1)] + '\t' + t[(col2 - 1)] + '\n')
            counts += 1

    ea.close()


def getBlockExonPositions():
    fn = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart65/ensembl/Mm/Mm_Ensembl_exon.txt'
    firstRow = True
    filename = fn[:-4] + '.block.txt'
    ea = export.ExportFile(filename)
    found = False
    lines = 0
    exon_db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        gene, exonid, chromosome, strand, start, stop, a, b, c, d = string.split(data, '\t')
        exonid = string.split(exonid, '.')[0]
        uid = gene + ':' + exonid
        if lines > 0:
            try:
                exon_db[(uid, strand)].append(int(start))
                exon_db[(uid, strand)].append(int(stop))
            except Exception:
                exon_db[(uid, strand)] = [
                 int(start)]
                exon_db[(uid, strand)].append(int(stop))

        lines += 1

    print len(exon_db)
    for uid, strand in exon_db:
        exon_db[(uid, strand)].sort()
        if strand == '-':
            exon_db[(uid, strand)].reverse()
        start = str(exon_db[(uid, strand)][0])
        stop = str(exon_db[(uid, strand)][1])
        coord = [start, stop]
        coord.sort()
        ea.write(uid + '\t' + strand + '\t' + coord[0] + '\t' + coord[1] + '\n')

    ea.close()


def combineVariants(fn):
    firstRow = True
    filename = fn[:-4] + '.gene-level.txt'
    ea = export.ExportFile(filename)
    found = False
    lines = 0
    gene_db = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        gene = t[9]
        if lines == 0:
            header = [
             'UID'] + t[16:]
            header = string.join(header, '\t') + '\n'
            ea.write(header)
            lines += 1
        else:
            var_calls = map(float, t[16:])
            if gene in gene_db:
                count_sum_array = gene_db[gene]
                count_sum_array = [ sum(value) for value in zip(*[count_sum_array, var_calls]) ]
                gene_db[gene] = count_sum_array
            else:
                gene_db[gene] = var_calls

    for gene in gene_db:
        var_calls = gene_db[gene]
        var_calls2 = []
        for i in var_calls:
            if i == 0:
                var_calls2.append('0')
            else:
                var_calls2.append('1')

        ea.write(gene + '\t' + string.join(var_calls2, '\t') + '\n')

    ea.close()


def compareFusions(fn):
    firstRow = True
    filename = fn[:-4] + '.matrix.txt'
    ea = export.ExportFile(filename)
    found = False
    lines = 0
    fusion_db = {}
    sample_list = []
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        if 'Gene_Fusion_Pair' in line:
            headers = string.split(data, '\t')[1:]
        try:
            sample, fusion = string.split(data, '\t')
            try:
                fusion_db[fusion].append(sample)
            except Exception:
                fusion_db[fusion] = [sample]

            if sample not in sample_list:
                sample_list.append(sample)
        except Exception:
            t = string.split(data, '\t')
            fusion = t[0]
            index = 0
            for i in t[1:]:
                if i == '1':
                    sample = headers[index]
                    try:
                        fusion_db[fusion].append(sample)
                    except Exception:
                        fusion_db[fusion] = [sample]

                    if sample not in sample_list:
                        sample_list.append(sample)
                index += 1

    fusion_db2 = []
    for fusion in fusion_db:
        samples = fusion_db[fusion]
        samples2 = []
        for s in sample_list:
            if s in samples:
                samples2.append('1')
            else:
                samples2.append('0')

        fusion_db[fusion] = samples2

    ea.write(string.join(['Fusion'] + sample_list, '\t') + '\n')
    for fusion in fusion_db:
        print [
         fusion]
        ea.write(fusion + '\t' + string.join(fusion_db[fusion], '\t') + '\n')

    ea.close()


def customCleanSupplemental(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        line = string.split(data, ', ')
        gene_data = []
        for gene in line:
            gene = string.replace(gene, ' ', '')
            if '/' in gene:
                genes = string.split(gene, '/')
                gene_data.append(genes[0])
                for i in genes[1:]:
                    gene_data.append(genes[0][:len(genes[1]) * -1] + i)

            elif '(' in gene:
                genes = string.split(gene[:-1], '(')
                gene_data += genes
            else:
                gene_data.append(gene)

        ea.write(string.join(gene_data, ' ') + '\n')

    ea.close()


def customCleanBinomial(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    from stats_scripts import statistics
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            headers = t
            firstRow = False
            ea.write(string.join(['uid'] + headers, '\t') + '\n')
        else:
            gene = t[0]
            values = map(float, t[1:])
            min_val = abs(min(values))
            values = map(lambda x: x + min_val, values)
            values = map(str, values)
            ea.write(string.join([gene] + values, '\t') + '\n')

    ea.close()


class MarkerFinderInfo():

    def __init__(self, gene, rho, tissue):
        self.gene = gene
        self.rho = rho
        self.tissue = tissue

    def Gene(self):
        return self.gene

    def Rho(self):
        return self.rho

    def Tissue(self):
        return self.tissue


def ReceptorLigandCellInteractions(species, lig_receptor_dir, cell_type_gene_dir):
    ligand_db = {}
    receptor_db = {}
    fn = filepath(lig_receptor_dir)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        ligand, receptor = string.split(data, '\t')
        if species == 'Mm':
            ligand = ligand[0] + string.lower(ligand[1:])
            receptor = receptor[0] + string.lower(receptor[1:])
        try:
            ligand_db[ligand].apepnd(receptor)
        except Exception:
            ligand_db[ligand] = [receptor]
        else:
            try:
                receptor_db[receptor].append(ligand)
            except Exception:
                receptor_db[receptor] = [ligand]

    firstRow = True
    filename = cell_type_gene_dir[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    found = False
    cell_specific_ligands = {}
    cell_specific_receptor = {}
    fn = filepath(cell_type_gene_dir)
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        gene, rho, tissue, notes, order = string.split(data, '\t')
        mf = MarkerFinderInfo(gene, rho, tissue)
        if gene in ligand_db:
            cell_specific_ligands[gene] = mf
        if gene in receptor_db:
            cell_specific_receptor[gene] = mf

    ligand_receptor_pairs = []
    for gene in cell_specific_ligands:
        receptors = ligand_db[gene]
        for receptor in receptors:
            if receptor in cell_specific_receptor:
                rmf = cell_specific_receptor[receptor]
                lmf = cell_specific_ligands[gene]
                gene_data = [gene, lmf.Tissue(), lmf.Rho(), receptor, rmf.Tissue(), rmf.Rho()]
                pair = (gene, receptor)
                if pair not in ligand_receptor_pairs:
                    ea.write(string.join(gene_data, '\t') + '\n')
                    ligand_receptor_pairs.append(pair)

    for receptor in cell_specific_receptor:
        ligands = receptor_db[receptor]
        for gene in ligands:
            if gene in cell_specific_ligands:
                rmf = cell_specific_receptor[receptor]
                lmf = cell_specific_ligands[gene]
                gene_data = [gene, lmf.Tissue(), lmf.Rho(), receptor, rmf.Tissue(), rmf.Rho()]
                pair = (gene, receptor)
                if pair not in ligand_receptor_pairs:
                    ea.write(string.join(gene_data, '\t') + '\n')
                    ligand_receptor_pairs.append(pair)

    ea.close()


def findReciprocal(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-filtered.txt'
    ea = export.ExportFile(filename)
    found = False
    gene_ko = {}
    gene_oe = {}
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            firstRow = False
            headers = {}
            TFs = {}
            i = 0
            for v in t[1:]:
                TF, direction = string.split(v, '-')
                headers[i] = (TF, direction, v)
                i += 1
                if v not in TFs:
                    f = filename[:-4] + '-' + v + '-up.txt'
                    tea = export.ExportFile(f)
                    TFs[v + '-up'] = tea
                    tea.write('GeneID\tEn\n')
                    f = filename[:-4] + '-' + v + '-down.txt'
                    tea = export.ExportFile(f)
                    TFs[v + '-down'] = tea
                    tea.write('GeneID\tEn\n')

        else:
            values = map(float, t[1:])
            gene = t[0]
            i = 0
            for v in values:
                TF, direction, name = headers[i]
                if 'KO' in direction:
                    if v > 1:
                        gene_ko[(gene, TF, 1)] = []
                        tea = TFs[(name + '-up')]
                        tea.write(gene + '\tEn\n')
                    else:
                        gene_ko[(gene, TF, -1)] = []
                        tea = TFs[(name + '-down')]
                        tea.write(gene + '\tEn\n')
                if 'OE' in direction:
                    if v > 1:
                        gene_oe[(gene, TF, 1)] = []
                        tea = TFs[(name + '-up')]
                        tea.write(gene + '\tEn\n')
                    else:
                        gene_oe[(gene, TF, -1)] = []
                        tea = TFs[(name + '-down')]
                        tea.write(gene + '\tEn\n')
                i += 1

    print len(gene_oe)
    for gene, TF, direction in gene_oe:
        alt_dir = direction * -1
        if (gene, TF, alt_dir) in gene_ko:
            ea.write(string.join([TF, gene, str(direction)], '\t') + '\n')

    ea.close()
    for TF in TFs:
        TFs[TF].close()


def effectsPrioritization(filename):
    fn = filepath(filename)
    firstRow = True
    filename = filename[:-4] + '-new.txt'
    ea = export.ExportFile(filename)
    from stats_scripts import statistics
    found = False
    for line in open(fn, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if firstRow:
            headers = t[1:]
            firstRow = False
        else:
            gene = t[0]
            values = map(float, t[1:])
            max_val = abs(max(values))
            max_header = headers[values.index(max_val)]
            ea.write(gene + '\t' + max_header + '\t' + str(max_val) + '\n')

    ea.close()


def simpleCombine(folder):
    filename = folder + '/combined/combined.txt'
    ea = export.ExportFile(filename)
    headers = ['UID']
    data_db = {}
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file:
            fn = filepath(folder + '/' + file)
            print fn
            firstRow = True
            for line in open(fn, 'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data, '\t')
                if firstRow:
                    for i in t[1:]:
                        headers.append(i + '.' + file[:-4])

                    firstRow = False
                else:
                    gene = t[0]
                    try:
                        data_db[gene] += t[1:]
                    except Exception:
                        data_db[gene] = t[1:]

    len_db = {}
    ea.write(string.join(headers, '\t') + '\n')
    for gene in data_db:
        if len(data_db[gene]) == len(headers) - 1:
            values = map(float, data_db[gene])
            count = 0
            for i in values:
                if i > 0.9:
                    count += 1

            if count > 7:
                ea.write(string.join([gene] + data_db[gene], '\t') + '\n')
        len_db[len(data_db[gene])] = []

    print len(len_db)
    for i in len_db:
        print i

    ea.close()


def simpleCombineFiles(folder, elite_output=True, uniqueOnly=False):
    filename = folder + '/combined/combined.txt'
    ea = export.ExportFile(filename)
    files = UI.read_directory(folder)
    unique_entries = []
    firstRow = True
    for file in files:
        if '.txt' in file:
            fn = filepath(folder + '/' + file)
            firstRow = True
            for line in open(fn, 'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data, '\t')
                if elite_output:
                    if firstRow:
                        t = 'UID\tSystemCode\tCategory'
                        ea.write(t + '\n')
                        firstRow = False
                    elif uniqueOnly:
                        if t[0] not in unique_entries:
                            ea.write(t[0] + '\t' + t[0] + '\t\n')
                        unique_entries.append(t[0])
                    else:
                        ea.write(string.join([string.split(t[0], '|')[0], 'Ae', string.replace(file[:-4], 'PSI.', '')], '\t') + '\n')
                elif firstRow:
                    t.append('Comparison')
                    ea.write(string.join(t, '\t') + '\n')
                    firstRow = False
                else:
                    t.append(file[:-4])
                    ea.write(string.join(t, '\t') + '\n')

    ea.close()


def combineDEGs(folder):
    filename = folder + '/combined_unique.txt'
    ea = export.ExportFile(filename)
    files = UI.read_directory(folder)
    for file in files:
        header = True
        if 'GE.' in file:
            for line in open(folder + '/' + file, 'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data, '\t')
                if header:
                    header = False
                else:
                    log_fold = abs(float(t[2]))
                    fold = math.pow(2, log_fold)
                    if fold > 1.2:
                        ea.write(t[0] + '\t' + t[1] + '\t' + file[3:-4] + '\n')

    ea.close()


def combinePSIFiles(folder, gene_file, PSI_annotations=None):
    filename = folder + '/combined/combined_unique.txt'
    psi_annotations = {}
    ensembl_symbol = {}
    if PSI_annotations != None:
        for line in open(PSI_annotations, 'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data, '\t')
            uid = t[0]
            symbol = t[2]
            altExons = t[4]
            protein_predictions = t[5]
            clusterID = t[(-3)]
            annotations = (symbol, clusterID, altExons, protein_predictions)
            psi_annotations[uid] = annotations
            try:
                a = string.split(uid, ':')
                ensembl_symbol[a[1]] = symbol
            except:
                pass

    cegs = {}
    for line in open(gene_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        geneID, counts, tissues = string.split(data, '\t')
        cegs[geneID] = tissues

    ea = export.ExportFile(filename)
    files = UI.read_directory(folder)
    cluster_results = {}
    firstRow = True
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            fn = filepath(folder + '/' + file)
            firstRow = True
            for line in open(fn, 'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data, '\t')
                uid = t[0]
                try:
                    symbol, clusterID, altExons, protein_predictions = psi_annotations[uid]
                except Exception:
                    uid_temp = string.split(uid, ':')
                    symbol = ensembl_symbol[uid_temp[1]]
                    uid = string.join([symbol] + uid_temp[1:], ':')
                    symbol, clusterID, altExons, protein_predictions = psi_annotations[uid]
                else:
                    if clusterID in cluster_results:
                        db = cluster_results[clusterID]
                        db[file] = (uid, symbol, altExons, protein_predictions)
                    else:
                        db = {}
                        db[file] = (
                         uid, symbol, altExons, protein_predictions)
                        cluster_results[clusterID] = db

    for clusterID in cluster_results:
        tissues = map(str, cluster_results[clusterID])
        tissues = string.join(tissues, '|')
        for tissue in cluster_results[clusterID]:
            uid, symbol, altExons, protein_predictions = cluster_results[clusterID][tissue]

        try:
            tissue_cegs = cegs[symbol]
        except:
            tissue_cegs = ''
        else:
            t = [
             uid, symbol, altExons, protein_predictions, tissues, tissue_cegs]
            ea.write(string.join(t, '\t') + '\n')

    ea.close()
    return


def evaluateMultiLinRegulatoryStructure(all_genes_TPM, MarkerFinder, SignatureGenes, state, query=None):
    """Predict multi-lineage cells and their associated coincident lineage-defining TFs"""
    ICGS_State_as_Row = True
    matrix, column_header, row_header, dataset_name, group_db = importData(all_genes_TPM)
    group_index = {}
    all_indexes = []
    for sampleName in group_db:
        ICGS_state = group_db[sampleName][0]
        try:
            group_index[ICGS_state].append(column_header.index(sampleName))
        except Exception:
            group_index[ICGS_state] = [column_header.index(sampleName)]
        else:
            all_indexes.append(column_header.index(sampleName))

    for ICGS_state in group_index:
        group_index[ICGS_state].sort()

    all_indexes.sort()

    def importGeneLists(fn):
        genes = {}
        for line in open(fn, 'rU').xreadlines():
            data = cleanUpLine(line)
            gene, cluster = string.split(data, '\t')[0:2]
            genes[gene] = cluster

        return genes

    def importMarkerFinderHits(fn):
        genes = {}
        skip = True
        for line in open(fn, 'rU').xreadlines():
            data = cleanUpLine(line)
            if skip:
                skip = False
            else:
                gene, symbol, rho, ICGS_State = string.split(data, '\t')
                if float(rho) > 0.0:
                    genes[gene] = (
                     float(rho), ICGS_State)
                    genes[symbol] = (float(rho), ICGS_State)

        return genes

    def importQueryDataset(fn):
        matrix, column_header, row_header, dataset_name, group_db = importData(fn)
        return (
         matrix, column_header, row_header, dataset_name, group_db)

    signatureGenes = importGeneLists(SignatureGenes)
    markerFinderGenes = importMarkerFinderHits(MarkerFinder)
    index = 0
    expressedGenesPerState = {}

    def freqCutoff(x, cutoff):
        if x > cutoff:
            return 1
        else:
            return 0

    for row in matrix:
        ICGS_state_gene_frq = {}
        gene = row_header[index]
        for ICGS_state in group_index:
            state_values = map(lambda i: row[i], group_index[ICGS_state])

            def freqCheck(x):
                if x > 1:
                    return 1
                else:
                    return 0

            expStateCells = sum(map(lambda x: freqCheck(x), state_values))
            statePercentage = float(expStateCells) / len(group_index[ICGS_state])
            ICGS_state_gene_frq[ICGS_state] = statePercentage

        multilin_frq = ICGS_state_gene_frq[state]
        datasets_values = map(lambda i: row[i], all_indexes)
        all_cells_frq = sum(map(lambda x: freqCheck(x), datasets_values)) / (len(datasets_values) * 1.0)
        all_states_frq = map(lambda x: ICGS_state_gene_frq[x], ICGS_state_gene_frq)
        all_states_frq.sort()
        rank = all_states_frq.index(multilin_frq)
        states_expressed = sum(map(lambda x: freqCutoff(x, 0.5), all_states_frq)) / (len(all_states_frq) * 1.0)
        if multilin_frq > 0.25 and rank > 0:
            if 'Rik' not in gene and 'Gm' not in gene:
                if gene in signatureGenes:
                    if ICGS_State_as_Row:
                        ICGS_State = signatureGenes[gene]
                    if gene in markerFinderGenes:
                        if ICGS_State_as_Row == False:
                            rho, ICGS_State = markerFinderGenes[gene]
                        else:
                            rho, ICGS_Cell_State = markerFinderGenes[gene]
                        score = int(rho * 100 * multilin_frq) * (float(rank) / len(all_states_frq))
                        try:
                            expressedGenesPerState[ICGS_State].append((score, gene))
                        except Exception:
                            expressedGenesPerState[ICGS_State] = [(score, gene)]

        index += 1

    if query != None:
        matrix, column_header, row_header, dataset_name, group_db = importQueryDataset(query)
    createPseudoCell = True
    representativeMarkers = {}
    for ICGS_State in expressedGenesPerState:
        expressedGenesPerState[ICGS_State].sort()
        expressedGenesPerState[ICGS_State].reverse()
        if '1Multi' not in ICGS_State:
            markers = expressedGenesPerState[ICGS_State][:5]
            print ICGS_State, ':', string.join(map(lambda x: x[1], list(markers)), ', ')
            if createPseudoCell:
                for gene in markers:

                    def getBinary(x):
                        if x > 1:
                            return 1
                        else:
                            return 0

                    if gene[1] in row_header:
                        row_index = row_header.index(gene[1])
                        binaryValues = map(lambda x: getBinary(x), matrix[row_index])
                        try:
                            representativeMarkers[ICGS_State].append(binaryValues)
                        except Exception:
                            representativeMarkers[ICGS_State] = [binaryValues]

            else:
                representativeMarkers[ICGS_State] = markers[0][(-1)]

    for ICGS_State in representativeMarkers:
        if createPseudoCell:
            signature_values = representativeMarkers[ICGS_State]
            signature_values = [ int(numpy.median(value)) for value in zip(*signature_values) ]
            representativeMarkers[ICGS_State] = signature_values
        else:
            gene = representativeMarkers[ICGS_State]
            row_index = row_header.index(gene)
            gene_values = matrix[row_index]
            representativeMarkers[ICGS_State] = gene_values

    expressedStatesPerCell = {}
    for ICGS_State in representativeMarkers:
        gene_values = representativeMarkers[ICGS_State]
        index = 0
        for cell in column_header:
            log2_tpm = gene_values[index]
            if log2_tpm >= 1:
                try:
                    expressedStatesPerCell[cell].append(ICGS_State)
                except Exception:
                    expressedStatesPerCell[cell] = [ICGS_State]

            index += 1

    cell_mutlilin_ranking = []
    for cell in expressedStatesPerCell:
        lineageCount = expressedStatesPerCell[cell]
        cell_mutlilin_ranking.append((len(lineageCount), cell))

    cell_mutlilin_ranking.sort()
    cell_mutlilin_ranking.reverse()
    for cell in cell_mutlilin_ranking:
        print cell[0], cell[1], string.join(expressedStatesPerCell[cell[1]], '|')

    return


def compareGenomicLocationAndICGSClusters():
    species = 'Mm'
    array_type = 'RNASeq'
    from build_scripts import EnsemblImport
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species, array_type, 'key_by_array')
    markerfinder = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/ExpressionOutput/MarkerFinder/AllCorrelationsAnnotated-ProteinCodingOnly.txt'
    eo = export.ExportFile(markerfinder[:-4] + '-bidirectional_promoters.txt')
    firstRow = True
    chr_cellTypeSpecific = {}
    for line in open(markerfinder, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        symbol = t[1]
        ensembl = t[0]
        try:
            rho = float(t[6])
        except Exception:
            pass
        else:
            cellType = t[7]
            if firstRow:
                firstRow = False
            elif ensembl in gene_location_db and rho > 0.2:
                chr, strand, start, end = gene_location_db[ensembl]
                start = int(start)
                end = int(end)
                try:
                    db = chr_cellTypeSpecific[(chr, cellType)]
                    try:
                        db[strand].append([start, end, symbol, ensembl])
                    except Exception:
                        db[strand] = [[start, end, symbol, ensembl]]

                except Exception:
                    db = {}
                    db[strand] = [
                     [
                      start, end, symbol, ensembl]]
                    chr_cellTypeSpecific[(chr, cellType)] = db

    bidirectional = {}
    eo.write(string.join(['CellType', 'Chr', 'Ensembl1', 'Symbol1', 'Start1', 'End1', 'Strand1', 'Ensembl2', 'Symbol2', 'Start2', 'End2', 'Strand2'], '\t') + '\n')
    for chr, cellType in chr_cellTypeSpecific:
        db = chr_cellTypeSpecific[(chr, cellType)]
        if len(db) > 1:
            for start, end, symbol, ens in db['+']:
                for start2, end2, symbol2, ens2 in db['-']:
                    if abs(start - end2) < 100000 and start > end2:
                        eo.write(string.join([cellType, chr, ens, symbol, str(start), str(end), '+', ens2, symbol2, str(end2), str(start2), '-'], '\t') + '\n')
                        try:
                            bidirectional[(chr, cellType)].append([start, end, symbol, ens, start2, end2, symbol2, ens2])
                        except Exception:
                            bidirectional[(chr, cellType)] = [[start, end, symbol, ens, start2, end2, symbol2, ens2]]

    eo.close()


def filterCountsFile(filename):
    fn = filepath(filename)
    firstRow = True

    def countif(value, cutoff=9):
        if float(value) > cutoff:
            return 1
        else:
            return 0

    header = True
    unique_genes = {}
    ea = export.ExportFile(filename[:-4] + '-filtered.txt')
    for line in open(fn, 'rU').xreadlines():
        if header:
            header = False
            ea.write(line)
        else:
            data = line.rstrip()
            t = string.split(data, '\t')
            gene = string.split(t[0], ':')[0]
            unique_genes[gene] = []
            expressedSamples = map(countif, t[1:])
            if sum(expressedSamples) > 2:
                ea.write(line)

    ea.close()
    print len(unique_genes), 'unique genes.'


def filterPSIValues(filename):
    fn = filepath(filename)
    firstRow = True
    header = True
    rows = 0
    filtered = 0
    new_file = filename[:-4] + '-75p.txt'
    new_file_clust = new_file[:-4] + '-clustID.txt'
    ea = export.ExportFile(new_file)
    eac = export.ExportFile(new_file_clust)
    added = []
    for line in open(fn, 'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data, '\t')
        if header:
            header = False
            t = [t[1]] + t[8:]
            header_length = len(t) - 1
            minimum_values_present = int(0.75 * int(header_length))
            not_detected = header_length - minimum_values_present
            new_line = string.join(t, '\t') + '\n'
            ea.write(new_line)
        else:
            cID = t[5]
            t = [t[1]] + t[8:]
            missing_values_at_the_end = header_length + 1 - len(t)
            missing = missing_values_at_the_end + t.count('')
            if missing < not_detected:
                added.append(cID)
                new_line = string.join(t, '\t') + '\n'
                ea.write(new_line)
                eac.write(t[0] + '\t' + cID + '\n')
                filtered += 1
        rows += 1

    print rows, filtered
    ea.close()
    eac.close()


def removeRedundantCluster(filename, clusterID_file):
    from scipy import stats
    import ExpressionBuilder
    sort_col = 0
    export_count = 0
    ExpressionBuilder.exportSorted(filename, sort_col, excludeHeader=True)
    new_file = filename[:-4] + '-unique.txt'
    ea = export.ExportFile(new_file)
    event_clusterID_db = {}
    for line in open(clusterID_file, 'rU').xreadlines():
        data = line.rstrip()
        eventID, clusterID = string.split(data, '\t')
        event_clusterID_db[eventID] = clusterID

        def compareEvents(events_to_compare, export_count):
            if len(events_to_compare) == 1:
                ea.write(events_to_compare[0][(-1)])
                export_count += 1
            else:
                exclude = {}
                compared = {}
                for event1 in events_to_compare:
                    if event1[0] not in exclude:
                        ea.write(event1[(-1)])
                        exclude[event1[0]] = []
                        export_count += 1
                    for event2 in events_to_compare:
                        if event2[0] not in exclude:
                            if event1[0] != event2[0] and (event1[0], event2[0]) not in compared:
                                uid1, values1, line1 = event1
                                uid2, values2, line2 = event2
                                coefr = numpy.ma.corrcoef(values1, values2)
                                rho = coefr[0][1]
                                if rho > 0.6 or rho < -0.6:
                                    exclude[event2[0]] = []
                                compared[(event1[0], event2[0])] = []
                                compared[(event2[0], event1[0])] = []

            for event in events_to_compare:
                if event[0] not in exclude:
                    ea.write(event[(-1)])
                    exclude.append(event[0])
                    export_count += 1

            return export_count

    header = True
    rows = 0
    filtered = 0
    prior_cID = 0
    events_to_compare = []
    for line in open(filename, 'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data, '\t')
        if header:
            ea.write(line)
            header_row = t
            header = False
        else:
            uid = t[0]
            cID = event_clusterID_db[uid]
            empty_offset = len(header_row) - len(t)
            t += [''] * empty_offset
            values = [ '0.000101' if x == '' else x for x in t[1:] ]
            values = map(float, values)
            values = numpy.ma.masked_values(values, 0.000101)
            if prior_cID == 0:
                prior_cID = cID
            if cID == prior_cID:
                events_to_compare.append((uid, values, line))
            else:
                export_count = compareEvents(events_to_compare, export_count)
                events_to_compare = [(uid, values, line)]
            prior_cID = cID

    if len(events_to_compare) > 0:
        export_count = compareEvents(events_to_compare, export_count)
    ea.close()
    print export_count, 'Non-redundant splice-events exported'


def convertToGOElite(folder):
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file:
            gene_count = 0
            up_count = 0
            down_count = 0
            new_filename = string.split(file[3:], '_')[0] + '.txt'
            ea = export.ExportFile(folder + '/GO-Elite/' + new_filename)
            fn = folder + '/' + file
            ea.write('GeneID\tSystemCode\n')
            firstLine = True
            for line in open(fn, 'rU').xreadlines():
                if firstLine:
                    firstLine = False
                    continue
                data = line.rstrip()
                t = string.split(data, '\t')
                if ':' in t[0]:
                    ea.write(string.split(t[0], ':')[0] + '\tSy\n')
                else:
                    gene_count += 1
                    if '-' in t[2]:
                        down_count += 1
                    else:
                        up_count += 1

            ea.close()
            print file, '\t', gene_count, '\t', up_count, '\t', down_count


def geneExpressionSummary(folder):
    import collections
    event_db = collections.OrderedDict()
    groups_list = ['']
    files = UI.read_directory(folder)
    for file in files:
        if '.txt' in file and 'GE.' in file:
            ls = []
            event_db[file[:-4]] = ls
            groups_list.append(file[:-4])
            fn = folder + '/' + file
            firstLine = True
            for line in open(fn, 'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data, '\t')
                if firstLine:
                    fold_index = t.index('LogFold')
                    firstLine = False
                    continue
                uid = t[0]
                if float(t[fold_index]) > 0:
                    fold_dir = 1
                else:
                    fold_dir = -1
                ls.append((uid, fold_dir))

    for file in event_db:
        print file, '\t', len(event_db[file])


def compareEventLists(folder):
    import collections
    event_db = collections.OrderedDict()
    groups_list = ['']
    files = UI.read_directory(folder)
    file_headers = {}
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            ls = {}
            event_db[file[:-4]] = ls
            groups_list.append(file[:-4])
            fn = folder + '/' + file
            firstLine = True
            for line in open(fn, 'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data, '\t')
                if firstLine:
                    file_headers[file[:-4]] = t
                    cid = t.index('ClusterID')
                    try:
                        event_index = t.index('Event-Direction')
                    except:
                        try:
                            event_index = t.index('Inclusion-Junction')
                        except:
                            print file, 'Event-Direction error'
                            sys.exit()

                    firstLine = False
                    continue
                uid = t[0]
                uid = string.split(uid, '|')[0]
                if 'U2AF1-l' in file or 'U2AF1-E' in file:
                    if t[2] == 'inclusion':
                        ls[(uid, t[event_index])] = t
                else:
                    ls[(uid, t[event_index])] = t

    def convertEvents(events):
        opposite_events = []
        for event, direction in events:
            if direction == 'exclusion':
                direction = 'inclusion'
            else:
                direction = 'exclusion'
            opposite_events.append((event, direction))

        return opposite_events

    ea1 = export.ExportFile(folder + '/overlaps-same-direction.txt')
    ea2 = export.ExportFile(folder + '/overlaps-opposite-direction.txt')
    ea3 = export.ExportFile(folder + '/concordance.txt')
    ea1.write(string.join(groups_list, '\t') + '\n')
    ea2.write(string.join(groups_list, '\t') + '\n')
    ea3.write(string.join(groups_list, '\t') + '\n')
    comparison_db = {}
    best_hits = {}
    for comparison1 in event_db:
        events1 = event_db[comparison1]
        hits1 = [comparison1]
        hits2 = [comparison1]
        hits3 = [comparison1]
        best_hits[comparison1] = []
        for comparison2 in event_db:
            events2 = event_db[comparison2]
            events3 = convertEvents(events2)
            overlapping_events = list(set(events1).intersection(events2))
            overlap = len(overlapping_events)
            inverse_overlap = len(set(events1).intersection(events3))
            min_events1 = min([len(events1), len(events2)])
            min_events2 = min([len(events1), len(events3)])
            denom = overlap + inverse_overlap
            if denom == 0:
                denom = 1e-05
            if min_events1 == 0:
                min_events1 = 1
            if overlap + inverse_overlap < 20:
                hits1.append('0.5')
                hits2.append('0.5')
                hits3.append('0.5|0.5')
            else:
                hits1.append(str(1.0 * overlap / min_events1))
                hits2.append(str(1.0 * inverse_overlap / min_events1))
                hits3.append(str(1.0 * overlap / denom) + '|' + str(1.0 * inverse_overlap / denom) + ':' + str(overlap + inverse_overlap))
                if 'Leu' not in comparison2:
                    comp_name = string.split(comparison2, '_vs')[0]
                    best_hits[comparison1].append([abs(1.0 * overlap / denom), 'cor', comp_name])
                    best_hits[comparison1].append([abs(1.0 * inverse_overlap / denom), 'anti', comp_name])
            if comparison1 != comparison2:
                if len(overlapping_events) > 0:
                    pass
                overlapping_events.sort()
                for event in overlapping_events:
                    vals = string.join([event[0], comparison1] + event_db[comparison1][event] + [comparison2] + event_db[comparison2][event], '\t')

        ea1.write(string.join(hits1, '\t') + '\n')
        ea2.write(string.join(hits2, '\t') + '\n')
        ea3.write(string.join(hits3, '\t') + '\n')

    ea1.close()
    ea2.close()
    ea3.close()
    for comparison in best_hits:
        best_hits[comparison].sort()
        best_hits[comparison].reverse()
        hits = best_hits[comparison][:10]
        hits2 = []
        for score, dir, comp in hits:
            h = str(score)[:4] + '|' + dir + '|' + comp
            hits2.append(h)

        print comparison, '\t', string.join(hits2, ', ')


def convertGroupsToBinaryMatrix(groups_file, sample_order, cellHarmony=False):
    eo = export.ExportFile(groups_file[:-4] + '-matrix.txt')
    print groups_file[:-4] + '-matrix.txt'
    firstRow = True
    samples = []
    for line in open(sample_order, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if 'row_clusters-flat' in t:
            samples = []
            samples1 = t[2:]
            for name in samples1:
                if ':' in name:
                    group, name = string.split(name, ':')
                samples.append(name)

            if cellHarmony == False:
                break
        elif 'column_clusters-flat' in t and cellHarmony:
            clusters = t[2:]
        elif groups_file == sample_order:
            samples.append(t[0])
        elif firstRow:
            samples = t[1:]
            firstRow = False

    import collections
    sample_groups = collections.OrderedDict()
    for line in open(groups_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        sample, groupNum, groupName = t[:3]
        if cellHarmony == False:
            if sample in samples:
                si = samples.index(sample)
                try:
                    sample_groups[groupName][si] = '1'
                except Exception:
                    sample_groups[groupName] = [
                     '0'] * len(samples)
                    sample_groups[groupName][si] = '1'

        else:
            sample_groups[groupNum] = groupName

    if cellHarmony:
        i = 0
        for sample in samples1:
            cluster = clusters[i]
            group_name = sample_groups[cluster]
            eo.write(sample + '\t' + cluster + '\t' + group_name + '\n')
            i += 1

        eo.close()
    else:
        eo.write(string.join(['GroupName'] + samples, '\t') + '\n')
        for group in sample_groups:
            eo.write(string.join([group] + sample_groups[group], '\t') + '\n')

        eo.close()


def returnIntronJunctionRatio(counts_file, species='Mm'):
    eo = export.ExportFile(counts_file[:-4] + '-intron-ratios.txt')
    header = True
    prior_gene = []
    exon_junction_values = []
    intron_junction_values = []
    eoi = export.ExportFile(counts_file[:-4] + '-intron-ratios-gene.txt')
    rows = 0
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species, ('hide', 'Ensembl-Symbol'))

    def logratio(list):
        try:
            return list[0] / list[1]
        except Exception:
            return 0

    for line in open(counts_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        junctionID = t[0]
        if header:
            eoi.write(line)
            samples = t[1:]
            global_intron_ratios = {}
            i = 0
            for val in samples:
                global_intron_ratios[i] = []
                i += 1

            header = False
            continue
        else:
            uid, coords = string.split(junctionID, '=')
            genes = string.split(uid, ':')
            if len(genes) > 2:
                trans_splicing = True
            else:
                trans_splicing = False
            coords = string.split(coords, ':')[1]
            coords = string.split(coords, '-')
            coords = map(int, coords)
            coord_diff = abs(coords[1] - coords[0])
        gene = string.split(junctionID, ':')[0]
        rows += 1
        if rows == 1:
            prior_gene = gene
        if gene != prior_gene:
            if len(intron_junction_values) == 0:
                pass
            else:
                intron_junction_values_original = list(intron_junction_values)
                exon_junction_values_original = list(exon_junction_values)
                intron_junction_values = [ sum(i) for i in zip(*intron_junction_values) ]
                exon_junction_values = [ sum(i) for i in zip(*exon_junction_values) ]
                intron_ratios = [ logratio(value) for value in zip(*[intron_junction_values, exon_junction_values]) ]
                intron_ratios2 = []
                if prior_gene in gene_to_symbol:
                    symbol = gene_to_symbol[prior_gene][0]
                else:
                    symbol = prior_gene
                i = 0
                if symbol == 'Pi4ka':
                    print samples[482:487]
                    for x in exon_junction_values_original:
                        print x[482:487]

                    print exon_junction_values[482:487]
                    print intron_ratios[482:487]
                for val in intron_ratios:
                    if exon_junction_values[i] > 9:
                        if val > 0:
                            if intron_junction_values[i] > 9:
                                intron_ratios2.append(val)
                            else:
                                intron_ratios2.append(0)
                        else:
                            intron_ratios2.append(0)
                    else:
                        intron_ratios2.append('')
                    i += 1

            eoi.write(string.join([symbol] + map(str, intron_ratios2), '\t') + '\n')
            i = 0
            for val in intron_ratios:
                if exon_junction_values[i] != 0:
                    global_intron_ratios[i].append(intron_ratios[i])
                i += 1

            exon_junction_values = []
            intron_junction_values = []
            prior_gene = gene
        values = map(float, t[1:])
        if 'I' in junctionID and '_' not in junctionID and coord_diff == 1 and trans_splicing == False:
            intron_junction_values.append(values)
            exon_junction_values.append(values)
        elif trans_splicing == False:
            exon_junction_values.append(values)

    print rows, 'processed'
    import numpy
    i = 0
    global_intron_ratios_values = []
    for val in samples:
        global_intron_ratios_values.append(100 * numpy.mean(global_intron_ratios[i]))
        i += 1

    eo.write(string.join(['UID'] + samples, '\t') + '\n')
    eo.write(string.join(['Global-Intron-Retention-Ratio'] + map(str, global_intron_ratios_values), '\t') + '\n')
    eo.close()
    eoi.close()


def convertSymbolLog(input_file, ensembl_symbol):
    gene_symbol_db = {}
    for line in open(ensembl_symbol, 'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl, symbol = string.split(data, '\t')
        gene_symbol_db[ensembl] = symbol

    convert = False
    eo = export.ExportFile(input_file[:-4] + '-log2.txt')
    header = 0
    added_symbols = []
    not_found = []
    for line in open(input_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data, '\t')
        gene = values[0]
        if header == 0:
            data = cleanUpLine(line)
            headers = []
            values = string.split(data, '\t')
            for v in values:
                if 'exp.' in v:
                    headers.append(string.split(v, '.exp.')[0])
                else:
                    headers.append(v)

            eo.write(string.join(headers, '\t') + '\n')
        header += 1
        if gene in gene_symbol_db:
            symbol = gene_symbol_db[gene]
            if symbol not in added_symbols:
                added_symbols.append(symbol)
                values = map(lambda x: math.log(float(x) + 1, 2), values[1:])
                if max(values) > 0.5:
                    values = map(lambda x: str(x)[:5], values)
                    eo.write(string.join([symbol] + values, '\t') + '\n')
        elif convert == False and header > 1:
            values = map(lambda x: math.log(float(x) + 1, 2), values[1:])
            if max(values) > 0.5:
                values = map(lambda x: str(x)[:5], values)
                eo.write(string.join([gene] + values, '\t') + '\n')
        else:
            not_found.append(gene)

    print len(not_found), not_found[:10]
    eo.close()


def convertXenaBrowserIsoformDataToStandardRatios(input_file):
    eo = open(input_file[:-4] + '-log2.txt', 'w')
    header = 0
    count = 0
    for line in open(input_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data, '\t')
        uid = string.split(values[0], '.')[0]
        isoform = values[0]
        if header == 0:
            eo.write(line)
            header += 1
        else:
            values = map(lambda x: math.pow(2, float(x)), values[1:])
            values = map(lambda x: math.log(float(x) + 1, 2), values)

            def percentExp(x):
                if x > 1:
                    return 1
                else:
                    return 0

            counts = map(lambda x: percentExp(x), values)
            if sum(counts) / (len(values) * 1.0) > 0.1:
                values = map(str, values)
                values = string.join([uid] + values, '\t')
                eo.write(values + '\n')
                count += 1

    eo.close()
    print count, 'genes written'


def outputForGOElite(folds_dir):
    matrix, column_header, row_header, dataset_name, group_db = importData(folds_dir, Normalize=False)
    matrix = zip(*matrix)
    ci = 0
    root_dir = findParentDir(folds_dir)
    for group_data in matrix:
        group_name = column_header[ci]
        eo = export.ExportFile(root_dir + '/folds/' + group_name + '.txt')
        gi = 0
        eo.write('geneID\tSy\tlog2-fold\n')
        for fold in group_data:
            gene = row_header[gi]
            if fold > 0:
                eo.write(gene + '\tSy\t' + str(fold) + '\n')
            gi += 1

        eo.close()
        ci += 1


def transposeMatrix(input_file):
    arrays = []
    eo = export.ExportFile(input_file[:-4] + '-transposed.txt')
    for line in open(input_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data, '\t')
        arrays.append(values)

    t_arrays = zip(*arrays)
    for t in t_arrays:
        eo.write(string.join(t, '\t') + '\n')

    eo.close()


def simpleStatsSummary(input_file):
    cluster_counts = {}
    header = True
    for line in open(input_file, 'rU').xreadlines():
        data = cleanUpLine(line)
        if header:
            header = False
        else:
            sample, cluster, counts = string.split(data, '\t')
            try:
                cluster_counts[cluster].append(float(counts))
            except Exception:
                cluster_counts[cluster] = [float(counts)]

    for cluster in cluster_counts:
        avg = statistics.avg(cluster_counts[cluster])
        stdev = statistics.stdev(cluster_counts[cluster])
        print cluster + '\t' + str(avg) + '\t' + str(stdev)


def latteralMerge(file1, file2):
    import collections
    cluster_db = collections.OrderedDict()
    eo = export.ExportFile(file2[:-4] + 'combined.txt')
    for line in open(file1, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        cluster_db[t[0]] = t

    for line in open(file2, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if t[0] in cluster_db:
            t1 = cluster_db[t[0]]
            eo.write(string.join(t1 + t[2:], '\t') + '\n')

    eo.close()


def removeMarkerFinderDoublets(heatmap_file, diff=1):
    matrix, column_header, row_header, dataset_name, group_db, priorColumnClusters, priorRowClusters = remoteImportData(heatmap_file)
    priorRowClusters.reverse()
    if len(priorColumnClusters) == 0:
        for c in column_header:
            cluster = string.split(c, ':')[0]
            priorColumnClusters.append(cluster)

        for r in row_header:
            cluster = string.split(r, ':')[0]
            priorRowClusters.append(cluster)

    import collections
    cluster_db = collections.OrderedDict()
    i = 0
    for cluster in priorRowClusters:
        try:
            cluster_db[cluster].append(matrix[i])
        except:
            cluster_db[cluster] = [matrix[i]]
        else:
            i += 1

    transposed_data_matrix = []
    clusters = []
    for cluster in cluster_db:
        cluster_cell_means = numpy.mean(cluster_db[cluster], axis=0)
        cluster_db[cluster] = cluster_cell_means
        transposed_data_matrix.append(cluster_cell_means)
        if cluster not in clusters:
            clusters.append(cluster)

    transposed_data_matrix = zip(*transposed_data_matrix)
    i = 0
    cell_max_scores = []
    cell_max_score_db = collections.OrderedDict()
    for cell_scores in transposed_data_matrix:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ci = clusters.index(cluster)
        cell_state_score = cell_scores[ci]
        alternate_state_scores = []
        for score in cell_scores:
            if score != cell_state_score:
                alternate_state_scores.append(score)

        alt_max_score = max(alternate_state_scores)
        alt_sum_score = sum(alternate_state_scores)
        cell_max_scores.append([cell_state_score, alt_max_score, alt_sum_score])
        try:
            cell_max_score_db[cluster].append([cell_state_score, alt_max_score, alt_sum_score])
        except:
            cell_max_score_db[cluster] = [[cell_state_score, alt_max_score, alt_sum_score]]
        else:
            i += 1

    for cluster in cell_max_score_db:
        cluster_cell_means = numpy.median(cell_max_score_db[cluster], axis=0)
        cell_max_score_db[cluster] = cluster_cell_means

    i = 0
    print len(cell_max_scores)
    keep = ['row_clusters-flat']
    keep_alt = ['row_clusters-flat']
    remove = ['row_clusters-flat']
    remove_alt = ['row_clusters-flat']
    min_val = 1000
    for cell_score, alt_score, alt_sum in cell_max_scores:
        cluster = priorColumnClusters[i]
        cell = column_header[i]
        ref_max, ref_alt, ref_sum = cell_max_score_db[cluster]
        ci = clusters.index(cluster)
        ref_diff = math.pow(2, ref_max - ref_alt) * diff
        ref_alt = math.pow(2, ref_alt)
        cell_diff = math.pow(2, cell_score - alt_score)
        cell_score = math.pow(2, cell_score)
        if cell_diff < min_val:
            min_val = cell_diff
        if cell_diff > ref_diff and cell_diff > diff:
            assignment = 0
            keep.append(cell)
            try:
                keep_alt.append(string.split(cell, ':')[1])
            except Exception:
                keep_alt.append(cell)

        else:
            remove.append(cell)
            try:
                remove_alt.append(string.split(cell, ':')[1])
            except Exception:
                remove_alt.append(cell)

            assignment = 1
        i += 1

    print min_val
    print len(keep), len(remove)
    from import_scripts import sampleIndexSelection
    input_file = heatmap_file
    output_file = heatmap_file[:-4] + '-Singlets.txt'
    try:
        sampleIndexSelection.filterFile(input_file, output_file, keep)
    except:
        sampleIndexSelection.filterFile(input_file, output_file, keep_alt)
    else:
        output_file = heatmap_file[:-4] + '-Multiplets.txt'
        try:
            sampleIndexSelection.filterFile(input_file, output_file, remove)
        except:
            sampleIndexSelection.filterFile(input_file, output_file, remove_alt)


def exportTFcorrelations(filename, TF_file, threshold, anticorrelation=False):
    eo = export.ExportFile(filename[:-4] + '-TF-correlations.txt')
    TFs = simpleListImport(TF_file)
    x, column_header, row_header, dataset_name, group_db = importData(filename)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        D1 = numpy.corrcoef(x)
    i = 0
    correlation_pairs = []
    for score_ls in D1:
        k = 0
        for v in score_ls:
            if str(v) != 'nan':
                if k != i:
                    if row_header[i] in TFs or row_header[k] in TFs:
                        if anticorrelation:
                            if v < -1 * threshold:
                                eo.write(row_header[i] + '\t' + row_header[k] + '\t' + str(v) + '\n')
                        elif v < -1 * threshold or v > threshold:
                            eo.write(row_header[i] + '\t' + row_header[k] + '\t' + str(v) + '\n')
            k += 1

        i += 1

    eo.close()


def TFisoformImport(filename):
    isoform_db = {}
    for line in open(filename, 'rU').xreadlines():
        data = line.rstrip()
        trans, prot, gene, symbol, uid, uid2, uid3 = string.split(data, '\t')
        isoform_db[trans] = (symbol, prot)

    return isoform_db

def exportIntraTFIsoformCorrelations(filename, TF_file, threshold, anticorrelation=False):
    eo = export.ExportFile(filename[:-4] + '-TF-correlations.txt')
    isoform_db = TFisoformImport(TF_file)
    x, column_header, row_header, dataset_name, group_db = importData(filename)
    ### For methylation data or other data with redundant signatures, remove these and only report the first one
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        D1 = numpy.corrcoef(x)
    i = 0
    correlation_pairs = []
    for score_ls in D1:
        k = 0
        for v in score_ls:
            if str(v) != 'nan':
                if k != i:
                    if row_header[i] in isoform_db or row_header[k] in isoform_db:
                        try:
                            gene1, prot1 = isoform_db[row_header[i]]
                            gene2, prot2 = isoform_db[row_header[k]]
                            if gene1 == gene2:
                                if anticorrelation:
                                    if v < -1 * threshold:
                                        eo.write(row_header[i] + '\t' + row_header[k] + '\t' + str(v) + '\n')
                                elif v < -1 * threshold or v > threshold:
                                    eo.write(row_header[i] + '\t' + row_header[k] + '\t' + str(v) + '\n')
                        except:
                            pass
            k += 1
        i += 1
    eo.close()

def PSIfilterAndImpute(folder):
    ### Filter a PSI file and impute missing values based on neighbors
    files = UI.read_directory(folder)
    for file in files:
        filename = folder + '/' + file
        if '.txt' in file:
            eo = export.ExportFile(filename[:-4] + '-impute.txt')
            header = True
            count = 0
            for line in open(filename, 'rU').xreadlines():
                data = cleanUpLine(line)
                values = string.split(data, '\t')
                t0 = values[1]
                tl = values[(-1)]
                vs = values[1:]
                if header:
                    header = False
                    eo.write(line)
                elif len(vs) == len(vs) - vs.count(''):
                    sum_val = sum(map(float, vs)) / len(vs)
                    if sum_val != 1 and sum_val != 0:
                        eo.write(line)
                        count += 1
                elif len(vs) - vs.count('') > len(vs) - 3:
                    new_values = []
                    i = 0
                    for v in vs:
                        if v=='':
                            if i==0: ### if the first element is null
                                try: new_values.append((float(vs[i+1])+float(tl))/2)
                                except: new_values.append(None) ### If two nulls occur in a row
                            elif i==len(vs)-1: ### if the last element is null
                                try: new_values.append((float(vs[i-1])+float(t0))/2)
                                except: new_values.append(None) ### If two nulls occur in a row
                            else: ### if the another element is null
                                try: new_values.append((float(vs[i-1])+float(vs[i+1]))/2)
                                except: new_values.append(None) ### If two nulls occur in a row
                        else:
                            new_values.append(v)
                        i += 1

                    if None not in new_values:
                        sum_val = sum(map(float, new_values)) / len(new_values)
                        if sum_val != 1 and sum_val != 0:
                            eo.write(string.join([values[0]] + map(str, new_values), '\t') + '\n')
                            count += 1
            eo.close()
            print count, '\t', fileg

def summarizePSIresults(folder, TF_file):
    TFs = simpleListImport(TF_file)
    ### Import PSI results and report number of impacted TFs
    files = UI.read_directory(folder)
    eo = export.ExportFile(folder + '/TF_events.txt')
    all_TFs = []
    for file in files:
        TFs_in_file = []
        filename = folder + '/' + file
        if '.txt' in file and 'PSI.' in file:
            header = True
            count = 0
            header = True
            for line in open(filename, 'rU').xreadlines():
                if header:
                    header = False
                else:
                    data = cleanUpLine(line)
                    t = string.split(data, '\t')
                    symbol = string.split(t[0], ':')[0]
                    dPSI = abs(float(t[(-5)]))
                    if symbol in TFs and symbol not in TFs_in_file and dPSI > 0.2:
                        eo.write(string.join(t + [file], '\t') + '\n')
                        TFs_in_file.append(symbol)
                        if symbol not in all_TFs:
                            all_TFs.append(symbol)
                        count += 1
            print file, count, len(all_TFs), string.join(TFs_in_file, ',')

    eo.close()

def convertPSICoordinatesToBED(folder):
    files = UI.read_directory(folder)
    eo = export.ExportFile(folder + '/combined.bed')
    all_TFs = []
    for file in files:
        TFs_in_file = []
        filename = folder + '/' + file
        if '.txt' in file:
            header = True
            count = 0
            header = True
            for line in open(filename, 'rU').xreadlines():
                if header:
                    header = False
                else:
                    data = cleanUpLine(line)
                    t = string.split(data, '\t')
                    symbol = string.split(t[0], ':')[0]
                    try:
                        coordinates = t[7]
                    except:
                        print t
                        sys.exit()
                    else:
                        j1, j2 = string.split(coordinates, '|')
                        c1a, c1b = map(int, string.split(j1.split(':')[1], '-'))
                        strand = '+'
                        if c1a > c1b:
                            c1a, c1b = c1b, c1a
                            strand = '-'
                        c2a, c2b = map(int, string.split(j2.split(':')[1], '-'))
                        if c2a > c2b:
                            c2a, c2b = c2b, c2a
                        chr = string.split(coordinates, ':')[0]
                        uid = string.replace(t[0], ':', '__')
                        eo.write(string.join([chr, str(c1a), str(c1b), uid + '--' + file, strand, str(c1a), str(c1b), '0'], '\t') + '\n')
                        eo.write(string.join([chr, str(c2a), str(c2b), uid + '--' + file, strand, str(c2a), str(c2b), '0'], '\t') + '\n')
    eo.close()

def convertPSIConservedCoordinatesToBED(Mm_Ba_coordinates, Ba_events):
    if 'Baboon' in Mm_Ba_coordinates:
            equivalencies={'Heart':['Heart'],
                   'Kidney':['Kidney-cortex','Kidney-medulla'],
                   'WFAT':['White-adipose-pericardial','White-adipose-mesenteric','White-adipose-subcutaneous','Omental-fat'],
                   'BFAT':['White-adipose-pericardial','White-adipose-mesenteric','White-adipose-subcutaneous','Omental-fat'],
                   'Lung':['Lungs'],
                   'Cere':['Cerebellum','Ventromedial-hypothalamus','Habenula','Pons','Pineal-gland','Visual-cortex','Lateral-globus-pallidus',
                           'Paraventricular-nuclei','Arcuate-nucleus','Suprachiasmatic-nuclei','Putamen','Optic-nerve-head', 'Medial-globus-pallidus',
                           'Amygdala','Prefontal-cortex','Dorsomedial-hypothalamus'],
                   'BS':['Cerebellum','Ventromedial-hypothalamus','Habenula','Pons','Pineal-gland','Visual-cortex','Lateral-globus-pallidus',
                           'Paraventricular-nuclei','Arcuate-nucleus','Suprachiasmatic-nuclei','Putamen','Optic-nerve-head', 'Medial-globus-pallidus',
                           'Amygdala','Prefontal-cortex','Dorsomedial-hypothalamus'],
                   'Hypo':['Cerebellum','Ventromedial-hypothalamus','Habenula','Pons','Pineal-gland','Visual-cortex','Lateral-globus-pallidus',
                           'Paraventricular-nuclei','Arcuate-nucleus','Suprachiasmatic-nuclei','Putamen','Optic-nerve-head', 'Medial-globus-pallidus',
                           'Amygdala','Prefontal-cortex','Dorsomedial-hypothalamus'],
                   'Adrenal':['Adrenal-cortex','Adrenal-medulla'],
                   'SM':['Muscle-gastrocnemian','Muscle-abdominal'],
                   'Liver':['Liver'],
                    }
    else:
            equivalencies={'Heart':['Heart'],
                   'Kidney':['Kidney','Kidney'],
                   'WFAT':['WFAT'],
                   'BFAT':['BFAT'],
                   'Lung':['Lungs'],
                   'Adrenal':['Adrenal'],
                   'Liver':['Liver'],
                    }
    eo = export.ExportFile(Mm_Ba_coordinates[:-4] + '-matched.txt')
    eo2 = export.ExportFile(Mm_Ba_coordinates[:-4] + '-matrix.txt')
    mouse_events = {}
    baboon_events = {}
    baboon_corridinates = {}
    ### This mouse circadian events file has been lifted over to baboon coordinates
    countX = 0
    for line in open(Mm_Ba_coordinates, 'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data, '\t')
        chr, c1, c2, event, strand, null, null, null = values
        event = string.replace(event, '__', ':')
        event, tissue = event.split('--')
        junctions = string.split(event, ':')[1:]
        junctions = string.join(junctions, ':')
        junctions = string.split(junctions, '|')
        junctions.sort() ### make a unique event
        junctions = string.join(junctions, '|')
        symbol = string.split(event, ':')[0]
        event = symbol + ':' + junctions
        countX += 1
        tissue = string.replace(tissue, '_event_annot_file.txt', '')
        tissue = string.replace(tissue, 'PSI.', '')
        tissue = string.replace(tissue, '_Mm', '')
        junction = chr + ':' + c2 + '-' + c1
        alt_junction1 = chr + ':' + str(int(c2) + 1) + '-' + str(int(c1) + 1)
        alt_junction2 = chr + ':' + str(int(c2) - 1) + '-' + str(int(c1) - 1)
        try:
            mouse_events[junction].append([event, tissue])
        except:
            mouse_events[junction] = [[event, tissue]]
        else:
            try:
                mouse_events[alt_junction1].append([event, tissue])
            except:
                mouse_events[alt_junction1] = [[event, tissue]]
            else:
                try:
                    mouse_events[alt_junction2].append([event, tissue])
                except:
                    mouse_events[alt_junction2] = [[event, tissue]]
                else:
                    junction = chr + ':' + c1 + '-' + c2
                    alt_junction1 = chr + ':' + str(int(c1) + 1) + '-' + str(int(c2) + 1)
                    alt_junction2 = chr + ':' + str(int(c1) - 1) + '-' + str(int(c2) - 1)
                    try:
                        mouse_events[junction].append([event, tissue])
                    except:
                        mouse_events[junction] = [[event, tissue]]

                try:
                    mouse_events[alt_junction1].append([event, tissue])
                except:
                    mouse_events[alt_junction1] = [[event, tissue]]

            try:
                mouse_events[alt_junction2].append([event, tissue])
            except:
                mouse_events[alt_junction2] = [[event, tissue]]

    for line in open(Ba_events, 'rU').xreadlines():
        data = cleanUpLine(line)
        values = string.split(data, '\t')
        event, tissue_num, tissues, coordinates = values
        junctions = string.split(event, ':')[1:]
        junctions = string.join(junctions, ':')
        junctions = string.split(junctions, '|')
        junctions.sort()
        junctions = string.join(junctions, '|')
        symbol = string.split(event, ':')[0]
        event = symbol + ':' + junctions
        baboon_corridinates[event] = coordinates
        try:
            j1, j2 = string.split(coordinates, '|')
        except:
            continue
        else:
            tissues = tissues.split('|')
            try:
                baboon_events[j1].append([event, tissues])
            except:
                baboon_events[j1] = [[event, tissues]]

            try:
                baboon_events[j2].append([event, tissues])
            except:
                baboon_events[j2] = [[event, tissues]]

    print len(mouse_events), len(baboon_events)
    common = 0
    matched_events = {}
    matched_mm_events = {}
    tissue_matrix = {}
    mm_single_tissue_counts = {}
    ba_single_tissue_counts = {}
    for junction in mouse_events:
        if junction in baboon_events:
            common += 1
            mm_events = {}
            for mm_event, mm_tissue in mouse_events[junction]:
                try:
                    mm_events[mm_event].append(mm_tissue)
                except:
                    mm_events[mm_event] = [mm_tissue]

            for mm_event in mm_events:
                mm_tissues = mm_events[mm_event]
                mm_tissues = unique.unique(mm_tissues)
                for ba_event, ba_tissues in baboon_events[junction]:
                    ba_tissues = unique.unique(ba_tissues)
                    matched_events[(mm_event, ba_event)] = (mm_tissues, ba_tissues)
                    matched_mm_events[mm_event] = []

    def matchingTissues(mouse, baboon):
        m_matches = []
        b_matches = []
        for m in mouse:
            for b in baboon:
                if m in equivalencies:
                    if b in equivalencies[m]:
                        m_matches.append(m)
                        b_matches.append(b)

        if len(m_matches) == 0:
            return ''
        m_matches = string.join(unique.unique(m_matches), ', ')
        b_matches = string.join(unique.unique(b_matches), ', ')
        return m_matches + ':' + b_matches

    for mm_event, ba_event in matched_events:
        mm_tissues, ba_tissues = matched_events[(mm_event, ba_event)]
        matching_tissues = matchingTissues(mm_tissues, ba_tissues)
        eo.write(string.join([mm_event, ba_event, string.join(mm_tissues, '|'), string.join(ba_tissues, '|'), str(len(mm_tissues)), str(len(ba_tissues)), matching_tissues], '\t') + '\n')
        for mt in mm_tissues:
            for bt in ba_tissues:
                try:
                    tissue_matrix[(mt, bt)] += 1
                except:
                    tissue_matrix[(mt, bt)] = 1
                else:
                    try:
                        mm_single_tissue_counts[mt] += 1
                    except:
                        mm_single_tissue_counts[mt] = 1

                    try:
                        ba_single_tissue_counts[bt] += 1
                    except:
                        ba_single_tissue_counts[bt] = 1

    print mm_single_tissue_counts['Heart']
    print tissue_matrix[('Heart', 'Heart')]
    tissue_matrix_table = []
    ba_tissues = ['Tissues']
    for bt in ba_single_tissue_counts:
        ba_tissues.append(bt)

    eo2.write(string.join(ba_tissues, '\t') + '\n')
    for mt in mm_single_tissue_counts:
        table = []
        for bt in ba_single_tissue_counts:
            if bt == 'Thyroid' and mt == 'Heart':
                print tissue_matrix[(mt, bt)]
                print tissue_matrix[(mt, bt)] / (1.0 * ba_single_tissue_counts[bt])
            try:
                table.append(str(tissue_matrix[(mt, bt)] / (1.0 * ba_single_tissue_counts[bt])))
            except:
                table.append('0')

        eo2.write(string.join([mt] + table, '\t') + '\n')

    print common, len(matched_events), len(matched_mm_events)
    eo.close()
    eo2.close()

def rankExpressionRescueFromCellHarmony(organized_diff_ref, repair1_folds, repair2_folds, reference_fold_dir, repair_dir1, repair_dir2):

    def importCellHarmonyDEGs(folder, repair=False):
        print folder
        files = os.listdir(folder)
        DEG_db = {}
        for file in files:
            filename = folder + '/' + file
            if '.txt' in file and 'GE.' in file:
                header = True
                count = 0
                header = True
                file = file[:-4]
                file = string.split(file[3:], '_')[0]
                if file[:2] == 'DM':
                    file = 'global'
                for line in open(filename, 'rU').xreadlines():
                    if header:
                        header = False
                    else:
                        data = cleanUpLine(line)
                        t = string.split(data, '\t')
                        GeneID, SystemCode, LogFold, rawp, adjp, Symbol, avg_g2, avg_g1 = t
                        rawp = float(rawp)
                        adjp = float(adjp)
                        if float(LogFold) > 0:
                            direction = 'positive'
                        else:
                            direction = 'negative'
                        if repair:
                            if float(LogFold) > 0:
                                fold = math.pow(2, float(LogFold))
                            else:
                                fold = -1 / math.pow(2, float(LogFold))
                            if Symbol == 'BC049762':
                                print 'BC049762', file, LogFold, fold
                            if abs(fold) > 1.5 and adjp < 0.05:
                                try:
                                    DEG_db[Symbol].append([file, direction])
                                except:
                                    DEG_db[Symbol] = [[file, direction]]

                        else:
                            try:
                                DEG_db[Symbol].append([file, direction])
                            except:
                                DEG_db[Symbol] = [[file, direction]]

        return DEG_db

    ref_DEGs = importCellHarmonyDEGs(reference_fold_dir)
    repaired_DEGs = importCellHarmonyDEGs(repair_dir1, repair=True)
    repaired2_DEGs = importCellHarmonyDEGs(repair_dir2, repair=True)

    def importCellHarmonyPseudoBulkFolds(filename):
        fold_db = {}
        header = True
        for line in open(filename, 'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data, '\t')
            if header:
                fold_db['header'] = t[1:]
                header = False
            else:
                uid = t[0]
                folds = t[1:]
                fold_db[uid] = folds

        return fold_db

    repaired_fold_db = importCellHarmonyPseudoBulkFolds(repair1_folds)
    repaired2_fold_db = importCellHarmonyPseudoBulkFolds(repair2_folds)
    import collections
    ordered_ref_degs = collections.OrderedDict()
    ordered_cluster_genes = collections.OrderedDict()
    repair_verified = collections.OrderedDict()
    repair2_verified = collections.OrderedDict()
    cluster_ordered_ref_db = collections.OrderedDict()
    header = True
    for line in open(organized_diff_ref, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if header:
            ref_header = t
            header = False
        else:
            cluster, geneID = string.split(t[0], ':')
            cluster = string.split(cluster, '_')[0]
            if cluster[:2] == 'DM':
                cluster = 'global'
            ordered_ref_degs[(geneID, cluster)] = t
            try:
                cluster_ordered_ref_db[cluster].append(geneID)
            except:
                cluster_ordered_ref_db[cluster] = [geneID]

    repaired_verified = {}
    verified = {}
    for geneID, ref_cluster in ordered_ref_degs:
        for cluster, ref_direction in ref_DEGs[geneID]:
            if geneID in repaired_DEGs:
                for repair_cluster, repair_direction in repaired_DEGs[geneID]:
                    if repair_cluster == cluster and ref_direction != repair_direction and ('Neu' in repair_cluster or 'global' in repair_cluster):
                        try:
                            repair_verified[repair_cluster].append(geneID)
                        except:
                            repair_verified[repair_cluster] = [geneID]
                        else:
                            print geneID + '\t' + repair_direction + '\t' + repair_cluster + '\tR412X-HMZ'
                            try:
                                verified[geneID].append('R412X-HMZ')
                            except:
                                verified[geneID] = ['R412X-HMZ']

            if geneID in repaired2_DEGs:
                for repair_cluster, repair_direction in repaired2_DEGs[geneID]:
                    if repair_cluster == cluster and ref_direction != repair_direction and ('Neu' in cluster or 'global' in cluster):
                        try:
                            repair2_verified[repair_cluster].append(geneID)
                        except:
                            repair2_verified[repair_cluster] = [geneID]
                        else:
                            print geneID + '\t' + repair_direction + '\t' + repair_cluster + '\t' + 'R412X-Irf8'
                            try:
                                verified[geneID].append('R412X-Irf8')
                            except:
                                verified[geneID] = ['R412X-Irf8']

    for gene in verified:
        verified[gene] = unique.unique(verified[gene])

    eo1 = export.ExportFile(organized_diff_ref[:-4] + '-Repair-Sorted.txt')
    eo2 = export.ExportFile(organized_diff_ref[:-4] + '-Repaired-Only.txt')
    header = ref_header + repaired_fold_db['header'] + repaired2_fold_db['header']
    eo1.write(string.join(header, '\t') + '\n')
    eo2.write(string.join(header, '\t') + '\n')
    print len(ordered_ref_degs)
    print len(repaired_fold_db)
    print len(repaired2_fold_db)
    print len(repair_verified)
    print len(repair2_verified)
    print len(verified)
    print len(ordered_ref_degs)
    prior_cluster = None
    added_genes = []
    for geneID, cluster in ordered_ref_degs:
        try:
            folds = ordered_ref_degs[(geneID, cluster)] + repaired_fold_db[geneID] + repaired2_fold_db[geneID]
        except:
            print '...Error in identifying match UID for:', geneID
            added_genes.append(geneID)
            continue
        else:
            if geneID not in verified:
                eo1.write(string.join(folds, '\t') + '\n')
            elif len(verified[geneID]) > 1:
                added_genes.append(geneID)
            elif 'R412X-HMZ' in verified[geneID]:
                added_genes.append(geneID)
            else:
                eo2.write(string.join(folds, '\t') + '\n')
                added_genes.append(geneID)

    eo1.close()
    eo2.close()
    return


if __name__ == '__main__':
    organized_diff_ref = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-4-Gfi1-R412X-ModGMP-1694-ADT/outs/filtered_gene_bc_matrices/Merged-Cells/centroid-revised/custom/cellHarmony/OrganizedDifferentials.txt'
    repair1_folds = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-5-Gfi1-R412X-R412X-ModGMP-1362-ADT/outs/filtered_gene_bc_matrices/Merged-Cells/hybrid/cellHarmony-vs-DM2-1.2-fold-adjp/OtherFiles/exp.ICGS-cellHarmony-reference__DM-5-Gfi1-R412X-R412X-ModGMP-1362-D7Cells-ADT-Merged_matrix_CPTT-AllCells-folds.txt'
    repair2_folds = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-6-Gfi1-R412X-Irf8-ModGMP-1499-ADT/outs/filtered_gene_bc_matrices/Merged-Cells-iseq/cellHarmony-centroid-revsied/hybrid/cellHarmony/OtherFiles/exp.ICGS-cellHarmony-reference__DM-6-Gfi1-R412X-Irf8-ModGMP-1499-ADT_matrix-3_matrix_CPTT-hybrid-AllCells-folds.txt'
    reference_fold_dir = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-4-Gfi1-R412X-ModGMP-1694-ADT/outs/filtered_gene_bc_matrices/Merged-Cells/centroid-revised/custom/cellHarmony/DifferentialExpression_Fold_1.2_adjp_0.05'
    repair_dir1 = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-5-Gfi1-R412X-R412X-ModGMP-1362-ADT/outs/filtered_gene_bc_matrices/Merged-Cells/hybrid/vs-R412X-het/cellHarmony/DifferentialExpression_Fold_1.2_adjp_0.05'
    repair_dir1 = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-5-Gfi1-R412X-R412X-ModGMP-1362-ADT/outs/filtered_gene_bc_matrices/Merged-Cells/hybrid/vs-R412X-het/cellHarmony/OtherFiles/DEGs-LogFold_0.0_rawp'
    repair_dir2 = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-6-Gfi1-R412X-Irf8-ModGMP-1499-ADT/outs/filtered_gene_bc_matrices/Merged-Cells-iseq/cellHarmony-centroid-revsied/hybrid/vs-R412X-Het/cellHarmony/DifferentialExpression_Fold_1.2_adjp_0.05'
    repair_dir2 = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/WuXi-David-Nature-Revision/PROJ-00584/fastqs/DM-6-Gfi1-R412X-Irf8-ModGMP-1499-ADT/outs/filtered_gene_bc_matrices/Merged-Cells-iseq/cellHarmony-centroid-revsied/hybrid/vs-R412X-Het/cellHarmony/OtherFiles/DEGs-LogFold_0.0_rawp'
    TF_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/NCI-R01/CCSB_TFIso_Clones.txt'
    PSI_dir = '/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/AltResults/AlternativeOutput/OncoSPlice-All-Samples-filtered-names/SubtypeAnalyses-Results/round1/Events-dPSI_0.1_adjp/'
    simpleCombineFiles('/Volumes/salomonis2/NCI-R01/Harvard/BRC_PacBio_Seq/metadataanalysis/PSICluster/TCGA/FilteredTF')
    sys.exit()
    filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Anukana/Breast-Cancer/TF-isoform/TF_ratio_correlation-analysis/tcga_rsem_isopct_filtered-filtered.2-filtered.txt'
    TF_file = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Anukana/Breast-Cancer/TF-isoform/Ensembl-isoform-key-CCSB.txt'
    input_file = '/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana/UO1analysis/xenabrowserFiles/tcga_rsem_isoform_tpm_filtered.txt'
    Mm_Ba_coordinates = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Krithika/Baboon-Mouse/mm10-circadian_liftOverTo_baboon.txt'
    Ba_events = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Krithika/Baboon-Mouse/Baboon_metacycle-significant-AS-coordinates.txt'
    Mm_Ba_coordinates = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Krithika/Human-Mouse/hg19-mm10-12-tissue-circadian.txt'
    Ba_events = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Krithika/Human-Mouse/Human_CYCLOPS-significant-AS-coordinates.txt'
    filename = '/Users/saljh8/Desktop/DemoData/Venetoclax/D4/cellHarmony-rawp-stringent/gene_summary.txt'
    filename = '/Volumes/salomonis2/LabFiles/Nathan/10x-PBMC-CD34+/AML-p27-pre-post/pre/cellHarmony-latest/gene_summary-p27.txt'
    filename = '/Volumes/salomonis2/LabFiles/Dan-Schnell/To_cellHarmony/MIToSham/Input/cellHarmony/cell-frequency-stats.txt'
    index1 = 2
    index2 = 3
    x_axis = 'Number of Differentially Expressed Genes'
    y_axis = 'Comparisons'
    title = 'Hippocampus - Number of Differentially Expressed Genes'
    index1 = 2
    index2 = 3
    x_axis = 'Number of DEGs'
    y_axis = 'Reference clusters'
    title = 'cellHarmony Differentially Expressed Genes'
    index1 = -2
    index2 = -1
    x_axis = 'Cell-State Percentage'
    y_axis = 'Reference clusters'
    title = 'Assigned Cell Frequencies'
    diff = 0.7
    print 'diff:', diff
    a = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/H202SC19091002-Davis-US/Merged-Human-data/MergedFiles-h34_38_T-h34_TA/ICGS-NMF_cosine/ADT-RNA/exp.CD34+CD38-combined-transposed.txt'
    b = '/Volumes/salomonis2/Immune-10x-data-Human-Atlas/Bone-Marrow/Stuart/Browser/ExpressionInput/HS-compatible_symbols.txt'
    transposeMatrix(a)
    sys.exit()
    b = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/Gly-CCO-Seq-Hs/10X-Grimes-Gly-CCO-Seq-20190925-3v3hg/10X-Grimes-Gly-CCO-Seq/outs/filtered_feature_bc_matrix/ICGS-NMF_cosine_cc/only1gRNA-per-cell.txt'
    a = '/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/H202SC19091002-Davis-US/Merged-Human-data/MergedFiles-h34_38_T-h34_TA/ICGS-NMF_cosine/ADT-RNA/exp.CD34+CD38-x.txt'
    convertGroupsToBinaryMatrix(b, b, cellHarmony=False)
    sys.exit()
    a = '/Users/saljh8/Desktop/temp/groups.TNBC.txt'
    b = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/tests/clusters.txt'
    a = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/SpliceICGS.R1.Depleted.12.27.17/all-depleted-and-KD'
    query_dataset = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/exp.GSE81682_HTSeq-cellHarmony-filtered.txt'
    all_tpm = '/Users/saljh8/Desktop/demo/BoneMarrow/ExpressionInput/exp.BoneMarrow-scRNASeq.txt'
    markerfinder = '/Users/saljh8/Desktop/demo/BoneMarrow/ExpressionOutput1/MarkerFinder/AllGenes_correlations-ReplicateBasedOriginal.txt'
    signature_genes = '/Users/saljh8/Desktop/Grimes/KashishNormalization/test/Panorama.txt'
    state = 'Multi-Lin'
    query_dataset = None
    all_tpm = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/MultiLin/Gottgens_HarmonizeReference.txt'
    signature_genes = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/MultiLin/Gottgens_HarmonizeReference.txt'
    markerfinder = '/Users/saljh8/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    state = 'Eryth_Multi-Lin'
    filename = '/Users/saljh8/Desktop/Grimes/GEC14078/MergedFiles.txt'
    filename = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart72/Mm/junction1/junction_critical-junction-seq.txt'
    filename = '/Users/saljh8/Downloads/CoexpressionAtlas.txt'
    folder = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart72/ensembl/Hs'
    try:
        files = UI.read_directory(folder)
        for file in files:
            if '.bed' in file:
                pass

    except Exception:
        pass

    countinp = '/Volumes/salomonis2/SinghLab/20150715_single_GCBCell/bams/ExpressionInput/counts.Bcells.txt'
    IGH_gene_file = '/Volumes/salomonis2/SinghLab/20150715_single_GCBCell/bams/ExpressionInput/IGH_genes.txt'
    import UI
    filename = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/genes.tpm_tracking-ordered.txt'
    TFs = '/Users/saljh8/Desktop/Grimes/KashishNormalization/3-25-2015/TF-by-gene_matrix/all-TFs2.txt'
    folder = '/Users/saljh8/Downloads/BLASTX2_Gecko.tab'
    genes = ['Gfi1', 'Irf8']
    gene_list = [
     'S100a8', 'Chd7', 'Ets1', 'Chd7', 'S100a8']
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
    gene_list_file = '/Volumes/salomonis2/PublicDatasets/GSE128423_RAW-Stroma/merged/Healthy-only/ExpressionInput/genes.txt'
    gene_list_file = '/Volumes/salomonis2-3/CCHMC-Collaborations/Theodosia-Kalfa/Combined-10X-CPTT/ExpressionInput/genes.txt'
    genesets = importGeneList(gene_list_file, n=35)
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
    filename = '/Volumes/salomonis2/PublicDatasets/GSE128423_RAW-Stroma/merged/Healthy-only/ExpressionInput/exp.BoneMarrow-10k.txt'
    filename = '/Volumes/salomonis2-3/CCHMC-Collaborations/Theodosia-Kalfa/Combined-10X-CPTT/ExpressionInput/exp.MergedFiles-ICGS.txt'
    print genesets
    for gene_list in genesets:
        multipleSubPlots(filename, gene_list, SubPlotType='column', n=35)

    sys.exit()
    plotHistogram(filename)
    sys.exit()
    filename = '/Users/saljh8/Desktop/Grimes/Expression_final_files/ExpressionInput/amplify-wt/DataPlots/Clustering-exp.myeloid-steady-state-PCA-all_wt_myeloid_SingleCell-Klhl7 Dusp7 Slc25a33 H6pd Bcorl1 Sdpr Ypel3 251000-hierarchical_cosine_cosine.cdt'
    openTreeView(filename)
    sys.exit()
    pdf1 = '/Users/saljh8/Desktop/Grimes/1.pdf'
    pdf2 = '/Users/saljh8/Desktop/Grimes/2.pdf'
    outPdf = '/Users/saljh8/Desktop/Grimes/3.pdf'
    merge_horizontal(outPdf, pdf1, pdf2)
    sys.exit()
    mergePDFs(pdf1, pdf2, outPdf)
    sys.exit()
    filename = '/Volumes/SEQ-DATA/CardiacRNASeq/BedFiles/ExpressionOutput/Clustering/SampleLogFolds-CardiacRNASeq.txt'
    ica(filename)
    sys.exit()
    features = 5
    matrix, column_header, row_header, dataset_name, group_db = importData(filename)
    Kmeans(features, column_header, row_header)
    sys.exit()
    filename = '/Users/saljh8/Desktop/delete.txt'
    filenames = [
     filename]
    outputClusters(filenames, [])
    sys.exit()
    pruned_folder = '/Users/nsalomonis/Desktop/CBD/LogTransformed/GO-Elite/GO-Elite_results/CompleteResults/ORA_pruned/'
    input_ora_folder = '/Users/nsalomonis/Desktop/CBD/LogTransformed/GO-Elite/input/'
    files = UI.read_directory(pruned_folder)
    for file in files:
        if '.sif' in file:
            input_file = string.join(string.split(file, '-')[:-1], '-') + '.txt'
            sif_file = pruned_folder + file
            input_file = input_ora_folder + input_file
            buildGraphFromSIF('Ensembl', 'Hs', sif_file, input_file)

    sys.exit()
    filenames = [filename]
    outputClusters(filenames, [])
