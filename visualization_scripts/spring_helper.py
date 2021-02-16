import os
import numpy as np
import scipy
import scipy.stats
import sklearn.cluster
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csc_matrix
import scipy.io
import pickle
from scipy.spatial.distance import pdist
from datetime import datetime
import json
import time
import matplotlib.pyplot as plt
import sys

########## USEFUL SPARSE FUNCTIONS

def sparse_var(E, axis=0):
    ''' variance across the specified axis '''
    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_multiply(E, a):
    ''' multiply each row of E by a scalar '''
    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def sparse_zscore(E):
    ''' z-score normalize each column of E '''
    mean_gene = E.mean(0)
    stdev_gene = np.sqrt(sparse_var(E))
    return sparse_multiply((E - mean_gene).T, 1/stdev_gene).T

def average_profile(E, all_genes, gene_set):
    gene_set = [g.upper() for g in gene_set]
    gix = np.array([i for i,g in enumerate(all_genes) if g.upper() in gene_set], dtype=int)
    if len(gix) == 0: 
        return np.zeros(E.shape[1])
    else:
        return sparse_zscore(E[:,gix]).mean(1).A.squeeze()

######### LOADING DATA

def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []
    gene_dict = {}

    with open(filename) as f:
        for iL in range(skip_rows):
            f.readline()
        for l in f:
            gene = l.strip('\n').split(delimiter)[column]
            gene = gene.replace("/", "\/")
            if gene in gene_dict:
                gene_dict[gene] += 1
                gene_list.append(gene + '__' + str(gene_dict[gene]))
                if gene_dict[gene] == 2:
                    i = gene_list.index(gene)
                    gene_list[i] = gene + '__1'
            else: 
                gene_dict[gene] = 1
                gene_list.append(gene)
    return gene_list


def load_pickle(fname):
    '''
    Load .pickle(.gz) data
    '''
    if fname.endswith('.gz'):
        tmpsuffix = str(np.random.randint(1e9))
        os.system('gunzip -c "' + fname + '" > tmp' + tmpsuffix)
        dat = pickle.load(open('tmp' + tmpsuffix, 'rb'))
        os.system('rm tmp' + tmpsuffix)
    else:
        dat = pickle.load(open(fname, 'rb'))
    return dat


### loading counts

def file_opener(filename):
    fileData = open(filename)
    if filename.endswith('.gz'):
        import gzip
        outData = gzip.GzipFile(fileobj = fileData, mode = 'rb')
    elif filename.endswith('.zip'):
        import zipfile
        zipData = zipfile.ZipFile(fileData, 'r')
        fnClean = filename.strip('/').split('/')[-1][:-4]
        outData = zipData.open(fnClean)
    else:
        outData = fileData
    return outData


def load_mtx(file_data):
    ''' Reads mtx file or, supposedly, an open file object
        Returns scipy.sparse.coo_matrix (if sparse)'''
    return scipy.io.mmread(file_data).tocsc()

def load_npz(file_data):
    return scipy.sparse.load_npz(file_data).tocsc()

def load_npy(file_data):
    return scipy.sparse.csc_matrix(np.load(file_data))

def load_text(file_data,delim='\t', load_cell_bcs=False):
    X_data = []
    X_row = []
    X_col = []
    cell_bcs = []

    start_column = -1
    start_row = -1
    for row_ix, dat in enumerate(file_data):
        dat = dat.strip('\n').split(delim)
        if start_row == -1:
            current_col = 0
            found_float = False
            while not found_float and current_col < len(dat):
                try: 
                    tmp = float(dat[current_col])
                    
                    try:
                        rowdat = np.array(map(float, dat[current_col:]))
                        ncol = len(rowdat)
                        col_ix = np.nonzero(rowdat)[0]

                        found_float = True
                        start_row = row_ix
                        start_column = current_col

                        X_col.extend(col_ix)
                        X_row.extend([row_ix - start_row] * len(col_ix))
                        X_data.extend(rowdat[col_ix])
                        if load_cell_bcs: cell_bcs.append(dat[0])
                    except:
                        current_col += 1

                except:
                    current_col += 1
        else:
            try:
                if load_cell_bcs: cell_bcs.append(dat[0])
                rowdat = np.array(map(float, dat[start_column:]))
                if len(rowdat) != ncol:
                    return 'ERROR: Rows have different numbers of numeric columns.'
                col_ix = np.nonzero(rowdat)[0]
                X_col.extend(col_ix)
                X_row.extend([row_ix - start_row] * len(col_ix))
                X_data.extend(rowdat[col_ix])
            except:
                return 'ERROR: Rows have different numbers of numeric columns.'

    if start_row == -1:
        return 'ERROR: no numeric values found'
    nrow = row_ix - start_row + 1
    E = scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=float, shape=(nrow, ncol)).tocsc()
    if load_cell_bcs: return E, np.array(cell_bcs)
    else: return E

def text_to_sparse(file_data,delim='\t',start_row=0,start_column=0,data_type=float):
    output = [[]]

    X_data = []
    X_row = []
    X_col = []

    for row_ix, dat in enumerate(file_data):
        dat = dat.strip('\n').split(delim)
        if row_ix >= start_row:
            rowdat = np.array(map(data_type, dat[start_column:]))
            col_ix = np.nonzero(rowdat)[0]
            X_col.extend(col_ix)
            X_row.extend([row_ix - start_row] * len(col_ix))
            X_data.extend(rowdat[col_ix])
    
    ncol = len(rowdat)
    nrow = row_ix - start_row + 1
    
    E = scipy.sparse.coo_matrix((X_data, (X_row, X_col)), dtype=data_type, shape=(nrow, ncol))
    
    return E

########## CELL FILTERING
def filter_dict(d, filt):
    for k,v in d.items():
        if k != 'meta':
            if len(v.shape) == 1:
                d[k] = v[filt]
            else:
                d[k] = v[filt,:]
    return d


########## GENE FILTERING

def runningquantile(x, y, p, nBins):

    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0]+dx/2, x[-1]-dx/2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i-1]
            else:
                yOut[i] = np.nan

    return xOut, yOut


def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    '''

    ncell = E.shape[0]

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:,gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
    b = b[:-1] + np.diff(b)/2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))
    errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func = errFun, x0=[b0], disp=False)
    a = c / (1+b) - 1


    v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
    CV_eff = np.sqrt((1+a)*(1+b) - 1);
    CV_input = np.sqrt(b);

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b

def filter_genes(E, base_ix = [], min_vscore_pctl = 85, min_counts = 3, min_cells = 3, show_vscore_plot = False, sample_name = ''):
    ''' 
    Filter genes by expression level and variability
    Return list of filtered gene indices
    '''

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(E[base_ix, :])
    ix2 = Vscores>0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = (((E[:,gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (Vscores >= min_vscore))
    
    if show_vscore_plot:
        import matplotlib.pyplot as plt
        x_min = 0.5*np.min(mu_gene)
        x_max = 2*np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max/x_min)*np.linspace(0,1,100))
        yTh = (1 + a)*(1+b) + b * xTh
        plt.figure(figsize=(8, 6));
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene), c = [.8,.8,.8], alpha = 0.3, edgecolors='');
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[ix], c = [0,0,0], alpha = 0.3, edgecolors='');
        plt.plot(np.log10(xTh),np.log10(yTh));
        plt.title(sample_name)
        plt.xlabel('log10(mean)');
        plt.ylabel('log10(Fano factor)');
        plt.show()

    return gene_ix[ix]


def remove_corr_genes(E, gene_list, exclude_corr_genes_list, test_gene_idx, min_corr = 0.1):
    seed_ix_list = []
    for l in exclude_corr_genes_list:
        seed_ix_list.append(np.array([i for i in range(len(gene_list)) if gene_list[i] in l], dtype=int))

    exclude_ix = []
    for iSet in range(len(seed_ix_list)):
        seed_ix = seed_ix_list[iSet][E[:,seed_ix_list[iSet]].sum(axis=0).A.squeeze() > 0]

        tmp = sparse_zscore(E[:,seed_ix])
        tmp = tmp.sum(1).A.squeeze()

        c = np.zeros(len(test_gene_idx))
        for iG in range(len(c)):
            c[iG],_ = scipy.stats.pearsonr(tmp, E[:,test_gene_idx[iG]].A.squeeze())

        exclude_ix.extend([test_gene_idx[i] for i in range(len(test_gene_idx)) if (c[i]) >= min_corr])
        #print len(exclude_ix)
    exclude_ix = np.array(exclude_ix)
    #print np.array(gene_list)[exclude_ix]

    return np.array([g for g in test_gene_idx if g not in exclude_ix], dtype=int)


########## CELL NORMALIZATION

def tot_counts_norm(E, exclude_dominant_frac = 1, included = [], target_mean = 0):
    ''' 
    Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
    Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts 
    '''

    E = E.tocsc()
    ncell = E.shape[0]
    if len(included) == 0:
        if exclude_dominant_frac == 1:
            tots_use = E.sum(axis=1)
        else:
            tots = E.sum(axis=1)
            wtmp = scipy.sparse.lil_matrix((ncell, ncell))
            wtmp.setdiag(1. / tots)
            included = np.asarray(~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0))[0,:]
            tots_use = E[:,included].sum(axis = 1)
            #print 'Excluded %i genes from normalization' %(np.sum(~included))
    else:
        tots_use = E[:,included].sum(axis = 1)

    if target_mean == 0:
        target_mean = np.mean(tots_use)

    w = scipy.sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_mean) / tots_use)
    Enorm = w * E

    return Enorm.tocsc(), target_mean, included

########## DIMENSIONALITY REDUCTION

def get_pca(E, base_ix=[], numpc=50, keep_sparse=False, normalize=True):
    '''
    Run PCA on the counts matrix E, gene-level normalizing if desired
    Return PCA coordinates
    '''
    # If keep_sparse is True, gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    if keep_sparse:
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(n_components=numpc)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc)

    pca.fit(Z[base_ix,:])
    return pca.transform(Z)


def preprocess_and_pca(E, total_counts_normalize=True, norm_exclude_abundant_gene_frac=1, min_counts=3, min_cells=5, min_vscore_pctl=85, gene_filter=None, num_pc=50, sparse_pca=False, show_vscore_plot=False):
    '''
    Total counts normalize, filter genes, run PCA
    Return PCA coordinates and filtered gene indices
    '''

    if total_counts_normalize:
        #print 'Total count normalizing'
        E = tot_counts_norm(E, exclude_dominant_frac = norm_exclude_abundant_gene_frac)[0]

    if gene_filter is None:
        #print 'Finding highly variable genes'
        gene_filter = filter_genes(E, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts, min_cells=min_cells, show_vscore_plot=show_vscore_plot)

    #print 'Using %i genes for PCA' %len(gene_filter)
    PCdat = get_pca(E[:,gene_filter], numpc=num_pc, keep_sparse=sparse_pca)

    return PCdat, gene_filter



########## GRAPH CONSTRUCTION

def get_knn_graph(X, k=5, dist_metric='euclidean', approx=False, return_edges=True):
    '''
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    '''

    t0 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except:
            approx = False
            #print 'Could not find library "annoy" for approx. nearest neighbor search'
    if approx:
        #print 'Using approximate nearest neighbor search'

        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)

        for i in xrange(ncell):
            annoy_index.add_item(i, list(X[i,:]))
        annoy_index.build(10) # 10 trees

        knn = []
        for iCell in xrange(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

    else:
        #print 'Using sklearn NearestNeighbors'

        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric, algorithm='brute').fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set([])
        for i in range(knn.shape[0]):
            for j in knn[i,:]:
                links.add(tuple(sorted((i,j))))

        t_elapse = time.time() - t0
        #print 'kNN graph built in %.3f sec' %(t_elapse)

        return links, knn
    return knn

def build_adj_mat(edges, n_nodes):
    A = scipy.sparse.lil_matrix((n_nodes, n_nodes))
    for e in edges:
        i, j = e
        A[i,j] = 1
        A[j,i] = 1
    return A.tocsc()


########## CLUSTERING

def get_spectral_clusters(A, k):
    from sklearn.cluster import SpectralClustering
    spec = SpectralClustering(n_clusters=k, random_state = 0, affinity = 'precomputed', assign_labels = 'discretize')
    return spec.fit_predict(A)


def get_louvain_clusters(nodes, edges):
    import networkx as nx
    import community
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    return np.array(community.best_partition(G).values())


########## EMBEDDING

def get_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
    from fa2 import ForceAtlas2
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
                  # Behavior alternatives
                  outboundAttractionDistribution=False,  # Dissuade hubs
                  linLogMode=False,  # NOT IMPLEMENTED
                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                  edgeWeightInfluence=edgeWeightInfluence,

                  # Performance
                  jitterTolerance=jitterTolerance,  # Tolerance
                  barnesHutOptimize=True,
                  barnesHutTheta=barnesHutTheta,
                  multiThreaded=False,  # NOT IMPLEMENTED

                  # Tuning
                  scalingRatio=scalingRatio,
                  strongGravityMode=False,
                  gravity=gravity,
                  # Log
                  verbose=verbose)

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    return positions

########## SPRING PREP

def save_hdf5_genes(E, gene_list, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_genes.hdf5"'''
    
    import h5py
    
    E = E.tocsc()
    
    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    cix_group = hf.create_group('cell_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iG, g in enumerate(gene_list):
        counts = E[:,iG].A.squeeze()
        cell_ix = np.nonzero(counts)[0]
        counts = counts[cell_ix]
        counts_group.create_dataset(g, data = counts)
        cix_group.create_dataset(g, data = cell_ix)

    hf.close()
    
def save_hdf5_cells(E, filename):
    '''SPRING standard: filename = main_spring_dir + "counts_norm_sparse_cells.hdf5" '''
    import h5py
    
    E = E.tocsr()
    
    hf = h5py.File(filename, 'w')
    counts_group = hf.create_group('counts')
    gix_group = hf.create_group('gene_ix')

    hf.attrs['ncells'] = E.shape[0]
    hf.attrs['ngenes'] = E.shape[1]

    for iC in range(E.shape[0]):
        counts = E[iC,:].A.squeeze()
        gene_ix = np.nonzero(counts)[0]
        counts = counts[gene_ix]
        counts_group.create_dataset(str(iC), data = counts)
        gix_group.create_dataset(str(iC), data = gene_ix)

    hf.close()
    
def save_sparse_npz(E, filename, compressed = False):
    ''' SPRING standard: filename = main_spring_dir + "/counts_norm.npz"'''
    E = E.tocsc()
    scipy.sparse.save_npz(filename, E, compressed = compressed)

def write_graph(filename, n_nodes, edges):
    nodes = [{'name':int(i),'number':int(i)} for i in range(n_nodes)]
    edges = [{'source':int(i), 'target':int(j), 'distance':0} for i,j in edges]
    out = {'nodes':nodes,'links':edges}
    open(filename,'w').write(json.dumps(out,indent=4, separators=(',', ': ')))

def write_edges(filename, edges):
    with open(filename, 'w') as f:
        for e in edges:
            f.write('%i;%i\n' %(e[0], e[1]))

def write_color_tracks(ctracks, fname):
    out = []
    for name,score in ctracks.items():
        line = name + ',' + ','.join(['%.3f' %x for x in score])
        out += [line]
    out = sorted(out,key=lambda x: x.split(',')[0])
    open(fname,'w').write('\n'.join(out))

def frac_to_hex(frac):
    rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3])*255,dtype=int))
    return '#%02x%02x%02x' % rgb

def get_color_stats_genes(color_stats, E, gene_list):
    means = E.mean(0).A.squeeze()
    stdevs = np.sqrt(sparse_var(E, 0))
    mins = E.min(0).todense().A1
    maxes = E.max(0).todense().A1
    
    pctl = 99.6
    pctl_n = (100-pctl) / 100. * E.shape[0]    
    pctls = np.zeros(E.shape[1], dtype=float)
    for iG in range(E.shape[1]):
        n_nonzero = E.indptr[iG+1] - E.indptr[iG]
        if n_nonzero > pctl_n:
            pctls[iG] = np.percentile(E.data[E.indptr[iG]:E.indptr[iG+1]], 100 - 100 * pctl_n / n_nonzero)
        else:
            pctls[iG] = 0
        color_stats[gene_list[iG]] = tuple(map(float, (means[iG], stdevs[iG], mins[iG], maxes[iG], pctls[iG])))
    return color_stats

def get_color_stats_custom(color_stats, custom_colors):
    for k,v in custom_colors.items():
        color_stats[k] = (np.mean(v),np.std(v),np.min(v),np.max(v),np.percentile(v,99))
    return color_stats

def save_color_stats(filename, color_stats):
    with open(filename,'w') as f:
        f.write(json.dumps(color_stats,indent=4, sort_keys=True).decode('utf-8'))

def build_categ_colors(categorical_coloring_data, cell_groupings):
    for k,labels in cell_groupings.items():
        label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
        categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
    return categorical_coloring_data

def save_cell_groupings(filename, categorical_coloring_data):
    with open(filename,'w') as f:
        f.write(json.dumps(categorical_coloring_data,indent=4, sort_keys=True).decode('utf-8'))

def save_spring_dir_sparse_hdf5(E,gene_list,project_directory, edges, custom_colors={}, cell_groupings={}):

    if not os.path.exists(project_directory):
        os.makedirs(project_directory)

    if not project_directory[-1] == '/': 
        project_directory += '/'

    # save custom colors
    custom_colors['Uniform'] = np.zeros(E.shape[0])
    write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')

    # create and save a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    color_stats = get_color_stats_genes(color_stats, E, gene_list)
    color_stats = get_color_stats_custom(color_stats, custom_colors)
    save_color_stats(project_directory + 'color_stats.json', color_stats)

    # save cell labels
    categorical_coloring_data = {}
    categorical_coloring_data = build_categ_colors(categorical_coloring_data, cell_groupings)
    save_cell_groupings(project_directory+'categorical_coloring_data.json', categorical_coloring_data)

    # write graph
    write_graph(project_directory + 'graph_data.json', E.shape[0], edges)
    write_edges(project_directory + 'edges.csv', edges)


#========================================================================================#

def make_spring_subplot(E, gene_list, save_path, base_ix = None, normalize = True, exclude_dominant_frac = 1.0, min_counts = 3, min_cells = 5, min_vscore_pctl = 75,show_vscore_plot = False, exclude_gene_names = None, num_pc = 30, sparse_pca = False, pca_norm = True, k_neigh = 4, cell_groupings = {}, num_force_iter = 100, output_spring = True, precomputed_pca = None, gene_filter = None, custom_colors = {}, exclude_corr_genes_list = None, exclude_corr_genes_minCorr = 0.2, dist_metric = 'euclidean', use_approxnn=False, run_doub_detector = False, dd_k=50, dd_frac=5, dd_approx=True, tot_counts_final = None):
    
    out = {}
    info_dict = {}
    info_dict['Date'] = '%s' %datetime.now()
    info_dict['Nodes'] = E.shape[0]
    info_dict['Num_Neighbors'] = k_neigh
    info_dict['Num_Force_Iter'] = num_force_iter
    info_dict['Gene_Var_Pctl'] = ''
    info_dict['Min_Cells'] = ''
    info_dict['Min_Counts'] = ''
    info_dict['Filtered_Genes'] = ''
    info_dict['Num_PCs'] = ''

    E = E.tocsc()
    if base_ix is None:
        base_ix = np.arange(E.shape[0])

    # total counts normalize
    if tot_counts_final is None:
        tot_counts_final = E.sum(1).A.squeeze()
    out['tot_counts_final'] = tot_counts_final

    if normalize:
        #print 'Normalizing'
        E = tot_counts_norm(E, exclude_dominant_frac = exclude_dominant_frac)[0]

    if precomputed_pca is None:
        if gene_filter is None:
            # Get gene stats (above Poisson noise, i.e. V-scores)
            #print 'Filtering genes'
            if (min_counts > 0) or (min_cells > 0) or (min_vscore_pctl > 0): 
                gene_filter = filter_genes(E, base_ix, min_vscore_pctl=min_vscore_pctl, min_counts=min_counts,min_cells=min_cells,show_vscore_plot = show_vscore_plot)

                info_dict['Gene_Var_Pctl'] = min_vscore_pctl
                info_dict['Min_Cells'] = min_cells
                info_dict['Min_Counts'] = min_counts
            else:
                gene_filter = np.arange(E.shape[1])


            if len(gene_filter) == 0:
                print 'Error: No genes passed filter'
                sys.exit(2)
                #print 'Error: All genes have mean expression < '+repr(min_exp) + ' or CV < '+repr(min_cv)
            #print 'Using %i genes' %(len(gene_filter))

            if not exclude_corr_genes_list is None:
                gene_filter = remove_corr_genes(E, gene_list, exclude_corr_genes_list, gene_filter, min_corr = exclude_corr_genes_minCorr)
                if len(gene_filter) == 0:
                    print 'Error: No genes passed filter'
                    sys.exit(2)

            # Remove user-excluded genes from consideration
            if not exclude_gene_names is None:
                keep_ix = np.array([ii for ii,gix in enumerate(gene_filter) if gene_list[gix] not in exclude_gene_names])
                #print 'Excluded %i user-provided genes' %(len(gene_filter)-len(keep_ix))
                gene_filter = gene_filter[keep_ix]
                if len(gene_filter) == 0:
                    print 'Error: No genes passed filter'
                    sys.exit(2)

        out['gene_filter'] = gene_filter
        info_dict['Filtered_Genes'] = len(gene_filter)
        # RUN PCA
        # if method == 'sparse': normalize by stdev
        # if method == anything else: z-score normalize
        #print 'Running PCA'
        num_pc = min(len(gene_filter), num_pc)
        Epca = get_pca(E[:,gene_filter], base_ix=base_ix, numpc=num_pc, keep_sparse=sparse_pca, normalize = pca_norm)
    else:
        Epca = precomputed_pca

    out['Epca'] = Epca
    out['num_pc'] = Epca.shape[1]
    info_dict['Num_PCs'] = Epca.shape[1]

    #print 'Building kNN graph'

    links, knn_graph = get_knn_graph(Epca, k=k_neigh, dist_metric = dist_metric, approx=use_approxnn)
    out['knn_graph'] = knn_graph

    if run_doub_detector:
        import doublet_detector as woublet
        #print 'Running woublet'
        doub_score, doub_score_full, doub_labels = woublet.detect_doublets([], counts=tot_counts_final, doub_frac=dd_frac, k=dd_k, use_approxnn=dd_approx, precomputed_pca=Epca)
        out['doub_score'] = doub_score
        out['doub_score_sim'] = doub_score_sim

    if output_spring:

        if not os.path.exists(save_path):
            os.makedirs(save_path)

        #print 'Saving SPRING files to %s' %save_path
        custom_colors['Total Counts'] = tot_counts_final
        np.savez_compressed(save_path + '/intermediates.npz', Epca = Epca, gene_filter = gene_filter, total_counts = tot_counts_final)

        if run_doub_detector:
            custom_colors['Doublet Score'] = doub_score

        if len(cell_groupings) > 0:
            save_spring_dir_sparse_hdf5(E, gene_list, save_path, list(links),
                            custom_colors = custom_colors,
                            cell_groupings = cell_groupings)
        else:
            save_spring_dir_sparse_hdf5(E, gene_list, save_path, list(links),
                            custom_colors = custom_colors)


    if num_force_iter > 0:
        positions = get_force_layout(links, Epca.shape[0], n_iter=num_force_iter, 
            edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, 
            jitterTolerance=1, verbose=False)
        positions = positions / 5.0
        positions = positions - np.min(positions, axis = 0) - np.ptp(positions, axis = 0) / 2.0
        positions[:,0] = positions[:,0]  + 750
        positions[:,1] = positions[:,1]  + 250
        out['coordinates'] = positions

    if output_spring:
        if num_force_iter > 0:
            np.savetxt(save_path + '/coordinates.txt',
                       np.hstack((np.arange(positions.shape[0])[:,None], positions)), fmt='%i,%.5f,%.5f')
         
        with open(save_path+'/run_info.json','w') as f:
            f.write(json.dumps(info_dict,indent=4, sort_keys=True).decode('utf-8'))
         
    return out

#========================================================================================#


############# PLOTTING

def gene_plot(x, y, E, gene_list, gene_name, col_range=(0,100), order_points=False, x_buffer=0, y_buffer=0,
        fig_size=(5,5), point_size=15, colormap='Reds', bg_color=[1,1,1], ax='', smooth_operator = []):
    '''
    Plot gene expression values on a scatter plot.

    Input
        x : x coordinates for scatter plot
        y : y coordinates for scatter plot
        E : gene expression counts matrix (cells x genes)
        gene_list (list of strings, length=n_cells): full list of gene names
        gene_name (string): name of gene to visualize
        col_range (float tuple, length=2): (color_floor, color_ceiling) percentiles
        order_points (boolean): if True, plot points with higher color values on top of points with lower values
        x_buffer (float): white space to add to x limits
        y_buffer (float): white space to add to y limits
        fig_size (float tuple, length=2): size of figure
        point_size (float): size of scatter plot points
        colormap: color scheme for coloring the scatter plot
        bg_color (RGB/HEX/color name): background color

    Output
        fig: figure handle
        ax: axis handle
        pl: scatter plot handle
    '''
    # get gene index and color data
    import matplotlib.pyplot as plt

    gene_ix = gene_list.index(gene_name)
    colordat = E[:,gene_ix].toarray()[:,0]

    if len(smooth_operator) > 0:
        colordat = np.dot(smooth_operator, colordat)

    # get min and max color values
    cmin = np.percentile(colordat, col_range[0])
    cmax = np.percentile(colordat, col_range[1])
    if cmax == 0:
        cmax = max(colordat)

    # order points by intensity, if desired
    if order_points:
        plot_ord = np.argsort(colordat)
    else:
        plot_ord = np.arange(len(colordat))

    # make the plot
    return_all = False
    if ax == '':
        return_all = True
        fig, ax = plt.subplots(1, 1, figsize = fig_size)

    pl = ax.scatter(x[plot_ord], y[plot_ord], c=colordat[plot_ord], s=point_size, edgecolor='none',
                    cmap=colormap, vmin=cmin, vmax=cmax)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim((min(x) - x_buffer, max(x) + x_buffer))
    ax.set_ylim((min(y) - y_buffer, max(y) + y_buffer))
    ax.patch.set_color(bg_color)

    if return_all:
        return fig, ax, pl
    else:
        return pl

def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii,0] = curcol[0] * .9
        cdat[ii,1] = curcol[1] * .9
        cdat[ii,2] = curcol[2] * .9
        cdat[ii,3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap

def custom_cmap(rgb_list):
    import matplotlib.pyplot as plt
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0],rgb_list)
    return cmap

def plot_groups(x, y, groups, lim_buffer = 50, saving = False, fig_dir = './', fig_name = 'fig', res = 300, close_after = False, title_size = 12, point_size = 3, ncol = 5):
    import matplotlib.pyplot as plt

    n_col = int(ncol)
    ngroup = len(np.unique(groups))
    nrow = int(np.ceil(ngroup / float(ncol)))
    fig = plt.figure(figsize = (14, 3 * nrow))
    for ii, c in enumerate(np.unique(groups)):
        ax = plt.subplot(nrow, ncol, ii+1)
        ix = groups == c

        ax.scatter(x[~ix], y[~ix], s = point_size, c = [.8,.8,.8], edgecolors = '')
        ax.scatter(x[ix], y[ix], s = point_size, c = [0,0,0], edgecolors = '')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([min(x) - lim_buffer, max(x) + lim_buffer])
        ax.set_ylim([min(y) - lim_buffer, max(y) + lim_buffer])

        ax.set_title(str(c), fontsize = title_size)

    fig.tight_layout()

    if saving:
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        plt.savefig(fig_dir + '/' + fig_name + '.png', dpi=res)

    if close_after:
        plt.close()


########## GENE ENRICHMENT

def rank_enriched_genes(E, gene_list, cell_mask, min_counts=3, min_cells=3):
    gix = (E[cell_mask,:]>=min_counts).sum(0).A.squeeze() >= min_cells
    print '%i cells in group' %(sum(cell_mask))
    print 'Considering %i genes' %(sum(gix))
    
    gene_list = gene_list[gix]
    
    z = sparse_zscore(E[:,gix])
    scores = z[cell_mask,:].mean(0).A.squeeze()
    o = np.argsort(-scores)
    
    return gene_list[o], scores[o]



