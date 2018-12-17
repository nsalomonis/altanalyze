
import sys,string,os
import traceback
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

try:
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        matplotlib.rcParams['backend'] = 'TkAgg'
        import matplotlib.pyplot as pylab
        matplotlib.rcParams['axes.linewidth'] = 0.5
        matplotlib.rcParams['pdf.fonttype'] = 42
        #matplotlib.rcParams['font.family'] = 'sans-serif'
        #matplotlib.rcParams['font.sans-serif'] = 'Arial'
        import numpy
except Exception:
    print traceback.format_exc()

import time
import random
import math
from stats_scripts import statistics
import ExpressionBuilder
import export

alphabet = map(chr, range(65, 91))+map(chr, range(97, 124)) ### Python magic

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    try:
        import unique ### local to AltAnalyze
        fn = unique.filepath(filename)
    except Exception:
        ### Should work fine when run as a script with this (AltAnalyze code is specific for packaging with AltAnalyze)
        dir=os.path.dirname(dirfile.__file__)
        try: dir_list = os.listdir(filename); fn = filename ### test to see if the path can be found (then it is the full path)
        except Exception: fn=os.path.join(dir,filename)
    return fn

def summarizeExpressionData(filename,qc_type):
    start_time = time.time()
    fn = filepath(filename)
    matrix=[]
    row_header=[]
    import RNASeq
    platform = RNASeq.checkExpressionFileFormat(fn,"3'array")
    x=0
    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0] =='#': x=0
        elif x==0:
            group_db, column_header, qc_db = assignGroupColors(t[1:],qc_type)
            x=1
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                if qc_type == 'distribution':
                    #values = map(lambda x: round(float(x), 1), t[1:]) ### report value to one decimal place
                    values = map(lambda x: float(x), t[1:])
                    i=0
                    for r in values:
                        if r!=0:
                            if 'counts' in dataset_name or platform == 'RNASeq':
                                r = round(math.log(r,2),1)
                            else:
                                r = round(r,1)
                            try:
                                qc_db[column_header[i]][r]+=1 ### count this rounded expression value once for this filename
                            except Exception:
                                qc_db[column_header[i]][r]=1
                        i+=1
                if qc_type == 'feature' or qc_type == 'totals':
                    if 'counts' in dataset_name:
                        feature_id = string.split(t[0],'=')[0]
                        if '-' in feature_id: feature = 'junction'
                        elif ':I' in feature_id: feature = 'intron'
                        elif ':E' in feature_id: feature = 'exon'
                        values = map(lambda x: float(x), t[1:])
                        i=0
                        for r in values:
                            if r!=0:
                                if qc_type == 'feature':
                                    r = round(math.log(r,2),1)
                                try:
                                    qc_db[column_header[i]][feature].append(r) ### add all expression values
                                except Exception:
                                    qc_db[column_header[i]][feature] = [r]
                            i+=1
            x+=1
            
    time_diff = str(round(time.time()-start_time,1))
    print 'Dataset import in %s seconds' % time_diff
    return qc_db,dataset_name

def reformatAltHeaders(headers):
    new_headers = []
    for i in headers:
        try: group, i = string.split(i,':')
        except Exception: pass
        new_headers.append(i)
    return new_headers
    
def importTableEntries(filename,filter_db,ensembl_exon_db,gene_db,root_dir,transpose,display,showIntrons,analysisType='plot'):
    import collections
    average_samples = True
    if showIntrons == 'yes': include_introns = True
    else: include_introns = False
    uid_db={} ### probeset or AltAnalyze RNA-Seq ID keyed
    uid_list={} ### ordered from first to last exon region
    uid_gene_db={} ### Lets us look at multiple genes
    try:
        import UI
        biotypes = UI.getBiotypes(filename)
    except Exception: biotypes={}
    for gene in ensembl_exon_db:
        uid_list[gene]=[]
        for (index,ed,id) in ensembl_exon_db[gene]:
            proceed = False
            if 'exp.' in filename:
                if include_introns:
                    proceed = True
                elif 'E' in ed.ExonID():
                    proceed = True
            else: ### Include introns for splicing index view
                if include_introns == True: proceed = True
                elif 'E' in ed.ExonID(): proceed = True
            if proceed:
                uid_db[id] = ed
                uid_list[gene].append(id)
            uid_gene_db[id]=gene

    if '_vs_' in filename: ### If one two groups, this is what will be output to the RawSplice folder - need to have this alternate way of getting the expression file location
        rootdir = string.split(filename, 'AltResults')[0]
        exp_dir = getValidExpFile(rootdir+'ExpressionInput')
        alt_groups_dir = string.split(exp_dir, 'ExpressionInput')[0]+'ExpressionInput/groups.'+findFilename(exp_dir)
        alt_groups_dir = string.replace(alt_groups_dir,'exp.','')
        
    start_time = time.time()
    fn = filepath(filename)
    matrix_gene_db={}
    stdev_gene_matrix_db={}
    row_header_gene={}
    ids={}
    x=0
    
    if 'heatmap' in analysisType:
        average_samples = False
        
    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
    for line in open(fn,'rU').xreadlines():         
        data = line.strip()
        t = string.split(data,'\t')
        if data[0]=='#': x=0
        elif x==0:
            if platform == 'RNASeq':
                removeExtension=True
            else:
                removeExtension=False
            group_db, column_header, sample_name_db = assignGroupColors(t[1:],'',removeExtension=removeExtension)
            x=1
            altresults = False
            if average_samples:
                if 'AltResults' in filename:
                    altresults=True
                    groups_dir = string.split(filename, 'AltResults')[0]+'ExpressionInput/groups.'+findFilename(filename)
                    if verifyFile(groups_dir)==False:
                        groups_dir = alt_groups_dir
                    new_column_header = reformatAltHeaders(t[3:])
                    start = 3
                else:
                    if 'exp.' in filename:
                        groups_dir = string.replace(filename,'exp.','groups.')
                    else:
                        groups_dir = string.replace(filename,'counts.','groups.')
                    new_column_header = column_header
                    start = 1 ### starting index with numeric values
                groups_dir = string.replace(groups_dir,'stats.','groups.')
                groups_dir = string.replace(groups_dir,'-steady-state.txt','.txt') ### groups is for the non-steady-state file
                
                try: group_index_db=collections.OrderedDict()
                except Exception:
                    import ordereddict
                    group_index_db = ordereddict.OrderedDict()
                ### use comps in the future to visualize group comparison changes
                sample_list,group_sample_db,group_db,group_name_sample_db,comp_groups,comps_name_db = ExpressionBuilder.simpleGroupImport(groups_dir)
                for item in sample_list:
                    group_name = group_db[item]
                    proceed=False
                    try: sample_index = new_column_header.index(item); proceed=True
                    except Exception:
                        try:
                            item = string.replace(item,'.bed','')
                            item = string.replace(item,'.CEL','') ### Probe-level analyses as RNA-Seq
                            item = string.replace(item,'.cel','')
                            item = string.replace(item,'.txt','')
                            item = string.replace(item,'.TXT','')
                            item = string.replace(item,'.TAB','')
                            item = string.replace(item,'.tab','')
                            sample_index = new_column_header.index(item)
                            proceed=True
                        except Exception:
                            pass
                            #print [item]
                            #print column_header
                            #print Error
                    if proceed:
                        try: group_index_db[group_name].append(sample_index)
                        except Exception:
                            try: group_index_db[group_name] = [sample_index] ### dictionary of group to input file sample indexes
                            except Exception: pass ### Occurs when analyzing splicing-index for two groups when more than two groups exist (error from 5 lines up)
                groups = map(str, group_index_db) ### store group names
                new_sample_list = map(lambda item: group_db[item], sample_list) ### lookup index of each sample in the ordered group sample list
                column_header = groups
            else:
                if 'AltResults' in filename: start = 3
                else: start = 1 ### starting index with numeric values
                column_header = t[start-1:]
            row_number=1   
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                uid = t[start-1]
                if ';' in uid:
                    uid = string.split(uid,';')[0]
                ids[uid]=None
                ens_geneID = string.split(uid,':')[0]
                #if ens_geneID in gene_db: print uid
                if uid in filter_db or ('heatmap' in analysisType and ens_geneID in gene_db):
                    try:
                        if len(biotypes)==1 and 'junction' in biotypes:
                            gene = ens_geneID
                        else:
                            gene = uid_gene_db[uid]
                        try: row_header_gene[gene].append(uid)
                        except Exception: row_header_gene[gene] = [uid]
                        if average_samples == False:
                            values = map(float,t[start:])
                            try: matrix_gene_db[gene].append(values)
                            except Exception: matrix_gene_db[gene]=[values]
                        else:
                            if platform == 'RNASeq' and altresults==False:
                                ### Convert to log2 RPKM values - or counts
                                values = map(lambda x: math.log(float(x),2), t[start:])
                            else:
                                values = map(float,t[start:])
                                
                            if 'AltResults' in filename: ### If splicing scores, normalize these to the mean values
                                mean = statistics.avg(values)
                                values = map(lambda x: x-mean, values)
                            avg_ls=[]; std_ls = []
                            for group_name in group_index_db:
                                group_values = map(lambda x: values[x], group_index_db[group_name]) ### simple and fast way to reorganize the samples
                                avg = statistics.avg(group_values)
                                try: st_err = statistics.stdev(group_values)/math.sqrt(len(group_values))
                                except Exception:
                                    ### Occurs if no replicates in the dataset
                                    st_err = 0
                                avg_ls.append(avg)
                                std_ls.append(st_err)
                            try: matrix_gene_db[gene].append(avg_ls)
                            except Exception: matrix_gene_db[gene]=[avg_ls]
                            try: stdev_gene_matrix_db[gene].append(std_ls)
                            except Exception: stdev_gene_matrix_db[gene]=[std_ls]
                    except Exception:
                        #print traceback.format_exc()
                        pass
            x+=1

    global colors
    original_column_header = list(column_header)
    if len(uid_list)==0:
        print 'No genes found in the exon expression database'; forceNoExonExpError
    successfully_output_genes=0
    display_count=0 ### Only display a certain number of genes
    
    for last_gene in uid_list: pass
    for gene in uid_list:
        fig = pylab.figure() ### Create this here - resulting in a single figure for memory purposes
        new_header = []
        new_matrix = []
        new_stdev = []
        annotation_list=[]
        gene_symbol = gene_db[gene]
        try: matrix = matrix_gene_db[gene]
        except Exception:
            #print gene_symbol, 'not in alternative expression database'
            continue ### go the next gene - no alt.expression for this gene
        row_header = row_header_gene[gene]

        try: stdev_matrix = stdev_gene_matrix_db[gene]
        except Exception: pass
        for uid in uid_list[gene]:
            #print row_header;sys.exit()
            try:
                i = row_header.index(uid) ### If the ID is in the filtered annotated exon list (not just core)
                new_header.append(uid)
                try: new_matrix.append(matrix[i])
                except Exception: print uid, i,len(matrix);sys.exit()
                ed = uid_db[uid]
                annotation_list.append(ed)
                try: new_stdev.append(stdev_matrix[i])
                except Exception: pass
            except Exception: pass

        if len(new_matrix)>0:
            matrix = new_matrix
        if len(new_header)>0:
            row_header = new_header
        if 'heatmap' in analysisType:
            export_dir = root_dir + gene_symbol + '-heatmap.txt'
            export_obj = export.ExportFile(export_dir)
            export_obj.write(string.join(column_header,'\t')+'\n')
            ki=0
            if len(annotation_list)>0:
                for ed in annotation_list:
                    if 'AltResults' not in filename and platform == 'RNASeq':
                        values = map(lambda x: math.log(x,2), matrix[ki])
                    else: values = matrix[ki]
                    export_obj.write(string.join([ed.ExonID()] + map(str,values),'\t')+'\n')
                    ki+=1
                row_metric = 'euclidean'; row_method = None
            else:
                ### Just junctions analyzed here... no sorted junctions yet
                ki=0
                for uid in row_header_gene[gene]:
                    if 'AltResults' not in filename and platform == 'RNASeq':
                        values = map(lambda x: math.log(x,2), matrix[ki])
                    else: values = matrix[ki]
                    export_obj.write(string.join([uid] + map(str,values),'\t')+'\n')
                    ki+=1
                row_metric = 'euclidean'; row_method = 'average'
            export_obj.close()
            from visualization_scripts import clustering
            
            column_metric = 'euclidean'; column_method = 'hopach'
            color_gradient = 'red_black_sky'; transpose = False; graphic_links=[]
            if ki>100: transpose = True
            if gene == last_gene: display = True
            else: display = False
            graphic_links = clustering.runHCexplicit(export_dir, graphic_links, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=display, Normalize=True, compressAxis = False, contrast = 2.5)
            successfully_output_genes+=1
        else:
            stdev_matrix = new_stdev
            time_diff = str(round(time.time()-start_time,1))
            #print '%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
            if transpose == True:
                matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
                column_header, row_header = row_header, original_column_header
                stdev_matrix = map(numpy.array, zip(*stdev_matrix))
            matrix = numpy.array(matrix)

            stdev_matrix = numpy.array(stdev_matrix)
            try:
                if len(uid_list)>10:
                    #if display_count==5: display=False
                    display=False
                if display_count==0:
                    ### store a consistent color palete to use
                    colors=[]
                    """
                    k=0
                    while k < len(row_header):
                        colors.append(tuple(rand(3)))
                        k+=1"""
                    #http://stackoverflow.com/questions/3016283/create-a-color-generator-from-given-colormap-in-matplotlib
                    cm = pylab.cm.get_cmap('gist_rainbow') #gist_ncar
                    for i in range(len(row_header)):
                        colors.append(cm(1.*i/len(row_header)))  # color will now be an RGBA tuple
        
                plotExonExpression(fig,matrix,stdev_matrix,row_header,column_header,dataset_name,annotation_list,gene_symbol,root_dir,display=display)
                successfully_output_genes+=1
                display_count+=1
            except Exception:
                print traceback.format_exc();sys.exit()
                print gene_symbol, 'failed'
        try: pylab.close()
        except Exception: pass
        if successfully_output_genes>0:
            #try: print 'Gene graphs exported to ExonPlots...'
            #except Exception: pass
            pass
        else:
            print '\nWARNING!!!! No genes with associated alternative exon evidence found\n'; forceNoExonExpError
        try:
            import gc
            fig.clf()
            pylab.close()
            gc.collect()
        except Exception:
            pass
        
def importDataSimple(filename,transpose):
    start_time = time.time()
    fn = filepath(filename)
    matrix=[]
    row_header=[]
    x=0
    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if data[0]=='#': x=0
        elif x==0:
            group_db, column_header, sample_name_db = assignGroupColors(t[1:],'')
            x=1
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                row_header.append(t[0])
                t = map(float,t[1:])
                if (abs(max(t)-min(t)))>0:
                    matrix.append(t)
            x+=1
            
    time_diff = str(round(time.time()-start_time,1))
    print '%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    if transpose == True:
        matrix = map(numpy.array, zip(*matrix)) ### coverts these to tuples
        column_header, row_header = row_header, column_header
    return numpy.array(matrix), column_header, row_header, dataset_name, group_db
  
def assignGroupColors(t,qc_type,removeExtension=True):
    """ Assign a unique color to each group """
    
    k = 0
    sample_name_db={}
    column_header=[]; group_db={}; color_db={}
    color_list = ['r', 'b', 'y', 'g', 'w', 'k', 'm']
    for i in t:
        if len(i)>4:
            if i[-4] == '.' and removeExtension:
                i = i[:-4] ### remove the .cel, .txt, .tab or .bed
            if ':' in i:
                group,i = string.split(i,':')
                try: color = color_db[group]
                except Exception:
                    try: color_db[group] = color_list[k]
                    except Exception:
                        ### If not listed in the standard color set add a new random color
                        rgb = tuple(rand(3)) ### random color
                        color_list.append(rgb)
                        color_db[group] = color_list[k]
                    color = color_db[group]
                    k+=1
                group_db[i] = group, color
        column_header.append(i)
        sample_name_db[i]=db={} ### Initialize a dictionary within the filename dictionary
    return group_db, column_header, sample_name_db

def plotNormalizationResults(matrix,row_headers,column_headers,dataset_name,plot_type):
    ### See - http://bib.oxfordjournals.org/content/early/2011/04/15/bib.bbq086.full Fig.1
    fig = pylab.figure()
    if plot_type == 'pm_mean':
        pylab.title('Average Raw Intensity Signal - %s' % dataset_name)
        pylab.ylabel('Signal Intensity')
    else:
        pylab.title('Deviation of Residuals from Median - %s' % dataset_name)
        pylab.ylabel('Mean absolute deviation')
    pylab.subplots_adjust(left=0.125, right=0.95, top=0.9, bottom=0.40)
    ax = fig.add_subplot(111)

    i = row_headers.index(plot_type)
    y_axis = matrix[i]
    pylab.plot(range(len(y_axis)),y_axis,'b*',marker='s',markersize=7)
    
    ### Increase the max of the y-axis to accomidate the legend
    pylab.ylim(ymin=0)
    new_max = increaseYaxis(max(y_axis))
    pylab.ylim(ymax=new_max)
    #ax.yaxis.set_ticks([min(all_totals),max(all_totals)])
    ax.xaxis.set_ticks([-1]+range(len(y_axis)+1))
    xtickNames = pylab.setp(pylab.gca(), xticklabels=['']+column_headers+[''])
    pylab.setp(xtickNames, rotation=90, fontsize=10)
    
    filename = 'QC-%s-%s.pdf' % (dataset_name,plot_type)
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename) #,dpi=200
    graphic_link.append(['QC - '+plot_type,root_dir+filename])
    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass


def plotExonExpression(fig,matrix,stdev_matrix,row_headers,column_headers,dataset_name,annotation_list,gene_symbol,root_dir,display=True):
    """ Display exon-level expression for splicing-index, RPKMs or Affymetrix probeset intensities """
    #print len(matrix);sys.exit()
    ax = fig.add_subplot(111)
    print '.',
    if 'exp.' in dataset_name:
        datatype = '-RawExonExp'
        pylab.ylabel('Exon Expression (log2)')
        pylab.title(gene_symbol+' Exon Expression - '+dataset_name)
    else:
        datatype = '-AltExonExp'
        pylab.ylabel('Splicing Index Fold Change')
        pylab.title(gene_symbol+' Exon Splicing Index - '+dataset_name)
    pylab.xlabel('Exon Regions')

    data_max=[]
    data_min=[]
    for values in matrix:
        data_max.append(max(values))
        data_max.append(abs(min(values)))
        data_min.append(min(values))

    ### Increase the max of the y-axis to accomidate the legend
    try: new_max = increaseYaxis(max(data_max))
    except Exception: print len(matrix);sys.exit()
    total_min = min(data_min)
    if max(data_max) < 3:
        new_max = 5
        if total_min > -3: total_min = -4
    else:
        new_max = new_max/1.3
    pylab.ylim(ymax=new_max,ymin=total_min)
    
    exon_annotations = []
    for ed in annotation_list:
        exon_annotations.append(ed.ExonID())
    if len(column_headers) != len(exon_annotations):
        exon_annotations = map(lambda x: string.split(x,':')[-1],column_headers)
    i=0
    color_list = ['r', '#ff9900', 'y', 'g', 'b', '#ff66cc', '#9900ff', '#996633', '#336666']
    font_ratio = 40.00/len(exon_annotations) ### at about 30 exons a font of 14 looks OK
    if font_ratio<1: fontsize = font_ratio*14
    else: fontsize = 12
    for sample_name in row_headers:
        #name = '(%s) %s' % (alphabet[i],sample_name)
        if len(row_headers)<10: color = color_list[i]
        else:
            color = colors[i]
        x_axis = list(range(0, len(column_headers), 1)) ### Numeric x-axis values (can be spaced differently)
        y_axis = matrix[i]
        err = stdev_matrix[i]
        if sum(err)==0:
            pylab.errorbar(x_axis,y_axis,color=color,linewidth=1.75,label=sample_name)
        else:
            #pylab.plot(x_axis,y_axis,color=color,linewidth=1.75)
            pylab.errorbar(x_axis,y_axis,yerr=err,color=color,linewidth=1.75,label=sample_name)
        #k = y_axis.index(max(y_axis))
        #pylab.text(x_axis[k],y_axis[k],alphabet[i],fontsize=12) ### Plot this on the maximum value
        a = [-1]+range(len(y_axis))+[len(y_axis)]
        ax.xaxis.set_ticks(a)
        xtickNames = pylab.setp(pylab.gca(), xticklabels=['']+exon_annotations)
        pylab.setp(xtickNames, rotation=90, fontsize=fontsize)
        i+=1
    loc = "upper center"
    size = 12
    ncol = 1
    if len(row_headers)>5:
        ncol = 2
        if len(row_headers)>20:
            loc = 'upper center'
            ncol = 3
            size = 8
        if len(row_headers)>30:
            size = 5.5
            ncol = 4
        if len(row_headers)>40:
            size = 5
            ncol = 4
    
    if loc == 'upper center':
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height*0.6])
        try: ax.legend(loc=loc, ncol=ncol, bbox_to_anchor=(0., 1.02,1.,.8),fontsize = size) ### move the legend over to the right of the plot
        except Exception:
            ### Older versions of Matplotlib
            try:
                pylab.legend(prop={fontsize: 'small'})
                ax.legend(loc=loc, ncol=ncol, bbox_to_anchor=(0., 1.02,1.,.8))
            except Exception:
                pass
            
    else:
        pylab.legend(loc=loc, ncol=ncol, prop={'size': size})

    filename = gene_symbol+datatype+'.pdf'
    pylab.savefig(root_dir + filename)
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename,dpi=120) #,dpi=200
    if display:
        pylab.show()
    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

def plotTotalExpression(qc_db,dataset_name,features):
    fig = pylab.figure()    
    #pylab.xlabel('Biological Sample Names')
    pylab.ylabel('Total Number of Reads')
    pylab.title('Total Expression for Transcript Features - %s' % dataset_name)
    pylab.subplots_adjust(left=0.125, right=0.95, top=0.9, bottom=0.40)
    ax = fig.add_subplot(111)
    
    color_db = {}
    color_db['exon']='r*'
    color_db['junction']='b*'
    color_db['intron']='y*'
    feature_summary_db={}
    samples=[]
    for sample_name in qc_db:
        samples.append(sample_name)
    samples.sort()
    
    all_totals=[]
    for sample_name in samples:
        ls=[]; x_ls=[]; y_ls=[]
        for feature_type in qc_db[sample_name]:
            total_exp = sum(qc_db[sample_name][feature_type])
            try:
                feature_summary_db[feature_type].append(total_exp)
            except Exception:
                feature_summary_db[feature_type] = [total_exp]
            all_totals.append(total_exp)
    for feature_type in feature_summary_db:
        y_axis = feature_summary_db[feature_type]
        pylab.plot(range(len(y_axis)),y_axis,color_db[feature_type],marker='o',markersize=7,label=feature_type)

    ### Increase the max of the y-axis to accomidate the legend
    new_max = increaseYaxis(max(all_totals))
    pylab.ylim(ymax=new_max)
    #ax.yaxis.set_ticks([min(all_totals),max(all_totals)])
    ax.xaxis.set_ticks([-1]+range(len(y_axis)+1))
    xtickNames = pylab.setp(pylab.gca(), xticklabels=['']+samples+[''])
    pylab.setp(xtickNames, rotation=90, fontsize=10)
    
    pylab.legend(loc="upper right", prop={'size': 11})

    filename = 'QC-%s-TotalFeatureExpression.pdf' % dataset_name
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename)
    graphic_link.append(['QC - Total Feature Expression',root_dir+filename])
    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass


def increaseYaxis(max):
    count = len(str(int(max)))-2
    round_down_max = round(max,-1*count) ### Rounds major unit down
    #new_max = round_down_max+(round_down_max/2)
    new_max = max*2
    return new_max
    
def plotFeatureBoxPlots(qc_db,dataset_name,feature_type):
    pylab.figure()    
    pylab.xlabel('Biological Sample Names')
    pylab.ylabel('Read Counts - Log2')
    pylab.title('Expression BoxPlots for %ss - %s' % (feature_type,dataset_name))
    #pylab.subplots_adjust(left=0.085, right=0.95, top=0.2, bottom=0.35)
    pylab.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.35)
    
    #axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    #pylab.axis(axes)

    boxplots=[]
    samples=[]
    
    sample_sorted_list=[]
    
    for sample_name in qc_db:
        try: qc = qc_db[sample_name][feature_type]
        except Exception:
            print 'No junction data found for at least one sample:',sample_name; forceExit
        sample_sorted_list.append([statistics.avg(qc),statistics.stdev(qc),sample_name])
    sample_sorted_list.sort()
    sample_sorted_list.reverse()
    
    filename = 'QC-%s-BoxPlot-%s.pdf' % (dataset_name,feature_type)
    export_obj = export.ExportFile(root_dir + filename[:-4]+'.txt')
    export_obj.write('SampleID\tAverage Expression\n')
    
    firstEntry=True
    for (mean,stdev,sample_name) in sample_sorted_list:
        ls=[]; x_ls=[]; y_ls=[]
        qc = qc_db[sample_name][feature_type]
        boxplots.append(qc)
        samples.append(sample_name)
        export_obj.write(sample_name+'\t'+str(mean)+'\n')
        if firstEntry:
            threshold=mean-2*stdev
            firstEntry=False
        else:
            if mean<threshold:
                print sample_name,'expression is considered very low (2 standard deviations away from the max).'
    pylab.boxplot(boxplots, notch=0, whis=1.5, positions=None, widths=None, patch_artist=False)
    #pylab.boxplot(boxplots, notch=0, sym='+', vert=1, whis=1.5, positions=None, widths=None, patch_artist=False)
    xtickNames = pylab.setp(pylab.gca(), xticklabels=samples)
    pylab.setp(xtickNames, rotation=90, fontsize=10)
    export_obj.close()

    #print 'Exporting:',filename
    pylab.savefig(root_dir + filename)
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename) #,dpi=200
    graphic_link.append(['QC - BoxPlot-'+feature_type+' Expression',root_dir+filename])
    try:
        import gc
        pylab.figure.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

def rand(i):
    ### There is a pylab and scipy rand which does not appear to be in numpy or matplotlib
    return map(lambda x: random.random(), ['']*i)
    
def plotExpressionDistribution(qc_db,dataset_name):
    pylab.figure()
    pylab.xlabel('Log2 Expression (x10-1 Precision)')
    pylab.ylabel('Number of Observations')
    pylab.title('Biological Sample Expression Distribution - '+dataset_name)
    
    #axes = getAxes(scores) ### adds buffer space to the end of each axis and creates room for a legend
    #pylab.axis(axes)

    i=0
    color_list = ['r', 'b', 'y', 'g', 'k', 'm']
    sample_list=[]
    for sample_name in qc_db:
        sample_list.append(sample_name)
    sample_list.sort()
    for sample_name in sample_list:
        ls=[]; x_ls=[]; y_ls=[]
        qc = qc_db[sample_name]
        try: code = alphabet[i]
        except Exception: code = str(i)
        name = '(%s) %s' % (code,sample_name)

        for x in qc:
            ls.append((x,qc[x]))
        ls.sort() ### Get all x,y values into a sorted list and then provide those two lists to plot
        for (x,y) in ls:
            x_ls.append(x); y_ls.append(y)
        if len(qc_db)<7: color = color_list[i]
        else: color = tuple(rand(3))
        pylab.plot(x_ls,y_ls,color=color,label=name,linewidth=1.75)
        try:
            k = y_ls.index(max(y_ls))
            pylab.text(x_ls[k],y_ls[k],code,fontsize=7) ### Plot this on the maximum value
        except Exception: pass
        i+=1
    if len(sample_list)<15:
        font_size = 11
    elif len(sample_list)<25:
        font_size = 8
    elif len(sample_list)<45:
        font_size = 5.5
    else:
        font_size = 5
    pylab.legend(loc="upper right", prop={'size': font_size})

    filename = 'QC-%s-distribution.pdf' % dataset_name
    pylab.savefig(root_dir + filename)
    #print 'Exporting:',filename
    filename = filename[:-3]+'png'
    pylab.savefig(root_dir + filename) #,dpi=200
    graphic_link.append(['QC - Expression Distribution',root_dir+filename])
    try:
        import gc
        pylab.figure.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

def getAxes(x_axis,y_axis):
    x_axis_min = min((x_axis)+min(x_axis)/2.5)
    x_axis_max = min((x_axis)+max(x_axis)/5)
    y_axis_min = min((y_axis)+min(y_axis)/5)
    y_axis_max = max((x_axis)+max(x_axis)/5)
    return [x_axis_min, x_axis_max, y_axis_min, y_axis_max]
    
def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]

def findFilename(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[x:]
    
def verifyFile(filename):
    fn=filepath(filename)
    try:
        for line in open(fn,'rU').xreadlines(): found = True; break
    except Exception: found = False
    return  found

def outputArrayQC(filename):
    """ QC plots for Affymetrix array analysis. The file is all probeset expression values (exp.DatasetName) """
    
    global root_dir
    global graphic_link
    graphic_link = []
    
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionInput','DataPlots')
    try: os.mkdir(root_dir)
    except Exception: null=[] ### dir exists
    
    try:
        ### If Affymetrix RMA result summaries available, use these
        if '/' in filename: delim = '/'
        else: delim = '\\'
        if 'ExpressionInput' in filename:
            apt_dir = string.split(filename,'ExpressionInput')[0]+'ExpressionInput/APT-output/rma-sketch.report.txt'                
        else:
            apt_dir = string.join(string.split(filename,delim)[:-1],delim)+'/ExpressionInput/APT-output/rma-sketch.report.txt'
        if verifyFile(apt_dir)==False:
            apt_dir = string.replace(apt_dir,'-sketch','')
        transpose = True
        matrix, column_header, row_header, dataset_name, group_db = importDataSimple(apt_dir,transpose)
        plotNormalizationResults(matrix,row_header,column_header,dataset_name,'pm_mean')
        plotNormalizationResults(matrix,row_header,column_header,dataset_name,'all_probeset_mad_residual_mean')
    except Exception:
        null=[]

    qc_type = 'distribution'
    qc_db,dataset_name = summarizeExpressionData(filename,qc_type)
    plotExpressionDistribution(qc_db,dataset_name)
    return graphic_link
    
def outputRNASeqQC(filename):
    """ QC plots for RNA-Seq analysis. The file should be exon-level read counts"""
    global root_dir
    global graphic_link
    graphic_link = []
    root_dir = findParentDir(filename)
    root_dir = string.replace(root_dir,'ExpressionInput','DataPlots')
    try: os.mkdir(root_dir)
    except Exception: null=[] ### dir exists
        
    qc_type = 'distribution'
    ### Distribution of each bin of expression values (log2 binned at the single decimal level)
    qc_db,dataset_name = summarizeExpressionData(filename,qc_type)
    plotExpressionDistribution(qc_db,dataset_name)
    
    qc_type = 'feature'
    qc_db,dataset_name = summarizeExpressionData(filename,qc_type)
    
    ### BoxPlot of expression values for each feature
    features=[]
    for s in qc_db:
        for feature in qc_db[s]:
            features.append(feature)
        break
    for feature in features:
        plotFeatureBoxPlots(qc_db,dataset_name,feature)
    
    qc_type = 'totals'
    ### Total expression (total reads) for each feature
    plotTotalExpression(qc_db,dataset_name,features)
    return graphic_link

def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>9: break
    except Exception: null=[]
    return count

def getValidExpFile(altanalyze_rawexp_dir):
    import unique
    dir_files = unique.read_directory(altanalyze_rawexp_dir)
    valid_file = ''
    for file in dir_files:
        if 'exp.' in file and 'state.txt' not in file and 'feature' not in file:
            valid_file = altanalyze_rawexp_dir+'/'+file
            break
    return valid_file

def displayExpressionGraph(species,Platform,exp_file,gene,transpose,display=True,showIntrons=False,analysisType='plot'):
    ### Get gene annotations (users can provide an Ensembl or symbol)
    print 'Importing exon-level expression data for visualization (be patient)...'
    from build_scripts import ExonAnalyze_module
    global platform
    platform = Platform
    if platform != 'AltMouse': gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
    else: gene_annotation_file = "AltDatabase/"+species+"/"+platform+"/"+platform+"_gene_annotations.txt"

    genes=[]
    gene=string.replace(gene,'|',',')
    gene=string.replace(gene,' ',',')
    if ',' in gene:
        genes += string.split(gene,',')
    else: genes.append(gene)    
    gene_db={}
    for gene in genes:
        try:
            if 'ENS' in gene:
                try: annotate_db ### If variable is defined
                except Exception:
                    annotate_db = ExonAnalyze_module.import_annotations(gene_annotation_file,platform,keyBySymbol=False) ### Make an SQLite call
                gene_symbol = annotate_db[gene].Symbol()
            else:
                try: annotate_db ### If variable is defined
                except Exception:
                    annotate_db = ExonAnalyze_module.import_annotations(gene_annotation_file,platform,keyBySymbol=True)
                gene_symbol = gene
                gene = annotate_db[gene].GeneID()
            gene_db[gene]=gene_symbol
        except Exception:
            #if len(gene)>0: print gene, 'not in database'
            pass
        
    if len(gene_db)==0:
        force_no_gene_found_error
    if 'AltResults' in exp_file:
        root_dir = string.split(exp_file,'AltResults')[0]+'ExonPlots/'
    else:
        root_dir = string.split(exp_file,'ExpressionInput')[0]+'ExonPlots/'
        
    from build_scripts import ExonAnalyze_module
    if platform == 'RNASeq': datatype = 'exons'
    else: datatype = 'probesets'
    export_exon_filename = 'AltDatabase/'+species+'/'+platform+'/'+species+'_Ensembl_'+datatype+'.txt'
    if verifyFileLength(export_exon_filename) == 0:
        rootdir = string.replace(root_dir,'ExonPlots/','')
        export_exon_filename = rootdir+'/'+export_exon_filename
    
    from build_scripts import ExonArrayEnsemblRules
    ensembl_exon_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'gene-probesets',gene_db) ### Make an SQLite call
    
    filter_db = {}
    for gene in ensembl_exon_db:
        ensembl_exon_db[gene].sort()
        
        for (index,ed,id) in ensembl_exon_db[gene]:
            filter_db[id] = []
            
    try: os.mkdir(root_dir)
    except Exception: None ### dir exists
    print 'Image results being saved to the folder "ExonPlots" in the AltAnalyze results directory.'
    importTableEntries(exp_file,filter_db,ensembl_exon_db,gene_db,root_dir,transpose,display,showIntrons,analysisType=analysisType) ### Make an SQLite call
 
if __name__ == '__main__':
    file = "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/Diabetes-Blood/ExpressionInput/exp.diabetes.txt"
    file = "/Users/nsalomonis/Desktop/User Diagnostics/Mm_spinal_cord_injury/AltResults/RawSpliceData/Mm/splicing-index/Mm_spinal.txt"
    #file = "/Users/nsalomonis/Desktop/User Diagnostics/Mm_spinal_cord_injury/ExpressionInput/exp.Mm_spinal.txt"
    file = "/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/AltResults/RawSpliceData/Mm/splicing-index/test.txt"
    file = '/Users/saljh8/Desktop/Archived/Desktop/dataAnalysis/CPMC/Liliana/CPMC_GB-samples/Cultured/AltResults/RawSpliceData/Hs/splicing-index/Hs_Exon_CBD_vs_Vehicle.p5_average.txt'
    file = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/Lattice/Full/AltResults/RawSpliceData/Mm/splicing-index/myeloblast.txt'
    file = '/Volumes/SEQ-DATA/Grimes/14018_gmp-pro/Lattice/Full/ExpressionInput/exp.myeloblast.txt'
    file = '/Volumes/salomonis1/projects/Grimes/GEC_14061/bams/ExpressionInput/counts.gec_14061.txt'
    #file = '/Volumes/SEQ-DATA/SingleCell-Churko/ExpressionInput/exp.CM.txt'
    outputRNASeqQC(file);sys.exit()
    species='Hs'
    Platform="RNASeq"
    #Platform="exon"
    ShowIntrons = 'yes'
    displayExpressionGraph(species,Platform,file,'ENSG00000081189',True,showIntrons=ShowIntrons);sys.exit() #ENSG00000140009 ENSG00000152284 ENSG00000133794 (ARNTL)
    outputRNASeqQC(file)
    