import sys, string
import export
import math
import random
import copy
import os
import os.path
import unique
R_present=True

try:
    ### If file is present use this location
    loc = unique.filepath('Config/R_location.txt')
    s = open(loc,'r')
    useStaticLocation=s.read()
    #print useStaticLocation
    #print 'Using the Config designated location
    for p in os.environ['PATH'].split(':'): ### For Unix cluster environments
        if os.path.exists(p + '/R'):
            #path = p + '/R'
            useStaticLocation = False
            
except Exception:
    #print 'NOT using the Config designated location'
    useStaticLocation = False
    
try:
    forceError ### This doesn't currently work with the compiled version of AltAnalyze
    import rpy2.robjects as robjects
    r = robjects.r
    print "\n---------Using RPY2---------\n"  
except Exception:
    from pyper import *
    #print "\n---------Using PypeR---------\n"
    ### Running the wrong one once is fine, but multiple times causes it to stall in a single session
    try:
        try:
            if 'Xdarwin' in sys.platform: ### Xdarwin is indicated since this if statement is invalid without a stand-alone Mac R package (ideal)
                #print 'Using AltAnalyze local version of R'
                #print 'A'
                path = unique.filepath("AltDatabase/tools/R/Mac/R")
                r = R(RCMD=path,use_numpy=True)
            elif os.name == 'nt':
                path = unique.filepath("AltDatabase/tools/R/PC/bin/x64/R.exe")
                r = R(RCMD=path,use_numpy=True)    
            else:
                #print 'B'
                if useStaticLocation == False or useStaticLocation=='no':
                    print 'NOT using static location'
                    r = R(use_numpy=True)
                else:
                    print 'Using static location'
                    path = '/usr/local/bin/R'
                    if os.path.exists(path): pass
                    else:
                        path = '/usr/bin/R'
                    if os.path.exists(path):
                        print 'Using the R path:',path
                        r = R(RCMD=path,use_numpy=True)
                    else:
                        r = None
                        R_present=False
                        print 'R does not appear to be installed... Please install first.'
        except Exception:
            #print 'C'
            r = R(use_numpy=True)

    except Exception:
        print traceback.format_exc()
        r = None
        R_present=False
        pass

LegacyMode = True
### Create a Directory for R packages in the AltAnalyze program directory (in non-existant)
r_package_path = string.replace(os.getcwd()+'/Config/R','\\','/') ### R doesn't link \\
r_package_path = unique.filepath(r_package_path) ### Remove the AltAnalyze.app location
try: os.mkdir(r_package_path)
except Exception: None

if R_present:
    ### Set an R-package installation path
    command = '.libPaths("'+r_package_path+'")'; r(command) ### doesn't work with %s for some reason
    #print_out = r('.libPaths()');print print_out; sys.exit()

def remoteMonocle(input_file,expPercent,pval,numGroups):
    #input_file="Altanalyze" 
    setWorkingDirectory(findParentDir(input_file)[:-1])
    try: os.mkdir(findParentDir(input_file)[:-1])
    except Exception: None
    z = RScripts(input_file)
    setWorkingDirectory(input_file)
    z.Monocle(input_file,expPercent,pval,numGroups)
    
def remoteHopach(input_file,cluster_method,metric_gene,metric_array,force_array='',force_gene=''):
    """ Run Hopach via a call from an external clustering and visualizaiton module """
    #input_file = input_file[1:] #not sure why, but the '\' needs to be there while reading initally but not while accessing the file late
    row_order = []
    column_order = []
    if 'ICGS-SubCluster' in input_file:
        force_array=2
    input_file = checkForDuplicateIDs(input_file) ### Duplicate IDs will cause R to exit when creating the data matrix
    z = RScripts(input_file)
    setWorkingDirectory(input_file)
    z.Hopach(cluster_method,metric_gene,force_gene,metric_array,force_array)

    if cluster_method == 'both' or cluster_method == 'gene':
        filename = findParentDir(input_file)+'/hopach/rows.'+findFileName(input_file)
        row_order = importHopachOutput(filename)
    if cluster_method == 'both' or cluster_method == 'array':
        filename = findParentDir(input_file)+'/hopach/columns.'+findFileName(input_file)
        column_order = importHopachOutput(filename)
    #print row_order; sys.exit()
    return input_file, row_order, column_order

def remoteAffyNormalization(input_file,normalization_method,probe_level,batch_effects):
    ### Input file is the path of the expression output from normalization
    setWorkingDirectory(findParentDir(input_file)[:-1])
    try: os.mkdir(findParentDir(input_file)[:-1])
    except Exception: None #Already exists
    z = RScripts(input_file)
    z.AffyNormalization(normalization_method,probe_level,batch_effects)
    
def checkForDuplicateIDs(input_file, useOrderedDict=True):
    
    if 'SamplePrediction' in input_file or '-Guide' in input_file:
        ### OrderedDict is prefered but will alter prior ICGS results
        useOrderedDict = False
    first_row = True
    import collections
    if useOrderedDict:
        try: key_db = collections.OrderedDict()
        except Exception:
            try:
                import ordereddict
                key_db = ordereddict.OrderedDict()
            except Exception:
                key_db={}
    else:
        key_db={}
    key_list=[]
    fn=filepath(input_file)
    offset=0
    nonNumericsPresent=False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if first_row == True:
            if ('row_clusters-flat' in t and 'row_clusters-flat' not in t[0]):
                headers = string.join(['uid']+t[2:],'\t')+'\n'
                offset = 1
            elif '-filtered.txt' in fn and ".R2." in t[1] and LegacyMode:
                headers = string.join(['uid']+t[2:],'\t')+'\n'
                offset = 1
            else:
                headers = line
            first_row = False
        else:
            key = t[0]
            try:
                k1,k2string.split(key,' ')
                print [k1, k2],
                if k1==k2: key = k1
                print key
            except Exception: pass
            if key!='column_clusters-flat':
                key_list.append(key)
                
                try: s = map(float,t[offset+1:])
                except Exception:
                    nonNumericsPresent=True
                key_db[key]=t
    
    if nonNumericsPresent:
        import numpy
        for key in key_db:
            t = key_db[key]
            s=[key]
            if offset ==1: s.append('')
            temp=[]
            for value in t[offset+1:]:
                try: temp.append(float(value))
                except Exception: pass
            avg=numpy.mean(temp)
            for value in t[offset+1:]:
                try: s.append(str(float(value)-avg))
                except Exception: s.append('0.000101')
            key_db[key]=s
            
    if len(key_db) != len(key_list) or offset>0 or nonNumericsPresent:
        print 'Writing a cleaned-up version of the input file:'
        ### Duplicate IDs present
        input_file = input_file[:-4]+'-clean.txt'
        export_text = export.ExportFile(input_file) ### create a new input file
        export_text.write(headers) ### Header is the same for each file
        for key in key_db:
            t = key_db[key]
            if offset > 0:
                t = [t[0]]+t[1+offset:]
            export_text.write(string.join(t,'\t')+'\n') ### Write z-score values and row names
        export_text.close()
        print 'File written...'
    return input_file

def importHopachOutput(filename):
    #print filename
    """ Import the ID order information """
    db={} ### Used to store the cluster data
    hopach_clusters=[]
    cluster_level=[]
    cluster_level2=[]
    cluster_level3=[]
    hopach_db={}
    cluster_db={}
    level2_level1={}
    firstLine = True
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if firstLine: firstLine = False
        else:
            t = string.split(data,'\t')
            final_level_order = int(t[-1])
            index, uid, cluster_number, cluster_label, cluster_level_order, final_label, final_level_order = string.split(data,'\t')
            try: l2 = str(int(round(float(cluster_label),0)))[:2]
            except Exception: l2 = int(cluster_label[0])
            try: l3 = str(int(round(float(cluster_label),0)))[:3]
            except Exception: l3 = int(cluster_label[0])
            hopach_clusters.append((int(final_level_order),int(index)-1)) ### Need to order according to the original index, sorted by the clustered order
            cluster_level.append(int(cluster_label[0])) ### This is the root cluster number
            cluster_level2.append(l2) ### Additional cluster levels
            cluster_level3.append(l3)
            hopach_db[uid] = cluster_label
            level2_level1[l2] = int(cluster_label[0])
            level2_level1[l3] = int(cluster_label[0])
            try: cluster_db[int(float(cluster_label[0]))].append(uid)
            except Exception: cluster_db[int(cluster_label[0])] = [uid]
            try: cluster_db[l2].append(uid)
            except Exception: cluster_db[l2] = [uid]
            try: cluster_db[l3].append(uid)
            except Exception: cluster_db[l3] = [uid]
        
    split_cluster=[]
    if 'column' in fn:
        cluster_limit = 50 ### typically less columns than rows
    else:
        cluster_limit = 75
    for cluster in cluster_db:
        #print cluster,len(cluster_db[cluster]),(float(len(cluster_db[cluster]))/len(hopach_db))
        if len(cluster_db[cluster])>cluster_limit and (float(len(cluster_db[cluster]))/len(hopach_db))>0.2:
            #print cluster
            if cluster<10:
                split_cluster.append(cluster)
    import unique
    levels1 = unique.unique(cluster_level)
    already_split={}
    updated_indexes={}
    if len(split_cluster)>0:
        print 'Splitting large hopach clusters:',split_cluster
        i=0
        for l2 in cluster_level2:
            l1 = level2_level1[l2]
            if l1 in split_cluster:
                cluster_level[i] = l2
                try:
                    l2_db = already_split[l1]
                    l2_db[l2]=[]
                except Exception: already_split[l1] = {l2:[]}
            i+=1
        ### Check and see if the l1 was split or not (might need 3 levels)
        i=0
        for l3 in cluster_level3:
            l1 = level2_level1[l3]
            if l1 in already_split:
                #l1_members = len(cluster_db[l1])
                l2_members = len(already_split[l1])
                #print l1, l3, l1_members, l2_members
                if l2_members == 1: ### Thus, not split
                    cluster_level[i] = l3
                    #print l1, l3, 'split'
            i+=1
    else:
        if len(cluster_level) > 50: ### Decide to use different hopach levels
            if len(levels1)<3:
                cluster_level = cluster_level2
        if len(cluster_level) > 200:
            if len(levels1)<4:
                cluster_level = cluster_level2
                
    hopach_clusters.sort()
    hopach_clusters = map(lambda x: x[1], hopach_clusters) ### Store the original file indexes in order based the cluster final order
    
    ### Change the cluster_levels from non-integers to integers for ICGS comparison group simplicity and better coloring of the color bar
    cluster_level2 = []
    ### Rename the sorted cluster IDs as integers
    cluster_level_sort = []
    for i in cluster_level:
        if str(i) not in cluster_level_sort:
            cluster_level_sort.append(str(i))
        cluster_level2.append(str(i))
    cluster_level_sort.sort()
    cluster_level = cluster_level2
    cluster_level2=[]
    i=1; cluster_conversion={}
    for c in cluster_level_sort:
        cluster_conversion[str(c)] = str(i)
        i+=1
        
    for c in cluster_level:
        cluster_level2.append(cluster_conversion[c])

    #print string.join(map(str,cluster_level2),'\t');sys.exit()
    db['leaves'] = hopach_clusters ### This mimics Scipy's cluster output data structure
    db['level'] = cluster_level2
    return db
            
class RScripts:
    def __init__(self,file):
        self._file = file
    def format_value_for_R(self,value):
        value = '"'+value+'"'
        return value
    def File(self):
        filename = self._file
        filename_list = string.split(filename,'/')
        filename = filename_list[-1]
        filename = self.format_value_for_R(filename)
        #root_dir = string.join(filename_list[:-1],'/')
        return filename
    
    def Monocle(self,samplelogfile,expPercent,p_val,numGroups):
        #samplelogfile='C:/Users/venz6v/Documents/Altanalyze R/data.txt' 
        #grp_list="C:/Users/venz6v/Documents/Altanalyze R/grous.txt"
        #gene_list="C:/Users/venz6v/Documents/Altanalyze R/gene.txt"
        filename=self.File()
        samplelogfile=findParentDir(filename)+'Monocle/expressionFile.txt"'
        grp_list=findParentDir(filename)+'Monocle/sampleGroups.txt"'
        gene_list=findParentDir(filename)+'Monocle/geneAnnotations.txt"'
        pseudo_tree=findParentDir(filename)+'Monocle/monoclePseudotime.pdf"'
        pseudo_txt=findParentDir(filename)+'Monocle/monoclePseudotime.txt"'
        #try: os.mkdir(findParentDir(samplelogfile)) ### create "hopach" dir if not present
        #except Exception: None
        #try: os.mkdir(findParentDir(grp_list)) ### create "hopach" dir if not present
        #except Exception: None
        #try: os.mkdir(findParentDir(gene_list)) ### create "hopach" dir if not present
        #except Exception: None
        #self._file = samplelogfile
        #samplelogfile = self.File()
        #self._file = grp_list
        #grp_list = self.File()
        #self._file = gene_list
        #gene_list = self.File()
        
        print 'Loading monocle package in R'
        print_out = r('library("monocle")')
        if "Error" in print_out:
            print 'Installing the R package "monocle" in Config/R'
            print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("monocle")')
            print print_out
            print_out = r('library("monocle")')
        if "Error" in print_out: print 'unable to download the package "monocle"';  
        print_out = r('library("monocle")')
        print "Reading Monocle data..."
        data_import = 'fpkm_matrix<-read.delim(%s,row.names=1,check.names=FALSE)' % samplelogfile
        #print [data_import]
        print_out = r(data_import);
        print print_out
    
        data_import = 'sample_sheet<-read.delim(%s,row.names=1,check.names=FALSE)' % grp_list
        #print [data_import]
        print_out = r(data_import);
        print print_out
        data_import = 'gene_ann<-read.delim(%s,row.names=1,check.names=FALSE)' % gene_list
        #print [data_import]
        print_out = r(data_import);
        print print_out
        print_out= r('pd <- new("AnnotatedDataFrame",data=sample_sheet)');
        print_out=r('fd <- new("AnnotatedDataFrame",data=gene_ann)');
        print_out=r('URMM <- newCellDataSet(as.matrix(fpkm_matrix),phenoData = pd,featureData =fd)');
        print print_out
        #colname(a) == colname(b)
        print_out=r('URMM<- detectGenes(URMM, min_expr = 0)')
        gene_exp='expressed_genes <- row.names(subset(fData(URMM), num_cells_expressed >=%s ))'% expPercent
        #print [gene_exp]
        try:print_out = r(gene_exp)
        except Exception:
                print "expression genes"
        print_out=r('length(expressed_genes)')
        print print_out

        # specify the grouping column for finding differential genes
        import multiprocessing
        cores = multiprocessing.cpu_count()
        print 'using', cores, 'cores'
        k =  'diff_test_res <- differentialGeneTest(URMM[expressed_genes, ], fullModelFormulaStr = "expression~Group",cores=%s)' % cores
        print [k]
        print_out=r(k)
        print print_out
        gene_ord='ordering_genes <- row.names(subset(diff_test_res, pval < %s))' %p_val
       
        print_out=r(gene_ord); print print_out
        print_out=r('write.table(ordering_genes,file="ordering_genes.txt")')  ### Writing out the informative genes used
        print print_out
        print_out=r('length(ordering_genes)'); print 'number or ordering genes',print_out
        
        print_out=r('ordering_genes <- intersect(ordering_genes, expressed_genes)'); print print_out
        print_out=r('URMM <- setOrderingFilter(URMM, ordering_genes)'); print print_out
        print_out=r('URMM <- reduceDimension(URMM, use_irlba = F)'); print print_out
        for i in range(numGroups,1,-1):
            span='URMM <- orderCells(URMM, num_paths = %s, reverse = F)'% i;
            
            print_out=r(span);
            print print_out
            if "Error" in print_out:
              continue
            else:
                print_out=r(span);print i
                print print_out
                break
        
        print_out=r('png("Monocle/monoclePseudotime.png")');
        print print_out
        print_out=r('plot_spanning_tree(URMM)'); print print_out
        print_out=r('dev.off()')
        print_out=r('pdf("Monocle/monoclePseudotime.pdf")');
        print print_out
        print_out=r('plot_spanning_tree(URMM)'); print print_out
        print_out=r('dev.off()')
        """
        print_out=r('pdf("Monocle/monoclePseudotimeOriginalGroups.pdf")');
        print print_out
        print_out=r('plot_spanning_tree(URMM), color_by = "originalGroups"'); print print_out
        print_out=r('dev.off()')
        """
        print_out=r('write.table(pData(URMM),file="Monocle/monoclePseudotime.txt")') 
        print " completed"
    
    def AffyNormalization(self,normalization_method,probe_level,batch_effects):
        print 'Loading affy package in R'
        print_out = r('library("affy")')
        if "Error" in print_out:
            #print_out = r('install.packages("ggplot2", repos="http://cran.us.r-project.org")')
            print 'Installing the R package "affy" in Config/R'
            print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("affy")')
            if "Error" in print_out: print 'unable to download the package "affy"'; forceError 
            print_out = r('library("affy")')
        if 'gcrma' in normalization_method:
            print 'Loading gcrma package in R'
            print_out = r('library("gcrma")')
            if "Error" in print_out:
                print 'Installing the R package "gcrma" in Config/R'
                print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("gcrma")')
                if "Error" in print_out: print 'unable to download the package "gcrma"'; forceError 
                print_out = r('library("gcrma")')
        if batch_effects == 'remove':
            ### Import or download support for SVA/Combat
            print 'Loading sva package in R'
            print_out = r('library("sva")')
            if "Error" in print_out:
                print 'Installing the R package "sva" in Config/R'
                print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("sva")')
                if "Error" in print_out: print 'unable to download the package "sva"'; forceError 
                print_out = r('library("sva")')
            
        print "Reading Affy files..."
        print_out = r('rawdata<-ReadAffy()')
        print print_out
        r('setwd("ExpressionInput")')
        if probe_level: ### normalize at the level of probes rahter than probeset (e.g., alt.exon analysis of 3' array)
            print_out = r('PM<-probes(rawdata,which="pm")'); print print_out
            print_out = r('AffyInfo<-dimnames(PM)[[1]]'); print print_out
            print_out = r('cutpos<-regexpr("\\d+$",AffyInfo,perl=T)'); print print_out
            print_out = r('AffyID<-substr(AffyInfo,1,cutpos-1)'); print print_out
            print_out = r('probe<-as.numeric(substr(AffyInfo,cutpos,nchar(AffyInfo)))'); print print_out
            print_out = r('data.bgc<-bg.correct(rawdata,method="rma")'); print print_out
            print_out = r('data.bgc.q<-normalize.AffyBatch.quantiles(data.bgc,type="pmonly")'); print print_out
            print_out = r('pm.bgc.q<-probes(data.bgc.q,which="pm")'); print print_out
            print_out = r('normalized<-cbind(AffyID,probe,pm.bgc.q)'); print print_out
            command = 'write.table(normalized,file='+self.File()+',sep="\t",row.names=FALSE, quote=FALSE)'
            print_out = r(command)
            print print_out
            print 'probe-level normalization complete'
        else:
            print "Begining %s normalization (will install array annotations if needed)... be patient" % normalization_method
            print_out = r('normalized<-%s(rawdata)') % normalization_method
            print print_out
            command = 'write.exprs(normalized,'+self.File()+')'; print_out = r(command)
            print print_out
        print self.File(), 'written...'
        if batch_effects == 'remove':
            ### Import data
            command = 'mod = model.matrix(~as.factor(cancer) + age, data=pheno)'
            print_out = r(command)
            command = 'cdata = ComBat(dat=normalized, batch=as.factor(pheno$batch), mod=mod, numCov=match("age", colnames(mod)))'
            print_out = r(command)
            command = 'write.table(cdata,file='+self.File()+',sep="\t",row.names=FALSE, quote=FALSE)'
            print_out = r(command)
        output_file = string.replace(self.File(),'exp.','stats.')
        print_out = r('calls<-mas5calls(rawdata)')
        #print_out = r('pvals<-se.exprs(calls)') ### outdated?
        print_out = r('pvals<-assayData(calls)[["se.exprs"]]')
        command = 'write.table(pvals,'+output_file+',sep = "\t", col.names = NA)'; print_out = r(command)
        print output_file, 'written...'

    def Limma(self,test_type):
        r('library("limma")')
        filename = self.File()
        try: output_file = string.replace(filename,'input','output-'+test_type)
        except ValueError: output_file = filename[0:-4]+'-output.txt'
        print "Begining to process",filename
        data_import = 'data<-read.table(%s,sep="\t",header=T,row.names=1,as.is=T)' % filename
        print_out = r(data_import)
        design_matrix_file = string.replace(filename,'input','design')
        design_import = 'design<-read.table(%s,sep="\t",header=T,row.names=1,as.is=T)' % design_matrix_file
        design_matrix = r(design_import)
        print_out = r('fit<-lmFit(data,design)')
        fit_data = r['fit']
        print_out = r('fit<-eBayes(fit)')
        fit_data = r['fit']
        contrast_matrix_file = string.replace(filename,'input','contrast')
        contrast_import = 'contrast<-read.table(%s,sep="\t",header=T,row.names=1,as.is=T)'  % contrast_matrix_file
        print_out = r(contrast_import)
        contrast_matrix = r['contrast']
        r('contrast<-as.matrix(contrast)')
        r('fit.contrast<-contrasts.fit(fit,contrast)')
        r('fit.contrast<-eBayes(fit.contrast)')
        r('nonadj<-fit.contrast$F.p.value')
        if test_type == 'fdr':
            print_out = r('results<-p.adjust(fit.contrast$F.p.value,method="fdr")')
        else:
            print_out = r('results<-nonadj')
        result = r['results']
        print 'test_type=',test_type
        print_out = r('sum(results<0.05)')
        summary = r['sum']
        print "Number of probeset with a p<0.05",summary,"using",test_type
        r('output<-cbind(data,results)')
        output = 'write.table(output,%s,sep="\t")' % output_file
        print_out = r(output)
        print output_file, 'written...'
        
    def Multtest(self,test_type):
        r('library("multtest")')
        filename = self.File()
        try: output_file = string.replace(filename,'input','output')
        except ValueError: output_file = filename[0:-4]+'-output.txt'
        print "Begining to process",filename
        parse_line = 'job<-read.table(%s,sep="\t", row.names=1, as.is=T)' % filename
        print_out = r(parse_line)
        print_out = r('matrix_size<-dim(job)')
        print_out = r('label<-job[1,2:matrix_size[2]]')
        print_out = r('jobdata<-job[2:matrix_size[1],2:matrix_size[2]]')
        if test_type == "f":
            print_out = r('ttest<-mt.maxT(jobdata,label, test="f", B=50000)')
        if test_type == "t":
            print_out = r('ttest<-mt.maxT(jobdata,label)')
        print_out = r('ttest2<-ttest[order(ttest[,1]),]')
        write_file = 'write.table(ttest2,%s,sep="\t")' % output_file
        print_out = r(write_file)
        print "Results written to:",output_file 
    def check_hopach_file_type(self):
        if 'hopach.input' in self.File():
            return 'continue'
        else: return 'break'
    def check_multtest_file_type(self):
        if 'output' not in self.File():
            return 'continue'
        else: return 'break'
    def check_limma_file_type(self):
        if 'input' in self.File():
            return 'continue'
        else: return 'break'
    def Hopach(self,cluster_method,metric_gene,force_gene,metric_array,force_array):
        if R_present==False:
            rNotPresent
        print_out = r('library("Biobase")')
        if "Error" in print_out:
            print 'Installing the R package "Biobase" in Config/R'
            print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("Biobase")')
            if "Error" in print_out: print 'unable to download the package "Biobase"'; forceError 
            print_out = r('library("Biobase")')
        print_out = r('library("hopach")')
        if "Error" in print_out:
            print 'Installing the R package "hopach" in Config/R'
            print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("hopach")')
            if "Error" in print_out: print 'unable to download the package "hopach"'; forceError 
            print_out = r('library("hopach")')
        filename = self.File()
        #r('memory.limit(2000)')
        print "Begining to process",filename,"using HOPACH"
        metric_g = self.format_value_for_R(metric_gene)
        metric_a = self.format_value_for_R(metric_array)
        parse_line = 'data<-read.table(%s,sep="\t",as.is=T,row.names=1,header=T)' % filename
        checklinelengths(self._file)
        print_out = r(parse_line)
        dat = r['data']
        print_out = r(parse_line)
        #print "Number of columns in input file:",len(dat)
        print_out = r('data<-as.matrix(data)')
        dat = r['data']
        #print "Number of columns in matrix:",len(dat)
        force1=''; force2=''; hopg='NULL'; hopa='NULL'; distmatg='NULL'; distmata = 'NULL' ### defaults for tree export
        if force_gene != '' and force_gene != 0: force1=',kmax='+str(force_gene)+', khigh='+str(force_gene)+', K='+str(force_array)
        if force_array != '' and force_array != 0: force2=',kmax='+str(force_array)+', khigh='+str(force_array)+', K='+str(force_array)
        if cluster_method == 'both' or cluster_method == 'gene':
            distance_matrix_line = 'distmatg<-distancematrix(data,d=%s)' % metric_g
            #print distance_matrix_line
            if len(dat) > 1:
                print_out1 = r(distance_matrix_line)
                print_out2 = r('hopg<-hopach(data,dmat=distmatg,ord="own"'+force1+')')
                #print 'hopg<-hopach(data,dmat=distmatg,ord="own"'+force1+')'
                try: hopach_run = r['hopg']
                except Exception:
                    print print_out1
                    print print_out2
                hopg = 'hopg'
                distmatg = 'distmatg'
                gene_output = self.HopachGeneOutputFilename(metric_gene,str(force_gene))
                output = 'out<-makeoutput(data,hopg,file=%s)' % gene_output
                #print output
                print_out = r(output)
                #print print_out
                output_file = r['out']
                status = 'stop'
                if 'clustering' in hopach_run:
                    if 'order' in hopach_run['clustering']:
                        try:
                            if len(hopach_run['clustering']['order']) > 10: status = 'continue'
                        except TypeError:
                            error = 'file: '+filename+": Hopach returned the array of cluster orders as blank while clustering GENES... can not process cluster... continuing with other files"
                            print error; errors.append(error)
                        if status == 'continue':
                            r(output_file); print 'hopach output written'
            else:
                error = 'file: '+filename+" Hopach returned data-matrix length zero...ARRAY clusters can not be generated"
                print error; errors.append(error)
        if cluster_method == 'both' or cluster_method == 'array':
            distance_matrix_line = 'distmata<-distancematrix(t(data),d=%s)' % metric_a
            if len(dat) > 1:
                dist = r(distance_matrix_line)
                #print distance_matrix_line
                print_out = r('hopa<-hopach(t(data),dmat=distmata,ord="own"'+force1+')') #,coll="all"
                #print ['hopa<-hopach(t(data),dmat=distmata,ord="own",'+force2+')']
                #print 'hopa<-hopach(t(data),dmat=distmata,ord="own"'+force2+')'
                hopach_run = r['hopa']
                hopa = 'hopa'
                distmata = 'distmata'
                array_output = self.HopachArrayOutputFilename(metric_array,str(force_array))
                output = 'out<-makeoutput(t(data),hopa,file=%s)' % array_output
                print_out = r(output)
                output_file = r['out']
                status = 'stop'
                if 'clustering' in hopach_run:
                    if 'order' in hopach_run['clustering']:
                        try:
                            if len(hopach_run['clustering']['order']) > 10: status = 'continue'
                        except TypeError:
                            error = 'file: '+filename+": Hopach returned the array of cluster orders as blank while clustering ARRAYS... can not process cluster"
                            print error; errors.append(error)
                        if status == 'continue':
                            r(output_file); print 'hopach output written'
            else:
                error = 'file: '+filename+"data-matrix length zero...ARRAY clusters can not be generated...continuing analysis"
                print error; errors.append(error)
        if len(metric_g)==0: metric_g = 'NULL'
        if len(metric_a)==0: metric_a = 'NULL'
        try:
            output_filename = string.replace(gene_output,'rows.','')
            cdt_output_line = 'hopach2tree(data, file = %s, hopach.genes = %s, hopach.arrays = %s, dist.genes = %s, dist.arrays = %s, d.genes = %s, d.arrays = %s, gene.wts = NULL, array.wts = NULL, gene.names = NULL)' % (output_filename,hopg,hopa,distmatg,distmata,metric_g,metric_a) ###7 values
        except Exception: None
        make_tree_line = 'makeTree(labels, ord, medoids, dist, side = "GENE")' ### Used internally by HOPACH
        #print cdt_output_line
        try: print_out = r(cdt_output_line)
        except Exception: None
        #print print_out
        
    def HopachGeneOutputFilename(self,value,force):
        filename = self.File() ### Relative to the set working directory
        if 'hopach.input' in filename: ### When running this module on it's own (requires nown filetypes)
            new_filename = string.replace(filename,'hopach.input','hopach.output')
            if len(value)>1: new_filename = string.replace(new_filename,'.txt','-'+value+'.txt')
            if len(force)>0: new_filename = string.replace(new_filename,'.txt','-'+'force_'+str(force)+'c.txt')
        else: ### When called from an external heatmap visualization module
            filename = self._file ### full path
            new_filename = findParentDir(filename)+'/hopach/rows.'+findFileName(filename)
            try: os.mkdir(findParentDir(new_filename)) ### create "hopach" dir if not present
            except Exception: None
            new_filename = '"'+new_filename+'"'
        return new_filename
    
    def HopachArrayOutputFilename(self,value,force):
        filename = self.File()
        if 'hopach.input' in filename: ### When running this module on it's own (requires nown filetypes)
            new_filename = string.replace(filename,'hopach.input','arrays.output')
            if len(value)>1: new_filename = string.replace(new_filename,'.txt','-'+value+'.txt')
            if len(force)>0: new_filename = string.replace(new_filename,'.txt','-'+'force_'+str(force)+'c.txt')
        else:
            filename = self._file ### full path
            filename = self._file ### full path
            new_filename = findParentDir(filename)+'/hopach/columns.'+findFileName(filename)
            try: os.mkdir(findParentDir(new_filename)) ### create "hopach" dir if not present
            except Exception: None
            new_filename = '"'+new_filename+'"'
        return new_filename
    
    def display(self):
        print self.data
    
class FormatData:
    def setdata(self,value):
        self.data = value
    def transform(self):
        self.data = checktype(self.data)
    def display(self):
        print self.data
    def returndata(self):
        return self.data
    
def checktype(object):
    ###Checks to see if item is a list or dictionary. If dictionary, convert to list
    import types
    if type(object) is types.DictType:
        object = converttolist(object)
    elif type(object) is types.ListType:
        object = object
    elif type(object) is types.TupleType:
        object = list(object)
    elif type(object) is types.StringType:
        object = importtable(object)
    return object

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def checklinelengths(filename):
    fn=filepath(filename); first_row='yes'; line_number=0
    for line in open(fn,'rU').xreadlines():
        try: data = cleanUpLine(line)
        except Exception: print 'error parsing the line:',[line], line_number
        t = string.split(data,'\t')
        if first_row == 'yes':
            elements = len(t)
            first_row = 'no'
        else:
            if len(t) != elements:
                print "Line number", line_number, "contains",len(t),"elements, when",elements,"expected...kill program"
                print filename; kill
        line_number+=1
        
def converttolist(dictionary):
    ###Converts dictionary to list by appending the dictionary key as the first item in the list
    converted_lists=[]
    for key in dictionary:
        dictionary_list = dictionary[key]
        dictionary_list.reverse(); dictionary_list.append(key); dictionary_list.reverse()
        converted_lists.append(dictionary_list)
    return converted_lists

############ IMPORT FILES BEGIN ############
def importtable(filename):
    fn=filepath(filename); tab_db = []
    for line in open(fn,'rU').readlines():
        data,null = string.split(line,'\n')
        t = string.split(data,'\t')
        tab_db.append(t)
    return tab_db

def filepath(filename):
    dir=os.path.dirname(__file__)  #directory file is input as a variable
    status = verifyFile(filename)
    if status:
        fn = filename
    else:
        fn=os.path.join(dir,filename)
    return fn

def verifyFile(filename):
    status = False
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = True;break
    except Exception: status = False
    return status

def findFileName(filename):
    filename = string.replace(filename,'\\','/')
    dataset_name = string.split(filename,'/')[-1]
    return dataset_name

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]
    
def setWorkingDirectory(filename):
    ### Set R's working directory when calling this module remotely
    working_dir = findParentDir(filename)
    setwd = 'setwd("%s")' % working_dir
    try: r(setwd)
    except Exception:
        print [filename]
        print [working_dir]
        print traceback.format_exc()
        kill
        
def read_directory(sub_dir):
    dir=os.path.dirname(__file__) 
    #print "Working Directory:", r('getwd()')
    working_dir = dir+'/'+sub_dir[1:]
    setwd = 'setwd("%s")' % working_dir
    r(setwd)
    #print "Working Directory:", r('getwd()')
    
    dir_list = os.listdir(dir +'/'+ sub_dir[1:]); dir_list2 = []
    for entry in dir_list: #add in code to prevent folder names from being included
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2
    
def CreateFilesMonocle(filename,rawExpressionFile,species='Hs'):
    first_row = True
    key_db={}
    key_list=[]
    fn=filepath(filename)
    offset=0
    nonNumericsPresent=False
    try:
        import gene_associations
        gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    except Exception:
        print "gene_symbols present"
        gene_to_symbol={}
    setWorkingDirectory(findParentDir(filename)[:-1])
    try: os.mkdir(findParentDir(filename)+'/Monocle')
    except Exception: None     
    #filename=self.File() 
    x = 0
    data_name=findParentDir(filename)+'/Monocle/expressionFile.txt'
    gene_name=findParentDir(filename)+'/Monocle/geneAnnotations.txt'
    sample_name=findParentDir(filename)+'/Monocle/sampleGroups.txt'
    gene_names = [];
    gene_list=[];
    dat=[];
    export_cdt = open(sample_name,'w')
    export_gene=open(gene_name,'w')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if first_row == True:
            if 'row_clusters-flat' in t and 'row_clusters-flat' not in t[0]:
                headers = string.join(t[2:],'\t')+'\n'
                offset = 1
            else:
                headers = string.join(t[1:],'\t')+'\n'
                
            first_row = False
        else:
            key = t[0]
            if key!='column_clusters-flat':
                key_list.append(key)
                
                try: s = map(float,t[offset+1:])
                except Exception:
                    nonNumericsPresent=True
                key_db[key]=t
            else:
                clusters = map(str,t[offset+1:])
    for key in key_list:
            t = key_db[key]
            s=[key]
            if offset ==1: s.append('')
            temp=[]
            for value in t[offset+1:]:
                try: temp.append(float(value))
                except Exception: pass
            min1=min(temp)
        
            for value in t[offset+1:]:
                try: s.append(str(float(value)-min1))
                except Exception: s.append('0.000101')
            key_db[key]=s
    export_object = open(data_name,'w') 
    export_object.write(''+'\t'+headers) ### Header is the same for each file
    for key in key_list:
            t = key_db[key]
            if offset > 0:
                t = [t[0]]+t[1+offset:]
            export_object.write(string.join(t,'\t')+'\n') ### Write z-score values and row names
    export_object.close()
    print 'File written...'
    #return input_file
    
    array_names = []; array_linker_db = {}; d = 0; i = 0
    for entry in headers.split('\t'):
        entry=cleanUpLine(entry)
        if '::' in entry:
            a = (entry.split("::"))
        elif ':' in entry:
            a = (entry.split(":"))
        else:
            a = (clusters[i],entry)
        #entry=string.join(a,'.')
              
        ent=entry+'\t'+a[0];
        #if(ent[0].isdigit()):
        #    ent='X'+ent[0:]     
        #if '-' in ent:
            #   ent=string.replace(ent,'-','.')
        #if '+' in ent:
            #   ent=string.replace(ent,'+','.')
            #print j
        array_names.append(ent);
        i+=1
        
    i=0
    eheader = string.join(['']+['Group'],'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eheader)
    for row in array_names:
        export_cdt.write(row+'\n')
        i+=1
    export_cdt.close()
    gheader = string.join(['']+ ['gene_short_name'],'\t')+'\n' ### format column-flat-clusters for export
    export_gene.write(gheader)
    for key in key_list:
        proceed=False
        ### The commented out code just introduces errors and is not needed - re-evaluate in the future if needed
        """
        if key in gene_to_symbol:
            symbol = gene_to_symbol[key][0]
            if symbol in gene_list:
                nid = symbol
                proceed = True
            if proceed:
                k=gene_list.index(nid)
                export_object.write(line)
                export_gene.write(key+'\n')
        else:
            export_gene.write(key+'\t'+key+'\n')"""
        export_gene.write(key+'\t'+key+'\n')
 
            
    export_object.close() 
    export_gene.close()

def reformatHeatmapFile(input_file):
    import unique
    export_file=string.replace(input_file,'Clustering-','Input-')
    eo = export.ExportFile(export_file)
    first_row = True
    fn=filepath(input_file)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if first_row == True:
            if 'column_clusters-flat' not in t:
                array_names = []
                for i in t[2:]:
                    array_names.append(string.replace(i,':','-'))
                    #print array_names;sys.exit()
                    #array_names.append(i)
            elif 'column_clusters-flat' in t:
                array_clusters = t[2:]
                unique_clusters = unique.unique(array_clusters)
                ind=0; headers=[]
                for c in array_clusters:
                    headers.append(c+'::'+array_names[ind])
                    ind+=1
                headers = string.join(['uid']+headers,'\t')+'\n'
                eo.write(headers)
                first_row = False
        else:
            values = string.join([t[0]]+t[2:],'\t')+'\n'
            eo.write(values)
    return export_file, len(unique_clusters)

def run_JTKcycle(expFile,annotFile,Time_range1, Time_range2,No_of_Timepoints,No_of_replicates,timepoint_difference):
    
    print 'Loading JTK-Cycle package in R'
    path='"'+r_package_path+'/JTK_CYCLE.R"'
    #print [path]
    line = 'source(%s)' %  path
    print_out = r(line)
    """
    if "Error" in print_out:
        print 'Installing the R package "JTK_CYCLE.R" in Config/R'
        print_out = r('install.packages("devtools")')
        
        print print_out
        print_out = r('library(devtools)')
        print print_out
        print_out = r('install_github("mfcovington/jtk-cycle")')
        #print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("jtk-cycle")')
        print print_out
        print_out = r('library("JTK_CYCLE.R")')
    sys,exit()
    print_out = r('source("/Users/ram5ge/Desktop/Krithika/JTK_Cycle/JTK_CYCLE.R")');print print_out
    if "Error" in print_out: print "JTK_CYCLE.R is missing"
    else: print 'Loading JTK Cycle'
    """
    
    print_out = r('project <- "JTK_output"')   
    print_out = r('options(stringsAsFactors=FALSE)');print print_out
  
    a = '"'+annotFile+'"'
    read_annot = 'annot <- read.delim(%s)' % a
    print [read_annot]
    print_out = r(read_annot);#print print_out
    v = '"'+expFile+'"'
    read_data = 'input_data <- read.delim(%s)' % v
    print [read_data]
    print_out = r(read_data);#print print_out
    
    print_out = r('rownames(input_data) <- input_data[,1]');#print print_out
    print_out = r('input_data <- input_data[,-1]');#print print_out
    #dist_calc = r('jtkdist(24,1)')
    
    dist_calc = 'jtkdist(%s,%s)' % (str(No_of_Timepoints), str(No_of_replicates))
    print [dist_calc]
    print_out = r(dist_calc);#print print_out
    
    period_calc = 'periods <- %s:%s' %(str(Time_range1), str(Time_range2))
    print [period_calc]

    print_out = r(period_calc);#print print_out
    j = str(timepoint_difference)
    jtk_calc = 'jtk.init(periods,%s)' % j
    print [jtk_calc]
    print_out = r(jtk_calc);#print print_out
    v = 'cat("JTK analysis started on",date(),"\n")'
    print [v]
    print_out = r(v);#print print_out
    print_out = r('flush.console()');#print print_out

    v = 'st <- system.time({res <- apply(data,1,function(z)'
    v+= ' {jtkx(z); c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)});'
    v+= ' res <- as.data.frame(t(res)); bhq <- p.adjust(unlist(res[,1]),"BH");'
    v+= ' res <- cbind(bhq,res); colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP");'
    v+= ' results <- cbind(annot,res,data); results <- results[order(res$ADJ.P,-res$AMP),]})'
    print [v]
    print_out = r(v); print print_out

    #print_out = r('dim(X)');print print_out

    print_out = r('print(st)');print #print_out
    
    print_out = r('save(results,file=paste("JTK",project,"rda",sep="."))');#print print_out
    print_out = r('write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")');#print print_out
    
    
def performMonocleAnalysisFromHeatmap(species,heatmap_output_dir,rawExpressionFile):
    numGroups=10
    if 'Clustering-' in heatmap_output_dir:
        export_file,numGroups = reformatHeatmapFile(heatmap_output_dir)
    #else:
    export_file = heatmap_output_dir;
    CreateFilesMonocle(export_file,rawExpressionFile,species=species)
    print 'Looking for',numGroups, 'Monocle groups in the input expression file.'
    remoteMonocle(export_file,expPercent=5,pval=0.05,numGroups=numGroups)
    
if __name__ == '__main__':
    expFile = '/Users/saljh8/Downloads/Liver_Smoothed_exp_steady_state.txt'
    annotFile = '/Users/saljh8/Downloads/Liver_annot.txt'
    Time_range1 = '10'
    Time_range2 = '12'
    No_of_Timepoints = '24'
    No_of_replicates = '1'
    timepoint_difference = '2'
    run_JTKcycle(expFile,annotFile,Time_range1, Time_range2,No_of_Timepoints,No_of_replicates,timepoint_difference);sys.exit()
 
    
    
    errors = []
    cluster_method='array';metric_gene="";force_gene='';metric_array="euclid";force_array=''
    analysis_method='hopach'; multtest_type = 'f'
    #Sample log File
    #Input-exp.MixedEffectsThanneer-DPF3%20DMRT3%20FOXA1%20SMAD6%20TBX3%20amplify%20monocle-hierarchical_cosine_correlated.txt
    filename='/Users/saljh8/Desktop/cardiacRNASeq/DataPlots/Clustering-additionalExpressionSingleCell-annotated-hierarchical_cosine_cosine2.txt'
    rawExpressionFile = filename
    #filename = "/Volumes/SEQ-DATA/Eric/embryonic_singlecell_kidney/ExpressionOutput/Clustering/SampleLogFolds-Kidney.txt"
    #filename = "/Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/NewCardiacMarkers1/FullDataset/ExpressionOutput/Clustering/SampleLogFolds-CM.txt"
    #rawExpressionFile = '/Volumes/SEQ-DATA/SingleCell-Churko/Filtered/Unsupervised-AllExons/NewCardiacMarkers1/FullDataset/ExpressionInput/exp.CM-steady-state.txt'
    #filename = '/Users/saljh8/Desktop/Stanford/ExpressionInput/amplify/DataPlots/Clustering-exp.EB-SingleCell-GPCR-hierarchical_cosine_correlation.txt'
    #rawExpressionFile = '/Users/saljh8/Desktop/Stanford/ExpressionInput/exp.EB-SingleCell.txt'
    performMonocleAnalysisFromHeatmap('Hs',filename,rawExpressionFile);sys.exit()
    CreateFilesMonocle(filename,rawExpressionFile)
    remoteMonocle(filename,expPercent=0,pval=0.01,numGroups=5);sys.exit()
    
    filename = '/Users/nsalomonis/Downloads/GSE9440_RAW/ExpressionInput/exp.differentiation.txt'
    remoteAffyNormalization(filename,'rma',True,'remove'); sys.exit()
    
    print "******Analysis Method*******"
    print "Options:"
    print "1) Multtest (permutation ftest/ttest)"
    print "2) HOPACH clustering"
    print "3) limma 2-way ANOVA"
    
    inp = sys.stdin.readline(); inp = inp.strip()
    analysis_method_val = int(inp)
    if analysis_method_val == 1: analysis_method = "multtest"
    if analysis_method_val == 2: analysis_method = "hopach"
    if analysis_method_val == 3: analysis_method = "limma"


    if analysis_method == "hopach":
        print "******Analysis Options*******"
        print "Cluster type:"
        print "1) genes only (cluster rows)"
        print "2) arrays only (cluster columns)"
        print "3) both"
      
        inp = sys.stdin.readline(); inp = inp.strip()
        cluster_type_call = int(inp)
        if cluster_type_call == 1: cluster_method = "gene"
        if cluster_type_call == 2: cluster_method = "array"
        if cluster_type_call == 3: cluster_method = "both"

        if cluster_method == "array" or cluster_method == "both":        
            print "******Analysis Options For Array Clustering*******"
            print "Cluster metrics:"
            print "1) euclidian distance (sensitive to magnitude)"
            print "2) cosine angle distance (not sensitive to magnitude)"
            print "3) correlation distance"
            inp = sys.stdin.readline(); inp = inp.strip()
        if cluster_method == "array" or cluster_method == "both":
            metric_array_call = int(inp)
            if metric_array_call == 1: metric_array = "euclid"
            if metric_array_call == 2: metric_array = "cosangle"
            if metric_array_call == 3: metric_array = "cor"
            
        if cluster_method == "gene" or cluster_method == "both":        
            print "******Analysis Options For Gene Clustering*******"
            print "Cluster metrics:"
            print "1) euclidian distance (sensitive to magnitude)"
            print "2) cosine angle distance (not sensitive to magnitude)"
            print "3) correlation distance"
            inp = sys.stdin.readline(); inp = inp.strip()
        if cluster_method == "gene" or cluster_method == "both":
            try: metric_gene_call = int(inp)
            except ValueError: print [inp], 'not a valid option'; sys.exit()
            if metric_gene_call == 1: metric_gene = "euclid"
            if metric_gene_call == 2: metric_gene = "cosangle"
            if metric_gene_call == 3: metric_gene = "cor"
            if metric_gene == "cosangle":
                print "******Analysis Options*******"
                print "Absolute Clustering:"
                print "1) yes"
                print "2) no"
                inp = sys.stdin.readline(); inp = inp.strip()
                if inp == "1": metric_gene = "abscosangle"
                
        print "Force Cluster Number for Arrays:"
        print "Enter 'n' if you don't want to "
        print "Enter number of clusters of arrays if you do"
      
        inp = sys.stdin.readline(); inp = inp.strip()
        if inp == 'n' or inp == 'N': force_array = ''
        else:force_array = int(inp)
        working_dir = '/hopach_input'
        
    if analysis_method == "multtest":
        print "******Analysis Options*******"
        print "Statistical test:"
        print "1) ftest (for multiple groups)"
        print "2) ttest (for two groups)"
      
        inp = sys.stdin.readline(); inp = inp.strip()
        multtest_type_call = int(inp)
        if multtest_type_call == 1: multtest_type = "f"
        if multtest_type_call == 2: multtest_type = "t"
        working_dir = '/multtest_input'
        
    if analysis_method == "limma":
        working_dir = '/limma_input'
        print "******Analysis Options*******"
        print "Statistical test:"
        print "1) Non-adjusted"
        print "2) FDR"
        
        inp = sys.stdin.readline(); inp = inp.strip()
        limma_type_call = int(inp)
        if limma_type_call == 1: limma_type = "nonadj"
        if limma_type_call == 2: limma_type = "fdr"

    dir_list = read_directory(working_dir)
    for input in dir_list:    #loop through each file in the directory to output results
        input_file = working_dir + "/"+ input
        input_file = input_file[1:] #not sure why, but the '\' needs to be there while reading initally but not while accessing the file late
        z = RScripts(input_file)
        if analysis_method == "hopach":
            status = z.check_hopach_file_type()
            if status == 'continue':
                z.Hopach(cluster_method,metric_gene,force_gene,metric_array,force_array)
        if analysis_method == "multtest":
            status = z.check_multtest_file_type()
            if status == 'continue':
                z.Multtest(multtest_type)
        if analysis_method == "limma":
            status = z.check_limma_file_type()
            if status == 'continue':
                design_matrix_file = string.replace(input,'input','design')
                contrast_matrix_file = string.replace(input,'input','contrast')
                if design_matrix_file in dir_list and contrast_matrix_file in dir_list:
                    z.Limma(limma_type)
    if analysis_method == "hopach":
        if len(errors)>0:
            print "**************ALL ERRORS**************"
            for entry in errors:
                print entry
        else: print 'Execution complete... check outputs for verification'
