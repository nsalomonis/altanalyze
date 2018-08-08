#!/usr/bin/env python
""" This script identifies correlations between splicing factors or other input gene lists and splicing events
in an pre-filtered PSI events file. The expression and PSI file must have the same sample order as input. """
    
import os, sys, string, getopt
import math, numpy
from scipy import stats
import warnings
import time

def importFile(filename,convertToFloat=False):
    db={}
    firstRow=True
    dataset_max=0
    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        uid = t[0]
        if firstRow: samples = t[1:]; firstRow=False
        if len(t)>1:
            values = t[1:]
            if convertToFloat:
                try:
                    values = map(float,values)
                    max_val = max(values)
                    if max_val > dataset_max:
                        dataset_max = max_val
                except:
                    continue ### Header row
            else:
                values = t[1] ### Just store the gene symbol
        else:
            values = uid
        db[uid]=values
    if convertToFloat:
        if dataset_max>100: ### Data is not log2
            print 'Converting gene expression data to log2 values'
            for uid in db:
                db[uid] = map(lambda x: math.log(x+1,2), db[uid])
        print 'Imported %d gene expression rows' % len(db)
        return db, samples
    else:
        print 'Imported %d splicing factors' % len(db) 
        return db
        
def importPSIData(PSI_dir,samples):
    """Import PSI data from either EventAnnotation or PSI value file"""
    firstRow=True
    PSI_data_db={}
    for line in open(PSI_dir,'rU').xreadlines():
        data = line.rstrip()
        PSI_data = string.split(data,'\t')
        if firstRow:
            data = string.replace(data,'.bed','')
            PSI_data = string.split(data,'\t')
            header_row = PSI_data
            if 'ProteinPredictions' in PSI_data:
                data_index = PSI_data.index('EventAnnotation')+1
                uid_index = PSI_data.index('UID')
            else:
                uid_index = 0
                data_index = 1
            psi_samples = PSI_data[data_index:]
            if psi_samples != samples:
                print 'Error: The gene expression sample order does not match the PSI. Exiting';sys.exit()
            else:
                print 'Confirmed: The sample order of the gene expression and splicing files match.'
            firstRow=False
        else:
            if len(PSI_data) != len(header_row):
                empty_offset = len(header_row)-len(PSI_data)
                PSI_data+=['']*empty_offset
            junctionID = PSI_data[uid_index]
            PSI_data = PSI_data[data_index:]
            try:
                values = map(lambda x: float(x), PSI_data)
            except Exception:
                values=[]
                for value in PSI_data:
                    try: values.append(float(value))
                    except:
                        values.append(0.000101) ### Missing value
                values = numpy.ma.masked_values(values,0.000101)
            PSI_data_db[junctionID]=values

    print 'Imported %d splicing event rows' % len(PSI_data_db)
    return PSI_data_db

def findcorrelations(SF_dir, PSI_dir, exp_dir, output_dir, PearsonCutoff):
    print ''
    ### Import the list of splicing factors or other genes of interest
    genesToExamine = importFile(SF_dir)
    ### Import the tab-delimited gene expression matrix
    geneExpression_db, samples = importFile(exp_dir,convertToFloat=True)
    ### Import the PSI data
    PSI_data_db = importPSIData(PSI_dir,samples)
    
    ### Create an export directory
    results_dir = output_dir+'/SFCorrelations_rho-'+str(PearsonCutoff)
    try: os.mkdir(results_dir)
    except: pass
    eo=open(output_dir+'/SFCorrelations/SF_correlations.txt','w')
    eo.write('Splicing Factor'+'\t'+'Events Count'+'\n')
    
    counter=0
    gene_correlation_time = []
    for gene in genesToExamine:
        gene_name = genesToExamine[gene]
        if gene in geneExpression_db:
            start_time = time.time()
            ### Hence, the gene is a splicing factor
            expression_values = geneExpression_db[gene]
            Corrsflist=[]
            count=0
            for junctionID in PSI_data_db:
                psi_values = PSI_data_db[junctionID]
                if 0.000101 in psi_values:
                    coefr=numpy.ma.corrcoef(expression_values,psi_values)
                    rho = coefr[0][1]
                else:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore",category=RuntimeWarning) 
                        rho,p = stats.pearsonr(expression_values,psi_values)
                if abs(rho)>PearsonCutoff:
                    count+=1
                    Corrsflist.append([junctionID,rho])
                gene_correlation_time.append(time.time()-start_time)
        
        eo.write(gene_name+'\t'+str(count)+'\n')
        filename=results_dir+"/"+gene_name+"_"+str(count)+".txt"
        if count>20:
            eg=open(filename,"w")
            eg.write("SplicingEvent\tSystemCode\tPearsonRho\n")
        
            for (junctionID,rho) in Corrsflist:
                eg.write(junctionID+"\t"+"Ae\t"+str(rho)+"\n")
            eg.close()
        counter+=1
        print '*',
    print '\n...Correlations obtained on average of %d seconds/gene' % numpy.mean(gene_correlation_time)

def performEventEnrichment(output_dir,eventDir,species):
    """Import significant splicing events from metaDataAnalysis.py comparisons and test for their
    statistical enrichmet relative to the Splicing Factor correlated events."""
    import collections
    import mappfinder
    event_db = collections.OrderedDict()
    import UI
    ### Import the splice-ICGS significant splicing events per signature
    files = UI.read_directory(eventDir)
    for file in files:
        if '.txt' in file and 'PSI.' in file:
            ls=[]
            event_db[file[:-4]]=ls ### This list is subsequently updated below
            fn = eventDir+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                if firstLine:
                    event_index = t.index('Event-Direction')
                    firstLine= False
                    continue
                uid = t[0]
                if 'U2AF1-like' in file:
                    if t[1] == "inclusion":
                        ls.append(uid) #ls.append((uid,t[event_index]))
                else:
                    ls.append(uid) #ls.append((uid,t[event_index]))
                    
    ### Import the splicing-factor correlated splicing events to identify associated signatures
    splicing_factor_correlated_scores={}
    gene_to_symbol=None
    files = UI.read_directory(output_dir)
    for file in files:
        if '.txt' in file and '_' in file:
            R_ls=[]
            if 'ENS' in file:
                splicing_factor = file[:-4]
                if gene_to_symbol==None: ### Import only once
                    import gene_associations
                    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
                sf = 'ENS'+string.split(splicing_factor,'ENS')[1]
                splicing_factor = string.split(sf,'_')[0]
                if splicing_factor in gene_to_symbol:
                    splicing_factor = gene_to_symbol[splicing_factor][0]
            else:
                splicing_factor = string.split(file[:-4],'_')[0]
            fn = output_dir+'/'+file
            firstLine = True
            for line in open(fn,'rU').xreadlines():
                data = line.rstrip()
                t = string.split(data,'\t')
                event = t[0]
                R_ls.append(event)
            R=len(R_ls)
            N=80000
            for signature in event_db:  
                n_ls=event_db[signature]
                n = len(n_ls)
                r_ls=set(R_ls).intersection(n_ls)
                r = len(r_ls)
                ### Calculate a Z-score
                try: z = Zscore(r,n,N,R)
                except ZeroDivisionError: z = 0.0000
                ### Calculate a Z-score assuming zero matching entries
                try: null_z = Zscore(0,n,N,R)
                except ZeroDivisionError: null_z = 0.000
                ### Calculate a Fischer's Exact P-value
                pval = mappfinder.FishersExactTest(r,n,R,N)
                 ### Store these data in an object
                zsd = mappfinder.ZScoreData(signature,r,n,z,null_z,n)
                zsd.SetP(pval)
                zsd.setAssociatedIDs(r_ls)
                #print splicing_factor,'\t', signature,'\t', z, pval;sys.exit()
                if splicing_factor in splicing_factor_correlated_scores:
                    signature_db = splicing_factor_correlated_scores[splicing_factor]
                    signature_db[signature]=zsd ### Necessary format for the permutation function
                else:
                    signature_db={signature:zsd}
                    splicing_factor_correlated_scores[splicing_factor] = signature_db
            
    results_dir = output_dir+'/SFEnrichmentResults'
    result_file = results_dir+'/SF-correlated_SignatureScores.txt'
    try: os.mkdir(results_dir)
    except: pass
    eo=open(result_file,'w')
    eo.write(string.join(['Splicing Factor','Signature', 'Number Changed', 'Number Measured', 'Z-score','FisherExactP','AdjustedP'],'\t')+'\n') #'Events'
    
    ### Perform a permutation analysis to get BH adjusted p-values
    for splicing_factor in splicing_factor_correlated_scores:
        sorted_results=[]
        signature_db = splicing_factor_correlated_scores[splicing_factor]
        ### Updates the adjusted p-value instances
        mappfinder.adjustPermuteStats(signature_db)
        for signature in signature_db:
            zsd = signature_db[signature]
            if float(zsd.ZScore())>1.96 and float(zsd.Changed())>2 and float(zsd.PermuteP())<0.05:
                enriched_SFs={}
                results = [splicing_factor,signature, zsd.Changed(), zsd.Measured(), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP()] #string.join(zsd.AssociatedIDs(),'|')
                sorted_results.append([float(zsd.PermuteP()),results])
        sorted_results.sort() ### Sort by p-value
        for (p,values) in sorted_results:
            eo.write(string.join(values,'\t')+'\n')
        if len(sorted_results)==0:
            eo.write(string.join([splicing_factor,'NONE','NONE','NONE','NONE','NONE','NONE'],'\t')+'\n')
    eo.close()
    
def Zscore(r,n,N,R):
    """where N is the total number of events measured: 
    R is the total number of events meeting the criterion:
    n is the total number of events in this specific reference gene-set: 
    r is the number of events meeting the criterion in the examined reference gene-set: """
    N=float(N) ### This bring all other values into float space
    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
    return z

if __name__ == '__main__':

    try:
        import multiprocessing as mlp
        mlp.freeze_support()
    except Exception:
        mpl = None

    ################  Default Variables ################
    species = 'Hs'
    platform = "RNASeq"
    useMulti = False
    output_dir = None
    eventDir = None
    PSI_dir = None
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print 'Example: SpliceEnricher.py --PSI "/Data/PSI_data.txt" --geneExp "/Data/GeneExp_data.txt" --geneList "/Data/SplicingFactors.txt" --rho 0.5'
    else:
        try:
            options, remainder = getopt.getopt(sys.argv[1:],'', ['PSI=','species=','o=','platform=','useMulti=',
                                                    'geneExp=','geneList=','rho=','eventDir='])
        except Exception,e:
            print "Error",e
        for opt, arg in options:
            if opt == '--PSI': PSI_dir=arg
            elif opt == '--geneExp': exp_dir=arg
            elif opt == '--geneList': SF_dir=arg
            elif opt == '--species': species=arg
            elif opt == '--o': output_dir=arg
            elif opt == '--platform': platform=arg
            elif opt == '--rho': PearsonCutoff=float(arg)
            elif opt == '--eventDir': eventDir=arg
    
    if output_dir==None:
        output_dir = string.replace(PSI_dir,'\\','/')
        output_dir = string.join(string.split(output_dir,'/')[:-1],'/')
    if PSI_dir !=None:
        findcorrelations(SF_dir, PSI_dir, exp_dir, output_dir, PearsonCutoff)
    if eventDir !=None:
        performEventEnrichment(output_dir,eventDir,species)
