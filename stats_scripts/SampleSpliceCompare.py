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
                ### If no annotations are present
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
    eo=open(results_dir+'/SF_correlations.txt','w')
    eo.write('Splicing Factor'+'\t'+'Events Count'+'\n')
    
    counter=0
    gene_correlation_time = []
    #print len(genesToExamine)
    count=0
    for i in genesToExamine:
        print [i]
        if count ==10: break
        count+=1
    count=0
    #print 'next'

    for gene in genesToExamine:
        gene_name = genesToExamine[gene]
        #print len(geneExpression_db),[gene]
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
        
        gene_name = string.replace(gene_name,':',"__")
        if '|' in gene_name:
            gene_name = string.split(gene_name,'|')[0]
        eo.write(gene_name+'\t'+str(count)+'\n')
        filename=results_dir+"/"+gene_name+"_"+str(count)+".txt"
        if count>20:
            eg=open(filename,"w")
            eg.write("SplicingEvent\tSystemCode\tPearsonRho\n")
        
            for (junctionID,rho) in Corrsflist:
                eg.write(junctionID+"\t"+"Ae\t"+str(rho)+"\n")
            eg.close()
        counter+=1
        print '.',
    print '\n...Correlations obtained on average of %d seconds/gene' % numpy.mean(gene_correlation_time)


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
        options, remainder = getopt.getopt(sys.argv[1:],'', ['reference=','output=','input='])

        for opt, arg in options:
            if opt == '--reference': PSI_dir=arg
            elif opt == '--input': dSPI_dir=arg
            elif opt == '--output': output_dir=arg
    
    findcorrelations(PSI_dir, dPSI_dir, output_dir)

