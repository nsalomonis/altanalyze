import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import export
import unique
import traceback

""" Intersecting Coordinate Files """

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def eCLIPimport(folder):
    eCLIP_dataset_peaks={}
    annotations=[]
    files = unique.read_directory(folder)
    for file in files:
        if '.bed' in file:
            peaks={}
            dataset = file[:-4]
            print dataset
            key_db={}
            fn = unique.filepath(folder+'/'+file)
            eo = export.ExportFile(folder+'/3-prime-peaks/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                chr = t[0]
                start = int(t[1])
                end = int(t[2])
                strand = t[5]
                gene = string.split(t[-3],'.')[0]
                annotation = string.split(t[-3],';')[-1]
                if 'three_prime_utrs' in annotation:
                    eo.write(line)
                if annotation not in annotations:
                    annotations.append(annotation)
                symbol = t[-2]
                key = chr,start,strand
                #"""
                if gene in coding_db:
                    coding_type = coding_db[gene][-1]
                    if 'protein_coding' in coding_type:
                        coding_type = 'protein_coding'
                ##"""
                        if key in key_db:
                            gene_db = key_db[key]
                            if gene in gene_db:
                                gene_db[gene].append(annotation)
                            else:
                                gene_db[gene]=[annotation]
                        else:
                            gene_db={}
                            gene_db[gene]=[annotation]
                            key_db[key]=gene_db
            for key in key_db:
                ranking=[]
                for gene in key_db[key]:
                    ranking.append((len(key_db[key][gene]),gene))
                ranking.sort()
                gene = ranking[-1][-1]
                for annotation in key_db[key][gene]:
                    if annotation in peaks:
                        peaks[annotation]+=1
                    else:
                        peaks[annotation]=1
            eCLIP_dataset_peaks[dataset]=peaks
    eo.close()
    
    annotations.sort()
    eo = export.ExportFile(folder+'/summary-annotations/summary.txt')
    header = string.join(['RBP']+map(str,annotations),'\t')+'\n'
    eo.write(header)
    for dataset in eCLIP_dataset_peaks:
        annot=[]
        peaks = eCLIP_dataset_peaks[dataset]
        for annotation in annotations:
            if annotation in peaks:
                annot.append(peaks[annotation])
            else:
                annot.append(0)
        annot = map(lambda x: (1.000*x/sum(annot)), annot)
        values = string.join([dataset]+map(str,annot),'\t')+'\n'
        eo.write(values)
    eo.close()

if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt
    CLIP_dir = None
    species = 'Hs'
    
    """ Usage:
    bedtools intersect -wb -a /Clip_merged_reproducible_ENCODE/K562/AARS-human.bed -b /annotations/combined/hg19_annotations-full.bed > /test.bed
    """

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['species=','clip='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--species':
                species = arg
            elif opt == '--clip':
                CLIP_dir = arg

    import ExpressionBuilder
    coding_db = ExpressionBuilder.importTranscriptBiotypeAnnotations(species)
    dataset_peaks = eCLIPimport(CLIP_dir)

            