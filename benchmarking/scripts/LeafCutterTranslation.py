#!/usr/bin/env python
import os,sys,string,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)
import unique

#Compare splicing annotations in LeafCutter and AltAnalyze's PSI EventAnnotation file

def verifyFile(filename):
    status = False
    try:
        fn=unique.filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = True;break
    except Exception: status = False
    return status

def importDatabaseEventAnnotations(species,platform):
    terminal_exons={}
    header=True
    count=0
    fn = 'AltDatabase/'+species+'/'+platform+'/'+species+'_Ensembl_exons.txt'
    fn = unique.filepath(fn)
    for line in open(fn,'rU'):
        line = line.rstrip('\n')
        values = string.split(line,'\t')
        if header:
            eI = values.index('splice_events')
            header=False
            continue
        
        exon = values[0]
        event = values[eI]
        if 'alt-N-term' in event or 'altPromoter' in event:
            if 'cassette' not in event:
                terminal_exons[exon] = 'altPromoter'
                count+=1    
        elif 'alt-C-term' in event:
            if 'cassette' not in event:
                terminal_exons[exon] = 'alt-C-term'
                count+=1
        """
        elif 'bleedingExon' in event or 'altFinish' in event:
            terminal_exons[exon] = 'bleedingExon'
            count+=1"""
    print count, 'terminal exon annotations stored'
    return terminal_exons

def formatFeatures(features):
    features2=[]
    for f in features:
        f = string.split(f,'|')
        direction = f[-1]
        annotation = f[0]
        f = '('+direction+')'+annotation
        features2.append(f)
    return features2

class SplicingAnnotations(object):
    def __init__(self, symbol, description,junc1,junc2,altExons,proteinPredictions,eventAnnotation,coordinates):
        self.symbol = symbol
        self.description = description
        self.junc1 = junc1
        self.junc2 = junc2
        self.altExons = altExons
        self.proteinPredictions = proteinPredictions
        self.eventAnnotation = eventAnnotation
        self.coordinates = coordinates
    def Symbol(self): return self.symbol
    def Description(self): return self.description
    def Junc1(self): return self.junc1
    def Junc2(self): return self.junc2
    def AltExons(self): return self.altExons
    def ProteinPredictions(self): return self.proteinPredictions
    def EventAnnotation(self): return self.eventAnnotation
    def Coordinates(self): return self.coordinates

def importPSIAnnotations(PSIpath):
    header=True
    count=0
    events={}
    for line in open(PSIpath,'rU').xreadlines():
        line = line.rstrip('\n')
        values = string.split(line,'\t')
        if header:
            sI = values.index('Symbol')
            dI = values.index('Description')
            eI = values.index('Examined-Junction')
            bI = values.index('Background-Major-Junction')
            aI = values.index('AltExons')
            pI = values.index('ProteinPredictions')
            vI = values.index('EventAnnotation')
            cI = values.index('Coordinates')
            rI = values.index('rawp')
            header=False
        else:
            symbol = values[sI]
            description = values[dI]
            junc1 = values[eI]
            junc2 = values[bI]
            altExons = values[aI]
            proteinPredictions = values[pI]
            eventAnnotation = values[vI]
            coordinates = values[cI]
            key = symbol+':'+junc1+"|"+junc2
            sa = SplicingAnnotations(symbol, description,junc1,junc2,altExons,proteinPredictions,eventAnnotation,coordinates)
            coord1,coord2 = string.split(coordinates,'|')
            events[coord1] = sa
            events[coord2] = sa
    return events

def importLeafCutterJunctions(leafcutter_clusters):
    count=0
    coordinate_clusters={}
    for line in open(leafcutter_clusters,'rU').xreadlines():
        line = line.rstrip('\n')
        if ':' in line:
            cluster = string.split(line,' ')[0]
            cluster = string.split(cluster,':')
            cluster_id = cluster[-1]
            coordinates = cluster[0]+':'+cluster[1]+'-'+cluster[1]
            try: coordinate_clusters[cluster_id].append(coordinates)
            except Exception: coordinate_clusters[cluster_id] = [coordinates]
    return coordinate_clusters

def importLeafSignificant(leafcutter_significant,coordinate_clusters):
    header=True
    count=0
    significant_coordinates={}
    unique_clusters=0
    for line in open(leafcutter_significant,'rU').xreadlines():
        values = line.rstrip('\n')
        values = string.split(values,'\t')
        if header:
            pI = values.index('p')
            header=False
        else:
            cluster = string.split(values[0],':')[1]
            try:
                p = float(values[pI])
            except Exception:
                p=1
            coordinates = coordinate_clusters[cluster]
            if p<0.05:
                unique_clusters+=1
                for coord in coordinates:
                    significant_coordinates[coord] = cluster
    print unique_clusters,'significant clusters'
    return significant_coordinates

def compare_algorithms(altanalyze_dir,leafcutter_significant,leafcutter_clusters,species,platform):
    """ Add splicing annotations for PSI results """
    
    ### Get all splice-junction pairs
    coordinate_clusters = importLeafCutterJunctions(leafcutter_clusters)
    significant_coordinates = importLeafSignificant(leafcutter_significant,coordinate_clusters)
    eventAnnotations = importPSIAnnotations(altanalyze_dir)
    unique_clusters={}
    for i in significant_coordinates:
        if i not in eventAnnotations:
            print i;sys.exit()
        else:
            print 'a',i;sys.exit()

    ### Get all de novo junction anntations (includes novel junctions) 
    psievents = importEventAnnotations(resultsDir,species,psievents,annotationType='de novo')
    ### Get all known junction annotations
    psievents = importEventAnnotations(resultsDir,species,psievents)
    ### Import the annotations that provide alternative terminal events
    terminal_exons = importDatabaseEventAnnotations(species,platform)
    
    ### Update our PSI annotation file with these improved predictions
    updatePSIAnnotations(PSIpath, species, psievents, terminal_exons, junctionPairFeatures)
    
if __name__ == '__main__':
    import multiprocessing as mlp
    import getopt
    #bam_dir = "H9.102.2.6.bam"
    #outputExonCoordinateRefBEDfile = 'H9.102.2.6__exon.bed'
    platform = 'RNASeq'
    species = 'Hs'
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []
        useMultiProcessing=False
        options, remainder = getopt.getopt(sys.argv[1:],'', ['altanalyze=','leafcutter=','leafclusters=','species=','platform=','array='])
        for opt, arg in options:
            if opt == '--altanalyze': altanalyze_dir=arg
            elif opt == '--leafcutter': leafcutter_significant=arg
            elif opt == '--leafclusters': leafcutter_clusters=arg
            elif opt == '--species': species=arg
            elif opt == '--platform': platform=arg
            elif opt == '--array': platform=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    compare_algorithms(altanalyze_dir,leafcutter_significant,leafcutter_clusters,species,platform)
