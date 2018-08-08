#!/usr/bin/env python
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import unique

#Annotate PSI file's splicing events using alternate_juntion  and alternate_junction_de-novo

def verifyFile(filename):
    status = False
    try:
        fn=unique.filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = True;break
    except Exception: status = False
    return status
 
class SplicingEventAnnotation:
    def __init__(self, event, critical_exon):
        self.event = event; self.critical_exon = critical_exon
    def Event(self): return self.event
    def CriticalExon(self): return self.critical_exon
    def setJunctions(self,junctions): self.junctions = junctions
    def Junctions(self): return self.junctions
    def ExonOnly(self):
        try: exon = string.split(self.critical_exon,':')[1]
        except Exception: exon = self.critical_exon
        return exon
    
def importPSIevents(PSIpath,species):
    header=True
    count=0;count1=0;
    initial_events={}
    psievents={}
    psijunctions={}
    # Create a dictionary for the psi event dict[min_isoform,major_isoform] and value will be the annotation that will assigned later in the code
    for line in open(PSIpath,'rU').xreadlines():
        line = line.rstrip('\n')
        #line=string.replace(line,"_",".")
        values=string.split(line,'\t')
        if header:
            aI = values.index('AltExons')
            header=False
            continue
        primary_junction = values[2]
        secondary_junction = values[3]
        critical_exon = values[aI]
        se = SplicingEventAnnotation(None,critical_exon)
        se.setJunctions([primary_junction,secondary_junction])
        psievents[primary_junction,secondary_junction] = se
        psijunctions[(primary_junction,)] = se ### format for importing protein annotations
        psijunctions[(secondary_junction,)] = se
    print"PSI events imported..."
    return psievents,psijunctions

def importEventAnnotations(resultsDir,species,psievents,annotationType=None):
    # Parse to the de-novo junction file and update the value of the psi events in the dictionary 
    header=True
    if annotationType == 'de novo':
        junction_file = string.split(resultsDir,'AltResults')[0]+'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions_de-novo.txt'
    else:
        junction_file = 'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions.txt'
    count=0
    fn = unique.filepath(junction_file)
    initial_events={}
    initial_critical_exons={}
    status = verifyFile(fn)
    if status:
        for line in open(fn,'rU').xreadlines():
            line = line.rstrip('\n')
            if header:
               header=False
               continue
            gene, critical_exon, j1, j2, event = string.split(line,'\t')
            critical_exon = gene+':'+critical_exon
            ### Fix this notation
            if len(j1)>17 and '_' not in j1:
                j1a,j1b = string.split(j1,'-')
                if len(j1a)>8:
                    a,b,c = string.split(j1a,'.')
                    j1a = a+'.'+b+'_'+c
                if len(j1b)>8:
                    a,b,c = string.split(j1b,'.')
                    j1b = a+'.'+b+'_'+c
                j1 = j1a+'-'+j1b
            if len(j2)>17 and '_' not in j2:
                j2a,j2b = string.split(j2,'-')
                if len(j2a)>8:
                    a,b,c = string.split(j2a,'.')
                    j2a = a+'.'+b+'_'+c
                if len(j2b)>8:
                    a,b,c = string.split(j2b,'.')
                    j2b = a+'.'+b+'_'+c
                j2 = j2a+'-'+j2b
            j1=gene+':'+j1
            j2=gene+':'+j2
            if ((j1,j2) in psievents) or ((j2,j1) in psievents):
                try: initial_events[j1,j2].append(event)
                except Exception: initial_events[j1,j2] = [event]
                try: initial_critical_exons[j1,j2].append(critical_exon)
                except Exception: initial_critical_exons[j1,j2] = [critical_exon]
              
        for junctions in initial_events:
            events = unique.unique(initial_events[junctions])
            events.sort()
            events = string.join(events,'|')
            critical_exons = unique.unique(initial_critical_exons[junctions])
            critical_exons.sort()
            critical_exons = string.join(critical_exons,'|')
            se = SplicingEventAnnotation(events,critical_exons)
            j1,j2 = junctions
            psievents[j1,j2] = se
            psievents[j2,j1] = se
            count+=1
        print count, "junction annotations imported..."

    return psievents

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

def inverseFeatureDirections(features):
    features2=[]
    for f in features:
        if '+' in f:
            f = string.replace(f,'+','-')
        else:
            f = string.replace(f,'-','+')
        features2.append(f)
    return features2

def formatFeatures(features):
    features2=[]
    for f in features:
        f = string.split(f,'|')
        direction = f[-1]
        annotation = f[0]
        f = '('+direction+')'+annotation
        features2.append(f)
    return features2

def importIsoformAnnotations(species,platform,psievents,annotType=None,junctionPairFeatures={},dataType='reciprocal'):
    count=0
    if annotType == 'domain':
        if dataType == 'reciprocal':
            fn = 'AltDatabase/'+species+'/'+platform+'/'+'probeset-domain-annotations-exoncomp.txt'
        else:
            fn = 'AltDatabase/'+species+'/'+platform+'/junction/'+'probeset-domain-annotations-exoncomp.txt'
    else:
        if dataType == 'reciprocal':
            fn = 'AltDatabase/'+species+'/'+platform+'/'+'probeset-protein-annotations-exoncomp.txt'
        else:
            fn = 'AltDatabase/'+species+'/'+platform+'/junction/'+'probeset-protein-annotations-exoncomp.txt'
    fn = unique.filepath(fn)
    for line in open(fn,'rU'):
        line = line.rstrip('\n')
        values = string.split(line,'\t')
        junctions = string.split(values[0],'|')
        features = formatFeatures(values[1:])
        antiFeatures = inverseFeatureDirections(features)
        if tuple(junctions) in psievents:
            try: junctionPairFeatures[tuple(junctions)].append(string.join(features,', '))
            except Exception: junctionPairFeatures[tuple(junctions)] = [string.join(features,', ')]
        if dataType == 'reciprocal':
            junctions.reverse()
            if tuple(junctions) in psievents:
                try: junctionPairFeatures[tuple(junctions)].append(string.join(antiFeatures,', '))
                except Exception: junctionPairFeatures[tuple(junctions)] = [string.join(antiFeatures,', ')]
        count+=1
    print count, 'protein predictions added'
    return junctionPairFeatures

def DetermineIntronRetention(coordinates):
    intronRetention = False
    coordinates1,coordinates2 = string.split(coordinates,'|')
    coordinates1 = string.split(coordinates1,':')[1]
    coordinate1a, coordinate1b = string.split(coordinates1,'-')
    coordinate1_diff = abs(float(coordinate1a)-float(coordinate1b))
    coordinates2 = string.split(coordinates2,':')[1]
    coordinate2a, coordinate2b = string.split(coordinates2,'-')
    coordinate2_diff = abs(float(coordinate2a)-float(coordinate2b))
    if coordinate1_diff==1 or coordinate2_diff==1:
        intronRetention = True
    return intronRetention    
    
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
    annotations={}
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
            annotations[key] = sa
    return annotations
            
def updatePSIAnnotations(PSIpath, species, psievents, terminal_exons, junctionPairFeatures, junctionFeatures):
    # write the updated psi file with the annotations into a new file and annotate events that have not been annotated by the junction files
    #print len(psievents)
    header=True
    export_path = PSIpath[:-4]+'_EventAnnotation.txt'
    export=open(export_path,'w')
    count=0
    for line in open(PSIpath,'rU').xreadlines():
        line = line.rstrip('\n')
        values = string.split(line,'\t')
        if header:
            try: fI = values.index('feature')
            except: fI = values.index('EventAnnotation')
            aI = values.index('AltExons')
            try: pI = values.index('PME')
            except: pI = values.index('ProteinPredictions')
            cI = values.index('Coordinates')
            values[fI] = 'EventAnnotation'
            values[pI] = 'ProteinPredictions'
            export.write(string.join(values,'\t')+'\n')
            header=False
            continue
        psiJunction=string.split(line,'\t')
        psiJunction_primary = psiJunction[2]
        psiJunction_secondary = psiJunction[3]
        key = psiJunction_primary,psiJunction_secondary
        #psiJunction_primary=string.replace(psiJunction[2],"_",".")
        #psiJunction_secondary=string.replace(psiJunction[3],"_",".")
        event = psievents[psiJunction_primary,psiJunction_secondary].Event()
        critical_exon = psievents[key].CriticalExon()
        if key in junctionPairFeatures:
            proteinAnnotation = string.join(junctionPairFeatures[key],'|')
        elif (psiJunction_primary,) in junctionFeatures:
            proteinAnnotation = string.join(junctionPairFeatures[(psiJunction_primary,)],'|')
        #elif (psiJunction_secondary,) in junctionFeatures:
            #proteinAnnotation = string.join(junctionPairFeatures[(psiJunction_secondary,)],'|')"""
        else:
            proteinAnnotation=''
        values[pI] = proteinAnnotation
        intronRetention = DetermineIntronRetention(values[cI])
        
        if critical_exon in terminal_exons:
            event = terminal_exons[critical_exon]
        if intronRetention:
            event = 'intron-retention'
        try:
            if event==None:
                primary_exons=string.split(psiJunction[2],":")
                secondary_exons=string.split(psiJunction[3],":")
                if len(primary_exons)>2 or len(secondary_exons)>2:
                    values[fI] = 'trans-splicing'
                    export.write(string.join(values,'\t')+'\n')
                    continue
                else:
                    primary_exonspos=string.split(primary_exons[1],"-")
                    secondary_exonspos=string.split(secondary_exons[1],"-")
                    if ('U0' in primary_exons[1]) or ('U0' in secondary_exons[1]):
                        if ('U0.' in primary_exonspos[0]) or ('U0.' in secondary_exonspos[0]):
                            values[fI] = 'altPromoter'
                            export.write(string.join(values,'\t')+'\n')
                            continue
                        else:
                            values[fI] = ''
                            export.write(string.join(values,'\t')+'\n')
                            continue
                    try: event = predictSplicingEventTypes(psiJunction_primary,psiJunction_secondary)
                    except Exception:
                        event = ''
                    values[fI] = event
                    export.write(string.join(values,'\t')+'\n')
                    continue

            else:
                values[fI] = event
                values[aI] = critical_exon
                count+=1
                export.write(string.join(values,'\t')+'\n')
        except Exception():
            values[fI] = ''
            export.write(string.join(values,'\t')+'\n')
    #print count
    return export_path

def predictSplicingEventTypes(junction1,junction2):
    if 'I' not in junction1 and '_' in junction1:
        junction1 = string.replace(junction1,'_','') ### allows this to be seen as an alternative splice site
    if 'I' not in junction2 and '_' in junction2:
        junction2 = string.replace(junction2,'_','') ### allows this to be seen as an alternative splice site
    if 'I' in junction1:
        forceError
    if 'I' in junction2:
        forceError

    j1a,j1b = string.split(junction1,'-')
    j2a,j2b = string.split(junction2,'-')
    j1a = string.split(j1a,':')[1][1:]
    j2a = string.split(j2a,':')[1][1:]  
    
    j1a,r1a = string.split(j1a,'.')
    j1b,r1b = string.split(j1b[1:],'.')

    j2a,r2a = string.split(j2a,'.')
    j2b,r2b = string.split(j2b[1:],'.')
    
    ### convert to integers
    j1a,r1a,j1b,r1b,j2a,r2a,j2b,r2b = map(lambda x: int(float(x)),[j1a,r1a,j1b,r1b,j2a,r2a,j2b,r2b])

    splice_event=[]
    if j1a == j2a and j1b==j2b: ### alt-splice site
        if r1a == r2a: splice_event.append("alt-3'")
        else: splice_event.append("alt-5'")
    elif j1a == j2a: splice_event.append("cassette-exon")
    elif j1b==j2b:
        if 'E1.' in junction1: splice_event.append("altPromoter")
        else: splice_event.append("cassette-exon")
    elif 'E1.' in junction1 or 'E1.1' in junction2:
        splice_event.append("altPromoter")
    else:
        splice_event.append("cassette-exon")
    splice_event = unique.unique(splice_event)
    splice_event.sort()
    splice_event = string.join(splice_event,'|')

    return splice_event
 
def parse_junctionfiles(resultsDir,species,platform):
    """ Add splicing annotations for PSI results """
    
    if 'top_alt' not in resultsDir:
        PSIpath = resultsDir+'/'+species+'_'+platform+'_top_alt_junctions-PSI.txt'
    else:
        PSIpath = resultsDir
    
    ### Get all splice-junction pairs
    psievents,psijunctions = importPSIevents(PSIpath,species)
    
    ### Get domain/protein predictions
    junctionPairFeatures = importIsoformAnnotations(species,platform,psievents)
    junctionPairFeatures = importIsoformAnnotations(species,platform,psievents,annotType='domain',junctionPairFeatures=junctionPairFeatures)

    junctionFeatures = importIsoformAnnotations(species,platform,psijunctions,dataType='junction')
    junctionFeatures = importIsoformAnnotations(species,platform,psijunctions,annotType='domain',junctionPairFeatures=junctionFeatures,dataType='junction')
    
    ### Get all de novo junction anntations (includes novel junctions) 
    psievents = importEventAnnotations(resultsDir,species,psievents,annotationType='de novo')
    ### Get all known junction annotations
    psievents = importEventAnnotations(resultsDir,species,psievents)
    ### Import the annotations that provide alternative terminal events
    terminal_exons = importDatabaseEventAnnotations(species,platform)
    
    ### Update our PSI annotation file with these improved predictions
    export_path = updatePSIAnnotations(PSIpath, species, psievents, terminal_exons, junctionPairFeatures, junctionFeatures)
    return export_path
    
if __name__ == '__main__':
    import multiprocessing as mlp
    import getopt
    #bam_dir = "H9.102.2.6.bam"
    #outputExonCoordinateRefBEDfile = 'H9.102.2.6__exon.bed'
    platform = 'RNASeq'
    species = 'Hs'
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a directory containing BAM files as input in the command-line"
        sys.exit()
    else:
        analysisType = []
        useMultiProcessing=False
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','platform=','array='])
        for opt, arg in options:
            if opt == '--i': resultsDir=arg
            elif opt == '--species': species=arg
            elif opt == '--platform': platform=arg
            elif opt == '--array': platform=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    parse_junctionfiles(resultsDir,species,platform)
