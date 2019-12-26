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

class EventInformation:
    def __init__(self, id, event_direction, clusterID, altExons, coordinates):
        self.id = id
        self.event_direction = event_direction
        self.clusterID = clusterID
        self.altExons = altExons
        self.coordinates = coordinates
    def ID(self):
        return self.id
    def GeneID(self):
        return string.split(self.id,':')[1]
    def Symbol(self):
        return string.split(self.id,':')[0]
    def EventDirection(self):
        return self.event_direction
    def ClusterID(self):
        return self.clusterID
    def AltExons(self):
        """ Can be multiple exons """
        altExons = string.split(self.altExons,'|')
        return altExons
    
    def JunctionCoordinateInterval(self):
        """ If the exon block is not in the database, return the full interval """
        try:
            j1,j2 = string.split(self.Coordinates(),'|')
            chr,coord1 = string.split(j1,':')
            chr,coord2 = string.split(j2,':')
            coords = map(int,string.split(coord1,'-'))
            coords += map(int,string.split(coord2,'-'))
            coords.sort()
            return chr,coords[0],coords[-1]
        except:
            """ See the class circInformation """
            return self.coordinates
    
    def AltExonBlockCoord(self):
        exon_block_coords=[]
        for altExon in self.AltExons():
            if 'ENSG' not in altExon:
                altExon = self.GeneID() + ':' + altExon
            altExon_block = string.split(altExon,'.')[0]
            try:
                chr, strand, start, end = exon_block_coordinates[altExon_block]
            except:
                chr, start, end = self.JunctionCoordinateInterval()
            exon_block_coords.append(chr+':'+str(start)+'-'+str(end))
        return string.join(exon_block_coords,'|')
    
    def FlankingBlockCoord(self):
        coords=[]
        for altExon in self.AltExons():
            if 'ENSG' in altExon:
                altExon = string.split(altExon,':')[1]
            altExon_block = string.split(altExon,'.')[0]
            block_num = int(altExon_block[1:])
            block_type = altExon_block[0]
            try:
                chr, strand, start, end = exon_block_coordinates[self.GeneID() + ':' + block_type+str(block_num)]
                coords+=[start, end]
            except:
                pass
            if block_type == 'E':
                if extendFlanking == False:
                    upstream_intron = self.GeneID() + ':' + "I"+str(block_num-1)
                else:
                    upstream_intron = self.GeneID() + ':' + "E"+str(block_num-1)
                    if upstream_intron not in exon_block_coordinates:
                        upstream_intron = self.GeneID() + ':' + "I"+str(block_num-1)
                try:
                    chr, strand, Ustart, Uend = exon_block_coordinates[upstream_intron]
                except:
                    chr, Ustart, Uend = self.JunctionCoordinateInterval()
                if extendFlanking == False:
                    downstream_intron = self.GeneID() + ':' + "I"+str(block_num)
                else:
                    downstream_intron = self.GeneID() + ':' + "E"+str(block_num+1)
                    if downstream_intron not in exon_block_coordinates:
                        downstream_intron = self.GeneID() + ':' + "I"+str(block_num)
                try:
                    chr, strand, Dstart, Dend = exon_block_coordinates[downstream_intron]
                except:
                    chr, Dstart, Dend = self.JunctionCoordinateInterval()
            else:
                upstream_exon = self.GeneID() + ':' + "E"+str(block_num)
                try:
                    chr, strand, Ustart, Uend = exon_block_coordinates[upstream_exon]
                except:
                    chr, Ustart, Uend = self.JunctionCoordinateInterval()
                downstream_exon = self.GeneID() + ':' + "E"+str(block_num+1)
                try:
                    chr, strand, Dstart, Dend = exon_block_coordinates[downstream_exon]
                except:
                    chr, Dstart, Dend = self.JunctionCoordinateInterval()
            coords += [Ustart, Uend, Dstart, Dend]
        coords.sort()
        start = coords[0]
        end = coords[-1]
        return start, end

    def Coordinates(self):
        return self.coordinates
    def Export(self):
        annotations = [self.Symbol(), self.ID(), self.EventDirection(), self.ClusterID(), self.Coordinates(), self.altExons, self.AltExonBlockCoord()]
        return annotations
    def __repr__(self):
        return self.ID(), self.EventDirection(), self.ClusterID(), self.Coordinates()
    
    
class circInformation(EventInformation):

    def __init__(self, id, pval, logFC):
        self.id = id
        self.pval = pval
        self.GeneID()
        self.Coordinates()
        self.FindAssociatedExonsBlocks()
        self.logFC = logFC
        self.AltExons()
        
    def pVal(self):
        return self.pval

    def EventDirection(self):
        if self.logFC > 0:
            return 'inclusion'
        else:
            return 'exclusion'
    
    def Symbol(self):
        symbol = string.split(self.ID(),'-')[0]
        self.symbol = symbol
        return self.symbol
    
    def GeneID(self):
        if self.Symbol() in symbol_to_gene:
            geneID=symbol_to_gene[self.Symbol()][0]
            for gene in symbol_to_gene[self.Symbol()]:
                if 'ENS' in gene:
                    geneID = gene
        else:
            geneID = self.Symbol()
        return geneID

    def Coordinates(self):
        c1,c2 = string.split(self.ID(),'-')[-2:]
        self.chr,c1 = string.split(c1,':')
        self.start = int(c1)
        self.end = int(c2)
        self.coordinates = self.chr, self.start, self.end
        return self.chr+':'+str(self.start)+'-'+str(self.end)
    
    def FindAssociatedExonsBlocks(self):
        """ Indentify matching exon/intron regions for the circRNA coordinates """
        
        chr, cstart, cend = self.coordinates
        circRNA_coords = [cstart, cend]
        circRNA_coords.sort()
        if self.GeneID() in gene_to_exons:
            search_blocks = gene_to_exons[self.GeneID()]
            aligned_blocks=[]
            for exon_block in search_blocks:
                chr, strand, start, end = exon_block_coordinates[exon_block]
                block_coords = [start, end]
                block_coords.sort()
                coords = circRNA_coords+block_coords
                coords.sort()
                if len(unique.unique(coords))==3:
                    if exon_block not in aligned_blocks:
                        aligned_blocks.append(exon_block)
                elif coords[:2]==circRNA_coords or coords[-2:]==circRNA_coords:
                    pass
                else:
                    if exon_block not in aligned_blocks:
                        aligned_blocks.append(exon_block)
            self.aligned_blocks = aligned_blocks
        else:
            self.aligned_blocks=[]
        return self.aligned_blocks
        
    def ClusterID(self):
        return self.GeneID()
    
    def AltExons(self):
        self.altExons = string.join(self.aligned_blocks,'|')
        return self.aligned_blocks

    def Chr(self):
        return self.chr
    def Start(self):
        return self.start
    def End(self):
        return self.end
    def Strand(self):
        return ''
    def Annotation(self):
        return ''

class PeakInformation:
    def __init__(self, chr, start, end, strand, annotation, gene, symbol):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        self.annotation = annotation
        self.symbol = symbol
        self.gene = gene
    def Chr(self):
        return self.chr
    def Start(self):
        return self.start
    def End(self):
        return self.end
    def Strand(self):
        return self.strand
    def GeneID(self):
        return self.gene
    def Annotation(self):
        return self.annotation
    def Symbol(self):
        return self.symbol
    def Coordinates(self):
        return self.chr+':'+str(self.Start())+'-'+str(self.End())
    def Export(self):
        annotations = [self.Coordinates(),self.Annotation()]
        return annotations
    def __repr__(self):
        return self.Symbol(), self.Annotation()

def importSplicingEvents(folder):
    dataset_events={}
    files = unique.read_directory(folder)
    for file in files:
        if 'PSI.' in file and '.txt' in file:
            events=[]
            dataset = file[:-4]
            fn = unique.filepath(folder+'/'+file)
            firstRow=True
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if firstRow:
                    index=0
                    """ Standard Fields from MultiPath-PSI """
                    for i in t:
                        if 'Event-Direction'==i:
                            ed = index
                        if 'ClusterID' == i:
                            ci = index
                        if 'AltExons' == i:
                            ae = index
                        if 'EventAnnotation' == i:
                            ea = index
                        if 'Coordinates' == i:
                            co = index
                        index+=1
                    firstRow = False
                else:
                    id = t[0]
                    event_direction = t[ed]
                    clusterID = t[ci]
                    altExons = t[ae]
                    coordinates = t[co]
                    ei = EventInformation(id, event_direction, clusterID, altExons, coordinates)
                    events.append(ei)
            dataset_events[dataset]=events
    return dataset_events

def eCLIPimport(folder):
    eCLIP_dataset_peaks={}
    files = unique.read_directory(folder)
    for file in files:
        if '.bed' in file:
            peaks=[]
            dataset = file[:-4]
            fn = unique.filepath(folder+'/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                chr = t[0]
                start = int(t[1])
                end = int(t[2])
                strand = t[5]
                annotation = t[6]
                gene = string.split(t[8],'.')[0]
                symbol = t[-2]
                pi = PeakInformation(chr, start, end, strand, annotation, gene, symbol)
                peaks.append(pi)
            eCLIP_dataset_peaks[dataset]=peaks
    return eCLIP_dataset_peaks

def importExonCoordinates(species):
    """ Import exon block, intron block and gene coordinates """
    
    firstRow=True
    exon_coordinate_path = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    fn = unique.filepath(exon_coordinate_path)
    gene_coordinates={}
    gene_to_exons={}
    exon_block_coordinates={}
    gene_chr_strand = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            firstRow = False
        else:
            gene, exonid, chr, strand, exon_region_starts, exon_region_ends, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            exon_region_starts  = map(int,string.split(exon_region_starts,'|'))
            exon_region_ends  = map(int,string.split(exon_region_ends,'|'))
            exon_block = gene+':'+string.split(exonid,'.')[0]
            gene_chr_strand[gene]=chr,strand
            if gene in gene_to_exons:
                gene_to_exons[gene].append(exon_block)
            else:
                gene_to_exons[gene] = [exon_block]
            if gene in gene_coordinates:
                gene_coordinates[gene]+=exon_region_starts+exon_region_ends
            else:
                gene_coordinates[gene]=exon_region_starts+exon_region_ends
            if exon_block in exon_block_coordinates:
                exon_block_coordinates[exon_block]+=exon_region_starts+exon_region_ends
            else:
                exon_block_coordinates[exon_block]=exon_region_starts+exon_region_ends
        
    for gene in gene_coordinates:
        gene_coordinates[gene].sort()
        start = gene_coordinates[gene][0]
        end = gene_coordinates[gene][-1]
        chr,strand = gene_chr_strand[gene]
        gene_coordinates[gene]= chr, strand, start, end
    
    for exon in exon_block_coordinates:
        exon_block_coordinates[exon].sort()
        start = exon_block_coordinates[exon][0]
        end = exon_block_coordinates[exon][-1]
        chr,strand = gene_chr_strand[string.split(exon,':')[0]]
        exon_block_coordinates[exon] = chr, strand, start, end
        
    print len(gene_coordinates), 'genes'
    print len(exon_block_coordinates), 'exons/introns'
    return gene_coordinates, exon_block_coordinates, gene_to_exons

def importGeneSymbols(species):
    import gene_associations
    gene_to_symbol = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
    
    from import_scripts import OBO_import
    symbol_to_gene = OBO_import.swapKeyValues(gene_to_symbol)
            
    return gene_to_symbol, symbol_to_gene

def importCircularRNAEvents(folder,circ_p):
    dataset_events={}
    files = unique.read_directory(folder)
    for file in files:
        if 'circRNA.' in file and '.txt' in file:
            events=[]
            dataset = file[:-4]
            fn = unique.filepath(folder+'/'+file)
            firstRow=True
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if firstRow:
                    index=0
                    """ Standard Fields from MultiPath-PSI """
                    for i in t:
                        if 'PValue'==i:
                            pv = index
                        if 'logFC'==i:
                            lf = index
                        index+=1
                    firstRow = False
                else:
                    id = t[0]
                    pval = float(t[pv])
                    logFC = float(t[lf])
                    ci = circInformation(id, pval, logFC)
                    if pval < circ_p:
                        events.append(ci)
            dataset_events[dataset]=events
    return dataset_events

def alignEventsAndPeaks(eCLIP, AS, eCLIP_peaks,AS_events,AS_event_dir):
    """ Compare genomic coordinates from these two datasets """
    eCLIP_gene_peaks = {}
    eCLIP_symbol_peaks = {}
    AS_gene_events = {}
    gene_to_symbol = {}
    
    """ Create gene indexes for both datasets, allowing alternative matches by symbol """
    for pi in eCLIP_peaks:
        if pi.GeneID() not in eCLIP_gene_peaks:
            eCLIP_gene_peaks[pi.GeneID()] = [pi]
        else:
            eCLIP_gene_peaks[pi.GeneID()].append(pi)
        if pi.Symbol() not in eCLIP_symbol_peaks:
            eCLIP_symbol_peaks[pi.Symbol()] = [pi]
        else:
            eCLIP_symbol_peaks[pi.Symbol()].append(pi)
            
    for ei in AS_events:
        if ei.GeneID() not in AS_gene_events:
            AS_gene_events[ei.GeneID()] = [ei]
        else:
            AS_gene_events[ei.GeneID()].append(ei)
        gene_to_symbol[ei.GeneID()] = ei.Symbol()
    
    """ Match peaks to splicing events based on annotated genes (could do by coordinates) """

    gene_export = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_gene.txt')    
    gene_export_incl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_gene-incl.txt')
    gene_export_excl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_gene-excl.txt')
    exon_export = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_exon.txt')
    exon_export_incl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_exon-incl.txt')
    exon_export_excl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_exon-excl.txt')
    flanking_export = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_flanking.txt')
    flanking_export_incl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_flanking-incl.txt')
    flanking_export_excl = export.ExportFile(AS_event_dir+'/eCLIP-overlaps/'+eCLIP+'_'+AS+'_flanking-excl.txt')
    
    header = ['Symbol', 'UID', 'EventDirection', 'ClusterID', 'Coordinates', 'altExons', 'AltExonBlockCoord','Peak-Coordinates','Peak-Annotations']
    header = string.join(header,'\t')+'\n'
    gene_export.write(header)
    gene_export_incl.write(header)
    gene_export_excl.write(header)
    exon_export.write(header)
    exon_export_incl.write(header)
    exon_export_excl.write(header)
    flanking_export.write(header)
    flanking_export_incl.write(header)
    flanking_export_excl.write(header)
    
    exported=[]
    for geneID in AS_gene_events:
        symbol = gene_to_symbol[geneID]
        pi_set = None
        if geneID in eCLIP_gene_peaks:
            pi_set = eCLIP_gene_peaks[geneID]
        elif symbol in eCLIP_symbol_peaks:
            pi_set = eCLIP_symbol_peaks[symbol]
        if pi_set != None:
            """ Matching peak and AS at the gene-level """
            for ei in AS_gene_events[geneID]:
                event_annotations = ei.Export()
                for altExon in ei.AltExons():
                    if 'ENSG' not in altExon:
                        altExon = ei.GeneID() + ':' + altExon
                    altExon_block = string.split(altExon,'.')[0]
                    try:
                        chr, strand, start, end = exon_block_coordinates[altExon_block]
                    except:
                        """ Can occur if the exon region is a novel 5' exon (not in database)
                        use the junction interval instead (less precise) """
                        chr, start, end = ei.JunctionCoordinateInterval()
                    for pi in pi_set:
                        peak_annotations = pi.Export()
                        overlaps =  string.join(event_annotations+peak_annotations,'\t')+'\n'
                        if overlaps in exported:
                            pass ### This comparison has already been exported
                        else:
                            exported.append(overlaps)
                            gene_export.write(overlaps)
                            if ei.EventDirection()=='inclusion':
                                gene_export_incl.write(overlaps)
                            else:
                                gene_export_excl.write(overlaps)
                            
                            """ Find direct exon overlaps """
                            AS_coords = [start,end]
                            AS_coords.sort()
                            Peak_coords = [pi.Start(),pi.End()]
                            Peak_coords.sort()
                            coords = AS_coords+Peak_coords
                            coords.sort()
                            if coords[:2]==AS_coords or coords[-2:]==AS_coords:
                                pass
                            else:
                                exon_export.write(overlaps)
                                if ei.EventDirection()=='inclusion':
                                    exon_export_incl.write(overlaps)
                                else:
                                    exon_export_excl.write(overlaps)
    
                            """ Find indirect flanking intron overlaps """
                            flank_start,flank_end = ei.FlankingBlockCoord()
                            AS_coords = [flank_start,flank_end]
                            AS_coords.sort()
                            Peak_coords = [pi.Start(),pi.End()]
                            Peak_coords.sort()
                            coords = AS_coords+Peak_coords
                            coords.sort()
                            if coords[:2]==AS_coords or coords[-2:]==AS_coords:
                                pass
                            else:
                                flanking_export.write(overlaps)
                                if ei.EventDirection()=='inclusion':
                                    flanking_export_incl.write(overlaps)
                                else:
                                    flanking_export_excl.write(overlaps)

if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt
    splicing_events_dir = None
    CLIP_dir = None
    circRNA_dir =  None
    species = 'Hs'
    extendFlanking = False
    circ_p = 0.05
    
    """ Usage:
    # For finding overlaps between annotated eCLIP and MultiPath-PSI splicing events from AltAnalyze
    import_scripts/peakOverlaps.py --species Hs --clip /Users/abcd/dataAnalysis/eCLIP --events /Users/abcd/dataAnalysis/AS-Events/ --extendFlanking True
    
    # For finding overlaps between annotated eCLIP and circRNA outputs from edgeR
    import_scripts/peakOverlaps.py --species Hs --clip /Users/abcd/dataAnalysis/eCLIP --events /Users/abcd/dataAnalysis/circRNA/ --extendFlanking True  --circ_p 0.1
    
    # For finding overlaps between annotated MultiPath-PSI splicing events and circRNA outputs from edgeR
    python import_scripts/peakOverlaps.py --species Hs --circ /Users/abcd/dataAnalysis/circRNA/ --circ_p 0.1 --events /Users/abcd/dataAnalysis/AS-Events/ --extendFlanking True 
    """

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['species=','clip=','events=', 'extendFlanking=', 'circ=', 'circ_p='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--species':
                species = arg
            elif opt == '--clip':
                CLIP_dir = arg
            elif opt == '--circ':
                circRNA_dir = arg
            elif opt == '--circ_p':
                circ_p = float(arg)
            elif opt == '--events':
                splicing_events_dir = arg
            elif opt == '--extendFlanking':
                if 'true' in string.lower(arg) or 'yes' in string.lower(arg):
                    extendFlanking = True
                    print 'Extend Flanking = TRUE'

    gene_to_symbol, symbol_to_gene = importGeneSymbols(species)
    gene_coordinates, exon_block_coordinates, gene_to_exons = importExonCoordinates(species)
    if CLIP_dir != None:
        dataset_peaks = eCLIPimport(CLIP_dir)
    if splicing_events_dir != None:
        AS_dataset_events = importSplicingEvents(splicing_events_dir)
        dataset_events = AS_dataset_events
        events_dir = splicing_events_dir
    if circRNA_dir != None:
        circRNA_dataset_events = importCircularRNAEvents(circRNA_dir, circ_p)
        dataset_events = circRNA_dataset_events
        events_dir = circRNA_dir
    if CLIP_dir == None:
        """ Comparison of circRNA coordinates as peaks to AS events """
        dataset_events = AS_dataset_events
        dataset_peaks = circRNA_dataset_events
    
    print len(dataset_peaks), 'peak datasets'
    print len(dataset_events), 'Alternative-event datasets\n'
    for eCLIP in dataset_peaks:
        for AE in dataset_events:
            print 'Aligning coordinates from', eCLIP, 'to', AE
            eCLIP_peaks = dataset_peaks[eCLIP]
            events = dataset_events[AE]
            alignEventsAndPeaks(eCLIP, AE, eCLIP_peaks, events,events_dir)
            