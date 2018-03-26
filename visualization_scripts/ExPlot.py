from __future__ import print_function
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy
import export
import traceback
import collections as c
import unique

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>0 and '--' in command_args: commandLine=True
else: commandLine=False

try:
    import math
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        #try: matplotlib.use('TkAgg')
        #except Exception: pass
        #if commandLine==False:
            #try: matplotlib.rcParams['backend'] = 'TkAgg'
            #except Exception: pass
        try:
            import matplotlib.pyplot as pylab
            import matplotlib.colors as mc
            import matplotlib.mlab as mlab
            import matplotlib.ticker as tic
            from matplotlib.patches import Circle
            from mpl_toolkits.mplot3d import Axes3D
            
            matplotlib.rcParams['axes.linewidth'] = 0.5
            matplotlib.rcParams['pdf.fonttype'] = 42
            matplotlib.rcParams['font.family'] = 'sans-serif'
            matplotlib.rcParams['font.sans-serif'] = 'Arial'
            from matplotlib.patches import Rectangle
            import matplotlib.patches
            import matplotlib.backend_bases as event_plot
            from mpldatacursor import datacursor
            from matplotlib.widgets import Slider as matplotSlider
        except Exception:
            print(traceback.format_exc())
            print('Matplotlib support not enabled')
except Exception:
    print(traceback.format_exc())

#os.chdir("/Users/saljh8/Desktop/Code/AltAnalyze/Config/ExonViewFiles")

class EnsemblRegionClass:
    def __init__(self,start,stop,ensembl_exon,exon_region,strand):
        self.start = start
        self.stop = stop
        self.ensembl_exon = ensembl_exon
        self.exon_region = exon_region
        self.strand = strand
    def Start(self): return int(self.start)
    def Stop(self): return int(self.stop)
    def EnsemblExon(self): return self.ensembl_exon
    def EnsemblRegion(self): return self.exon_region
    def ExonBlock(self):
        return string.split(self.EnsemblRegion(),'.')[0]
    def Strand(self): return self.strand
    def Length(self): return abs(int(self.Start())-int(self.Stop()))
    def setChr(self, chr): self.chr = chr
    def Chr(self): return self.chr

class SplicingIndexClass:
    def __init__(self, reg_call, splicing_index, p_val, midas):
        self.reg_call = reg_call
        self.splicing_index = splicing_index
        self.p_val = p_val
        self.midas = midas
    def RegCall(self): return self.reg_call
    def SplicingIndex(self): return self.splicing_index
    def PVal(self): return self.p_val
    def Midas(self): return self.midas

class MicroRNAClass:
    def __init__(self, exon_name, description, basepairs, algorithms):
        self.exon_name = exon_name
        self.description = description
        self.basepairs = basepairs
        self.algorithms = algorithms
    def ExonBlock(self): return self.exon_name
    def Description(self): return self.description
    def BP(self): return self.basepairs
    def Algorithms(self): return self.algorithms
    

def ProteinCentricIsoformView(Selected_Gene):
    Transcript_List = []
    Transcript_db = {}
    Exon_db = {}
    for line in open(Transcript_Annotations_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        if(line[0] == Selected_Gene):
            transcriptID = line[-1]
            exonID = line[5]
            start = line[3]
            stop = line[4]
            strand = line[2]
            chr = line[1]
            if 'chr' not in chr:
                chr = 'chr'+chr
            exon_data = EnsemblRegionClass(start,stop,exonID,None,strand)
            exon_data.setChr(chr)
            Transcript_List.append((transcriptID, exonID))
            try:
                Transcript_db[transcriptID].append(exon_data)
            except Exception:
                Transcript_db[transcriptID]=[exon_data]
            try:
                Exon_db[exonID].append(transcriptID)
            except Exception:
                Exon_db[exonID]=[transcriptID]
    
    Transcript_Protein_db = {}
    Protein_Transcript_db = {}
    Protein_List = []
    count = 0
    for line in open(Prt_Trans_File, "rU").xreadlines():
        if(count == 0):
            count = 1
            continue
        line = line.rstrip()
        line = line.split("\t")
        if(len(line) != 3):
            continue
        geneID = line[0]
        transcriptID = line[1]
        proteinID = line[2]
        if Selected_Gene == geneID:
            Transcript_Protein_db[transcriptID] = proteinID
            Protein_Transcript_db[proteinID] = transcriptID
            Protein_List.append(proteinID)

    #MicroRNA File
    microRNA_db = {}
    for line in open(microRNA_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        try:
            gene_and_exon_id = line[0].split(":")
            current_gene_id = gene_and_exon_id[0]
            current_exon_id = gene_and_exon_id[1]   
        except Exception:
            continue
        #print([current_gene_id,current_exon_id,Selected_Gene]);break
        current_description = line[1]
        current_base_pairs = line[2]
        algorithms = line[3]
        if(current_gene_id == Selected_Gene):
            m = MicroRNAClass(current_exon_id, current_description, current_base_pairs, algorithms)           
            try:
                if(len(microRNA_db[current_exon_id]) > 6):
                    continue
                microRNA_db[current_exon_id].append(m)
                #print("ADDED!")
                
            except:
                microRNA_db[current_exon_id] = [m]
            
    Transcript_ExonRegion_db={}
    geneExonRegion_db={}
    exon_coord_db={}
    exonRegion_db={}
    AllBlocks = [("E", []), ("I", [])]
    # Store the exon region positions and later link them to the Ensembl exons
    for line in open(ExonRegion_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        geneID = line[0]
        exon_region = line[1]
        chr = line[2]
        exonID = line[1]
        strand = line[3]
        start = line[4]
        stop = line[5]
        er = EnsemblRegionClass(start,stop,exonID,exon_region,strand)
        if(geneID == Selected_Gene):
                Block_Num = exon_region[1:]
                I_E_id = exon_region[0]
                if(I_E_id == "E"):
                    AllBlocks[0][1].append(Block_Num)
                if(I_E_id == "I"):
                    AllBlocks[1][1].append(Block_Num)
                    continue
                exon_added = False
                #Exon_List = line[7].split("|")
                exon_coord_db[chr,int(start),'start'] = exon_region
                exon_coord_db[chr,int(stop),'stop'] = exon_region
                exonRegion_db[Selected_Gene,exon_region] = er
                #print chr,start,'start'
    probeset_to_ExonID={}
    if platform != 'RNASeq':
        for line in open(unique.filepath('AltDatabase/'+species+'/'+string.lower(platform)+'/'+species+'_Ensembl_probesets.txt'), "rU").xreadlines():
            line = line.rstrip()
            line = line.split("\t")
            gene = line[2]
            if gene == Selected_Gene:
                probeset = line[0]
                exon_region = line[12]
                if '.' not in exon_region:
                    exon_region = string.replace(exon_region,'-','.')
                probeset_to_ExonID[probeset] = exon_region

    ETC_List = []
    for line in open(SplicingIndex_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        if ':' in line[0]:
            GeneLine = line[0].split(":")
            FeatureID = GeneLine[1]
        else:
            FeatureID = line[0]
        Gene = line[1]
        regcall = line[2]
        spl_index = line[3]
        pval = line[4]
        midas = line[5]
        S_I_data = SplicingIndexClass(regcall, spl_index, pval, midas)
        if(Gene == Selected_Gene):
            if platform != 'RNASeq':
                if FeatureID in probeset_to_ExonID:
                    FeatureID = probeset_to_ExonID[FeatureID]
                    #print(FeatureID)
                    ETC_List.append((FeatureID, S_I_data))
            else:
                try:
                    FeatureID = FeatureID.split("_")
                    FeatureID = FeatureID[0]         
                    ETC_List.append((FeatureID, S_I_data))
                except:
                    pass

    ETC_dict = {}
       
    # Link the exon regions to the Ensembl exons
    for transcriptID in Transcript_db:
        for exon_data in Transcript_db[transcriptID]:
            start = exon_data.Start()
            stop = exon_data.Stop()
            chr = exon_data.Chr()
            strand = exon_data.Strand()
            try:
                start_exon_region = exon_coord_db[chr,start,'start']
                stop_exon_region = exon_coord_db[chr,stop,'stop']
                proceed = True
            except Exception: ### Not clear why this error occurs. Erroring region was found to be an intron region start position (I7.2 ENSMUSG00000020385)
                proceed = False
            if proceed:
                if '-' in strand:
                    stop_exon_region,start_exon_region = start_exon_region,stop_exon_region
                regions = [start_exon_region]
                block,start_region = start_exon_region.split('.')
                start_region = int(float(start_region))
                block,stop_region = stop_exon_region.split('.')
                stop_region = int(float(stop_region))
                region = start_region+1
                while region<stop_region:
                    er = block+'.'+str(region)
                    regions.append(er)
                    region+=1
                if stop_region != start_region:
                    regions.append(stop_exon_region)
                for region in regions:
                    er = exonRegion_db[Selected_Gene,region]
                    try:
                        Transcript_ExonRegion_db[transcriptID].append(er)
                    except:
                        Transcript_ExonRegion_db[transcriptID] = [er]
    
    exon_virtualToRealPos= c.OrderedDict()
    junction_transcript_db = {}
    for transcriptID in Transcript_ExonRegion_db:
            #print('transcripts:',transcriptID)
            position=0
            Buffer=15
            for exon_object in Transcript_ExonRegion_db[transcriptID]:
                if position!=0:
                    if last_exon != exon_object.ExonBlock():
                        #print last_exon_region+'-'+exon_object.EnsemblRegion(),position,Buffer
                        junctionID = last_exon_region+'-'+exon_object.EnsemblRegion()
                        try: junction_transcript_db[transcriptID].append((position,position+Buffer, junctionID)) ### virtual junction positions
                        except: junction_transcript_db[transcriptID] = [(position,position+Buffer, junctionID)]
                        position+=Buffer
                        
                virtualStart = position
                virtualStop = virtualStart + exon_object.Length()
                position = virtualStop
                try:                    
                    exon_virtualToRealPos[transcriptID].append(([virtualStart,virtualStop],[exon_object.Start(), exon_object.Stop()],exon_object))
                except Exception:                    
                    exon_virtualToRealPos[transcriptID]=[([virtualStart,virtualStop],[exon_object.Start(), exon_object.Stop()],exon_object)]
                #print transcriptID,exon_object.ExonBlock(),exon_object.EnsemblExon(),exon_object.EnsemblRegion(),exon_object.Start(),exon_object.Stop(),virtualStart,virtualStop,"\n"
                last_exon = exon_object.ExonBlock()
                last_exon_region = exon_object.EnsemblRegion()

    for i in ETC_List:
        Region = i[0]
        S_I = i[1]
        Region = Region.split("-")
        if(len(Region) > 1):

            #Delete underscores from values.

            R_Start = Region[0]
            R_End = Region[1]
            R_Start = R_Start.split("_")
            R_End = R_End.split("_")
            R_Start = R_Start[0]
            R_End = R_End[0]
            R_Final = R_Start + "-" + R_End
            R_Type = R_Final[0]
            #print(R_Final)
            ETC_dict[R_Final] = S_I
            
        else:
            Region = Region[0]
            Region = Region.split("_")
            Region = Region[0]
            Region_type = Region[0]
            ETC_dict[Region] = S_I
            #if(Region_type == "E"):
            #    for entry in AllBlocks[0][1]:
            #        if(Region[1:] == entry):
            #            ETC_dict[("E" + entry)] = S_I
            #if(Region_type == "I"):
            #    for entry in AllBlocks[1][1]:    
            #        if(Region[1:] == entry):
            #            ETC_dict[("I" + entry)] = S_I

    #for a in ETC_dict:
    #        print(ETC_dict[a].RegCall(), a)
          
    #for i in junction_transcript_db:
    #    print i, junction_transcript_db[i], "\n"
    
    Protein_Pos_Db = {}
    last_protein=None
    stored_stop=None
    for line in open(Prt_Boundaries_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        proteinID = line[0]
        if(proteinID in Protein_List):
            Stop = int(line[-1])
            Start = int(line[-2])
            if(proteinID != last_protein):
                if stored_stop !=None:
                    #print proteinID,stored_start,stored_stop
                    Protein_Pos_Db[last_protein] = [[stored_start,stored_stop,None]]
                stored_start = int(Start)
            if(proteinID == last_protein):
                stored_stop = int(Stop)
            last_protein = str(proteinID)
    
    Protein_Pos_Db[last_protein] = [(stored_start,stored_stop,None)]
    Protein_virtualPos = RealToVirtual(Protein_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db)
    
    Domain_Pos_Db={}
    domainAnnotation_db={}
    #"""
    for line in open(Prt_Regions_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        proteinID = line[0]
        if proteinID in Protein_Pos_Db:
            domain_start = int(float(line[3]))
            domain_stop = int(float(line[4]))
            domainID = line[-2]
            domainName = line[-1]
            try:
                Domain_Pos_Db[proteinID].append((domain_start,domain_stop,domainID))
            except:
                Domain_Pos_Db[proteinID] = [(domain_start,domain_stop,domainID)]
            domainAnnotation_db[domainID] = domainName

    #"""
    for line in open(UniPrt_Regions_File, "rU").xreadlines():
        line = line.rstrip()
        line = line.split("\t")
        proteinID = line[0]
        if proteinID in Protein_Pos_Db:
            domain_start = int(float(line[3]))
            domain_stop = int(float(line[4]))
            domainID = line[-1]
            domainName = line[-1]
            try:
                Domain_Pos_Db[proteinID].append((domain_start,domain_stop,domainID))
            except:
                Domain_Pos_Db[proteinID] = [(domain_start,domain_stop,domainID)]
            domainAnnotation_db[domainID] = domainName
            #print('--',domainName,domain_start,domain_stop)

    # Do the same for domain coordinates
    Domain_virtualPos = RealToVirtual(Domain_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db)

    return_val = ((junction_transcript_db, Protein_virtualPos, Domain_virtualPos, Transcript_db, exon_virtualToRealPos, ETC_dict, microRNA_db, domainAnnotation_db))
    return return_val
    
def RealToVirtual(Protein_Pos_Db, exon_virtualToRealPos, Protein_Transcript_db,Transcript_ExonRegion_db):
    Transcript_to_Protein_Coords = {}
    for proteinID in Protein_Pos_Db:
        transcript = Protein_Transcript_db[proteinID]
        strand = Transcript_ExonRegion_db[transcript][0].Strand()
        e_coords = exon_virtualToRealPos[transcript]
        if proteinID in Protein_Pos_Db:
            for q in Protein_Pos_Db[proteinID]:
                #print("Protein: ", proteinID)
                p_start = q[0]
                p_stop = q[1]
                annotation = q[2]
                if '-' not in strand:
                    p_start +=1
                    p_stop -= 1 ### seems to be off by 1
                virtual_p_start = None
                virtual_p_stop = None
                #print("E", len(e_coords))
                #print("Protein: ", proteinID)
                
                for i in range(len(e_coords)):
                    e_virtual_start, e_virtual_stop = e_coords[i][0]
                    #print("Sub-E", e_virtual_start, e_virtual_stop)
                    e_real_start,e_real_stop = e_coords[i][1]
                    #print e_real_start,e_real_stop
                    e = [e_real_start,e_real_stop]
                    e.sort()
                    p = [p_start,p_stop]
                    p.sort()
                    coord = e+p
                    coord.sort()
                    if (p_start<e[1] and p_start>e[0]) or p_start==e[1] or p_start==e[0]:
                        if '-' in strand:
                            offset = e_real_stop-p_start
                        else:
                            offset = p_start-e_real_start
                        virtual_p_start = offset+e_virtual_start
                        #print("Final_Val", proteinID, virtual_p_start)
                    if (p_stop<e[1] and p_stop>e[0]) or p_stop==e[1] or p_stop==e[0]:
                        if '-' in strand:
                            offset = e_real_stop-p_stop
                        else:
                            offset = p_stop-e_real_start
                        virtual_p_stop = offset+e_virtual_start
                if annotation != None:
                    #print("Entered", proteinID, virtual_p_start)
                    try:
                        Transcript_to_Protein_Coords[transcript].append((proteinID, annotation, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1]))
                    except Exception:
                        Transcript_to_Protein_Coords[transcript] = [(proteinID, annotation, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1])]
                else:
                    #print("Entered_2", proteinID, virtual_p_start)
                    Transcript_to_Protein_Coords[transcript] = proteinID, virtual_p_start, virtual_p_stop, e_coords[0][0][0],e_coords[-1][0][1]
                #print transcript, proteinID, virtual_p_start, virtual_p_stop, p_start,p_stop, e_coords[0][0][0],e_coords[-1][0][1],annotation
    return Transcript_to_Protein_Coords
     
def searchDirectory(directory,var,secondary=None):
    directory = unique.filepath(directory)

    files = unique.read_directory(directory)
    for file in files:
        if var in file:
            if secondary== None:
                return directory+'/'+file
                break
            elif secondary in file:
                return directory+'/'+file
                break
            
    ### if all else fails
    return directory+'/'+file 
    
def getPlatform(filename):
    prefix = string.split(export.findFilename(filename),'.')[0]
    array_type = string.split(prefix,'_')[1]
    if array_type != 'RNASeq':
        array_type = string.lower(array_type)
    return array_type

def remoteGene(gene,Species,root_dir,comparison_file):
    global Transcript_Annotations_File
    global ExonRegion_File
    global Selected_Gene
    global Prt_Trans_File
    global Prt_Regions_File
    global Prt_Boundaries_File
    global SplicingIndex_File
    global UniPrt_Regions_File
    global microRNA_File
    global domainAnnotation_db
    global platform
    global species

    Selected_Gene = str(gene)
    species = Species
    
    comparison_name = string.split(export.findFilename(comparison_file),'.')[0]
    ExonRegion_File = unique.filepath("AltDatabase/ensembl/"+species+"/"+species+"_Ensembl_exon.txt")
    Transcript_Annotations_File = unique.filepath("AltDatabase/ensembl/"+species+"/"+species+"_Ensembl_transcript-annotations.txt")
    Prt_Trans_File = searchDirectory("AltDatabase/ensembl/"+species+"/",'Ensembl_Protein')
    Prt_Regions_File = searchDirectory("AltDatabase/ensembl/"+species+"/",'ProteinFeatures')
    Prt_Boundaries_File = searchDirectory("AltDatabase/ensembl/"+species+"/",'ProteinCoordinates')
    UniPrt_Regions_File = searchDirectory("AltDatabase/uniprot/"+species+"/",'FeatureCoordinate')
    SplicingIndex_File = searchDirectory(root_dir+'/AltResults/ProcessedSpliceData/','splicing-index',secondary=comparison_name)
    platform = getPlatform(SplicingIndex_File)
    microRNA_File = searchDirectory("AltDatabase/"+species+"/"+platform,'microRNAs_multiple')
    #print(SplicingIndex_File)

    total_val = ProteinCentricIsoformView(Selected_Gene)
    junctions = total_val[0]
    p_boundaries = total_val[1]
    p_domains = total_val[2]
    transcript_db = total_val[3]
    exon_db = total_val[4]
    splice_db = total_val[5]
    microRNA_db = total_val[6]
    domainAnnotation_db = total_val[7]

    #for i in exon_db:
    #    print("THE", i, exon_db[i], "\n")

    #for i in microRNA_db:
    #        m_test = microRNA_db[i]
    #    print(len(m_test))
    #    for q in m_test:
    #        print("microRNA", q.ExonBlock(), q.Description(), q.BP(), "\n")

    #for i in exon_db["ENST00000349238"]:
    #    print(i[2].EnsemblRegion())
    
    domain_color_list = []
    for i in p_domains:
        ploy = p_domains[i]
        for a in ploy:
            domain_color_list.append(a[1])

    domain_color_list = list(set(domain_color_list))
    domain_color_key = {}
    c_color1 = [0.8, 0.6, 0.1]
    c_color2 = [0.1, 0.6, 0.8]
    c_color3 = [0.6, 0.1, 0.8]
    c_color4 = [0.95, 0.6, 0.3]
    c_color5 = [0.3, 0.6, 0.95]
    c_color6 = [0.6, 0.3, 0.95]
    FLAG = 1

    for item in domain_color_list:
        if(FLAG == 1):
            domain_color_key[item] = c_color1
            FLAG = FLAG + 1
            continue
        if(FLAG == 2):
            domain_color_key[item] = c_color2
            FLAG = FLAG + 1
            continue
        if(FLAG == 3):
            domain_color_key[item] = c_color3
            FLAG = FLAG + 1
            continue
        if(FLAG == 4):
            domain_color_key[item] = c_color4
            FLAG = FLAG + 1
            continue
        if(FLAG == 5):
            domain_color_key[item] = c_color5
            FLAG = FLAG + 1
            continue
        if(FLAG == 6):
            domain_color_key[item] = c_color6
            FLAG = 1
            continue

    #for i in domain_color_key:
        #print(i, domain_color_key[i], "\n")
    
    Y = 100
    Transcript_to_Y = {}
    for transcript in transcript_db:
        Transcript_to_Y[transcript] = Y
        Y = Y + 300
    import traceback

    def onpick(event):
        #ind = event.ind
        print(event.artist.get_label())

    #for i in domainAnnotation_db: print(i,len(domainAnnotation_db));break
    
    fig = pylab.figure()
    
    ylim = Y + 200
    currentAxis = pylab.gca()
    #ax = pylab.axes()
    ax = fig.add_subplot(111)
    X_Pos_List = []
    CoordsBank = []
    
    for transcript in transcript_db:
        try:
            Junc_List = junctions[transcript]
            y_pos = Transcript_to_Y[transcript]
            Gene_List = exon_db[transcript]
            color_flag = 1
            for entry in Gene_List:
                G_start = entry[0][0]
                G_end = entry[0][1]
                Exon_Object = entry[2]
                try:
                    LabelClass = splice_db[Exon_Object.EnsemblRegion()]
                    ExonName = Exon_Object.EnsemblExon()
                    RegCall = LabelClass.RegCall()
                    SplicingIndex = LabelClass.SplicingIndex()
                    PVal = LabelClass.PVal()
                    Midas = LabelClass.Midas()
                    Label = "\n" + "Exon: " + str(ExonName) + "\n" + "RegCall: "  + str(RegCall) + "\n" + "Splicing Index: " + str(SplicingIndex) + "\n" + "P-Value: " + str(PVal) + "\n" + "Midas Value: " + str(Midas) + "\n"
                    Label = string.replace(Label,"\n"," ")
                    if(RegCall == "UC"):
                        color_choice = "Grey"
                    else:
                        S_Int = float(SplicingIndex)
                        if(S_Int > 0):
                            #color_choice = (0.7, 0.7, 0.99)
                            color_choice = 'blue'
                        if(S_Int < 0):
                            #color_choice = (0.8, 0.4, 0.4)
                            color_choice = 'red'
                                            
                except:
                    #print(traceback.format_exc());sys.exit()
                    Label = ""
                    color_choice = "Grey"
                #print("Start", G_start, "end", G_end, "Region", entry[2].EnsemblRegion())
                if((color_flag % 2) == 0):
                    currentAxis.add_patch(Rectangle((G_start, y_pos), (G_end - G_start), 50, color = color_choice, label = (entry[2].EnsemblRegion() + Label), picker = True))
                    y_end = y_pos + 50
                    try: CoordsBank.append((G_start, G_end, y_pos, y_end, 'Exon: '+entry[2].EnsemblRegion()+' '+ 'SI: '+str(SplicingIndex)[:4]+' Pval: '+str(Midas)[:4]))
                    except Exception:
                        CoordsBank.append((G_start, G_end, y_pos, y_end, 'Exon: '+entry[2].EnsemblRegion()))
                    #print(entry[2].EnsemblRegion(),y_pos,y_end)
                if((color_flag % 2) != 0):                   
                    currentAxis.add_patch(Rectangle((G_start, y_pos), (G_end - G_start), 50, color = color_choice, label = (entry[2].EnsemblRegion() + Label), picker = True))
                    y_end = y_pos + 50
                    try: CoordsBank.append((G_start, G_end, y_pos, y_end, 'Exon: '+entry[2].EnsemblRegion()+' '+ 'SI: '+str(SplicingIndex)[:4]+' p-value: '+str(Midas)[:4]))
                    except Exception:
                        CoordsBank.append((G_start, G_end, y_pos, y_end, 'Exon: '+entry[2].EnsemblRegion()))
                    #print(entry[2].EnsemblRegion(),y_pos,y_end)
                color_flag = color_flag + 1
                if(entry[2].EnsemblRegion() in microRNA_db):
                    microRNA_object = microRNA_db[entry[2].EnsemblRegion()]
                    mr_label = "MICRORNA MATCHES" + "\n"
                    for class_object in microRNA_object:
                        mr_exonname = class_object.ExonBlock()
                        mr_desc = class_object.Description() + " " + class_object.Algorithms()
                        #print(mr_desc)
                        mr_label = mr_label + mr_desc + "\n"
                    
                    currentAxis.add_patch(Rectangle((G_start, (y_pos - 75)), (G_end - G_start), 40, color = "Green", label = (mr_label), picker = True))
                    y_start = y_pos - 75
                    y_end = y_pos - 35
                    CoordsBank.append((G_start, G_end, y_start, y_end, mr_desc))
                
            for entry in Junc_List:
                junctionID = entry[-1]
                try:
                    LabelClass = splice_db[entry[2]]
                    RegCall = LabelClass.RegCall()
                    SplicingIndex = LabelClass.SplicingIndex()
                    PVal = LabelClass.PVal()
                    Midas = LabelClass.Midas()
                    Label = "\n" + "RegCall: " + str(RegCall) + "\n" + "Splicing Index: " + str(SplicingIndex) + "\n" + "P-Value: " + str(PVal) + "\n" + "Midas Value: " + str(Midas) + "\n"
                    if(float(SplicingIndex) > 0):
                        color_junc = "blue"
                    if(float(SplicingIndex) < 0):
                        color_junc = "red"
                    if(RegCall == "UC"):
                        color_junc = "grey"
                except:
                    Label = ""
                    color_junc = "grey"
                currentAxis.add_patch(Rectangle((entry[0], y_pos), (entry[1] - entry[0]), 50, color = "White", label = (str(entry[2]) + Label), picker = True))
                ax.arrow(entry[0], (y_pos+50), 8, 40, label = (str(entry[2]) + Label), color = color_junc, picker = True)
                ax.arrow((entry[0] + 8), (y_pos+90), 11, -40, label = (str(entry[2]) + Label), color = color_junc, picker = True)
                y_start = y_pos
                y_end = y_pos + 30
                #print(junctionID,y_start,y_end)
                CoordsBank.append((G_start, G_end, y_start, y_end, junctionID))

            try:
                P_Bound_List = p_boundaries[transcript]
                E_Start = P_Bound_List[-2]
                E_End = P_Bound_List[-1]
                P_Start = P_Bound_List[1]
                P_End = P_Bound_List[2]
                #print("Boundaries: ", P_Start, P_End)
                X_Pos_List.append(int(E_End))
                #currentAxis.add_patch(Rectangle((E_Start, y_pos), E_End, 50, color = "Blue"))
                try:
                    currentAxis.add_patch(Rectangle((P_Start, (y_pos + 120)), (P_End - P_Start), 10))
                except:
                    pass
                p_label_list = ["DEF"]
                #CoordsBank.append((P_Start, P_End, y_pos, P_End - P_Start, transcript)) ### Added by NS - needs work
                try: P_Domain_List = p_domains[transcript]
                except Exception: P_Domain_List=[]
                for entry in P_Domain_List:
                    #print("Domain", entry)
                    color_domain_choice = domain_color_key[entry[1]]
                    domain_annotation = domainAnnotation_db[entry[1]]
                    #domain_annotation = string.replace(domain_annotation,'REGION-','')
                    p_label = (str(entry[0]) +  " " + str(domain_annotation))
                    #print(entry[0], entry[2], entry[3], P_Start, P_End, domain_annotation, )
                    Repeat_Flag = 0
                    for i in p_label_list:
                        if(p_label == i):
                            Repeat_Flag = 1
                    if(Repeat_Flag == 1):
                        continue
                    p_label_list.append(p_label)               
                    currentAxis.add_patch(Rectangle((entry[2], y_pos + 100), (entry[3] - entry[2]), 50, color = color_domain_choice, label= p_label, picker = True))
                    y_start = y_pos + 100
                    y_end = y_pos + 150
                    CoordsBank.append((entry[2], entry[3], y_start, y_end, p_label))
            except Exception:
                pass
                #print(traceback.format_exc())
        except:
            #print(traceback.format_exc())
            pass
    pylab.ylim([0.0, ylim])
    try:
        max_x = max(X_Pos_List)
    except:
        max_x = 5000
    try:
        pylab.xlim([0.0, max_x])
    except:
        pylab.xlim([0.0, 3000])
    fig.canvas.mpl_connect('pick_event', onpick)
    def format_coord(x, y):
        for m in CoordsBank:
            if(x >= m[0] and x <= m[1] and y >= m[2] and y <= m[3]):
                string_display = m[4]
                return string_display
        string_display = "  "
        return string_display

    ax.format_coord = format_coord
    #datacursor(hover=True, formatter='{label}'.format, bbox=dict(fc='yellow', alpha=1), arrowprops=None)
    pylab.show()
    
if __name__ == "__main__":
    #Selected_Gene = sys.argv[1]
    Selected_Gene = 'ENSG00000132906'
    Species = 'Hs'
    root_dir = '/Volumes/salomonis2/Leucegene_project/STAR_TOP-SRSF2-U2AF1-like/'
    comparison_file = '/Volumes/salomonis2/Leucegene_project/STAR_TOP-SRSF2-U2AF1-like/AltResults/RawSpliceData/Hs/splicing-index/Hs_RNASeq_U2AF1-like_vs_SRSF2-like.ExpCutoff-5.0_average.txt'
    remoteGene(Selected_Gene,Species,root_dir,comparison_file)