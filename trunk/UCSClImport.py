import sys, string
import os.path
import unique
import EnsemblImport
import copy
import time
import alignToKnownAlt

dirfile = unique
py2app_adj = '/AltAnalyze.app/Contents/Resources/Python/site-packages.zip'

def filepath(filename):
    dir=os.path.dirname(dirfile.__file__)       #directory file is input as a variable under the main            
    fn=os.path.join(dir,filename)
    fn = string.replace(fn,py2app_adj,'')
    fn = string.replace(fn,'\\library.zip','') ###py2exe on some systems, searches for all files in the library file, eroneously
    return fn

def read_directory(sub_dir):
    dirfile = unique
    dir=os.path.dirname(dirfile.__file__)
    dir = string.replace(dir,py2app_adj,'')
    dir = string.replace(dir,'\\library.zip','')
    dir_list = os.listdir(dir + sub_dir)
    return dir_list

################### Import exon coordinate/transcript data from BIOMART
def eliminate_redundant_dict_values(database):
    db1={}
    for key in database:
        list = unique.unique(database[key])
        list.sort()
        db1[key] = list
    return db1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def covertStrToListInt(str_data):
    str_data = string.replace(str_data,'"','')
    list = string.split(str_data,',');list2=[]
    try:
        for i in list[:-1]: list2.append(int(i))
    except ValueError: print list;kill
    return list2

def importEnsExonStructureData(species):
    start_time = time.time(); global ensembl_annotations; global ensembl_exon_data; ensembl_exon_data={}; global ensembl_gene_coordinates
    global ensembl_gene_exon_db;  ensembl_gene_exon_db={}; global ensembl_exon_pairs; global coordinate_gene_count_db
    global ensembl_transcript_structure_db; ensembl_transcript_structure_db={}; global ens_transcript_gene_db; ens_transcript_gene_db={}
    global ensembl_const_exon_db; ensembl_const_exon_db = {}
    ###Simple function to import and organize exon/transcript data
    if round == 1: filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    if round == 2: filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_NonGene_structure_'+mRNA_type+'.txt'
    fn=filepath(filename); ensembl_gene_coordinates={}; ensembl_annotations={}; x=0
    transcript_exon_db = {}; initial_junction_db = {}; ensembl_transcript_exons={}; coordinate_gene_count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            if 'Chromosome' in t[0]:  type = 'old'###when using older builds of EnsMart versus BioMart
            else: type = 'current'
            x=1
        else:
            if type == 'old': chr, strand, gene, ens_transcriptid, ens_exonid, exon_start, exon_end, constitutive_exon = t
            else: gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            ###switch exon-start and stop if in the reverse orientation
            if strand == '-1': strand = '-'
            else: strand = '+'
            if '_' in chr: c = string.split(chr,'_'); chr = c[0]
            exon_end = int(exon_end); exon_start = int(exon_start)
            ens_exonid_data = string.split(ens_exonid,'.'); ens_exonid = ens_exonid_data[0]
            if test == 'yes':
                if gene in test_gene: proceed = 'yes'
                else: proceed = 'no'
            else: proceed = 'yes'
            if abs(exon_end-exon_start)>0 and proceed == 'yes':
                
                ###Create temporary databases storing just exon and just trascript and then combine in the next block of code
                try: ensembl_gene_coordinates[gene].append(exon_start)
                except KeyError: ensembl_gene_coordinates[gene] = [exon_start]
                try: coordinate_gene_count_db[gene].append([exon_start,exon_end])
                except KeyError: coordinate_gene_count_db[gene] = [[exon_start,exon_end]]
                ensembl_gene_coordinates[gene].append(exon_end)
                ensembl_annotations[gene] = chr,strand
                if analysis_type == 'mRNA': ### This conserves energy
                    ensembl_exon_data[(chr,exon_start,exon_end)] = ens_exonid
                    try: ensembl_gene_exon_db[gene]+= [ens_exonid]
                    except KeyError: ensembl_gene_exon_db[gene] = [ens_exonid]
                    try: ensembl_transcript_exons[ens_transcriptid,strand].append((exon_start,ens_exonid))
                    except KeyError: ensembl_transcript_exons[ens_transcriptid,strand] = [(exon_start,ens_exonid)]
                    try: ensembl_transcript_structure_db[ens_transcriptid].append((chr,exon_start,exon_end))
                    except KeyError: ensembl_transcript_structure_db[ens_transcriptid] = [(chr,exon_start,exon_end)]
                    ens_transcript_gene_db[ens_transcriptid] = gene
                    ensembl_const_exon_db[ens_exonid] = constitutive_exon
                
    end_time = time.time(); time_diff = int(end_time-start_time)
    print filename,"parsed in %d seconds" % time_diff

    if analysis_type == 'ncRNA':
        ### Delete dbase entries where there are multiple exons
        coordinate_gene_count_db = delMultipleEntries(coordinate_gene_count_db)
        ensembl_transcript_exons = delMultipleEntries(ensembl_transcript_exons)
        ensembl_transcript_structure_db = delMultipleEntries(ensembl_transcript_structure_db)
        
        ensembl_gene_coordinates = filterMultipleEntries(coordinate_gene_count_db,ensembl_gene_coordinates)
        ensembl_gene_exon_db = filterMultipleEntries(coordinate_gene_count_db,ensembl_gene_exon_db)

    ###Sort exon info in the transcript to store pairs of Ensembl exons to find unique pairs of exons from UCSC
    ensembl_exon_pairs={}
    for (transcript,strand) in ensembl_transcript_exons:
        a = ensembl_transcript_exons[(transcript,strand)]; a.sort()
        if strand == '-': a.reverse()
        index=0
        try:
            while index<len(a):
                exon_pair = a[index][1],a[index+1][1]
                ensembl_exon_pairs[exon_pair]=[]; index+=1
        except IndexError: index+=1
        
    ensembl_chr_coordinate_db={}
    for gene in ensembl_gene_coordinates:
        a = ensembl_gene_coordinates[gene]; a.sort()
        gene_start = a[0]; gene_stop = a[-1]
        chr,strand = ensembl_annotations[gene]
        if chr in ensembl_chr_coordinate_db:
            ensembl_gene_coordinates2 = ensembl_chr_coordinate_db[chr]
            ensembl_gene_coordinates2[(gene_start,gene_stop)] = gene,strand
        else:
            try: ensembl_gene_coordinates2={}; ensembl_gene_coordinates2[(gene_start,gene_stop)] = gene,strand
            except TypeError: print gene,ensembl_gene_coordinates[gene];kill
            ensembl_chr_coordinate_db[chr]=ensembl_gene_coordinates2
    return ensembl_chr_coordinate_db

def filterMultipleEntries(db1,db2):
    db3 = {}
    for i in db1:
        if i in db2: db3[i]=db2[i]
    return db3

def delMultipleEntries(db):
    db1 = {}
    for key in db:
        if len(db[key])==1: db1[key]=db[key]

    db2 = {}
    for key in db1:
        vals = db1[key][0]
        if len(vals) == 2:
            try:
                exon_length = int(vals[1]) - int(vals[0])
                if exon_length < 1350: db2[key]=db1[key]
            except ValueError: db2[key]=db1[key]
        else: db2[key]=db1[key]
    return db2

def makeUnique(item):
    db1={}; list1=[]
    for i in item: db1[i]=[]
    for i in db1: list1.append(i)
    list1.sort()
    return list1

def mergeFragmentedExons(exon_coordiantes,strand):
    ###If the intron between two exons is < 10bp, make one exon
    index=0; merged='no'
    #print len(exon_coordiantes),exon_coordiantes,strand,'\n'
    while index<(len(exon_coordiantes)-1):
        deleted = 'no'
        if strand == '-':
            (stop1,start1) = exon_coordiantes[index]
            (stop2,start2) = exon_coordiantes[index+1]
        else:
            (start1,stop1) = exon_coordiantes[index]
            (start2,stop2) = exon_coordiantes[index+1]
        if abs(start2-stop1)<9:
            if strand == '+':  new_exon = (start1,stop2)
            else:  new_exon = (stop2,start1)###for neg_strand
            exon_coordiantes[index+1] = new_exon; del exon_coordiantes[index]
            deleted = 'yes'; merged = 'yes'
        index+=1
        if deleted == 'yes':
            if index == (len(exon_coordiantes)):break
            else: index-=1 ###reset this since the number of elements has changed
    exon_coordiantes = makeUnique(exon_coordiantes)
    if strand == '-': exon_coordiantes.reverse()
    #print len(exon_coordiantes),exon_coordiantes,strand;kill
    return exon_coordiantes,merged
    
def importUCSCExonStructureData(species):
    start_time = time.time(); global ucsc_annotations; global input_gene_file; global ucsc_genomic_accession_db
    ###Simple function to import and organize exon/transcript data
    if analysis_type == 'mRNA': filename = 'AltDatabase/ucsc/'+species+'/all_'+mRNA_type+'.txt'; input_gene_file = filename; merged_ac_count=0;ac_count=0
    elif analysis_type == 'ncRNA': filename = 'AltDatabase/ucsc/'+species+'/all_ncRNA.txt'; input_gene_file = filename; merged_ac_count=0;ac_count=0
    fn=filepath(filename); ucsc_gene_coordinates={}; transcript_structure_db={}; ucsc_interim_gene_transcript_db={}
    ucsc_annotations={}; ucsc_transcript_coordinates={}; temp_data_db={}; accession_count={}; clusterid_count={}
    ucsc_genomic_accession_db = {}; last_exon_stop=0; last_exon_start=0
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if analysis_type == 'mRNA':
            coordinates = t[-1]; coordinates = covertStrToListInt(coordinates)
            exon_size_list = t[-3]; exon_size_list = covertStrToListInt(exon_size_list) ##coordinates and are the same length
            strand = t[9]; accession = t[10]; clusterid = t[0]; chr = t[14][3:]
            index=0; exon_coordiantes = []
            
            while index<len(coordinates):
                exon_start = coordinates[index]; exon_stop = coordinates[index]+exon_size_list[index]
                exon_coordiantes.append((exon_start+1,exon_stop))
                try: ###Below code is available if we want to join exons that have very small introns (really one exon)... however, EnsemblImport will agglomerate these exons and ignore the splice event (if in the same exon)
                    if (coordinates[index+1]-exon_stop)<gap_length: null = []#print accession, exon_stop,coordinates;kill
                except IndexError: null=[]
                index+=1
            if strand == '-': exon_coordiantes.reverse()
            exon_coordiantes,merged = mergeFragmentedExons(exon_coordiantes,strand)
            merged='no'
            if merged=='yes': merged_ac_count+=1
            unique_geneid = clusterid,strand,chr  ###can have strand specific isoforms
        elif analysis_type == 'ncRNA':
            exon_coordiantes = []
            try: clusterid,chr,exon_start,exon_stop,name,score,strand,null,null,nc_class = t
            except ValueError: print t;kill
            chr = string.replace(chr,'chr','')
            accession = name+"_"+str(exon_start)
            exon_start = int(exon_start); exon_stop = int(exon_stop)
            exon_coordiantes.append((exon_start+1,exon_stop))
            unique_geneid = accession,strand,chr  ###can have strand specific isoforms
        if '_' in chr: c = string.split(chr,'_'); chr = c[0]
        if exon_start == last_exon_start or last_exon_stop == exon_stop: null=[] ### Don't record entry - only applies to ncRNA data
        elif len(exon_coordiantes)>1: null=[]
        elif (exon_stop - exon_start)>1500: null=[]
        else:   
            #ucsc_annotations[chr,start,end]=[]
            temp_data_db[accession] = unique_geneid,strand,exon_coordiantes,chr
            ucsc_genomic_accession_db[accession]=chr,strand
            try: accession_count[accession]+=1
            except KeyError: accession_count[accession]=1
            clusterid_count[unique_geneid] = []
            ac_count+=1
        last_exon_start = exon_start; last_exon_stop = exon_stop
    print len(clusterid_count), 'clusters imported'
    print merged_ac_count, "transcripts had at least two exons merged into one, out of",ac_count

    for accession in temp_data_db:
        if accession_count[accession]==1: ###Why would an accession have multiple entries: shouldn't but does
            unique_geneid,strand,exon_coordiantes,chr = temp_data_db[accession]
            ###Add the first and list positions in the transcript
            ucsc_annotations[unique_geneid] = chr ###don't include strand, since transcripts in a cluster can be on different strands
            try: ucsc_gene_coordinates[unique_geneid]+=[exon_coordiantes[0][0]]
            except KeyError: ucsc_gene_coordinates[unique_geneid]=[exon_coordiantes[0][0]]
            ucsc_gene_coordinates[unique_geneid]+=[exon_coordiantes[-1][1]]

            try: ucsc_transcript_coordinates[unique_geneid]+=[exon_coordiantes[0][0]]
            except KeyError: ucsc_transcript_coordinates[accession]=[exon_coordiantes[0][0]]
            ucsc_transcript_coordinates[accession]+=[exon_coordiantes[-1][1]]
            
            try: ucsc_interim_gene_transcript_db[unique_geneid].append(accession)
            except KeyError: ucsc_interim_gene_transcript_db[unique_geneid] = [accession]
            transcript_structure_db[accession] = exon_coordiantes
            temp_data_db[accession]
    end_time = time.time(); time_diff = int(end_time-start_time)
    print filename,"parsed in %d seconds" % time_diff
    
    ###Build a gene cluster level database (start and stop coordinates) for UCSC clusters
    ucsc_chr_coordinate_db={}
    for geneid in ucsc_gene_coordinates:
        strand = geneid[1]
        a = ucsc_gene_coordinates[geneid]; a.sort()
        gene_start = a[0]; gene_stop = a[-1]
        chr = ucsc_annotations[geneid]
        if chr in ucsc_chr_coordinate_db:
            ucsc_gene_coordinates2 = ucsc_chr_coordinate_db[chr]
            ucsc_gene_coordinates2[(gene_start,gene_stop)] = geneid,strand
        else:
            ucsc_gene_coordinates2={}; ucsc_gene_coordinates2[(gene_start,gene_stop)] = geneid,strand
            ucsc_chr_coordinate_db[chr] = ucsc_gene_coordinates2
    
    return ucsc_chr_coordinate_db,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates

def getChromosomalOveralap(ucsc_chr_db,ensembl_chr_db):
    print len(ucsc_chr_db),len(ensembl_chr_db); start_time = time.time()
    """Find transcript_clusters that have overlapping start positions with Ensembl gene start and end (based on first and last exons)"""
    ###exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
    y = 0; l =0; ensembl_transcript_clusters={}; no_match_list=[]
    ###(bp1,ep1) = (47211632,47869699); (bp2,ep2)  =  (47216942, 47240877)
    for chr in ucsc_chr_db:
        ucsc_db = ucsc_chr_db[chr]
        try:
         for (bp1,ep1) in ucsc_db:
            x = 0
            gene_clusterid,ucsc_strand = ucsc_db[(bp1,ep1)]
            accession = gene_clusterid[0]
            try:
                ensembl_db = ensembl_chr_db[chr]
                delta = 0
                for (bp2,ep2) in ensembl_db:
                    y += 1; ensembl,ens_strand = ensembl_db[(bp2,ep2)]
                    """if accession == 'n3578(miRNA)_20337643':
                        print ensembl,bp1,ep1,bp2,ep2
                        delta = 1"""
                    if ucsc_strand == ens_strand:
                        ###if the two gene location ranges overlapping
                        ##########FORCE UCSC mRNA TO EXIST WITHIN THE SPACE OF ENSEMBL TO PREVENT TRANSCRIPT CLUSTER EXCLUSION IN ExonArrayEnsemblRules
                        add = 0
                        if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((ep1 >= bp2)  and (ep2 >= ep1)): add = 1 ###if the start or stop of the UCSC region is inside the Ensembl start and stop
                        elif ((bp2 >= bp1)  and (ep1 >= bp2)) or ((ep2 >= bp1)  and (ep1 >= ep2)): add = 1 ###opposite
                        if add == 1:
                            #if ((bp1 >= (bp2-bp_offset))  and ((ep2+bp_offset) >= ep1)): a = ''
                            #else: print gene_clusterid,ensembl,bp1,bp2,ep1,ep2;kill
                            x = 1
                            try: ensembl_transcript_clusters[gene_clusterid].append(ensembl)
                            except KeyError: ensembl_transcript_clusters[gene_clusterid] = [ensembl]
                            l += 1
            except KeyError: null=[]; #print chr, 'not found'
            if x == 0: no_match_list.append(gene_clusterid)
        except ValueError:
            for y in ucsc_db: print y;kill
    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC genes matched up to Ensembl in %d seconds" % time_diff            
    print "UCSC Transcript Clusters (or accession numbers) overlapping with Ensembl:",len(ensembl_transcript_clusters)
    print "With NO overlapp",len(no_match_list)
    return ensembl_transcript_clusters,no_match_list

def getChromosomalOveralapSpecific(ucsc_db,ensembl_db):
    print len(ucsc_db),len(ensembl_db); start_time = time.time()
    """Find transcript_clusters that have overlapping start positions with Ensembl gene start and end (based on first and last exons)"""
    ###exon_location[transcript_cluster_id,chr,strand] = [(start,stop,exon_type,probeset_id)]
    y = 0; l =0; ensembl_transcript_clusters={}; no_match_list=[]
    ###(bp1,ep1) = (47211632,47869699); (bp2,ep2)  =  (47216942, 47240877)
    for (chr,bp1,ep1) in ucsc_db:
        x = 0
        accession,ucsc_chr,ucsc_strand,ensembls = ucsc_db[(chr,bp1,ep1)]
        status = 'go'
        for ensembl in ensembls:
                y += 1; bp2,ep2 = ensembl_db[ensembl]
                add = 0
                if ((bp1 >= bp2)  and (ep2 >= bp1)) or ((ep1 >= bp2)  and (ep2 >= ep1)): add = 1 ###if the start or stop of the UCSC region is inside the Ensembl start and stop
                elif ((bp2 >= bp1)  and (ep1 >= bp2)) or ((ep2 >= bp1)  and (ep1 >= ep2)): add = 1 ###opposite
                if add == 1:
                    #if ((bp1 >= bp2)  and (ep2 >= ep1)):      
                    x = 1
                    try: ensembl_transcript_clusters[accession].append(ensembl)
                    except KeyError: ensembl_transcript_clusters[accession] = [ensembl]
                    l += 1; status = 'break'
        if x == 0: no_match_list.append(accession)

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "UCSC genes matched up to Ensembl in %d seconds" % time_diff            
    print "UCSC mRNA accession numbers overlapping with Ensembl:",len(ensembl_transcript_clusters)
    print "With NO overlapp",len(no_match_list)
    return ensembl_transcript_clusters,no_match_list

def customDeepCopy(db):
    db2={}
    for i in db:
        for e in db[i]:
            try: db2[i].append(e)
            except KeyError: db2[i]=[e]
    return db2

def identifyNewExonsForAnalysis(ensembl_transcript_clusters,no_match_list,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates):
    ###There's currently a lot of data here we are not using: Ensembl exons and transcript data
    ###Currently deemed as "over-kill" to use this.

    all_accession_gene_associations = {} ###record this for export: needed to align all transcripts to genes for probeset seqeunce alignment
    key = ('n3578(miRNA)_20337643', '+', '22')
    ucsc_multiple_alignments = {}; ensembl_gene_accession_structures={}
    for clusterid in ensembl_transcript_clusters:
        ens_geneids = ensembl_transcript_clusters[clusterid]
        ucsc_transcripts = ucsc_interim_gene_transcript_db[clusterid]
        
        for accession in ucsc_transcripts:
            try: all_accession_gene_associations[accession] += ens_geneids
            except KeyError: all_accession_gene_associations[accession] = ens_geneids
        if len(ens_geneids)>1: ###If a cluster ID associates with multiple Ensembl IDs
            chr = ucsc_annotations[clusterid]
            strand = clusterid[1]
            for accession in ucsc_transcripts:
                if test == 'yes': print "Multiple Ensembls for",accession
                a = ucsc_transcript_coordinates[accession]; a.sort(); trans_start = a[0];trans_stop = a[-1]
                ucsc_multiple_alignments[(chr,trans_start,trans_stop)] = accession,chr,strand,ens_geneids
        else:
            ens_geneid = ens_geneids[0]
            transcripts = ucsc_interim_gene_transcript_db[clusterid]
            for accession in transcripts:
                if test == 'yes': print "One Ensembl for",accession
                exon_structure_list = transcript_structure_db[accession]
                try: ensembl_gene_accession_structures[ens_geneid].append((accession,exon_structure_list))
                except KeyError: ensembl_gene_accession_structures[ens_geneid]= [(accession,exon_structure_list)]

    ###Create a new database for ensembl gene boundaries indexed by ensembl id for quicker reference in the faster lookup chr overlap function
    ensembl_gene_coordinates2={}
    for chr in ensembl_chr_coordinate_db:
        ensembl_gene_coordinates_data = ensembl_chr_coordinate_db[chr]
        for (bp,ep) in ensembl_gene_coordinates_data:
            ensembl,strand = ensembl_gene_coordinates_data[(bp,ep)]
            ensembl_gene_coordinates2[ensembl] = (bp,ep)

    ###Add all Accession #'s, for which the cluster ID did not correspond to a gene and re-run the chr overlap function
    accession_chr_coordinate_db={}
    for clusterid in no_match_list:
        ucsc_transcripts = ucsc_interim_gene_transcript_db[clusterid]
        chr = ucsc_annotations[clusterid]
        strand = clusterid[1]
        for accession in ucsc_transcripts:
            a = ucsc_transcript_coordinates[accession]; a.sort(); trans_start = a[0];trans_stop = a[-1]
            if chr in accession_chr_coordinate_db:
                accession_gene_coordinate_db = accession_chr_coordinate_db[chr]
                accession_gene_coordinate_db[(trans_start,trans_stop)] = accession,strand
            else:
                accession_gene_coordinate_db={}
                accession_gene_coordinate_db[(trans_start,trans_stop)] = accession,strand
                accession_chr_coordinate_db[chr] = accession_gene_coordinate_db
                
    ###Re-run this query with the accession numbers rather than the cluster number (may take a while)
    new_ensembl_transcript_clusters,new_no_match_list = getChromosomalOveralap(accession_chr_coordinate_db,ensembl_chr_coordinate_db)
    
    ###Add the single gene accession associations to ensembl_gene_accession_structures
    for accession in new_ensembl_transcript_clusters:
        ens_geneids = new_ensembl_transcript_clusters[accession]
        try: all_accession_gene_associations[accession] += ens_geneids
        except KeyError: all_accession_gene_associations[accession] = ens_geneids
        if len(ens_geneids)>1: ###If a cluster ID associates with multiple Ensembl IDs
            null = []###don't do anything with these
        else:
            ens_geneid = ens_geneids[0]
            exon_structure_list = transcript_structure_db[accession]
            try: ensembl_gene_accession_structures[ens_geneid].append((accession,exon_structure_list))
            except KeyError: ensembl_gene_accession_structures[ens_geneid]= [(accession,exon_structure_list)]
            
    ###Re-run the chromosomal overlap analysis specifically on transcripts, where the gene overlapped with multiple ensembls
    ensembl_transcript_clusters2,no_match_list2 = getChromosomalOveralapSpecific(ucsc_multiple_alignments,ensembl_gene_coordinates2)

    """At this point we have no_match_list and no_match_list2, both of which are theoretically non-coding RNAs that do not
    align to a known gene (within the region of a gene)"""
    
    for accession in ensembl_transcript_clusters2:
        ens_geneids = ensembl_transcript_clusters2[accession]
        all_accession_gene_associations[accession] = ens_geneids ###Write over existing if a more specific set of gene associations found
        if len(ens_geneids)==1: ###Otherwise there are multiple associations
            ens_geneid = ens_geneids[0]
            exon_structure_list = transcript_structure_db[accession]
            ###This is the list of Ensembl genes to GenBank accessions and exon coordiantes
            try: ensembl_gene_accession_structures[ens_geneid].append((accession,exon_structure_list))
            except KeyError: ensembl_gene_accession_structures[ens_geneid]= [(accession,exon_structure_list)]

    if analysis_type != 'ncRNA':
        ###Verify accession to gene associations for multiple associations or pick the one propper gene call among several incorrect
        """A problem is that if an Ensembl pseudo-transcript (not a real gene), with no exons overlapping with UCSC transcript exists, the accession
        could not be annotated with a gene, but this is not ideal, since the exons in the transcript may just overlap with one gene"""
        all_accession_gene_associations2 = []; number_of_associated_exons={}; removed=0; ensembl_gene_accession_structures_deleted={}; exon_annotation_del_db={}
        for accession in all_accession_gene_associations:
            exon_structure = transcript_structure_db[accession] ###coordinates for 'exons' provided by UCSC
            unique_genes = unique.unique(all_accession_gene_associations[accession])
            ensembl_gene_exons_temp={}
            for gene in unique_genes:
                chr,strand = ensembl_annotations[gene]
                for exon_coordinates in exon_structure:
                    exons = ensembl_gene_exon_db[gene]  ###Ensembl exonids for this gene
                    new_exon_coordinates = chr,exon_coordinates[0],exon_coordinates[1] ###create a exon coordiante tuple analagous to the one created for Ensembl
                    if new_exon_coordinates in ensembl_exon_data:  ###Therefore this is an Ensembl aligning exon (same start and stop)
                        ensembl_exon = ensembl_exon_data[new_exon_coordinates]  ###Get the id for this exon
                        if ensembl_exon in exons:  ###Since this exon could be in any gene, check to make sure it's specific for just this one
                            try: ensembl_gene_exons_temp[gene].append(ensembl_exon)
                            except KeyError: ensembl_gene_exons_temp[gene] = [ensembl_exon]
            #if accession == 'X97298': print accession, unique_genes, ensembl_gene_exons_temp, len(ensembl_gene_exons_temp);kill
            #if strand == '+': print accession, unique_genes, ensembl_gene_exons_temp, len(ensembl_gene_exons_temp);kill 
            if len(ensembl_gene_exons_temp) == 1: ###Therefore, only one Ensembl gene contained overlapping exons
                for gene in ensembl_gene_exons_temp:
                    all_accession_gene_associations[accession] = [gene]
                    number_of_associated_exons[gene,accession] = len(ensembl_gene_exons_temp[gene])
                    if len(unique_genes)>1: ###If multiple genes, then this accession number has not been updated in our main accession to Ensembl Dbase
                        exon_structure_list = transcript_structure_db[accession]
                        try: ensembl_gene_accession_structures[gene].append((accession,exon_structure_list))
                        except KeyError: ensembl_gene_accession_structures[gene]= [(accession,exon_structure_list)]
            elif len(ensembl_gene_exons_temp) == 0 or len(ensembl_gene_exons_temp) > 1:
                ###Therfore, no Ensembl exon overlaps with the transcript or exons overlapp with several genes.  If in the main accession to Ensembl Dbase, delete it
                for gene in unique_genes:
                    if gene in ensembl_gene_accession_structures:
                        accession_data_list = ensembl_gene_accession_structures[gene]
                        index = 0
                        for accession_data in accession_data_list:
                            if accession in accession_data:
                                del accession_data_list[index]; removed +=1
                                ### add all of the gene accession info to a new db to look for overlap with UCSC annotated alt events
                                if len(ensembl_gene_exons_temp) == 0 and len(unique_genes) == 1: ### This occurs if a transcript has no overlaping ensembl exons, but does contain an annotated event according to UCSC
                                    all_accession_gene_associations[accession] = [gene]  ###although the entry is deleted, probably no issue with exporting the data to LinkEST
                                    if len(unique_genes)==1: ###If multiple genes, then this accession number has not been updated in our main accession to Ensembl Dbase
                                        exon_structure_list = transcript_structure_db[accession]
                                        try: ensembl_gene_accession_structures_deleted[gene].append((accession,exon_structure_list))
                                        except KeyError: ensembl_gene_accession_structures_deleted[gene]= [(accession,exon_structure_list)]
                                        chr,strand = ensembl_annotations[gene]
                                        exon_annotation_del_db[(gene,chr,strand)] = exon_structure_list ###This should mimic the Ensembl database used for alignToKnownAlt
                            index += 1
                            
        ###Check to see if any of the unique accession-gene associations that didn't have an ensembl exon overlap with a known UCSC alt-event
        ###first update gene coordiantes with approved UCSC mRNAs (Ensembl frak's up sometmes and actually makes mRNAs too short)
        for gene in ensembl_gene_accession_structures:
            for (accession,exon_structure_list) in ensembl_gene_accession_structures[gene]:
                for exon_info in exon_structure_list:
                    exon_start = exon_info[0]; exon_stop = exon_info[1]
                    ensembl_gene_coordinates[gene].append(exon_start); ensembl_gene_coordinates[gene].append(exon_stop)
        ucsc_splicing_annot_db = alignToKnownAlt.importEnsExonStructureData(species,ensembl_gene_coordinates,ensembl_annotations,exon_annotation_del_db)
        # ucsc_splicing_annot_db[ensembl].append((start,stop,annotation_str))

        ### Only retain genes that have at least one example of UCSC alternative splicing    
        for gene in ensembl_gene_accession_structures_deleted:
            if gene in ucsc_splicing_annot_db:
                for (accession,exon_structure_list) in ensembl_gene_accession_structures_deleted[gene]:
                    add = []
                    for ucsc_exon_info in ucsc_splicing_annot_db[gene]:
                        bp1 = ucsc_exon_info[0]; ep1 = ucsc_exon_info[1]; annotation = ucsc_exon_info[2]
                        for exon_info in exon_structure_list:
                            bp2 = exon_info[0]; ep2 = exon_info[1]
                            #if accession == 'BC092441':
                            if ((bp1 >= bp2) and (ep2 >= bp1)) or ((ep1 >= bp2) and (ep2 >= ep1)): add.append(annotation) ###if the start or stop of the UCSC Alt is inside the UCSC mRNA exon start and stop
                            elif ((bp2 >= bp1) and (ep1 >= bp2)) or ((ep2 >= bp1) and (ep1 >= ep2)): add.append(annotation) ###opposite
                    if len(add)>0:
                        try: ensembl_gene_accession_structures[gene].append((accession,exon_structure_list))
                        except KeyError: ensembl_gene_accession_structures[gene]= [(accession,exon_structure_list)]

        print removed, "accessions removed from analysis"
        ###Export all possible transcript to ensembl anntotations
        """This is used for LinkESTSequenceToArrays, for figuring out which UCSC mRNAs can be specifically aligned to which Ensembl genes"""
        export_file = string.replace(input_gene_file,'all',species+'_UCSC-accession-to-gene')
        fn=filepath(export_file); data = open(fn,'w')
        for accession in all_accession_gene_associations:
            unique_genes = unique.unique(all_accession_gene_associations[accession])
            unique_genes = string.join(unique_genes,'|')
            if (unique_genes,accession) in number_of_associated_exons:
                number_of_ens_exons = number_of_associated_exons[(unique_genes,accession)]
            else: number_of_ens_exons = 'NA'
            values = accession +'\t'+ unique_genes +'\t'+ str(number_of_ens_exons)+'\n'
            data.write(values)
        data.close()
    return ensembl_gene_accession_structures

def matchUCSCExonsToEnsembl(ensembl_gene_accession_structures):
    global distal_ens_exon_pos_db; distal_ens_exon_pos_db = {}
    ###Use this database to determine which start and stop exon from Ensembl correspond to UCSC start and stop exons (with different distal positions and same junctions)
    for transcript in ensembl_transcript_structure_db:
        ens_geneid = ens_transcript_gene_db[transcript]
        chr,strand = ensembl_annotations[ens_geneid]
        """
        ###Performs the same operation as below but just for perifpheral ensembl exons, not all of them (less conservative)
        a = ensembl_transcript_structure_db[transcript]; a.sort()
        if strand == '+':
            start_pos = a[0]; start_ensembl_exon = ensembl_exon_data[a[0]]
            end_pos = a[-1]; stop_ensembl_exon = ensembl_exon_data[a[-1]]
        else:
            start_pos = a[-1]; start_ensembl_exon = ensembl_exon_data[a[-1]]
            end_pos = a[0]; stop_ensembl_exon = ensembl_exon_data[a[0]]
        v = [(start_pos,start_ensembl_exon),(end_pos,stop_ensembl_exon)]
        try: distal_ens_exon_pos_db[ens_geneid] += v
        except KeyError: distal_ens_exon_pos_db[ens_geneid] = v"""
        trans_data = ensembl_transcript_structure_db[transcript]; trans_data.sort()
        for a in trans_data:
            ###Do this for all exons, since we don't care if an mRNA from UCSC starts in the middl of the the 3rd exon
            #ensembl_exon_data[(chr,exon_start,exon_end)] = ens_exonid
            pos = a; ensembl_exon = ensembl_exon_data[a]
            v = [(pos,ensembl_exon)]
            try: distal_ens_exon_pos_db[ens_geneid] += v
            except KeyError: distal_ens_exon_pos_db[ens_geneid] = v

    distal_ens_exon_pos_db = eliminate_redundant_dict_values(distal_ens_exon_pos_db)        
    ###Determine which exons are constitutive, by simply counting their number in each transcript
    constitutive_gene_db={}; coordinate_to_ens_exon = {}
    for ens_geneid in ensembl_gene_accession_structures:
        a = ensembl_gene_accession_structures[ens_geneid]
        coordinate_count_db={}
        chr,strand = ensembl_annotations[ens_geneid]
        for (accession,exon_structure) in a:
            index=1
            for exon_coordinates in exon_structure:
                ###check whether the exon is an ensembl exon
                exon_coordinates = list(exon_coordinates)
                new_exon_coordinates = chr,exon_coordinates[0],exon_coordinates[1]
                ensembl_exon = 'null'
                #if new_exon_coordinates == ('1', 112109882, 112111042): print index, accession,len(exon_structure);kill
                try:
                    ensembl_exon = ensembl_exon_data[new_exon_coordinates]
                    #if ensembl_exon == 'ENSE00001473402': print new_exon_coordinates;kill
                    coordinate_to_ens_exon[new_exon_coordinates] = ensembl_exon
                except KeyError:
                  if len(exon_structure)>1: ###Ensures we don't do this for single exon transcripts, where this would be in-appropriate
                    if index == 1 and strand == '+': #or index == len(exon_structure)
                        for (pos,exon) in distal_ens_exon_pos_db[ens_geneid]:
                            cr,start,stop = pos
                            if exon_coordinates[1] == stop: ensembl_exon = exon; exon_coordinates[0] = start #; print 'a'
                    elif index == 1 and strand == '-': #or index == len(exon_structure)
                        for (pos,exon) in distal_ens_exon_pos_db[ens_geneid]:
                            cr,start,stop = pos
                            if exon_coordinates[0] == start: ensembl_exon = exon; exon_coordinates[1] = stop #; print 'b'
                    elif index == len(exon_structure) and strand == '+': #or index == len(exon_structure)
                        for (pos,exon) in distal_ens_exon_pos_db[ens_geneid]:
                            cr,start,stop = pos
                            if exon_coordinates[0] == start: ensembl_exon = exon; exon_coordinates[1] = stop #; print 'c'
                    elif index == len(exon_structure) and strand == '-': #or index == len(exon_structure)
                        for (pos,exon) in distal_ens_exon_pos_db[ens_geneid]:
                            cr,start,stop = pos
                            if exon_coordinates[1] == stop: ensembl_exon = exon; exon_coordinates[0] = start #; print 'd',exon_coordinates,start,stop,ensembl_exon;kill
                    if ensembl_exon != 'null':
                        #print accession,ensembl_exon,exon_coordinates;kill
                        coordinate_to_ens_exon[new_exon_coordinates] = ensembl_exon
                index+=1
                ###count each exon in all transcripts
                exon_coordinates = tuple(exon_coordinates)
                try: coordinate_count_db[exon_coordinates]+=1
                except KeyError: coordinate_count_db[exon_coordinates]=1
        count_list=[]
        ###process for identifying putative constitutive IDs.  Found cases were it eliminate real constitutive exons (check out SLK and EF139853...other transcripts eliminated due to ncRNA overlap)
        """
        print coordinate_count_db,'\n'
        ###Now do this for Ensembl gene exons
        for exon_coordinates in coordinate_gene_count_db[ens_geneid]:
            try: coordinate_count_db[tuple(exon_coordinates)]+=1
            except KeyError: coordinate_count_db[tuple(exon_coordinates)]=1
        print coordinate_count_db;kill
        for exon_coordinates in coordinate_count_db:
            count = coordinate_count_db[exon_coordinates]
            count_list.append((count,exon_coordinates))
        count_list.sort(); count_list.reverse(); max_count = count_list[0][0]; constitutive=[]
        for (count,coor) in count_list:
            if count == max_count: constitutive.append(coor)
        constitutive_gene_db[ens_geneid] = constitutive"""
    
    return ensembl_gene_accession_structures,constitutive_gene_db,coordinate_to_ens_exon

def exportNonGeneExonStructures(ensembl_gene_accession_structures,transcript_structure_db):
    export_file = string.replace(input_gene_file,'all',species+'_UCSC_NonGene_structure'); ensembl_linked_accession_db={}
    fn=filepath(export_file); data = open(fn,'w')
    title = ['Ensembl Gene ID','Chromosome','Strand','Exon Start (bp)','Exon End (bp)','Custom Exon ID','Constitutive Exon','NCBI Accession']
    title = string.join(title,'\t')+'\n'
    data.write(title)
    for ens_geneid in ensembl_gene_accession_structures:
        for (accession,exon_structure) in ensembl_gene_accession_structures[ens_geneid]:
            ensembl_linked_accession_db[accession] = []
    non_gene_acs = 0
    for accession in transcript_structure_db:
        if accession not in ensembl_linked_accession_db:
            non_gene_acs+=1
            chr,strand = ucsc_genomic_accession_db[accession]
            common_values = [accession,chr,strand]
            exon_structure = transcript_structure_db[accession]             
            index=1
            ens_exon_count=[]###See if we should export this transcript (if it contains no unique exons)
            for (exon_start,exon_stop) in exon_structure:
                exonid = accession+'-'+str(index) ###custom ID designating the relative exon position in the exon_structure
                constitutive_call = '0'
                values = common_values+[str(exon_start),str(exon_stop),exonid,constitutive_call,accession]
                index+=1
                values = string.join(values,'\t')+'\n'
                data.write(values)
    data.close()
    print non_gene_acs, "NON-Ensembl linked accession numbers exported"
    
def exportExonClusters(ensembl_gene_accession_structures,constitutive_gene_db,coordinate_to_ens_exon,species):
    export_file = string.replace(input_gene_file,'all',species+'_UCSC_transcript_structure'); accessions_exlcuded=0; accessions_included=0
    if analysis_type == 'ncRNA': export_file = string.replace(export_file,'structure','structure_'+str(round))
    fn=filepath(export_file); data = open(fn,'w')
    title = ['Ensembl Gene ID','Chromosome','Strand','Exon Start (bp)','Exon End (bp)','Custom Exon ID','Constitutive Exon','NCBI Accession']
    title = string.join(title,'\t')+'\n'
    data.write(title)
    for ens_geneid in ensembl_gene_accession_structures:
        chr,strand = ensembl_annotations[ens_geneid]
        common_values = [ens_geneid,chr,strand]
        for (accession,exon_structure) in ensembl_gene_accession_structures[ens_geneid]:
            index=1
            #constitutive_coordinates = constitutive_gene_db[ens_geneid]
            ens_exon_count=[]###See if we should export this transcript (if it contains no unique exons)
            for (exon_start,exon_stop) in exon_structure:
                try:
                    exonid = coordinate_to_ens_exon[(chr,exon_start,exon_stop)]
                    ###Verify that the exon corresponds to that gene (some Ensembl exon regions belong to more than one gene)
                    if exonid in ensembl_gene_exon_db[ens_geneid]: ens_exon_count.append(exonid)
                except KeyError: null = []
            if (len(ens_exon_count) != len(exon_structure)) or (analysis_type == 'ncRNA'):  ###REMOVE TRANSCRIPTS FOR WHICH THERE ARE ONLY ENSEMBL EXONS: This is an issue since we really want to get rid of transcripts with only Ensembl Junctions (too strict)
                accessions_included+=1
                for (exon_start,exon_stop) in exon_structure:
                    ###used to try to emperically determine which exons are constitutive... instead, just trust Ensembl
                    #if (exon_start,exon_stop) in constitutive_coordinates: constitutive_call = '1'
                    #else: constitutive_call = '0'
                    try:
                        exonid = coordinate_to_ens_exon[(chr,exon_start,exon_stop)]
                        if exonid in ensembl_gene_exon_db[ens_geneid]: exonid = exonid ###Perform the same check as above
                        else: exonid = accession+'-'+str(index) 
                    except KeyError: exonid = accession+'-'+str(index) ###custom ID designating the relative exon position in the exon_structure
                    if exonid in ensembl_const_exon_db: constitutive_call = ensembl_const_exon_db[exonid]
                    else: constitutive_call = '0'
                    if analysis_type == 'ncRNA':
                        ###For ncRNA analyses, the ncRNA coordinates are small and many can overlap with a single AC, so use the primary AC coordinates
                        gene_coordinates = coordinate_gene_count_db[ens_geneid]
                        exon_start = gene_coordinates[0][0]; exon_start = gene_coordinates[-1][1]
                    values = common_values+[str(exon_start),str(exon_stop),exonid,constitutive_call,accession]
                    index+=1
                    values = string.join(values,'\t')+'\n'
                    data.write(values)
            else: accessions_exlcuded+=1
    data.close()
    print accessions_exlcuded, "Accession numbers excluded (same as Ensembl transcript), out of",(accessions_included+accessions_exlcuded)
    print 'data written to:',export_file

def exportSimple(filtered_data,title,input_gene_file):
    export_file = string.replace(input_gene_file,'all',species+'_UCSC_transcript_structure_filtered')
    fn=filepath(export_file); data = open(fn,'w')
    data.write(title)
    for index,line in filtered_data: data.write(line)
    data.close()
    print 'data written to:',export_file
    
def filterBuiltAssociations(species):
    if analysis_type == 'mRNA': input_gene_file = 'AltDatabase/ucsc/'+species+'/'+mRNA_type+'.txt'
    elif analysis_type == 'ncRNA': input_gene_file = 'AltDatabase/ucsc/'+species+'/all_ncRNA.txt'
    filename = string.replace(input_gene_file,'all',species+'_UCSC_transcript_structure')

    ###Re-import the data and filter it to remove junctions that should be present in the Ensembl database (ensembl-ensembl junction)
    fn=filepath(filename); x=0
    raw_data = {}; global accession_db; accession_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: title = line; x+=1 ###first line
        else:
            gene, chr, strand, exon_start, exon_end, exonid, constitutive_exon, accession = t
            try: raw_data[accession].append([exonid,x,line])
            except KeyError: raw_data[accession] = [[exonid,x,line]]
            #y = EnsemblImport.ExonStructureData(gene,chr,strand,exon_start,exon_end,constitutive_exon,exonid,accession)
            x+=1
            
    keep=[]; k=0
    ###Remove the terminal portions of the transcripts that are comprised ONLY of continuous Ensembl exons
    for accession in raw_data:
        #accession = 'AK027479'
        index=0; x = raw_data[accession]
        while index<(len(x)-1): ###for each exon in the accession
            try:
                y,n,n = x[index]
                s,n,n = x[index+1] ###Next exon in transcript
                #elif 'ENS' in y and 'ENS' in s: k+=1 #; print (y,s),accession;kill
                if (y,s) in ensembl_exon_pairs or (s,y) in ensembl_exon_pairs: x = x[index+1:]; index = -1
                else: break
            except IndexError: break
            index+=1
        index=0; x.reverse()
        while index<(len(x)-1): ###for each exon in the accession
            try:
                y,n,n = x[index]
                s,n,n = x[index+1] ###Next exon in transcript
                #elif 'ENS' in y and 'ENS' in s: k+=1#; print (y,s),accession;kill
                if (y,s) in ensembl_exon_pairs or (s,y) in ensembl_exon_pairs: x = x[index+1:]; index = -1
                else: x.reverse(); raw_data[accession] = x; break
            except IndexError: break
            index+=1
        for (exonid,i,line) in raw_data[accession]: keep.append((i,line))
    keep = unique.unique(keep); keep.sort()
    #print k, "exon pairs that contained two ensembl exons not paired in an Ensembl transcript";kill
    exportSimple(keep,title,input_gene_file)


def importExonStructureDataGeneric(filename,filter_db,raw_data):
    ###Re-import the data and filter it to remove junctions that should be present in the Ensembl database (ensembl-ensembl junction)
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: title = line; x+=1 ###first line
        else:
            gene, chr, strand, exon_start, exon_end, exonid, constitutive_exon, accession = t
            proceed = 'no'
            if len(filter_db)>0:
                if gene in filter_db: proceed = 'yes' ### for Ensembl, only include predicted ncRNAs
            elif accession not in raw_data: proceed = 'yes' ### only include ncRNAs not already added based on accession overlapp
            if proceed == 'yes':
                try: raw_data[accession].append([exonid,x,line])
                except KeyError: raw_data[accession] = [[exonid,x,line]]
            x+=1
    return raw_data

def combineNCRNAPredictions(species):

    ### Include ncRNAs NOT inlcuded in the lists from http://www.noncode.org/
    ### These include missing microRNAs, tRNAs, rRNAs, snoRNAs, etc.
    ensembl_annotation_db = getEnsemblAnnotations(species)
    additional_ncRNAs={}
    for symbol in ensembl_annotation_db:
        description,ens_geneid = ensembl_annotation_db[symbol]
        if 'hsa-mir' in description: additional_ncRNAs[ens_geneid]=[]
        #elif ('rna' in description or 'RNA' in description) and 'mRNA' not in description: additional_ncRNAs[ens_geneid]=[]
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    raw_data = importExonStructureDataGeneric(filename,additional_ncRNAs,{})
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_1_ncRNA.txt'
    raw_data = importExonStructureDataGeneric(filename,{},raw_data)
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_2_ncRNA.txt'
    raw_data = importExonStructureDataGeneric(filename,{},raw_data)
    print len(raw_data);kill
    #exportSimple(keep,title,input_gene_file)

def getEnsemblAnnotations(species):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    fn=filepath(filename); ensembl_annotation_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if 'Description' in data: ensembl_gene_id,description,symbol = string.split(data,'\t')
        else: ensembl_gene_id,symbol,description = string.split(data,'\t')
        ensembl_annotation_db[symbol] = description, ensembl_gene_id
    return ensembl_annotation_db

def getUCSCAssociations(Species):
    global species; species = Species
    global ensembl_gene_coordinates
    ensembl_chr_coordinate_db = importEnsExonStructureData(species)
    ucsc_gene_coordinates,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates = importUCSCExonStructureData(species)
    ensembl_transcript_clusters,no_match_list = getChromosomalOveralap(ucsc_gene_coordinates,ensembl_chr_coordinate_db)
    ensembl_gene_accession_structures,constitutive_gene_db = identifyNewExonsForAnalysis(ensembl_transcript_clusters,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates)
    exportExonClusters(ensembl_gene_accession_structures,constitutive_gene_db,species)

def runUCSCEnsemblAssociations(Species,Analysis_type,Round,mRNA_Type,run_from_scratch):
    global species; species = Species; global round; round = Round; global mRNA_type; mRNA_type = mRNA_Type
    global test; global bp_offset; global gap_length; global test_gene; global analysis_type; analysis_type = Analysis_type
    bp_offset = 100 ###allowed excesive base pairs added to the distal ends of the Ensembl genes
    gap_length = 12 ###maximum allowed gap length to qualify for merging exons
    test = 'no'
    test_gene = ['ENSMUSG00000022194']#,'ENSG00000154889','ENSG00000156026','ENSG00000148584','ENSG00000063176','ENSG00000126860'] #['ENSG00000105968']
    
    if run_from_scratch == 'yes':
        global ensembl_chr_coordinate_db
        ensembl_chr_coordinate_db = importEnsExonStructureData(species)
        ucsc_chr_coordinate_db,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates = importUCSCExonStructureData(species)
        ensembl_transcript_clusters,no_match_list = getChromosomalOveralap(ucsc_chr_coordinate_db,ensembl_chr_coordinate_db)
        ensembl_gene_accession_structures = identifyNewExonsForAnalysis(ensembl_transcript_clusters,no_match_list,transcript_structure_db,ucsc_interim_gene_transcript_db,ucsc_transcript_coordinates)
        
        ensembl_gene_accession_structures,constitutive_gene_db,coordinate_to_ens_exon = matchUCSCExonsToEnsembl(ensembl_gene_accession_structures)
        exportExonClusters(ensembl_gene_accession_structures,constitutive_gene_db,coordinate_to_ens_exon,species)
        exportNonGeneExonStructures(ensembl_gene_accession_structures,transcript_structure_db)
        if (analysis_type == 'mRNA'): filterBuiltAssociations(species)
    else:
        if (analysis_type == 'mRNA'):
            ensembl_chr_coordinate_db = importEnsExonStructureData(species) ###need this to get the unique Ensembl pairs
            filterBuiltAssociations(species)

def updateUCSCEnsemblAssociations(Species,Analysis_type,mRNA_Type,run_from_scratch):
    #Round = 1; runUCSCEnsemblAssociations(Species,Analysis_type,Round,mRNA_Type,run_from_scratch)
    #if Analysis_type != 'mRNA':
    Round = 2; runUCSCEnsemblAssociations(Species,Analysis_type,Round,mRNA_Type,run_from_scratch)
        
if __name__ == '__main__':
    dirfile = unique
    run_from_scratch = 'yes'
    Species = 'Hs'
    mRNA_Type = 'est'
    #Analysis_type = 'mRNA'
    #updateUCSCEnsemblAssociations(Species,Analysis_type,mRNA_Type,run_from_scratch)
    Analysis_type = 'ncRNA'
    #combineNCRNAPredictions(Species);kill
    updateUCSCEnsemblAssociations(Species,Analysis_type,mRNA_Type,run_from_scratch)