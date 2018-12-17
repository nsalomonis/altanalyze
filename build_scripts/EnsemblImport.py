###EnsemblImport
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
from build_scripts import GO_parsing
import copy
import time
try: from build_scripts import alignToKnownAlt
except Exception: pass ### circular import error
import export
import traceback

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def findAligningIntronBlock(ed,intron_block_db,exon_block_db,last_region_db):
    aligned_region = 'no'
    for block in intron_block_db:
        for rd in intron_block_db[block]:
            intron_pos = [rd.ExonStart(), rd.ExonStop()]; intron_pos.sort()
            intronid = rd.IntronRegionID() ### The intron block (relative to the exon block) is assigned in the function exon_clustering
            retained_pos = [ed.ExonStart(), ed.ExonStop()]; retained_pos.sort()
            aligned_region = 'no'
            #print ed.GeneID(), rd.GeneID(), retained_pos,intron_pos,intronid
            if intron_pos == retained_pos: aligned_region = 'yes'
            elif intron_pos[0] == retained_pos[0] or intron_pos[1] == retained_pos[1]: aligned_region = 'yes'
            elif (intron_pos[0]<retained_pos[0] and intron_pos[1]>retained_pos[0]) or (intron_pos[0]<retained_pos[1] and intron_pos[1]>retained_pos[1]): aligned_region = 'yes'
            elif (retained_pos[0]<intron_pos[0] and retained_pos[1]>intron_pos[0]) or (retained_pos[0]<intron_pos[1] and retained_pos[1]>intron_pos[1]): aligned_region = 'yes'
            if aligned_region == 'yes':
                #print 'intron-aligned',intronid
                ed.setAssociatedSplicingEvent('intron-retention')
                rd.updateDistalIntronRegion()
                intronid = intronid[:-1]+str(rd.DistalIntronRegion())
                intronid = string.replace(intronid,'-','.')
                ed.setNewIntronRegion(intronid); break
        if aligned_region == 'yes': break
        
    if aligned_region == 'no':                
        for block in exon_block_db:
            for rd in exon_block_db[block]:
                exonregion_pos = [rd.ExonStart(), rd.ExonStop()]; exonregion_pos.sort()
                exonid = rd.ExonRegionID() ### The exon block (relative to the exon block) is assigned in the function exon_clustering
                retained_pos = [ed.ExonStart(), ed.ExonStop()]; retained_pos.sort()
                aligned_region = 'no'
                #print ed.GeneID(),rd.GeneID(), retained_pos,exonregion_pos,exonid
                if exonregion_pos == retained_pos: aligned_region = 'yes'
                elif exonregion_pos[0] == retained_pos[0] or exonregion_pos[1] == retained_pos[1]: aligned_region = 'yes'
                elif (exonregion_pos[0]<retained_pos[0] and exonregion_pos[1]>retained_pos[0]) or (exonregion_pos[0]<retained_pos[1] and exonregion_pos[1]>retained_pos[1]): aligned_region = 'yes'
                elif (retained_pos[0]<exonregion_pos[0] and retained_pos[1]>exonregion_pos[0]) or (retained_pos[0]<exonregion_pos[1] and retained_pos[1]>exonregion_pos[1]): aligned_region = 'yes'
                if aligned_region == 'yes':
                    #print 'exon aligned',exonid
                    #print len(last_region_db),last_region_db[rd.ExonNumber()]
                    last_region_db[rd.ExonNumber()].sort(); last_region = last_region_db[rd.ExonNumber()][-1]
                    last_region_db[rd.ExonNumber()].append(last_region+1)
                    #print len(last_region_db),last_region_db[rd.ExonNumber()],'E'+str(rd.ExonNumber())+'.'+str(last_region+1)
                    ed.setAssociatedSplicingEvent('exon-region-exclusion')
                    exonid = 'E'+str(rd.ExonNumber())+'.'+str(last_region+1)
                    ed.setNewIntronRegion(exonid); break
            if aligned_region == 'yes': break
    return ed,last_region_db

def getUCSCSplicingAnnotations(ucsc_events,splice_events,start,stop):
    splice_events = string.split(splice_events,'|')    
    for (r_start,r_stop,splice_event) in ucsc_events:
        if ((start >= r_start) and (start < r_stop)) or ((stop > r_start) and (stop <= r_stop)):
            splice_events.append(splice_event)
        elif ((r_start >= start) and (r_start <= stop)) or ((r_stop >= start) and (r_stop <= stop)): ### Applicable to polyA annotations
            splice_events.append(splice_event)
    unique.unique(splice_events)
    splice_events = string.join(splice_events,'|')
    try:
        if splice_events[0] == '|': splice_events=splice_events[1:]
    except IndexError: null=[]
    return splice_events 

def exportSubGeneViewerData(exon_regions,exon_annotation_db2,critical_gene_junction_db,intron_region_db,intron_retention_db,full_junction_db,excluded_intronic_junctions,ucsc_splicing_annot_db):
    intron_retention_db2={}
    for (gene,chr,strand) in intron_retention_db:
        for intron_info in intron_retention_db[(gene,chr,strand)]:
            pos1,pos2,ed = intron_info; pos_list=[pos1,pos2]; pos_list.sort()
            try: intron_retention_db2[gene].append(pos_list)
            except KeyError: intron_retention_db2[gene] = [pos_list]
            
    exon_annotation_export = 'AltDatabase/ensembl/'+species+'/'+species+'_SubGeneViewer_exon-structure-data.txt'
    print 'Writing the file:',exon_annotation_export
    fn=filepath(exon_annotation_export); sgvdata = open(fn,'w')
    title = ['gene','exon-id','type','block','region','constitutive','start-exon','annotation']
    title = string.join(title,'\t')+'\n'; sgvdata.write(title)

    full_exon_structure_export = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    full_junction_structure_export = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'
    exondata = export.ExportFile(full_exon_structure_export); junctiondata = export.ExportFile(full_junction_structure_export)
    exontitle = ['gene', 'exon-id', 'chromosome', 'strand', 'exon-region-start(s)', 'exon-region-stop(s)', 'constitutive_call', 'ens_exon_ids', 'splice_events', 'splice_junctions']
    exontitle = string.join(exontitle,'\t')+'\n'; exondata.write(exontitle); junctiondata.write(exontitle)

    export_annotation = '_intronic'
    alt_junction_export = 'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions'+export_annotation+'.txt'
    print 'Writing the file:',alt_junction_export
    fn=filepath(alt_junction_export); reciprocol_junction_data = open(fn,'w')
    title = ['gene','critical-exon-id','junction1','junction2']
    title = string.join(title,'\t')+'\n'; reciprocol_junction_data.write(title)
    
    ### exon_annotation_db2 contains coordinates and ensembl exon IDs for each exon - extract out just this data
    exon_coordinate_db={}
    for key in exon_annotation_db2:
        gene=key[0]; exon_coords={}
        for exon_data in exon_annotation_db2[key]:
            exon_start = exon_data[1][0]; exon_stop = exon_data[1][1]; ed = exon_data[1][2]
            exon_coords[exon_start,exon_stop] = ed.ExonID()
        exon_coordinate_db[gene] = exon_coords

    intron_coordinate_db={}; exon_region_annotations={}; intron_junction_db={}
    for key in intron_retention_db:
        gene=key[0]; exon_coords={}
        for exon_data in intron_retention_db[key]:
            exon_start = exon_data[0]; exon_stop = exon_data[1]; ed = exon_data[2]
            exon_coords[exon_start,exon_stop] = ed.ExonID()
        intron_coordinate_db[gene] = exon_coords

    for gene in exon_regions:
        previous_exonid=''; previous_intronid=''
        block_db = exon_regions[gene]
        try:intron_block_db = intron_region_db[gene]; introns = 'yes'
        except KeyError: introns = 'no'
        if gene in ucsc_splicing_annot_db: ucsc_events = ucsc_splicing_annot_db[gene]
        else: ucsc_events = []
        utr_data = [gene,'U0.1','u','0','1','n','n','']
        values = string.join(utr_data,'\t')+'\n'; sgvdata.write(values)
        index=1
        for block in block_db:
            for rd in block_db[block]:
                splice_event = rd.AssociatedSplicingEvent();exon_pos = [rd.ExonStart(), rd.ExonStop()]; exon_pos.sort()
                if gene in intron_retention_db2:
                    for retained_pos in intron_retention_db2[gene]:
                        if exon_pos == retained_pos: splice_event = 'exon-region-exclusion'
                        elif exon_pos[0] == retained_pos[0] or exon_pos[1] == retained_pos[1]: splice_event = 'exon-region-exclusion'
                        elif exon_pos[0]>retained_pos[0] and exon_pos[0]<retained_pos[1] and exon_pos[1]>retained_pos[0] and exon_pos[1]<retained_pos[1]: splice_event = 'exon-region-exclusion'
                        if 'exon-region-exclusion' in splice_event: rd.setAssociatedSplicingEvent(splice_event); id = rd
                if len(splice_event)>0: constitutive_call = 'no' ### If a splice-event is associated with a recommended constitutive region, over-ride it
                else: constitutive_call = rd.ConstitutiveCallAbrev()
                values = [gene,rd.ExonRegionID2(),'e',str(index),str(rd.RegionNumber()),constitutive_call[0],'n',splice_event]#,str(exon_pos[0]),str(exon_pos[1])]
                values = string.join(values,'\t')+'\n'; sgvdata.write(values)

                ens_exons = getMatchingEnsExons(rd.ExonStart(),rd.ExonStop(),exon_coordinate_db[gene])
                start_stop = [rd.ExonStart(),rd.ExonStop()]; start_stop.sort(); start,stop = start_stop
                splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,start,stop)
                if len(splice_event)>0: constitutive_call = 'no' ### If a splice-event is associated with a recommended constitutive region, over-ride it
                else: constitutive_call = rd.ConstitutiveCallAbrev()
                values = [gene,rd.ExonRegionID2(),'chr'+rd.Chr(),rd.Strand(),str(start),str(stop),constitutive_call,ens_exons,splice_event,rd.AssociatedSplicingJunctions()]
                values = string.join(values,'\t')+'\n'; exondata.write(values)
                exon_region_annotations[gene,rd.ExonRegionID2()]=ens_exons
                if len(previous_intronid)>0:
                    ### Intron junctions are not stored and are thus not analyzed as recipricol junctions (thus add them hear for RNASeq)
                    try: intron_junction_db[gene].append((previous_exonid,previous_intronid))
                    except Exception: intron_junction_db[gene] = [(previous_exonid,previous_intronid)]
                    intron_junction_db[gene].append((previous_intronid,rd.ExonRegionID2()))
                    values = [gene,previous_intronid,previous_exonid+'-'+previous_intronid,previous_exonid+'-'+rd.ExonRegionID2(),id.AssociatedSplicingEvent()]
                    values = string.join(values,'\t')+'\n'; reciprocol_junction_data.write(values)
                    values = [gene,previous_intronid,previous_intronid+'-'+rd.ExonRegionID2(),previous_exonid+'-'+rd.ExonRegionID2(),id.AssociatedSplicingEvent()]
                    values = string.join(values,'\t')+'\n'; reciprocol_junction_data.write(values)
                    id.setAssociatedSplicingJunctions(previous_exonid+'-'+previous_intronid+'|'+previous_intronid+'-'+rd.ExonRegionID2())
                    previous_intronid=''
                if splice_event == 'exon-region-exclusion':
                    previous_intronid = rd.ExonRegionID2()
                else: previous_exonid=rd.ExonRegionID2()
                    
            index+=1
            if introns == 'yes':
                try:
                    intronid = rd.IntronRegionID() ### The intron block (relative to the exon block) is assigned in the function exon_clustering
                    for rd in intron_block_db[block]:
                        intron_pos = [rd.ExonStart(), rd.ExonStop()]; intron_pos.sort()
                        splice_event = rd.AssociatedSplicingEvent()
                        if gene in intron_retention_db2:
                            for retained_pos in intron_retention_db2[gene]:
                                #if '15' in intronid: print intron_pos,retained_pos;kill
                                if intron_pos == retained_pos: splice_event = 'intron-retention'
                                elif intron_pos[0] == retained_pos[0] or intron_pos[1] == retained_pos[1]: splice_event = 'intron-retention'
                                elif intron_pos[0]>retained_pos[0] and intron_pos[0]<retained_pos[1] and intron_pos[1]>retained_pos[0] and intron_pos[1]<retained_pos[1]: splice_event = 'intron-retention'
                                if 'intron' in splice_event: rd.setAssociatedSplicingEvent(splice_event); id = rd
                        if len(splice_event)>0: constitutive_call = 'no' ### If a splice-event is associated with a recommended constitutive region, over-ride it
                        else: constitutive_call = rd.ConstitutiveCallAbrev()
                        values = [gene,intronid,'i',str(index),str(rd.RegionNumber()),constitutive_call[0],'n',splice_event]#,str(intron_pos[0]),str(intron_pos[1])]
                        values = string.join(values,'\t')+'\n'; sgvdata.write(values); ens_exons=''
                        
                        if len(splice_event)>0:
                            ens_exons = getMatchingEnsExons(rd.ExonStart(),rd.ExonStop(),intron_coordinate_db[gene])
                            exon_region_annotations[gene,intronid]=ens_exons
                            previous_intronid=intronid
                        start_stop = [rd.ExonStart(),rd.ExonStop()]; start_stop.sort(); start,stop = start_stop
                        splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,start,stop)
                        if len(splice_event)>0: constitutive_call = 'no' ### If a splice-event is associated with a recommended constitutive region, over-ride it
                        else: constitutive_call = rd.ConstitutiveCallAbrev()
                        values = [gene,intronid,'chr'+rd.Chr(),rd.Strand(),str(start),str(stop),constitutive_call,ens_exons,splice_event,rd.AssociatedSplicingJunctions()]
                        values = string.join(values,'\t')+'\n'; exondata.write(values)
                    index+=1
                except KeyError: null=[]
        last_exon_region,null = string.split(rd.ExonRegionID2(),'.') ### e.g. E13.1	becomes E13, 1
        
        utr_data = [gene,'U'+last_exon_region[1:]+'.1','u',str(index),'1','n','n','']
        values = string.join(utr_data,'\t')+'\n'; sgvdata.write(values)
    sgvdata.close(); reciprocol_junction_data.close()

    for gene in excluded_intronic_junctions:
        excluded_junction_exons=[]
        last_region_db={}
        exon_block_db=exon_regions[gene]
        for block in exon_block_db:
            for rd in exon_block_db[block]:
                try: last_region_db[rd.ExonNumber()].append(rd.RegionNumber())
                except KeyError: last_region_db[rd.ExonNumber()] = [rd.RegionNumber()]
        
        ###excluded_intronic_junctions These are Ensembl exons reclassified as belonging to a retained intron
        ### Nonetheless, we need to store the original junctions since they are not novel
        for (ed1,ed2) in excluded_intronic_junctions[gene]:
            #print ed1.ExonStop(),ed2.ExonStart()
            ### Store all exons aligning to retained introns and sort by position (needed for ordering)
            if ed1.IntronDeletionStatus()== 'yes':
                if ed1.Strand() == '-': start = ed1.ExonStop(); stop = ed1.ExonStart()
                else: start = ed1.ExonStart(); stop = ed1.ExonStop()
                excluded_junction_exons.append((start,stop,ed1))
            if ed2.IntronDeletionStatus()== 'yes':
                if ed2.Strand() == '-': start = ed2.ExonStop(); stop = ed2.ExonStart()
                else: start = ed2.ExonStart(); stop = ed2.ExonStop()
                excluded_junction_exons.append((start,stop,ed2))
        excluded_junction_exons = unique.unique(excluded_junction_exons); excluded_junction_exons.sort()
        for (start,stop,ed) in excluded_junction_exons:
            ### update the IntronIDs and annotations for each exon (could be present in multiple junctions)
            #print gene, ed.GeneID(), len(intron_region_db[gene]), len(exon_regions[gene]), len(last_region_db)
            try: intron_regions = intron_region_db[gene]
            except KeyError: intron_regions = [] ### In very rare cases where we have merged exons with <500bp introns, no introns for the gene will 'exist'
            ed,last_region_db=findAligningIntronBlock(ed,intron_regions,exon_regions[gene],last_region_db)
            #print [[ed.ExonID(), ed.NewIntronRegion(), start, stop]]

            ens_exons = getMatchingEnsExons(ed.ExonStart(),ed.ExonStop(),exon_coordinate_db[gene])
            start_stop = [ed.ExonStart(),ed.ExonStop()]; start_stop.sort(); start,stop = start_stop
            splice_event = getUCSCSplicingAnnotations(ucsc_events,ed.AssociatedSplicingEvent(),start,stop)
            try: values = [gene,ed.NewIntronRegion(),'chr'+ed.Chr(),ed.Strand(),str(start),str(stop),'no',ens_exons,splice_event,ed.AssociatedSplicingJunctions()]
            except Exception:
                print gene, 'chr'+ed.Chr(),ed.Strand(),str(start),str(stop),'no',ens_exons,splice_event,ed.AssociatedSplicingJunctions()
                print len(exon_regions[gene]), len(intron_regions), len(exon_regions[gene]), len(last_region_db); kill
            values = string.join(values,'\t')+'\n'; exondata.write(values)
                
    ###Export the two individual exon regions for each exon junction
    critical_gene_junction_db = eliminate_redundant_dict_values(critical_gene_junction_db)
    exon_annotation_export = 'AltDatabase/ensembl/' +species+'/'+species+ '_SubGeneViewer_junction-data.txt'
    fn=filepath(exon_annotation_export); sgvjdata = open(fn,'w')
    title = ['gene',"5'exon-region","3'exon-region"]
    title = string.join(title,'\t')+'\n'; sgvjdata.write(title)
    alt_junction={}
    for gene in critical_gene_junction_db:
        for junction_ls in critical_gene_junction_db[gene]:
            values = [gene,junction_ls[0],junction_ls[1]]
            values = string.join(values,'\t')+'\n'; sgvjdata.write(values)
            try: alt_junction[gene].append((junction_ls[0],junction_ls[1]))
            except KeyError: alt_junction[gene]=[(junction_ls[0],junction_ls[1])]
        ### Include junctions for intron retention and exon-exclusion
        if gene in intron_junction_db:
            for junction_ls in intron_junction_db[gene]:
                values = [gene,junction_ls[0],junction_ls[1]]
                values = string.join(values,'\t')+'\n'; sgvjdata.write(values)
                try: full_junction_db[gene].append((junction_ls[0],junction_ls[1]))
                except KeyError: full_junction_db[gene] = [(junction_ls[0],junction_ls[1])]
                try: alt_junction[gene].append((junction_ls[0],junction_ls[1]))
                except KeyError: alt_junction[gene]=[(junction_ls[0],junction_ls[1])]
    for gene in full_junction_db:
        ### Add junctions to the exon database
        block_db = exon_regions[gene]
        try: intron_block_db = intron_region_db[gene]
        except KeyError: intron_block_db={}  
        for (exon1,exon2) in full_junction_db[gene]:
            found1='no'; found2='no'
            for block in block_db:
                if found1 == 'no' or found2 == 'no':
                    for rd in block_db[block]:
                        if rd.ExonRegionID2() == exon1: le = rd; found1 = 'yes' ### Left exon found
                        if rd.ExonRegionID2() == exon2: re = rd; found2 = 'yes' ### Right exon found

            if found1 == 'no' or found2 == 'no':
                for block in intron_block_db:
                    for rd in intron_block_db[block]:
                        if rd.IntronRegionID() == exon1: le = rd; found1 = 'yes' ### Left exon found
                        if rd.IntronRegionID() == exon2: re = rd; found2 = 'yes' ### Right exon found                
            if found1 == 'yes' and found2 == 'yes':
                ens_exons1 = exon_region_annotations[gene,exon1]
                ens_exons2 = exon_region_annotations[gene,exon2]
                const_call = 'no'
                if gene in alt_junction:
                    if (exon1,exon2) in alt_junction[gene]: const_call = 'no'
                    elif le.ConstitutiveCall() == 'yes' and re.ConstitutiveCall() == 'yes': const_call = 'yes'
                elif le.ConstitutiveCall() == 'yes' and re.ConstitutiveCall() == 'yes': const_call = 'yes'
                ens_exons=combineAnnotations([ens_exons1,ens_exons2])
                splice_event=combineAnnotations([le.AssociatedSplicingEvent(),re.AssociatedSplicingEvent()])
                splice_junctions=combineAnnotations([le.AssociatedSplicingJunctions(),re.AssociatedSplicingJunctions()])

                le_start_stop = [le.ExonStart(),le.ExonStop()]; le_start_stop.sort(); le_start,le_stop = le_start_stop
                re_start_stop = [re.ExonStart(),re.ExonStop()]; re_start_stop.sort(); re_start,re_stop = re_start_stop
                splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,le_start,le_stop)
                splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,re_start,re_stop)
                if len(splice_event)>0: const_call = 'no'
                values = [gene,exon1+'-'+exon2,'chr'+le.Chr(),le.Strand(),str(le_start)+'|'+str(le_stop),str(re_start)+'|'+str(re_stop),const_call,ens_exons,splice_event,splice_junctions]
                values = string.join(values,'\t')+'\n'; junctiondata.write(values)
                #print exon1+'-'+exon2, le_start_stop,re_start_stop
        if gene in excluded_intronic_junctions:
            ### Repeat for junctions that occur in AltAnalyze determined retained introns
            for (le,re) in excluded_intronic_junctions[gene]:
                if le.IntronDeletionStatus()== 'yes': exon1 = le.NewIntronRegion()
                else: exon1 = le.ExonRegionID2()
                if re.IntronDeletionStatus()== 'yes': exon2 = re.NewIntronRegion()
                else: exon2 = re.ExonRegionID2()

                ens_exons1 = le.ExonID()
                ens_exons2 = re.ExonID()
                const_call = 'no'
                ens_exons=combineAnnotations([ens_exons1,ens_exons2])
                splice_event=combineAnnotations([le.AssociatedSplicingEvent(),re.AssociatedSplicingEvent()])
                splice_junctions=combineAnnotations([le.AssociatedSplicingJunctions(),re.AssociatedSplicingJunctions()])

                le_start_stop = [le.ExonStart(),le.ExonStop()]; le_start_stop.sort(); le_start,le_stop = le_start_stop
                re_start_stop = [re.ExonStart(),re.ExonStop()]; re_start_stop.sort(); re_start,re_stop = re_start_stop
                splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,le_start,le_stop)
                splice_event = getUCSCSplicingAnnotations(ucsc_events,splice_event,re_start,re_stop)
                values = [gene,exon1+'-'+exon2,'chr'+le.Chr(),le.Strand(),str(le_start)+'|'+str(le_stop),str(re_start)+'|'+str(re_stop),const_call,ens_exons,splice_event,splice_junctions]
                values = string.join(values,'\t')+'\n'; junctiondata.write(values)
                #print exon1, le_start_stop, exon2, re_start_stop           
    sgvjdata.close(); exondata.close(); junctiondata.close()

def combineAnnotations(annotation_list):
    annotation_list2=[]
    for annotations in annotation_list:
        annotations = string.split(annotations,'|')
        for annotation in annotations:
            if len(annotation)>0: annotation_list2.append(annotation)
    annotation_list = unique.unique(annotation_list2)
    annotation_list = string.join(annotation_list,'|')
    return annotation_list

def getMatchingEnsExons(region_start,region_stop,exon_coord_db):
    region_coord=[region_start,region_stop]; ens_exons=[]
    region_coord.sort(); region_start,region_stop = region_coord
    for exon_coord in exon_coord_db:
        combined=region_coord+list(exon_coord)
        combined.sort()
        if region_start==combined[1] and region_stop==combined[-2]:
            ens_exons.append(exon_coord_db[exon_coord])
    ens_exons = string.join(ens_exons,'|')
    return ens_exons

################# Begin Analysis from parsing files
class EnsemblInformation:
    def __init__(self, chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, constitutive_exon, new_exon_start, new_exon_stop, new_gene_start, new_gene_stop):
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._genestart = gene_start; self._genestop = gene_stop
        self._exonstart = exon_start; self._exonstop = exon_stop
        self._constitutive_exon = constitutive_exon; self._strand = strand
        self._newgenestart = new_gene_start;  self._newgenestop = new_gene_stop
        self._newexonstart = new_exon_start; self._newexonstop = new_exon_stop
        self._del_status = 'no' #default value
        self._distal_intron_region = 1 #default value
        self.mx='no'
    def GeneID(self): return self._geneid
    def ExonID(self): return self._exonid
    def reSetExonID(self,exonid): self._exonid = exonid
    def Chr(self): return self._chr
    def Strand(self): return self._strand
    def GeneStart(self): return self._genestart
    def GeneStop(self): return self._genestop
    def ExonStart(self): return self._exonstart
    def ExonStop(self): return self._exonstop
    def NewGeneStart(self): return self._newgenestart
    def NewGeneStop(self): return self._newgenestop
    def NewExonStart(self): return self._newexonstart
    def NewExonStop(self): return self._newexonstop
    def setConstitutive(self,constitutive_exon): self._constitutive_exon = constitutive_exon
    def Constitutive(self): return self._constitutive_exon
    def ConstitutiveCall(self):
        if self.Constitutive() == '1': call = 'yes'
        else: call = 'no'
        return call
    def ConstitutiveCallAbrev(self):
        if self.Constitutive() == 'yes': call = 'yes'
        elif self.Constitutive() == '1': call = 'yes'
        else: call = 'no'
        return call
    def setSpliceData(self,splice_event,splice_junctions):
        self._splice_event = splice_event; self._splice_junctions = splice_junctions
    def setMutuallyExclusive(self): self.mx='yes'
    def setIntronDeletionStatus(self,del_status):
        self._del_status = del_status
    def setAssociatedSplicingEvent(self,splice_event): self._splice_event = splice_event
    def AssociatedSplicingEvent(self):
        if self.mx == 'yes':
            try:
                self._splice_event = string.replace(self._splice_event,'cassette-exon','mutually-exclusive-exon')
                self._splice_event = string.join(unique.unique(string.split(self._splice_event,'|')),'|')
            except AttributeError: return 'mutually-exclusive-exon'
        try: return self._splice_event
        except AttributeError: return ''
    def updateDistalIntronRegion(self): self._distal_intron_region+=1
    def DistalIntronRegion(self): return self._distal_intron_region
    def setNewIntronRegion(self,new_intron_region): self.new_intron_region = new_intron_region
    def NewIntronRegion(self): return self.new_intron_region
    def IntronDeletionStatus(self):
        try: return self._del_status
        except AttributeError: return 'no'
    def setAssociatedSplicingJunctions(self,splice_junctions): self._splice_junctions = splice_junctions
    def AssociatedSplicingJunctions(self):
        try: return self._splice_junctions
        except AttributeError: return ''
    def setExonRegionIDs(self,id): self._exon_region_ids = id
    def ExonRegionIDs(self): return self._exon_region_ids
    def setJunctionCoordinates(self,start1,stop1,start2,stop2):
        self.start1=start1;self.stop1=stop1;self.start2=start2;self.stop2=stop2
    def JunctionCoordinates(self): return [self.start1,self.stop1,self.start2,self.stop2]
    def JunctionDistance(self):
        jc = self.JunctionCoordinates(); jc.sort(); distance = int(jc[2])-int(jc[1])
        return distance
    def ExonNumber(self): return self._exon_num
    def RegionNumber(self): return self._region_num
    def ExonRegionNumbers(self): return (self._exon_num,self._region_num)
    def ExonRegionID(self): return 'E'+str(self._exon_num)+'-'+str(self._region_num)
    def ExonRegionID2(self): return 'E'+str(self._exon_num)+'.'+str(self._region_num)
    def IntronRegionID(self): return 'I'+str(self._exon_num)+'.1'
    def setExonSeq(self,exon_seq): self._exon_seq = exon_seq
    def ExonSeq(self): return self._exon_seq
    def setPrevExonSeq(self,prev_exon_seq): self._prev_exon_seq = prev_exon_seq
    def PrevExonSeq(self):
        try: return self._prev_exon_seq
        except Exception: return ''
    def setNextExonSeq(self,next_exon_seq): self._next_exon_seq = next_exon_seq
    def NextExonSeq(self):
        try: return self._next_exon_seq
        except Exception: return ''
    def setPrevIntronSeq(self,prev_intron_seq): self._prev_intron_seq = prev_intron_seq
    def PrevIntronSeq(self): return self._prev_intron_seq
    def setNextIntronSeq(self,next_intron_seq): self._next_intron_seq = next_intron_seq
    def NextIntronSeq(self): return self._next_intron_seq
    def setPromoterSeq(self,promoter_seq): self._promoter_seq = promoter_seq
    def PromoterSeq(self): return self._promoter_seq
    
    def AllGeneValues(self):
        output = str(self.ExonID())
        return output
    def __repr__(self): return self.AllGeneValues()

class ExonStructureData(EnsemblInformation):
    def __init__(self, ensembl_gene_id, chr, strand, exon_start, exon_stop, constitutive_exon, ensembl_exon_id, transcriptid):
        self._transcriptid = transcriptid
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._exonstart = exon_start; self._exonstop = exon_stop
        self._constitutive_exon = constitutive_exon; self._strand = strand
        self._distal_intron_region = 1; self.mx='no'
    def TranscriptID(self): return self._transcriptid

class ExonRegionData(EnsemblInformation):
    def __init__(self, ensembl_gene_id, chr, strand, exon_start, exon_stop, ensembl_exon_id, exon_region_id, exon_num, region_num, constitutive_exon):
        self._exon_region_id = exon_region_id; self._constitutive_exon = constitutive_exon
        self._geneid = ensembl_gene_id; self._exonid = ensembl_exon_id; self._chr = chr
        self._exonstart = exon_start; self._exonstop = exon_stop; self._distal_intron_region = 1; self.mx='no'
        self._exon_num = exon_num; self._region_num = region_num; self._strand = strand

class ProbesetAnnotation(EnsemblInformation):
    def __init__(self, ensembl_exon_id, constitutive_exon, exon_region_id, splice_event, splice_junctions, exon_start, exon_stop):
        self._region_num = exon_region_id; self._constitutive_exon = constitutive_exon;self._splice_event = splice_event
        self._splice_junctions = splice_junctions; self._exonid = ensembl_exon_id
        self._exonstart = exon_start; self._exonstop = exon_stop; self.mx='no'
        
class ExonAnnotationsSimple(EnsemblInformation):
    def __init__(self, chr, strand, exon_start, exon_stop, ensembl_gene_id, ensembl_exon_id ,constitutive_exon, exon_region_id, splice_event, splice_junctions):
        self._exon_region_ids = exon_region_id; self._constitutive_exon = constitutive_exon;self._splice_event =splice_event
        self._splice_junctions = splice_junctions; self._exonid = ensembl_exon_id; self._geneid = ensembl_gene_id
        self._chr = chr; self._exonstart = exon_start; self._exonstop = exon_stop; self._strand = strand; self.mx='no'

class CriticalExonInfo:
    def __init__(self,geneid,critical_exon,splice_type,junctions):
        self._geneid = geneid; self._junctions = junctions  
        self._critical_exon = critical_exon; self._splice_type = splice_type
    def GeneID(self): return self._geneid
    def Junctions(self): return self._junctions
    def CriticalExonRegion(self): return self._critical_exon
    def SpliceType(self): return self._splice_type

class RelativeExonLocations:
    def __init__(self,exonid,pes,pee,nes,nee):
        self._exonid = exonid; self._pes = pes; self._pee = pee
        self._nes = nes; self._nee = nee
    def ExonID(self): return self._exonid
    def PrevExonCoor(self): return (self._pes,self._pee)
    def NextExonCoor(self): return (self._nes,self._nee)
    def __repr__(self): return self.ExonID()
    
################### Import exon coordinate/transcript data from BIOMART
def importEnsExonStructureDataSimple(species,type,gene_strand_db,exon_location_db,adjacent_exon_locations):
    if type == 'ensembl': filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    elif type == 'ucsc': filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_mrna.txt'
    elif type == 'ncRNA': filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_ncRNA.txt'
    start_time = time.time()
    fn=filepath(filename); x=0; k=[]; relative_exon_locations={}
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
            if strand == '-1' or strand == '-': strand = '-'#; exon_start2 = int(exon_end); exon_end2 = int(exon_start); exon_start=exon_start2; exon_end=exon_end2
            else: strand = '+'; exon_end = int(exon_end)#; exon_start = int(exon_start)
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            exon_end = int(exon_end); exon_start = int(exon_start)
            exon_data = (exon_start,exon_end,ens_exonid)
            try: relative_exon_locations[ens_transcriptid,gene,strand].append(exon_data)
            except KeyError: relative_exon_locations[ens_transcriptid,gene,strand] = [exon_data]
            gene_strand_db[gene] = strand
            exon_location_db[ens_exonid] = exon_start,exon_end

    ###Generate a list of exon possitions for adjacent exons for all exons
    first_exon_dbase={}
    for (transcript,gene,strand) in relative_exon_locations:
        relative_exon_locations[(transcript,gene,strand)].sort()
        if strand == '-': relative_exon_locations[(transcript,gene,strand)].reverse()
        i = 0
        ex_ls = relative_exon_locations[(transcript,gene,strand)]
        for exon_data in ex_ls:
            exonid = exon_data[-1]
            if i == 0: ### first exon
                pes = -1; pee = -1 ###Thus, Index should be out of range, but since -1 is valid, it won't be
                if strand == '-': ces = ex_ls[i][1]
                else: ces = ex_ls[i][0]
                try: first_exon_dbase[gene].append([ces,exonid])
                except KeyError: first_exon_dbase[gene] = [[ces,exonid]]
            else: pes = ex_ls[i-1][0]; pee = ex_ls[i-1][1] ###pes: previous exon start, pee: previous exon end
            try:  nes = ex_ls[i+1][0]; nee = ex_ls[i+1][1]
            except IndexError: nes = -1; nee = -1
            rel = RelativeExonLocations(exonid,pes,pee,nes,nee)
            """if exonid in adjacent_exon_locations:
                rel1 = adjacent_exon_locations[exonid]
                prev_exon_start,prev_exon_stop = rel1.NextExonCoor()
                next_exon_start,next_exon_stop = rel1.PrevExonCoor()
                if prev_exon_start == -1 or next_exon_start == -1:
                    adjacent_exon_locations[exonid] = rel ###Don't over-ride the exisitng entry if no exon is proceeding or following
            else: adjacent_exon_locations[exonid] = rel"""
            adjacent_exon_locations[exonid] = rel
            i+=1

    for gene in first_exon_dbase:
        first_exon_dbase[gene].sort()
        strand = gene_strand_db[gene]
        if strand == '-': first_exon_dbase[gene].reverse()
        first_exon_dbase[gene] = first_exon_dbase[gene][0][1] ### select the most 5' of the start exons for the gene
        
    #print relative_exon_locations['ENSMUST00000025142','ENSMUSG00000024293','-']; kill
    end_time = time.time(); time_diff = int(end_time-start_time)
    print filename,"parsed in %d seconds" % time_diff
    print len(gene_strand_db),'genes imported'
    return gene_strand_db,exon_location_db,adjacent_exon_locations,first_exon_dbase

def importEnsExonStructureDataSimpler(species,type,relative_exon_locations):
    if type == 'ensembl': filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    elif type == 'ucsc': filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_COMPLETE-mrna.txt'
    start_time = time.time()
    fn=filepath(filename); x=0; k=[]
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
            #if gene == 'ENSG00000143776'
            ###switch exon-start and stop if in the reverse orientation
            if strand == '-1' or strand == '-': strand = '-'#; exon_start2 = int(exon_end); exon_end2 = int(exon_start); exon_start=exon_start2; exon_end=exon_end2
            else: strand = '+'; exon_end = int(exon_end)#; exon_start = int(exon_start)
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            exon_end = int(exon_end); exon_start = int(exon_start)
            exon_data = (exon_start,exon_end,ens_exonid)
            try: relative_exon_locations[ens_transcriptid,gene,strand].append(exon_data)
            except KeyError: relative_exon_locations[ens_transcriptid,gene,strand] = [exon_data]

    return relative_exon_locations

def importEnsGeneData(species):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    global ens_gene_db; ens_gene_db={}
    importEnsExonStructureData(filename,species,'gene')
    return ens_gene_db
    
def importEnsExonStructureData(filename,species,data2process):
    start_time = time.time()
    fn=filepath(filename); x=0; k=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if data2process == 'all':
                ###switch exon-start and stop if in the reverse orientation
                if strand == '-1' or strand == '-': strand = '-'; exon_start2 = int(exon_end); exon_end2 = int(exon_start); exon_start=exon_start2; exon_end=exon_end2
                else: strand = '+'; exon_end = int(exon_end); exon_start = int(exon_start)
                if 'ENS' in gene: ens_exonid_data = string.split(ens_exonid,'.'); ens_exonid = ens_exonid_data[0]
                if abs(exon_end-exon_start)>-1:
                    if abs(exon_end-exon_start)<1:
                        try: too_short[gene].append(exon_start)
                        except KeyError: too_short[gene] = [exon_start]
                    exon_coordinates = [exon_start,exon_end]
                    continue_analysis = 'no'
                    if test=='yes': ###used to test the program for a single gene
                        if gene in test_gene: continue_analysis='yes' 
                    else: continue_analysis='yes'
                    if  continue_analysis=='yes':
                        ###Create temporary databases storing just exon and just trascript and then combine in the next block of code
                        initial_exon_annotation_db[ens_exonid] = gene,chr,strand,exon_start,exon_end,constitutive_exon                    
                        try: exon_transcript_db[ens_exonid].append(ens_transcriptid)
                        except KeyError: exon_transcript_db[ens_exonid] = [ens_transcriptid]
                        ###Use this database to figure out which ensembl exons represent intron retention as a special case down-stream
                        try: transcript_exon_db[gene,chr,strand,ens_transcriptid].append(exon_coordinates)
                        except KeyError: transcript_exon_db[gene,chr,strand,ens_transcriptid] = [exon_coordinates]
                        ###Store transcript data for downstream analyses
                        transcript_gene_db[ens_transcriptid] = gene,chr,strand
                        try: gene_transcript[gene].append(ens_transcriptid)
                        except KeyError: gene_transcript[gene] = [ens_transcriptid]
                        #print ens_exonid, exon_end
            elif data2process == 'exon-transcript':
                continue_analysis = 'yes'
                if test=='yes': ###used to test the program for a single gene
                    if gene not in test_gene: continue_analysis='no'
                if continue_analysis == 'yes':
                    try: exon_transcript_db[ens_exonid].append(ens_transcriptid)
                    except KeyError: exon_transcript_db[ens_exonid] = [ens_transcriptid]
            elif data2process == 'gene': ens_gene_db[gene] = chr,strand
    end_time = time.time(); time_diff = int(end_time-start_time)
    try: print len(transcript_gene_db), "number of transcripts included"
    except Exception: null=[]
    print filename,"parsed in %d seconds" % time_diff

def getEnsExonStructureData(species,data_type):
    start_time = time.time()
    ###Simple function to import and organize exon/transcript data
    filename1 = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    filename2 = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_mrna.txt'
    filename3 = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_filtered_ncRNA.txt'
    data2process = 'all'
    global initial_exon_annotation_db; initial_exon_annotation_db={}
    global exon_transcript_db; exon_transcript_db={}
    global transcript_gene_db; transcript_gene_db={}
    global gene_transcript; gene_transcript={}
    global transcript_exon_db; transcript_exon_db={}
    global initial_junction_db; initial_junction_db = {}; global too_short; too_short={}
    global ensembl_annotations; ensembl_annotations={}; global ensembl_gene_coordinates; ensembl_gene_coordinates = {}
    if data_type == 'mRNA' or data_type == 'gene':
        importEnsExonStructureData(filename2,species,data2process)
        importEnsExonStructureData(filename1,species,data2process)
    elif data_type == 'ncRNA':  ###Builds database based on a mix of Ensembl, GenBank and UID ncRNA IDs
        importEnsExonStructureData(filename3,species,data2process)
        
    global exon_annotation_db; exon_annotation_db={} ###Agglomerate the data from above and store as an instance of the ExonStructureData class
    for ens_exonid in initial_exon_annotation_db:
        gene,chr,strand,exon_start,exon_end,constitutive_exon = initial_exon_annotation_db[ens_exonid]
        ens_transcript_list = exon_transcript_db[ens_exonid]
        ###Record the coordiantes for matching up UCSC exon annotations with Ensembl genes
        try: ensembl_gene_coordinates[gene].append(exon_start)
        except KeyError: ensembl_gene_coordinates[gene] = [exon_start]
        ensembl_gene_coordinates[gene].append(exon_end)
        ensembl_annotations[gene] = chr,strand
        
        y = ExonStructureData(gene, chr, strand, exon_start, exon_end, constitutive_exon, ens_exonid, ens_transcript_list)
        exon_info = [exon_start,exon_end,y]
        if gene in too_short: too_short_list = too_short[gene]
        else: too_short_list=[]
        if exon_start in too_short_list:# and exon_end in too_short_list: pass ###Ensembl exons that are one bp (start and stop the same) - fixed minor issue in 2.0.9
            pass
        else:
            try: exon_annotation_db[(gene,chr,strand)].append(exon_info)
            except KeyError: exon_annotation_db[(gene,chr,strand)] = [exon_info]

    initial_exon_annotation_db={}; exon_transcript_db={}
    print 'Exon and transcript data obtained for %d genes' % len(exon_annotation_db)
    ###Grab the junction data
    for (gene,chr,strand,transcript) in transcript_exon_db:
        exon_list_data = transcript_exon_db[(gene,chr,strand,transcript)];exon_list_data.sort()
        if strand == '-': exon_list_data.reverse()
        index=0
        if gene in too_short: too_short_list = too_short[gene]
        else: too_short_list=[]
        while index+1 <len(exon_list_data):
            junction_data = exon_list_data[index],exon_list_data[index+1]
            if junction_data[0][0] not in too_short_list and junction_data[0][1] not in too_short_list and junction_data[1][0] not in too_short_list and junction_data[1][1] not in too_short_list: ###Ensembl exons that are one bp (start and stop the same)
                try: initial_junction_db[(gene,chr,strand)].append(junction_data)
                except KeyError: initial_junction_db[(gene,chr,strand)] = [junction_data]
            index+=1
    print 'Exon-junction data obtained for %d genes' % len(initial_junction_db)
    
    ###Find exons that suggest intron retention: Within a junction database, find exons that overlapp with junction boundaries
    intron_retention_db={}; retained_intron_exons={}; delete_db={}
    for key in initial_junction_db:
        for junction_info in initial_junction_db[key]:
            e5_pos = junction_info[0][1]; e3_pos = junction_info[1][0] ###grab the junction coordiantes
            for exon_info in exon_annotation_db[key]:
                exon_start,exon_stop,ed = exon_info
                loc = [e5_pos,e3_pos]; loc.sort() ###The downstream functions need these two sorted
                new_exon_info = loc[0],loc[1],ed
                retained = compareExonLocations(e5_pos,e3_pos,exon_start,exon_stop)
                exonid = ed.ExonID()
                if retained == 'yes':
                    #print exonid,e5_pos,e3_pos,exon_start,exon_stop
                    intron_length = abs(e5_pos-e3_pos)
                    if intron_length>500:
                        ed.setIntronDeletionStatus('yes')
                        try: delete_db[key].append(exon_info)
                        except KeyError: delete_db[key] = [exon_info]
                    try: intron_retention_db[key].append(new_exon_info)
                    except KeyError: intron_retention_db[key]=[new_exon_info]
                    retained_intron_exons[exonid]=[]
                    #print key,ed.ExonID(),len(exon_annotation_db[key])

    k=0 ### Below code removes exons that have been classified as retained introns - not needed when we can selective remove these exons with ed.IntronDeletionStatus()
    #exon_annotation_db2={}
    for key in exon_annotation_db:
        for exon_info in exon_annotation_db[key]:
            if key in delete_db:
                delete_info = delete_db[key] ### coordinates and objects... looks like you can match up based on object memory locations
                if exon_info in delete_info: k+=1
                """else:
                    try: exon_annotation_db2[key].append(exon_info)
                    except KeyError: exon_annotation_db2[key]=[exon_info]"""
            """else:
                try: exon_annotation_db2[key].append(exon_info)
                except KeyError: exon_annotation_db2[key]=[exon_info]"""
    #exon_annotation_db = exon_annotation_db2; exon_annotation_db2=[]
    transcript_exon_db=[]

    print k, 'exon entries removed from primary exon structure, which occur in predicted retained introns'    
    initial_junction_db={}
    print len(retained_intron_exons),"ensembl exons in %d genes, show evidence of being retained introns (only sequences > 500bp are removed from the database" % len(intron_retention_db)

    end_time = time.time(); time_diff = int(end_time-start_time)
    print "Primary databases built in %d seconds" % time_diff

    ucsc_splicing_annot_db = alignToKnownAlt.importEnsExonStructureData(species,ensembl_gene_coordinates,ensembl_annotations,exon_annotation_db) ### Should be able to exclude
    #except Exception: ucsc_splicing_annot_db={}

    return exon_annotation_db,transcript_gene_db,gene_transcript,intron_retention_db,ucsc_splicing_annot_db

def compareExonLocations(e5_pos,e3_pos,exon_start,exon_stop):
    sort_list = [e5_pos,e3_pos,exon_start,exon_stop]; sort_list.sort()
    new_sort = sort_list[1:-1]
    if e5_pos in new_sort and e3_pos in new_sort: retained = 'yes'
    else: retained = 'no'
    return retained
    
################### Import exon sequence data from BIOMART (more flexible alternative function to above)
def getSeqLocations(sequence,ed,strand,start,stop,original_start,original_stop):
    cd = seqSearch(sequence,ed.ExonSeq())
    if cd != -1:
        if strand == '-': exon_stop = stop - cd; exon_start = exon_stop - len(ed.ExonSeq()) + 1
        else: exon_start = start + cd; exon_stop = exon_start + len(ed.ExonSeq())-1
        ed.setExonStart(exon_start); ed.setExonStop(exon_stop); ed.setGeneStart(original_start); ed.setGeneStop(original_stop)
        #if ed.ExonID() == 'E10' and ed.ArrayGeneID() == 'G7225860':
        #print exon_start, exon_stop,len(ed.ExonSeq()),ed.ExonSeq();kill
    else:
        cd = seqSearch(sequence,ed.ExonSeq()[:15])
        #print ed.ExonSeq()[:15],ed.ExonSeq();kill
        if cd == -1: cd = seqSearch(sequence,ed.ExonSeq()[-15:])
        if cd != -1:
            if strand == '-': exon_stop = stop - cd; exon_start = exon_stop - len(ed.ExonSeq()) + 1
            else: exon_start = start + cd; exon_stop = exon_start + len(ed.ExonSeq())-1
            ed.setExonStart(exon_start); ed.setExonStop(exon_stop); ed.setGeneStart(original_start); ed.setGeneStop(original_stop)
        else:
            ed.setGeneStart(original_start); ed.setGeneStop(original_stop) ### set these at a minimum for HTA arrays so that the pre-set exon coordinates are reported
    return ed

def import_sequence_data(filename,filter_db,species,analysis_type):
    """Note: A current bug in this module is that the last gene is not analyzed"""
    print "Begining generic fasta import of",filename
    fn=filepath(filename);fasta = {}; exon_db = {}; gene_db = {}; cDNA_db = {};sequence = '';count = 0; seq_assigned=0
    global temp_seq; temp_seq=''; damned =0; global failed; failed={}
    addition_seq_len = 2000; var = 1000
    if len(analysis_type)==2: analysis_parameter,analysis_type = analysis_type
    else: analysis_parameter = 'null'
    if 'gene' in fn or 'chromosome' in fn:
        gene_strand_db,exon_location_db,adjacent_exon_locations,null = importEnsExonStructureDataSimple(species,'ucsc',{},{},{})
        gene_strand_db,exon_location_db,adjacent_exon_locations,first_exon_db = importEnsExonStructureDataSimple(species,'ensembl',gene_strand_db,exon_location_db,adjacent_exon_locations)
        null=[]
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        try:
            if data[0] == '>':
                    if len(sequence) > 0:
                        if 'gene' in fn or 'chromosome' in fn:
                            start = int(start); stop = int(stop)
                        else:
                            try: exon_start = int(exon_start); exon_stop = int(exon_stop)
                            except ValueError: exon_start = exon_start
                            if strand == '-1': strand = '-'
                            else: strand = '+'
                            if strand == '-': new_exon_start = exon_stop; new_exon_stop = exon_start; exon_start = new_exon_start; exon_stop = new_exon_stop #"""
                        if 'exon' in fn:
                            exon_info = [exon_start,exon_stop,exon_id,exon_annot]
                            try: exon_db[(gene,chr,strand)].append(exon_info)
                            except KeyError: exon_db[(gene,chr,strand)] = [exon_info] #exon_info = (exon_start,exon_stop,exon_id,exon_annot)
                            fasta[geneid]=description
                        if 'cDNA' in fn:
                            cDNA_info = [transid,sequence]
                            try: cDNA_db[(gene,strand)].append(cDNA_info)
                            except KeyError: cDNA_db[(gene,strand)] = [cDNA_info]                   
                        if 'gene' in fn or 'chromosome' in fn:
                            temp_seq = sequence
                            if analysis_type == 'gene_count': fasta[gene]=[]
                            if gene in filter_db:
                                count += 1
                                if count == var: print var,'genes examined...'; var+=1000
                                #print gene, filter_db[gene][0].ArrayGeneID();kill
                                if (len(sequence) -(stop-start)) != ((addition_seq_len*2) +1):
                                    ###multiple issues can occur with trying to grab sequence up and downstream - for now, exlcude these
                                    ###some can be accounted for by being at the end of a chromosome, but not all.
                                    damned +=1
                                    try: failed[chr]+=1
                                    except KeyError: failed[chr]=1
                                    """if chr in failed:
                                        gene_len = (stop-start); new_start = start - addition_seq_len
                                        null = sequence[(gene_len+2000):]"""
                                else:
                                    original_start = start; original_stop = stop
                                    start = start - addition_seq_len
                                    stop = stop + addition_seq_len
                                    strand = gene_strand_db[gene]
                                    first_exonid = first_exon_db[gene]
                                    fexon_start,fexon_stop = exon_location_db[first_exonid]
                                    for ed in filter_db[gene]:
                                        if analysis_type == 'get_locations':
                                            ed = getSeqLocations(sequence,ed,strand,start,stop,original_start,original_stop)
                                        if analysis_type == 'get_sequence':
                                            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = ed
                                            if analysis_parameter == 'region_only': ens_exon_list = [exon_id]
                                            else: ens_exon_list = ed.ExonID()
                                            for ens_exon in ens_exon_list:
                                                if len(ens_exon)>0:
                                                    if analysis_parameter == 'region_only':
                                                        ### Only extract the specific region exon sequence
                                                        exon_start,exon_stop = int(probe_start),int(probe_stop)
                                                        #if ':440657' in probeset_id: print ens_exon_list,probe_start,probe_stop
                                                        #probe_coord = [int(probe_start),int(probe_stop)]; probe_coord.sort()
                                                        #exon_start,exon_stop = probe_coord
                                                    else: exon_start,exon_stop = exon_location_db[ens_exon]
                                                    #if ':440657' in probeset_id: print [exon_start,exon_stop]
                                                    exon_sequence = grabSeq(sequence,strand,start,stop,exon_start,exon_stop,'exon')
                                                    #print [exon_id,exon_sequence, start,stop,exon_start,exon_stop];kill
                                                    """Could repeat if we build another dictionary with exon->adjacent exon positions (store in a class where
                                                    you designate last and next exon posiitons for each transcript relative to that exon), to grab downstream, upsteam exon and intron sequences"""
                                                    try:
                                                        rel = adjacent_exon_locations[ens_exon]
                                                        prev_exon_start,prev_exon_stop = rel.PrevExonCoor()
                                                        next_exon_start,next_exon_stop = rel.NextExonCoor()
                                                        prev_exon_sequence = grabSeq(sequence,strand,start,stop,prev_exon_start,prev_exon_stop,'exon')
                                                        next_exon_sequence = grabSeq(sequence,strand,start,stop,next_exon_start,next_exon_stop,'exon')
                                                        seq_type = 'intron'
                                                        if strand == '-':
                                                            if 'alt-N-term' in ed.AssociatedSplicingEvent() or 'altPromoter' in ed.AssociatedSplicingEvent(): seq_type = 'promoter' ### Thus prev_intron_seq is used to designate an alternative promoter sequence
                                                            prev_intron_sequence = grabSeq(sequence,strand,start,stop,exon_stop,prev_exon_start,seq_type)
                                                            promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                            next_intron_sequence = grabSeq(sequence,strand,start,stop,next_exon_stop,exon_start,'intron')
                                                        else:
                                                            prev_intron_sequence = grabSeq(sequence,strand,start,stop,prev_exon_stop,exon_start,seq_type)
                                                            promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                            next_intron_sequence = grabSeq(sequence,strand,start,stop,exon_stop,next_exon_start,'intron')
                                                        """if 'ENS' in ens_exon:
                                                            print ens_exon, strand
                                                            print '1)',exon_sequence
                                                            print '2)',prev_intron_sequence[:20],prev_intron_sequence[-20:], len(prev_intron_sequence), strand,start,stop,prev_exon_stop,exon_start,seq_type, ed.AssociatedSplicingEvent()
                                                            print '3)',next_intron_sequence[:20],next_intron_sequence[-20:], len(next_intron_sequence)
                                                            print '4)',promoter_sequence[:20],promoter_sequence[-20:], len(promoter_sequence);kill"""
                                                        ###Intron sequences can be extreemly long so just include the first and last 1kb
                                                        if len(prev_intron_sequence)>2001 and seq_type == 'promoter': prev_intron_sequence = prev_intron_sequence[-2001:]
                                                        elif len(prev_intron_sequence)>2001: prev_intron_sequence = prev_intron_sequence[:1000]+'|'+prev_intron_sequence[-1000:]
                                                        if len(next_intron_sequence)>2001: next_intron_sequence = next_intron_sequence[:1000]+'|'+next_intron_sequence[-1000:]
                                                        if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]
                                                    except Exception:
                                                        ### When analysis_parameter == 'region_only' an exception is desired since we only want exon_sequence
                                                        prev_exon_sequence=''; next_intron_sequence=''; next_exon_sequence=''; promoter_sequence=''
                                                        if analysis_parameter != 'region_only':
                                                            if strand == '-': promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                            else: promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                            if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]                                                     
                                                    ed.setExonSeq(exon_sequence) ### Use to replace the previous probeset/critical exon sequence with sequence corresponding to the full exon
                                                    seq_assigned+=1
                                                    if len(prev_exon_sequence)>0:
                                                        ### Use to output sequence for ESE/ISE type motif searches
                                                        ed.setPrevExonSeq(prev_exon_sequence);
                                                    if len(next_exon_sequence)>0: ed.setNextExonSeq(next_exon_sequence)
                                                    if analysis_parameter != 'region_only':
                                                        ed.setPrevIntronSeq(prev_intron_sequence[1:-1]); ed.setNextIntronSeq(next_intron_sequence[1:-1])
                                                        ed.setPromoterSeq(promoter_sequence[1:-1])
                                                else:
                                                    if strand == '-': promoter_sequence = grabSeq(sequence,strand,start,stop,fexon_stop,-1,"promoter")
                                                    else: promoter_sequence = grabSeq(sequence,strand,start,stop,-1,fexon_start,"promoter")
                                                    if len(promoter_sequence)>2001: promoter_sequence = promoter_sequence[-2001:]
                                                    ed.setPromoterSeq(promoter_sequence[1:-1])
                            sequence = ''; data2 = data[1:]; t= string.split(data2,'|'); gene,chr,start,stop = t
                    else:
                        data2 = data[1:]; t= string.split(data2,'|'); gene,chr,start,stop = t
        except IndexError: continue
        try:
            if data[0] != '>': sequence = sequence + data
        except IndexError: continue

    if analysis_type == 'get_locations': ### Applies to the last gene sequence read
        if ('gene' in fn or 'chromosome' in fn) and len(sequence) > 0:
            start = int(start); stop = int(stop)
            if analysis_type == 'gene_count': fasta[gene]=[]
            if gene in filter_db:
                    original_start = start; original_stop = stop
                    start = start - addition_seq_len
                    stop = stop + addition_seq_len
                    strand = gene_strand_db[gene]
                    first_exonid = first_exon_db[gene]
                    fexon_start,fexon_stop = exon_location_db[first_exonid]
                    for ed in filter_db[gene]:
                        try: print ed.ExonStart(),'1'
                        except Exception: null=[]
                        ed = getSeqLocations(sequence,ed,strand,start,stop,original_start,original_stop)
                    for ed in filter_db[gene]:
                        print ed.ExonStart(),'2'

    probesets_analyzed=0
    for gene in filter_db:
        for probe_data in filter_db[gene]: probesets_analyzed+=1
        
    print "Number of imported sequences:", len(fasta),count
    print "Number of assigned probeset sequences:",seq_assigned,"out of",probesets_analyzed
    if len(exon_db) > 0: return exon_db,fasta
    elif len(cDNA_db) > 0: return cDNA_db
    elif len(fasta) > 0: return fasta
    else: return filter_db

def grabSeq(sequence,strand,start,stop,exon_start,exon_stop,type):
    proceed = 'yes'
    if type != 'promoter' and (exon_start == -1 or exon_stop == -1): proceed = 'no' ###Thus no preceeding or succedding exons and thus no seq reported
    if proceed == 'yes':
        if strand == '-':
            if exon_stop == -1: exon_stop = stop
            exon_sequence = sequence[(stop-exon_stop):(stop-exon_stop)+(exon_stop-exon_start)+1]
        else:
            #print type, exon_start,start,exon_stop
            if exon_start == -1: exon_start = start ### For last intron 
            exon_sequence = sequence[(exon_start-start):(exon_stop-start+1)]
    else: exon_sequence = ''
    return exon_sequence

def seqSearch(sequence,exon_seq):
    cd = string.find(sequence,exon_seq)
    if cd == -1:
        rev_seq = reverse_orientation(exon_seq); cd = string.find(sequence,rev_seq)
    return cd

def reverse_string(astring):
    "http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/65225"
    revchars = list(astring)        # string -> list of chars
    revchars.reverse()              # inplace reverse the list
    revchars = ''.join(revchars)    # list of strings -> string
    return revchars

def reverse_orientation(sequence):
    """reverse the orientation of a sequence (opposite strand)"""
    exchange = []
    for nucleotide in sequence:
        if nucleotide == 'A': nucleotide = 'T'
        elif nucleotide == 'T': nucleotide = 'A'
        elif nucleotide == 'G': nucleotide = 'C'
        elif nucleotide == 'C': nucleotide = 'G'
        exchange.append(nucleotide)
    complementary_sequence = reverse_string(exchange)
    return complementary_sequence

############# First pass for annotating exons into destict, ordered regions for further annotation
def annotate_exons(exon_location):
    print "Begining to assign intial exon block and region annotations"
    ### Sort and reverse exon orientation for transcript_cluster exons
    original_neg_strand_coord={}
    ###make negative strand coordinates look like positive strand to identify overlapping exons
    for (geneid,chr,strand) in exon_location:
        exon_location[(geneid,chr,strand)].sort()
        if strand == '-':
            exon_location[(geneid,chr,strand)].reverse()
            denominator = exon_location[(geneid,chr,strand)][0][0] ###numerical coordiantes to subtract from to normalize negative strand data
            for exon_info in exon_location[(geneid,chr,strand)]:
                start,stop,ed = exon_info; ens_exon = ed.ExonID()
                coordinates = stop,start; coordinates = copy.deepcopy(coordinates)###format these in reverse for analysis in other modules
                original_neg_strand_coord[ens_exon] = coordinates
                exon_info[0] = abs(start - denominator);exon_info[1] = abs(stop - denominator)
 
    #alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z']
    exon_location2={}; exon_temp_list=[]
    for key in exon_location:
        index = 1; index2 = 1; exon_list=[]; y = 0 #if key[-1] == '-':
        for exon in exon_location[key]: #print exon[0],exon[1],len(exon_temp_list)
            if exon[-1].IntronDeletionStatus() == 'yes': null=[] ### retained intron (don't include)
            else:
                if y == 0:
                    exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                    exon_list.append(exon_info); y = 1; last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list =[]; exon_temp_list.append(last_start); exon_temp_list.append(last_stop)
                elif y == 1:
                    current_start = exon[0];  current_stop = exon[1]
                    if ((current_start >= last_start) and (current_start <= last_stop)) or ((last_start >= current_start) and (last_start <= current_stop)):
                        exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)] #+alphabet[index2]
                        exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                    elif (abs(current_start - last_stop) < 1) or (abs(last_start - current_stop) < 1):
                        exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)] #+alphabet[index2]
                        exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                    elif len(exon_temp_list)>3:
                        exon_temp_list.sort()
                        if (((current_start-1) > exon_temp_list[-1]) and ((current_stop-1) > exon_temp_list[-1])) or (((current_start+1) < exon_temp_list[0]) and ((current_stop+1) < exon_temp_list[0])):
                            ###Thus an overlapp with atleast one exon DOESN'T exists
                            exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                            exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list=[]; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                        else:
                            exon_info = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)]  #+alphabet[index2]
                            exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index2 += 1; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
                    else:
                        exon_info = ['E'+str(index)+'-1',exon,(index,1)]
                        exon_list.append(exon_info); last_start = exon[0]; last_stop = exon[1]; index += 1; index2 = 2; exon_temp_list=[]; exon_temp_list.append(current_start); exon_temp_list.append(current_stop)
        exon_location2[key] = exon_list

    for key in exon_location2:
        ###Re-assign exon coordiantes back to the actual coordiantes for negative strand genes
        strand = key[-1]
        if strand == '-':
            for exon_data in exon_location2[key]:
                ed = exon_data[1][2]; ens_exon = ed.ExonID()
                ###re-assing the coordiantes
                start,stop = original_neg_strand_coord[ens_exon]
                exon_data[1][0] = start; exon_data[1][1] = stop
    """
    for key in exon_location2:
        if key[0] == 'ENSG00000129566':
            for item in exon_location2[key]: print key, item"""
    return exon_location2

def exon_clustering(exon_location):
    """for i in exon_location[('ENSG00000129566','14','-')]:print i"""
    #exon_info = exon_start,exon_end,y
    #try: exon_annotation_db[(geneid,chr,strand)].append(exon_info)
    #except KeyError: exon_annotation_db[(geneid,chr,strand)] = [exon_info]
    
    exon_clusters={}; region_locations={}; region_gene_db={}
    for key in exon_location:
        chr = 'chr'+key[1];entries = exon_location[key];strand = key[-1]; gene = key[0]
        temp_exon_db={}; temp_exon_db2={}; exon_block_db={}
        for exon in entries:
            a = string.find(exon[0],'-')  ###outdated code: if all have a '-'
            if a == -1: exon_number = int(exon[0][1:])
            else: exon_number = int(exon[0][1:a])
            ###Group overlapping exons into exon clusters
            try: temp_exon_db[exon_number].append(exon[1])
            except KeyError: temp_exon_db[exon_number] = [exon[1]]
            ###add all start and stop values to the temp database
            try: temp_exon_db2[exon_number].append(exon[1][0])
            except KeyError: temp_exon_db2[exon_number] = [exon[1][0]]
            temp_exon_db2[exon_number].append(exon[1][1])
            try: exon_block_db[exon_number].append(exon)
            except KeyError: exon_block_db[exon_number] = [exon]
        for exon in temp_exon_db:
            #intron = 'I'+str(exon)+'-1'
            exon_info = unique.unique(temp_exon_db2[exon])
            exon_info.sort();start = exon_info[0];stop = exon_info[-1];type=[]; accession=[]
            for (exon_start,exon_stop,ed) in temp_exon_db[exon]:
                exon_type = ed.Constitutive()
                exon_id = ed.ExonID()
                type.append(exon_type); accession.append(exon_id)
                #if exon_id == 'ENSE00000888906': print exon_info,temp_exon_db[exon]
            type=unique.unique(type); accession=unique.unique(accession)
            exon_data = exon,(start,stop),accession,type
            key1 = key[0],chr,strand
            try: exon_clusters[key1].append(exon_data)
            except KeyError: exon_clusters[key1] = [exon_data]
            #if len(exon_info)-len(temp_exon_db[exon])>2: print key1,len(exon_info),len(temp_exon_db[exon]);kill
            if len(exon_info)>2:
                if strand == '-': exon_info.reverse()
                index=0; exon_data_list=[]
                while index < (len(exon_info)-1):
                    region = str(index+1)#;region_locations[key,'E'+str(exon)+'-'+region] = exon_info[index:index+2]
                    if strand == '-': new_stop,new_start = exon_info[index:index+2]
                    else: new_start,new_stop = exon_info[index:index+2]
                    ned = ExonStructureData(gene, key[1], strand, new_start, new_stop, '', '', [])
                    new_exon_info = ['E'+str(exon)+'-'+region,[new_start,new_stop,ned],(exon,index+1)]
                    exon_data_list.append(new_exon_info)
                    index+=1
                    #if gene == 'ENSG00000171735':print new_exon_info, 0
                region_locations[key,exon] = exon_data_list
            else:
                exon_data_list = [exon_block_db[exon][0]]  ###There can be multiples that occur - 2 exons 1 set of coordinates
                #if gene == 'ENSG00000171735':print exon_data_list, 1
                region_locations[key,exon] = exon_data_list
                
    ###Resort and re-populated the new exon_location database where we've re-defined the region entries    
    interim_location_db={}
    for (key,exon) in region_locations:
        exon_data_list = region_locations[(key,exon)]
        for exon_data in exon_data_list:
            try: interim_location_db[key].append((exon,exon_data))
            except KeyError: interim_location_db[key] = [(exon,exon_data)]
    for key in interim_location_db:
        interim_location_db[key].sort(); new_exon_list=[]
        for (e,i) in interim_location_db[key]: new_exon_list.append(i)
        exon_location[key] = new_exon_list

    #for i in exon_location[('ENSG00000171735', '1', '+')]:print i

    """
    for i in region_locations:
        if 'ENSG00000129566' in i[0]:print i, region_locations[i]
    
    ###Transform coordinates from the source Ensembl exon to discrete regions (regions may not be biological relevant in the 3' or 5' exon of the gene due to EST assemblies).
    for (key,exon_id) in region_locations:
        gene,chr,strand = key; id_added='no'; exon_annot=[]; ens_exonids=[]; ens_transcripts=[]
        if strand == '-': stop,start = region_locations[(key,exon_id)]
        else: start,stop = region_locations[(key,exon_id)]
        ###If the number of old regions is greater than the new, delete the old
        previous_region_number = len(exon_location[key])
        new_region_number = region_gene_db[key]
        if previous_region_number>new_region_number: 
        for exon_data in exon_location[key]:
            exon_start,exon_stop,ed = exon_data[1]; ens_exon_id = ed.ExonID(); exon_annotation = ed.Constitutive()
            #if exon_id == exon_data[0]: exon_data[1][0] = start; exon_data[1][1] = stop; id_added = 'yes'
            if exon_start == start or exon_stop == stop: exon_annot.append(exon_annotation); ens_exonids.append(ens_exon_id); ens_transcripts+=ed.TranscriptID()
            if exon_stop == start or exon_start == stop: exon_annot.append(exon_annotation); ens_exonids.append(ens_exon_id)
        exon_annot = unique.unique(exon_annot); ens_exonids = unique.unique(ens_exonids); ens_transcripts = unique.unique(ens_transcripts)
        exon_annot = string.join(exon_annot,'|'); ens_exonids = string.join(ens_exonids,'|')
        for exon_data in exon_location[key]:
            if exon_id == exon_data[0]:
                ###Replace exsting entries (we're basically replacing these with completely new data, just like below, but this is easier than deleting the existing)
                y = ExonStructureData(gene, chr, strand, start, stop, exon_annot, ens_exonids, ens_transcripts)
                exon_data[1][0] = start; exon_data[1][1] = stop; exon_data[1][2] = y; id_added = 'yes'
                #start is the lower number, with either strand
                break
        if id_added == 'no': ###This occurs when a new region must be introduced from a large now broken large exon (e.g. E1-1 pos: 1-300, E1-2 pos: 20-200, now need a 3rd E1-3 200-300)
            indeces = string.split(exon_id[1:],'-')
            index1 = int(indeces[0]); index2 = int(indeces[1])
            #new_entry = ['E'+str(index-1) +'-'+ str(index2),exon,(index-1,index2)]
            #exon_info = [exon_start,exon_stop,exon_id,exon_annot]
            ###Can include inappopriate exon IDs and annotations, but not worth specializing the code
            y = ExonStructureData(gene, chr, strand, start, stop, exon_annot, ens_exonids, ens_transcripts)
            exon_info = [start,stop,y]
            new_entry = [exon_id,exon_info,(index1,index2)]
            #if key == ('ENSG00000129566', '14', '-'): print key,new_entry;kill
            exon_location[key].append(new_entry)"""
                
    exon_regions = {}
    for (gene,chr,strand) in exon_location:
        for exon_data in exon_location[(gene,chr,strand)]:
            try: exon_region_id,exon_detailed,(exon_num,region_num) = exon_data
            except ValueError: print exon_data;kill
            start,stop,ed = exon_detailed
            if strand == '+': start,stop,ed = exon_detailed
            else: stop,start,ed = exon_detailed
            y = ExonRegionData(gene, chr, strand, start, stop, ed.ExonID(), exon_region_id, exon_num, region_num, ed.Constitutive())
            try: exon_regions[gene].append(y)
            except KeyError: exon_regions[gene] = [y]
    """
    for key in exon_location:
        if key[0] == 'ENSG00000075413':
            print key
            for item in exon_location[key]: print item"""

    ###Create a corresponding database of intron and locations for the clustered block exon database
    intron_clusters={}; intron_region_db={}
    for key in exon_clusters:
        gene,chr,strand = key; chr = string.replace(chr,'chr','')
        exon_clusters[key].sort(); index = 0
        for exon_data in exon_clusters[key]:
            try: 
                exon_num,(start,stop),null,null = exon_data
                next_exon_data = exon_clusters[key][index+1]
                en,(st,sp),null,null = next_exon_data
                intron_num = exon_num
                if strand == '+': intron_start = stop+1; intron_stop = st-1; intron_start_stop = (intron_start,intron_stop)
                else: intron_start = start-1; intron_stop = sp+1; intron_start_stop = (intron_stop,intron_start)
                index+=1
                intron_data = intron_num,intron_start_stop,'','no'
                try: intron_clusters[key].append(intron_data)
                except KeyError: intron_clusters[key] = [intron_data]
                ###This database is used for SubGeneViewer and is analagous to region_db
                intron_region_id = 'I'+str(exon)+'-1'
                rd = ExonRegionData(gene, chr, strand, intron_start, intron_stop, ed.ExonID(), intron_region_id, intron_num, 1, 0)
                #rd = ExonRegionData(gene, chr, strand, intron_start, intron_stop, ed.ExonID(), intron_region_id, intron_num, 1, 0)
                if gene in intron_region_db:
                    block_db = intron_region_db[gene]
                    block_db[intron_num] = [rd]
                else:
                    block_db={}; block_db[intron_num] = [rd]
                    intron_region_db[gene] = block_db
            except IndexError: continue ### If no gene is added to intron_region_db, can be due to complete merger of all exons (multi-exon gene) when exon-exclusion occurs

    return exon_clusters,intron_clusters,exon_regions,intron_region_db

def eliminate_redundant_dict_values(database):
    for key in database:
        list = makeUnique(database[key])
        list.sort()
        database[key] = list
    return database

def getEnsemblAnnotations(filename,rna_processing_ensembl):
    fn=filepath(filename)
    ensembl_annotation_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if 'Description' in data: ensembl_gene_id,description,symbol = string.split(data,'\t')
        else: ensembl_gene_id,symbol,description = string.split(data,'\t')
        ensembl_annotation_db[ensembl_gene_id] = description, symbol

    for gene in ensembl_annotation_db:
        if gene in rna_processing_ensembl: mRNA_processing = 'RNA_processing/binding'
        else: mRNA_processing = ''
        index = ensembl_annotation_db[gene] 
        ensembl_annotation_db[gene] = index[0],index[1],mRNA_processing
        
    exportEnsemblAnnotations(ensembl_annotation_db)
    return ensembl_annotation_db

def exportEnsemblAnnotations(ensembl_annotation_db):
    exon_annotation_export = 'AltDatabase/ensembl/' +species+'/'+species+ '_Ensembl-annotations.txt'
    fn=filepath(exon_annotation_export); data = open(fn,'w')
    for ensembl_gene in ensembl_annotation_db:
        a = ensembl_annotation_db[ensembl_gene]
        values = ensembl_gene +'\t'+ a[0] +'\t'+ a[1] +'\t'+ a[2] + '\n'
        data.write(values)
    data.close()

def reimportEnsemblAnnotations(species,symbolKey=False):
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations.txt'
    fn=filepath(filename)
    ensembl_annotation_db = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        ensembl_gene_id,symbol,description,mRNA_processing = string.split(data,'\t')
        if symbolKey:
            ensembl_annotation_db[symbol] = ensembl_gene_id
        else:
            ensembl_annotation_db[ensembl_gene_id] = symbol,description,mRNA_processing
    return ensembl_annotation_db
    
def importPreviousBuildTranslation(filename,use_class_structures):
    ###When previous genome build coordinates a provided for an array - create a translation
    ###file between ensembl exon ids
    fn=filepath(filename)
    exon_coord_translation_db = {}; gene_coor_translation_db = {}; exon_db = {}
    gene_coor_translation_list = []; x=0
    for line in open(fn,'r').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x+=1 ###ignore the first line
        else:
            chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, null, null, null,null, constitutive_exon = string.split(data,'\t')
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            gene_start = int(gene_start); gene_stop = int(gene_stop)
            exon_start = int(exon_start); exon_stop = int(exon_stop)
            if strand == '-1': strand = '-'
            if strand == '1': strand = '+'
            if strand == '-': new_gene_start = gene_stop; new_gene_stop = gene_start; new_exon_start = exon_stop; new_exon_stop = exon_start 
            else: new_gene_start = gene_start; new_gene_stop = gene_stop; new_exon_start = exon_start; new_exon_stop = exon_stop 
            new_exon_start = abs(new_gene_start - new_exon_start)
            new_exon_stop = abs(new_gene_start - new_exon_stop)
            ###constitutive_exon column contains a 0 or 1: ensembl_exon_id is versioned
            ensembl_exon_id,null = string.split(ensembl_exon_id,'.')

            if use_class_structures == 'yes':
                exon_coord_info = ensembl_gene_id,(exon_start,exon_stop),(gene_start,gene_stop),int(constitutive_exon)
                ei = EnsemblInformation(chr, gene_start, gene_stop, strand, ensembl_gene_id, ensembl_exon_id, exon_start, exon_stop, constitutive_exon, new_exon_start, new_exon_stop, new_gene_start, new_gene_stop)    
                ###Also independently determine exon clusters for previous build exon data
                exon_annot=''; exon_info = (exon_start,exon_stop,ensembl_exon_id,exon_annot)
                try: exon_db[(ensembl_gene_id,chr,strand)].append(exon_info)
                except KeyError: exon_db[(ensembl_gene_id,chr,strand)] = [exon_info]
                y = ei
            else:
                ei = [chr,(gene_start,gene_stop),ensembl_gene_id,(new_gene_start,new_gene_stop)]
                y = [chr,(exon_start,exon_stop),ensembl_exon_id,(new_exon_start,new_exon_stop)]
            try: exon_coord_translation_db[ensembl_gene_id].append(y)
            except KeyError: exon_coord_translation_db[ensembl_gene_id] = [y]
            gene_coor_translation_db[(chr,strand),gene_start,gene_stop] = ei
    for key in gene_coor_translation_db:
        ei = gene_coor_translation_db[key]
        gene_coor_translation_list.append([key,ei])
    gene_coor_translation_list.sort()
    gene_coor_translation_list2={}
    for key in gene_coor_translation_list:
        chr_strand = key[0][0]
        try: gene_coor_translation_list2[chr_strand].append(key[-1])
        except KeyError: gene_coor_translation_list2[chr_strand] = [key[-1]]
    gene_coor_translation_list = gene_coor_translation_list2
    ###Determine Exon Clusters for current Ensembl genes with poor linkage properties (can't be converted to old coordiantes) 
    exon_db2 = annotate_exons(exon_db)
    exon_clusters = exon_clustering(exon_db2); exon_db2={}
    return exon_coord_translation_db, exon_db, exon_clusters

def importExonTranscriptAnnotations(filename):
    fn=filepath(filename)
    exon_trans_association_db = {}; x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 0: x+=1
        else:
            ensembl_gene_id,ensembl_trans_id,ensembl_exon_id,constitutive = string.split(data,'\t')
            ensembl_exon_id,null = string.split(ensembl_exon_id,'.')
            try: exon_trans_association_db[ensembl_exon_id].append([ensembl_trans_id,constitutive])
            except KeyError: exon_trans_association_db[ensembl_exon_id] = [[ensembl_trans_id,constitutive]]
    return exon_trans_association_db 

def importEnsemblDomainData(filename):
    fn=filepath(filename); x = 0; ensembl_ft_db = {}; ensembl_ft_summary_db = {} # Use the last database for summary statistics
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x == 1:
            try: ensembl_gene, chr, mgi, uniprot, ensembl_prot, seq_data, position_info = string.split(data,'\t')
            except ValueError: continue
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            ft_info_list = string.split(position_info,' | ')
            for entry in ft_info_list:
                try: peptide_start_end, gene_start_end, feature_source, interprot_id, description = string.split(entry,' ')
                except ValueError: continue
                ###142-180 3015022-3015156 Pfam IPR002050 Env_polyprotein
                ft_start_pos, ft_end_pos = string.split(peptide_start_end,'-')
                pos1 = int(ft_start_pos); pos2 = int(ft_end_pos)
                
                #sequence_fragment = seq_data[pos1:pos2]
                if len(description)>1 or len(interprot_id)>1:
                    ft_info = [description,sequence_fragment,interprot_id]
                    ft_info2 = description,interprot_id
                    ###uniprot_ft_db[id].append([ft,pos1,pos2,annotation])
                    try: ensembl_ft_db[id].append(ft_info)
                    except KeyError: ensembl_ft_db[id] = [ft_info]
                    try: ensembl_ft_summary_db[id].append(ft_info2)
                    except KeyError: ensembl_ft_summary_db[id] = [ft_info2]         
        elif data[0:6] == 'GeneID': x = 1
        
    ensembl_ft_db = eliminate_redundant_dict_values(ensembl_ft_db)
    ensembl_ft_summary_db = eliminate_redundant_dict_values(ensembl_ft_summary_db)
    domain_gene_counts = {}
    ###Count the number of domains present in all genes (count a domain only once per gene)
    for gene in ensembl_ft_summary_db:
        for domain_info in ensembl_ft_summary_db[gene]:
            try: domain_gene_counts[domain_info] += 1
            except KeyError: domain_gene_counts[domain_info] = 1
    print "Number of Ensembl genes, linked to array genes with domain annotations:",len(ensembl_ft_db)
    print "Number of Ensembl domains:",len(domain_gene_counts)
    return ensembl_ft_db,domain_gene_counts
       
def getEnsemblAssociations(Species,data_type,test_status):
    global species; species = Species
    global test; test = test_status
    global test_gene
    meta_test = ["ENSG00000224972","ENSG00000107077"]
    test_gene = ['ENSG00000163132'] #'ENSG00000215305
    #test_gene = ['ENSMUSG00000000037'] #'test Mouse - ENSMUSG00000000037
    test_gene = ['ENSMUSG00000065005'] #'ENSG00000215305
    test_gene = ['ENSRNOE00000194194']
    test_gene = ['ENSG00000229611','ENSG00000107077','ENSG00000107077','ENSG00000107077','ENSG00000107077','ENSG00000107077','ENSG00000163132', 'ENSG00000115295']
    test_gene = ['ENSG00000110955']
    #test_gene = ['ENSMUSG00000059857'] ### for JunctionArrayEnsemblRules
    #test_gene = meta_test
    exon_annotation_db,transcript_gene_db,gene_transcript,intron_retention_db,ucsc_splicing_annot_db = getEnsExonStructureData(species,data_type)
    exon_annotation_db2 = annotate_exons(exon_annotation_db); ensembl_descriptions={}
    exon_db = customDBDeepCopy(exon_annotation_db2) ##having problems with re-writting contents of this db when I don't want to
    exon_clusters,intron_clusters,exon_regions,intron_region_db = exon_clustering(exon_db); exon_db={}
    exon_junction_db,putative_as_junction_db,exon_junction_db,full_junction_db,excluded_intronic_junctions = processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,intron_retention_db)
    exon_regions,critical_gene_junction_db = compareJunctions(species,putative_as_junction_db,exon_regions)

    exportSubGeneViewerData(exon_regions,exon_annotation_db2,critical_gene_junction_db,intron_region_db,intron_retention_db,full_junction_db,excluded_intronic_junctions,ucsc_splicing_annot_db)
    full_junction_db=[]
    
    ###Grab rna_processing Ensembl associations
    use_exon_data='no';get_splicing_factors = 'yes'
    try: rna_processing_ensembl = GO_parsing.parseAffyGO(use_exon_data,get_splicing_factors,species)
    except Exception: rna_processing_ensembl={}
    
    ensembl_annot_file = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = getEnsemblAnnotations(ensembl_annot_file,rna_processing_ensembl)
    #exportExonClusters(exon_clusters,species)

    return exon_annotation_db2,ensembl_annotation_db,exon_clusters,intron_clusters,exon_regions,intron_retention_db,ucsc_splicing_annot_db,transcript_gene_db
  
def getExonTranscriptDomainAssociations(Species):
    global species; species = Species
    import_dir = '/AltDatabase/ensembl/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for file_name in dir_list:    #loop through each file in the directory to output results
        dir_file = 'AltDatabase/ensembl/'+species+'/'+file_name
        if 'Exon_cDNA' in dir_file: exon_trans_file = dir_file
        elif 'Domain' in dir_file: domain_file = dir_file
    exon_trans_association_db = importExonTranscriptAnnotations(exon_trans_file)
    return exon_trans_association_db
    
def exportExonClusters(exon_clusters,species):
    exon_cluster_export = 'AltDatabase/ensembl/'+species+'/'+species+'-Ensembl-exon-clusters.txt'
    fn=filepath(exon_cluster_export); data = open(fn,'w')
    for key in exon_clusters:
        ensembl_gene = key[0]
        chr = key[1]
        strand = key[2]
        for entry in exon_clusters[key]:
            exon_cluster_num = str(entry[0])
            exon_start = str(entry[1][0])
            exon_stop = str(entry[1][1])
            exon_id_list = string.join(entry[2],'|')
            annotation_list = string.join(entry[3],'|')
            values = ensembl_gene +'\t'+ chr +'\t'+ strand +'\t'+ exon_cluster_num +'\t'+ exon_start +'\t'+ exon_stop +'\t'+ exon_id_list +'\t'+ annotation_list +'\n'
            data.write(values)
    data.close()

def checkforEnsemblExons(trans_exon_data):
    proceed_status = 0
    for (start_pos,ste,spe) in trans_exon_data:
        #print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
        if ste.ExonID() == '' or spe.ExonID() == '': print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
        if ('-' in ste.ExonID()) and ('-' in spe.ExonID()): null=[]
        else: proceed_status +=1
        #else: print ste.ExonRegionNumbers(),spe.ExonRegionNumbers(), [ste.ExonID(),spe.ExonID()];kill
    if proceed_status>0: proceed_status = 'yes'
    else: proceed_status = 'no'
    return proceed_status

def getEnsemblGeneLocations(species,array_type,key):
    if key == 'key_by_chromosome':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    else:
        if array_type == 'RNASeq':
            filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
        else:
            filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_Ensembl_probesets.txt'
    fn=filepath(filename); x=0; gene_strand_db={}; gene_location_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            try: gene, chr, strand, exon_start, exon_end, ens_exonid, constitutive_exon, ens_transcriptid = t
            except Exception:
                if array_type == 'RNASeq': gene = t[0]; chr=t[2]; strand=t[3]; exon_start=t[4]; exon_end=t[5] ### for Ensembl_exon.txt
                else: gene = t[2]; chr=t[4]; strand=t[5]; exon_start=t[6]; exon_end=t[7] ### for Ensembl_probesets.txt
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if strand == '-1' or strand == '-': strand = '-'
            else: strand = '+'
            exon_end = int(exon_end); exon_start = int(exon_start)
            gene_strand_db[gene] = strand, chr
            try: gene_location_db[gene]+=[exon_start,exon_end]
            except Exception: gene_location_db[gene]=[exon_start,exon_end]
    if key == 'key_by_chromosome':
        location_gene_db={}; chr_gene_db={}
        for gene in gene_location_db:
            gene_location_db[gene].sort()
            start = gene_location_db[gene][0]
            end = gene_location_db[gene][-1]
            strand,chr = gene_strand_db[gene]
            location_gene_db[chr,start,end]=gene,strand
            try: chr_gene_db[chr].append([start,end])
            except Exception: chr_gene_db[chr]=[[start,end]]
        return chr_gene_db,location_gene_db
    else:
        for gene in gene_location_db:
            gene_location_db[gene].sort()
            start = gene_location_db[gene][0]
            end = gene_location_db[gene][-1]
            strand,chr = gene_strand_db[gene]
            gene_location_db[gene]=chr,strand,str(start),str(end)
        return gene_location_db
    
def getAllEnsemblUCSCTranscriptExons():
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_transcript-annotations.txt'
    importEnsExonStructureData(filename,species,'exon-transcript')
    filename = 'AltDatabase/ucsc/'+species+'/'+species+'_UCSC_transcript_structure_mrna.txt'
    importEnsExonStructureData(filename,species,'exon-transcript')
    exon_transcripts = eliminate_redundant_dict_values(exon_transcript_db)
    return exon_transcripts

def processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,intron_retention_db):
    ###Parses transcript to exon relationships and links these to previously described distinct exon regions.
    ###The exon blocks and region numbers provide a simple semantic for determining where in the transcript and what event is occuring
    ###when directly comparing junctions to each other
    exon_transcripts = getAllEnsemblUCSCTranscriptExons()
    transcript_exon_db={}; constitutive_region_gene_db={}; null=[]; k=[]
    for (gene,chr,strand) in exon_annotation_db:
        constitutive_region_count_db={}; constitutive_region_count_list=[]; count_db={}
        for exon_info in exon_annotation_db[(gene,chr,strand)]:
            if strand == '+': exon_start,exon_end,y = exon_info ###Link the original exon-start/stop data back to the region based data
            else: exon_end,exon_start,y = exon_info ###This shouldn't be necessary, but the coordinates get reversed in annotate_exons and copying the original via deepcopy takes way too much time
            transcripts = exon_transcripts[y.ExonID()] ### ERROR WILL OCCUR HERE IF YOU DON'T REBUILD THE UCSC DATABASE!!! (different version still left over from Protein analysis methods)
            if gene in exon_regions:
                ste=''; spe=''
                #print exon_start, y.ExonID(); y.TranscriptID(); y.IntronDeletionStatus()
                for ed in exon_regions[gene]:
                    #print ed.ExonRegionID(),transcripts, [ed.ExonStart(), ed.ExonStop(), exon_start,exon_end]
                    
                    ### Since we know which exon region each Ensembl/UCSC exon begins and ends, we can just count the number of transcripts
                    ### that contain each of the Ensembl/UCSC exon ID, which will give us our constitutive exon count for each gene
                    proceed = 'no'
                    if ((ed.ExonStart() >= exon_start) and (ed.ExonStart() <= exon_end)) and ((ed.ExonStop() >= exon_start) and (ed.ExonStop() <= exon_end)): proceed = 'yes'
                    elif ((ed.ExonStart() <= exon_start) and (ed.ExonStart() >= exon_end)) and ((ed.ExonStop() <= exon_start) and (ed.ExonStop() >= exon_end)): proceed = 'yes'
                    if proceed == 'yes': ### Ensures that each examined exon region lies directly within the examined exon (begining and end) -> version 2.0.6
                        try: constitutive_region_count_db[ed.ExonRegionID()]+=list(transcripts) #rd.ExonNumber() rd.ExonRegionID()
                        except KeyError: constitutive_region_count_db[ed.ExonRegionID()]=list(transcripts) ### If you don't convert to list, the existing list object will be updated non-specifically
                    ### Identify the exon regions that associate with the begining and end-positions of each exon to later determine the junction in exon region space
                    if exon_start == ed.ExonStart(): ste = ed  ### Two exon regions can not have the same start position in this database
                    if exon_end == ed.ExonStop(): spe = ed ### Two exon regions can not have the same end position in this database
                    if spe!='' and ste!='': break 
                    #print ed.ExonStart(),ed.ExonStop()
                    
                if spe!='' and ste!='':
                    values = exon_start,ste,spe #,y
                    ste.reSetExonID(y.ExonID()); spe.reSetExonID(y.ExonID())  ###Need to represent the original source ExonID to eliminate UCSC transcripts without Ensembl exons
                    for ens_transcriptid in y.TranscriptID():
                        #print exon_start, ste.ExonID(), spe.ExonID(),y.TranscriptID()
                        try: transcript_exon_db[ens_transcriptid].append(values)
                        except KeyError: transcript_exon_db[ens_transcriptid] = [values]
                else:
                    ### Indicates this exon has been classified as a retained intron
                    ### Keep it so we know where this exon is and to prevent inclusion of faulty flanking junctions
                    ### This is needed so we can retain junctions that are unique to this transcript but not next to the retained intron
                    values = y.ExonStart(),y,y
                    for ens_transcriptid in y.TranscriptID():
                        try: transcript_exon_db[ens_transcriptid].append(values)
                        except KeyError: transcript_exon_db[ens_transcriptid] = [values]    
            else: k.append(gene)
        ### Get the number of transcripts associated with each region
        for exon_region in constitutive_region_count_db:
            #print exon_region, unique.unique(constitutive_region_count_db[exon_region])
            transcript_count = len(unique.unique(constitutive_region_count_db[exon_region]))
            constitutive_region_count_list.append(transcript_count)
            try: count_db[transcript_count].append(exon_region)
            except KeyError: count_db[transcript_count] = [exon_region]
        constitutive_region_count_list = unique.unique(constitutive_region_count_list)
        constitutive_region_count_list.sort()
        try: max_count = constitutive_region_count_list[-1]
        except Exception:
            print gene, constitutive_region_count_list, constitutive_region_count_db, count_db;kill
        #print count_db
        cs_exon_region_ids = list(count_db[max_count])
        ### If there is only one constitutive region and one set of strong runner ups, include these to improve the constitutive expression prediction
        if (len(cs_exon_region_ids)==1 and len(constitutive_region_count_list)>2) or len(constitutive_region_count_list)>3: ### Decided to add another heuristic that will add the second highest scoring exons always (if at least 4 frequencies present)
            second_highest_count = constitutive_region_count_list[-2]
            if (max_count-second_highest_count)<3: ### Indicates that the highest and second highest common exons be similiar in terms of their frequency (different by two transcripts at most)
                #print cs_exon_region_ids
                #print second_highest_count
                if second_highest_count != 1: ### Don't inlcude if there is only one transcript contributing
                    cs_exon_region_ids += list(count_db[second_highest_count])
        constitutive_region_gene_db[gene] = cs_exon_region_ids
        #print gene, cs_exon_region_ids;sys.exit()

    exon_transcripts=[]; del exon_transcripts
    ### Reset the constitutive exon, previously assigned by Ensembl - ours should be more informative since it uses more transcript data and specific region info
    #"""
    for gene in constitutive_region_gene_db:
        if gene in exon_regions:
            for ed in exon_regions[gene]:
                if ed.ExonRegionID() in constitutive_region_gene_db[gene]: ed.setConstitutive('1')
                else: ed.setConstitutive('0')
                #print ed.ExonRegionID(), ed.Constitutive()
    #"""
    null = unique.unique(null); #print len(null)
    print len(k), 'genes, not matched up to region database'
    print len(transcript_gene_db), 'transcripts from Ensembl being analyzed'

    transcript_exon_db = eliminate_redundant_dict_values(transcript_exon_db)
    gene_transcript = eliminate_redundant_dict_values(gene_transcript)
    
    gene_transcript_multiple={}; tc=0
    for gene in gene_transcript:
        if len(gene_transcript[gene])>1: tc+=1; gene_transcript_multiple[gene]=len(gene_transcript[gene])
    print tc,"genes with multiple transcripts associated from Ensembl"

    """
    ###If a retained intron is present we must ignore all exons downstream of where we DELETED that exon information (otherwise there is false junciton information introduced)    
    ###Here we simply delete the information for that transcript all together
    td=0
    for key in intron_retention_db:
        transcripts_to_delete = {}
        for (s1,s2,y) in intron_retention_db[key]:
            if y.IntronDeletionStatus() == 'yes': ###only do this for exon's deleted upon import
                for transcript in y.TranscriptID(): transcripts_to_delete[transcript] = []
        for transcript in transcripts_to_delete:
            ###may not be present if the exon deleted constituted the whole transcript
            if transcript in transcript_exon_db:
                for (start_pos,ste,spe) in transcript_exon_db[transcript]:
                    print ste.ExonRegionID(),spe.ExonRegionID()
                del transcript_exon_db[transcript]; td+=1
    print td, "transcripts deleted with intron retention. Required for down-stream junction analysis"
    """

    exon_junction_db={}; full_junction_db={}; excluded_intronic_junctions={}; rt=0
    mx_detection_db={}; transcript_exon_region_db={}; #junction_transcript_db={}
    ###Sort and filter the junction data
    for transcript in transcript_exon_db:
        gene,chr,strand = transcript_gene_db[transcript]
        transcript_exon_db[transcript].sort()
        if strand == '-': transcript_exon_db[transcript].reverse()
        index=0
        #if 'BC' in transcript: print transcript, transcript_exon_db[transcript]
        ###Introduced a filter to remove transcripts from UCSC with no supporting Ensembl exons (distinct unknown transcript type)
        ###Since these are already exon regions, we exclude them from alt. splicing/promoter assignment.
        proceed_status = checkforEnsemblExons(transcript_exon_db[transcript])
        ###Loop through the exons in each transcript
        if proceed_status == 'yes':
            #print transcript
            #print transcript_exon_db[transcript]
            for (start_pos,ste,spe) in transcript_exon_db[transcript]:
                if (index+1) != len(transcript_exon_db[transcript]): ###Don't grab past the last exon in the transcript
                    start_pos2,ste2,spe2 = transcript_exon_db[transcript][index+1]
                    #print spe.ExonID(),ste.ExonID(),spe2.ExonID(),ste2.ExonID(), spe.ExonRegionID2(),ste.ExonRegionID2(),spe2.ExonRegionID2(),ste2.ExonRegionID2(),ste.IntronDeletionStatus(),ste2.IntronDeletionStatus(),spe.ExonStop(),ste2.ExonStart()
                    #print transcript,spe.ExonStop(),ste2.ExonStart(), ste.IntronDeletionStatus(),ste2.IntronDeletionStatus()
                    if ste.IntronDeletionStatus() == 'no' and ste2.IntronDeletionStatus() == 'no':
                        ### Don't include junctions where the current or next junction was a removed retained intron (but keep other junctions in the transcript)
                        exon_junction = (ste.ExonRegionNumbers(),spe.ExonRegionNumbers()),(ste2.ExonRegionNumbers(),spe2.ExonRegionNumbers())
                        try: exon_junction_db[gene].append(exon_junction)
                        except KeyError: exon_junction_db[gene] = [exon_junction]
                        #try: junction_transcript_db[gene,exon_junction].append(transcript)
                        #except KeyError: junction_transcript_db[gene,exon_junction] = [transcript]
                        #try: transcript_exon_region_db[transcript]+=[ste.ExonRegionNumbers(),spe.ExonRegionNumbers(),ste2.ExonRegionNumbers(),spe2.ExonRegionNumbers()]
                        #except Exception: transcript_exon_region_db[transcript] = [ste.ExonRegionNumbers(),spe.ExonRegionNumbers(),ste2.ExonRegionNumbers(),spe2.ExonRegionNumbers()]
                        try: full_junction_db[gene].append((spe.ExonRegionID2(),ste2.ExonRegionID2()))
                        except KeyError: full_junction_db[gene] = [(spe.ExonRegionID2(),ste2.ExonRegionID2())]
                        ### Look for potential mutually-exclusive splicing events
                        if (index+2) != len(transcript_exon_db[transcript]):
                            start_pos3,ste3,spe3 = transcript_exon_db[transcript][index+2]
                            if ste3.IntronDeletionStatus() == 'no':
                                try: mx_detection_db[gene,spe.ExonRegionID2(),ste3.ExonRegionID2()].append(spe2)
                                except KeyError: mx_detection_db[gene,spe.ExonRegionID2(),ste3.ExonRegionID2()]=[spe2]
                    else:
                        ### Retain this information when exporting all known junctions
                        #print transcript,spe.ExonStop(),ste2.ExonStart()
                        try: excluded_intronic_junctions[gene].append((spe,ste2))
                        except KeyError: excluded_intronic_junctions[gene]=[(spe,ste2)]

                index+=1
        else: rt +=1

    mx_exons=[]
    for key in mx_detection_db:
        gene = key[0]
        if len(mx_detection_db[key])>1:
            cassette_block_ids=[]; cassette_block_db={}
            for spe2 in mx_detection_db[key]:
                cassette_block_ids.append(spe2.ExonNumber())
                cassette_block_db[spe2.ExonNumber()]=spe2
            cassette_block_ids = unique.unique(cassette_block_ids)
            if len(cassette_block_ids)>1:
                for exon in cassette_block_ids:
                    spe2 = cassette_block_db[exon]
                    spe2.setMutuallyExclusive()
                    #print key, spe2.ExonRegionID(), spe2.ExonID()
                    #try: mx_exons[gene].append(spe2)
                    #except KeyError: mx_exons[gene]=[sp2]
                    
    print rt, "transcripts removed from analysis with no Ensembl exon evidence. Results in more informative splicing annotations downstream"
    #print len(junction_transcript_db), 'length of junction_transcript_db'
    print len(exon_junction_db),'length of exon_junction_db'

    full_junction_db = eliminate_redundant_dict_values(full_junction_db)
    
    ###Stringent count, since it requires all region information for each exon and common splice events occur for just one region to another     
    ###example: (((8, 1), (8, 1)), ((9, 1), (9, 1))), (((8, 1), (8, 1)), ((9, 1), (9, 2)))
    putative_as_junction_db={}
    for gene in exon_junction_db:
        junctions = exon_junction_db[gene]
        junction_count={}
        for junction in exon_junction_db[gene]:
            try: junction_count[junction]+=1
            except KeyError: junction_count[junction]=1
        count_junction={}; count_list=[]
        for junction in junction_count:
            count = junction_count[junction]
            try: count_junction[count].append(junction)
            except KeyError: count_junction[count] = [junction]
            count_list.append(count)
        count_list = unique.unique(count_list); count_list.sort() ###number of unique counts
        if len(count_list)>1 and gene in gene_transcript_multiple: ###Otherwise, there is no variation in junction number between transcripts
            transcript_number = gene_transcript_multiple[gene]
            max_count = count_list[-1] ###most common junction - not alternatively spliced
            if max_count == transcript_number: ###Ensures we are grabbing everything but constitutive exons (max_count can include AS if no Ensembl constitutive).
                for count in count_junction:
                    if count != max_count:
                        junctions = count_junction[count]
                        junctions.sort()
                        try: putative_as_junction_db[gene]+=junctions
                        except KeyError: putative_as_junction_db[gene]=junctions
            else:
                try: putative_as_junction_db[gene]+=junctions
                except KeyError: putative_as_junction_db[gene]=junctions
        elif gene in gene_transcript_multiple:
            ###If there are multiple transcripts, descriminating these is difficult, just include all junctions for that gene
            try: putative_as_junction_db[gene]+=junctions
            except KeyError: putative_as_junction_db[gene]=junctions
    return exon_junction_db,putative_as_junction_db,exon_junction_db,full_junction_db,excluded_intronic_junctions

def reformatJunctions(exons,type):
    exons2=[]
    for (b,i) in exons:
        exons2.append('E'+str(b)+'.'+str(i))
    if type == 'junction': exons2 = string.join(exons2,'-')
    else: exons2 = string.join(exons2,'|')
    return exons2
    
def compareJunctions(species,putative_as_junction_db,exon_regions,rootdir=None,searchChr=None):
    ### Add polyA site information and mutually-exclusive splicing site
    if len(exon_regions)==0:
        export_annotation = '_de-novo'
        alt_junction_export = rootdir+'/AltDatabase/ensembl/'+species+'/denovo/'+species+'_alternative_junctions'+export_annotation+'.'+searchChr+'.txt'
        import export
        data = export.ExportFile(alt_junction_export)
    else:
        export_annotation = ''
        alt_junction_export = 'AltDatabase/ensembl/'+species+'/'+species+'_alternative_junctions'+export_annotation+'.txt'
    if export_annotation != '_de-novo':
        print 'Writing the file:',alt_junction_export
        print len(putative_as_junction_db),'genes being examined for AS/alt-promoters in Ensembl'
        fn=filepath(alt_junction_export); data = open(fn,'w')
    title = ['gene','critical-exon-id','junction1','junction2']
    title = string.join(title,'\t')+'\n'; data.write(title)
    
    ###Find splice events based on structure based evidence
    
    critical_exon_db={}; j=0; global add_to_for_terminal_exons; add_to_for_terminal_exons={}
    complex3prime_event_db={}; complex5prime_event_db={}; cassette_exon_record={}
    for gene in putative_as_junction_db:
        #if gene == 'ENSMUSG00000000028': print putative_as_junction_db[gene];kill
        for j1 in putative_as_junction_db[gene]:
            for j2 in putative_as_junction_db[gene]: ### O^n squared query
                if j1 != j2:
                    temp_junctions = [j1,j2]; temp_junctions.sort(); junction1,junction2 = temp_junctions
                    splice_junctions=[]; critical_exon=[]; splice_type=''
                    e1a,e2a = junction1; e1b,e2b = junction2 ###((8, 2), (8, 2)) break down the exon into exon_block,region tubles.
                    e1a3,e1a5 = e1a; e2a3,e2a5 = e2a; e1b3,e1b5 = e1b; e2b3,e2b5 = e2b ###(8, 2) break down the exons into single tuples (designating 5' and 3' ends of the exon): 
                    e1a3_block,e1a3_reg = e1a3; e1a5_block,e1a5_reg = e1a5; e2a3_block,e2a3_reg = e2a3; e2a5_block,e2a5_reg = e1a5
                    e1b3_block,e1b3_reg = e1b3; e1b5_block,e1b5_reg = e1b5 ;e2b3_block,e2b3_reg = e2b3; e2b5_block,e2b5_reg = e1b5
                    splice_junctions = [(e1a5,e2a3),(e1b5,e2b3)] ###three junctions make up the cassette event, record the two evidenced by this comparison and agglomerate after all comps
                    splice_junction_str = reformatJunctions(splice_junctions[0],'junction')+'\t'+reformatJunctions(splice_junctions[1],'junction')
                    ###IMPORTANT NOTE: The temp_junctions are sorted, but doesn't mean that splice_junctions is sorted correctly... must account for this
                    splice_junctions2 = list(splice_junctions); splice_junctions2.sort()
                    if splice_junctions2 != splice_junctions: ###Then the sorting is wrong and the down-stream method won't work
                        ###Must re-do the above assingments
                        junction2,junction1 = temp_junctions
                        e1a,e2a = junction1; e1b,e2b = junction2 ###((8, 2), (8, 2)) break down the exon into exon_block,region tubles.
                        e1a3,e1a5 = e1a; e2a3,e2a5 = e2a; e1b3,e1b5 = e1b; e2b3,e2b5 = e2b ###(8, 2) break down the exons into single tuples (designating 5' and 3' ends of the exon): 
                        e1a3_block,e1a3_reg = e1a3; e1a5_block,e1a5_reg = e1a5; e2a3_block,e2a3_reg = e2a3; e2a5_block,e2a5_reg = e1a5
                        e1b3_block,e1b3_reg = e1b3; e1b5_block,e1b5_reg = e1b5 ;e2b3_block,e2b3_reg = e2b3; e2b5_block,e2b5_reg = e1b5
                        splice_junctions = [(e1a5,e2a3),(e1b5,e2b3)]
                        splice_junction_str = reformatJunctions(splice_junctions[0],'junction')+'\t'+reformatJunctions(splice_junctions[1],'junction')
                    if e1a5_block == e2a3_block or e1b5_block == e2b3_block: continue ###suggests splicing within a block... we won't deal with these
                    if e1a5 == e1b5:  ###If 5'exons in the junction are the same
                        ###make sure the difference isn't in the 5' side of the next splice junction (or exon end)
                        if e2a3 != e2b3: #(((5, 1), (5, 1)*), ((6, 1)*, (6, 1)))  --  (((5, 1), (5, 1)*), ((7, 1)*, (7, 1)))
                            if e2a3_block == e2b3_block: #[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 1)*), ((2, 3)*, (2, 3)))]
                                splice_type = "alt-3'"; critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')   
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                            else:
                                critical_exon = [e2a3]; splice_type = 'cassette-exon' #[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 1)*), ((3, 1)*, (3, 1)))]    
                                try: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                                except KeyError: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]
                                try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                                except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                #print critical_exon,splice_type,splice_junction_str
                    if splice_type =='' and e2a3 == e2b3: ###If 3'exons in the junction are the same
                        if e1a5 != e1b5:
                            if e1a5_block == e1b5_block: ###simple alt 5' splice site encountered
                                splice_type = "alt-5'"; critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)#[(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 1), (1, 3)*), ((2, 1)*, (2, 1)))]
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                            else:
                                splice_type = 'cassette-exon'; critical_exon = [e1b5] #[(((1, 1), (1, 1)*), ((3, 1)*, (3, 1))) ----(((2, 1), (2, 1)*), ((3, 1)*, (3, 1)))]       
                                try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                                except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                                try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                                except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]                                
                                #if gene == 'ENSG00000128606' and critical_exon == [(4, 1)] : print junction1,junction2;kill
                                y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                #print critical_exon,splice_type,splice_junction_str
                    if splice_type =='' and e2a3_block == e2b3_block and e1a5_block != e2a3_block and e1b5_block != e2b3_block: ###Begin looking at complex examples: If 3'exon blocks in the junction are the same
                        if e1a5_block == e1b5_block: #alt5'-alt3' [(((1, 1), (1, 1)*), ((2, 1)*, (2, 1))) ----(((1, 3), (1, 3)*), ((2, 3)*, (2, 3)))]
                            critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon); splice_type = "alt5'-alt3'"
                            if len(critical_exon)>0:
                                alt5_exon = [critical_exon[0]]; alt3_exon = [critical_exon[1]]
                                #print alt5_exon,critical_exon,critical_exon;kill
                                y = CriticalExonInfo(gene,alt5_exon,"alt-5'",splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y]
                                y = CriticalExonInfo(gene,alt3_exon,"alt-3'",splice_junctions)
                                data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                                try: critical_exon_db[gene].append(y)
                                except KeyError: critical_exon_db[gene] = [y] 
                        else: #cassette-alt3' [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((2, 1), (2, 3)*), ((4, 3)*, (4, 3)))]
                            critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                            splice_type = "alt-3'"
                            y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                            data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                            try: critical_exon_db[gene].append(y)
                            except KeyError: critical_exon_db[gene] = [y]
                            if e1a5_block < e1b5_block:
                                critical_exon = [e1b5]
                                try: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                                except KeyError: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                                complex3prime_event_db[gene,e1b5_block] = e1a5
                                try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                                except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]
                            elif e1a5_block != e1b5_block:
                                critical_exon = [e1a5]
                                try: add_to_for_terminal_exons[gene,e1a5_block] = [e1b5]
                                except KeyError: add_to_for_terminal_exons[gene,e1a5_block].append(e1b5)
                                complex3prime_event_db[gene,e1a5_block] = e1b5
                                try: cassette_exon_record[gene,e1a5_block].append(e2a3_block)
                                except KeyError: cassette_exon_record[gene,e1a5_block] = [e2a3_block]
                            splice_type = "cassette-exon"
                            y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                            data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                            try: critical_exon_db[gene].append(y)
                            except KeyError: critical_exon_db[gene] = [y] 
                    if splice_type =='' and e1a5_block == e1b5_block and e1a5_block != e2a3_block and e1b5_block != e2b3_block:
                        #alt5'-cassette' [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((1, 1), (1, 3)*), ((5, 1)*, (5, 1)))]
                        critical_exon = pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon)
                        splice_type = "alt-5'"
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                            
                        if e2a3_block < e2b3_block:
                            critical_exon = [e2a3]
                            try: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]
                            except KeyError: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                            complex5prime_event_db[gene,e2a3_block] = e2b3
                            try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                            except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                        elif e2a3_block != e2b3_block:
                            critical_exon = [e2b3]
                            try: add_to_for_terminal_exons[gene,e2b3_block] = [e2a3]
                            except KeyError: add_to_for_terminal_exons[gene,e2b3_block].append(e2a3)
                            complex5prime_event_db[gene,e2b3_block] = e2a3
                            try: cassette_exon_record[gene,e2b3_block].append(e1b5_block)
                            except KeyError: cassette_exon_record[gene,e2b3_block] = [e1b5_block]
                        splice_type = "cassette-exon"
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                    if splice_type =='' and e1a5_block<e1b5_block and e2b3_block>e2a3_block and e2a3_block>e1b5_block:
                        #mx-mix [(((1, 1), (1, 1)*), ((4, 1)*, (4, 1))) ----(((2, 1), (2, 1)*), ((5, 1)*, (5, 1)))]
                        critical_exon = [e2a3,e1b5]#; mx_event_db[gene,e2a3] = e1b5; mx_event_db[gene,e1b5] = e2a3
                        splice_type = 'cassette-exon'
                        #"""
                        try: add_to_for_terminal_exons[gene,e2a3_block].append(e2b3)
                        except KeyError: add_to_for_terminal_exons[gene,e2a3_block] = [e2b3]                        
                        try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                        except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                        #"""       
                        try: cassette_exon_record[gene,e2a3_block].append(e1a5_block)
                        except KeyError: cassette_exon_record[gene,e2a3_block] = [e1a5_block]
                        try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                        except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]   
                        #print 'mx-mx',critical_exon, splice_junctions
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
                        #print splice_type,critical_exon, gene 
                    if splice_type =='' and e1a5_block<e1b5_block and e2a3_block>e2b3_block:
                        #(((2, 1), (2, 1)), ((9, 6), (9, 9))) (((4, 2), (4, 4)), ((7, 1), (7, 1))) ###one junction inside another
                        splice_type = "cassette-exon"; critical_exon = [e1b5,e2b3]
                        #"""

                        try: add_to_for_terminal_exons[gene,e2b3_block].append(e2a3)
                        except KeyError: add_to_for_terminal_exons[gene,e2b3_block] = [e2a3]                        
                        try: add_to_for_terminal_exons[gene,e1b5_block].append(e1a5)
                        except KeyError: add_to_for_terminal_exons[gene,e1b5_block] = [e1a5]
                        
                        try: cassette_exon_record[gene,e2b3_block].append(e1b5_block)
                        except KeyError: cassette_exon_record[gene,e2b3_block] = [e1b5_block]
                        try: cassette_exon_record[gene,e1b5_block].append(e2b3_block)
                        except KeyError: cassette_exon_record[gene,e1b5_block] = [e2b3_block]
                        #"""
                        y = CriticalExonInfo(gene,critical_exon,splice_type,splice_junctions)
                        data.write(gene+'\t'+reformatJunctions(critical_exon,'exon')+'\t'+splice_junction_str+'\t'+splice_type+'\n')
                        try: critical_exon_db[gene].append(y)
                        except KeyError: critical_exon_db[gene] = [y]
    data.close()
    ###Determine unique splice events and improve the annotations
    if export_annotation != '_de-novo':
        print len(critical_exon_db), 'genes identified from Ensembl, with alternatively regulated junctions'
    cassette_exon_record = eliminate_redundant_dict_values(cassette_exon_record)
    """for i in cassette_exon_record:
        print i, cassette_exon_record[i]"""
    add_to_for_terminal_exons = eliminate_redundant_dict_values(add_to_for_terminal_exons)
    ###Store all region information in a dictionary for efficient recall
    global region_db; region_db={}
    for gene in exon_regions:
        for rd in exon_regions[gene]:
            try: region_db[gene,rd.ExonRegionNumbers()]=rd
            except AttributeError: print gene, rd;kill
                
    if len(exon_regions) == 0:
        critical_exon_db_original = copy.deepcopy(critical_exon_db) ### get's modified somehow below
        #critical_exon_db_original = manualDeepCopy(critical_exon_db) ### won't work because it is the object that is chagned
    
    alternative_exon_db={}; critical_junction_db={}; critical_gene_junction_db={}
    alternative_terminal_exon={}
    for gene in critical_exon_db:
        critical_exon_junctions={}; critical_exon_splice_type={}
        for sd in critical_exon_db[gene]:
            for critical_exon in sd.CriticalExonRegion():
                try: critical_exon_junctions[critical_exon]+=sd.Junctions()
                except KeyError: critical_exon_junctions[critical_exon]=sd.Junctions()
                try: critical_exon_splice_type[critical_exon].append(sd.SpliceType())
                except KeyError: critical_exon_splice_type[critical_exon]=[sd.SpliceType()]                
                for junction in sd.Junctions():
                    try: critical_junction_db[tuple(junction)].append(junction)
                    except KeyError: critical_junction_db[tuple(junction)]=[sd.SpliceType()]
        critical_exon_junctions = eliminate_redundant_dict_values(critical_exon_junctions)
        critical_exon_splice_type = eliminate_redundant_dict_values(critical_exon_splice_type)
        for critical_exon in critical_exon_junctions:
            cj = critical_exon_junctions[critical_exon]
            splice_events = critical_exon_splice_type[critical_exon]
            status = 'stop'
            #print splice_events,critical_exon
            if splice_events == ['cassette-exon'] or ((gene,critical_exon[0]) in complex3prime_event_db) or ((gene,critical_exon[0]) in complex5prime_event_db):
                exons_blocks_joined_to_critical = cassette_exon_record[gene,critical_exon[0]]
                cassette_status = check_exon_polarity(critical_exon[0],exons_blocks_joined_to_critical)
                if len(critical_exon_junctions[critical_exon])<3 or cassette_status == 'no': ###Thus, it is not supported by 3 independent junctions
                    if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                        if cj[0][1] == cj[1][1]:
                            splice_events = ['alt-N-term']
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]; status = 'add_another'
                            alternative_terminal_exon[gene,critical_exon] = 'alt-N-term'
                        elif cj[0][0] == cj[1][0]:
                            splice_events = ['alt-C-term']
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]; status = 'add_another'
                            alternative_terminal_exon[gene,critical_exon] = 'alt-C-term'
                        else:
                            if critical_exon == cj[0][1]:
                                splice_events = ['alt-C-term'] ###this should be the only alt-exon
                                alternative_terminal_exon[gene,critical_exon] = 'alt-C-term'
                elif (gene,critical_exon[0]) in complex3prime_event_db:
                    #print '3prime',splice_events,critical_exon
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        #print critical_exon,len(cassette_exon_record[gene,critical_exon[0]]),cassette_exon_record[gene,critical_exon[0]];kill
                        if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                            splice_events = ['alt-N-term']; status = 'add_another'
                            alternative_terminal_exon[gene,critical_exon] = 'alt-N-term'
                elif (gene,critical_exon[0]) in complex5prime_event_db:
                    #print '5prime',splice_events,critical_exon
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        #print critical_exon,len(cassette_exon_record[gene,critical_exon[0]]),cassette_exon_record[gene,critical_exon[0]];kill
                        if len(exons_blocks_joined_to_critical)<2 or cassette_status == 'no':
                            second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                            splice_events = ['alt-C-term']; status = 'add_another'
                            alternative_terminal_exon[gene,critical_exon] = 'alt-C-term'
                """if 'mx-mx' in splice_events and (gene,critical_exon) in mx_event_db:
                    ###if one exon is a true cassette exon, then the mx-mx is not valid
                    if (gene,critical_exon[0]) in add_to_for_terminal_exons:
                        second_critical_exon = add_to_for_terminal_exons[gene,critical_exon[0]]
                        #print gene,critical_exon,second_critical_exon;kill"""
            splice_events = string.join(splice_events,'|'); exon_junction_str_list=[]
            splice_events = string.replace(splice_events, 'cassette-exons','cassette-exon(s)')
            ###if the junctions comprising the splice event for an alt-cassette show evidence of multiple exons, annotate as such
            if "alt5'-cassette" in splice_events: #or "cassette-alt3'"
                for junction in cj:
                    splicing_annotations = critical_junction_db[tuple(junction)]
                    if 'cassette-exons' in splicing_annotations:
                        splice_events = string.replace(splice_events, "alt5'-cassette","alt5'-cassette(s)"); break
            if "cassette-alt3'" in splice_events: #or "cassette-alt3'"
                for junction in cj:
                    splicing_annotations = critical_junction_db[tuple(junction)]
                    if 'cassette-exons' in splicing_annotations:
                        splice_events = string.replace(splice_events, "cassette-alt3'","cassette(s)-alt3'"); break
            ###Currently, 'cassette-exon(s)' is redundant with "alt5'-cassette(s)" or "cassette(s)-alt3'", so simplify
            if "alt5'-cassette(s)" in splice_events and 'cassette-exon(s)' in splice_events:
                splice_events = string.replace(splice_events, 'cassette-exon(s)','')
            if "cassette(s)-alt3'" in splice_events and 'cassette-exon(s)' in splice_events:
                splice_events = string.replace(splice_events, 'cassette-exon(s)','')
            splice_events = string.replace(splice_events, '||','|')
            for j in cj:
                nj=[]
                for exon in j: e = 'E'+str(exon[0])+'.'+str(exon[1]); nj.append(e)
                try: critical_gene_junction_db[gene].append(nj)
                except KeyError: critical_gene_junction_db[gene] = [nj]
                nj = string.join(nj,'-')
                exon_junction_str_list.append(nj)
            exon_junction_str = string.join(exon_junction_str_list,'|')
            try:
                rd = region_db[gene,critical_exon] ###('ENSG00000213588', (26, 1)) and ('ENSG00000097007', (79, 1)) didn't work
            except KeyError:
                ###Occurs as a results of either exons or transcripts with sketchy or complex assignments
                null = []
            try:
                se = rd.AssociatedSplicingEvent()
                if len(se)>1:
                    if splice_events not in se: se = se+'|'+ splice_events
                else: se = splice_events
                rd.setSpliceData(se,exon_junction_str)
        
                #print critical_exon,se,exon_junction_str,gene
                if status == 'add_another':
                    for critical_exon in second_critical_exon: 
                        rd = region_db[gene,critical_exon]
                        se = rd.AssociatedSplicingEvent()
                        if len(se)>1:
                            if splice_events not in se:
                                if se != 'cassette-exon': se = se+'|'+ splice_events
                        else: se = splice_events
                        rd.setSpliceData(se,exon_junction_str)
                        #print critical_exon, se, exon_junction_str, gene,'second'
                """
                ###create an index for easy searching of exon content in downstream modules
                critical_exon_block = critical_exon[0]
                if gene in alternative_exon_db:
                    block_db = alternative_exon_db[gene]
                    try: block_db[critical_exon_block].append(rd)
                    except KeyError: block_db[critical_exon_block] = [rd]
                else:
                    block_db = {}; block_db[critical_exon_block]=[rd]
                    alternative_exon_db[gene]=block_db"""
            except Exception: null=[] ### Occurs when analyzing novel junctions, rather than Ensembl
            
    #alternative_terminal_exon[gene,critical_exon] = 'alt-C-term'
    ###Since setSpliceData will update the existing instance, we can just re-roder the region db for easy searching in downstream modules
    ### (note: the commented out code above could be useful for exon-structure output)
    for gene in exon_regions:
        block_db = {}
        for rd in exon_regions[gene]:
            try: block_db[rd.ExonNumber()].append(rd)
            except KeyError: block_db[rd.ExonNumber()] = [rd]
        exon_regions[gene] = block_db ###Replace the existing list with a dictionary for faster look-ups
    
    if len(exon_regions)==0: exon_regions = critical_exon_db_original ### See JunctionArray.inferJunctionComps()
    
    return exon_regions,critical_gene_junction_db

def manualDeepCopy(db):
    ### Same as deep copy, possibly less memory intensive
    db_copy={}
    for i in db:
        db_copy[i] = list(tuple(db[i]))
    return db_copy

def check_exon_polarity(critical_exon_block,exons_blocks_joined_to_critical):
    g=0;l=0
    for joined_exon_blocks in exons_blocks_joined_to_critical:
        if joined_exon_blocks>critical_exon_block: g+=1
        if joined_exon_blocks<critical_exon_block: l+=1
    if g>0 and l>0: cassette_status = 'yes'
    else: cassette_status = 'no'
    return cassette_status

def pickOptimalCriticalExons(e1b5,e1a5,e2a3,e2b3,critical_exon):
    e1a5_block,e1a5_reg = e1a5
    e2a3_block,e2a3_reg = e2a3
    e1b5_block,e1b5_reg = e1b5
    e2b3_block,e2b3_reg = e2b3

    if e1a5_block == e1b5_block:                    
        if e1a5_reg < e1b5_reg: critical_exon += [e1b5]
        elif e1a5_reg != e1b5_reg: critical_exon += [e1a5]
    if e2a3_block == e2b3_block:
        if e2a3_reg < e2b3_reg: critical_exon += [e2a3]
        elif e2a3_reg != e2b3_reg: critical_exon += [e2b3]
    return critical_exon

def customDBDeepCopy(db):
    db2={}
    for i in db:
        for e in db[i]:
            try: db2[i].append(e)
            except KeyError: db2[i]=[e]
    return db2

def customLSDeepCopy(ls):
    ls2=[]
    for i in ls: ls2.append(i)
    return ls2

def importExonRegionCoordinates(species):
    """ Used to export Transcript-ExonRegionIDs """
    filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    fn=filepath(filename); x=0
    exon_region_coord_db={}
    all_coord_db={}
    exon_region_db={}
    strand_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x+=1
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            start_end = [start, stop]; start_end.sort()
            all_coord_db[gene,exonid] = start_end
            try:
                ### start and stop should be unique
                db = exon_region_coord_db[gene]
                db['s',float(start)] = exonid ### start be the same as another region stop - designate as start
                db['e',float(stop)] = exonid ### start be the same as another region stop - designate as end
            except Exception:
                db={}
                db['s',float(start)] = exonid
                db['e',float(stop)] = exonid
                exon_region_coord_db[gene] = db
            #if 'I' not in exonid: ### if it spans an intron, include it
            try: exon_region_db[gene].append((exonid,start))
            except Exception: exon_region_db[gene] = [(exonid,start)]
            strand_db[gene] = strand
    return exon_region_coord_db, all_coord_db, exon_region_db, strand_db

def exportTranscriptExonIDAssociations(species):
    """ Determine the exon region ID composition of all analyzed mRNAs """
    try: relative_exon_locations = importEnsExonStructureDataSimpler(species,'ucsc',{})
    except Exception: relative_exon_locations={}
    relative_exon_locations = importEnsExonStructureDataSimpler(species,'ensembl',relative_exon_locations)
    from build_scripts import IdentifyAltIsoforms
    seq_files, transcript_protein_db = IdentifyAltIsoforms.importProteinSequences(species,just_get_ids=True) ### get mRNA-protein IDs
    
    exon_region_coord_db, all_coord_db, exon_region_db, strand_db = importExonRegionCoordinates(species)
    export_file = 'AltDatabase/ensembl/'+species+'/mRNA-ExonIDs.txt'
    export_data = export.ExportFile(export_file)
    transcript_region_db={}
    transcripts_with_missing_regions={}
    errorCount=0
    retainedIntrons=0
    for key in relative_exon_locations:  ### From UCSC or Ensembl transcript coordinates only
        all_exon_intron_regions = {}
        ens_transcriptid,gene,strand = key
        #if ens_transcriptid != 'AK185721': continue
        region_db={} ### track the coordinates for cleaning up the order
        coord_db = exon_region_coord_db[gene]
        for (region,start) in exon_region_db[gene]:
            all_exon_intron_regions[region] = all_coord_db[gene,region]
            region_db[region] = start
        for exon_data in relative_exon_locations[key]: ### each transcript exon
            regions=[]
            exon_start,exon_end,ens_exonid = exon_data
            start_end = [exon_start,exon_end]; start_end.sort(); start,end = start_end
            partial_matches=[]
            added=False
            for region in all_exon_intron_regions: ### search each transcript exon against every annotated region
                e5,e3 = all_exon_intron_regions[region]
                annotated = [int(e5),int(e3)]
                annotated.sort()
                coords = [start_end[0],start_end[1]]+annotated
                coords.sort()
                if coords[0] == start_end[0] and coords[-1] == start_end[-1]:
                    ### Hence, the exon/intron region is contained within or is equal to the transcript exon
                    if (gene,ens_transcriptid) in transcript_region_db:
                        if region not in transcript_region_db[gene,ens_transcriptid]:
                            transcript_region_db[gene,ens_transcriptid].append((region))
                    else:
                        transcript_region_db[gene,ens_transcriptid] = [region]
                    added=True
                    retainedIntrons+=1
                else:
                    if annotated[0] == start_end[0] or annotated[-1] == start_end[-1]:
                        #print exon_start,exon_end, start_end, region, ens_transcriptid
                        partial_matches.append(region)
            if added ==False:
                if len(partial_matches)>0:
                    ### Occurs typically when a UCSC exon has a longer or shorter 3' or 5'UTR exon boundaries
                    for region in partial_matches:
                        if (gene,ens_transcriptid) in transcript_region_db:
                            if region not in transcript_region_db[gene,ens_transcriptid]:
                                transcript_region_db[gene,ens_transcriptid].append((region))
                        else:
                            transcript_region_db[gene,ens_transcriptid] = [region]
                else:
                    transcripts_with_missing_regions[ens_transcriptid]=None
                    errorCount+=1
                    #if errorCount<100: print 'Error:',key, exon_start,exon_end,ens_exonid
        if (gene,ens_transcriptid) in transcript_region_db:
            regions = transcript_region_db[gene,ens_transcriptid]
            regions_sort=[]
            for region in regions:
                try: start = region_db[region]
                except Exception,e: print e, coord_db;sys.exit()
                regions_sort.append([start,region])
            regions_sort.sort()
            if strand_db[gene] == '-':regions_sort.reverse()
            transcript_region_db[gene,ens_transcriptid] = map(lambda (s,r): r, regions_sort)
            #print gene, ens_transcriptid, transcript_region_db[gene,ens_transcriptid];sys.exit()
    print len(transcripts_with_missing_regions), 'transcripts with missing exon regions out of', len(transcript_region_db)+len(transcripts_with_missing_regions)

    t1=[]
    for i in transcripts_with_missing_regions:
        t1.append(i)
    print 'missing:',t1[:15]

    exon_db={}
    gene_region_db={}
    for (gene,transcript) in transcript_region_db:
        try: proteinAC = transcript_protein_db[transcript]
        except Exception: proteinAC = 'None'
        try: regions = string.join(transcript_region_db[gene,transcript],'|')
        except Exception: print gene, transcript, transcript_region_db[gene,transcript];sys.exit()
        gene_region_db[gene,regions]=[]
        export_data.write(string.join([gene,transcript,proteinAC,regions],'\t')+'\n')
        for exonID in transcript_region_db[gene,transcript]:
            exon_db[gene+':'+exonID]=None
    export_data.close()

    print 'Unique-transcripts by regionID makeup:',len(gene_region_db) ### Complete UCSC gives 272071 versus 237607 for the non-Complete
    filterExonRegionSeqeunces(exon_db,species)

def filterExonRegionSeqeunces(exon_db,species):
    filename = 'AltDatabase/'+species+'/RNASeq/RNASeq_critical-exon-seq_updated.txt'
    export_file = 'AltDatabase/'+species+'/RNASeq/RNASeq_critical-exon-seq_filtered.txt'
    print 'importing', filename
    fn=filepath(filename)
    export_data = export.ExportFile(export_file)
    for line in open(fn,'r').xreadlines():
        data = line.strip()
        t = string.split(data,'\t')
        exonID = t[0]; sequence = t[-1]
        try:
            y = exon_db[exonID]
            export_data.write(line)
        except Exception: null=[] ### Occurs if there is no Ensembl for the critical exon or the sequence is too short to analyze
    export_data.close()

def createExonRegionSequenceDB(species,platform):
    """ Store the filtered exon sequence data in an SQL database for faster retreival """
    start=time.time()
    import SQLInterface
    DBname = 'ExonSequence'
    schema_text ='''-- Schema for species specific AltAnalyze transcript data.

-- Genes store general information on each Ensembl gene ID
create table ExonSeq (
    uid          text primary key,
    gene         text,
    sequence     text
);
'''
    conn = SQLInterface.populateSQLite(species,platform,DBname,schema_text=schema_text) ### conn is the database connnection interface
    
    ### Populate the database
    filename = 'AltDatabase/'+species+'/RNASeq/RNASeq_critical-exon-seq_filtered.txt'
    print 'importing', filename
    fn=filepath(filename)
    for line in open(fn,'r').xreadlines():
        data = line.strip()
        t = string.split(data,'\t')
        exonID = t[0]; sequence = t[-1]
        gene,region = string.split(exonID,':')
        #print exonID,gene,sequence
        ### Store this data in the SQL database
        command = """insert into ExonSeq (uid, gene, sequence)
        values ('%s', '%s','%s')""" % (exonID,gene,sequence)
        conn.execute(command)
        
    conn.commit() ### Needed to commit changes
    conn.close()
    time_diff = str(round(time.time()-start,1))
    print 'Exon Region Sequences added to SQLite database in %s seconds' % time_diff

def importTranscriptExonIDs(species):
    start=time.time()
    filename = 'AltDatabase/ensembl/'+species+'/mRNA-ExonIDs.txt'
    fn=filepath(filename)
    gene_transcript_structure={}
    protein_db = {}
    for line in open(fn,'r').xreadlines():
        data = line.strip()
        gene,transcript,proteinAC,regions = string.split(data,'\t')
        if gene in gene_transcript_structure:
            tdb=gene_transcript_structure[gene]
            tdb[transcript] = regions
        else:
            tdb={}
            tdb[transcript] = regions+'|'
            gene_transcript_structure[gene] = tdb
        protein_db[proteinAC] = transcript
            
    time_diff = str(round(time.time()-start,1))
    #print 'Transcript-ExonID associations imported in %s seconds' % time_diff
    return gene_transcript_structure, protein_db
        
def identifyPCRregions(species,platform,uid,inclusion_junction,exclusion_junction,isoform1,isoform2):
    #print uid,inclusion_junction,exclusion_junction,isoform1,isoform2;sys.exit()
    try:
        x = len(gene_transcript_structure)
        print_outs = False
    except Exception:
        gene_transcript_structure, protein_db = importTranscriptExonIDs(species)
        print_outs = True

    gene,region = string.split(uid,':')
    isoform_db = copy.deepcopy(gene_transcript_structure[gene])
    
    import SQLInterface
    conn = SQLInterface.connectToDB(species,platform,'ExonSequence')
    ids = [gene]
    query = "select uid, sequence from ExonSeq where gene = ?"
    uid_sequence_list = SQLInterface.retreiveDatabaseFields(conn,ids,query)
    
    ex1,ex2 = string.split(exclusion_junction,'-')
    ex1b,ex2b = string.split(inclusion_junction,'-')

    exon_seq_db={}
    for (uid1,seq) in uid_sequence_list:
        exon_seq_db[uid1] = seq
    try: 
        print '('+ex1+')'+exon_seq_db[gene+':'+ex1]
        print '('+ex2+')'+exon_seq_db[gene+':'+ex2]
        print '('+ex1b+')'+exon_seq_db[gene+':'+ex1b]
        print '('+ex2b+')'+exon_seq_db[gene+':'+ex2b]
    except Exception:
        pass
    
    if print_outs == True:
        #"""
        for (uid1,seq) in uid_sequence_list:
            if uid1 == uid:
                print seq
        #"""
    try:
        mRNA1 = protein_db[isoform1]
        mRNA1_s = isoform_db[mRNA1]
    except Exception:
        mRNA1_s = string.replace(inclusion_junction,'-','|')
    try:
        mRNA2 = protein_db[isoform2]
        mRNA2_s = isoform_db[mRNA2]
    except Exception:
        mRNA2_s = string.replace(exclusion_junction,'-','|')

    print ex1,ex2
    print [mRNA1_s]
    print [mRNA2_s]
    if mRNA1_s != None:
        ex1_pos = string.find(mRNA1_s,ex1+'|') ### This is the location in the string where the exclusion junctions starts
        ex1_pos = ex1_pos+1+string.find(mRNA1_s[ex1_pos:],'|') ### This is the location in the string where the inclusion exon starts
        ex2_pos = string.find(mRNA1_s,ex2+'|')-1 ### This is the location in the string where the inclusion exon ends
    print ex2_pos, ex1_pos
    if ex2_pos<ex1_pos:
        if mRNA2_s != None:
            mRNA1_s = mRNA2_s
            #mRNA1 = mRNA2
            ex1_pos = string.find(mRNA1_s,ex1+'|') ### This is the location in the string where the exclusion junctions starts
            ex1_pos = ex1_pos+1+string.find(mRNA1_s[ex1_pos:],'|') ### This is the location in the string where the inclusion exon starts
            ex2_pos = string.find(mRNA1_s,ex2+'|')-1 ### This is the location in the string where the inclusion exon ends
    if abs(ex1_pos-ex2_pos)<2:
        ### Incorrect isoform assignments resulting in faulty primer design
        if ex1b+'|' in mRNA2_s and ex2b+'|' in mRNA2_s:
            ex1 = ex1b
            ex2 = ex2b
            ex1_pos = string.find(mRNA1_s,ex1+'|') ### This is the location in the string where the exclusion junctions starts
            ex1_pos = ex1_pos+1+string.find(mRNA1_s[ex1_pos:],'|') ### This is the location in the string where the inclusion exon starts
            ex2_pos = string.find(mRNA1_s,ex2+'|')-1 ### This is the location in the string where the inclusion exon ends      

    inclusion_exons = string.split(mRNA1_s[ex1_pos:ex2_pos],'|') ### These are the missing exons from the exclusion junction (in between)
    #if '-' in mRNA1_s:
    common_exons5p = string.split(mRNA1_s[:ex1_pos-1],'|') ### Flanking full 5' region (not just the 5' exclusion exon region)
    if (ex2_pos+1) == -1: ### Hence the 3' exon is the last exon in the mRNA
        common_exons3p = [string.split(mRNA1_s,'|')[-1]]
        inclusion_exons = string.split(mRNA1_s[ex1_pos:],'|')[:-1]
    else:
        common_exons3p = string.split(mRNA1_s[ex2_pos+1:],'|')  ### Flanking full 3' region (not just the 3' exclusion exon region)
    if gene == 'E1NSG00000205423':
        #print mRNA1_s, ex2;sys.exit()
        #print mRNA1_s;sys.exit()
        print uid,inclusion_junction,exclusion_junction,isoform1,isoform2
        print inclusion_exons
        print common_exons5p, common_exons3p
        print mRNA1_s, ex2_pos, ex1, ex2
        sys.exit()

    inclusion_exons = map(lambda x: gene+':'+x, inclusion_exons) ### add geneID prefix
    common_exons5p = map(lambda x: gene+':'+x, common_exons5p) ### add geneID prefix
    common_exons3p = map(lambda x: gene+':'+x, common_exons3p) ### add geneID prefix
    inclusion_junction = string.replace(inclusion_junction,'-','|')+'|'
    exclusion_junction = string.replace(exclusion_junction,'-','|')+'|'
    
    #if inclusion_junction in mRNA1_s: print '1 true'
    #if exclusion_junction in mRNA2_s: print '2 true'

    #sys.exit()

    uid_seq_db={} ### convert list to dictionary
    for (uid,seq) in uid_sequence_list:
        uid_seq_db[uid] = seq
        #E213.1-E214.1 vs. E178.1-E223.1
        if uid == 'ENSG00000145349:E13.1': print uid, seq, len(seq)
        if uid == 'ENSG00000145349:E12.1': print uid, seq, len(seq)
        if uid == 'ENSG00000145349:E19.1': print uid, seq, len(seq)
        if uid == 'ENSG00000145349:E15.2': print uid, seq, len(seq)
        if uid == 'ENSG00000145349:E18.1': print uid, seq, len(seq)
        
    #sys.exit()
    ### Get the common flanking and inclusion transcript sequence
    #print common_exons5p,common_exons3p;sys.exit()
    print 1
    common_5p_seq = grabTranscriptSeq(common_exons5p,uid_seq_db)
    print 2
    common_3p_seq = grabTranscriptSeq(common_exons3p,uid_seq_db)
    print 3
    inclusion_seq = grabTranscriptSeq(inclusion_exons,uid_seq_db)
    print 'common_5p_seq:',[common_5p_seq]
    print 'common_3p_seq:',[common_3p_seq]
    print 'inclusion_seq:',[inclusion_seq], inclusion_exons
    incl_isoform_seq = common_5p_seq+inclusion_seq+common_3p_seq
    incl_isoform_seq_formatted = common_5p_seq[-100:]+'['+inclusion_seq+']'+common_3p_seq[:100]
    excl_isoform_seq = common_5p_seq+common_3p_seq

    c5pl=len(common_5p_seq)
    c3pl=len(common_3p_seq)
    ipl1=len(inclusion_seq)
    il=len(incl_isoform_seq)
    
    # Metrics to designate minimal search regions for primer3 > release 2.3
    if c5pl>100: s1 = c5pl-100; e1 = 100 ### start looking for primers in the region (forward)
    else: s1 = 0; e1 = c5pl
    #if c3pl>200: s2 = c5pl+ipl1; e2 = 200
    if c3pl>100: s2 = c5pl+ipl1; e2 = 100
    else: s2 = c5pl+ipl1; e2 = c3pl
    include_region = [s1,s2+e2]
    include_region = [s1,e1+ipl1+e2]
    include_region = map(lambda x: str(x), include_region)
    target_region = [c5pl,ipl1]
    target_region = map(lambda x: str(x), target_region)

    input_dir = filepath('AltDatabase/primer3/temporary-files/temp1.txt')
    output_file = filepath('AltDatabase/primer3/temporary-files/output1.txt')
    try: os.remove(input_dir)
    except Exception: pass
    try: os.remove(output_file)
    except Exception: pass
    #print incl_isoform_seq
    #print include_region
    #print target_region;sys.exit()
    
    input_dir = exportPrimerInputSeq(incl_isoform_seq,include_region,target_region)
    primer3_file = getPrimer3Location()
    #for Primer3 release 2.3 and greater: SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100
    #Left primer in the 50 bp region starting at position 100 AND right primer in the 50 bp region starting at position 300

    ### Run Primer3
    command_line = [primer3_file, "<",'"'+input_dir+'"', ">",'"'+output_file+'"'] #-format_output
    #command_line = [primer3_file, "-format_output <",'"'+input_dir+'"', ">",'"'+output_file+'"']
    command_line = string.join(command_line,' ')
    #print [command_line];sys.exit()
    retcode = os.popen(command_line)
    time.sleep(4)
    """
    commandFinshed = False
    while commandFinshed == False:
        try: commandFinshed = checkFileCompletion(output_file)
        except Exception,e: pass
    """
        
    left_primer,right_primer,amplicon_size = importPrimer3Output(output_file,gene)
    primers = 'F:'+left_primer+' R:'+right_primer+' sizes:'+str(amplicon_size)+'|'+str(amplicon_size-ipl1) + '  (inclusion-isoform:%s)' % incl_isoform_seq_formatted
    right_primer_sense = reverse_orientation(right_primer)
    if left_primer in excl_isoform_seq:
        print 'Left',
    else: kill
    if right_primer_sense in excl_isoform_seq:
        print 'Right'
    else: kill
    #if print_outs:
    print primers
    return primers
   
def checkFileCompletion(fn):
    complete = False
    for line in open(fn,'r').xreadlines():
        if 'PRIMER_PRODUCT_SIZE' in line: complete = True
    return complete

def grabTranscriptSeq(exons,uid_seq_db):
    seq=''
    #print exons
    for uid in exons:
        try: seq += uid_seq_db[uid]
        except Exception: break ### MISSING SEQUENCE FROM THE DATABASE - OFTEN OCCURS IN THE LAST FEW EXONS - UNCLEAR WHY
    return seq

def exportPrimerInputSeq(sequence,include_region,target_region):
    tempdir = filepath('AltDatabase/primer3/temporary-files/temp1.txt')
    export_obj = export.ExportFile(tempdir)
    export_obj.write('PRIMER_SEQUENCE_ID=temp1\nSEQUENCE=')
    export_obj.write(sequence+'\n')
    export_obj.write('INCLUDED_REGION=')
    export_obj.write(string.join(include_region,',')+'\n')
    export_obj.write('PRIMER_PRODUCT_SIZE_RANGE=75-1000\n')
    
    export_obj.write('TARGET=')
    export_obj.write(string.join(target_region,',')+'\n')
    export_obj.write('=')
    export_obj.close()
    return tempdir

def importPrimer3Output(fn,gene):
    gene_transcript_structure={}
    protein_db = {}
    for line in open(fn,'r').xreadlines():
        data = line.strip()
        if 'PRIMER_LEFT_SEQUENCE=' in data:
            left_primer = string.split(data,'PRIMER_LEFT_SEQUENCE=')[-1]
            #print left_primer;
            #print fn; sys.exit()
        if 'PRIMER_RIGHT_SEQUENCE=' in data:
            right_primer = string.split(data,'PRIMER_RIGHT_SEQUENCE=')[-1]
        if 'PRIMER_PRODUCT_SIZE=' in data:
            amplicon_size = int(string.split(data,'PRIMER_PRODUCT_SIZE=')[-1])
            break
    return left_primer,right_primer,amplicon_size

def getPrimer3Location():
    primer3_dir = 'AltDatabase/primer3/'
    if os.name == 'nt':
        if '32bit' in architecture: primer3_file = primer3_dir + '/PC/32bit/primer3_core'; plat = 'Windows'
        elif '64bit' in architecture: primer3_file = primer3_dir + '/PC/64bit/primer3_core'; plat = 'Windows'
    elif 'darwin' in sys.platform: primer3_file = primer3_dir + '/Mac/primer3_core'; plat = 'MacOSX'
    elif 'linux' in sys.platform:
        if '32bit' in platform.architecture(): primer3_file = primer3_dir + '/Linux/32bit/primer3_core'; plat = 'linux32bit'
        elif '64bit' in platform.architecture(): primer3_file = primer3_dir + '/Linux/64bit/primer3_core'; plat = 'linux64bit'
    primer3_file = filepath(primer3_file)
    return primer3_file

def importComparisonSplicingData4Primers(filename,species):
    fn=filepath(filename)
    stringent_regulated_exons = {}
    firstLine = True
    global gene_transcript_structure
    global protein_db
    gene_transcript_structure, protein_db = importTranscriptExonIDs(species)
    
    ei = export.ExportFile(filename[:-4]+'-primers.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if firstLine:
            header = t
            firstLine = False
            ei.write(line)
        else:
            try:
                if 'comparison' in header:
                    t = t[:-1]
                if 'PSI' not in filename:
                    exonid = t[0]; symbol = t[2]; confirmed = t[5]; fold_change = abs(float(t[-2])); percent_exp = float(t[-3])
                    junctions = t[4]; isoforms = t[9]; splice_type = t[8]
                else:
                    """ Update the code to work with PSI results files from the metaDataAnalysis script """
                    #UID,InclusionNotation,ClusterID,UpdatedClusterID,AltExons,EventAnnotation,Coordinates,ProteinPredictions,dPSI,rawp,adjp,avg1,avg2
                    uid = t[0]
                    print uid
                    uid_objects = string.split(uid,':')
                    symbol = uid_objects[0]
                    junctions = string.join(uid_objects[1:],':')
                    junctions = string.split(junctions,'|')
                    fold_change = abs(float(t[9]))
                    isoforms = t[8]
                    splice_type = t[6]
                    exonid = t[5]
                #print symbol, junctions,fold_change,percent_exp,confirmed;sys.exit()
                #if fold_change<2 and percent_exp>0.25 and (confirmed == 'yes'):
                if fold_change < 50:
                    if 'alternative_polyA' not in splice_type and 'altPromoter' not in splice_type:
                        if len(junctions)==2:
                            j1, j2 = junctions
                            if 'ENS' in j1:
                                j1 = string.split(j1,':')[1]
                                j2 = string.split(j2,':')[1]
                        else:
                            j1, j2 = string.split(string.split(junctions,'|')[0],' vs. ')
                        #(-)alt-C-terminus,(-)AA:188(ENSP00000397452)->238(ENSP00000410667),(-)microRNA-target(hsa-miR-599:miRanda,hsa-miR-186:miRanda)
                        try:
                            iso1, iso2 = string.split(string.split(isoforms,'AA:')[1],')->')
                            iso1 = string.split(iso1,'(')[1]
                            iso2 = string.split(string.split(iso2,'(')[1],')')[0]
                            #print iso1, iso2
                        except:
                            iso1 = ''
                            iso2 = ''
                        #print j1, j2
                        #print symbol
                        try:
                            primer = identifyPCRregions(species,'RNASeq',exonid,j1,j2,iso1,iso2)
                            #print exonid, j1, j2, iso1, iso2
                            ei.write(string.join(t+[primer],'\t')+'\n')
                            #print primer, symbol, exonid
                            #sys.exit()
                        except Exception:
                            pass
                            print traceback.format_exc(),'\n'; #sys.exit()
                        #sys.exit()
            except Exception:
                #print traceback.format_exc(),'\n'#;sys.exit()
                pass
    ei.close()
                
if __name__ == '__main__':
    ###KNOWN PROBLEMS: the junction analysis program calls exons as cassette-exons if there a new C-terminal exon occurs downstream of that exon in a different transcript (ENSG00000197991).
    Species = 'Hs'
    test = 'yes'
    Data_type = 'ncRNA'
    Data_type = 'mRNA'
    #E6.1-E8.2 vs. E5.1-E8.3
    filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Ichi/August.11.2017/Events-dPSI_0.0_rawp/PSI.R636S_Homo_vs_WTC-limma-updated-Domains2-filtered.txt'
    #filename = '/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/July-2017/PSI/Events-dPSI_0.1_adjp/PSI.U2AF1-like_vs_OthersQPCR.txt'
    #filename = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Ichi/Combined-junction-exon-evidence.txt'
    
    #exportTranscriptExonIDAssociations(Species);sys.exit()
    #createExonRegionSequenceDB(Species,'RNASeq'); sys.exit()
    importComparisonSplicingData4Primers(filename,Species); sys.exit()
    
    #ENSG00000154556:E8.2-E9.7|ENSG00000154556:E5.12-E9.7	alt-N-terminus|-	alt-C-terminus|-	AA:41(ENST00000464975-PEP)->43(ENST00000493709-PEP)|-
    #ENSG00000161999:E2.3-E2.5	nonsense_mediated_decay|+	retained_intron|+	alt-N-terminus|+	alt-C-terminus|+	AA:72(ENST00000564436-PEP)->81(ENSP00000454700)|+
    #ENSG00000122591:E12.2-E13.1|ENSG00000122591:E12.1-E13.1	alt-N-terminus|+	alt-C-terminus|+	AA:116(ENST00000498833-PEP)->471(ENSP00000397168)|+

    #ENSG00000128891:E1.11-E2.1|ENSG00000128891:E1.10-E2.1	alt-N-terminus|+	AA:185(ENSP00000350695)->194(ENSP00000452773)|+
    identifyPCRregions(Species,'RNASeq','ENSG00000128891:E1.12','E1.12-E2.1','E1.11-E2.1','ENSP00000350695','ENSP00000350695'); sys.exit()
    #identifyPCRregions(Species,'RNASeq','ENSG00000133226:E9.1','E8.1-E9.1','E8.1-E11.3','ENST00000564436-PEP','ENSP00000454700'); sys.exit()
    #createExonRegionSequenceDB(Species,'RNASeq'); sys.exit()
    
    getEnsemblAssociations(Species,Data_type,test); sys.exit()
    
    test_gene = ['ENSG00000143776']#,'ENSG00000154889','ENSG00000156026','ENSG00000148584','ENSG00000063176','ENSG00000126860'] #['ENSG00000105968']
    meta_test = ["ENSG00000215305","ENSG00000179676","ENSG00000170484","ENSG00000138180","ENSG00000100258","ENSG00000132170","ENSG00000105767","ENSG00000105865","ENSG00000108523","ENSG00000150045","ENSG00000156026"]
    #test_gene = meta_test
    gene_seq_filename = 'AltDatabase/ensembl/'+Species+'/'+Species+'_gene-seq-2000_flank'
    gene_db = import_sequence_data(gene_seq_filename,{},Species,'gene_count'); print len(gene_db); sys.exit()
    import_dir = '/AltDatabase/ensembl/'+species
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for file_name in dir_list:    #loop through each file in the directory to output results
        dir_file = 'AltDatabase/ensembl/'+species+'/'+file_name
        if 'exon' in dir_file: exon_file = dir_file
        elif 'transcript' in dir_file: trans_file = dir_file
        elif 'gene' in dir_file: gene_seq_file = dir_file
        elif 'Exon_cDNA' in dir_file: exon_trans_file = dir_file
        elif 'Domain' in dir_file: domain_file = dir_file
    #"""
    exon_annotation_db,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db,ucsc_splicing_annot_db = getEnsExonStructureData(species,data_type)

    #exon_db = customDBDeepCopy(exon_annotation_db)
    #"""
    exon_annotation_db2 = annotate_exons(exon_annotation_db)
    #kill
    exon_db2 = customDBDeepCopy(exon_annotation_db2) ##having problems with re-writting contents of this db when I don't want to
    
    exon_clusters,intron_clusters,exon_regions,intron_region_db = exon_clustering(exon_db2); exon_db2={}
    #"""
    exon_junction_db,putative_as_junction_db,exon_junction_db = processEnsExonStructureData(exon_annotation_db,exon_regions,transcript_gene_db,gene_transcript,transcript_exon_db,intron_retention_db)
    s
    #ej = {};ej['ENSG00000124104'] = exon_junction_db['ENSG00000124104']
    #exon_regions,critical_gene_junction_db = compareJunctions(putative_as_junction_db,exon_regions)
    #exportSubGeneViewerData(exon_regions,critical_gene_junction_db,intron_region_db,intron_retention_db)        
    #ENSG00000149792 possible retained 3' intron... check out
    kill
    use_exon_data='no';get_splicing_factors = 'yes'
    rna_processing_ensembl = GO_parsing.parseAffyGO(use_exon_data,get_splicing_factors,species)
    ensembl_annot_file = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl-annotations_simple.txt'
    ensembl_annotation_db = getEnsemblAnnotations(ensembl_annot_file,rna_processing_ensembl)    
    
    ###Print ovelap statistics for exon-blocks
    x=0;y=0; m=0; l=0
    for key in exon_clusters:
        for exon_block in exon_clusters[key]:
            if len(exon_block[2])>1: y += 1; x += 1; m += len(exon_block[2]); l += len(exon_block[2])
            else: x += 1; m += 1
            #if x < 50:
                #print key[0], exon_block[2], len(exon_block[2]),x,y,m,l
    print 'x',x,'y',y,'m',m,'l',l
    """
    for gene in exon_regions:
        db = exon_regions[gene]
        for block in db:
            for rd in db[block]:
                try: print rd.AssociatedSplicingEvent(),rd.AssociatedSplicingJunctions();kill
                except AttributeError: continue"""
