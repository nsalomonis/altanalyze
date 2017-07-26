###ResultsExport_module
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
from stats_scripts import statistics
import unique
import export
dirfile = unique

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    #add in code to prevent folder names from being included
    dir_list2 = [] 
    for entry in dir_list:
        #if entry[-4:] == ".txt" or entry[-4:] == ".all" or entry[-5:] == ".data" or entry[-3:] == ".fa":
        dir_list2.append(entry)
    return dir_list2

def returnDirectories(sub_dir):
    dir=os.path.dirname(dirfile.__file__)
    dir_list = os.listdir(dir + sub_dir)
    ###Below code used to prevent FILE names from being included
    dir_list2 = []
    for entry in dir_list:
        if "." not in entry: dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

class GrabFiles:
    def setdirectory(self,value): self.data = value
    def display(self): print self.data
    def searchdirectory(self,search_term):
        #self is an instance while self.data is the value of the instance
        files = getDirectoryFiles(self.data,search_term)
        if len(files)<1: print 'files not found'
        return files
    def returndirectory(self):
        dir_list = getAllDirectoryFiles(self.data)
        return dir_list

def getAllDirectoryFiles(import_dir):
    all_files = []
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        all_files.append(data_dir)
    return all_files

def getDirectoryFiles(import_dir,search_term):
    dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
    matches=[]
    for data in dir_list:    #loop through each file in the directory to output results
        data_dir = import_dir[1:]+'/'+data
        if search_term not in data_dir: matches.append(data_dir)
    return matches

############ Result Export Functions #############

def outputSummaryResults(summary_results_db,name,analysis_method,root_dir):
    #summary_results_db[dataset_name] = udI,udI-up_diff,ddI,ddI-down_diff,udI_mx,udI_mx-mx_diff,up_dI_genes,down_gene, annotation_list
    annotation_db = {}
    for dataset in summary_results_db:
        for entry in summary_results_db[dataset][-1]:
            annotation = entry[0]
            count = entry[1]
            if 'AA:' not in annotation:
                try: annotation_db[annotation].append((dataset,count))
                except KeyError: annotation_db[annotation] = [(dataset,count)]
    annotation_ls = []

    for annotation in annotation_db: annotation_ls.append(annotation)
    annotation_ls.sort()
    annotation_db2={}
    for annotation in annotation_ls:
        for dataset in summary_results_db:
            y=0
            for entry in summary_results_db[dataset][-1]:
                annotation2 = entry[0]
                count = entry[1]
                if annotation2 == annotation:
                    y=1; new_count = count
            if y == 1:
                try: annotation_db2[dataset].append((annotation,new_count))
                except KeyError: annotation_db2[dataset] = [(annotation,new_count)]
            else:
                try: annotation_db2[dataset].append((annotation,0))
                except KeyError: annotation_db2[dataset] = [(annotation,0)]
      
    summary_output = root_dir+'AltResults/AlternativeOutput/'+analysis_method+'-summary-results'+name+'.txt'
    fn=filepath(summary_output)
    data = export.createExportFile(summary_output,'AltResults/AlternativeOutput')
    if analysis_method == 'splicing-index' or analysis_method == 'FIRMA':
        event_type1 = 'inclusion-events'; event_type2 = 'exclusion-events'; event_type3 = 'alternative-exons'
    else:
        event_type1 = 'inclusion-events'; event_type2 = 'exclusion-events'; event_type3 = 'mutually-exlusive-events'
    title = 'Dataset-name' +'\t'+ event_type1+'\t'+event_type2 +'\t'+ event_type3 +'\t'+ 'up-deltaI-genes' +'\t'+ 'down-deltaI-genes' +'\t'+ 'total-'+analysis_method+'-genes'
    title = title +'\t' + 'upregulated_genes' +'\t'+ 'downregulated_genes' +'\t'+ analysis_method+'-genes-differentially-exp'+'\t'+ 'RNA_processing/binding-factors-upregulated' +'\t'+ 'RNA_processing/binding-factors-downregulated' +'\t'+ analysis_method+'_RNA_processing/binding-factors'
    title = title +'\t'+ 'avg-downregulated-peptide-length' +'\t'+ 'std-downregulated-peptide-length' +'\t'+ 'avg-upregulated-peptide-length' +'\t'+ 'std-upregulated-peptide-length' +'\t'+ 'ttest-peptide-length' +'\t'+ 'median-peptide-length-fold-change'

    for entry in annotation_ls: title = title +'\t'+ entry
    data.write(title+'\n')
    for dataset in summary_results_db:
        values = dataset
        for entry in summary_results_db[dataset][0:-1]: values = values +'\t'+ str(entry)
        if dataset in annotation_db2:
            for entry in annotation_db2[dataset]: values = values +'\t'+ str(entry[1])
        data.write(values+'\n')
    data.close()

def compareAltAnalyzeResults(aspire_output_list,annotate_db,number_events_analyzed,analyzing_genes,analysis_method,array_type,root_dir):
    aspire_gene_db = {}; aspire_event_db = {}; event_annotation_db = {}; dataset_name_list = []
    #annotate_db[affygene] = name, symbol,ll_id,splicing_annotation

    include_all_other_genes = 'yes'    
    for filename in aspire_output_list:
        x = 0
        fn=filepath(filename)
        if '\\' in filename: names = string.split(filename,'\\')  #grab file name
        else: names = string.split(filename,'/')
        try: names = string.split(names[-1],'-'+analysis_method)
        except ValueError: print names;kill
        name = names[0]
        dataset_name_list.append(name)
        for line in open(fn,'rU').xreadlines():             
            data = cleanUpLine(line)
            data = string.split(data,'\t')  #remove endline
            y = 0
            if x == 0: x=1
            else:
                if analyzing_genes == 'no':
                    if (array_type == 'exon' or array_type == 'gene') and analysis_method in filename:
                        lowest_pvalue = float(data[8]);
                        try: si_p = float(data[20])
                        except Exception: si_p = 1
                        try: midas_p = float(data[9])
                        except ValueError: midas_p = 0
                        #print si_p,midas_p;kill
                        #if lowest_pvalue < 0.05:
                        y = 1
                        affygene = data[0]; dI = float(data[1])
                        symbol = data[2]; description = data[3]
                        exon_set1 = data[4]; exon_set2 = ''
                        event_call = data[27]; functional_attribute = data[14]
                        uniprot_attribute = data[15]; gene_expression_change = data[22]
                        dI = dI*(-1)
                    elif analysis_method in filename:
                        y = 1
                        affygene = data[0]; dI = float(data[1])
                        symbol = data[2]; description = data[3]
                        exon_set1 = data[4]+'('+data[8]+')'; exon_set2 = data[5]+'('+data[10]+')'
                        event_call = data[27]; functional_attribute = data[14]
                        uniprot_attribute = data[15]; gene_expression_change = data[22]
                        if analysis_method == 'linearregres' or analysis_method == 'ASPIRE':
                            functional_attribute = data[19]; uniprot_attribute = data[20]
                        #print exon_set1, exon_set2, data[:5];kill
                else:
                    if (array_type == 'exon' or array_type == 'gene') and analysis_method in filename:
                        y = 1
                        affygene = data[0]; dI = float(data[1])
                        symbol = data[3]; description = data[5]
                        dI_direction = data[6]; locus_link = affygene
                        exon_set1 = ''; exon_set2 = ''
                        event_call = data[-4]; functional_attribute = data[-9]
                        uniprot_attribute = data[-8]; gene_expression_change = data[-5]
                        if dI_direction == 'upregulated': dI = dI*(-1)
                    elif analysis_method in filename:
                        y = 1
                        affygene = data[0]; dI = float(data[1])
                        symbol = data[3]; description = data[5]
                        dI_direction = data[6]; locus_link = data[4]
                        exon_set1 = ''; exon_set2 = ''
                        event_call = data[-4]; functional_attribute = data[-9]
                        uniprot_attribute = data[-8]; gene_expression_change = data[-5]
                        if dI_direction == 'downregulated': dI = dI*(-1)
                        #print affygene,data[-10:];kill
                if y == 1:
                    data_tuple = [name,functional_attribute,uniprot_attribute,gene_expression_change,dI]
                    try: aspire_event_db[affygene,exon_set1,exon_set2].append(data_tuple)
                    except KeyError: aspire_event_db[affygene,exon_set1,exon_set2] = [data_tuple]
                    event_annotation_db[affygene,exon_set1,exon_set2] = event_call,symbol,description

    aspire_event_db2 = {}; splice_gene_db = {}; dataset_name_list.sort()
    for name in dataset_name_list:
        for key in event_annotation_db:
            ###record all genes in the event_annotation_db
            splice_gene_db[key[0]] = key[0]
            if key in aspire_event_db:
                x = 0
                for entry in aspire_event_db[key]:
                    if entry[0] == name:
                        x = 1
                        dI = entry[1],entry[2],entry[3],entry[4]
                        try: aspire_event_db2[key].append(dI)
                        except KeyError: aspire_event_db2[key] = [dI]
                if x ==0:
                    try: aspire_event_db2[key].append(('','','',0))
                    except KeyError: aspire_event_db2[key] = [('','','',0)]
            else:
                try: aspire_event_db2[key].append(('','','',0))
                except KeyError: aspire_event_db2[key] = [('','','',0)]

    for key in aspire_event_db2:
        dataset_size = len(aspire_event_db2[key])
        break

    ###Add all other Affygene's
    temp=[]; x = 0
    while x < dataset_size:
        temp.append(('','','',0))
        x +=1
    
    for affygene in annotate_db:
        if affygene not in splice_gene_db:
            aspire_event_db2[affygene,'',''] = temp
            
    if include_all_other_genes == 'yes': analysis_method+= '-all-genes'
    if analyzing_genes == 'no': summary_output = root_dir+'AltResults/AlternativeOutput/'+analysis_method+'-comparisons-events.txt'
    else: summary_output = root_dir+'AltResults/AlternativeOutput/'+analysis_method+'-'+ 'GENE-' +'comparisons-events.txt'
    
    fn=filepath(summary_output)
    data = open(fn,'w')
    title = 'GeneID' +'\t'+ 'symbol'+'\t'+'description' +'\t'+ 'exon_set1' +'\t'+ 'exon_set2' +'\t'+ 'event_call' +'\t'+ 'splicing_factor_call'
    for entry in dataset_name_list:
        title = title +'\t'+ entry + '-functional-attribute' +'\t'+ entry + '-uniprot-attribute' +'\t'+ entry +'-GE-change' +'\t'+ entry +'-dI'
    data.write(title +'\t'+ 'common-hits' + '\n')
    for key in aspire_event_db2:
        affygene = key[0]; exon_set1 = key[1]; exon_set2 = key[2]
        if affygene in annotate_db: splicing_factor_call = annotate_db[affygene].RNAProcessing()
        else: splicing_factor_call = ''
        try:
            event_call = event_annotation_db[key][0]
            symbol = event_annotation_db[key][1]
            description = event_annotation_db[key][2]
        except KeyError:
            event_call = ''; symbol = ''; description = ''
        values = affygene +'\t'+ symbol +'\t'+ description +'\t'+ exon_set1 +'\t'+ exon_set2 +'\t'+ event_call +'\t'+ splicing_factor_call
        x=0
        for entry in aspire_event_db2[key]:
            for info in entry: values = values +'\t'+ str(info)
            if entry[-1] != 0: x +=1
        values = values +'\t'+ str(x) + '\n'
        if include_all_other_genes == 'no':
            if x>0: data.write(values)
        else: data.write(values)
    data.close()

def exportTransitResults(array_group_list,array_raw_group_values,array_group_db,avg_const_exp_db,adj_fold_dbase,exon_db,dataset_name,apt_dir):
    """Export processed raw expression values (e.g. add global fudge factor or eliminate probe sets based on filters) to txt files
    for analysis with MiDAS"""
    #array_group_list contains group names in order of analysis
    #array_raw_group_values contains expression values for the x number of groups in above list
    #array_group_db key is the group name and values are the list of array names
    #avg_const_exp_db contains the average expression values for all arrays for all constitutive probesets, with gene as the key

    ordered_array_header_list=[]
    for group in array_group_list: ###contains the correct order for each group
        for array_id in array_group_db[group]:
            ordered_array_header_list.append(str(array_id))
    ordered_exp_val_db = {} ###new dictionary containing all expression values together, but organized based on group
    probeset_affygene_db = {} ###lists all altsplice probesets and corresponding affygenes
    for probeset in array_raw_group_values:
        try:
            include_probeset = 'yes'
            ###Examines user input parameters for inclusion of probeset types in the analysis
            if include_probeset == 'yes':
                if probeset in adj_fold_dbase: ###indicates that this probeset is analyzed for splicing (e.g. has a constitutive probeset)
                    for group_val_list in array_raw_group_values[probeset]:
                        non_log_group_exp_vals = statistics.log_fold_conversion(group_val_list)
                        for val in non_log_group_exp_vals:
                            try: ordered_exp_val_db[probeset].append(str(val))
                            except KeyError: ordered_exp_val_db[probeset] = [str(val)]
                    affygene = exon_db[probeset].GeneID()
                    try: probeset_affygene_db[affygene].append(probeset)
                    except KeyError: probeset_affygene_db[affygene] = [probeset]
        except KeyError:
            ###Indicates that the expression dataset file was not filtered for whether annotations exist in the exon annotation file
            ###In that case, just ignore the entry
            null = ''

    gene_count = 0
    ordered_gene_val_db={}
    for affygene in avg_const_exp_db: ###now, add all constitutive gene level expression values (only per anlayzed gene)
        if affygene in probeset_affygene_db: ###ensures we only include gene data where there are altsplice examined probesets
            non_log_ordered_exp_const_val = statistics.log_fold_conversion(avg_const_exp_db[affygene])
            gene_count+=1
            for val in non_log_ordered_exp_const_val:
                try: ordered_gene_val_db[affygene].append(str(val))
                except KeyError: ordered_gene_val_db[affygene] = [str(val)]

    convert_probesets_to_numbers={}
    convert_affygene_to_numbers={}; array_type = 'junction'
    probeset_affygene_number_db={}; x=0; y=0
    for affygene in probeset_affygene_db:
        x+=1; y = x  ###each affygene has a unique number, from other affygenes and probesets and probesets count up from each affygene
        x_copy = x
        example_gene_probeset = probeset_affygene_db[affygene][0]
        #if exon_db[example_gene_probeset].ArrayType() == 'exon': x_copy = exon_db[example_gene_probeset].SecondaryGeneID()
        if x_copy not in exon_db:
            convert_affygene_to_numbers[affygene] = str(x_copy)
        else: print affygene, x_copy,'new numeric for MIDAS already exists as a probeset ID number'; kill
        for probeset in probeset_affygene_db[affygene]:
            y = y+1; y_copy = y
            if exon_db[probeset].ArrayType() == 'exon':
                y_copy = probeset ### Only appropriate when the probeset ID is a number
                array_type = 'exon'
            convert_probesets_to_numbers[probeset] = str(y_copy)
            try: probeset_affygene_number_db[str(x_copy)].append(str(y_copy))
            except KeyError: probeset_affygene_number_db[str(x_copy)] = [str(y_copy)]
        x=y
    
    metafile = 'AltResults/MIDAS/meta-'+dataset_name[0:-1]+'.txt'
    data1 = export.createExportFile(metafile,'AltResults/MIDAS')
    title = 'probeset_id\ttranscript_cluster_id\tprobeset_list\tprobe_count\n'    
    data1.write(title)
    for affygene in probeset_affygene_number_db:
        probeset_list = probeset_affygene_number_db[affygene]; probe_number = str(len(probeset_list)*6)
        probeset_list = [string.join(probeset_list,' ')]
        probeset_list.append(affygene); probeset_list.append(affygene); probeset_list.reverse(); probeset_list.append(probe_number)
        probeset_list = string.join(probeset_list,'\t'); probeset_list=probeset_list+'\n'
        data1.write(probeset_list)
    data1.close()

    junction_exp_file = 'AltResults/MIDAS/'+array_type+'-exp-'+dataset_name[0:-1]+'.txt'
    fn2=filepath(junction_exp_file)
    data2 = open(fn2,'w')
    ordered_array_header_list.reverse(); ordered_array_header_list.append('probeset_id'); ordered_array_header_list.reverse()
    title = string.join(ordered_array_header_list,'\t')
    data2.write(title+'\n')
    for probeset in ordered_exp_val_db:
        probeset_number = convert_probesets_to_numbers[probeset]
        exp_values = ordered_exp_val_db[probeset]; exp_values.reverse(); exp_values.append(probeset_number); exp_values.reverse()
        exp_values = string.join(exp_values,'\t'); exp_values = exp_values +'\n'
        data2.write(exp_values)
    data2.close()

    gene_exp_file = 'AltResults/MIDAS/gene-exp-'+dataset_name[0:-1]+'.txt'
    fn3=filepath(gene_exp_file)
    data3 = open(fn3,'w')
    title = string.join(ordered_array_header_list,'\t')
    data3.write(title+'\n')
    for affygene in ordered_gene_val_db:
        try: affygene_number = convert_affygene_to_numbers[affygene]
        except KeyError: print len(convert_affygene_to_numbers), len(ordered_gene_val_db); kill
        exp_values = ordered_gene_val_db[affygene]; exp_values.reverse(); exp_values.append(affygene_number); exp_values.reverse()
        exp_values = string.join(exp_values,'\t'); exp_values = exp_values +'\n'
        data3.write(exp_values)
    data3.close()

    exportMiDASArrayNames(array_group_list,array_group_db,dataset_name,'new')
    
    coversionfile = 'AltResults/MIDAS/probeset-conversion-'+dataset_name[0:-1]+'.txt'
    fn5=filepath(coversionfile)
    data5 = open(fn5,'w')
    title = 'probeset\tprobeset_number\n'; data5.write(title)
    for probeset in convert_probesets_to_numbers: ###contains the correct order for each group
        probeset_number = convert_probesets_to_numbers[probeset]
        values = probeset+'\t'+probeset_number+'\n'
        data5.write(values)
    data5.close()

    """
    ### This code is obsolete... used before AltAnalyze could connect to APT directly.
    commands = 'AltResults/MIDAS/commands-'+dataset_name[0:-1]+'.txt'
    data = export.createExportFile(commands,'AltResults/MIDAS')
    path = filepath('AltResults/MIDAS'); path = string.replace(path,'\\','/'); path = 'cd '+path+'\n\n'
    metafile = 'meta-'+dataset_name[0:-1]+'.txt'
    junction_exp_file = array_type+'-exp-'+dataset_name[0:-1]+'.txt'
    gene_exp_file = 'gene-exp-'+dataset_name[0:-1]+'.txt'
    celfiles = 'celfiles-'+dataset_name[0:-1]+'.txt'
    command_line = 'apt-midas -c '+celfiles+' -g '+gene_exp_file+' -e '+junction_exp_file+' -m '+metafile+' -o '+dataset_name[0:-1]+'-output'
    data.write(path); data.write(command_line); data.close()
    """

    status = runMiDAS(apt_dir,array_type,dataset_name,array_group_list,array_group_db)    
    return status

def exportMiDASArrayNames(array_group_list,array_group_db,dataset_name,type):
    celfiles = 'AltResults/MIDAS/celfiles-'+dataset_name[0:-1]+'.txt'
    fn4=filepath(celfiles)
    data4 = open(fn4,'w')
    if type == 'old': cel_files = 'cel_file'
    elif type == 'new': cel_files = 'cel_files'
    title = cel_files+'\tgroup_id\n'; data4.write(title)
    for group in array_group_list: ###contains the correct order for each group
        for array_id in array_group_db[group]:
            values = str(array_id) +'\t'+ str(group) +'\n'
            data4.write(values)
    data4.close()
    
def getAPTDir(apt_fp):
    ###Determine if APT has been set to the right directory and add the analysis_type to the filename
    if 'bin' not in apt_fp: ###This directory contains the C+ programs that we wish to call
        if 'apt' in apt_fp: ###This directory is the parent directory to 'bin'
            apt_fp = apt_fp+'/'+'bin'
        elif 'Affymetrix Power Tools' in apt_fp: ###The user selected the parent directory
            dir_list = read_directory(apt_fp); versions = [] ###See what folders are in this directory (e.g., specific APT versions)
            for folder in dir_list: ###If there are multiple versions
                if 'apt' in folder:
                    version = string.replace(folder,'apt-','')
                    version = string.split(version,'.'); version_int_list = []
                    try:
                        for val in version: version_int_list.append(int(val))
                        versions.append([version_int_list,folder])
                    except Exception:
                        versions.append([folder,folder]) ### arbitrarily choose an APT version if multiple exist and the folder name does not conform to the above
            if len(versions)>0:
                ###By making the versions indexed by the integer list value of the version, we can grab the latest version
                versions.sort(); apt_fp = apt_fp+'/'+versions[-1][1]+'/'+'bin' ###Add the full path to the most recent version
    return apt_fp

def runMiDAS(apt_dir,array_type,dataset_name,array_group_list,array_group_db):
    if '/bin' in apt_dir: apt_file = apt_dir +'/apt-midas' ### if the user selects an APT directory
    elif os.name == 'nt':
        import platform
        if '32bit' in platform.architecture(): apt_file = apt_dir + '/PC/32bit/apt-midas'
        elif '64bit' in platform.architecture(): apt_file = apt_dir + '/PC/64bit/apt-midas' 
    elif 'darwin' in sys.platform: apt_file = apt_dir + '/Mac/apt-midas'
    elif 'linux' in sys.platform:
        import platform
        if '32bit' in platform.architecture(): apt_file = apt_dir + '/Linux/32bit/apt-midas'
        elif '64bit' in platform.architecture(): apt_file = apt_dir + '/Linux/64bit/apt-midas' 
    apt_file = filepath(apt_file)

    ### Each input file for MiDAS requires a full file path, so get the parent path
    midas_input_dir = 'AltResults/MIDAS/'
    path=filepath(midas_input_dir)

    ### Remotely connect to the previously verified APT C+ midas program and run analysis
    metafile = path + 'meta-'+dataset_name[0:-1]+'.txt'
    exon_or_junction_file = path + array_type+'-exp-'+dataset_name[0:-1]+'.txt'
    gene_exp_file = path + 'gene-exp-'+dataset_name[0:-1]+'.txt'
    celfiles = path + 'celfiles-'+dataset_name[0:-1]+'.txt'
    output_file = path + dataset_name[0:-1]+'-output'

    ### Delete the output folder if it already exists (may cause APT problems)
    delete_status = export.deleteFolder(output_file)
    try:
        import subprocess
        retcode = subprocess.call([
        apt_file, "--cel-files", celfiles,"-g", gene_exp_file, "-e", exon_or_junction_file,
        "-m", metafile, "-o", output_file])
        if retcode: status = 'failed'
        else: status = 'run'
    except NameError: status = 'failed'
    if status == 'failed':
        try:
            ### Try running the analysis with old MiDAS file headers and command
            exportMiDASArrayNames(array_group_list,array_group_db,dataset_name,'old')
            import subprocess
            retcode = subprocess.call([
            apt_file, "-c", celfiles,"-g", gene_exp_file, "-e", exon_or_junction_file,
            "-m", metafile, "-o", output_file])
            if retcode: status = 'failed'
            else: status = 'run'        
        except Exception: status = 'failed'
        if status == 'failed': print "apt-midas failed"
    else: print "apt-midas run successfully"
    return status

def importMidasOutput(dataset_name):
    coversionfile = 'AltResults/MIDAS/probeset-conversion-'+dataset_name[0:-1]+'.txt'
    #print "Looking for", coversionfile
 
    fn=filepath(coversionfile); x=0; probeset_conversion_db={}
    for line in open(fn,'rU').xreadlines():         
        data = cleanUpLine(line)
        if x==0: x=1 ###Occurs for the header line
        else:
            probeset,probeset_number = string.split(data,'\t')
            probeset_conversion_db[probeset_number] = probeset
            
    midas_results = 'AltResults/MIDAS/'+dataset_name[:-1]+'-output'+'/midas.pvalues.txt'
    fn=filepath(midas_results); x=0; midas_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line) 
        if data[0] == '#': continue
        elif x==0: x=1 ###Occurs for the header line
        else:
            t = string.split(data,'\t')
            try: probeset_number,geneid,p = t
            except ValueError: print t;kill
            try: p = float(p)
            except ValueError: p = 1.000 ### "-1.#IND" can occur, when the constituitive and probeset are the same
            probeset = probeset_conversion_db[probeset_number]
            midas_db[probeset] = p
    return midas_db

def combineRawSpliceResults(species,analysis_method):
    import_dir = '/AltResults/RawSpliceData/'+species+'/'+analysis_method
    g = GrabFiles(); g.setdirectory(import_dir)
    files_to_merge = g.searchdirectory('combined') ###made this a term to excluded

    headers =[]; combined_data={}
    for filename in files_to_merge:
        fn=filepath(filename); x=0
        for line in open(fn,'rU').xreadlines():         
            data = cleanUpLine(line) 
            t = string.split(data,'\t')
            if x==0:
                headers += t[1:]
                x=1 ###Occurs for the header line
            else:
                values = t; values = values[1:]; key = t[0]
                try: combined_data[key]+=values
                except KeyError: combined_data[key]=values

    max_len=0 ###Some files will contain data that others won't... normalize for this, so we only include raws where there is data in all files examined
    for key in combined_data:
        if len(combined_data[key])>max_len: max_len = len(combined_data[key])

    combined_data2 = {}; k=0; j=0
    for key in combined_data:
        #print combined_data[key];kill
        ### '1' in the list, then there was only one constitutive probeset that was 'expressed' in that dataset: thus comparisons are not likely valid
        count = list.count(combined_data[key],'1.0')
        if len(combined_data[key])==max_len and count <3: combined_data2[key] = combined_data[key]
        elif len(combined_data[key])!=max_len: k+=1#; print key,max_len, len(combined_data[key]),combined_data[key]; kill
        elif count >2: j+=1
    combined_data = combined_data2

    #print k,j
    export_file = import_dir[1:]+'/combined.txt'
    fn=filepath(export_file);data = open(fn,'w')
    title = string.join(['gene-probeset']+headers,'\t')+'\n'; data.write(title)
    for key in combined_data:
        values = string.join([key]+combined_data[key],'\t')+'\n'; data.write(values)
    data.close()
    print "exported",len(combined_data),"to",export_file

def import_annotations(filename,array_type):
    from build_scripts import ExonAnalyze_module
    fn=filepath(filename); annotate_db = {}; x = 0
    if array_type == 'AltMouse':
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line) 
            if x == 0: x = 1
            else:
                try: affygene, description, ll_id, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError: affygene, description, ll_id, symbol = string.split(data,'\t'); splicing_annotation = ''
                if '"' in description: null,description,null = string.split(description,'"')
                rna_processing_annot =''
                y = ExonAnalyze_module.GeneAnnotationData(affygene, description, symbol, ll_id, rna_processing_annot)
                annotate_db[affygene] = y
    else:
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line) 
            if x == 0: x = 1
            else:
                rna_processing_annot=''
                try: ensembl, description, symbol, rna_processing_annot = string.split(data,'\t')
                except ValueError: ensembl, description, symbol = string.split(data,'\t')
                y = ExonAnalyze_module.GeneAnnotationData(ensembl, description, symbol, ensembl, rna_processing_annot)
                annotate_db[ensembl] = y
    return annotate_db

if __name__ == '__main__':
    array_type = 'exon'
    a = 'Mm'; b = 'Hs'
    e = 'ASPIRE'; f = 'linearregres'; g = 'ANOVA'; h = 'splicing-index'
    analysis_method = h
    species = b ### edit this

    if array_type != 'AltMouse': gene_annotation_file = "AltDatabase/ensembl/"+species+"/"+species+"_Ensembl-annotations.txt"
    annotate_db = import_annotations(gene_annotation_file,array_type)
    number_events_analyzed = 0
    analyzing_genes = 'no'
    root_dir = 'C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/Combined-GSE14588_RAW/junction/' ### edit this
    a = root_dir+'AltResults/AlternativeOutput/Hs_Junction_d14_vs_d7.p5_average-ASPIRE-exon-inclusion-results.txt' ### edit this
    b = root_dir+'AltResults/AlternativeOutput/Hs_Junction_d14_vs_d7.p5_average-splicing-index-exon-inclusion-results.txt' ### edit this
    aspire_output_list = [a,b]
    compareAltAnalyzeResults(aspire_output_list,annotate_db,number_events_analyzed,analyzing_genes,analysis_method,array_type,root_dir)

    #dataset_name = 'test.'; apt_dir = 'AltDatabase/affymetrix/APT'
    #aspire_output_gene_list = ['AltResults/AlternativeOutput/Hs_Exon_CS-d40_vs_hESC-d0.p5_average-splicng_index-exon-inclusion-GENE-results.txt', 'AltResults/AlternativeOutput/Hs_Exon_Cyt-NP_vs_Cyt-ES.p5_average-splicing_index-exon-inclusion-GENE-results.txt', 'AltResults/AlternativeOutput/Hs_Exon_HUES6-NP_vs_HUES6-ES.p5_average-splicing_index-exon-inclusion-GENE-results.txt']
    #runMiDAS(apt_dir,array_type,dataset_name,{},()); sys.exit()
    #midas_db = importMidasOutput(dataset_name)