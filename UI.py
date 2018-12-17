###UI
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

import math
from stats_scripts import statistics
import sys, string
import shutil
import os.path
import unique
import update; reload(update)
import export
import ExpressionBuilder
import time
import webbrowser
import traceback
import AltAnalyze
from sys import argv

"""
import numpy
import scipy
from PIL import Image as PIL_Image
import ImageTk
import matplotlib
import matplotlib.pyplot as pylab
"""

try:
    try:
        from visualization_scripts import WikiPathways_webservice
    except Exception:
        #print traceback.format_exc()
        if 'URLError' in traceback.format_exc():
            print 'No internet connection found'
        else:
            print 'WikiPathways visualization not supported (requires installation of suds)'
    try:
        from PIL import Image as PIL_Image
        try: import ImageTk
        except Exception: from PIL import ImageTk
        import PIL._imaging
        import PIL._imagingft
    except Exception:
        print traceback.format_exc()
        #print 'Python Imaging Library not installed... using default PNG viewer'
        None

    try:
        ### Only used to test if matplotlib is installed
        #import matplotlib
        #import matplotlib.pyplot as pylab
        None
    except Exception:
        #print traceback.format_exc()
        print 'Graphical output mode disabled (requires matplotlib, numpy and scipy)'
        None

except Exception:
    None

command_args = string.join(sys.argv,' ')
if len(sys.argv[1:])>1 and '-' in command_args and '--GUI' not in command_args:
    runningCommandLine = True
else:
    runningCommandLine = False
    try:
        import Tkinter 
        #import bwidget; from bwidget import *
        from Tkinter import *
        from visualization_scripts import PmwFreeze
        from Tkconstants import LEFT
        import tkMessageBox
        import tkFileDialog
    except Exception: print "\nPmw or Tkinter not found... proceeding with manual input"
    
mac_print_mode = 'no' 
if os.name == 'posix':  mac_print_mode = 'yes' #os.name is  'posix', 'nt', 'os2', 'mac', 'ce' or 'riscos'
debug_mode = 'no'
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def osfilepath(filename):
    fn = filepath(filename)
    fn = string.replace(fn,'\\','/')
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def getFolders(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        if entry[-4:] != ".txt" and entry[-4:] != ".csv" and ".zip" not in entry: dir_list2.append(entry)
    return dir_list2

def returnDirectoriesNoReplace(dir):
    dir_list = unique.returnDirectoriesNoReplace(dir); dir_list2 = []
    for entry in dir_list:
        if '.' not in entry and 'affymetrix' not in entry:
            if 'EnsMart' in entry: dir_list2.append(entry)
    return dir_list2

def returnFilesNoReplace(dir):
    dir_list = unique.returnDirectoriesNoReplace(dir); dir_list2 = []
    for entry in dir_list:
        if '.' in entry: dir_list2.append(entry)
    return dir_list2
    
def identifyCELfiles(dir,array_type,vendor):
    dir_list = read_directory(dir); dir_list2=[]; full_dir_list=[]
    datatype = 'arrays'
    types={}

    for file in dir_list:
        original_file = file
        file_lower = string.lower(file); proceed = 'no'
        ### "._" indicates a mac alias
        if ('.cel' in file_lower[-4:] and '.cel.' not in file_lower) and file_lower[:2] != '._':
            proceed = 'yes'
        elif ('.bed' in file_lower[-4:] or '.tab' in file_lower or '.junction_quantification.txt' in file_lower or '.bam' in file_lower) and file_lower[:2] != '._' and '.bai' not in file_lower:
            proceed = 'yes'
            datatype = 'RNASeq'
        elif array_type == "3'array" and '.cel' not in  file_lower[-4:] and '.txt' in file_lower[-4:] and vendor != 'Affymetrix':
            proceed = 'yes'
        if proceed == 'yes':
            if '__' in file and '.cel' not in string.lower(file):
                #print file,string.split(file,'__'),file[-4:]
                file=string.split(file,'__')[0]+file[-4:]
                if '.tab' in original_file: file = string.replace(file,'.txt','.tab')
                elif '.bed' in original_file: file = string.replace(file,'.txt','.bed')
                if '.TAB' in original_file: file = string.replace(file,'.txt','.TAB')
                elif '.BED' in original_file: file = string.replace(file,'.txt','.BED')
            dir_list2.append(file)
            file = dir+'/'+file
            full_dir_list.append(file)
    dir_list2 = unique.unique(dir_list2)
    full_dir_list = unique.unique(full_dir_list)
    dir_list2.sort(); full_dir_list.sort()
    
    if datatype == 'RNASeq':
        checkBEDFileFormat(dir) ### Make sure the names are wonky
        dir_list3=[]
        c = string.lower(string.join(dir_list2,''))
        if '.bam' in c and '.bed' in c: #If bed present use bed and not bam
            for i in dir_list2:
                if '.bam' not in i:
                    dir_list3.append(i)
            dir_list2 = dir_list3
        elif '.bam' in c:
            for i in dir_list2:
                if '.bam' in i:
                    dir_list3.append(string.replace(i,'.bam','.bed'))
                elif '.BAM' in i:
                    dir_list3.append(string.replace(i,'.BAM','.bed'))
            dir_list2 = dir_list3         
        
    return dir_list2,full_dir_list

def checkBEDFileFormat(bed_dir):
    """ This checks to see if some files have two underscores and one has none or if double underscores are missing from all."""
    dir_list = read_directory(bed_dir)
    condition_db={}
    for filename in dir_list:
        if '.tab' in string.lower(filename) or '.bed' in string.lower(filename) or '.junction_quantification.txt' in string.lower(filename):
            condition_db[filename]=[]

    if len(condition_db)==0: ### Occurs if BAMs present but not .bed files
        for filename in dir_list:
            if '.bam' in string.lower(filename):
                condition_db[filename]=[]
                 
    ### Check to see if exon.bed and junction.bed file names are propper or faulty (which will result in downstream errors)
    double_underscores=[]
    no_doubles=[]
    for condition in condition_db:
        if '__' in condition:
            double_underscores.append(condition)
        else:
            no_doubles.append(condition)
    
    exon_beds=[]
    junctions_beds=[] 
    if len(double_underscores)>0 and len(no_doubles)>0:
        ### Hence, a problem is likely due to inconsistent naming
        print_out = 'The input files appear to have inconsistent naming. If both exon and\njunction sample data are present, make sure they are named propperly.\n\n'
        print_out += 'For example: cancer1__exon.bed, cancer1__junction.bed\n(double underscore required to match these samples up)!\n\n'
        print_out += 'Exiting AltAnalyze'
        IndicatorWindowSimple(print_out,'Quit')
        sys.exit()
    elif len(no_doubles)>0:
        for condition in no_doubles:
            condition = string.lower(condition)
            if 'exon' in condition:
                exon_beds.append(condition)
            if 'junction' in condition:
                junctions_beds.append(condition)
        if len(exon_beds)>0 and len(junctions_beds)>0:
            print_out = 'The input files appear to have inconsistent naming. If both exon and\njunction sample data are present, make sure they are named propperly.\n\n'
            print_out += 'For example: cancer1__exon.bed, cancer1__junction.bed\n(double underscore required to match these samples up)!\n\n'
            print_out += 'Exiting AltAnalyze'
            IndicatorWindowSimple(print_out,'Quit')
            sys.exit()
            
def identifyArrayType(full_dir_list):
    #import re
    arrays={}; array_type=None ### Determine the type of unique arrays in each directory
    for filename in full_dir_list:
        fn=filepath(filename); ln=0
        for line in open(fn,'rU').xreadlines():
            if '\x00' in line: ### Simple way of determining if it is a version 4 file with encoding
                line = string.replace(line,'\x00\x00',' ') ### retains spaces
                line = string.replace(line,'\x00','') ### returns human readable line
            if ln<150:
                data = cleanUpLine(line); ln+=1
                if 'sq' in data:
                    try:
                        #fileencoding = "iso-8859-1"
                        #txt = line.decode(fileencoding); print [txt];kill ### This works but so does the above
                        array_info,null = string.split(data,'sq')
                        array_info = string.split(array_info,' ')
                        array_type = array_info[-1]
                        if '.' in array_type: array_type,null = string.split(array_type,'.')
                        #array_type = string.join(re.findall(r"\w",array_type),'') ### should force only alphanumeric but doesn't seem to always work
                        arrays[array_type]=[]
                        #print array_type+'\t'+filename
                        break
                    except Exception: pass
                elif 'affymetrix-array-type' in data:
                    null, array_type = string.split(data,'affymetrix-array-type')
                    if '.' in array_type: array_type,null = string.split(array_type,'.')
                    arrays[array_type]=[]
                """else: ### some CEL file versions are encoded
                    fileencoding = "iso-8859-1"
                    txt = line.decode(fileencoding)
                    print txt;kill"""
            else: break
    
    array_ls = []
    for array in arrays:
        if len(array)<50: array_ls.append(array) ### Occurs with version 4 encoding (bad entries added)
    return array_ls, array_type

def getAffyFilesRemote(array_name,arraytype,species):
    global backSelect; global array_type; global debug_mode
    debug_mode = 'yes'
    backSelect = 'yes'
    array_type = arraytype
    library_dir, annotation_dir, bgp_file, clf_file = getAffyFiles(array_name,species)
    return library_dir, annotation_dir, bgp_file, clf_file
    
def getAffyFiles(array_name,species):#('AltDatabase/affymetrix/LibraryFiles/'+library_file,species)
    sa = supproted_array_db[array_name]; library_file = sa.LibraryFile(); annot_file = sa.AnnotationFile(); original_library_file = library_file
    filename = 'AltDatabase/affymetrix/LibraryFiles/'+library_file
    fn=filepath(filename); library_dir=filename; bgp_file = ''; clf_file = ''
    local_lib_files_present = False
    if backSelect == 'yes': warn = 'no'
    else: warn = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
        ### Hence, the library file was found!!!
        local_lib_files_present = True
        input_cdf_file = filename
        if '.pgf' in input_cdf_file:
            ###Check to see if the clf and bgp files are present in this directory 
            icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
            clf_short = string.replace(cdf_short,'.pgf','.clf')
            if array_type == 'exon' or array_type == 'junction':
                bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
            else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
            try: dir_list = read_directory(parent_dir)
            except Exception: dir_list = read_directory('/'+parent_dir)
            if clf_short in dir_list and bgp_short in dir_list:
                pgf_file = input_cdf_file; clf_file = string.replace(pgf_file,'.pgf','.clf')
                if array_type == 'exon' or array_type == 'junction': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
            else:
                try:
                    print_out = "The directory;\n"+parent_dir+"\ndoes not contain either a .clf or antigenomic.bgp\nfile, required for probeset summarization."
                    IndicatorWindow(print_out,'Continue')
                except Exception: print print_out; sys.exit()
                        
    except Exception:
        print_out = "AltAnalyze was not able to find a library file\nfor your arrays. Would you like AltAnalyze to\nautomatically download these files?"
        try:
            dw = DownloadWindow(print_out,'Download by AltAnalyze','Select Local Files')
            warn = 'no' ### If already downloading the library, don't warn to download the csv too
            dw_results = dw.Results(); option = dw_results['selected_option']
        except Exception: option = 1 ### Occurs when Tkinter is not present - used by CommandLine mode
        if option == 1:
            library_file = string.replace(library_file,'.cdf','.zip')
            filename = 'AltDatabase/affymetrix/LibraryFiles/'+library_file
            input_cdf_file = filename
            if '.pgf' in input_cdf_file:
                pgf_file = input_cdf_file; clf_file = string.replace(pgf_file,'.pgf','.clf')
                if array_type == 'exon' or array_type == 'junction': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                filenames = [pgf_file+'.gz',clf_file+'.gz',bgp_file+'.gz']
                if 'Glue' in pgf_file:
                    kil_file = string.replace(pgf_file,'.pgf','.kil') ### Only applies to the Glue array
                    filenames.append(kil_file+'.gz')
            else: filenames = [input_cdf_file]
            for filename in filenames:
                var_list = filename,'LibraryFiles'
                if debug_mode == 'no': StatusWindow(var_list,'download')
                else:
                    for filename in filenames:
                        continue_analysis = update.downloadCurrentVersion(filename,'LibraryFiles','')
                try: os.remove(filepath(filename)) ### Not sure why this works now and not before
                except Exception: pass
        else: library_dir = ''
    
    filename = 'AltDatabase/affymetrix/'+species+'/'+annot_file
    fn=filepath(filename); annotation_dir = filename
    
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception:
        if warn == 'yes' and local_lib_files_present == False:
            ### Indicates that library file wasn't present to prior to this method
            print_out = "AltAnalyze was not able to find a CSV annotation file\nfor your arrays. Would you like AltAnalyze to\nautomatically download these files?"
            try:
                dw = DownloadWindow(print_out,'Download by AltAnalyze','Select Local Files'); warn = 'no'
                dw_results = dw.Results(); option = dw_results['selected_option']
            except OSError: option = 1 ### Occurs when Tkinter is not present - used by CommandLine mode
        else:
            try: option = option
            except Exception: option = 2
        if option == 1 or debug_mode=='yes':
            annot_file += '.zip'
            filenames = ['AltDatabase/affymetrix/'+species+'/'+annot_file]
            for filename in filenames:
                var_list = filename,'AnnotationFiles'
                if debug_mode == 'no': StatusWindow(var_list,'download')
                else:
                    for filename in filenames:
                        try: update.downloadCurrentVersionUI(filename,'AnnotationFiles','',Tk())
                        except Exception:
                            try: update.downloadCurrentVersion(filename,'AnnotationFiles',None)
                            except Exception: pass ### Don't actually need Affy's annotations in most cases - GO-Elite used instead
                try: os.remove(filepath(filename)) ### Not sure why this works now and not before
                except Exception: pass
        else: annotation_dir = ''
    return library_dir, annotation_dir, bgp_file, clf_file

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

########### Status Window Functions ###########
def copyFiles(file1,file2,root):
    print 'Copying files from:\n',file1
    data = export.ExportFile(file2) ### Ensures the directory exists
    try: shutil.copyfile(file1,file2)
    except Exception: print "This file already exists in the destination directory."
    root.destroy()
    
class StatusWindow:
    def __init__(self,info_list,analysis_type,windowType='parent'):
        try:
            if windowType == 'child':
                root = Toplevel()
            else:
                root = Tk()
            self._parent = root
            root.title('AltAnalyze version 2.1.2')
            statusVar = StringVar() ### Class method for Tkinter. Description: "Value holder for strings variables."

            height = 300; width = 700
            if os.name != 'nt': height+=100; width+=50
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = 'Download File Status Window',
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            
            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Output')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                
            Label(group.interior(),width=180,height=1000,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)

            status = StringVarFile(statusVar,root) ### Captures the stdout (or print) to the GUI instead of to the terminal
            self.original_sys_out = sys.stdout ### Save the original stdout mechanism
            #ProgressBar('Download',self._parent)
        except Exception: None
        if analysis_type == 'download':
            filename,dir = info_list
            try: sys.stdout = status; root.after(100,update.downloadCurrentVersionUI(filename,dir,None,self._parent))
            except Exception:
                update.downloadCurrentVersion(filename,dir,None)
        if analysis_type == 'copy':
            file1,file2 = info_list
            try: sys.stdout = status; root.after(100,copyFiles(file1,file2,self._parent))
            except Exception: copyFiles(file1,file2,None)
        if analysis_type == 'getOnlineDBConfig':
            file_location_defaults = info_list
            try: sys.stdout = status; root.after(100,getOnlineDBConfig(file_location_defaults,self._parent))
            except Exception,e: getOnlineDBConfig(file_location_defaults,None)
        if analysis_type == 'getOnlineEliteDatabase':
            file_location_defaults,db_version,new_species_codes,update_goelite_resources = info_list
            try: sys.stdout = status; root.after(100,getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,update_goelite_resources,self._parent))
            except Exception,e: getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,update_goelite_resources,None)
        if analysis_type == 'getAdditionalOnlineResources':
            species_code,additional_resources = info_list
            try: sys.stdout = status; root.after(100,getAdditionalOnlineResources(species_code,additional_resources,self._parent))
            except Exception,e: getAdditionalOnlineResources(species_code,additional_resources,None)
        if analysis_type == 'createHeatMap':
            filename, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast = info_list
            try: sys.stdout = status; root.after(100,createHeatMap(filename, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast, self._parent))
            except Exception,e: createHeatMap(filename, row_method, row_metric, column_method, column_metric, color_gradient, transpose,contrast,None)
        if analysis_type == 'performPCA':
            filename, pca_labels, dimensions, pca_algorithm, transpose, geneSetName, species, zscore, colorByGene, reimportModelScores = info_list
            try: sys.stdout = status; root.after(100,performPCA(filename, pca_labels, pca_algorithm, transpose, self._parent, plotType = dimensions, geneSetName=geneSetName, species=species, zscore=zscore, colorByGene=colorByGene, reimportModelScores=reimportModelScores))
            except Exception,e: performPCA(filename, pca_labels, pca_algorithm, transpose, None, plotType = dimensions, geneSetName=geneSetName, species=species, zscore=zscore, colorByGene=colorByGene, reimportModelScores=reimportModelScores)
        if analysis_type == 'runLineageProfiler':
            fl, filename, vendor, custom_markerFinder, geneModel_file, modelDiscovery = info_list
            try: sys.stdout = status; root.after(100,runLineageProfiler(fl, filename, vendor, custom_markerFinder, geneModel_file, self._parent, modelSize=modelDiscovery))
            except Exception,e: runLineageProfiler(fl, filename, vendor, custom_markerFinder, geneModel_file, None, modelSize=modelDiscovery)
        if analysis_type == 'MergeFiles':
            files_to_merge, join_option, ID_option, output_merge_dir = info_list
            try: sys.stdout = status; root.after(100,MergeFiles(files_to_merge, join_option, ID_option, output_merge_dir, self._parent))
            except Exception,e: MergeFiles(files_to_merge, join_option, ID_option, output_merge_dir, None)
        if analysis_type == 'VennDiagram':
            files_to_merge, output_venn_dir = info_list
            try: sys.stdout = status; root.after(100,vennDiagram(files_to_merge, output_venn_dir, self._parent))
            except Exception,e: vennDiagram(files_to_merge, output_venn_dir, None)
        if analysis_type == 'AltExonViewer':
            species,platform,exp_file,gene,show_introns,analysisType = info_list
            try: sys.stdout = status; root.after(100,altExonViewer(species,platform,exp_file,gene,show_introns,analysisType,self._parent))
            except Exception,e: altExonViewer(species,platform,exp_file,gene,show_introns,analysisType,None)
        if analysis_type == 'network':
            inputDir,inputType,outputdir,interactionDirs,degrees,expressionFile,gsp = info_list
            try: sys.stdout = status; root.after(100,networkBuilder(inputDir,inputType,outputdir,interactionDirs,degrees,expressionFile,gsp, self._parent))
            except Exception,e: networkBuilder(inputDir,inputType,outputdir,interactionDirs,degrees,expressionFile,gsp, None)      
        if analysis_type == 'IDConverter':
            filename, species_code, input_source, output_source = info_list
            try: sys.stdout = status; root.after(100,IDconverter(filename, species_code, input_source, output_source, self._parent))
            except Exception,e: IDconverter(filename, species_code, input_source, output_source, None)
        if analysis_type == 'predictGroups':
            try: expFile, mlp_instance, gsp, reportOnly = info_list
            except Exception: expFile, mlp_instance, gsp, reportOnly = info_list
            try: sys.stdout = status; root.after(100,predictSampleExpGroups(expFile, mlp_instance, gsp, reportOnly, self._parent))
            except Exception,e: predictSampleExpGroups(expFile, mlp_instance, gsp, reportOnly, None)
        if analysis_type == 'preProcessRNASeq':
            species,exp_file_location_db,dataset,mlp_instance  = info_list
            try: sys.stdout = status; root.after(100,preProcessRNASeq(species,exp_file_location_db,dataset,mlp_instance, self._parent))
            except Exception,e: preProcessRNASeq(species,exp_file_location_db,dataset,mlp_instance, None)

        try:
            self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
            self._parent.mainloop()
            self._parent.destroy()
        except Exception: None  ### This is what typically get's called
        try:
            sys.stdout = self.original_sys_out ### Has to be last to work!!!
        except Exception: None

    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()
    def quit(self):
        try: self._parent.destroy(); sys.exit() #self._parent.quit(); 
        except Exception: sys.exit() #self._parent.quit();
    def SysOut(self):
        return self.original_sys_out
    
def preProcessRNASeq(species,exp_file_location_db,dataset,mlp_instance,root):
    for dataset in exp_file_location_db:
        flx = exp_file_location_db[dataset]
    if root == None: display=False
    else: display=True
    runKallisto = False
    try:
        import RNASeq
        from build_scripts import ExonArray
        expFile = flx.ExpFile()
        count = verifyFileLength(expFile)
        try: fastq_folder = flx.RunKallisto()
        except Exception: fastq_folder = []
        try: customFASTA = flx.CustomFASTA()
        except Exception: customFASTA = None
        try: matrix_file = flx.ChromiumSparseMatrix()
        except Exception: matrix_file = []
        if len(matrix_file)>0:
            print 'Exporting Chromium sparse matrix file to tab-delimited-text'
            try:
                #print expFile, 'expFile'
                #print flx.RootDir(), 'root_dir'
                output_dir = export.findParentDir(expFile)
                try: os.mkdir(output_dir)
                except Exception: pass
                matrix_dir = export.findParentDir(matrix_file)
                genome = export.findFilename(matrix_dir[:-1])
                parent_dir = export.findParentDir(matrix_dir[:-1])
                from import_scripts import ChromiumProcessing
                ChromiumProcessing.import10XSparseMatrix(parent_dir,genome,dataset,expFile=expFile)
            except Exception:
                print 'Chromium export failed due to:',traceback.format_exc()
            try: root.destroy()
            except Exception: pass
            return None
        elif len(fastq_folder)>0 and count<2:
            print 'Pre-processing input files'
            try:
                parent_dir = export.findParentDir(expFile)
                flx.setRootDir(parent_dir)
                RNASeq.runKallisto(species,dataset,flx.RootDir(),fastq_folder,mlp_instance,returnSampleNames=False,customFASTA=customFASTA)   
            except Exception:
                print 'Kallisto failed due to:',traceback.format_exc()
            try: root.destroy()
            except Exception: pass
            return None
        elif len(fastq_folder)>0 and count>1:
            try: root.destroy()
            except Exception: pass
            return None ### Already run
        elif count<2:
            print 'Pre-processing input BED/BAM files\n'
            analyzeBAMs=False
            bedFilesPresent=False
            dir_list = unique.read_directory(flx.BEDFileDir())
            for file in dir_list:
                if '.bam' in string.lower(file):
                    analyzeBAMs=True
                if '.bed' in string.lower(file):
                    bedFilesPresent=True
            if analyzeBAMs and bedFilesPresent==False:
                print 'No bed files present, deriving from BAM files'
                from import_scripts import multiBAMtoBED
                bam_dir = flx.BEDFileDir()
                refExonCoordinateFile = filepath('AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt')
                outputExonCoordinateRefBEDfile = bam_dir+'/BedRef/'+species+'_'+string.replace(dataset,'exp.','')
                analysisType = ['exon','junction','reference']
                #analysisType = ['junction']
                try: multiBAMtoBED.parallelBAMProcessing(bam_dir,refExonCoordinateFile,outputExonCoordinateRefBEDfile,analysisType=analysisType,useMultiProcessing=flx.multiThreading(),MLP=mlp_instance,root=root)
                except Exception:
                    print traceback.format_exc()
            try: biotypes = RNASeq.alignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset,Multi=mlp_instance)
            except Exception:
                print traceback.format_exc()
            biotypes = getBiotypes(expFile)
        else:
            biotypes = getBiotypes(expFile)
        
        array_linker_db,array_names = ExonArray.remoteExonProbesetData(expFile,{},'arraynames',flx.ArrayType())
        steady_state_export = expFile[:-4]+'-steady-state.txt'
        normalize_feature_exp = flx.FeatureNormalization()
        try: excludeLowExpressionExons = flx.excludeLowExpressionExons()
        except Exception: excludeLowExpressionExons = True
 
        if flx.useJunctionsForGeneExpression():
            if 'junction' in biotypes:
                feature = 'junction'
            else:
                feature = 'exon'
        else:
            ### Use all exons either way at this step since more specific parameters will apply to the next iteration
            if 'exon' in biotypes:
                feature = 'exon'
            else:
                feature = 'junction'
        probeset_db = getAllKnownFeatures(feature,species,flx.ArrayType(),flx.Vendor(),flx)
        print 'Calculating gene-level expression values from',feature+'s'
        RNASeq.calculateGeneLevelStatistics(steady_state_export,species,probeset_db,normalize_feature_exp,array_names,flx,excludeLowExp=excludeLowExpressionExons,exportRPKMs=True)
        
        #if display == False: print print_out
        #try: InfoWindow(print_out, 'Continue')
        #except Exception: None
        try: root.destroy()
        except Exception: pass
    except Exception:
        error = traceback.format_exc()
        #print error
        try:
            logfile = filepath(fl.RootDir()+'Error.log')
            log_report = open(logfile,'a')
            log_report.write(traceback.format_exc())
        except Exception:
            None
        print_out = 'Expression quantification failed..\n',error
        if runningCommandLine==False:
            try: print print_out
            except Exception: pass ### Windows issue with the Tk status window stalling after pylab.show is called
            try: WarningWindow(print_out,'Continue')
            except Exception: pass
        try: root.destroy()
        except Exception: pass

def getBiotypes(filename):
    biotypes={}
    firstRow=True
    if 'RawSpliceData' in filename: index = 2
    else: index = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            t = string.split(line,'\t')
            if firstRow:
                firstRow = False
            else:
                if '-' in t[index]:
                    biotypes['junction']=[]
                else:
                    biotypes['exon']=[]
    except Exception: pass
    return biotypes

def getAllKnownFeatures(feature,species,array_type,vendor,fl):
    ### Simple method to extract gene features of interest
    from build_scripts import ExonArrayEnsemblRules

    source_biotype = 'mRNA'
    if array_type == 'gene': source_biotype = 'gene'
    elif array_type == 'junction': source_biotype = 'junction'

    if array_type == 'AltMouse':
        import ExpressionBuilder
        probeset_db,constitutive_gene_db = ExpressionBuilder.importAltMerge('full')
        source_biotype = 'AltMouse'
    elif vendor == 'Affymetrix' or array_type == 'RNASeq':
        if array_type == 'RNASeq':
            source_biotype = array_type, fl.RootDir()
            
    dbs = ExonArrayEnsemblRules.getAnnotations('no','Ensembl',source_biotype,species)
    probeset_db = dbs[0]; del dbs
    probeset_gene_db={}
    for probeset in probeset_db:
        probe_data = probeset_db[probeset]
        gene = probe_data[0]; external_exonid = probe_data[-2]
        if len(external_exonid)>2: ### These are known exon only (e.g., 'E' probesets)
            proceed = True
            if feature == 'exon': ### Restrict the analysis to exon RPKM or count data for constitutive calculation
                if '-' in probeset and '_' not in probeset: proceed = False
            else:
                if '-' not in probeset and '_' not in probeset: proceed = False ### Use this option to override 
            if proceed:
                try: probeset_gene_db[gene].append(probeset)
                except Exception: probeset_gene_db[gene] = [probeset]
    return probeset_gene_db
                        
def RemotePredictSampleExpGroups(expFile, mlp_instance, gsp, globalVars):
    global species
    global array_type
    species, array_type = globalVars
    predictSampleExpGroups(expFile, mlp_instance, gsp, False, None, exportAdditionalResults=False)
    return graphic_links

def getICGSNMFOutput(root,filename):
    if 'exp.' in filename:
        filename = string.replace(filename,'exp.','')
    umap_scores_file=''
    for file in unique.read_directory(root):
        if '-UMAP_coordinates.txt' in file and filename in file:
            umap_scores_file = file
    return root+'/'+umap_scores_file

def exportAdditionalICGSOutputs(expFile,group_selected,outputTSNE=True):
    
    ### Remove OutlierRemoved files (they will otherwise clutter up the directory)
    if 'OutliersRemoved' in group_selected and 'OutliersRemoved' not in expFile:
        try: os.remove(expFile[:-4]+'-OutliersRemoved.txt')
        except Exception: pass
        try: os.remove(string.replace(expFile[:-4]+'-OutliersRemoved.txt','exp.','groups.'))
        except Exception: pass
        try: os.remove(string.replace(expFile[:-4]+'-OutliersRemoved.txt','exp.','comps.'))
        except Exception: pass
        
    ### Create the new groups file but don't over-write the old
    import RNASeq
    new_groups_dir = RNASeq.exportGroupsFromClusters(group_selected,expFile,array_type,suffix='ICGS')
    from import_scripts import sampleIndexSelection

    ### Look to see if a UMAP file exists in ICGS-NMF
    if 'ICGS-NMF' in group_selected or 'NMF-SVM' in group_selected:
        root = export.findParentDir(export.findParentDir(expFile)[:-1])+'/ICGS-NMF'
        umap_scores_file = getICGSNMFOutput(root,export.findFilename(expFile)[:-4])
        tSNE_score_file = umap_scores_file
        extension = '-UMAP_scores.txt'
    elif outputTSNE:
        try:
            ### Build-tSNE plot from the selected ICGS output (maybe different than Guide-3)
        
            tSNE_graphical_links = performPCA(group_selected, 'no', 'UMAP', False, None, plotType='2D',
                display=False, geneSetName=None, species=species, zscore=True, reimportModelScores=False,
                separateGenePlots=False,returnImageLoc=True)
            if 'SNE' in tSNE_graphical_links[-1][-1]:
                tSNE_score_file =tSNE_graphical_links[-1][-1][:-10]+'-t-SNE_scores.txt'
                extension = '-t-SNE_scores.txt'
            else:
                tSNE_score_file =tSNE_graphical_links[-1][-1][:-9]+'-UMAP_scores.txt'
                extension = '-UMAP_scores.txt'
        except Exception:
            print traceback.format_exc()
            pass
    
    if '-steady-state' in expFile:
        newExpFile = string.replace(expFile,'-steady-state','-ICGS-steady-state')
        ICGS_order = sampleIndexSelection.getFilters(new_groups_dir)
        sampleIndexSelection.filterFile(expFile,newExpFile,ICGS_order)
        
        ### Copy the steady-state files for ICGS downstream-specific analyses
        ssCountsFile = string.replace(expFile,'exp.','counts.')
        newExpFile = string.replace(expFile,'-steady-state','-ICGS-steady-state')
        newssCountsFile = string.replace(newExpFile,'exp.','counts.')                    
        exonExpFile = string.replace(expFile,'-steady-state','')
        exonCountFile = string.replace(exonExpFile,'exp.','counts.')
        newExonExpFile = string.replace(newExpFile,'-steady-state','')
        newExonCountsFile = string.replace(newExonExpFile,'exp.','counts.')

        sampleIndexSelection.filterFile(ssCountsFile,newssCountsFile,ICGS_order)
        sampleIndexSelection.filterFile(exonExpFile,newExonExpFile,ICGS_order)
        sampleIndexSelection.filterFile(exonCountFile,newExonCountsFile,ICGS_order)
        exonExpFile = exonExpFile[:-4]+'-ICGS.txt'
    else:
        newExpFile = expFile[:-4]+'-ICGS.txt'
        ICGS_order = sampleIndexSelection.getFilters(new_groups_dir)
        sampleIndexSelection.filterFile(expFile,newExpFile,ICGS_order)
        exonExpFile = newExpFile
    
    if outputTSNE:
        try:
            status = verifyFile(tSNE_score_file)
            if status=='no':
                tSNE_score_file = string.replace(tSNE_score_file,'Clustering-','')
            ### Copy the t-SNE scores to use it for gene expression analyses
            if 'DataPlots' not in tSNE_score_file:
                outdir = export.findParentDir(export.findParentDir(tSNE_score_file)[:-1])+'/DataPlots'
            else:
                outdir = export.findParentDir(tSNE_score_file)
            exp_tSNE_score_file = outdir+'/'+export.findFilename(exonExpFile)[:-4]+extension
            import shutil
            shutil.copyfile(tSNE_score_file,exp_tSNE_score_file)
        except Exception:
            print traceback.format_exc()
            pass
    
    return exonExpFile,newExpFile,new_groups_dir 
        
def predictSampleExpGroups(expFile, mlp_instance, gsp, reportOnly, root, exportAdditionalResults=True):
    global graphic_links; graphic_links=[];
    if root == None: display=False
    else: display=True
    
    import RNASeq,ExpressionBuilder; reload(RNASeq) ### allows for GUI testing with restarting
    try:
        if gsp.FeaturestoEvaluate() != 'AltExon':
            from stats_scripts import ICGS_NMF
            reload(ICGS_NMF)
            scaling = True ### Performs pagerank downsampling if over 2,500 cells - currently set as a hard coded default
            dynamicCorrelation=True
            graphic_links=ICGS_NMF.runICGS_NMF(expFile,scaling,array_type,species,gsp,enrichmentInput='',dynamicCorrelation=True)
            #graphic_links = RNASeq.singleCellRNASeqWorkflow(species, array_type, expFile, mlp_instance, parameters=gsp, reportOnly=reportOnly)
        if gsp.FeaturestoEvaluate() != 'Genes':
            ### For splice-ICGS (needs to be updated in a future version to ICGS_NMF updated code)
            graphic_links2,cluster_input_file=ExpressionBuilder.unbiasedComparisonSpliceProfiles(fl.RootDir(),species,array_type,expFile=fl.CountsFile(),min_events=gsp.MinEvents(),med_events=gsp.MedEvents())
            gsp.setCountsCutoff(0);gsp.setExpressionCutoff(0)  
            graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', cluster_input_file, mlp_instance, parameters=gsp, reportOnly=reportOnly)
            graphic_links+=graphic_links2+graphic_links3
        print_out  = 'Predicted sample groups saved.'
        
        if exportAdditionalResults:
            ### Optionally automatically generate t-SNE and MarkerFinder Results
            guide3_results = graphic_links[-1][-1][:-4]+'.txt'
            exportAdditionalICGSOutputs(expFile,guide3_results)

            
        if len(graphic_links)==0:
            print_out  = 'No predicted sample groups identified. Try different parameters.'
        if display == False: print print_out
        try: InfoWindow(print_out, 'Continue')
        except Exception: None
        try: root.destroy()
        except Exception: pass
    except Exception:
        error = traceback.format_exc()
        if 'score_ls' in error:
            error = 'Unknown error likely due to too few genes resulting from the filtering options.'
        if 'options_result_in_no_genes' in error:
            error = 'ERROR: No genes differentially expressed with the input criterion'
        print_out = 'Predicted sample export failed..\n'+error
        try: print print_out
        except Exception: pass ### Windows issue with the Tk status window stalling after pylab.show is called
        try: WarningWindow(print_out,'Continue')
        except Exception: pass
        try: root.destroy()
        except Exception: pass
    try: print error
    except Exception: pass
        
def openDirectory(output_dir):
    if runningCommandLine:
        pass
    elif os.name == 'nt':
        try: os.startfile('"'+output_dir+'"')
        except Exception:  os.system('open "'+output_dir+'"')
    elif 'darwin' in sys.platform: os.system('open "'+output_dir+'"')
    elif 'linux' in sys.platform: os.system('xdg-open "'+output_dir+'"')
        
def networkBuilder(inputDir,inputType,outputdir,interactionDirs_short,degrees,expressionFile,gsp,root):
    species = gsp.Species()
    Genes = gsp.GeneSelection()
    PathwaySelect = gsp.PathwaySelect()
    OntologyID = gsp.OntologyID()
    GeneSet = gsp.GeneSet()
    IncludeExpIDs = gsp.IncludeExpIDs()
    if 'Ontology' in GeneSet: directory = 'nested'
    else: directory = 'gene-mapp'
    interactionDirs=[]
    obligatorySet=[] ### Always include interactions from these if associated with any input ID period
    secondarySet=[]
    print 'Species:',species, '| Algorithm:',degrees, ' | InputType:',inputType, ' | IncludeExpIDs:',IncludeExpIDs
    print 'Genes:',Genes
    print 'OntologyID:',gsp.OntologyID(), gsp.PathwaySelect(), GeneSet
    print ''
    if interactionDirs_short == None or len(interactionDirs_short)==0:
        interactionDirs_short = ['WikiPathways']
    for i in interactionDirs_short:
        if i == None: None
        else:
            if 'common-' in i:
                i = string.replace(i,'common-','')
                secondarySet.append(i)
            if 'all-' in i:
                i = string.replace(i,'all-','')
                obligatorySet.append(i)
            fn = filepath('AltDatabase/goelite/'+species+'/gene-interactions/Ensembl-'+i+'.txt')
            interactionDirs.append(fn)
    print "Interaction Files:",string.join(interactionDirs_short,' ')
    import InteractionBuilder
    try:
        output_filename = InteractionBuilder.buildInteractions(species,degrees,inputType,inputDir,outputdir,interactionDirs,Genes=Genes,
                      geneSetType=GeneSet,PathwayFilter=PathwaySelect,OntologyID=OntologyID,directory=directory,expressionFile=expressionFile,
                      obligatorySet=obligatorySet,secondarySet=secondarySet,IncludeExpIDs=IncludeExpIDs)
        if output_filename==None:
            print_out = 'Network creation/visualization failed..\nNo outputs produced... try different options.\n'
            print_out += traceback.format_exc()
            if root != None and root != '':
                try: InfoWindow(print_out, 'Continue')
                except Exception: None
        else:
            if root != None and root != '':
                try: openDirectory(outputdir)
                except Exception: None
            else:
                print 'Results saved to:',output_filename
        if root != None and root != '':
            GUI(root,'ViewPNG',[],output_filename) ### The last is default attributes (should be stored as defaults in the option_db var)
    except Exception:
        error = traceback.format_exc()
        if 'queryGeneError' in error:
            print_out = 'No valid gene IDs present in the input text search\n(valid IDs = FOXP1,SOX2,NANOG,TCF7L1)'
        else: print_out = 'Network creation/visualization failed..\n',error
        if root != None and root != '':
            try: InfoWindow(print_out, 'Continue')
            except Exception: None
    try: root.destroy()
    except Exception: None
    
def vennDiagram(files_to_merge, output_venn_dir, root, display=True):
    from visualization_scripts import VennDiagram
    if root == None and display==False: display=False
    else: display=True
    try:
        VennDiagram.compareInputFiles(files_to_merge,output_venn_dir,display=display)
        if display == False: print 'VennDiagrams saved to:',output_venn_dir
    except Exception:
        error = traceback.format_exc()
        print_out = 'Venn Diagram export failed..\n',error
        if root != None and root != '':
            try: InfoWindow(print_out, 'Continue')
            except Exception: None
    try: root.destroy()
    except Exception: None
    
def altExonViewer(species,platform,exp_file,gene,show_introns,analysisType,root):
    from visualization_scripts import QC
    transpose=True
    if root == None: display = False
    else: display = True
    if analysisType == 'Sashimi-Plot':
        showEvent = False
        try:
            ### Create sashimi plot index
            from visualization_scripts import SashimiIndex
            print 'Indexing splicing-events'
            SashimiIndex.remoteIndexing(species,exp_file)
            from visualization_scripts import SashimiPlot
            #reload(SashimiPlot)
            print 'Running Sashimi-Plot...'
            genes=None
            if '.txt' in gene:
                events_file = gene
                events = None
            else:
                gene = string.replace(gene,',',' ')
                genes = string.split(gene,' ')
                events_file = None
                if len(genes)==1:
                    showEvent = True
            SashimiPlot.remoteSashimiPlot(species,exp_file,exp_file,events_file,events=genes,show=showEvent) ### assuming the bam files are in the root-dir
            if root != None and root != '':
                print_out = 'Sashimi-Plot results saved to:\n'+exp_file+'/SashimiPlots'
                try: InfoWindow(print_out, 'Continue')
                except Exception: None
        except Exception:
            error = traceback.format_exc()
            print_out = 'AltExon Viewer failed..\n',error
            if root != None and root != '':
                try: WarningWindow(print_out, 'Continue')
                except Exception: None
        try: root.destroy()
        except Exception: None
    else:
        #print [analysisType, species,platform,exp_file,gene,transpose,display,show_introns]
        try: QC.displayExpressionGraph(species,platform,exp_file,gene,transpose,display=display,showIntrons=show_introns,analysisType=analysisType)
        except Exception:
            error = traceback.format_exc()
            print_out = 'AltExon Viewer failed..\n',error
            if root != None and root != '':
                try: WarningWindow(print_out, 'Continue')
                except Exception: None
        try: root.destroy()
        except Exception: None
    
def MergeFiles(files_to_merge, join_option, ID_option, output_merge_dir, root):
    from import_scripts import mergeFiles
    try: outputfile = mergeFiles.joinFiles(files_to_merge, join_option, ID_option, output_merge_dir)
    except Exception:
        outputfile = 'failed'
        error = traceback.format_exc()
    if outputfile == 'failed':
        print_out = 'File merge failed due to:\n',error
    else:
        print_out = 'File merge complete. See the new file:\n'+outputfile
    if root != None and root!= '':
        try: InfoWindow(print_out, 'Continue')
        except Exception: None
        try: root.destroy()
        except Exception: None
        if outputfile != 'failed': ### Open the folder
            try: openDirectory(output_merge_dir)
            except Exception: None
    
def IDconverter(filename, species_code, input_source, output_source, root):
    import gene_associations
    try: outputfile = gene_associations.IDconverter(filename, species_code, input_source, output_source)
    except Exception:
        outputfile = 'failed'
        error = traceback.format_exc()
    if outputfile == 'failed':
        print_out = 'Translation failed due to:\n',error
        print print_out
    else:
        print_out = 'ID translation complete. See the new file:\n'+outputfile
    if root != None and root!= '':
        try: InfoWindow(print_out, 'Continue')
        except Exception: None
        try: root.destroy()
        except Exception: None
        if outputfile != 'failed': ### Open the folder
            try: openDirectory(export.findParentDir(filename))
            except Exception: None

def remoteLP(fl, expr_input_dir, vendor, custom_markerFinder, geneModel, root, modelSize=None,CenterMethod='centroid'):
    global species; global array_type
    species = fl.Species()
    array_type = fl.PlatformType()
    runLineageProfiler(fl, expr_input_dir, vendor, custom_markerFinder, geneModel, root, modelSize=modelSize,CenterMethod=CenterMethod)

def runLineageProfiler(fl, expr_input_dir, vendor, custom_markerFinder, geneModel, root, modelSize=None,CenterMethod='centroid'):
    try:
        global logfile
        root_dir = export.findParentDir(expr_input_dir)
        time_stamp = AltAnalyze.timestamp()    
        logfile = filepath(root_dir+'/AltAnalyze_report-'+time_stamp+'.log')
    except: pass
    
    if custom_markerFinder == '': custom_markerFinder = False
    if modelSize != None and modelSize != 'no':
        try: modelSize = int(modelSize)
        except Exception: modelSize = 'optimize'
    try:
        classificationAnalysis = fl.ClassificationAnalysis()
        if custom_markerFinder == False and classificationAnalysis == 'cellHarmony':
            classificationAnalysis = 'LineageProfiler'
    except Exception:
        if custom_markerFinder == False:
            classificationAnalysis = 'LineageProfiler'
        else:  
            classificationAnalysis = 'cellHarmony'

    if ((geneModel == None or geneModel == False) and (modelSize == None or modelSize == 'no')) and classificationAnalysis != 'cellHarmony':
        print 'LineageProfiler'
        import ExpressionBuilder; reload(ExpressionBuilder)
        compendium_type = fl.CompendiumType()
        compendium_platform = fl.CompendiumPlatform()
        if 'exp.' in expr_input_dir:
            ### Correct the input file to be the gene-expression version
            if array_type != "3'array" and 'AltExon' not in compendium_type:
                if 'steady' not in expr_input_dir:
                    expr_input_dir = string.replace(expr_input_dir,'.txt','-steady-state.txt')
    
        print '****Running LineageProfiler****'
        graphic_links = ExpressionBuilder.remoteLineageProfiler(fl,expr_input_dir,array_type,species,vendor,customMarkers=custom_markerFinder,specificPlatform=True,visualizeNetworks=False)
        if len(graphic_links)>0:
            print_out = 'Lineage profiles and images saved to the folder "DataPlots" in the input file folder.'
            try: InfoWindow(print_out, 'Continue')
            except Exception: None
        else:
            print_out = 'Analysis error occured...\nplease see warning printouts.'
            try: print print_out
            except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
            try: WarningWindow(print_out,'Continue')
            except Exception: None
        try: root.destroy()
        except Exception: None
    else:
        import LineageProfilerIterate
        reload(LineageProfilerIterate)
        print '****Running cellHarmony****'
        codingtype = 'exon'; compendium_platform = 'exon'
        platform = array_type,vendor
        try: LineageProfilerIterate.runLineageProfiler(species,platform,expr_input_dir,expr_input_dir,codingtype,compendium_platform,customMarkers=custom_markerFinder,geneModels=geneModel,modelSize=modelSize,fl=fl)
        except Exception:
            print_out = traceback.format_exc()
            try: InfoWindow(print_out, 'Continue') ### Causes an error when peforming heatmap visualizaiton
            except Exception: None
        print_out = 'LineageProfiler classification results saved to the folder "CellClassification".'
        if root!=None and root!='':
            try: openDirectory(export.findParentDir(expr_input_dir)+'/cellHarmony')
            except Exception: None
            try: InfoWindow(print_out, 'Continue') ### Causes an error when peforming heatmap visualizaiton
            except Exception: None
        else:
            print print_out

        try: root.destroy()
        except Exception: None

    
def performPCA(filename, pca_labels, pca_algorithm, transpose, root, plotType='3D',display=True,
            geneSetName=None, species=None, zscore=True, colorByGene=None, reimportModelScores=True,
            separateGenePlots=False, returnImageLoc=False):
    from visualization_scripts import clustering; reload(clustering)
    graphics = []
    if pca_labels=='yes' or pca_labels=='true'or pca_labels=='TRUE': pca_labels=True
    else: pca_labels=False
    if zscore=='yes': zscore = True
    elif zscore=='no': zscore = False
    pca_graphical_links=[]
    try:
        pca_graphical_links = clustering.runPCAonly(filename, graphics, transpose, showLabels=pca_labels,
                    plotType=plotType,display=display, algorithm=pca_algorithm, geneSetName=geneSetName,
                    species=species, zscore=zscore, colorByGene=colorByGene, reimportModelScores=reimportModelScores,
                    separateGenePlots=separateGenePlots)
        try: print'Finished building exporting plot.'
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
    except Exception:
        if 'importData' in traceback.format_exc():
            try: print traceback.format_exc(),'\n'
            except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
            print_out = 'Bad input file! Should be a tab-delimited text file with a single\nannotation column and row and the remaining as numeric values.'
        else:
            try: print traceback.format_exc(),'\n'
            except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
            print_out = 'Analysis error occured...\nplease try again with different parameters.'
        try: print print_out
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
        try: WarningWindow(print_out,'Continue')
        except Exception: None
    try: root.destroy()
    except Exception: pass
    if returnImageLoc:
        try: return pca_graphical_links
        except Exception: pass
    
def createHeatMap(filename, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast, root, display=True):
    graphics = []
    try:
        from visualization_scripts import clustering; reload(clustering)
        clustering.runHCexplicit(filename, graphics, row_method, row_metric, column_method, column_metric, color_gradient, transpose, display=display, contrast = contrast)
        print_out = 'Finished building heatmap.'
        try: print print_out
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
        try: root.destroy()
        except Exception: pass ### DO NOT PRINT HERE... CONFLICTS WITH THE STOUT
    except Exception:
        if 'importData' in traceback.format_exc():
            try: print traceback.format_exc(),'\n'
            except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
            print_out = 'Bad input file! Should be a tab-delimited text file with a single\nannotation column and row and the remaining as numeric values.'
        else:
            try: print traceback.format_exc(),'\n'
            except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
            print_out = 'Analysis error occured...\nplease try again with different parameters.'
        try: print print_out
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
        try: WarningWindow(print_out,'Continue')
        except Exception: None
        try: root.destroy()
        except Exception: pass
    
def getAdditionalOnlineResources(species_code,additional_resources,root):
    if additional_resources[0] == 'customSet':
        additional_resources = additional_resources[1]
    elif additional_resources == 'All Resources':
        additional_resources = importResourceList()
    else: additional_resources = [additional_resources]
    try:
        print 'Adding supplemental GeneSet and Ontology Collections'
        from build_scripts import GeneSetDownloader; force = 'yes'
        GeneSetDownloader.buildAccessoryPathwayDatabases([species_code],additional_resources,force)
        try: print'Finished incorporating additional resources.'
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
    except Exception:
        print_out = 'Download error encountered for additional ontologies and gene-sets...\nplease try again later.'
        try: print print_out
        except Exception: None ### Windows issue with the Tk status window stalling after pylab.show is called
        try: WarningWindow(print_out,'Continue')
        except Exception: None
    try: root.destroy()
    except Exception: pass
    
class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s): ### Write is called by python when any new print statement is called
        new = self.__stringvar.get()
        for c in s:
            #if c == '\n': self.__newline = 1
            if c == '\k': self.__newline = 1 ### This should not be found and thus results in a continous feed rather than replacing a single line
            else:
                if self.__newline: new = ""; self.__newline = 0
                new = new+c
        try: self.set(new)
        except Exception: pass
        #except Exception: None ### Not sure why this occurs
        try:
            log_report = open(logfile,'a')
            log_report.write(s); log_report.close() ### Variable to record each print statement
        except Exception: pass
    def set(self,s):
        try: self.__stringvar.set(s); self.__window.update()
        except Exception: pass
    def get(self): return self.__stringvar.get()
    def flush(self): pass
    

################# GUI #################

class ProgressBar:
	def __init__(self,method,t):
		#http://tkinter.unpythonic.net/wiki/ProgressBar
		self.progval = IntVar(t)
		self.progmsg = StringVar(t); self.progmsg.set(method+" in Progress...")
		#b = Button(t, relief=LINK, text="Quit (using bwidget)", command=t.destroy); b.pack()
		self.c = ProgressDialog(t, title="Please wait...",
			type="infinite",
			width=30,
			textvariable=self.progmsg,
			variable=self.progval,
			command=lambda: self.c.destroy()
			)
		self.update_progress()
	def update_progress(self):
		self.progval.set(2)
		self.c.after(10, self.update_progress)
	
class ImageFiles:
    def __init__(self,shortname,fullpath,return_gif=False):
        self.shortname = shortname
        self.fullpath = fullpath
        self.return_gif = return_gif
    def ShortName(self): return self.shortname
    def FullPath(self): return self.fullpath
    def returnGIF(self): return self.return_gif
    def Thumbnail(self):
        if self.returnGIF():
            gif_path = string.replace(self.FullPath(),'.png','.gif')
            return gif_path
        else:
            png_path = string.replace(self.FullPath(),'.png','_small.png')
            return png_path
        
class GUI:
    def PredictGroups(self):
        self.button_flag = True
        self.graphic_link = {}
        try: import ImageTk
        except Exception:
            from PIL import ImageTk
            from PIL import Image
        
        self.toplevel_list=[] ### Keep track to kill later
        
        self.filename_db={}
        filenames=[]
        i=1
        for (name,file) in graphic_links:
            self.filename_db['clusters '+str(i)]=file
            filenames.append('clusters '+str(i))
            i+=1
        
        filenames_ls = list(filenames)
        filenames_ls.reverse()
        self.title = 'Select cluster groups for further analysis'
        self.option = 'group_select' ### choose a variable name here
        self.options = filenames_ls
        self.default_option = 0
        self.comboBox() ### This is where the cluster group gets selected and stored
  
        # create a frame and pack it
        frame1 = Tkinter.Frame(self.parent_type)
        frame1.pack(side=Tkinter.TOP, fill=Tkinter.X)
        
        ### Convert PNG to GIF and re-size
        assigned_index=1
        #print filenames
        for image_file in filenames:
            file_dir = self.filename_db[image_file]
            iF = ImageFiles(image_file,file_dir)
            im = Image.open(file_dir) 
            #im.save('Gfi1.gif')
            size = 128, 128
            im.thumbnail(size, Image.ANTIALIAS)
            im.save(iF.Thumbnail()) ### write out the small gif file
            option = 'imageView'
            self.option=option
            #photo1 = Tkinter.PhotoImage(file=iF.Thumbnail())
            photo1 = ImageTk.PhotoImage(file=iF.Thumbnail()) ### specifically compatible with png files
            # create the image button, image is above (top) the optional text
            def view_FullImageOnClick(image_name):
                tl = Toplevel() #### This is the critical location to allow multiple TopLevel instances that don't clash, that are created on demand (by click)
                self.toplevel_list.append(tl)
                self.graphic_link['WP'] = self.filename_db[image_name]
                try: self.viewPNGFile(tl) ### ImageTK PNG viewer
                except Exception:
                    print traceback.format_exc()
                    try: self.openPNGImage() ### OS default PNG viewer
                    except Exception: pass
            #print assigned_index
            if assigned_index == 1:
                image_file1 = image_file; #tl1 = Toplevel() ### not good to create here if we have to destroy it, because then we can't re-invoke
                button1 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file1, bg='green', command=lambda:view_FullImageOnClick(image_file1)) ### without lamda, the command is called before being clicked
                button1.image = photo1; button1.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 2:
                image_file2 = image_file; #tl2 = Toplevel()
                button2 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file2, bg='green', command=lambda:view_FullImageOnClick(image_file2))
                button2.image = photo1; button2.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 3:
                image_file3 = image_file; #tl3 = Toplevel()
                button3 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file3, bg='green', command=lambda:view_FullImageOnClick(image_file3))
                button3.image = photo1; button3.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 4:
                image_file4 = image_file; #tl4 = Toplevel()
                button4 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file4, bg='green', command=lambda:view_FullImageOnClick(image_file4))
                button4.image = photo1; button4.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 5:
                image_file5 = image_file; #tl5 = Toplevel()
                button5 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file5, bg='green', command=lambda:view_FullImageOnClick(image_file5))
                button5.image = photo1; button5.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 6:
                image_file6 = image_file; #tl4 = Toplevel()
                button6 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file6, bg='green', command=lambda:view_FullImageOnClick(image_file6))
                button6.image = photo1; button6.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 7:
                image_file7 = image_file; #tl5 = Toplevel()
                button7 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file7, bg='green', command=lambda:view_FullImageOnClick(image_file7))
                button7.image = photo1; button7.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 8:
                image_file8 = image_file; #tl5 = Toplevel()
                button8 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file8, bg='green', command=lambda:view_FullImageOnClick(image_file8))
                button8.image = photo1; button8.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 9:
                image_file9 = image_file; #tl4 = Toplevel()
                button9 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file9, bg='green', command=lambda:view_FullImageOnClick(image_file9))
                button9.image = photo1; button9.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 10:
                image_file10 = image_file; #tl5 = Toplevel()
                button10 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file10, bg='green', command=lambda:view_FullImageOnClick(image_file10))
                button10.image = photo1; button10.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 11:
                image_file11 = image_file; #tl5 = Toplevel()
                button11 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file11, bg='green', command=lambda:view_FullImageOnClick(image_file11))
                button11.image = photo1; button11.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 12:
                image_file12 = image_file; #tl5 = Toplevel()
                button12 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file12, bg='green', command=lambda:view_FullImageOnClick(image_file12))
                button12.image = photo1; button12.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 13:
                image_file13 = image_file; #tl5 = Toplevel()
                button13 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file13, bg='green', command=lambda:view_FullImageOnClick(image_file13))
                button13.image = photo1; button13.pack(side=Tkinter.TOP, padx=2, pady=2)
            elif assigned_index == 14:
                image_file14 = image_file; #tl5 = Toplevel()
                button14 = Tkinter.Button(frame1, compound=Tkinter.TOP, image=photo1,
                    text=image_file14, bg='green', command=lambda:view_FullImageOnClick(image_file14))
                button14.image = photo1; button14.pack(side=Tkinter.TOP, padx=2, pady=2)
            assigned_index+=1

        # start the event loop
        use_selected_button = Button(self._parent, text="Use Selected", command=self.UseSelected) 
        use_selected_button.pack(side = 'right', padx = 10, pady = 5)
        
        recluster_button = Button(self._parent, text="Re-Cluster", command=self.ReCluster) 
        recluster_button.pack(side = 'right', padx = 10, pady = 5)

        quit_button = Button(self._parent, text="Quit", command=self.quit) 
        quit_button.pack(side = 'right', padx = 10, pady = 5)

        try: help_button = Button(self._parent, text='Help', command=self.GetHelpTopLevel); help_button.pack(side = 'left', padx = 5, pady = 5)
        except Exception: help_button = Button(self._parent, text='Help', command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()
        
    def UseSelected(self):
        status = self.checkAllTopLevelInstances()
        if status:
            self.checkAllTopLevelInstances()
            self._user_variables['next'] = 'UseSelected'
            try: self._parent.quit(); self._parent.destroy()
            except Exception: self._parent.quit()

    def ReCluster(self):
        status = self.checkAllTopLevelInstances()
        if status:
            self._user_variables['next'] = 'ReCluster'
            try: self._parent.quit(); self._parent.destroy()
            except Exception:
                try: self._parent.destroy()
                except Exception: pass

    def checkAllTopLevelInstances(self):
        ### Ideally, we would just kill any open toplevel instances, but this was causing a "ghost" process
        ### to continue running even after all of the tls and roots were destroyed
        if len(self.toplevel_list)>0:
            removed=[]
            for tl in self.toplevel_list:
                try:
                    if 'normal' == tl.state():
                        InfoWindow('Please close all cluster windows before proceeding.', 'Continue')
                        break
                except Exception:
                    removed.append(tl)
            for tl in removed:
                self.toplevel_list.remove(tl)
        if len(self.toplevel_list)==0:
            return True
        else:
            return False
                
    def killAllTopLevelInstances(self):
        ### destroy's any live TopLevel instances
        removed=[]
        for tl in self.toplevel_list:
            try: tl.quit(); tl.destroy(); removed.append(tl)
            except Exception: pass
        for tl in removed:
            self.toplevel_list.remove(tl)
        
    def ViewWikiPathways(self):
        """ Canvas is already drawn at this point from __init__ """
        global pathway_db
        pathway_db={}
        button_text = 'Help'

        ### Create a species drop-down option that can be updated
        current_species_names,manufacturers_list = getSpeciesList('') ### pass the variable vendor to getSpeciesList (none in this case) --- different than the GO-Elite UI call
        self.title = 'Select species to search for WikiPathways '
        self.option = 'species_wp'
        self.options = ['---']+current_species_names #species_list
        self.default_option = 0
        self.comboBox()
        
        ### Create a label that can be updated below the dropdown menu
        self.label_name = StringVar()
        self.label_name.set('Pathway species list may take several seconds to load')
        self.invokeLabel() ### Invoke a new label indicating that the database is loading
                    
        ### Create a MOD selection drop-down list
        system_list,mod_list = importSystemInfo() ### --- different than the GO-Elite UI call
        self.title = 'Select the ID system to translate to (MOD)'
        self.option = 'mod_wp'
        self.options = mod_list
        try: self.default_option = mod_list.index('Ensembl') ### Get the Ensembl index number
        except Exception: self.default_option = 0
        self.dropDown()
        
        ### Create a file selection option
        self.title = 'Select GO-Elite input ID text file'
        self.notes = 'note: ID file must have a header row and at least three columns:\n'
        self.notes += '(1) Identifier, (2) System Code, (3) Value to map (- OR +)\n'
        self.file_option = 'goelite_input_file'
        self.directory_type = 'file'
        self.FileSelectionMenu()
        
        dispaly_pathway = Button(text = 'Display Pathway', command = self.displayPathway)
        dispaly_pathway.pack(side = 'right', padx = 10, pady = 10)

        back_button = Button(self._parent, text="Back", command=self.goBack) 
        back_button.pack(side = 'right', padx =10, pady = 5)
        
        quit_win = Button(self._parent, text="Quit", command=self.quit) 
        quit_win.pack(side = 'right', padx =10, pady = 5)

        try: help_button = Button(self._parent, text=button_text, command=self.GetHelpTopLevel); help_button.pack(side = 'left', padx = 5, pady = 5)
        except Exception: help_button = Button(self._parent, text=button_text, command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()

    def FileSelectionMenu(self):
        option = self.file_option
        group = PmwFreeze.Group(self.parent_type,tag_text = self.title)
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
        
        def filecallback(callback=self.callback,option=option): self.getPath(option)
        default_option=''
        entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
        entrytxt.set(default_option)
        self.pathdb[option] = entrytxt
        self._user_variables[option] = default_option
        entry = Entry(group.interior(),textvariable=self.pathdb[option]);
        entry.pack(side='left',fill = 'both', expand = 0.7, padx = 10, pady = 2)
        button = Button(group.interior(), text="select "+self.directory_type, width = 10, fg="black", command=filecallback)
        button.pack(side=LEFT, padx = 2,pady = 2)
        if len(self.notes)>0: ln = Label(self.parent_type, text=self.notes,fg="blue"); ln.pack(padx = 10)
        
    def dropDown(self):
        def comp_callback(tag,callback=self.callbackWP,option=self.option):
            callback(tag,option)
        self.comp = PmwFreeze.OptionMenu(self.parent_type,
            labelpos = 'w', label_text = self.title, items = self.options, command = comp_callback)
        if self.option == 'wp_id_selection':
            self.wp_dropdown = self.comp ### update this variable later (optional)
        self.comp.pack(anchor = 'w', padx = 10, pady = 0, fill = 'x')
        self.comp.invoke(self.default_option) ###Just pick the first option

    def comboBox(self):
        """ Alternative, more sophisticated UI than dropDown (OptionMenu).
        Although it behaves similiar it requires different parameters, can not be
        as easily updated with new lists (different method) and requires explict
        invokation of callback when a default is set rather than selected. """
        
        def comp_callback(tag,callback=self.callbackWP,option=self.option):
            callback(tag,option)
        self.comp = PmwFreeze.ComboBox(self.parent_type,
            labelpos = 'w', dropdown=1, label_text = self.title,
            unique = 0, history = 0,
            scrolledlist_items = self.options, selectioncommand = comp_callback)

        try: self.comp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.comp.invoke())
        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
        
        if self.option == 'wp_id_selection':
            self.wp_dropdown = self.comp ### update this variable later (optional)
        self.comp.pack(anchor = 'w', padx = 10, pady = 0)
        try: self.comp.selectitem(self.default_option) ###Just pick the first option
        except Exception: pass
        try: self.callbackWP(self.options[0],self.option)  ### Explicitly, invoke first option (not automatic)
        except Exception: pass
        
    def invokeLabel(self):
        self.label_object = Label(self.parent_type, textvariable=self.label_name,fg="blue"); self.label_object.pack(padx = 10)
        
    def enterMenu(self):
        if len(self.notes)>0:
            lb = Label(self.parent_type, text=self.notes,fg="black"); lb.pack(pady = 5)
        ### Create and pack a horizontal RadioSelect widget
        def custom_validate(tag,custom_validate=self.custom_validate,option=self.option):
            validate = custom_validate(tag,self.option)
        self.entry_field = PmwFreeze.EntryField(self.parent_type,
                labelpos = 'w', label_text = self.title, validate = custom_validate, 
                value = self.default_option, hull_borderwidth = 2)
        self.entry_field.pack(fill = 'x', expand = 0.7, padx = 10, pady = 5)

    def displayAnyPNG(self,png_file):
        self.graphic_link={}
        self.graphic_link['WP'] = png_file
        self.graphic_link['quit']=None
        try: tl = Toplevel()
        except Exception:
            import Tkinter
            tl = Tkinter.Toplevel()
        try: self.viewPNGFile(tl) ### ImageTK PNG viewer
        except Exception:
            try: self.openPNGImage() ### OS default PNG viewer
            except Exception:
                print 'Unable to open PNG file for unknown reasons'
                        
    def displayPathway(self):
        filename = self._user_variables['goelite_input_file']
        mod_type = self._user_variables['mod_wp']
        species = self._user_variables['species_wp']
        pathway_name = self._user_variables['wp_id_selection']
        wpid_selected = self._user_variables['wp_id_enter']
        species_code = species_codes[species].SpeciesCode()
        wpid = None
        if len(wpid_selected)>0:
            wpid = wpid_selected
        elif len(self.pathway_db)>0:
            for wpid in self.pathway_db:
                if pathway_name == self.pathway_db[wpid].WPName():
                    break
        if len(filename)==0:
            print_out = 'Select an input ID file with values first'
            WarningWindow(print_out,'Error Encountered!')
        else:
            try:
                self.graphic_link = WikiPathways_webservice.visualizePathwayAssociations(filename,species_code,mod_type,wpid)
                if len(self.graphic_link)==0:
                    force_no_matching_error
                self.wp_status = 'Pathway images colored and saved to disk by webservice\n(see image title for location)'
                self.label_status_name.set(self.wp_status)
                tl = Toplevel()
                try: self.viewPNGFile(tl) ### ImageTK PNG viewer
                except Exception:
                    try: self.openPNGImage() ### OS default PNG viewer
                    except Exception:
                        self.wp_status = 'Unable to open PNG file using operating system'
                        self.label_status_name.set(self.wp_status)
            except Exception,e:
                try:
                    wp_logfile = filepath('webservice.log')
                    wp_report = open(wp_logfile,'a')
                    wp_report.write(traceback.format_exc())
                except Exception:
                    None
                try:
                    print traceback.format_exc()
                except Exception:
                    null=None ### Occurs when transitioning back from the Official Database download window (not sure why) -- should be fixed in 1.2.4 (sys.stdout not re-routed)
                if 'force_no_matching_error' in traceback.format_exc():
                    print_out = 'None of the input IDs mapped to this pathway'
                elif 'force_invalid_pathway' in traceback.format_exc():
                    print_out = 'Invalid pathway selected'
                elif 'IndexError' in traceback.format_exc():
                    print_out = 'Input ID file does not have at least 3 columns, with the second column being system code'
                elif 'ValueError' in traceback.format_exc():
                    print_out = 'Input ID file error. Please check that you do not have extra rows with no data' 
                elif 'source_data' in traceback.format_exc():
                    print_out = 'Input ID file does not contain a valid system code' 
                else:
                    print_out = 'Error generating the pathway "%s"' % pathway_name
                WarningWindow(print_out,'Error Encountered!')
        
    def getSpeciesPathways(self,species_full):
        pathway_list=[]
        self.pathway_db = WikiPathways_webservice.getAllSpeciesPathways(species_full)
        for wpid in self.pathway_db:
            if self.pathway_db[wpid].WPName() != None: ### Not sure where the None comes from but will break the UI if not exlcuded
                pathway_list.append(self.pathway_db[wpid].WPName())
        pathway_list = unique.unique(pathway_list)
        pathway_list.sort()
        return pathway_list

    def callbackWP(self, tag, option):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables[option] = tag
        if option == 'group_select':
            ### set group_select equal to the filename
            self._user_variables[option] = self.filename_db[tag]                
            #print option, tag
            #print option, self._user_variables[option], self.filename_db[tag]
        if option == 'species_wp':
            ### Add additional menu options based on user selection
            if tag != '---':
                ### If this already exists from an earlier iteration
                hault = False
                self.label_name.set('Loading available WikiPathways')
                try:
                    self.pathway_list=self.getSpeciesPathways(tag)
                    traceback_printout = ''
                except Exception,e:
                    if 'not supported' in traceback.format_exc():
                        print_out = 'Species not available at WikiPathways'
                        WarningWindow(print_out,'Species Not Found!')
                        traceback_printout=''
                        hault = True
                    elif 'URLError' in traceback.format_exc():
                        print_out = 'Internet connection could not be established'
                        WarningWindow(print_out,'Internet Error')
                        traceback_printout=''
                        hault = True
                    else:
                        traceback_printout = traceback.format_exc()
                    try: 
                        if len(self.pathway_list)>0: ### When true, a valid species was selected in a prior interation invoking the WP fields (need to repopulate)
                            hault = False
                    except Exception: None
                    self.pathway_list = ['None']; self.pathway_db={}
                self.label_name.set('')
                if hault == False:
                    try:
                        ### If the species specific wikipathways drop down exists, just update it
                        self.wp_dropdown._list.setlist(self.pathway_list)
                        self.wp_dropdown.selectitem(self.pathway_list[0])
                        self.callbackWP(self.pathway_list[0],'wp_id_selection')
                    except Exception:
                        ### Create a species specific wikipathways drop down
                        self.option = 'wp_id_selection'
                        self.title = 'Select WikiPathways to visualize your data'
                        if len(traceback_printout)>0:
                            self.title += traceback_printout ### Display the actual problem in the GUI (sloppy but efficient way for users to indicate the missing driver)
                        self.options = self.pathway_list
                        self.default_option = 0
                        self.comboBox() ### Better UI for longer lists of items (dropDown can't scroll on Linux)
                        
                        ### Create a species specific wikipathways ID enter option
                        self.notes = 'OR'
                        self.option = 'wp_id_enter'
                        self.title = 'Enter the WPID (example: WP254) '
                        self.default_option = ''
                        self.enterMenu()
                    try:
                        ### Create a label that can be updated below the dropdown menu
                        
                        self.wp_status = 'Pathway image may take several seconds to a minute to load...\n'
                        self.wp_status += '(images saved to "WikiPathways" folder in input directory)'
                        try: self.label_status_name.set(self.wp_status)
                        except Exception:
                            self.label_status_name = StringVar()
                            self.label_status_name.set(self.wp_status)
                            self.invokeStatusLabel() ### Invoke a new label indicating that the database is loading
                    except Exception:
                        None
                        
        if option == 'wp_id_selection':
            ### Reset any manually input WPID if a new pathway is selected from dropdown
            try: self.entry_field.setentry('') 
            except Exception: pass
         
    def ShowImageMPL(self):
        png_file_dir = self.graphic_link['WP']
        fig = pylab.figure()
        pylab.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.00) ### Fill the plot area left to right
        ax = fig.add_subplot(111)
        ax.set_xticks([]) ### Hides ticks
        ax.set_yticks([])
        img= pylab.imread(png_file_dir)
        imgplot = pylab.imshow(img)
        pylab.show()
        
    def viewPNGFile(self,tl):
        """ View PNG file within a PMW Tkinter frame """
        try: import ImageTk ### HAVE TO CALL HERE TO TRIGGER AN ERROR - DON'T WANT THE TopLevel to open otherwise
        except Exception:
            from PIL import ImageTk
        png_file_dir = self.graphic_link['WP']
        img = ImageTk.PhotoImage(file=png_file_dir)
        
        sf = PmwFreeze.ScrolledFrame(tl, labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 800, hull_height = 550)
        sf.pack(padx = 0, pady = 0, fill = 'both', expand = 1)
        frame = sf.interior()

        tl.title(png_file_dir)
        can = Canvas(frame)
        can.pack(fill=BOTH, padx = 0, pady = 0)
        w = img.width()
        h = height=img.height()
        
        can.config(width=w, height=h)        
        can.create_image(2, 2, image=img, anchor=NW)
        if 'quit' in self.graphic_link:
            tl.protocol("WM_DELETE_WINDOW", lambda: self.tldeleteWindow(tl))
            tl.mainloop()
        else:
            tl.protocol("WM_DELETE_WINDOW", lambda: self.tldeleteWindow(tl))
            tl.mainloop()

    def openPNGImage(self):
        png_file_dir = self.graphic_link['WP']
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+png_file_dir+'"')
            except Exception:  os.system('open "'+png_file_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')   

    def centerPage(self):
            # Example of how to use the yview() method of Pmw.ScrolledFrame.
            top, bottom = self.sf.yview()
            size = bottom - top
            middle = 0.5 - size / 2
            self.sf.yview('moveto', middle)

    def __init__(self, parent, option_db, option_list, defaults):
        if option_db == 'ViewPNG':
            output_filename = defaults
            self.displayAnyPNG(output_filename)
            return None
        self._parent = parent; self._option_list = option_list; self._option_db = option_db
        self._user_variables = user_variables; self.pathdb={}; i = -1
        enter_index=0; radio_index=0; dropdown_index=0; check_index=0 ### used to keep track of how many enter boxes we have
        self.default_dir = PathDir; self.default_file = PathFile
        self.defaults = defaults
        
        filename = 'Config/icon.gif'; orient_type = 'left'

        if 'input_cel_dir' in option_list:
            filename = 'Config/aa_0.gif'
            if array_type == 'RNASeq': filename = 'Config/aa_0_rs.gif'
            if '10X' in vendor: filename = 'Config/aa_0_rs.gif'
        if 'include_raw_data' in option_list:
            filename = 'Config/aa_1.gif'; orient_type = 'top'
            if array_type == 'RNASeq': filename = 'Config/aa_1_rs.gif'
        if 'filter_for_AS' in option_list:
            filename = 'Config/aa_2.gif'; orient_type = 'top'
            if array_type == 'RNASeq': filename = 'Config/aa_2_rs.gif'
        if 'pathway_permutations' in option_list: filename = 'Config/goelite.gif'
        if 'GeneSelectionPredict' in option_list:
            filename = 'Config/aa_3.gif'
            if array_type == 'RNASeq': filename = 'Config/aa_3_rs.gif'
        
        fn=filepath(filename); img = PhotoImage(file=fn)
        self.can = Canvas(parent); self.can.pack(side='top'); self.can.config(width=img.width(), height=img.height())        
        try: self.can.create_image(2, 2, image=img, anchor=NW)
        except Exception:
            try: self.can.delete("all")
            except Exception: pass
        #except Exception: print filename; 'what?';kill
        self.pathdb={}; use_scroll = 'no'

        #if defaults == 'groups' or defaults == 'comps' or 'filter_for_AS' in option_list:
        if defaults != 'null':
            height = 350; width = 400
            if defaults == 'groups':
                notes = "For each sample, type in a name for the group it belongs to\n(e.g., 24hrs, 48hrs, 4days, etc.)."
                Label(self._parent,text=notes).pack(); label_text_str = 'AltAnalyze Group Names'
                if len(option_list)<15: height = 320; width = 400

            elif defaults == 'batch':
                notes = "For each sample, type in a name for the BATCH it belongs to\n(e.g., batch1, batch2, batch3 etc.)."
                Label(self._parent,text=notes).pack(); label_text_str = 'AltAnalyze Group Names'
                if len(option_list)<15: height = 320; width = 400
            elif defaults == 'comps':
                notes = "Experimental Group\t\t\tBaseline Group     "
                label_text_str = 'AltAnalyze Pairwise Group Comparisons'
                if len(option_list)<5: height = 250; width = 400
            elif 'filter_for_AS' in option_list:
                label_text_str = 'AltAnalyze Alternative Exon Analysis Parameters'
                height = 350; width = 400; use_scroll = 'yes'
                if os.name != 'nt': width+=100
            elif 'pathway_permutations' in option_list:
                label_text_str = 'GO-Elite Parameters'
                height = 350; width = 425; use_scroll = 'yes'
            elif 'expression_data_format' in option_list:
                label_text_str = "AltAnalyze Expression Dataset Parameters"
                height = 350; width = 400; use_scroll = 'yes'
                if os.name != 'nt': width+=100
            elif 'Genes_network' in option_list:
                label_text_str = "Network Analysis Parameters"
                height = 350; width = 400; use_scroll = 'yes'
                #if os.name != 'nt': width+=50
            elif 'GeneSelectionPredict' in option_list:
                notes = "Perform an unsupervised or supervised analysis to identify the\npredominant sample groups via expression clustering"
                Label(self._parent,text=notes).pack()
                label_text_str = "AltAnalyze Prediction Sample Group Parameters"
                height = 310; width = 400; use_scroll = 'yes'
            elif 'join_option' in option_list:
                label_text_str = "AltAnalyze Merge Files Parameters"
                height = 310; width = 400; use_scroll = 'yes'
            else:
                label_text_str = "AltAnalyze Main Dataset Parameters"
                height = 310; width = 400; use_scroll = 'yes'
                
            if os.name != 'nt':height+=75; width+=150
            if os.name== 'nt':height+=25; width+=50
            if 'linux' in sys.platform: offset = 25
            else: offset=0
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = label_text_str,
                    usehullsize = 1, hull_width = width-offset, hull_height = height-offset)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            if defaults == 'comps':
                Label(self.frame,text=notes).pack()

        create_group = 'yes'
        if 'pathway_permutations' in option_list or 'expression_data_format' in option_list or 'filter_probe_types' in option_list:
            if 'ge_ptype' in option_list:
                self.group_tag = 'GO-Elite Gene Expression Analysis Filters'
            elif 'pathway_permutations' in option_list:
                self.group_tag = 'GO-Elite Over-Representation and Filtering Parameters'
            if 'expression_threshold' in option_list:
                self.group_tag = 'Exon/Junction Filtering Options'
                od = option_db['expression_threshold']
                if od.ArrayOptions() == ['NA']: create_group = 'no'
                if ('rpkm_threshold' in option_list and create_group== 'no'):
                    create_group='yes'
                    self.group_tag = 'Gene Expression Filtering Options'
                    od = option_db['rpkm_threshold']
                    if od.ArrayOptions() == ['NA']: create_group = 'no'
            elif 'expression_data_format' in option_list and 'rpkm_threshold' not in option_list:
                self.group_tag = 'Gene Expression Analysis Options'
            if 'filter_probe_types' in option_list:
                self.group_tag = 'Primary Alternative Exon Parameters'
            if create_group == 'yes': 
                custom_group = PmwFreeze.Group(self.sf.interior(),tag_text = self.group_tag)
                custom_group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
                insert_into_group = 'yes'
            else: insert_into_group = 'no'
        else: insert_into_group = 'no'

        object_directions = ['top','bottom','up','down']

        if option_db == 'ViewWikiPathways':
            width = 520
            self.parent_type = self.sf.interior()
            self.ViewWikiPathways()
        if option_db == 'PredictGroups':
            width = 520
            self.parent_type = self.sf.interior()
            self.PredictGroups()
        for option in option_list:
          i+=1 ####Keep track of index - if options are deleted, count these to select the appropriate default from defaults
          if option in option_db:
            od = option_db[option]; self.title = od.Display(); notes = od.Notes()
            self.display_options = od.ArrayOptions()
            try: override_default = od.DefaultOption()
            except Exception: override_default = ''
            if 'radio' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if 'pathway_permutations' in option_list or 'new_run' in option_list: orient_type = 'top'
                if insert_into_group == 'yes': parent_type = custom_group.interior(); radio_index+=1

                ### Create and pack a RadioSelect widget, with radiobuttons.
                self._option = option
                def radiocallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                radiobuttons = PmwFreeze.RadioSelect(parent_type,                       
                        buttontype = 'radiobutton', orient = 'vertical',
                        labelpos = 'w', command = radiocallback, label_text = self.title,
                        hull_borderwidth = 2, hull_relief = 'ridge')
                if insert_into_group == 'no': radiobuttons.pack(side = orient_type, expand = 1, padx = 10, pady = 10)
                elif radio_index == 1: radiobuttons1 = radiobuttons
                elif radio_index == 2: radiobuttons2 = radiobuttons
                ### print self.display_options
                ### Add some buttons to the radiobutton RadioSelect.
                for text in self.display_options:
                    if text != ['NA']: radiobuttons.add(text)
                if len(override_default)>0: self.default_option = override_default
                elif len(defaults) <1:
                    try: self.default_option = self.display_options[0]
                    except Exception: print option; kill
                else: self.default_option = defaults[i]
                radiobuttons.invoke(self.default_option)
                if len(notes)>0: Label(self._parent, text=notes).pack()
            if 'button' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                
                self._option = option
                if mac_print_mode == 'yes' or 'radbutton' in od.DisplayObject(): button_type = 'radiobutton'
                else: button_type = 'button'
                ### Create and pack a horizontal RadioSelect widget.
                if len(override_default)>0: self.default_option = override_default
                elif len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                def buttoncallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                    #self.sf.update_idletasks()
                    #self.centerPage()
                orientation = 'vertical'
                #if 'pathway_permutations' in option_list or 'new_run' in option_list: orientation = 'vertical'
                #elif 'run_from_scratch' in option_list: orientation = 'vertical'
                #else: orientation = 'vertical'
                horiz = PmwFreeze.RadioSelect(parent_type, buttontype = button_type, orient = orientation,
                        labelpos = 'w', command = buttoncallback,
                        label_text = self.title, frame_borderwidth = 2,
                        frame_relief = 'ridge'
                ); horiz.pack(fill = 'x',padx = 10, pady = 10)

                ### Add some buttons to the horizontal RadioSelect
                for text in self.display_options:
                    if text != ['NA']: horiz.add(text)
                horiz.invoke(self.default_option)
                if len(notes)>0: Label(self._parent, text=notes).pack()

            if ('folder' in od.DisplayObject() or 'file' in od.DisplayObject()) and self.display_options != ['NA']:
              if use_scroll == 'yes': parent_type = self.sf.interior()
              else: parent_type = self._parent
              if 'sparse-matrix' in self.title:
                od.setDisplayObject('file') ### set the object equal to a file
              proceed = 'yes'
              #if option == 'raw_input': proceed = 'no'
              if proceed == 'yes':
                self._option = option
                group = PmwFreeze.Group(parent_type,tag_text = self.title)
                group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
                def filecallback(callback=self.callback,option=option): self.getPath(option)
                entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
                try: default_option = string.replace(override_default,'---','')
                except Exception: default_option = ''
                entrytxt.set(default_option)
                self.pathdb[option] = entrytxt
                self._user_variables[option] = default_option
                if option == 'input_cel_dir' and '10X' in vendor:
                    od.setDisplayObject('file')
                #l = Label(group.interior(), text=self.title); l.pack(side=LEFT)        
                entry = Entry(group.interior(),textvariable=self.pathdb[option]);
                entry.pack(side='left',fill = 'both', expand = 1, padx = 10, pady = 2)
                button = Button(group.interior(), text="select "+od.DisplayObject(), width = 10, fg="black", command=filecallback); button.pack(side=LEFT, padx = 2,pady = 2)                    

                #print option,run_mappfinder, self.title, self.default_option
                if len(notes)>0: ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)

            if ('update-entry' in od.DisplayObject()) and self.display_options != ['NA']:
              if use_scroll == 'yes': parent_type = self.sf.interior()
              else: parent_type = self._parent
              proceed = 'yes'
              #if option == 'raw_input': proceed = 'no'
              if proceed == 'yes':
                self._option = option

                #group = PmwFreeze.Group(parent_type,tag_text = self.title)
                #group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
                        
                entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
                try: default_option = defaults[i]
                except Exception: default_option = ''
                entrytxt.set(default_option)
                self.pathdb[option] = entrytxt
                self._user_variables[option] = default_option
                
                #l = Label(parent_type, text=self.title); l.pack(side=LEFT)         
                #entry = Entry(parent_type,textvariable=self.pathdb[option]);
                #entry.pack(side='left',fill = 'both', expand = 1, padx = 10, pady = 2)
                
                l = Label(self.sf.interior(), text=self.title); l.pack()         
                entry = Entry(self.sf.interior(),textvariable=self.pathdb[option]);
                entry.pack()

                #print option,run_mappfinder, self.title, self.default_option
                #if len(notes)>0: ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)
                
            if 'drop-down' in od.DisplayObject() and self.display_options != ['NA']:
                #print option, defaults,  self.display_options
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if insert_into_group == 'yes': parent_type = custom_group.interior(); dropdown_index+=1

                self._option = option
                self.default_option = self.display_options

                def comp_callback1(tag,callback=self.callback,option=option):
                    callback(tag,option)

                self.comp = PmwFreeze.OptionMenu(parent_type,
                    labelpos = 'w', label_text = self.title,                                                 
                    items = self.default_option, command = comp_callback1)

                try: selected_default = od.DefaultOption()
                except Exception:
                    if len(defaults)>0: selected_default = defaults[i]
                    else: selected_default = self.default_option[0] ###Just pick the first option
                #print option, dropdown_index
                if 'species' in option:
                    if 'selected_species2' in option:
                        self.speciescomp2 = self.comp; self.speciescomp2.pack(anchor = 'e', padx = 10, pady = 0, expand = 1, fill = 'both')
                    elif 'selected_species3' in option:
                        self.speciescomp3 = self.comp; self.speciescomp3.pack(anchor = 'e', padx = 10, pady = 0, expand = 1, fill = 'both')
                    else: self.speciescomp = self.comp; self.speciescomp.pack(anchor = 'e', padx = 10, pady = 0, expand = 1, fill = 'both')
                    self.speciescomp.invoke(selected_default) 
                elif 'array_type' in option:
                    self.arraycomp = self.comp; self.arraycomp.pack(anchor = 'e', padx = 10, pady = 0, expand = 1, fill = 'both')
                    self.arraycomp.invoke(selected_default)
                elif 'manufacturer_selection' in option:
                    self.vendorcomp = self.comp; self.vendorcomp.pack(anchor = 'e', padx = 10, pady = 0, expand = 1, fill = 'both')
                    self.vendorcomp.invoke(selected_default)
                else:
                    if insert_into_group == 'no':
                        if 'version' in option: pady_int = 0
                        else: pady_int = 1
                        self.comp.pack(anchor = 'w', padx = 10, pady = pady_int, expand = 1, fill = 'both')
                    elif dropdown_index == 1: comp1 = self.comp
                    elif dropdown_index == 2: comp2 = self.comp
                    elif dropdown_index == 3: comp3 = self.comp
                    elif dropdown_index == 4: comp4 = self.comp
                    elif dropdown_index == 5: comp5 = self.comp
                    elif dropdown_index == 6: comp6 = self.comp
                    elif dropdown_index == 7: comp7 = self.comp
                    elif dropdown_index == 8: comp8 = self.comp
                    elif dropdown_index == 9: comp9 = self.comp
                    elif dropdown_index == 10: comp10 = self.comp
                    elif dropdown_index == 11: comp11 = self.comp
                    elif dropdown_index == 12: comp12 = self.comp
                    elif dropdown_index == 13: comp13 = self.comp
                    try: self.comp.invoke(selected_default)
                    except Exception:
                        #self.comp.invoke(self.display_options[0]) # better to know the variable incase their is a conflict
                        print self.display_options, selected_default, option, option_list;kill
                    if option == 'selected_version':
                        notes = 'Note: Available species may vary based on database selection. Also,\n'
                        notes += 'different Ensembl versions will relate to different genome builds\n'
                        notes += '(e.g., EnsMart54-74 for hg19 and EnsMart75-current for hg38).\n'
                        ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)
                    if option == 'probability_algorithm':
                        notes = 'Note: Moderated tests only run for gene-expression analyses      \n'
                        ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 3)

            if 'comboBox' in od.DisplayObject() and self.display_options != ['NA']:
                
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if insert_into_group == 'yes': parent_type = custom_group.interior(); dropdown_index+=1

                self._option = option
                self.default_option = self.display_options
                
                try: selected_default = od.DefaultOption()
                except Exception:
                    if len(defaults)>0: selected_default = defaults[i]
                    else: selected_default = self.default_option[0] ###Just pick the first option
                    
                listbox_selectmode = 'single'
                if 'multiple' in od.DisplayObject():
                    listbox_selectmode = 'multiple'
                
                def comp_callback1(tag,callback=self.callbackComboBox,option=option):
                    callback(tag,option)
                    
                def mult_callback(tag,callback=self.callbackComboBox,option=option):
                    if 'PathwaySelection' in option:
                        tag = self.pathwayselect.getcurselection() ### there is a conflict otherwise with another multi-comboBox multcomp object
                    elif 'HeatmapAdvanced' in option:
                        tag = self.HeatmapAdvanced.getcurselection() ### there is a conflict otherwise with another multi-comboBox multcomp object
                    else:
                        tag = self.multcomp.getcurselection() ### get the multiple item selection
                    callback(tag,option)
                    
                if 'selected_version' not in option_list: ### For clustering UI
                    label_pos = 'e' ### Orients to the text left -> east
                    entrywidth = 20 ### Width of entry
                    #entry_foreground = 'black'
                    hullsize = 1 #http://pmw.sourceforge.net/doc/ScrolledListBox.html -> doesn't seem to work here
                else:
                    label_pos = 'w' ### Orients to the text right -> west
                    entrywidth = 20 ### Width of entry
                    hullsize = 1
                    
                if listbox_selectmode == 'multiple':
                    self.comp = PmwFreeze.ComboBox(parent_type,
                        labelpos = label_pos, dropdown=1, label_text = self.title,
                        unique = 0, history = 0, entry_background="light gray", entry_width=entrywidth,
                        scrolledlist_usehullsize=1,listbox_selectmode=listbox_selectmode,
                        scrolledlist_items = self.default_option,
                        selectioncommand = mult_callback)
                    self.multcomp = self.comp     
                else:  
                    self.comp = PmwFreeze.ComboBox(parent_type,
                        labelpos = label_pos, dropdown=1, label_text = self.title,
                        unique = 0, history = 0, entry_background="light gray", entry_width=entrywidth,
                        scrolledlist_usehullsize=1,listbox_selectmode=listbox_selectmode,
                        scrolledlist_items = self.default_option,
                        selectioncommand = comp_callback1)
                if 'HeatmapAdvanced' in option:
                    self.HeatmapAdvanced = self.multcomp
                if 'PathwaySelection' in option or 'PathwaySelection_network' in option:
                    if 'network' in option:
                        geneset_param = 'GeneSetSelection_network' ### for network visualization
                    else:
                        geneset_param = 'GeneSetSelection' ### for heatmap visualization
                    self.pathwayselect = self.multcomp; self.pathwayselect.pack(anchor = 'w', padx = 10, pady = 0)
                    try: self.pathwayselect.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.pathwayselect.invoke())
                    except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                    try:
                        ### The next several lines are for a second iteration of this analysis to re-select the previously selected parameters
                        tag = self._user_variables[geneset_param]
                        if 'Ontology' in tag: directory = 'gene-go'
                        else: directory = 'gene-mapp'
                        supported_genesets = self._user_variables[tag]
                        #print 'loading pathways from memory A1'
                        #supported_genesets = listAllGeneSetCategories(species,tag,directory)
                        self.pathwayselect._list.setlist(supported_genesets)
                        self.pathwayselect.selectitem(selected_default) 
                        self.callbackComboBox(selected_default,option)
                    except Exception:
                        try:
                            self.pathwayselect.selectitem(self.default_option[-1]) ###Just pick the first option
                            self.callbackComboBox(self.default_option[-1],option)
                        except Exception: pass
                            
                if 'species' in option:
                    if 'selected_species2' in option:
                        self.speciescomp2 = self.comp; self.speciescomp2.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp2.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp2.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp2.selectitem(selected_default) 
                            self.callbackComboBox(selected_default,option) 
                        except Exception:
                            self.speciescomp2.selectitem(self.default_option[0]) ###Just pick the first option
                            self.callbackComboBox(self.default_option[0],option)
                    elif 'selected_species3' in option:
                        self.speciescomp3 = self.comp; self.speciescomp3.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp3.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp3.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp3.selectitem(selected_default) ###Just pick the first option
                            self.callbackComboBox(selected_default,option) 
                        except Exception:
                            self.speciescomp3.selectitem(self.default_option[0])
                            self.callbackComboBox(self.default_option[0],option)
                    else:
                        self.speciescomp = self.comp; self.speciescomp.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp.selectitem(selected_default)
                            self.callbackComboBox(selected_default,option) 
                        except Exception:
                            self.speciescomp.selectitem(self.default_option[0])
                            self.callbackComboBox(self.default_option[0],option)
                elif 'array_type' in option:
                    self.arraycomp = self.comp; self.arraycomp.pack(anchor = 'w', padx = 10, pady = 0)
                    try: self.arraycomp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.arraycomp.invoke())
                    except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                    try:
                        self.arraycomp.selectitem(selected_default)
                        self.callbackComboBox(selected_default,option) 
                    except Exception:
                        self.arraycomp.selectitem(self.default_option[0])
                        self.callbackComboBox(self.default_option[0],option)
                elif 'manufacturer_selection' in option:
                    self.vendorcomp = self.comp; self.vendorcomp.pack(anchor = 'w', padx = 10, pady = 0)
                    try: self.vendorcomp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.vendorcomp.invoke())
                    except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                    try:
                        self.vendorcomp.selectitem(selected_default)
                        self.callbackComboBox(selected_default,option) 
                    except Exception:
                        self.vendorcomp.selectitem(self.default_option[0])
                        self.callbackComboBox(self.default_option[0],option)
                else:
                    self.combo = self.comp ### has to be a unique combo box to refer to itself in the component call below
                    self.combo.pack(anchor = 'w', padx = 10, pady = 1)
                    try: self.combo.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.combo.invoke())
                    except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                    """
                    if listbox_selectmode == 'multiple':
                        if len(od.DefaultOption()[0])>1: ###Hence it is a list
                            self.combo.ApplyTypeSelections(od.DefaultOption()) 
                            for opt in od.DefaultOption():
                                self.combo.invoke(opt)
                            self.callbackComboBox(tuple(od.DefaultOption()),option) 
                            #self.combo.selectitem(opt)
                            #self.callbackComboBox(opt,option)
                    """
                    #print selected_default
                    try:
                        if len(selected_default[0])>1: ###Hence it is a list
                            for opt in selected_default:
                                self.combo.selectitem(opt)
                                self.callbackComboBox(opt,option); break
                        else:
                            ### This is where the default for the combobox is actually selected for GeneSets
                            self.combo.selectitem(selected_default)
                            self.callbackComboBox(selected_default,option) 
                    except Exception:
                        try:
                            self.combo.selectitem(self.default_option[0])
                            self.callbackComboBox(self.default_option[0],option)
                        except Exception:
                            None
                    if option == 'selected_version':
                        notes = 'Note: Available species may vary based on database selection. Also,\n'
                        notes += 'different Ensembl versions will relate to different genome builds\n'
                        notes += '(e.g., EnsMart54-74 for hg19 and EnsMart75-current for hg38).\n'
                        ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)
                        
            if 'pulldown_comps' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = self.display_options
                ###From the option, create two new options, one for each group in the comparison
                option1 = option+'-1'; option2 = option+'-2'
                ### Pack these into a groups to maintain organization
                group = PmwFreeze.Group(self.sf.interior(),tag_text = self.title)
                group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)

                if check_index == -1:
                    check_option = 'analyze_all_conditions'
                    def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=check_option):
                        #print tag,state,option
                        checkbuttoncallback(tag,state,option)                
                    ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                    self.checkbuttons = PmwFreeze.RadioSelect(self._parent,
                            buttontype = 'checkbutton', command = checkbuttoncallback)
                    self.checkbuttons.pack(side = 'top', expand = 1, padx = 0, pady = 0)
                    ### Add some buttons to the checkbutton RadioSelect.
                    self.checkbuttons.add('Analyze ALL GROUPS in addition to specifying comparisons')
                    self._user_variables[check_option] = 'no'
                    check_index+=1
                    
                def comp_callback1(tag,callback=self.callback,option1=option1):
                    callback(tag,option1)
                def comp_callback2(tag,callback=self.callback,option2=option2):
                    callback(tag,option2)

                #labelpos = 'w', label_text = self.title,  -inside of OptionMenu
                self.comp1 = PmwFreeze.OptionMenu(group.interior(),
                    items = self.default_option, menubutton_width = 20, command = comp_callback1)
                self.comp1.pack(side = LEFT, anchor = 'w', padx = 10, pady = 0)

                self.comp2 = PmwFreeze.OptionMenu (group.interior(), 
                    items = self.default_option, menubutton_width = 20, command = comp_callback2,
                ); self.comp2.pack(side = LEFT, anchor = 'w', padx = 10, pady = 0)

                try: self.comp1.invoke(notes[0])
                except Exception: pass
                try: self.comp2.invoke(notes[1])
                except Exception: pass
                
            if 'simple_entry' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                if len(override_default)>0: self.default_option = override_default
                else: self.default_option = self.display_options[0]
                def enter_callback(tag,enter_callback=self.enter_callback,option=option):
                    enter_callback(tag,option)
                #self.title = self.title + '\t ' #entry_width=entrywidth
                self.entry_field = PmwFreeze.EntryField(self.sf.interior(),
                        labelpos = 'e', label_text = self.title,
                        validate = enter_callback,
                        value = self.default_option
                ); self.entry_field.pack(anchor='w',padx = 10, pady = 1)

            if 'enter' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if insert_into_group == 'yes': parent_type = custom_group.interior(); enter_index+=1
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                if len(override_default)>0: self.default_option = override_default
                elif len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                #print self.default_option, self.title; kill
                ### entrytxt object for alt_exon_fold_cutoff in option
                def custom_validate(tag,custom_validate=self.custom_validate,option=option):
                    validate = custom_validate(tag,option)
                def custom_validate_p(tag,custom_validate_p=self.custom_validate_p,option=option):
                    validate = custom_validate_p(tag,option)
                    #print [validate], tag, option
                if 'Genes_network' in option_list:
                    label_pos = 'e' ### Orients to the text left -> east
                    entrywidth = 20 ### Width of entry entry_width=entrywidth          
                elif 'GeneSelection' in option_list or 'GeneSelectionPredict' in option_list: ### For clustering UI
                    label_pos = 'e' ### Orients to the text left -> east
                    entrywidth = 20 ### Width of entry entry_width=entrywidth
                elif 'JustShowTheseIDs' in option_list or 'JustShowTheseIDsPredict' in option_list: ### For clustering UI
                    label_pos = 'e' ### Orients to the text left -> east
                    entrywidth = 20 ### Width of entry entry_width=entrywidth
                else:
                    label_pos = 'e'
                try:
                    if float(self.default_option) <= 1: use_method = 'p'
                    else: use_method = 'i'
                except ValueError:
                    #self.default_option = 'CHANGE TO A NUMERIC VALUE'; use_method = 'i'
                    self.default_option = string.replace(self.default_option,'---','')
                    use_method = 'i'

                if use_method == 'p':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = label_pos, label_text = self.title, validate = custom_validate_p, 
                            value = self.default_option, hull_borderwidth = 1)                   
                if use_method == 'i':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = label_pos, label_text = self.title, validate = custom_validate,
                            value = self.default_option, hull_borderwidth = 1)
                
                #if 'GeneSelection' in option_list:
                #self.entry_field.component("entry").configure(width=5)
                
                if insert_into_group == 'no': self.entry_field.pack(anchor = 'w', padx = 10, pady = 0)	
                elif enter_index == 1: self.entry_field1 = self.entry_field
                elif enter_index == 2: self.entry_field2 = self.entry_field
                elif enter_index == 3: self.entry_field3 = self.entry_field
                elif enter_index == 4: self.entry_field4 = self.entry_field
                elif enter_index == 5: self.entry_field5 = self.entry_field
                if len(notes)>0: Label(self._parent, text=notes).pack()

            if 'multiple-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                self._option = option
                if len(override_default)>0: self.default_option = override_default
                elif len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=option):
                    checkbuttoncallback(tag,state,option)                 
                ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                        buttontype = 'checkbutton', orient = 'vertical',
                        labelpos = 'w', command = self.checkbuttoncallback,
                        label_text = self.title, hull_borderwidth = 2)
                self.checkbuttons.pack(padx = 10, pady = 0)

                ### Add some buttons to the checkbutton RadioSelect.
                for text in self.display_options:
                     if text != ['NA']:
                        self.checkbuttons.add(text)
                        #if 'common-' not in text and 'all-' not in text:
                        #self.checkbuttons.invoke(text)
                if len(notes)>0: Label(self._parent, text=notes).pack()
            
            if 'single-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if defaults == 'comps': parent_type = self._parent; orient_type = 'top'
                if insert_into_group == 'yes': parent_type = custom_group.interior(); check_index+=1
                if defaults == 'groups': parent_type = self.sf.interior(); orient_type = 'top'
                
                self._option = option
                proceed = 'yes'
                """if option == 'export_splice_index_values':
                    if analysis_method != 'splicing-index': proceed = 'no' ### only export corrected constitutive ratios if splicing index method chosen"""
                if proceed == 'yes':
                    if len(override_default)>0: self.default_option = override_default
                    elif len(defaults) <1: self.default_option = self.display_options[0]
                    else: self.default_option = defaults[i]
                    if self.default_option != 'NA':
                        def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=option):
                            checkbuttoncallback(tag,state,option)                
                        ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                        self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                                buttontype = 'checkbutton', command = checkbuttoncallback)
                                #hull_borderwidth = 2, hull_relief = 'ridge')
                        if insert_into_group == 'no': self.checkbuttons.pack(anchor = 'w',side = orient_type, padx = 10, pady = 1)
                        elif check_index == 1:  checkbuttons1 = self.checkbuttons
                        elif check_index == 2:  checkbuttons2 = self.checkbuttons
                        elif check_index == 3:  checkbuttons3 = self.checkbuttons
                        elif check_index == 4:  checkbuttons4 = self.checkbuttons
                        ### Add some buttons to the checkbutton RadioSelect.
                        self.checkbuttons.add(self.title)
                        if self.default_option == 'yes': self.checkbuttons.invoke(self.title)
                        else: self._user_variables[option] = 'no'

            custom_group_endpoints = ['ge_ptype', 'get_additional', 'expression_threshold', 'run_goelite', 'gene_expression_cutoff', 'microRNA_prediction_method']
            try:
                eod = option_db['expression_threshold']
                if eod.ArrayOptions() == ['NA']:
                    custom_group_endpoints.append('rpkm_threshold') ### Ensures that if analyzing pre-compiled gene expression values, only certain items are shown and in a frame
                    custom_group_endpoints.remove('expression_threshold')
                    #insert_into_group = 'yes'
            except Exception:  pass 
            if option in custom_group_endpoints and insert_into_group == 'yes':
                ### This is employed when we want to place several items into a group frame together.
                ### Since this is a generic class, we need to setup special cases to do this, however,
                ### this same code could be used in other instances as well
                reorganize = 'no'
                self.group_tag = 'GO-Elite Over-Representation and Filtering Parameters'; pady_int = 5
                if 'run_goelite' in option_list: self.group_tag = 'Gene Expression Analysis Options'; pady_int = 1                
                if 'microRNA_prediction_method' in option_list: self.group_tag = 'Advanced Options'; pady_int = 1; reorganize = 'yes'

                try: checkbuttons1.pack(anchor = 'w', side = 'top', padx = 9, pady = 0)
                except Exception: pass
                try: checkbuttons2.pack(anchor = 'w', side = 'top', padx = 9, pady = 0)
                except Exception: pass
                try: checkbuttons3.pack(anchor = 'w', side = 'top', padx = 9, pady = 0)
                except Exception: pass
                try: checkbuttons4.pack(anchor = 'w', side = 'top', expand = 1, padx = 9, pady = 0)
                except Exception: pass
                try: radiobuttons2.pack(side = orient_type, expand = 1, padx = 10, pady = 5)
                except Exception: pass
                
                try: comp1.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: radiobuttons1.pack(side = orient_type, expand = 1, padx = 10, pady = 5)
                except Exception: pass
                if reorganize == 'yes':
                    try: comp2.pack(anchor = 'w', padx = 10, pady = pady_int)
                    except Exception: pass
                try: self.entry_field1.pack(anchor = 'w', padx = 10, pady = 0)
                except Exception: pass
                try: self.entry_field2.pack(anchor = 'w', padx = 10, pady = 0);
                except Exception: pass
                try: self.entry_field3.pack(anchor = 'w', padx = 10, pady = 0)
                except Exception: pass
                try: self.entry_field4.pack(anchor = 'w', padx = 10, pady = 0)
                except Exception: pass
                try: self.entry_field5.pack(anchor = 'w', padx = 10, pady = 0)
                except Exception: pass
                if reorganize == 'no':
                    try: comp2.pack(anchor = 'w', padx = 10, pady = pady_int)
                    except Exception: pass
                try: comp3.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp4.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp5.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp6.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp7.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp8.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp9.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp10.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp11.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp12.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                try: comp13.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: pass
                enter_index=0; radio_index=0; dropdown_index=0
                if 'ge_ptype' in option or 'expression_threshold' in option or 'gene_expression_cutoff' in option or 'rpkm_threshold' in option:
                    custom_group = PmwFreeze.Group(self.sf.interior(),tag_text = self.group_tag)
                    custom_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
                    insert_into_group = 'yes'
                
            #i+=1 ####Keep track of index
            
        if len(option_list)>0: ### Used when visualizing WikiPathways (no option_list supplied - all parameters hard coded)
            #def quitcommand(): parent.destroy; sys.exit()
            #self.button = Button(text="   Quit  ", command=quitcommand)
            #self.button.pack(side = 'bottom', padx = 10, pady = 10)
    
            if 'input_cdf_file' in option_list: ### For the CEL file selection window, provide a link to get Library files
                button_text = 'Download Library Files'; d_url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
                self.d_url = d_url; text_button = Button(self._parent, text=button_text, command=self.Dlinkout); text_button.pack(side = 'left', padx = 5, pady = 5)

            if 'GeneSelectionPredict' in option_list:
                run_button = Button(self._parent, text='Run Analysis', command=self.runPredictGroups)
                run_button.pack(side = 'right', padx = 10, pady = 10)
            else:
                continue_to_next_win = Button(text = 'Continue', command = self._parent.destroy)
                continue_to_next_win.pack(side = 'right', padx = 10, pady = 10) 
    
            if 'input_annotation_file' in option_list:
                skip_win = Button(text = 'Skip', command = self._parent.destroy)
                skip_win.pack(side = 'right', padx = 10, pady = 10)
                
            back_button = Button(self._parent, text="Back", command=self.goBack) 
            back_button.pack(side = 'right', padx =10, pady = 5)
            
            quit_win = Button(self._parent, text="Quit", command=self.quit) 
            quit_win.pack(side = 'right', padx =10, pady = 5)
            
            button_text = 'Help'
            url = 'http://www.altanalyze.org/help_main.htm'; self.url = url
            pdf_help_file = 'Documentation/AltAnalyze-Manual.pdf'; pdf_help_file = filepath(pdf_help_file); self.pdf_help_file = pdf_help_file
            
            try: help_button = Button(self._parent, text=button_text, command=self.GetHelpTopLevel); help_button.pack(side = 'left', padx = 5, pady = 5)
            except Exception: help_button = Button(self._parent, text=button_text, command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)
 
            if 'species' in option_list:
                new_species_button = Button(self._parent, text='Add New Species', command=self.newSpecies)
                new_species_button.pack(side = 'left', padx = 5, pady = 5)
    
            def runPredictGroupsTest():
                self.runPredictGroups(reportOnly=True)
            
            if 'GeneSelectionPredict' in option_list:
                expFilePresent = self.verifyExpressionFile()
                if expFilePresent:
                    button_instance = Button(self._parent, text='Test Settings', command=runPredictGroupsTest)
                    button_instance.pack(side = 'left', padx = 5, pady = 5)
                    
            if 'build_exon_bedfile' in option_list and array_type == 'RNASeq':
                self.pdf_help_file = filepath('AltDatabase/kallisto/license.txt')
                button_instance = Button(self._parent, text='Kallisto License', command=self.openPDFHelp)
                button_instance.pack(side = 'left', padx = 5, pady = 5)
                
            self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
            self._parent.mainloop()
        
    def verifyExpressionFile(self):
        continue_analysis = False ### See if the input file is already present
        try:
            expFile = fl.ExpFile()
            count = verifyFileLength(expFile[:-4]+'-steady-state.txt')
            if count>1: continue_analysis = True
            else:
                count = verifyFileLength(expFile)
                if count>1: continue_analysis = True
        except Exception: pass
        return continue_analysis
    
    def goBack(self):
        self._parent.destroy()
        selected_options = selected_parameters; selected_options2=[] ### If we don't do this we get variable errors
        if 'Library' == selected_options[-1]: selected_options[-1] = ''
        for i in selected_options:
            if i!='Library': selected_options2.append(i)
        selected_options = selected_options2
        try:
            while selected_options[-2]==selected_options[-1]:
                selected_options = selected_options[:-1] ### When clicking back on the next loop of a back, makes sure you don't get looped back to the same spot
        except Exception: selected_options = selected_options
        if len(selected_options)<3: run_parameter = 'no'
        else: run_parameter = selected_options[:-1], self._user_variables
        AltAnalyze.AltAnalyzeSetup(run_parameter); sys.exit()
        
    def newSpecies(self):
        self._user_variables['species'] = 'Add Species'
        self._parent.destroy()
        
    def runPredictGroups(self,reportOnly=False):
        column_metric = self.Results()['column_metric_predict']
        column_method = self.Results()['column_method_predict']
        GeneSetSelection = self.Results()['GeneSetSelectionPredict']
        try: PathwaySelection = self.Results()['PathwaySelectionPredict']
        except Exception: PathwaySelection = 'None Selected'
        GeneSelection = self.Results()['GeneSelectionPredict']
        JustShowTheseIDs = self.Results()['JustShowTheseIDsPredict']
        ExpressionCutoff = self.Results()['ExpressionCutoff']
        CountsCutoff = self.Results()['CountsCutoff']
        rho_cutoff = self.Results()['rho_cutoff']
        FoldDiff = self.Results()['FoldDiff']
        SamplesDiffering = self.Results()['SamplesDiffering']
        try: featurestoEvaluate = self.Results()['featuresToEvaluate']
        except Exception: featurestoEvaluate = 'Genes'
        removeOutliers = self.Results()['removeOutliers']
        dynamicCorrelation = self.Results()['dynamicCorrelation']
        restrictBy = self.Results()['restrictBy']
        excludeCellCycle = self.Results()['excludeCellCycle']
        gsp = GeneSelectionParameters(species,array_type,vendor)
        gsp.setGeneSet(GeneSetSelection)
        gsp.setPathwaySelect(PathwaySelection)
        gsp.setGeneSelection(GeneSelection)
        gsp.setJustShowTheseIDs(JustShowTheseIDs)
        gsp.setNormalize('median')
        gsp.setSampleDiscoveryParameters(ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering,dynamicCorrelation,
            removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff)
        self._user_variables['gsp'] = gsp
        
        import RNASeq  
        expFile = fl.ExpFile()
        mlp_instance = fl.MLP()
        count = verifyFileLength(expFile[:-4]+'-steady-state.txt')
        if count>1: expFile = expFile[:-4]+'-steady-state.txt'
        
        if reportOnly:
            ### Only used to report back what the number of regulated genes are if the gene expression file is present
            reload(RNASeq)
            try: report = RNASeq.singleCellRNASeqWorkflow(species, array_type, expFile, mlp_instance, parameters=gsp, reportOnly=reportOnly)
            except Exception: report = traceback.format_exc()
            if 'options_result_in_no_genes' in report:
                report = 'Options are too stringent. Try relaxing the thresholds.'
            try: InfoWindow(report, 'Continue')
            except Exception: print report
        else:
            ### If the parameters look OK, or user wishes to run, collapse this GUI adn proceed (once exited it will run)
            self._parent.quit()
            self._parent.destroy()
            """
            values = expFile, mlp_instance, gsp, reportOnly
            StatusWindow(values,'predictGroups') ### display an window with download status
            root = Tk()
            root.title('AltAnalyze: Evaluate Sampled Groupings')
            gu = GUI(root,'PredictGroups',[],'')
            nextStep = gu.Results()['next']
            group_selected = gu.Results()['group_select']
            if nextStep == 'UseSelected':
                print group_selected;sys.exit()
                group_selected = group_selected
                ### When nothing returned here, the full analysis will run
            else:
                #print 're-initializing window'
                AltAnalyze.AltAnalyzeSetup((selected_parameters,user_variables)); sys.exit()
                
            """
            
    def GetHelpTopLevel(self):
        message = ''
        self.message = message; self.online_help = 'Online Documentation'; self.pdf_help = 'Local PDF File'
        tl = Toplevel(); self._tl = tl; nulls = '\t\t\t\t'; tl.title('Please select one of the options')

        self.sf = PmwFreeze.ScrolledFrame(self._tl,
                labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 320, hull_height = 200)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()

        group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Options')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)

        filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(group.interior()); can.pack(side='left',padx = 10, pady = 20); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)

        l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
        text_button2 = Button(group.interior(), text=self.online_help, command=self.openOnlineHelp); text_button2.pack(side = 'top', padx = 5, pady = 5) 
        try: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        except Exception: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        tl.mainloop()

    def openPDFHelp(self):
        if runningCommandLine:
            pass
        elif os.name == 'nt':
            try: os.startfile('"'+self.pdf_help_file+'"')
            except Exception:  os.system('open "'+self.pdf_help_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.pdf_help_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.pdf_help_file+'"')
        if 'license' not in self.pdf_help_file:
            try: self._tl.destroy()
            except Exception: pass

    def openOnlineHelp(self):
        try: webbrowser.open(self.url)
        except Exception: pass
        self._tl.destroy()
    
    def linkout(self):
        try: webbrowser.open(self.url)
        except Exception: pass
        
    def Dlinkout(self):
        try: webbrowser.open(self.d_url)
        except Exception: pass
            
    def setvscrollmode(self, tag):
        self.sf.configure(vscrollmode = tag)

    def info(self):
        tkMessageBox.showinfo("title","message",parent=self._parent)

    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()

    def tldeleteWindow(self,tl):
        try: tl.quit(); tl.destroy()#; print 1
        except Exception: tl.destroy()#; print 2
    
    def quit(self):
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: self._parent.quit(); sys.exit()

    def continue_win(self):
        ### Currently not used - can be used to check data integrity before closing a window
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: self._parent.quit(); sys.exit()
        
    def chooseDirectory(self,option):
        tag = tkFileDialog.askdirectory(parent=self._parent)
        self._user_variables[option] = tag

    def chooseFile(self,option):
        tag = tkFileDialog.askopenfile(parent=self._parent)
        self._user_variables[option] = tag.name
        
    def getPath(self,option):
        
        try: ### Below is used to change a designated folder path to a filepath
            if option == 'input_cel_dir' and '10X' in vendor:
                processFile = True
            else:
                forceException
        except Exception:
            #print traceback.format_exc()
            if 'dir' in option or 'folder' in option:
                processFile = False
            else:
                processFile = True
        #print option, processFile
        if ('dir' in option or 'folder' in option) and processFile ==False:
            try: dirPath = tkFileDialog.askdirectory(parent=self._parent,initialdir=self.default_dir)
            except Exception: 
                self.default_dir = ''
                try: dirPath = tkFileDialog.askdirectory(parent=self._parent,initialdir=self.default_dir)
                except Exception: 
                    try: dirPath = tkFileDialog.askdirectory(parent=self._parent)
                    except Exception:
                        print_out = "AltAnalyze is unable to initialize directory opening.\nThis error may be due to one of the following issues:\n"
                        print_out += "1) (if running directly from source code) Tkinter/PMW components are not installed or are incompatible.\n"
                        print_out += "2) (if running directly from source code) Python version is untested with AltAnalyze.\n"
                        print_out += "3) There is a conflict between the AltAnalyze packaged version of python on the OS.\n\n"
                        print_out += "Contact altanalyze@gmail.com if this error persists with your system information.\n"
                        print_out += "Alternatively, try AltAnalyze using command-line options (http://www.AltAnalyze.org)."
                        try: InfoWindow(print_out,'Continue')
                        except Exception: print print_out
                        try: self._parent.destroy(); sys.exit()
                        except Exception: sys.exit()
            self.default_dir = dirPath
            entrytxt = self.pathdb[option]; entrytxt.set(dirPath)
            self._user_variables[option] = dirPath
            try: file_location_defaults['PathDir'].SetLocation(dirPath) ### Possible unknown exception here... may need to correct before deployment
            except Exception:
                try:
                    ### Entry was deleted from Config file - re-create it
                    fl = FileLocationData('local', dirPath, 'all')
                    file_location_defaults['PathDir'] = fl
                except Exception: pass
            try: exportDefaultFileLocations(file_location_defaults)
            except Exception: pass

        if 'file' in option or processFile:
            try: tag = tkFileDialog.askopenfile(parent=self._parent,initialdir=self.default_file)
            except Exception: 
                self.default_file = ''
                try: tag = tkFileDialog.askopenfile(parent=self._parent,initialdir=self.default_file)
                except Exception: 
                    try: tag = tkFileDialog.askopenfile(parent=self._parent)
                    except Exception:
                        print_out = "AltAnalyze is unable to initialize directory opening.\nThis error may be due to one of the following issues:\n"
                        print_out += "1) (if running directly from source code) Tkinter/PMW components are not installed or are incompatible.\n"
                        print_out += "2) (if running directly from source code) Python version is untested with AltAnalyze.\n"
                        print_out += "3) There is a conflict between the AltAnalyze packaged version of python on the OS.\n\n"
                        print_out += "Contact altanalyze@gmail.com if this error persists with your system information.\n"
                        print_out += "Alternatively, try AltAnalyze using command-line options (see documentation at http://AltAnalyze.org)."
                        try: InfoWindow(print_out,'Continue')
                        except Exception: print print_out
                        try: self._parent.destroy(); sys.exit()
                        except Exception: sys.exit()
            try: filePath = tag.name #initialdir=self.default_dir
            except AttributeError: filePath = ''
            filePath_dir = string.join(string.split(filePath,'/')[:-1],'/')
            self.default_file = filePath_dir
            entrytxt = self.pathdb[option]
            entrytxt.set(filePath)
            self._user_variables[option] = filePath
            try: file_location_defaults['PathFile'].SetLocation(filePath_dir)
            except Exception:
                try:
                    ### Entry was deleted from Config file - re-create it
                    fl = FileLocationData('local', filePath_dir, 'all')
                    file_location_defaults['PathFile'] = fl
                except Exception: null = None
            try: exportDefaultFileLocations(file_location_defaults)
            except Exception: pass
        
    def Report(self,tag,option):
        output = tag
        return output
    def __repr__(self,tag,option): return self.Report(tag,option)
    
    def Results(self):
        for i in self._user_variables:
            user_variables[i]=self._user_variables[i]
        return self._user_variables
    
    def custom_validate(self, text, option):
        self._user_variables[option] = text
        #try: text = float(text);return 1
        #except ValueError: return -1

    def enter_callback(self, tag, option):
        if self.defaults == 'batch':
            ### Bath removal array annotation UI
            self._user_variables[option,'batch'] = tag
        else:
            self._user_variables[option] = tag

    def custom_validate_p(self, text, option):
        #print [option],'text:', text
        self._user_variables[option] = text
        try:
            text = float(text)
            if text <1:return 1
            else:return -1
        except ValueError:return -1
        
    def callback(self, tag, option):
        #print 'Button',[option], tag,'was pressed.'
        change_var = ''
        self._user_variables[option] = tag
        if option == 'dbase_version':
            ###Export new species info
            exportDBversion(tag); change_var = 'all'
            try: self.changeVendorSelection(); self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: pass
        elif option == 'species':
            try: self.changeArraySelection()
            except Exception: pass
        elif option == 'manufacturer_selection':
            try: self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: pass
        #elif option == 'array_type':
            #self.checkSpeciesArraySelection(array_type)
        elif option == 'analysis_method':
            if tag == 'ASPIRE':
                try: self.entry_field2.setentry('0.2')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '0.2'
            elif tag == 'linearregres':
                try: self.entry_field2.setentry('2')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '2'
            elif tag == 'MultiPath-PSI':
                try: self.entry_field2.setentry('0.1')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '0.1'
        elif option == 'selected_version':
            current_species_names = db_versions[tag]
            current_species_names.sort()
            try: self.speciescomp.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp2.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp3.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        """
        ### Doesn't work right now because self.entry_field only has one object instance and complicated to get multiple
        elif option == 'ORA_algorithm':
            if tag == 'Permute p-value':
                try: self.entry_field.setentry('2000')
                except Exception: pass
                self._user_variables['permutation'] = '2000'
            elif tag == 'Fisher Exact Test':
                try: self.entry_field.setentry('NA')
                except Exception: pass
                self._user_variables['permutation'] = '0'
        """
    def callbackComboBox(self, tag, option):
        """ Similiar to the above, callback, but ComboBox uses unique methods """
        #print 'Button',[option], tag,'was pressed.'
        if option == 'interactionDirs' or 'PathwaySelection' in option or 'HeatmapAdvanced' in option: ### Allow multiple selections
            if len(tag)==0:
                self._user_variables[option] = None ### no options selected
            else:
                if isinstance(tag, tuple) or isinstance(tag, list):
                    pass
                else:
                    try: tag = self._user_variables[option] ### This indicates that this option was previously set and in the new window was not explicitly set, suggesting we should re-apply the original settings
                    except Exception: None
                try: ### Occurs when no items selected
                    if len(tag[0])==1: ### Hence, just one item selected
                        self._user_variables[option] = [tag]
                except Exception:
                    pass
                if len(list(tag)[0]) == 1:
                    tag_list = [tag]
                else: tag_list = list(tag)
                self._user_variables[option] = tag_list
        else:
            self._user_variables[option] = tag
        if option == 'selected_version':
            current_species_names = db_versions[tag]
            current_species_names.sort()
            current_species_names = ['---']+current_species_names
            species_option = current_species_names[0]
            try:
                self.speciescomp._list.setlist(current_species_names) ### This is the way we set a new list for ComboBox
                ### Select the best default option to display (keep existing or re-set)
                if 'selected_species1' in self._user_variables: ### If this is the species downloader
                    species_option = 'selected_species1'
                else:
                    for i in self._user_variables:
                        if 'species' in i: species_option = i
                default = self.getBestDefaultSelection(species_option,current_species_names)
                self.speciescomp.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try:
                self.speciescomp2._list.setlist(current_species_names)
                default = self.getBestDefaultSelection('selected_species2',current_species_names)
                self.speciescomp2.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try:
                self.speciescomp3._list.setlist(current_species_names)
                default = self.getBestDefaultSelection('selected_species3',current_species_names)
                self.speciescomp3.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        elif option == 'dbase_version':
            ###Export new species info
            exportDBversion(tag); change_var = 'all'
            try: self.changeVendorSelection(); self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: pass
        elif option == 'species':
            try: self.changeArraySelection()
            except Exception: pass
        elif option == 'manufacturer_selection':
            try: self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: pass
        #elif option == 'array_type':
            #self.checkSpeciesArraySelection(array_type)
        elif option == 'analysis_method':
            if tag == 'ASPIRE':
                try: self.entry_field2.setentry('0.2')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '0.2'
            elif tag == 'linearregres':
                try: self.entry_field2.setentry('2')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '2'
            elif tag == 'MultiPath-PSI':
                try: self.entry_field2.setentry('0.1')
                except Exception: pass
                self._user_variables['alt_exon_fold_cutoff'] = '0.1'
        elif 'GeneSetSelection' in option or 'GeneSetSelectionPredict' in option:
            #print option,tag
            if 'network' in option: suffix='_network'
            else: suffix=''
            #species = self._user_variables['species']
            try:
                if 'Ontology' in tag: directory = 'gene-go'
                else: directory = 'gene-mapp'
                if tag in self._user_variables and 'StoredGeneSets' not in tag: ### Want to reload StoredGeneSets each time
                    supported_genesets = self._user_variables[tag]
                    #print 'loading pathways from memory'
                else:
                    #print 'importing all pathways from scratch'
                    supported_genesets = listAllGeneSetCategories(species,tag,directory)
                    self._user_variables[tag] = supported_genesets ### Store this so we don't waste time reloading it the next time
                self.pathwayselect._list.setlist(supported_genesets)
                ##### self.pathwayselect.selectitem(supported_genesets[0]) # This sets the default for multi- or single-combo boxes... DON'T SELECT UNLESS YOU WANT TO HAVE TO DE-SELECT IT
                ##### self._user_variables['PathwaySelection'+suffix] = supported_genesets[0] ### store this default
                ##### self._user_variables['PathwaySelectionPredict'+suffix] = supported_genesets[0] ### store this default
                self.pathwayselect.selectitem(0,setentry = 1) # Select the item but then re-set the list to deselect it
                self.pathwayselect._list.setlist(supported_genesets)
            except Exception, e:
                #print e
                pass ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            
    def getBestDefaultSelection(self,option,option_list):
        default = option_list[0] ### set the default to the first option listed
        if option in self._user_variables:
            selected = self._user_variables[option]
            if selected in option_list: ### If selected species exists in the new selected version of EnsMart
                default = selected
            else:
                self._user_variables[option] = default ### Hence, the default has changed, so re-set it
        return default
    
    def changeSpeciesSelection(self):
        vendor = self._user_variables['manufacturer_selection'] ### Get vendor (stored as global)
        current_species_names = getSpeciesList(vendor) ### Need to change species, manufacturers and array_type
        for i in self._option_list:
            if 'species' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                last_selected_species = self._user_variables[i]
                if last_selected_species not in current_species_names:
                    try: self._user_variables[i] = current_species_names[0]
                    except Exception: null = []
                    try:
                        self.speciescomp._list.setlist(current_species_names)
                        self.speciescomp.selectitem(current_species_names[0])
                    except Exception:
                        print traceback.format_exc()
                        null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated

    def checkSpeciesArraySelection(self,array_type):
        current_species_names = getSpeciesForArray(array_type)
        try:
            self.speciescomp._list.setlist(current_species_names)
            self.speciescomp.selectitem(current_species_names[0])
        except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'species' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                try: self._user_variables[i] = current_species_names[0]
                except Exception: null = []
        
    def changeArraySelection(self):
        species_name = self._user_variables['species'] ### Get species (stored as global)
        vendor = self._user_variables['manufacturer_selection'] ### Get vendor (stored as global)
        species = species_codes[species_name].SpeciesCode()
        current_array_types, manufacturer_list = getArraysAndVendors(species,vendor)
        if 'Other ID'==vendor: ### Populate the current_array_types as all Ensembl linked systems
            current_array_types = getSupportedGeneSystems(species,'uid-gene')         
        try:
            self.arraycomp._list.setlist(current_array_types)
            self.arraycomp.selectitem(current_array_types[0])
        except Exception:
            pass ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'array_type' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                if self._user_variables[i] not in current_array_types: ### If the current array type is supported by the new species selection, keep it the same
                    try: self._user_variables[i] = current_array_types[0]
                    except Exception: null = []

    def changeVendorSelection(self):
        species_name = self._user_variables['species'] ### Get species (stored as global)
        vendor = self._user_variables['manufacturer_selection']
        current_array_types, manufacturer_list = getArraysAndVendors(species,'')
        try:
            self.vendorcomp._list.setlist(manufacturer_list)
            self.vendorcomp.selectitem(manufacturer_list[0])
        except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'manufacturer_selection' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                if vendor in manufacturer_list: new_vendor = vendor
                else: new_vendor = manufacturer_list[0]
                try: self._user_variables[i] = new_vendor
                except Exception: null = []
                
    def multcallback(self, tag, state):
        if state: action = 'pressed.'
        else: action = 'released.'
        """print 'Button', tag, 'was', action, \
                'Selection:', self.multiple.getcurselection()"""
        self._user_variables[option] = tag

    def checkbuttoncallback(self, tag, state, option):
        if state: action = 'pressed.'
        else: action = 'released.'
        """print 'Button',[option], tag, 'was', action, \
                'Selection:', self.checkbuttons.getcurselection()"""
        if state==0: tag2 = 'no'
        else: tag2 = 'yes'
        #print '---blahh', [option], [tag], [state], [action], [self.checkbuttons.getcurselection()]
        self._user_variables[option] = tag2

 ################# Database Version Handling ##################

class PreviousResults:
    def __init__(self, user_variables): 
        self._user_variables = user_variables
    def Results(self): return self._user_variables
        
def getSpeciesList(vendor):
    try: current_species_dirs = unique.read_directory('/AltDatabase')
    except Exception: ### Occurs when the version file gets over-written with a bad directory name
        try:
            ### Remove the version file and wipe the species file
            os.remove(filepath('Config/version.txt'))
            #raw = export.ExportFile('Config/species.txt'); raw.close()
            os.mkdir(filepath('AltDatabase'))
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        except Exception:
            #print traceback.format_exc()
            print 'Cannot write Config/version.txt to the Config directory (likely Permissions Error)'
        try: db_versions = returnDirectoriesNoReplace('/AltDatabase')
        except Exception:
            try: os.mkdir(filepath('AltDatabase'))
            except Exception: pass
            db_versions = returnDirectoriesNoReplace('/AltDatabase')
        for db_version in db_versions:
            if 'EnsMart' in db_version:
                try: exportDBversion(db_version)
                except Exception: exportDBversion('')
                break
        current_species_dirs = unique.read_directory('/AltDatabase')

    current_species_names=[]; manufacturers_list=[]
    for species in species_codes:
        species_code = species_codes[species].SpeciesCode()
        if species_code in current_species_dirs:
            if len(vendor)>0:
                proceed = 'no'
                for array_name in array_codes:
                    manufacturer = array_codes[array_name].Manufacturer()
                    if manufacturer == vendor:
                        if species_code in array_codes[array_name].SpeciesCodes(): proceed = 'yes'
            else:
                for array_name in array_codes:
                    manufacturer = array_codes[array_name].Manufacturer()
                    if species_code in array_codes[array_name].SpeciesCodes():
                        manufacturers_list.append(manufacturer)
                proceed = 'yes'
            if proceed == 'yes': current_species_names.append(species)
    current_species_names.sort(); manufacturers_list = unique.unique(manufacturers_list); manufacturers_list.sort()
    if len(vendor)>0:
        return current_species_names
    else: return current_species_names, manufacturers_list

def exportDBversion(db_version):
    import datetime
    db_version = string.replace(db_version,'Plant','')
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    try: exportVersionData(db_version,today,'Config/')
    except Exception:
        print traceback.format_exc()
        print 'Cannot write Config/version.txt to the Config directory (likely Permissions Error)'

def exportVersionData(version,version_date,dir):
    new_file = dir+'version.txt'
    new_file_default = unique.filepath(new_file,force='application-path') ### can use user directory local or application local
    try:
        data.write(str(version)+'\t'+str(version_date)+'\n'); data.close()
    except:
        data = export.ExportFile(new_file)
        data.write(str(version)+'\t'+str(version_date)+'\n'); data.close()

def importResourceList():
    filename = 'Config/resource_list.txt'
    fn=filepath(filename); resource_list=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        resource = data
        resource_list.append(resource)
    return resource_list

def importGeneList(filename,limit=None):
    ### Optionally limit the number of results imported
    gene_list=[]
    fn=filepath(filename); resource_list=[]; count=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        gene = string.split(data,'\t')[0]
        if ' ' in gene:
            gene = string.split(gene,' ')[0]
        if ':' in gene:
            gene = string.split(gene,':')[0]
        if gene not in gene_list and gene != 'GeneID' and gene != 'UID' and gene != 'probesetID':
            gene_list.append(gene)
            count+=1
            if limit != None:
                if limit==count: break
    gene_list = string.join(gene_list,',')
    return gene_list

def exportJunctionList(filename,limit=None):
    ### Optionally limit the number of results imported
    parent = export.findParentDir(filename)
    export_file = parent+'/top'+str(limit)+'/MultiPath-PSI.txt'#+file
    #export_file = filename[:-4]+'-top-'+str(limit)+'.txt'
    eo = export.ExportFile(export_file)
    fn=filepath(filename); count=0; firstLine=True
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstLine:
            firstLine = False
        elif '-' in t[0]:
            junctions = string.split(data,'\t')[0]
            junctions = string.replace(junctions,'|',' ')
            junctions = string.join(string.split(junctions,':')[1:],':')
            eo.write(junctions+'\n')
            count+=1
            if limit==count: break
        else:
            junctions = t[1] #Atg9a:ENSMUSG00000033124:E1.1-E3.1|ENSMUSG00000033124:E1.1-E3.2	
            junctions = string.split(junctions,'|') #ENSMUSG00000032314:I11.1_55475101-E13.1-ENSMUSG00000032314:E11.1-E13.1|ENSMUSG00000032314:I11.1_55475153;I11.1_55475101
            for junction_pair in junctions:
                if '-' in junction_pair:
                    try:
                        a,b = string.split(junction_pair,'-ENS')
                        b = 'ENS'+b
                        eo.write(a+' '+b+'\n')
                        count+=1
                        if limit==count: break      
                    except Exception:
                        pass
            if count>limit: break
        if count>limit: break
    eo.close()
    return export_file

def importConfigFile():
    #print "importing config file"
    filename = 'Config/config.txt'
    fn=filepath(filename); config_db={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        config_type,options = string.split(data,'\t')
        config_db[config_type] = options
    return config_db

def exportConfigFile(config_db):
    #print "exporting config file"
    new_file = 'Config/config.txt'
    data = export.ExportFile(new_file)
    for config in config_db:
        try: data.write(config+'\t'+str(config_db[config])+'\n'); data.close()
        except Exception:
            print 'Cannot write Config/config.txt to the Config directory (like Permissions Error)'

def remoteOnlineDatabaseVersions():
    db_versions = importOnlineDatabaseVersions()
    return db_versions_vendors,db_versions

def importOnlineDatabaseVersions():
    filename = 'Config/array_versions.txt'
    fn=filepath(filename); global db_versions; db_versions={}; global db_versions_vendors; db_versions_vendors={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        species,version,vendors = string.split(data,'\t')
        vendors  = string.split(vendors,'|')
        ad = ArrayData('','',vendors,'',species)
        version = string.replace(version,'Plus','') ### The user won't understand the Plus which relates to the GO-Elite version (AltAnalyze does not have a plus but we want the Plus for goelite)
        try: db_versions[version].append(species)
        except KeyError: db_versions[version] = [species]
        try: db_versions_vendors[version].append(ad)
        except KeyError: db_versions_vendors[version] = [ad]
    return db_versions

def getOnlineDBConfig(file_location_defaults,root):
    base_url = file_location_defaults['url'].Location()
    try:
        fln1,status1 = update.download(base_url+'Config/species_all.txt','Config/','')
        fln2,status2 = update.download(base_url+'Config/source_data.txt','Config/','')
        fln3,status3 = update.download(base_url+'Config/array_versions.txt','Config/','')
    except Exception:
        print 'Could not download the latest online configuration files (likely Permissions Error)'

    try:
        if 'Internet' not in status3:
            print 'Finished downloading the latest configuration files.'; root.destroy()
        else:
            try: WarningWindow(status3,'Error Encountered!'); root.destroy(); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
            except Exception: print status3; root.destroy(); sys.exit()
    except Exception: pass

def getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,update_goelite_resources,root):
    base_url = file_location_defaults['url'].Location()
    goelite_url = file_location_defaults['goelite'].Location()
    dbs_added = 0

    AltAnalyze_folders = read_directory(''); Cytoscape_found = 'no'
    for dir in AltAnalyze_folders:
        if 'Cytoscape_' in dir: Cytoscape_found='yes'
    if Cytoscape_found == 'no':
        fln,status = update.download(goelite_url+'Cytoscape/cytoscape.tar.gz','','')
        if 'Internet' not in status: print "Cytoscape program folder downloaded."

    count = verifyFileLength('AltDatabase/TreeView/TreeView.jar')
    if count==0:
        fln,status = update.download(goelite_url+'TreeView.zip','AltDatabase/NoVersion','')
        if 'Internet' not in status: print "TreeView program downloaded."
        
    fln,status = update.download(goelite_url+'Databases/'+db_version+'Plus/OBO.zip','AltDatabase/goelite/','')
    if 'Internet' not in status: print "Gene Ontology structure files downloaded."
    
    for species_code in new_species_codes:
        #print [base_url+'AltDatabase/'+db_version+'/'+species_code+'.zip']
        if species_code == 'Mm' or species_code == 'Hs' or species_code == 'Rn': specific_extension=''
        else: specific_extension='_RNASeq'
        fln,status = update.download(base_url+'AltDatabase/updated/'+db_version+'/'+species_code+specific_extension+'.zip','AltDatabaseNoVersion/','')
        if 'Internet' not in status:
            print 'Finished downloading the latest species database files.'
            dbs_added+=1
            #print goelite_url+'Databases/'+db_version+'Plus/'+species_code+'.zip'
        try: fln,status = update.download(goelite_url+'Databases/'+db_version+'Plus/'+species_code+'.zip','AltDatabase/goelite/','')
        except Exception: print "No species GO-Elite database found."
        
        if update_goelite_resources == 'yes': ### Get all additional GeneSet database types (can be lengthy download times)
            try: getAdditionalOnlineResources(species_code, 'All Resources',None)
            except Exception: print "Unable to update additional GO-Elite resources."
        
        if 'Internet' not in status: print "GO-Elite database installed." ; dbs_added+=1
        else: print "No species GO-Elite database found."
        try: os.mkdir(filepath('AltDatabase/'+species_code))
        except Exception: pass
    if dbs_added>0:
        print_out = "New species data successfully added to database."
        if root !='' and root !=None:
            try: InfoWindow(print_out,'Continue')
            except Exception: print print_out
        else: print print_out
        try: root.destroy()
        except Exception: pass
    else:
        if root !='' and root !=None: WarningWindow(status,'Error Encountered!'); root.destroy(); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        else: print status; root.destroy(); sys.exit()
    
def filterExternalDBs(all_external_ids,externalDBName_list,external_ids,array_db):
    filtered_external_list=[]
    for name in externalDBName_list:
        if name in external_ids:
            id = external_ids[name]
            if id in all_external_ids:
                if name != 'GO': filtered_external_list.append(name)
    for array in array_db:
        if '\\N_' not in array: filtered_external_list.append(array)
    return filtered_external_list

def updateOBOfiles(file_location_defaults,update_OBO,OBO_url,root):
    run_parameter = "Create/Modify Databases"
    if update_OBO == 'yes':
        from import_scripts import OBO_import
        c = OBO_import.GrabFiles()
        c.setdirectory('/OBO'); file_dirs = c.searchdirectory('.ontology')+c.searchdirectory('.obo')
        if len(OBO_url)>0: obo = OBO_url
        else: ### If not present, get the gene-ontology default OBO file
            obo = file_location_defaults['OBO'].Location()
        fln,status = update.download(obo,'OBO/','')
        run_parameter='Create/Modify Databases'
        if 'Internet' not in status:
            OBO_import.moveOntologyToArchiveDir()

            print_out = 'Finished downloading the latest Ontology OBO files.'
            print print_out

            try: system_codes,source_types,mod_types = GO_Elite.getSourceData()
            except Exception: null=[]
            
            if root !='' and root !=None:
                InfoWindow(print_out,'Update Complete!')
                continue_to_next_win = Button(text = 'Continue', command = root.destroy)
                continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()            
                GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            else: null=[]
        else:
            if root !='' and root !=None: WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            else: print status
    else:
        print_out = 'Download Aborted.'
        if root !='' and root !=None: WarningWindow(print_out,print_out); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        else: print print_out
        
def importExternalDBs(species_full):
    filename = 'Config/EnsExternalDBs.txt'
    fn=filepath(filename); x = 0; external_dbs=[]; external_system={}; all_databases={}; external_ids={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            id, database, species_specific, exclude, system_code = string.split(data,'\t')
            external_ids[database] = int(id)
            if database != 'GO':
                all_databases[database]=system_code
                if (species_full == species_specific) or len(species_specific)<2:
                    if len(exclude)<2:
                        external_system[database] = system_code

    filename = 'Config/external_db.txt'; external_system2={}
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            t = string.split(data,'\t'); id = int(t[0]); database = t[1]
            external_ids[database] = id
            if database in external_system:
                external_system2[database] = external_system[database]
            elif database not in all_databases: ### Add it if it's new
                try:
                    try: system = database[:3]
                    except Exception: system = database[:2]
                    external_system2[database] = system
                except Exception: null=[]
        except Exception: null=[] ### Occurs when a bad end of line is present

    filename = 'Config/array.txt'
    global array_db; array_db={}
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try:
            array = t[1]; vendor = t[3]
            database = vendor+'_'+array; array_db[database]=[]
            if database in external_system:
                external_system2[database] = external_system[database]
            if database in all_databases:
                external_system2[database] = all_databases[database]
            elif database not in all_databases: ### Add it if it's new
                try:
                    if vendor == 'AFFY': system = 'X'
                    if vendor == 'ILLUMINA': system = 'Il'
                    if vendor == 'CODELINK': system = 'Co'
                    if vendor == 'AGILENT': system = 'Ag'
                    else: system = 'Ma'  ###Miscelaneous Array type
                    external_system2[database] = system
                except Exception: null=[]
        except Exception: null=[]
    external_system = external_system2
    #try: del external_system['GO']
    #except Exception: null=[]
    for database in external_system: external_dbs.append(database)
    external_dbs.append(' '); external_dbs = unique.unique(external_dbs); external_dbs.sort()
    return external_dbs,external_system,array_db,external_ids
    
class SupprotedArrays:
    def __init__(self, array_name, library_file, annotation_file, species, array_type):
        self.array_name = array_name; self.library_file = library_file; self.annotation_file = annotation_file
        self.species = species; self.array_type = array_type
    def ArrayName(self): return self.array_name
    def LibraryFile(self): return self.library_file
    def AnnotationFile(self): return self.annotation_file
    def Species(self): return self.species
    def ArrayType(self): return self.array_type
    def __repr__(self): return self.ArrayName()
    
def importSupportedArrayInfo():
    filename = 'Config/ArrayFileInfo.txt'; x=0
    fn=filepath(filename); global supproted_array_db; supproted_array_db={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        array_name,library_file,annotation_file,species,array_type = string.split(data,'\t')
        if x==0: x=1
        else:
            sd = SupprotedArrays(array_name,library_file,annotation_file,species,array_type)
            supproted_array_db[array_name] = sd
    return supproted_array_db

def exportSupportedArrayInfo():
    fn=filepath('Config/ArrayFileInfo.txt'); data = open(fn,'w'); x=0
    header = string.join(['ArrayName','LibraryFile','AnnotationFile','Species','ArrayType'],'\t')+'\n'
    data.write(header)
    for array_name in supproted_array_db:
        sd = supproted_array_db[array_name]
        values = [array_name,sd.LibraryFile(),sd.AnnotationFile(),sd.Species(),sd.ArrayType()]
        values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()

class SystemData:
    def __init__(self, syscode, sysname, mod):
        self._syscode = syscode; self._sysname = sysname; self._mod = mod
    def SystemCode(self): return self._syscode
    def SystemName(self): return self._sysname
    def MOD(self): return self._mod
    def __repr__(self): return self.SystemCode()+'|'+self.SystemName()+'|'+self.MOD()

def remoteSystemInfo():
    system_codes,system_list,mod_list = importSystemInfo(returnSystemCode=True)
    return system_codes,system_list,mod_list

def getSystemInfo():
    importSystemInfo()
    return system_codes

def importSystemInfo(returnSystemCode=False):
    filename = 'Config/source_data.txt'; x=0
    fn=filepath(filename); global system_list; system_list=[]; global system_codes; system_codes={}; mod_list=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if '!DOCTYPE' in data:
            fn2 = string.replace(fn,'.txt','_archive.txt')
            import shutil; shutil.copyfile(fn2,fn) ### Bad file was downloaded (with warning)
            importSystemInfo(); break

	elif '<html>' in data:
	    print_out = "WARNING!!! Connection Error. Proxy may not be allowed from this location."
	    try: WarningWindow(print_out,' Continue ')
	    except NameError: print print_out
	    importSystemInfo(); break
	else:
	    try: sysname=t[0];syscode=t[1]
	    except Exception: sysname=''            
	try: mod = t[2]
	except Exception: mod = ''
	if x==0: x=1
	else:
	    system_list.append(sysname)
	    ad = SystemData(syscode,sysname,mod)
	    if len(mod)>1: mod_list.append(sysname)
	    system_codes[sysname] = ad

    if returnSystemCode:
        return system_codes,system_list,mod_list
    else:
        return system_list,mod_list

def exportSystemInfoRemote(system_code_db):
    global system_codes; system_codes = system_code_db
    exportSystemInfo()
    
def exportSystemInfo():
    if len(system_codes)>0:
        filename = 'Config/source_data.txt'
        fn=filepath(filename); data = open(fn,'w')
        header = string.join(['System','SystemCode','MOD_status'],'\t')+'\n'
        data.write(header)
        for sysname in system_codes:
            ad = system_codes[sysname]
            values = string.join([sysname,ad.SystemCode(),ad.MOD()],'\t')+'\n'
            data.write(values)
        data.close()

class SpeciesData:
    def __init__(self, abrev, species, algorithms):
        self._abrev = abrev; self._species = species; self._algorithms = algorithms
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def Algorithms(self): return self._algorithms
    def __repr__(self): return self.Report()

def getSpeciesInfo():
    ### Used by AltAnalyze
    global integrate_online_species; integrate_online_species = 'yes'
    importSpeciesInfo(); species_names={}
    for species_full in species_codes:
        sc = species_codes[species_full]; abrev = sc.SpeciesCode()
        species_names[abrev] = species_full
    return species_names

def remoteSpeciesInfo():
    global integrate_online_species; integrate_online_species = 'yes'
    importSpeciesInfo()
    return species_codes
    
def remoteSpeciesAlt():
    ### Replicates the output of GO-Elite's species importer
    global integrate_online_species; integrate_online_species = 'yes'
    importSpeciesInfo()
    species_names={}
    for species in species_codes:
        sd = species_codes[species]
        species_names[sd.SpeciesCode()] = sd
    return species_names

def importSpeciesInfo():
    try:
        if integrate_online_species == 'yes': filename = 'Config/species_all.txt'
        else: filename = 'Config/species.txt'
    except Exception: filename = 'Config/species.txt'
    
    fn=filepath(filename)
    global species_list
    species_list=[]
    global species_codes
    species_codes={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            try: abrev,species,algorithms = string.split(data,'\t')
            except Exception:
                try: abrev,species = string.split(data,'\t'); algorithms = ''
                except Exception:
                    abrev,species,taxid,compatible_mods = string.split(data,'\t')
                    algorithms = ''
        except Exception:
            if '!DOCTYPE': print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
            else: print_out = "Unknown file error encountered."
            IndicatorWindow(print_out,'Continue')
            raw = export.ExportFile(fn); raw.close(); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        if x==0: x=1
        else:
            algorithms = string.split(algorithms,'|')
            species_list.append(species)
            sd = SpeciesData(abrev,species,algorithms)
            species_codes[species] = sd
    return species_codes

def exportSpeciesInfo(species_codes):
    fn=filepath('Config/species.txt'); data = open(fn,'w'); x=0
    header = string.join(['species_code','species_name','compatible_algorithms'],'\t')+'\n'
    data.write(header)
    for species in species_codes:
        sd = species_codes[species]; algorithms = string.join(sd.Algorithms(),'|')
        values = [sd.SpeciesCode(),species,algorithms]
        values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()
            
class ArrayGroupData:
    def __init__(self, array_header, group, group_name):
        self._array_header = array_header; self._group = group; self._group_name = group_name
    def Array(self): return self._array_header
    def Group(self): return self._group
    def setGroup(self,group): self._group = group
    def GroupName(self): return self._group_name
    def setGroupName(self,group_name): self._group_name = group_name
    def Report(self): return self.Array()
    def __repr__(self): return self.Report()

def importArrayGroupsSimple(expr_group_dir,cel_files):
    array_group_list = []; group_db={}
    fn=filepath(expr_group_dir)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        array_header,group,group_name = string.split(data,'\t')
        if group_name == 'NA': group_name = 'None'
        #print [array_header],cel_files
        if (array_header in cel_files) or len(cel_files)==0: ### restrict import to array files listed in the groups file
            try: group = int(group); group_db[group]=group_name
            except ValueError: print group, group_name;kill
            agd = ArrayGroupData(array_header,group,group_name)
            array_group_list.append(agd)
    if len(cel_files)>0:
        if len(cel_files)!=len(array_group_list):
            #print len(cel_files),len(array_group_list)
            #print cel_files
            array_group_list2=[]
            for i in array_group_list:
                if i.Array() not in cel_files:
                    print [i.Array()], 'not in CEL file dir (in groups file)'
                array_group_list2.append(i.Array())
            for i in cel_files:
                if i not in array_group_list2:
                    print [i], 'not in groups file (in CEL file dir)'
            raise NameError('Samples In Groups Not Found In Dir')
    return array_group_list,group_db

class ArrayData:
    def __init__(self, abrev, array, manufacturer, constitutive_source, species):
        self._abrev = abrev; self._array = array; self._manufacturer = manufacturer; self._species = species
        self._constitutive_source = constitutive_source
    def ArrayCode(self): return self._abrev
    def ArrayName(self): return self._array
    def Manufacturer(self): return self._manufacturer
    def ConstitutiveSource(self): return self._constitutive_source
    def SpeciesCodes(self): return self._species
    def setSpeciesCodes(self,species): species = self._species
    def __repr__(self): return self.ArrayCode()+'|'+str(self.SpeciesCodes())+'|'+str(self.Manufacturer())

def remoteArrayInfo():
    importArrayInfo()
    return array_codes

def importArrayInfo():
    filename = 'Config/arrays.txt'; x=0
    fn=filepath(filename); global array_list; array_list=[]; global array_codes; array_codes={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        abrev,array,manufacturer,constitutive_source,species = string.split(data,'\t')
        if x==0: x=1
        else:
            species = string.split(species,'|')
            array_list.append(array)
            ad = ArrayData(abrev,array,manufacturer,constitutive_source,species)
            array_codes[array] = ad
    return array_list

def exportArrayInfo(array_codes):
    fn=filepath('Config/arrays.txt'); data = open(fn,'w'); x=0
    header = string.join(['array_type','array_name','manufacturer','constitutive_source','compatible_species'],'\t')+'\n'
    data.write(header)
    for array in array_codes:
        ad = array_codes[array]; species = string.join(ad.SpeciesCodes(),'|')
        values = [ad.ArrayCode(),array,ad.Manufacturer(),ad.ConstitutiveSource(),species]
        values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()
            
class FileLocationData:
    def __init__(self, status, location, species):
        self._status = status; self._location = location; self._species = species
    def Status(self): return self._status
    def Location(self): return self._location
    def SetLocation(self,location): self._location = location
    def Species(self): return self._species
    def __repr__(self): return self.Report()
    
def importDefaultFileLocations():
    filename = 'Config/default-files.csv'; x=0
    fn=filepath(filename); file_location_defaults={}
    for line in open(fn,'rU').readlines():
        line = string.replace(line,',','\t') ### Make tab-delimited (had to make CSV since Excel would impoperly parse otherwise)
        data = cleanUpLine(line)
        ###Species can be multiple species - still keep in one field
        try: app,status,location,species = string.split(data,'\t')
        except Exception:
            try:
                t = string.split(data,'\t')
                app=t[0]; status=t[1]; location=t[2]; species=t[3]
            except Exception:
                continue
        fl = FileLocationData(status, location, species)
        if species == 'all': file_location_defaults[app] = fl
        else:
            try: file_location_defaults[app].append(fl)
            except KeyError: file_location_defaults[app] = [fl]
    return file_location_defaults
    
def exportDefaultFileLocations(file_location_defaults):
    ### If the user supplies new defaults, over-write the existing
    fn=filepath('Config/default-files.csv'); data = open(fn,'w')
    for app in file_location_defaults:
        fl_list = file_location_defaults[app]
        try:
            for fl in fl_list:
                values = [app,fl.Status(),fl.Location(),fl.Species()]
                values = '"'+string.join(values,'","')+'"'+'\n'
                data.write(values)
        except Exception:
            fl = fl_list
            values = [app,fl.Status(),fl.Location(),fl.Species()]
            values = '"'+string.join(values,'","')+'"'+'\n'
            data.write(values)
    data.close()

def exportGroups(exp_file_location_db,array_group_list,filetype='Groups'):
    ### If the user supplies new defaults, over-write the existing
    for dataset_name in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]
        groups_file = fl.GroupsFile()
        if filetype =='Batch':
            groups_file = string.replace(groups_file,'groups.','batch.')
        fn=filepath(groups_file); data = open(fn,'w')
        value_list = [] ### Sort grouped results based on group number
        for agd in array_group_list:
            values = [agd.Array(), str(agd.Group()), agd.GroupName()]
            values = string.join(values,'\t')+'\n'; value_list.append(((agd.Group(),agd.Array()),values))
        value_list.sort()
        for values in value_list: data.write(values[-1])
        data.close() 

def exportComps(exp_file_location_db,comp_group_list):
    ### If the user supplies new defaults, over-write the existing
    for dataset_name in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]; comps_file = fl.CompsFile()
        fn=filepath(comps_file); data = open(fn,'w')
        for comp_num, groups in comp_group_list:
            group1, group2 = groups
            values = [str(group1), str(group2)]
            values = string.join(values,'\t')+'\n'; data.write(values)
        data.close()
        
class Defaults:
    def __init__(self, abrev, array, species):
        self._abrev = abrev; self._array = array; self._species = species
    def ArrayCode(self): return self._abrev
    def ArrayName(self): return self._array
    def Species(self): return self._species
    def __repr__(self): return self.Report()

def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found

def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>9: break
    except Exception: pass
    return count

def getGeneSystem(filename):
    firstRow=True
    count=0
    system = 'Symbol'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            if firstRow: firstRow=False
            else:
                id = string.split(line,'\t')[0]
                if 'ENS' in id: system = 'Ensembl'
                count+=1
            if count>9: break
    except Exception: pass
    return system

def determinePlatform(filename):
    platform = ''
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            #print [line]
            if len(line)>0: platform = line
    except Exception: pass
    return platform
         
def importDefaults(array_type,species):
    filename = 'Config/defaults-expr.txt'
    expr_defaults = importDefaultInfo(filename,array_type)
    #perform_alt_analysis, expression_data_format, dabg_p, expression_threshold, avg_all_for_ss, include_raw_data
    
    filename = 'Config/defaults-alt_exon.txt'
    alt_exon_defaults = importDefaultInfo(filename,array_type)
    #analysis_method, alt_exon_fold_variable, p_threshold, filter_probeset_types, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold,run_MiDAS, export_splice_index_values = values
    
    filename = 'Config/defaults-funct.txt'
    functional_analysis_defaults = importDefaultInfo(filename,array_type)
    #analyze_functional_attributes,microRNA_prediction_method = functional_analysis_defaults

    filename = 'Config/defaults-goelite.txt'
    goelite_defaults = importDefaultInfo(filename,array_type)
    return expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults

def importDefaultInfo(filename,array_type):
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if '-expr' in filename:
            array_abrev, dabg_p, rpkm_threshold, gene_exp_threshold, exon_exp_threshold, exon_rpkm_threshold, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, normalize_feature_exp, normalize_gene_data, avg_all_for_ss, include_raw_data, probability_algorithm, FDR_statistic, batch_effects, marker_finder, visualize_results, run_lineage_profiler, run_goelite = string.split(data,'\t')
            if array_type == array_abrev:
                return dabg_p, rpkm_threshold, gene_exp_threshold, exon_exp_threshold, exon_rpkm_threshold, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, normalize_feature_exp, normalize_gene_data, avg_all_for_ss, include_raw_data, probability_algorithm, FDR_statistic, batch_effects, marker_finder, visualize_results, run_lineage_profiler, run_goelite
            
        if '-alt' in filename:
            array_abrev, analysis_method, additional_algorithms, filter_probeset_types, analyze_all_conditions, p_threshold, alt_exon_fold_variable, additional_score, permute_p_threshold, gene_expression_cutoff, remove_intronic_junctions, perform_permutation_analysis, export_splice_index_values, run_MiDAS, calculate_splicing_index_p, filter_for_AS = string.split(data,'\t')
            if array_type == array_abrev:
                return  [analysis_method, additional_algorithms, filter_probeset_types, analyze_all_conditions, p_threshold, alt_exon_fold_variable, additional_score, permute_p_threshold, gene_expression_cutoff, remove_intronic_junctions, perform_permutation_analysis, export_splice_index_values, run_MiDAS, calculate_splicing_index_p, filter_for_AS]
            
        if '-funct' in filename:
            array_abrev, analyze_functional_attributes, microRNA_prediction_method = string.split(data,'\t')
            if array_type == array_abrev:
                return [analyze_functional_attributes,microRNA_prediction_method]

        if '-goelite' in filename:
            array_abrev, ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, ORA_algorithm, resources_to_analyze, pathway_permutations, mod, returnPathways, get_additional = string.split(data,'\t')
            if array_type == array_abrev:
                return [ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, ORA_algorithm, resources_to_analyze, pathway_permutations, mod, returnPathways, get_additional]
        
class OptionData:
    def __init__(self,option,displayed_title,display_object,notes,array_options,global_default):
        self._option = option; self._displayed_title = displayed_title; self._notes = notes
        self._array_options = array_options; self._display_object = display_object
        if len(global_default)>0:
            if '|' in global_default:
                global_default = string.split(global_default,'|') ### store as a list
            self._default_option = global_default
    def Option(self): return self._option
    def VariableName(self): return self._option
    def Display(self): return self._displayed_title
    def setDisplay(self,display_title): self._displayed_title = display_title
    def setDisplayObject(self,display_object): self._display_object = display_object
    def DisplayObject(self): return self._display_object
    def Notes(self): return self._notes
    def setNotes(self,notes): self._notes = notes
    def DefaultOption(self): return self._default_option
    def setDefaultOption(self,default_option): self._default_option = default_option
    def setNotes(self,notes): self._notes = notes
    def ArrayOptions(self): return self._array_options
    def setArrayOptions(self,array_options): self._array_options = array_options
    def Options(self): return self._array_options
    def __repr__(self): return self.Option()+'|'+self.Display()

def importUserOptions(array_type,vendor=None):
    filename = 'Config/options.txt'; option_db={}; option_list_db={}
    fn=filepath(filename); x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        data = string.replace(data,'\k','\n') ###Used \k in the file instead of \n, since these are removed above
        if array_type == 'RNASeq':
            data = string.replace(data,'probeset','junction')
            data = string.replace(data,'probe set','junction, exon or gene')
            data = string.replace(data,'CEL file','BED, BAM, TAB or TCGA junction file')
        if vendor != 'Affymetrix':
            data = string.replace(data,'probe set','gene')      
        if vendor == 'Agilent':
            if 'CEL file' in data:
                data = string.replace(data,'CEL file','Feature Extraction file')
                data = string.replace(data,' (required)','')
        if array_type == '10XGenomics':
            data = string.replace(data,'CEL file containing folder','Chromium filtered matrix.mtx file')
        try:
            if '10X' in vendor:
                data = string.replace(data,'CEL file containing folder','Chromium filtered matrix.mtx file')
        except Exception: pass
        
        t = string.split(data,'\t')
        #option,mac_displayed_title,pc_displayed_title,pc_display2,linux_displayed_title,display_object,group,notes,description,global_default = t[:10]
        option,displayed_title,display_object,group,notes,description,global_default = t[:7]
        """
        if os.name == 'nt':
            import platform
            if '64' in platform.machine(): displayed_title = pc_display2
            elif '32' in platform.machine(): displayed_title = pc_display2
            elif '64bit' in platform.architecture(): displayed_title = pc_display2
            else: displayed_title = pc_displayed_title
        elif 'darwin' in sys.platform: displayed_title = mac_displayed_title
        elif 'linux' in sys.platform: displayed_title = linux_displayed_title
        else: displayed_title = linux_displayed_title
        """
        """
        try:
            if option == 'rho_cutoff' and '10X' in vendor:
                global_default = '0.3'
            if option == 'restrictBy' and '10X' in vendor:
                global_default = 'yes'
            if option == 'column_metric_predict' and '10X' in vendor:
                global_default = 'euclidean'    
        except Exception:
            pass
        """
 
        if 'junction' in displayed_title: displayed_title+=' '
        """if array_type == 'RNASeq':
            if option == 'dabg_p': ### substitute the text for the alternatitve text in notes
                displayed_title = notes"""
        if x == 0:
            i = t.index(array_type) ### Index position of the name of the array_type selected by user (or arbitrary to begin with)
            x = 1
        else:
            array_options = t[i]
            if array_type == "3'array":
                """
                if 'normalize_gene_data' in data and vendor != 'Agilent':
                    array_options = 'NA' ### only applies currently to Agilent arrays """
                if 'channel_to_extract' in data and vendor != 'Agilent':
                    array_options = 'NA' ### only applies currently to Agilent arrays
            array_options = string.split(array_options,'|')
            od = OptionData(option,displayed_title,display_object,notes,array_options,global_default)
            option_db[option] = od
            try: option_list_db[group].append(option) ###group is the name of the GUI menu group
            except KeyError: option_list_db[group] = [option]
    return option_list_db,option_db

class SummaryResults:
    def __init__(self):
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        LINKS=('http://www.altanalyze.org','')
        self.LINKS = LINKS
        try: tl = Toplevel()
        except Exception: tl = Tkinter.Toplevel()
        tl.title('AltAnalyze')
        
        filename = 'Config/icon.gif'
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(tl); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW); use_scroll = 'no'

        label_text_str = 'AltAnalyze Result Summary'; height = 250; width = 700      
        self.sf = PmwFreeze.ScrolledFrame(tl,
            labelpos = 'n', label_text = label_text_str,
            usehullsize = 1, hull_width = width, hull_height = height)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        tl.mainloop()
        
class FeedbackWindow:
    def __init__(self,message,button_text,button_text2):
        self.message = message; self.button_text = button_text; self.button_text2 = button_text2
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')
        self._user_variables={}
        
        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
            
        try:
            can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())
            can.create_image(2, 2, image=img, anchor=NW)
        except Exception: pass
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()
        text_button = Button(parent, text=self.button_text, command=self.button1); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button2 = Button(parent, text=self.button_text2, command=self.button2); text_button2.pack(side = 'bottom', padx = 5, pady = 5)
        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def button1(self): self._user_variables['button']=self.button_text; self._parent.destroy()
    def button2(self): self._user_variables['button']=self.button_text2; self._parent.destroy()
    def ButtonSelection(self): return self._user_variables
    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()

class IndicatorWindowSimple:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        try: 
            can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
            can.create_image(2, 2, image=img, anchor=NW)
        except Exception: pass
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
        
class IndicatorWindow:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        try:
            can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
            can.create_image(2, 2, image=img, anchor=NW)
        except Exception: pass
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        quit_button = Button(parent, text='Quit', command=self.quit); quit_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
    def quit(self):
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: self._parent.quit(); sys.exit()

class DownloadWindow:
    def __init__(self,message,option1,option2):
        self._user_variables = user_variables
        if len(option2)==2: option2,option3 = option2; num_options = 3; self.option3 = option3
        else: num_options = 2
        self.message = message; self.option1 = option1; self.option2 = option2
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        text_button = Button(parent, text=self.option1, command=self.selected1); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button2 = Button(parent, text=self.option2, command=self.selected2); text_button2.pack(side = 'bottom', padx = 5, pady = 5)
        if num_options == 3:
            text_button3 = Button(parent, text=self.option3, command=self.selected3); text_button3.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
    def selected1(self):
        self._user_variables['selected_option']=1; self._parent.destroy()
    def selected2(self):
        self._user_variables['selected_option']=2; self._parent.destroy()
    def selected3(self):
        self._user_variables['selected_option']=3; self._parent.destroy()
    def Results(self): return self._user_variables

class IndicatorLinkOutWindow:
    def __init__(self,message,button_text,url):
        self.message = message; self.button_text = button_text; nulls = '\t\t\t\t\t\t\t';
        parent = Tk(); self._parent = parent; parent.title('Attention!!!'); self.url = url

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        continue_button = Button(parent, text='Continue', command=parent.destroy); continue_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button = Button(parent, text=self.button_text, command=self.linkout); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
    def linkout(self):
        webbrowser.open(self.url)
            
class IndicatorChooseWindow:
    def __init__(self,message,button_text):
        self.message = message; self               .button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        #text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        option=''
        def foldercallback(callback=self.callback,option=option): self.chooseDirectory(option)
        choose_win = Button(self._parent, text=self.button_text,command=foldercallback); choose_win.pack(padx =  3, pady = 3)
        quit_button = Button(parent, text='Quit', command=self.quit); quit_button.pack(padx = 3, pady = 3)
        parent.mainloop()
    def quit(self):
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: self._parent.quit(); sys.exit()
    def callback(self, tag, option): null = ''
    def chooseDirectory(self,option):
        tag = tkFileDialog.askdirectory(parent=self._parent)
        ### Below is code specific for grabbing the APT location
        from import_scripts import ResultsExport_module
        apt_location = ResultsExport_module.getAPTDir(tag)
        if 'bin' not in apt_location:
            print_out = "WARNING!!! Unable to find a valid Affymetrix Power Tools directory."
            try: WarningWindow(print_out,' Continue ')
            except NameError: print print_out
            self._tag = ''
        else: self._tag = apt_location
        self.destroy_win()
    def destroy_win(self):
        try: self._parent.quit(); self._parent.destroy()
        except Exception: self._parent.quit(); sys.exit()
    def Folder(self): return self._tag
    
class WarningWindow:
    def __init__(self,warning,window_name):
        try: tkMessageBox.showerror(window_name, warning)
        except Exception:
            print warning
            #print window_name; sys.exit()
            kill

class InfoWindow:
    def __init__(self,dialogue,header):
        try: tkMessageBox.showinfo(header, dialogue)
        except Exception:
            print dialogue
            #print 'Attempted to open a GUI that is not accessible...exiting program';sys.exit()
            #print "Analysis finished...exiting AltAnalyze."; sys.exit()
       
class MacConsiderations:
    def __init__(self):
        parent = Tk()
        self._parent = parent
        parent.title('AltAnalyze: Considerations for Mac OSX')
        self._user_variables={}
        filename = 'Config/MacOSX.png'
        try:
            import ImageTk
            img = ImageTk.PhotoImage(file=filepath(filename))
        except Exception:
            try:
                from PIL import ImageTk
                img = ImageTk.PhotoImage(file=filepath(filename))
            except Exception: 
                filename = 'Config/MacOSX.gif'
                fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='top',fill=BOTH); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)

        ### Add some buttons to the horizontal RadioSelect
        continue_to_next_win = Tkinter.Button(text = 'Continue', command = parent.destroy)
        continue_to_next_win.pack(side = 'right', padx = 5, pady = 5);

        info_win = Button(self._parent, text="Online Help", command=self.Linkout)
        info_win.pack(side = 'left', padx = 5, pady = 5)
        self.url = 'http://www.altanalyze.org/MacOSX_help.html'

        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def Linkout(self):
        try: webbrowser.open(self.url)
        except Exception,e: print e
        
    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()

    def callback(self, tag):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables['continue'] = tag   
        
class MainMenu:
    def __init__(self):
        parent = Tk()
        self._parent = parent
        parent.title('AltAnalyze: Introduction')
        self._user_variables={}
        filename = 'Config/logo.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='top',fill=BOTH); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)

        """      
        ### Create and pack a horizontal RadioSelect widget.
        def buttoncallback(tag,callback=self.callback):
            callback(tag)
        horiz = PmwFreeze.RadioSelect(parent,
                labelpos = 'w', command = buttoncallback,
                label_text = 'AltAnalyze version 1.155 Main', frame_borderwidth = 2,
                frame_relief = 'ridge'
        ); horiz.pack(fill = 'x', padx = 10, pady = 10)
        for text in ['Continue']: horiz.add(text)
        """
        ### Add some buttons to the horizontal RadioSelect
        continue_to_next_win = Tkinter.Button(text = 'Begin Analysis', command = parent.destroy)
        continue_to_next_win.pack(side = 'bottom', padx = 5, pady = 5);

        info_win = Button(self._parent, text="About AltAnalyze", command=self.info)
        info_win.pack(side = 'bottom', padx = 5, pady = 5)

        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def info(self):
        
        """
        ###Display the information using a messagebox
        about = 'AltAnalyze version 2.1.2.\n'
        about+= 'AltAnalyze is an open-source, freely available application covered under the\n'
        about+= 'Apache open-source license. Additional information can be found at:\n'
        about+= 'http://www.altanalyze.org\n'
        about+= '\nDeveloped by:\nNathan Salomonis Lab\nCincinnati Childrens Hospital Medical Center 2008-2016'
        tkMessageBox.showinfo("About AltAnalyze",about,parent=self._parent)
        """
        
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        LINKS=('http://www.altanalyze.org','')
        self.LINKS = LINKS
        tl = Toplevel() ### Create a top-level window separate than the parent        
        txt=Text(tl)

        #filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        #can = Canvas(tl); can.pack(side='left'); can.config(width=img.width(), height=img.height())        
        #can.create_image(2, 2, image=img, anchor=NW)
        
        txt.pack(expand=True, fill="both")
        txt.insert(END, 'AltAnalyze version 2.1.2.\n')
        txt.insert(END, 'AltAnalyze is an open-source, freely available application covered under the\n')
        txt.insert(END, 'Apache open-source license. Additional information can be found at:\n')
        txt.insert(END, "http://www.altanalyze.org\n", ('link', str(0)))
        txt.insert(END, '\nDeveloped by:\nDr. Nathan Salomonis Research Group\nCincinnati Childrens Hospital Medical Center 2008-2016')
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)
        
    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()

    def callback(self, tag):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables['continue'] = tag        

class LinkOutWindow:
    def __init__(self,output):
        ### Text window with link included
        url,text_list = output
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        LINKS=(url,'')
        self.LINKS = LINKS
        tl = Toplevel() ### Create a top-level window separate than the parent        
        txt=Text(tl)
        
        txt.pack(expand=True, fill="both")
        for str_item in text_list:
            txt.insert(END, str_item+'\n')
            txt.insert(END, "http://www.altanalyze.org\n", ('link', str(0)))
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()

def exportCELFileList(cel_files,cel_file_dir):
    fn=cel_file_dir+'/cel_files.txt'; data = open(fn,'w')
    data.write('cel_files'+'\n') ###header
    for cel_file in cel_files:
        data.write(cel_file+'\n')
    data.close()
    return fn

def predictGroupsAndComps(cel_files,output_dir,exp_name):
    fn1=output_dir+'/ExpressionInput/groups.'+exp_name+'.txt'; gdata = export.ExportFile(fn1)
    fn2=output_dir+'/ExpressionInput/comps.'+exp_name+'.txt'; cdata = export.ExportFile(fn2)
    fn3=output_dir+'/ExpressionInput/exp.'+exp_name+'.txt'
    delimited_db={}; delim_type={}; files_exported = 'no'

    for cel_file in cel_files:
        cel_name = cel_file
        cel_file = string.replace(cel_file,'.CEL','')
        cel_file = string.replace(cel_file,'.cel','')
        dashed_delim = string.split(cel_file,'-')
        dot_delim = string.split(cel_file,'.')
        under_delim = string.split(cel_file,'_')
        if len(dashed_delim) == 2:
            delim_type[1]=None
            try: delimited_db[dashed_delim[0]].append(cel_name)
            except KeyError: delimited_db[dashed_delim[0]] = [cel_name]
        elif len(under_delim) == 2:
            delim_type[2]=None
            try: delimited_db[under_delim[0]].append(cel_name)
            except KeyError: delimited_db[under_delim[0]] = [cel_name]
        elif len(dot_delim) == 2:
            delim_type[3]=None
            try: delimited_db[dot_delim[0]].append(cel_name)
            except KeyError: delimited_db[dot_delim[0]] = [cel_name]
    if len(delim_type)==1 and len(delimited_db)>1: ###only 1 type of delimiter used and at least 2 groups present
        group_index=0; group_db={}; files_exported = 'yes'
        for group in delimited_db:
            group_index+=1; group_db[str(group_index)]=None
            for array in delimited_db[group]:
                gdata.write(string.join([array,str(group_index),group],'\t')+'\n')
        for index1 in group_db: ### Create a comps file for all possible comps
            for index2 in group_db:
                if index1 != index2:
                    cdata.write(string.join([index1,index2],'\t')+'\n')
    gdata.close(); cdata.close()
    if files_exported == 'no':
        os.remove(fn1); os.remove(fn2)
        try: ExpressionBuilder.checkArrayHeaders(fn3,fn1) ### Create just the groups template file
        except Exception: pass ### This error will more likely occur since no expression file has been created 
    return files_exported

def formatArrayGroupsForGUI(array_group_list, category = 'GroupArrays'):
        ### Format input for GUI like the imported options.txt Config file, except allow for custom fields in the GUI class
        option_db={}; option_list={}
        
        if category != 'BatchArrays':
            ### Add a checkbox at the top to allow for automatic assignment of groups (e.g., Single Cell Data)
            option='PredictGroups';displayed_title='Run de novo cluster prediction (ICGS) to discover groups, instead';display_object='single-checkbox';notes='';array_options=['---']
            od = OptionData(option,displayed_title,display_object,notes,array_options,'')
            option_db[option] = od
            option_list[category] = [option]

        for agd in array_group_list:
            option = agd.Array(); array_options = [agd.GroupName()]; displayed_title=option; display_object='simple_entry'; notes=''
            od = OptionData(option,displayed_title,display_object,notes,array_options,'')
            option_db[option] = od
            try: option_list[category].append(option) ###group is the name of the GUI menu group
            except KeyError: option_list[category] = [option]
        return option_db,option_list
    
def importExpressionFiles():
    exp_file_location_db={}; exp_files=[]; parent_dir = 'ExpressionInput'+'/'+array_type
    fn =filepath(parent_dir+'/'); dir_files = read_directory('/'+parent_dir)

    stats_file_dir='' 
    for file in dir_files:
        if 'exp.' in file: exp_files.append(file)
    for file in exp_files:
        stats_file = string.replace(file,'exp.','stats.')
        groups_file = string.replace(file,'exp.','groups.')
        comps_file = string.replace(file,'exp.','comps.')
        if stats_file in dir_files: stats_file_dir = fn+stats_file
        if groups_file in dir_files and comps_file in dir_files:
            groups_file_dir = fn+groups_file; comps_file_dir = fn+comps_file
            exp_file_dir = fn+file
            fl = ExpressionFileLocationData(exp_file_dir,stats_file_dir,groups_file_dir,comps_file_dir)
            exp_file_location_db[file] = fl
    return exp_file_location_db

class ExpressionFileLocationData:
    def __init__(self, exp_file, stats_file, groups_file, comps_file):
        self._exp_file = exp_file; self._stats_file = stats_file; self._groups_file = groups_file
        self._comps_file = comps_file; self.biotypes='NA'
        import platform; self.architecture = platform.architecture()[0]
        self.normalize_feature_exp = 'NA'
        self.normalize_gene_data = 'NA'
        self.runKallisto =''
    def setExpFile(self, exp_file):self._exp_file=exp_file
    def ExpFile(self): return self._exp_file
    def StatsFile(self): return self._stats_file
    def CountsFile(self):
        import AltAnalyze
        counts_file = string.replace(self.ExpFile(),'exp.','counts.')
        file_length = AltAnalyze.verifyFileLength(counts_file)
        if file_length>0:
            return counts_file
        else:
            return self.ExpFile()
    def GroupsFile(self): return self._groups_file
    def setGroupsFile(self, groups_file):self._groups_file=groups_file
    def CompsFile(self): return self._comps_file
    def setCompsFile(self, comps_file):self._comps_file=comps_file
    def setArchitecture(self,architecture): self.architecture = architecture
    def setAPTLocation(self,apt_location): self._apt_location = osfilepath(apt_location)
    def setInputCDFFile(self,cdf_file): self._cdf_file = osfilepath(cdf_file)
    def setCLFFile(self,clf_file): self._clf_file = osfilepath(clf_file)
    def setBGPFile(self,bgp_file): self._bgp_file = osfilepath(bgp_file)
    def setCELFileDir(self,cel_file_dir): self._cel_file_dir = osfilepath(cel_file_dir)
    def setBEDFileDir(self,cel_file_dir): self._cel_file_dir = osfilepath(cel_file_dir)
    def setFeatureNormalization(self,normalize_feature_exp): self.normalize_feature_exp = normalize_feature_exp
    def setExcludeLowExpressionExons(self, excludeNonExpExons): self.excludeNonExpExons = excludeNonExpExons
    def setNormMatrix(self,normalize_gene_data): self.normalize_gene_data = normalize_gene_data
    def setProbabilityStatistic(self,probability_statistic): self.probability_statistic = probability_statistic
    def setFDRStatistic(self, FDR_statistic): self.FDR_statistic = FDR_statistic
    def setBatchEffectRemoval(self,batch_effects): self.batch_effects = batch_effects
    def setProducePlots(self,visualize_results): self.visualize_results = visualize_results
    def setPerformLineageProfiler(self, run_lineage_profiler): self.run_lineage_profiler = run_lineage_profiler
    def setCompendiumType(self,compendiumType): self.compendiumType = compendiumType
    def setCompendiumPlatform(self,compendiumPlatform): self.compendiumPlatform = compendiumPlatform
    def setClassificationAnalysis(self, classificationAnalysis): self.classificationAnalysis = classificationAnalysis
    def setReturnCentroids(self,returnCentroids): self.returnCentroids = returnCentroids
    def setMultiThreading(self, multithreading): self.multithreading = multithreading
    def setVendor(self,vendor): self.vendor = vendor
    def setPredictGroups(self, predictGroups): self.predictGroups = predictGroups
    def setPredictGroupsParams(self, predictGroupsObjects): self.predictGroupsObjects = predictGroupsObjects
    def setGraphicLinks(self,graphic_links): self.graphic_links = graphic_links ### file location of image files
    def setSTDOUT(self, stdout): self.stdout = stdout
    def setExonExpThreshold(self,exon_exp_threshold):
        try: exon_exp_threshold = float(exon_exp_threshold)
        except Exception: exon_exp_threshold = exon_exp_threshold
        self.exon_exp_threshold = exon_exp_threshold
    def setExonRPKMThreshold(self,exon_rpkm_threshold):
        try: exon_rpkm_threshold = float(exon_rpkm_threshold)
        except Exception: exon_rpkm_threshold = exon_rpkm_threshold
        self.exon_rpkm_threshold = exon_rpkm_threshold
    def setGeneExpThreshold(self,gene_exp_threshold):
        try: gene_exp_threshold = float(gene_exp_threshold)
        except Exception: gene_exp_threshold = gene_exp_threshold
        self.gene_exp_threshold = gene_exp_threshold
    def setJunctionExpThreshold(self,junction_exp_threshold):
        try: junction_exp_threshold = float(junction_exp_threshold)
        except Exception: junction_exp_threshold = junction_exp_threshold
        self.junction_exp_threshold = junction_exp_threshold
    def setRPKMThreshold(self,rpkm_threshold):
        try: rpkm_threshold = float(rpkm_threshold)
        except Exception: rpkm_threshold = rpkm_threshold
        self.rpkm_threshold = rpkm_threshold
    def setMarkerFinder(self,marker_finder): self.marker_finder = marker_finder
    def ReturnCentroids(self): return self.returnCentroids
    def FDRStatistic(self): return self.FDR_statistic
    def multiThreading(self): return self.multithreading
    def STDOUT(self): return self.stdout
    def ExonExpThreshold(self): return self.exon_exp_threshold
    def BatchEffectRemoval(self): return self.batch_effects
    def MarkerFinder(self): return self.marker_finder
    def PredictGroups(self): self.predictGroups
    def PredictGroupsObjects(self): self.predictGroupsObjects
    def ExonRPKMThreshold(self): return self.exon_rpkm_threshold
    def GeneExpThreshold(self): return self.gene_exp_threshold
    def JunctionExpThreshold(self): return self.junction_exp_threshold
    def RPKMThreshold(self): return self.rpkm_threshold
    def ProbabilityStatistic(self): return self.probability_statistic
    def ProducePlots(self): return self.visualize_results
    def PerformLineageProfiler(self): return self.run_lineage_profiler
    def CompendiumType(self): return self.compendiumType
    def CompendiumPlatform(self): return self.compendiumPlatform
    def ClassificationAnalysis(self): return self.classificationAnalysis
    def GraphicLinks(self): return self.graphic_links
    def setArrayType(self,array_type): self._array_type = array_type
    def setOutputDir(self,output_dir): self._output_dir = output_dir
    def setBiotypes(self,biotypes): self.biotypes = biotypes
    def setRootDir(self,parent_dir):
        ### Get directory above ExpressionInput
        split_dirs = string.split(parent_dir,'ExpressionInput')
        root_dir = split_dirs[0]
        self._root_dir = root_dir + '/'
    def setXHybRemoval(self,xhyb): self._xhyb = xhyb
    def XHybRemoval(self): return self._xhyb
    def setExonBedBuildStatus(self,bed_build_status): self.bed_build_status = bed_build_status
    def setRunKallisto(self, runKallisto): self.runKallisto = runKallisto
    def RunKallisto(self): return self.runKallisto
    def setCustomFASTA(self, customFASTA): self.customFASTA = customFASTA
    def CustomFASTA(self): return self.customFASTA
    def setChromiumSparseMatrix(self, chromiumSparseMatrix): self.chromiumSparseMatrix = chromiumSparseMatrix
    def ChromiumSparseMatrix(self): return self.chromiumSparseMatrix
    def setChannelToExtract(self,channel_to_extract): self.channel_to_extract = channel_to_extract
    def ExonBedBuildStatus(self): return self.bed_build_status
    def ChannelToExtract(self): return self.channel_to_extract
    def FeatureNormalization(self): return self.normalize_feature_exp
    def setUseJunctionsForGeneExpression(self, use_junctions_for_geneexpression): self.use_junctions_for_geneexpression = use_junctions_for_geneexpression
    def useJunctionsForGeneExpression(self):
        try: return self.use_junctions_for_geneexpression
        except Exception: return False
    def excludeLowExpressionExons(self): return self.excludeNonExpExons
    def NormMatrix(self): return self.normalize_gene_data
    def RootDir(self): return self._root_dir
    def APTLocation(self): return self._apt_location
    def InputCDFFile(self): return self._cdf_file
    def CLFFile(self): return self._clf_file
    def BGPFile(self): return self._bgp_file
    def CELFileDir(self): return self._cel_file_dir
    def BEDFileDir(self): return self._cel_file_dir+'/'
    def ArrayType(self):
        try: return self._array_type
        except Exception: return 'RNASeq'
    def OutputDir(self): return self._output_dir
    def Vendor(self):
        try: return self.vendor
        except Exception: return 'RNASeq'
    def setSpecies(self, species): self.species = species
    def Species(self): return self.species
    def setPlatformType(self, platformType): self.platformType = platformType
    def setAnalysisMode(self, analysis_mode): self.analysis_mode = analysis_mode
    def setMLP(self,mlpr): self.mlp = mlpr
    def setExonMapFile(self, exonMapFile): self.exonMapFile = exonMapFile
    def ExonMapFile(self): return self.exonMapFile
    def setCorrelationDirection(self, correlationDirection): self.correlationDirection = correlationDirection
    def CorrelationDirection(self): return self.correlationDirection
    def setPearsonThreshold(self, pearsonThreshold): self.pearsonThreshold = pearsonThreshold
    def PearsonThreshold(self): return self.pearsonThreshold
    def setUseAdjPvalue(self, useAdjPval): self.useAdjPval = useAdjPval
    def UseAdjPvalue(self):
        if string.lower(self.useAdjPval)=='false' or string.lower(self.useAdjPval)=='no' or self.useAdjPval==False:
            return False
        else:
            return True
    def setFoldCutoff(self, foldCutoff): self.foldCutoff = foldCutoff
    def FoldCutoff(self): return self.foldCutoff
    def setPvalThreshold(self, pvalThreshold): self.pvalThreshold = pvalThreshold
    def PvalThreshold(self): return self.pvalThreshold
    def setPeformDiffExpAnalysis(self, peformDiffExpAnalysis): self.peformDiffExpAnalysis = peformDiffExpAnalysis
    def PeformDiffExpAnalysis(self):
        if string.lower(self.peformDiffExpAnalysis)=='false' or string.lower(self.peformDiffExpAnalysis)=='no' or self.peformDiffExpAnalysis==False:
            return False
        else:
            return True
    def MLP(self): return self.mlp
    def PlatformType(self): return self.platformType
    def AnalysisMode(self): return self.analysis_mode
    def DatasetFile(self):
        if 'exp.' in self.ExpFile():
            dataset_dir = string.replace(self.ExpFile(),'exp.','DATASET-')
        else:
            parent = export.findParentDir(self.ExpFile())
            file = export.findFilename(self.ExpFile())
            if 'DATASET-' not in file:
                dataset_dir = parent + 'DATASET-'+file
            else:
                dataset_dir = self.ExpFile()
        dataset_dir = string.replace(dataset_dir,'ExpressionInput','ExpressionOutput')
        return dataset_dir
    def Architecture(self): return self.architecture
    def BioTypes(self): return self.biotypes
    def Report(self): return 'fl printout'+self.ExpFile()+'|'+str(len(self.StatsFile()))+'|'+str(len(self.GroupsFile()))+'|'+str(len(self.CompsFile()))
    def __repr__(self): return self.Report()

class AdditionalAlgorithms:
    def __init__(self, additional_algorithm):
        self._additional_algorithm = additional_algorithm
    def Algorithm(self): return self._additional_algorithm
    def setScore(self,score): self._score = score
    def Score(self): return self._score
    def __repr__(self): return self.Algorithm()
    
def getDirectoryFiles():
    status = 'repeat'
    while status == 'repeat':
        if backSelect == 'no' or 'InputExpFiles' == selected_parameters[-1]:
            root = Tk(); root.title('AltAnalyze: Select Expression File for Filtering')
            selected_parameters.append('InputExpFiles'); backSelect = 'no'
            gu = GUI(root,option_db,option_list['InputExpFiles'],'')
        else: gu = PreviousResults(old_options)
        try: input_exp_file = gu.Results()['input_exp_file']
        except KeyError: input_exp_file = '' ### Leave this blank so that the default directory is used
        try: input_stats_file = gu.Results()['input_stats_file']
        except KeyError: input_stats_file = '' ### Leave this blank so that the default directory is used
        #if array_type == 'exon':
        if 'steady-state' in input_exp_file or 'steady-state' in input_stats_file:
            print_out = "Do not select steady-state expression files.."
            IndicatorWindow(print_out,'Continue'); output_dir=''
        elif len(input_exp_file)>0:
            try: output_dir = gu.Results()['output_dir']
            except KeyError: output_dir = '' ### Leave this blank so that the default directory is used
            try: cel_files, array_linker_db = ExpressionBuilder.getArrayHeaders(input_exp_file)
            except Exception:
                print_out = "Input Expression file does not have a valid format."
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()  
            if len(cel_files)>0: status = 'continue'
            else:
                print_out = "The expression file:\n"+input_exp_file+"\ndoes not appear to be a valid expression file. Check to see that\nthis is the correct tab-delimited text file."
                IndicatorWindow(print_out,'Continue') 
        else:
            print_out = "No input expression file selected."
            IndicatorWindow(print_out,'Continue')
            
        if len(output_dir)<1:
            ### Set output to the same directory or parent if none selected
            if 'ExpressionInput' in input_exp_file: i = -2
            else: i = -1
            output_dir = string.join(string.split(input_exp_file,'/')[:i],'/')
            

def getUpdatedParameters(array_type,species,run_from_scratch,file_dirs):
    ### Get default options for ExpressionBuilder and AltAnalyze
    na = 'NA'; log = 'log'; no = 'no'
    global user_variables; user_variables={}; global selected_parameters; selected_parameters = []

    run_goelite=no; change_threshold=na;pathway_permutations=na;mod=na; ge_ptype = 'rawp';resources_to_analyze = na
    ge_fold_cutoffs=2;ge_pvalue_cutoffs=0.05;filter_method=na;z_threshold=1.96;p_val_threshold=0.05
    returnPathways = 'no'
    
    option_list,option_db = importUserOptions(array_type)
                                              
    global root
    if run_from_scratch != 'Prefiltered': ### This is when AltAnalyze has finished an analysis
        root = Tk()
        root.title('AltAnalyze: Perform Additional Analyses')
        selected_parameters.append('AdditionalOptions'); backSelect = 'no'
        gu = GUI(root,option_db,option_list['AdditionalOptions'],'')
        new_run = gu.Results()['new_run']
    else: new_run = None
    
    if new_run == 'Change Parameters and Re-Run': AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
    else:
        expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = importDefaults(array_type,species)
        option_db['get_additional'].setArrayOptions(['---']+importResourceList())
        option_db['get_additional'].setDefaultOption('---')
        default_resources = option_db['resources_to_analyze'].ArrayOptions()
        import_dir1 = '/AltDatabase/goelite/'+species+'/gene-mapp'
        import_dir2 = '/AltDatabase/goelite/'+species+'/gene-go'
        try:
            gene_mapp_list = read_directory(import_dir1)
            gene_mapp_list.sort()
            for file in gene_mapp_list:
                resource = string.split(file,'-')[-1][:-4]
                if resource != 'MAPP' and resource not in default_resources and '.txt' in file:
                    default_resources.append(resource)
        except Exception: pass
        try:
            gene_go_list = read_directory(import_dir2)
            gene_go_list.sort()
            for file in gene_go_list:
                resource = string.split(file,'-')[-1][:-4]
                if resource != 'GeneOntology' and resource not in default_resources and 'version' not in resource and '.txt' in file:
                    default_resources.append(resource)
        except Exception: pass
        option_db['resources_to_analyze'].setArrayOptions(default_resources)

        proceed = 'no'
        while proceed == 'no':
            root = Tk(); root.title('AltAnalyze: Pathway Analysis Parameters')
            if 'filtered' in run_from_scratch: ### Not relevant for 'Process AltAnalyze filtered'
                option_list['GOElite'] = option_list['GOElite'][3:]; goelite_defaults = goelite_defaults[3:]
            selected_parameters.append('GOElite'); backSelect = 'no'
            gu = GUI(root,option_db,option_list['GOElite'],goelite_defaults)
            if 'filtered' not in run_from_scratch: ### Not relevant for 'Process AltAnalyze filtered'
                ge_fold_cutoffs = gu.Results()['ge_fold_cutoffs']
                ge_pvalue_cutoffs = gu.Results()['ge_pvalue_cutoffs']
                ge_ptype = gu.Results()['ge_ptype']
            filter_method = gu.Results()['filter_method']
            z_threshold = gu.Results()['z_threshold']
            returnPathways = gu.Results()['returnPathways']
            p_val_threshold = gu.Results()['p_val_threshold']
            change_threshold = gu.Results()['change_threshold']
            resources_to_analyze = gu.Results()['resources_to_analyze']
            pathway_permutations = gu.Results()['pathway_permutations']
            ORA_algorithm = gu.Results()['ORA_algorithm']
            mod = gu.Results()['mod']
            get_additional = gu.Results()['get_additional']
            try:
                z_threshold = float(z_threshold)
                change_threshold = float(change_threshold)-1 ### This reflects the > statement in the GO-Elite filtering
                p_val_threshold = float(p_val_threshold)
                pathway_permutations = int(pathway_permutations)
                if run_from_scratch != 'Process AltAnalyze filtered':
                    ge_fold_cutoffs = float(ge_fold_cutoffs)
                    ge_pvalue_cutoffs = float(ge_pvalue_cutoffs)
                proceed = 'yes'
            except Exception:
                print_out = "Invalid numerical entry. Try again."
                IndicatorWindow(print_out,'Continue')
    if get_additional != '---':
        analysis = 'getAdditionalOnlineResources'
        values = species,get_additional
        StatusWindow(values,analysis) ### display an window with download status
    try:
        criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs
        import GO_Elite
        if run_from_scratch != 'Prefiltered': ### Only applies to AltAnalyze generated GO-Elite input
            ###Export dataset criterion using user-defined filters
            ExpressionBuilder.buildCriterion(ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, main_output_folder, 'goelite')
            #except Exception: null = []; # print 'No expression files to summarize'
        if ORA_algorithm == 'Fisher Exact Test':
            pathway_permutations = 'FisherExactTest'
        goelite_var = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,returnPathways,file_dirs,''
        GO_Elite.remoteAnalysis(goelite_var,'UI',Multi=mlp)
        AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
    except Exception:
        print traceback.format_exc()
        print_out = "Unexpected error encountered. Please see log file."
        IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
   
def addOnlineSpeciesDatabases(backSelect):
    StatusWindow(file_location_defaults,'getOnlineDBConfig')
    #except Exception,e: print [e]; null = []

    importSystemInfo()
    try: exportSystemInfo() ### By re-importing we incorporate new source data from the downloaded file
    except Exception:
        print 'Cannot write Config/source_data.txt to the Config directory (like Permissions Error)'
    existing_species_codes = species_codes
    importSpeciesInfo(); online_species = ['']
    for species in species_codes: online_species.append(species)
    online_species.sort()
    
    importOnlineDatabaseVersions(); db_version_list=[]
    for version in db_versions: db_version_list.append(version)
    db_version_list.sort(); db_version_list.reverse(); select_version = db_version_list[0]
    db_versions[select_version].sort()
    option_db['selected_species1'].setArrayOptions(['---']+db_versions[select_version])
    option_db['selected_species2'].setArrayOptions(['---']+db_versions[select_version])
    option_db['selected_species3'].setArrayOptions(['---']+db_versions[select_version])
    option_db['selected_version'].setArrayOptions(db_version_list)
    proceed = 'no'
    
    while proceed == 'no':
        if backSelect == 'no' or 'OnlineDatabases' == selected_parameters[-1]:
            selected_parameters.append('OnlineDatabases'); backSelect = 'no'
            root = Tk(); root.title('AltAnalyze: Species Databases Available for Download')
            gu = GUI(root,option_db,option_list['OnlineDatabases'],'')
        else: gu = PreviousResults(old_options); print 'alpha'
        db_version = gu.Results()['selected_version']
        exportDBversion(db_version)
        try: species1 = gu.Results()['selected_species1']
        except Exception: species1='---'
        try: species2 = gu.Results()['selected_species2']
        except Exception: species2='---'
        try: species3 = gu.Results()['selected_species3']
        except Exception: species3='---'
        try: species_full = gu.Results()['species']
        except Exception: species_full = ''
        try: update_goelite_resources = gu.Results()['update_goelite_resources']
        except Exception: update_goelite_resources = ''
        #if species_full == 'Add Species': AltAnalyze.AltAnalyzeSetup(species_full); sys.exit()
        new_species_list = [species1,species2,species3]; new_species_codes={}
        for species in new_species_list:
            if '---' not in species:
                #try:
                ### Export basic species information
                sc = species_codes[species].SpeciesCode()
                existing_species_codes[species] = species_codes[species]
                new_species_codes[sc]=[]
                #except Exception: sc = None
            if sc != None:
                for ad in db_versions_vendors[db_version]:
                    if ad.SpeciesCodes() == species:
                        for array_system in array_codes:
                            ac = array_codes[array_system]
                            compatible_species = ac.SpeciesCodes()
                            if ac.Manufacturer() in ad.Manufacturer() and ('expression' in ac.ArrayName() or 'RNASeq' in ac.ArrayName() or 'RNA-seq' in ac.ArrayName()):
                                if sc not in compatible_species: compatible_species.append(sc)
                            ac.setSpeciesCodes(compatible_species)
                try: exportArrayInfo(array_codes)
                except Exception:
                    print 'Cannot write Config/arrays.txt to the Config directory (like Permissions Error)'
        if len(new_species_codes) > 0:
            analysis = 'getOnlineEliteDatabase'
            values = file_location_defaults,db_version,new_species_codes,update_goelite_resources ### Download the online databases
            StatusWindow(values,analysis)
            proceed = 'yes'
        else:
            print_out = "Please select a species before continuing."
            IndicatorWindow(print_out,'Try Again')
    
    #db_versions_vendors
    try: exportSpeciesInfo(existing_species_codes)
    except Exception:
        print 'Cannot write Config/species.txt to the Config directory (like Permissions Error)'
    integrate_online_species = 'no'

def getArraysAndVendors(species,vendor):
    array_list2=[]; manufacturer_list=[]
    compatible_species,manufacturer_list_all = getSpeciesList('')
    for array_name in array_list:
        manufacturer = array_codes[array_name].Manufacturer()
        if species in array_codes[array_name].SpeciesCodes():
            manufacturer_list.append(manufacturer)
            if len(vendor)>0:
                if vendor == manufacturer: proceed = 'yes'
                else: proceed = 'no'
            else: proceed = 'yes'
            if proceed == 'yes':
                array_list2.append(array_name)
                
    manufacturer_list = unique.unique(manufacturer_list)
    array_list2 = unique.unique(array_list2) ### Filtered based on compatible species arrays
    array_list2.sort(); manufacturer_list.sort()
    if vendor == 'RNASeq':
        array_list2.reverse()
    return array_list2, manufacturer_list

def getSpeciesForArray(array_type):
    array_list2=[]; manufacturer_list=[]; manufacturer_list_all=[]
    for array_name in array_list:
        current_species_codes = array_codes[array_type].SpeciesCodes()

    try: current_species_dirs = unique.read_directory('/AltDatabase')
    except Exception: current_species_dirs = current_species_codes
    
    current_species_names=[]
    for species in species_codes:
        species_code = species_codes[species].SpeciesCode()
        if species_code in current_species_codes:
            if species_code in current_species_dirs: current_species_names.append(species)
    current_species_names.sort()
    return current_species_names

def verifyLineageProfilerDatabases(species,run_mode):
    import AltAnalyze
    installed = False
    download_species = species
    try: gene_database = unique.getCurrentGeneDatabaseVersion()
    except Exception: gene_database='00'
    try:
        if int(gene_database[-2:]) < 62:
            print_out = 'LineageProfiler is not supported in this database version (EnsMart62 and higher required).'
            print print_out
            return False
        else:
            if species == 'Hs':
                source_file = 'AltDatabase/ensembl/'+species+'/'+species+'_exon_tissue-specific_protein_coding.txt'
                download_species = 'Hs'
            elif species == 'Mm':
                source_file = 'AltDatabase/ensembl/'+species+'/'+species+'_gene_tissue-specific_protein_coding.txt'
                download_species = 'Mm'
            else: ### Use the mouse version instead - less variable data
                source_file = 'AltDatabase/ensembl/'+species+'/'+species+'_gene_tissue-specific_protein_coding.txt'
                download_species = 'Mm'
            file_length = AltAnalyze.verifyFileLength(source_file)
            if file_length>0:
                installed = True
            else:
                print_out = 'To perform a LineageProfiler analysis AltAnalyze must\nfirst download the appropriate database.'
                if run_mode == 'GUI':
                    IndicatorWindow(print_out,'Download')
                else:
                    print print_out  ### Occurs in command-line mode
                filename = 'AltDatabase/ensembl/'+download_species+'_LineageProfiler.zip'
                dir = 'AltDatabase/updated/'+gene_database ### Directory at altanalyze.org
                var_list = filename,dir
                if debug_mode == 'no' and run_mode == 'GUI': StatusWindow(var_list,'download')
                else: update.downloadCurrentVersion(filename,dir,None)
                file_length = AltAnalyze.verifyFileLength(source_file)
                if file_length>0: installed = True
                else:
                    try:
                        from build_scripts import GeneSetDownloader
                        GeneSetDownloader.translateBioMarkersBetweenSpecies('AltDatabase/ensembl/'+download_species,species)
                    except Exception:
                        None
    except Exception: installed = False
    return installed

def checkForLocalArraySupport(species,array_type,specific_arraytype,run_mode):
    specific_arraytype = string.lower(specific_arraytype) ### Full array name
    if array_type == 'junction' or array_type == 'RNASeq':
        try: gene_database = unique.getCurrentGeneDatabaseVersion()
        except Exception: gene_database='00'
        if int(gene_database[-2:]) < 0:
            print_out = 'The AltAnalyze database indicated for '+array_type+' analysis\n is not supported for alternative exon analysis.\nPlease update to EnsMart55 or greater before\nproceeding.'
            if run_mode == 'GUI': IndicatorWindow(print_out,'Continue')
            else: print print_out ### Occurs in command-line mode
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        downloaded_junction_db = 'no'; file_problem='no'; wrong_junction_db = 'no'
        while downloaded_junction_db == 'no': ### Used as validation in case internet connection is unavailable
            try: dirs = read_directory('/AltDatabase/'+species)
            except Exception: dirs=[]
            if wrong_junction_db == 'yes':
                print_out = 'Another junction database is installed. Select "Contine" to overwrite or manually change the name of this folder:\n'+filepath('AltDatabase/'+species+'/'+array_type)
                if run_mode == 'GUI': IndicatorWindow(print_out,'Continue')
                else: print print_out  ### Occurs in command-line mode
            if array_type not in dirs or file_problem == 'yes' or wrong_junction_db == 'yes':
                if file_problem == 'yes':
                    print_out = 'Unknown installation error occured.\nPlease try again.'
                else:
                    print_out = 'To perform an '+array_type+' analysis allow AltAnalyze \nto download the appropriate database now.'
                if run_mode == 'GUI': IndicatorWindow(print_out,'Download')
                else: print print_out  ### Occurs in command-line mode
                if array_type == 'RNASeq': filename = 'AltDatabase/'+species+'_'+array_type+'.zip'
                elif 'glue' in specific_arraytype: filename = 'AltDatabase/'+species+'/'+species+'_'+array_type+'_Glue.zip'
                elif 'hta 2.0' in specific_arraytype: filename = 'AltDatabase/'+species+'/'+species+'_'+array_type+'_HTA-2_0.zip'
                elif 'mta 1.0' in specific_arraytype: filename = 'AltDatabase/'+species+'/'+species+'_'+array_type+'_MTA-1_0.zip'
                else: filename = 'AltDatabase/'+species+'/'+species+'_'+array_type+'.zip'
                dir = 'AltDatabase/updated/'+gene_database; var_list = filename,dir
                if debug_mode == 'no' and run_mode == 'GUI':
                    StatusWindow(var_list,'download')
                else: update.downloadCurrentVersion(filename,dir,None)
            try: dirs = read_directory('/AltDatabase/'+species)
            except Exception: dirs=[]
            if array_type in dirs:
                import AltAnalyze
                file_length = AltAnalyze.verifyFileLength('AltDatabase/'+species+'/'+array_type+'/probeset-domain-annotations-exoncomp.txt')
                if file_length>0: downloaded_junction_db = 'yes'
                elif species == 'Mm' or species == 'Hs' or species == 'Rn': file_problem = 'yes'
                else: downloaded_junction_db = 'yes' ### Occurs when no alternative exons present for species
                if array_type == 'junction':
                    specific_platform = determinePlatform('AltDatabase/'+species+'/'+array_type+'/platform.txt')
                    if 'glue' in specific_arraytype and 'Glue' not in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    elif 'glue' not in specific_arraytype and 'Glue' in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    elif 'hta 2.0' in specific_arraytype and 'HTA-2_0' not in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    elif 'hta 2.0' not in specific_arraytype and 'HTA-2_0' in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    elif 'mta 1.0' in specific_arraytype and 'MTA-1_0' not in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    elif 'mta 1.0' not in specific_arraytype and 'MTA-1_0' in specific_platform: wrong_junction_db = 'yes'; downloaded_junction_db = 'no'
                    #print [specific_arraytype], [specific_platform], wrong_junction_db, downloaded_junction_db
      
def exportGeneList(gene_list,outputFolder):
    filename = string.join(gene_list,' ')[:25]
    eo = export.ExportFile(outputFolder+'/GO-Elite_input/'+filename+'.txt')
    eo.write('Symbol\tSytemCode\n')
    for i in gene_list:
        eo.write(i+'\tSy\n')
    return outputFolder+'/GO-Elite_input'

def getUserParameters(run_parameter,Multi=None):
    global AltAnalyze; import AltAnalyze; global mlp; mlp=Multi ### multiprocessing support
    if run_parameter == 'yes':
        try: MainMenu()
        except Exception:
            print traceback.format_exc()
            print_out = "\nCritical error encountered!!! This machine does not have either:\n"
            print_out += "1) Have the required Tcl/Tk components installed.\n"
            print_out += "2) Is being run from a compiled version that has critical incompatibilities your OS or hardware or\n"
            print_out += "3) Is being run from source-code in the same-directory as executable code resulting in a conflict\n"
            print_out += "\nIf any of these apply, we recommend downloading the Python source-code version of AltAnalyze "
            print_out += "(installing necessary dependencies - see our Wiki or Documentation)."
            print_out += "Otherwise, please contact AltAnalyze support (http://code.google.com/p/altanalyze/wiki/ContactUs).\n\n"
            print_out += "Installation Wiki: http://code.google.com/p/altanalyze/wiki/Installation\n\n"
            print print_out
            try:
                ### Create a log report of this
                try: log_file = filepath('AltAnalyze_error-report.log')
                except Exception: log_file = filepath('/AltAnalyze_error-report.log')
                log_report = open(log_file,'w');
                log_report.write(print_out)
                log_report.write(traceback.format_exc())
                log_report.close()
                ### Open this file
                if os.name == 'nt':
                    try: os.startfile('"'+log_file+'"')
                    except Exception:  os.system('open "'+log_file+'"')
                elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
                elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
            except Exception: None
            sys.exit()
    global species; species=''; global user_variables; user_variables={}; global analysis_method; global array_type; global vendor
    global PathDir; global PathFile; global file_location_defaults; global integrate_online_species; integrate_online_species = 'no'
    global option_db; global option_list; global analysis_status; analysis_status = 'continue'; global selected_parameters; selected_parameters=[]
    global backSelect; global fl; predictGroups = False

    if os.name == 'posix' and run_parameter == 'yes':
        try: MacConsiderations()
        except Exception:
            print traceback.format_exc()
            sys.exit()
    ### Get default options for ExpressionBuilder and AltAnalyze

    na = 'NA'; log = 'log'; no = 'no'
    run_from_scratch=na; expression_threshold=na; perform_alt_analysis=na; expression_data_format=log
    include_raw_data=na; avg_all_for_ss=no; dabg_p=na; normalize_feature_exp=na; normalize_gene_data = na
    analysis_method=na; p_threshold=na; filter_probeset_types=na; alt_exon_fold_cutoff=na
    permute_p_threshold=na; perform_permutation_analysis=na; export_splice_index_values=no
    run_MiDAS=no; analyze_functional_attributes=no; microRNA_prediction_method=na
    gene_expression_cutoff=na; cel_file_dir=na; input_exp_file=na; input_stats_file=na; filter_for_AS=no
    remove_intronic_junctions=na; build_exon_bedfile=no; input_cdf_file = na; bgp_file = na
    clf_file = na; remove_xhyb = na; multiThreading = True; input_fastq_dir = ''; sparse_matrix_file=''
    
    compendiumType = 'protein_coding'; compendiumPlatform = 'gene'
    calculate_splicing_index_p=no; run_goelite=no; ge_ptype = 'rawp'; probability_algorithm = na
    ge_fold_cutoffs=2;ge_pvalue_cutoffs=0.05;filter_method=na;z_threshold=1.96;p_val_threshold=0.05
    change_threshold=2;pathway_permutations=na;mod=na; analyze_all_conditions=no; resources_to_analyze=na
    additional_algorithms = na; rpkm_threshold = na; exon_exp_threshold = na; run_lineage_profiler = no
    gene_exp_threshold = na; exon_rpkm_threshold = na; visualize_results = no; returnPathways = 'no'
    batch_effects = na; marker_finder = na
                
    try: option_list,option_db = importUserOptions('exon')  ##Initially used to just get the info for species and array_type
    except IOError:
        ### Occurs if Config folder is absent or when the source code is run outside AltAnalyze root
        print_out = '\nWarning! The Config folder in the AltAnalyze program directory cannot be found. The likely cause is:\n'
        print_out +='   A): The AltAnalyze source-code is being run outside the root AltAnalyze directory or   \n'
        print_out +='   B): AltAnalyze was zip extracted/installed in a weird way (incommpatible zip extractor)\n'
        print_out +='\nIf you beleive (B) is possible, unzip with another unzip program (e.g., default Windows unzip program).'
        print_out +='\nIf neither applies, we recommend contacting our help desk (http://code.google.com/p/altanalyze/wiki/ContactUs).'

        try: IndicatorWindow(print_out,'Exit')
        except Exception: print printout
        sys.exit()
        
    importSpeciesInfo()
    file_location_defaults = importDefaultFileLocations()
    importArrayInfo()

    try: elite_db_versions = returnDirectoriesNoReplace('/AltDatabase')
    except Exception:
        try: elite_db_versions=[]; os.mkdir(filepath('AltDatabase'))
        except Exception: pass ### directory already exists      
    try: gene_database_dir = unique.getCurrentGeneDatabaseVersion()
    except Exception: gene_database_dir=''
    if len(elite_db_versions)>0 and gene_database_dir == '':
        for db_version in elite_db_versions:
            if 'EnsMart' in db_version:
                gene_database_dir = db_version; exportDBversion(db_version)
    
    current_species_names,manufacturer_list_all = getSpeciesList('')
    option_db['species'].setArrayOptions(current_species_names)

    try: PathDir = file_location_defaults['PathDir'].Location()
    except Exception:
        try:
            ### Entry was deleted from Config file - re-create it
            fl = FileLocationData('local', '', 'all')
            file_location_defaults['PathDir'] = fl
        except Exception: null = None   
        PathDir = ''
    try: PathFile = file_location_defaults['PathFile'].Location()
    except Exception:
        try:
            ### Entry was deleted from Config file - re-create it
            fl = FileLocationData('local', '', 'all')
            file_location_defaults['PathFile'] = fl
        except Exception: null = None        
        PathFile = ''
    
    old_options = []
    try:
        #### Get information from previous loop
        if len(run_parameter) == 2 and run_parameter != 'no': ### Occurs when selecting "Back" from Elite parameter window
            old_options = run_parameter[1]; selected_parameters = run_parameter[0]
            try:
                if selected_parameters[-2]==selected_parameters[-1]: selected_parameters = selected_parameters[:-1] 
            except Exception: selected_parameters = selected_parameters
            backSelect = 'yes'
            #print selected_parameters
            #print old_options,'\n'
            for option in old_options: ### Set options to user selected
                try: option_db[option].setDefaultOption(old_options[option])
                except Exception: pass
                user_variables[option] = old_options[option]
            if 'array_type' in old_options:
                specific_array = old_options['array_type']
                vendor = old_options['manufacturer_selection']
                species_full = old_options['species']
                species = species_codes[species_full].SpeciesCode()
            if selected_parameters == []: backSelect = 'no'
        else: backSelect = 'no'; old_options=[]
        
        ###Update this informatin in option_db which will be over-written after the user selects a species and array_type
        option_db['species'].setArrayOptions(current_species_names)
        if len(current_species_names)==0 and run_parameter != 'Add Species':
            print_out = "No species databases found. Select\ncontinue to proceed with species download."
            IndicatorWindow(print_out,'Continue')
            integrate_online_species = 'yes'
            addOnlineSpeciesDatabases(backSelect)
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        
        ### Set defaults based on avialable species
        
        #default_vendor = 'RNASeq'
        #default_specific_array = 'RNA-seq aligned read counts'
        default_vendor = 'RNASeq'
        default_specific_array='Raw sequence or processed'
        
        """
        try: ### If the users have already analyzed Affy data, make this the default
            affymetrix_library_dir = 'AltDatabase/affymetrix/LibraryFiles'
            affy_dir_list = read_directory(filepath(affymetrix_library_dir))
            if len(affy_dir_list)>0:
                default_vendor = 'Affymetrix'
                default_specific_array='Affymetrix expression array'
        except Exception:
            None ### Occurs if this directory is missing (possible in future versions)
        """
        
        if run_parameter == 'Add Species':
            species_full = 'Homo sapiens'; species = 'Hs'; vendor = 'Affymetrix'; specific_array = 'Exon 1.0 ST array'
        if backSelect == 'yes' and 'array_type' in old_options:
            pass
        elif 'Homo sapiens' in current_species_names:
            species_full = 'Homo sapiens'; species = 'Hs'; vendor = default_vendor; specific_array = default_specific_array
        elif 'Mus musculus' in current_species_names:
            species_full = 'Mus musculus'; species = 'Mm'; vendor = default_vendor; specific_array = default_specific_array
        elif 'Rattus norvegicus' in current_species_names:
            species_full = 'Rattus norvegicus'; species = 'Rn'; vendor = default_vendor; specific_array = default_specific_array
        else:
            for species_full in current_species_names:
                species = species_codes[species_full].SpeciesCode()
                for array_name in array_list:
                    vendor = array_codes[array_name].Manufacturer()
                    if species in array_codes[array_name].SpeciesCodes(): specific_array = array_name; break
        array_list2, manufacturer_list = getArraysAndVendors(species,vendor)
        #print [[array_list2]], species, vendor
        option_db['species'].setDefaultOption(species_full)
        option_db['array_type'].setArrayOptions(array_list2)
        option_db['array_type'].setDefaultOption(specific_array)
        option_db['manufacturer_selection'].setArrayOptions(manufacturer_list_all)
        option_db['manufacturer_selection'].setDefaultOption(vendor)
        
        manufacturer_list_all_possible=[]        
        for array_name in array_list:
            manufacturer = array_codes[array_name].Manufacturer(); manufacturer_list_all_possible.append(manufacturer) 
        manufacturer_list_all_possible = unique.unique(manufacturer_list_all_possible); manufacturer_list_all_possible.sort()
        
        if len(elite_db_versions)>1:
            option_db['dbase_version'].setArrayOptions(elite_db_versions)
            option_db['dbase_version'].setDefaultOption(gene_database_dir)
        else:
            ### Otherwise, remove this option
            del option_db['dbase_version']
        
        ### Get user array and species selections
        if run_parameter != 'Add Species':
            if backSelect == 'no' or 'ArrayType' == selected_parameters[-1]:
                selected_parameters.append('ArrayType'); backSelect = 'no'
                root = Tk(); root.title('AltAnalyze: Select Species and Experimental Platform')
                gu = GUI(root,option_db,option_list['ArrayType'],'')
            else: gu = PreviousResults(old_options)
            species_full = gu.Results()['species']

        new_analysis_options=[]
        try: update_dbs = gu.Results()['update_dbs']
        except Exception: update_dbs = 'no'
        
        try:
            selected_parameters[-1]
        except Exception:
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
            
        if update_dbs == 'yes' or species_full == 'Add Species' or 'NewSpecies' == selected_parameters[-1]:
            integrate_online_species = 'yes'
            addOnlineSpeciesDatabases(backSelect)
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()

        elif species_full == 'Add Species' or 'NewSpecies' == selected_parameters[-1]: ### outdated code - bypassed by the above
            species_added = 'no'
            option_db['new_manufacturer'].setArrayOptions(manufacturer_list_all_possible)
            while species_added == 'no':
                if backSelect == 'no' or 'NewSpecies' == selected_parameters[-1]:
                    selected_parameters.append('NewSpecies'); backSelect = 'no'
                    root = Tk(); root.title('AltAnalyze: Add New Species Support')
                    gu = GUI(root,option_db,option_list['NewSpecies'],'')
                else: gu = PreviousResults(old_options)
                new_species_code = gu.Results()['new_species_code']
                new_species_name = gu.Results()['new_species_name']
                new_manufacturer = gu.Results()['new_manufacturer']

                if len(new_species_code)==2 and len(new_species_name)>0 and len(new_manufacturer)>0:
                    species_added = 'yes'
                    sd = SpeciesData(new_species_code,new_species_name,[''])
                    species_codes[new_species_name] = sd
                    try: exportSpeciesInfo(species_codes)
                    except Exception:
                        print 'Cannot write Config/species.txt to the Config directory (like Permissions Error)'
                    try: os.mkdir(filepath('AltDatabase/'+new_species_code))
                    except Exception: pass
                    for array_system in array_codes:
                        ac = array_codes[array_system]
                        manufacturer=ac.Manufacturer()
                        compatible_species = ac.SpeciesCodes()
                        if manufacturer == new_manufacturer and 'expression array' in ac.ArrayName():
                            if new_species_code not in compatible_species: compatible_species.append(new_species_code)
                        ac.setSpeciesCodes(compatible_species)
                        try: exportArrayInfo(array_codes)
                        except Exception:
                            print 'Cannot write Config/arrays.txt to the Config directory (like Permissions Error)'
                    fn = filepath('AltDatabase/affymetrix/'+ new_species_code)
                    try: os.mkdir(fn)
                    except OSError: null = [] ### Directory already exists
                    AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
                else:
                    print_out = "Valid species data was not added. You must\nindicate a two letter species code and full species name."
                    IndicatorWindow(print_out,'Continue')  
        else: species = species_codes[species_full].SpeciesCode()    
        array_full = gu.Results()['array_type'] ### Can be 10X Genomics
        vendor = gu.Results()['manufacturer_selection']
        if '10X' in array_full:
            vendor = '10XGenomics'

        try:
            array_type = array_codes[array_full].ArrayCode() ### Here, 10X Genomics would be converted to RNASeq
        except Exception:
            if vendor == 'Other ID':
                #"""
                ### An error occurs because this is a system name for the Other ID option
                array_type = "3'array"
                if array_full == "3'array" and vendor == 'RNASeq':
                    ### Occurs when hitting the back button
                    ### When RNASeq is selected as the platform but change to "3'array" when normalized data is imported.
                    vendor = 'other:Symbol'
                else:
                    vendor = 'other:'+array_full ### Ensembl linked system name

        if array_type == 'gene':
            try: gene_database = unique.getCurrentGeneDatabaseVersion()
            except Exception: gene_database='00'
            if int(gene_database[-2:]) < 54:
                print_out = 'The AltAnalyze database indicated for Gene 1.0 ST\narray analysis is not supported for alternative exon\nanalysis. Please update to EnsMart54 or greater\nbefore proceeding.'
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()

        ### Examine the AltDatabase folder for directories required for specific array analyses      
        checkForLocalArraySupport(species,array_type,array_full,'GUI')
        
        if array_type == 'exon' or array_type == 'AltMouse' or array_type == 'gene' or array_type == 'junction':
            try: dirs = read_directory('/AltDatabase/'+species)
            except Exception: dirs=[]
            if len(dirs)==0:
                print_out = 'Valid database directories were not found for this array.\nPlease re-install database.'
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()

        if '10X' in vendor:
            ### Needed when the the back button is selected for the 10X platform
            array_type = '10XGenomics'
            array_full == '10X Genomics sparse matrix'
        option_list,option_db = importUserOptions(array_type,vendor=vendor)  ##Initially used to just get the info for species and array_type
        
        if array_type == "3'array" and '10X' not in vendor:
            if species == 'Hs': compendiumPlatform = "3'array"
            for i in option_db['run_from_scratch'].ArrayOptions():
                if 'AltAnalyze' not in i:
                    if array_type == "3'array":
                        if 'CEL' in i and vendor != 'Affymetrix': proceed = 'no'
                        else: proceed = 'yes'
                    else: proceed = 'yes'
                    if proceed == 'yes': new_analysis_options.append(i)
            option_db['run_from_scratch'].setArrayOptions(new_analysis_options)

        proceed = 'no'
        if len(new_analysis_options)!=1:
            if backSelect == 'no' or 'AnalysisType' == selected_parameters[-1]:
                selected_parameters.append('AnalysisType'); backSelect = 'no'
                root = Tk(); root.title('AltAnalyze: Select Analysis Method')
                gu = GUI(root,option_db,option_list['AnalysisType'],'')
            else: gu = PreviousResults(old_options)
            run_from_scratch = gu.Results()['run_from_scratch']
        else: run_from_scratch = 'Process Expression file'
        try: vendor = array_codes[array_full].Manufacturer()
        except Exception: None ### Key the existing vendor
        try: constitutive_source = array_codes[array_full].ConstitutiveSource()
        except Exception: constitutive_source = vendor

        if '10X' in array_full or '10X' in vendor:
            array_type = "3'array"
            vendor = 'other:'+array_full ### Ensembl linked system name
            #option_list,option_db = importUserOptions(array_type,vendor=vendor)  ##Initially used to just get the info for species and array_type

        if backSelect == 'yes':
            for option in old_options: ### Set options to user selected
                try: option_db[option].setDefaultOption(old_options[option])
                except Exception: pass
                
        if run_from_scratch == 'Interactive Result Viewer':
            AltAnalyze.AltAnalyzeSetup('remoteViewer');sys.exit()
            
        def rebootAltAnalyzeGUI(selected_parameters,user_variables):
            commandline_args = ['--selected_parameters',selected_parameters[-1]]
            for uv in user_variables:
                if isinstance(user_variables[uv], list):
                    commandline_args += ['--'+uv,user_variables[uv][0]]
                else:
                    try:
                        if len(user_variables[uv])>0:
                            commandline_args += ['--'+uv,user_variables[uv]]
                    except Exception: pass 
            commandline_args = map(lambda x: string.replace(x,' ','__'),commandline_args)
            commandline_args = str(string.join(commandline_args,' '))
            if os.name == 'posix' or os.name == 'nt':
                try:
                    package_path = filepath('python')
                    if os.name == 'posix':
                        package_path = string.replace(package_path,'python','AltAnalyze.app/Contents/MacOS/AltAnalyze')
                    else:
                        package_path = string.replace(package_path,'python','AltAnalyze.exe')
                        package_path = 'AltAnalyze.exe'
                    #print [package_path+' --GUI yes '+commandline_args]
                    os.system(package_path+' --GUI yes '+commandline_args);sys.exit()
                except Exception:   
                    package_path = filepath('python')
                    package_path = string.replace(package_path,'python','AltAnalyze.py')
                    package_path = 'python '+package_path
                    if os.name == 'nt':
                        package_path = 'python AltAnalyze.py'
                    os.system(package_path+' --GUI yes '+commandline_args);sys.exit()
            else:
                AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                                
        if run_from_scratch == 'Additional Analyses':

            if backSelect == 'no' or 'Additional Analyses' == selected_parameters[-1]:
                selected_parameters.append('Additional Analyses'); backSelect = 'no'
                root = Tk()
                root.title('AltAnalyze: Additional Analysis Options')
                gu = GUI(root,option_db,option_list['Additional Analyses'],'')
                ### Venn Diagram Error here with _tkinter.TclError: image "pyimage36" doesn't exist
            else: gu = PreviousResults(old_options)
            additional_analyses = gu.Results()['additional_analyses']

            if 'nrichment' in additional_analyses:
                status = 'repeat'
                while status == 'repeat':
                    if backSelect == 'no' or 'InputGOEliteDirs' == selected_parameters[-1]:
                        root = Tk(); root.title('AltAnalyze: Select Expression File for Filtering')
                        selected_parameters.append('InputGOEliteDirs'); backSelect = 'no'
                        gu = GUI(root,option_db,option_list['InputGOEliteDirs'],'')
                    else: gu = PreviousResults(old_options)
                    try: criterion_input_folder = gu.Results()['criterion_input_folder']
                    except KeyError: criterion_input_folder = '' ### Leave this blank so that the default directory is used
                    try: criterion_denom_folder = gu.Results()['criterion_denom_folder']
                    except KeyError: criterion_denom_folder = '' ### Leave this blank so that the default directory is used
                    try:
                        try: main_output_folder = gu.Results()['main_output_folder']
                        except KeyError: main_output_folder = 'GO-Elite/input/' ### Leave this blank so that the default directory is
                        inputIDs = gu.Results()['inputIDs']
                        if len(inputIDs)>0:
                            inputIDs = string.replace(inputIDs, '\r',' ')
                            inputIDs = string.replace(inputIDs, '\n',' ')
                            inputIDs = string.split(inputIDs, ' ')
                            criterion_input_folder = exportGeneList(inputIDs,main_output_folder)
                    except Exception: inputIDs=[]
                    
                    if len(criterion_input_folder)>0:# and len(criterion_denom_folder)>0:
                        try: main_output_folder = gu.Results()['main_output_folder']
                        except KeyError: main_output_folder = '' ### Leave this blank so that the default directory is
                        if len(main_output_folder)<1:
                            ### Set output to the same directory or parent if none selected
                            i = -1 ### 1 directory up
                            main_output_folder = string.join(string.split(criterion_input_folder,'/')[:i],'/')
                        status = 'continue'
                    else:
                        print_out = "No GO-Elite input or denominator folder(s) selected."
                        IndicatorWindow(print_out,'Continue')
                        
                file_dirs = criterion_input_folder, criterion_denom_folder, main_output_folder
                #print file_dirs
                ### Get GO-Elite Input Parameters
                getUpdatedParameters(array_type,species,'Prefiltered',file_dirs)

            if additional_analyses == 'Pathway Visualization':
                root = Tk()
                root.title('AltAnalyze: Visualize Data on WikiPathways')
                selected_parameters.append('Pathway Visualization')
                GUI(root,'ViewWikiPathways',[],'') ### The last is default attributes (should be stored as defaults in the option_db var)

            if additional_analyses == 'Identifier Translation':

                try:
                    selected_parameters.append('Identifier Translation')
                    
                    supported_geneid_types = getSupportedGeneSystems(species,'uid-gene')
                    option_db['input_source'].setArrayOptions(['None Selected']+supported_geneid_types)
                    option_db['output_source'].setArrayOptions(['None Selected']+supported_geneid_types)
                    #option_db['PathwaySelection'].setArrayOptions(supported_genesets)
                except Exception,e:
                    print traceback.format_exc()
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Translate Input File Identifiers to Another System')
                    gu = GUI(root,option_db,option_list['IDConverter'],'')
                    try: input_cluster_file = gu.Results()['input_cluster_file']
                    except Exception: input_cluster_file = ''
                    input_data_file = gu.Results()['input_data_file']
                    input_source = gu.Results()['input_source']
                    output_source = gu.Results()['output_source']
                    
                    if len(input_data_file)>0 and input_source != 'None Selected' and output_source != 'None Selected':
                        analysis = 'IDConverter'
                        values = input_data_file, species, input_source, output_source
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')

            if additional_analyses == 'Merge Files':
                selected_parameters.append('Merge Files')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Merge Multiple Text Files Containing Common IDs')
                    gu = GUI(root,option_db,option_list['MergeFiles'],'')
                    input_file1 = gu.Results()['input_file1']
                    input_file2 = gu.Results()['input_file2']
                    input_file3 = gu.Results()['input_file3']
                    input_file4 = gu.Results()['input_file4']
                    join_option = gu.Results()['join_option']
                    ID_option = gu.Results()['ID_option']
                    output_merge_dir = gu.Results()['output_merge_dir']
                    
                    if len(input_file1)>0 and len(input_file2)>0 and len(output_merge_dir)>0:
                        if ID_option == 'False': ID_option = False
                        if ID_option == 'True': ID_option = True
                        analysis = 'MergeFiles'
                        files_to_merge = [input_file1, input_file2]
                        if len(input_file3)>0: files_to_merge.append(input_file3)
                        if len(input_file4)>0: files_to_merge.append(input_file4)
                        values = files_to_merge, join_option, ID_option, output_merge_dir
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')

            if additional_analyses == 'Venn Diagram':
                selected_parameters.append('Venn Diagram')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: View Venn Diagram from AltAnalyze or Input Files')
                    gu = GUI(root,option_db,option_list['VennDiagram'],'')
                    input_file1 = gu.Results()['venn_input_file1']
                    input_file2 = gu.Results()['venn_input_file2']
                    input_file3 = gu.Results()['venn_input_file3']
                    input_file4 = gu.Results()['venn_input_file4']
                    venn_output_dir = gu.Results()['venn_output_dir']
                    
                    if len(input_file1)>0 and len(input_file2)>0 and len(venn_output_dir)>0:
                        analysis = 'VennDiagram'
                        files_to_merge = [input_file1, input_file2]
                        if len(input_file3)>0: files_to_merge.append(input_file3)
                        if len(input_file4)>0: files_to_merge.append(input_file4)
                        values = files_to_merge, venn_output_dir
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')

            if additional_analyses == 'AltExon Viewer':
                selected_parameters.append('AltExon Viewer')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Visualize Exon-Level Expression Results')
                    gu = GUI(root,option_db,option_list['AltExonViewer'],'')
                    altanalyze_results_folder = gu.Results()['altanalyze_results_folder']
                    data_type = gu.Results()['data_type']
                    show_introns = gu.Results()['show_introns']
                    gene_symbol = gu.Results()['gene_symbol']
                    altgenes_file = gu.Results()['altgenes_file']
                    analysisType = gu.Results()['analysisType']
                    if len(altgenes_file)>0 and analysisType != 'Sashimi-Plot':
                        gene_symbol = importGeneList(altgenes_file) ### list of gene IDs or symbols
                    if analysisType == 'Sashimi-Plot':
                        altanalyze_results_folder = string.split(altanalyze_results_folder,'AltResults')[0]
                        exp_file = altanalyze_results_folder
                        if len(gene_symbol)<1:
                            gene_symbol = altgenes_file
                    elif data_type == 'raw expression': ### Switch directories if expression
                        altanalyze_results_folder = string.replace(altanalyze_results_folder,'AltResults','ExpressionInput')
                        exp_file = getValidExpFile(altanalyze_results_folder)
                    else:
                        altanalyze_results_folder += '/RawSpliceData/'+species
                        try: exp_file = getValidSplicingScoreFile(altanalyze_results_folder)
                        except Exception,e:
                            print_out = "No files found in: "+altanalyze_results_folder
                            IndicatorWindow(print_out,'Continue')

                    if len(exp_file)>0 or ((len(exp_file)>0 or len(gene_symbol)>0) and analysisType == 'Sashimi-Plot'):
                        analysis = 'AltExonViewer'
                        values = species,array_type,exp_file,gene_symbol,show_introns,analysisType
                        try: StatusWindow(values,analysis) ### display an window with download status
                        except Exception: pass
                        #if len(altgenes_file)>0 or ' ' in gene_symbol or ((len(exp_file)>0 or len(gene_symbol)>0) and analysisType == 'Sashimi-Plot'):
                        if len(analysisType)>0:
                            ### Typically have a Tkinter related error
                            rebootAltAnalyzeGUI(selected_parameters,user_variables)
                        else:
                            AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "Either no gene or no AltResults folder selected."
                        IndicatorWindow(print_out,'Continue')

            if additional_analyses == 'Network Visualization':
                selected_parameters.append('Network Visualization')

                supported_interaction_types = getSupportedGeneSetTypes(species,'gene-interactions')                
                supported_geneset_types = getSupportedGeneSetTypes(species,'gene-mapp')
                supported_geneset_types += getSupportedGeneSetTypes(species,'gene-go')
                option_db['GeneSetSelection_network'].setArrayOptions(['None Selected']+supported_geneset_types)
                option_db['PathwaySelection_network'].setArrayOptions(['None Selected'])
                #option_db['PathwaySelection'].setArrayOptions(supported_genesets)
            
                status = 'repeat'
                while status == 'repeat':
                    ### If no databases present download and populate gene-interactions folder
                    if len(supported_interaction_types)==0:
                        print_out = 'No interaction databases available.\nPress Continue to download interaction\ndatabases for this species.' 
                        IndicatorWindow(print_out,'Continue')
                        downloadInteractionDBs(species,'parent')
                        
                    ### Get present interaction databases (including custom added)
                    updated_list=[]
                    if 'WikiPathways' in supported_interaction_types: updated_list.append('WikiPathways')
                    if 'KEGG' in supported_interaction_types: updated_list.append('KEGG')
                    if 'TFTargets' in supported_interaction_types: updated_list.append('TFTargets')
                    if 'BioGRID' in supported_interaction_types: updated_list.append('BioGRID')
                    for db in supported_interaction_types:
                        if 'microRNATargets' in db:
                            updated_list.append('common-microRNATargets'); updated_list.append('all-microRNATargets')
                        elif 'DrugBank' in db:
                            updated_list.append('common-DrugBank'); updated_list.append('all-DrugBank')
                        elif db not in updated_list: updated_list.append(db)
                    option_db['interactionDirs'].setArrayOptions(updated_list)
                
                    root = Tk()
                    root.title('AltAnalyze: Create and Visualize Interaction Networks')
                    gu = GUI(root,option_db,option_list['network'],'')
                    Genes_network = gu.Results()['Genes_network']
                    inputDir_network = gu.Results()['input_ID_file']
                    GeneSetSelection_network = gu.Results()['GeneSetSelection_network']
                    inputType_network = gu.Results()['inputType_network']
                    PathwaySelection_network = gu.Results()['PathwaySelection_network']
                    OntologyID_network = gu.Results()['OntologyID_network']
                    interactionDirs = gu.Results()['interactionDirs']
                    degrees = gu.Results()['degrees']
                    update_interactions = gu.Results()['update_interactions']
                    expressionFile_network = gu.Results()['elite_exp_file']
                    outputDir_network = gu.Results()['output_net_folder']
                    includeExpIDs_network = gu.Results()['includeExpIDs_network']

                    ### Set the below variables to the appropriate object types
                    if update_interactions == 'yes': update_interactions = True
                    else: update_interactions = False
                    if len(inputDir_network) == 0: inputDir_network = None
                    if len(expressionFile_network) == 0: expressionFile_network = None
                    if len(Genes_network) == 0: Genes_network = None
                    if len(outputDir_network) == 0: outputDir_network = None
                    if len(GeneSetSelection_network) == 'None Selected': GeneSetSelection_network = None
                    if includeExpIDs_network=='yes': includeExpIDs_network = True
                    else: includeExpIDs_network = False
                    
                    ### Save these as instances of GeneSelectionParameters (easier this way to add more object types in the future)
                    gsp = GeneSelectionParameters(species,array_type,vendor) ### only species currently neeed
                    gsp.setGeneSet(GeneSetSelection_network)
                    gsp.setPathwaySelect(PathwaySelection_network)
                    gsp.setGeneSelection(Genes_network)
                    gsp.setOntologyID(OntologyID_network)
                    gsp.setIncludeExpIDs(includeExpIDs_network)
                    if update_interactions:
                        downloadInteractionDBs(species,'parent')
                    if outputDir_network==None:
                        print_out = "No output directory selected."
                        IndicatorWindow(print_out,'Continue')
                    elif inputDir_network != None or GeneSetSelection_network != None or Genes_network != None:
                        analysis = 'network'
                        values = inputDir_network,inputType_network,outputDir_network,interactionDirs,degrees,expressionFile_network,gsp
                        StatusWindow(values,analysis,windowType='parent') ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input gene IDs, expression file or GeneSet selected."
                        IndicatorWindow(print_out,'Continue')
                        
            if additional_analyses == 'Hierarchical Clustering':
                selected_parameters.append('Hierarchical Clustering')
                
                supported_geneset_types = getSupportedGeneSetTypes(species,'gene-mapp')
                supported_geneset_types += getSupportedGeneSetTypes(species,'gene-go')
                option_db['GeneSetSelection'].setArrayOptions(['None Selected']+supported_geneset_types)
                option_db['PathwaySelection'].setArrayOptions(['None Selected'])
                option_db['ClusterGOElite'].setArrayOptions(['None Selected','all']+supported_geneset_types)
                #option_db['PathwaySelection'].setArrayOptions(supported_genesets)
                
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Create a Heatmap from an Expression Matrix')
                    gu = GUI(root,option_db,option_list['heatmap'],'')
                    try: input_cluster_file = gu.Results()['input_cluster_file']
                    except Exception: input_cluster_file = ''
                    column_metric = gu.Results()['column_metric']
                    column_method = gu.Results()['column_method']
                    row_metric = gu.Results()['row_metric']
                    row_method = gu.Results()['row_method']
                    color_selection = gu.Results()['color_selection']
                    cluster_rows = gu.Results()['cluster_rows']
                    cluster_columns = gu.Results()['cluster_columns']
                    GeneSetSelection = gu.Results()['GeneSetSelection']
                    PathwaySelection = gu.Results()['PathwaySelection']
                    GeneSelection = gu.Results()['GeneSelection']
                    ClusterGOElite = gu.Results()['ClusterGOElite']
                    HeatmapAdvanced = gu.Results()['HeatmapAdvanced']
                    JustShowTheseIDs = gu.Results()['JustShowTheseIDs']
                    geneSetName = gu.Results()['heatmapGeneSets']
                    try: CorrelationCutoff = float(gu.Results()['CorrelationCutoff'])
                    except Exception: CorrelationCutoff=None
                    OntologyID = gu.Results()['OntologyID']
                    transpose = gu.Results()['transpose']
                    normalization = gu.Results()['normalization']
                    contrast = gu.Results()['contrast']
                    if transpose == 'yes': transpose = True
                    else: transpose = False
                    translate={'None Selected':'','Exclude Cell Cycle Effects':'excludeCellCycle',
                               'Top Correlated Only':'top','Positive Correlations Only':'positive',
                               'Perform Iterative Discovery':'guide', 'Intra-Correlated Only':'IntraCorrelatedOnly',
                               'Correlation Only to Guides':'GuideOnlyCorrelation','Perform Monocle':'monocle'}
                    try:
                        if 'None Selected' in HeatmapAdvanced: pass
                    except Exception: HeatmapAdvanced = ('None Selected')
                    if ('None Selected' in HeatmapAdvanced and len(HeatmapAdvanced)==1) or 'None Selected' == HeatmapAdvanced: pass
                    else:
                        try:
                            GeneSelection += ' '+string.join(list(HeatmapAdvanced),' ')
                            for name in translate:
                                GeneSelection = string.replace(GeneSelection,name,translate[name])
                            GeneSelection = string.replace(GeneSelection,'  ',' ')
                            if 'top' in GeneSelection or 'positive' in GeneSelection or 'IntraCorrelatedOnly' in GeneSelection: #or 'guide' in GeneSelection or 'excludeCellCycle' in GeneSelection - will force correlation to selected genes
                                GeneSelection+=' amplify'
                        except Exception: pass
                    ### This variable isn't needed later, just now to indicate not to correlate in the first round
                    if GeneSetSelection  != 'None Selected' and PathwaySelection == ['None Selected']:
                        PathwaySelection = [gu.Results()[GeneSetSelection][0]] ### Default this to the first selection

                    GeneSetSelection = string.replace(GeneSetSelection,'\n',' ')
                    GeneSetSelection = string.replace(GeneSetSelection,'\r',' ')
                    #print [GeneSetSelection, JustShowTheseIDs, GeneSelection,ClusterGOElite,normalization]
                    if GeneSetSelection != 'None Selected' or GeneSelection != '' or normalization != 'NA' or JustShowTheseIDs != '' or JustShowTheseIDs != 'None Selected':
                        gsp = GeneSelectionParameters(species,array_type,vendor)
                        if CorrelationCutoff!=None: #len(GeneSelection)>0 and 
                            gsp.setRhoCutoff(CorrelationCutoff)
                            GeneSelection = 'amplify '+GeneSelection
                            if 'GuideOnlyCorrelation' in GeneSelection:
                                ### Save the correlation cutoff for ICGS but don't get expanded correlation sets in the first round
                                GeneSelection = string.replace(GeneSelection,'GuideOnlyCorrelation','')
                                GeneSelection = string.replace(GeneSelection,'amplify','')
                                GeneSelection = string.replace(GeneSelection,'  ','')
                        gsp.setGeneSet(GeneSetSelection)
                        gsp.setPathwaySelect(PathwaySelection)
                        gsp.setGeneSelection(GeneSelection)
                        gsp.setOntologyID(OntologyID)
                        gsp.setTranspose(transpose)
                        gsp.setNormalize(normalization)
                        gsp.setJustShowTheseIDs(JustShowTheseIDs)
                        gsp.setClusterGOElite(ClusterGOElite)
                        gsp.setStoreGeneSetName(geneSetName)
                        transpose = gsp ### this allows methods that don't transmit this object to also work
                    if len(input_cluster_file)>0:
                        analysis = 'createHeatMap'
                        color_selection=string.replace(color_selection, '-','_')
                        if cluster_rows == 'no': row_method = None
                        if cluster_columns == 'no': column_method = None
                        values = input_cluster_file, row_method, row_metric, column_method, column_metric, color_selection, transpose, contrast
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')
                        
            if additional_analyses == 'Dimensionality Reduction':
                selected_parameters.append('Dimensionality Reduction')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Perform Dimensionality Reduction from an Expression Matrix')
                    gu = GUI(root,option_db,option_list['PCA'],'')
                    try: input_cluster_file = gu.Results()['input_cluster_file']
                    except Exception: input_cluster_file = ''
                    dimensions = gu.Results()['dimensions']
                    pca_labels = gu.Results()['pca_labels']
                    pca_algorithm = gu.Results()['pca_algorithm']
                    zscore = gu.Results()['zscore']
                    transpose = gu.Results()['transpose']
                    geneSetName = gu.Results()['pcaGeneSets']
                    reimportModelScores = gu.Results()['reimportModelScores']
                    if reimportModelScores == 'yes':
                        reimportModelScores = True
                    else:
                        reimportModelScores = False
                    try:
                        colorByGene = gu.Results()['colorByGene']
                        colorByGene_temp = string.replace(colorByGene,' ','')
                        if len(colorByGene_temp)==0:
                            colorByGene = None
                        else:
                            #Standardize the delimiter
                            colorByGene = string.replace(colorByGene,'|',' ')
                            colorByGene = string.replace(colorByGene,',',' ')
                            colorByGene = string.replace(colorByGene,'\r',' ')
                            colorByGene = string.replace(colorByGene,'\n',' ')
                            colorByGene = string.replace(colorByGene,'  ',' ')
                            if colorByGene[0] == ' ': colorByGene=colorByGene[1:]
                            if colorByGene[-1] == ' ': colorByGene=colorByGene[:-1]
                    except Exception: colorByGene = None
                    if len(geneSetName)==0:
                        geneSetName = None
                    if len(input_cluster_file)>0:
                        analysis = 'performPCA'
                        if transpose == 'yes': transpose = True
                        else: transpose = False
                        values = input_cluster_file, pca_labels, dimensions, pca_algorithm, transpose, geneSetName, species, zscore, colorByGene, reimportModelScores
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')

            if additional_analyses == 'Lineage Analysis' or additional_analyses == 'Cell Classification':
                selected_parameters.append('Lineage Analysis')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    if species == 'Mm':
                        option_db['compendiumPlatform'].setDefaultOption('gene')
                    if species == 'Hs':
                        option_db['compendiumPlatform'].setDefaultOption('exon')
                    if array_type == "3'array":
                        option_db['compendiumType'].setArrayOptions(["protein_coding"])
                    root.title('AltAnalyze: Perform CellHarmony and LineageProfiler Analysis')
                    gu = GUI(root,option_db,option_list['LineageProfiler'],'')
                    input_exp_file = gu.Results()['input_lineage_file']
                    compendiumPlatform = gu.Results()['compendiumPlatform']
                    try: classificationAnalysis = gu.Results()['classificationAnalysis']
                    except: classificationAnalysis = 'cellHarmony'
                    compendiumType = gu.Results()['compendiumType']
                    markerFinder_file = gu.Results()['markerFinder_file']
                    geneModel_file = gu.Results()['geneModel_file']
                    modelDiscovery = gu.Results()['modelDiscovery']
                    pearsonThreshold = gu.Results()['PearsonThreshold']
                    returnCentroids = gu.Results()['returnCentroids']
                    performDiffExp = gu.Results()['performDiffExp']
                    useAdjPval = gu.Results()['UseAdjPval']
                    pvalThreshold = gu.Results()['pvalThreshold']
                    foldCutoff = gu.Results()['FoldCutoff']
                    if '.png' in markerFinder_file or '.pdf' in markerFinder_file:
                        markerFinder_file=markerFinder_file[:-4]+'.txt'
                    if len(geneModel_file) == 0: geneModel_file = None
                    if len(modelDiscovery) == 0: modelDiscovery = None
                    if len(input_exp_file)>0:
                        analysis = 'runLineageProfiler'
                        fl = ExpressionFileLocationData('','','','') ### Create this object to store additional parameters for LineageProfiler
                        fl.setCompendiumType(compendiumType)
                        fl.setCompendiumPlatform(compendiumPlatform)
                        fl.setClassificationAnalysis(classificationAnalysis)
                        fl.setPearsonThreshold(float(pearsonThreshold))
                        fl.setReturnCentroids(returnCentroids)
                        fl.setPeformDiffExpAnalysis(performDiffExp)
                        fl.setUseAdjPvalue(useAdjPval)
                        fl.setPvalThreshold(pvalThreshold)
                        fl.setFoldCutoff(foldCutoff)
                        """
                        print fl.PeformDiffExpAnalysis()
                        print fl.CompendiumType()
                        print fl.CompendiumPlatform()
                        print fl.ClassificationAnalysis()
                        print fl.PearsonThreshold()
                        print fl.ReturnCentroids()
                        print fl.PeformDiffExpAnalysis()
                        print fl.PvalThreshold()
                        print fl.FoldCutoff()
                        print fl.UseAdjPvalue()
                        print fl.PearsonThreshold()"""
                        values = fl, input_exp_file, vendor, markerFinder_file, geneModel_file, modelDiscovery
                        StatusWindow(values,analysis) ### display an window with download status
                        ### Typically have a Tkinter related error
                        try: rebootAltAnalyzeGUI(selected_parameters[:-1],user_variables)
                        except: 
                            AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                        
                        #else:
                        #AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')
                    print 'here'
                    
            if additional_analyses == 'MarkerFinder Analysis':
                selected_parameters.append('MarkerFinder Analysis')
                status = 'repeat'
                while status == 'repeat':
                    root = Tk()
                    root.title('AltAnalyze: Perform MarkerFinder Analysis from the Input in ExpressionInput')
                    gu = GUI(root,option_db,option_list['MarkerFinder'],'')
                    input_exp_file = gu.Results()['input_markerfinder_file']
                    genes_to_output = gu.Results()['compendiumPlatform']
                    if len(geneModel_file) == 0: geneModel_file = None
                    if len(modelDiscovery) == 0: modelDiscovery = None
                    if len(input_exp_file)>0:
                        analysis = 'runLineageProfiler'
                        fl = ExpressionFileLocationData('','','','') ### Create this object to store additional parameters for LineageProfiler
                        fl.setCompendiumType(compendiumType)
                        fl.setCompendiumPlatform(compendiumPlatform)
                        values = fl, input_exp_file, vendor, markerFinder_file, geneModel_file, modelDiscovery
                        StatusWindow(values,analysis) ### display an window with download status
                        AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()
                    else:
                        print_out = "No input expression file selected."
                        IndicatorWindow(print_out,'Continue')

        if 'CEL files' in run_from_scratch or 'RNA-seq reads' in run_from_scratch or 'Feature Extraction' in run_from_scratch or 'Chromium' in run_from_scratch:
            """Designate CEL, Agilent or BED file directory, Dataset Name and Output Directory"""
            assinged = 'no'
            while assinged == 'no': ### Assigned indicates whether or not the CEL directory and CDF files are defined
                if species == 'Rn' or array_type == 'RNASeq': del option_list['InputCELFiles'][-1] ### Don't examine xyb
                #print (((backSelect,selected_parameters)))
                if backSelect == 'no' or 'InputCELFiles' == selected_parameters[-1]:
                    selected_parameters.append('InputCELFiles'); backSelect = 'no'
                    root = Tk()
                    if array_type == 'RNASeq':
                        root.title('AltAnalyze: Select Exon and/or Junction files to analyze'); import_file = 'BED, BAM, TAB or TCGA'
                    elif '10X' in vendor:
                        root.title('AltAnalyze: Select Chromium Sparse Matrix Filtered Matrix'); import_file = 'Filtered Matrix'
                    elif vendor == 'Agilent':
                        root.title('AltAnalyze: Select Agilent Feature Extraction text files to analyze'); import_file = '.txt'
                    else:
                        root.title('AltAnalyze: Select CEL files for APT'); import_file = '.CEL'
                    gu = GUI(root,option_db,option_list['InputCELFiles'],'')
                else: gu = PreviousResults(old_options)
                dataset_name = gu.Results()['dataset_name']
                try: remove_xhyb = gu.Results()['remove_xhyb']
                except KeyError: remove_xhyb = 'no'
                try:
                    multiThreading = gu.Results()['multithreading']
                    if multiThreading == 'yes': multiThreading = True
                    else: multiThreading = False
                except KeyError: multiThreading = True
                try:
                    build_exon_bedfile = gu.Results()['build_exon_bedfile']
                    try: normalize_feature_exp = 'RPKM'
                    except Exception: pass
                except KeyError: build_exon_bedfile = 'no'
                try:
                    input_fastq_dir = gu.Results()['input_fastq_dir']
                except Exception: pass
                try: channel_to_extract = gu.Results()['channel_to_extract']
                except Exception: channel_to_extract = 'no'
                
                if build_exon_bedfile == 'yes' and len(input_fastq_dir)==0:
                    print_out = 'Please note: AltAnalyze will exit immediately after\nimporting your junction results to allow you to build\nyour exon count files and reload this data.' 
                    IndicatorWindowSimple(print_out,'Continue')
                    run_from_scratch = 'buildExonExportFiles'
                if len(dataset_name)<1:
                    print_out = "Please provide a name for the dataset before proceeding."
                    IndicatorWindow(print_out,'Continue')
                elif 'input_cel_dir' in gu.Results() or 'input_fastq_dir' in gu.Results():
                    if len(input_fastq_dir)>0:
                        import RNASeq
                        cel_files = RNASeq.runKallisto(species,'',input_fastq_dir,input_fastq_dir,mlp,returnSampleNames=True)
                        try: output_dir = gu.Results()['output_CEL_dir']
                        except KeyError: output_dir = input_fastq_dir
                        """ ### Change made in version 2.1.2
                        option_db['perform_alt_analysis'].setArrayOptions(['NA'])
                        option_db['exon_exp_threshold'].setArrayOptions(['NA'])
                        option_db['exon_rpkm_threshold'].setArrayOptions(['NA'])
                        option_db['expression_threshold'].setArrayOptions(['NA'])
                        option_db['gene_exp_threshold'].setArrayOptions(['NA'])
                        """
                        assinged = 'yes'
                    else:
                        cel_file_dir = gu.Results()['input_cel_dir']
                        if '10X' in vendor:
                            sparse_matrix_file = gu.Results()['input_cel_dir'] # 'filtered_gene_bc_matrices'
                            def import10XSparseMatrixHeaders(matrix_file):
                                import csv
                                barcodes_path = string.replace(matrix_file,'matrix.mtx','barcodes.tsv' )
                                barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
                                barcodes = map(lambda x: string.replace(x,'-1',''), barcodes)
                                return barcodes
                            barcodes = import10XSparseMatrixHeaders(sparse_matrix_file)
                            cel_files = barcodes
                        else:
                            cel_files,cel_files_fn=identifyCELfiles(cel_file_dir,array_type,vendor)
                        try: output_dir = gu.Results()['output_CEL_dir']
                        except KeyError: output_dir = cel_file_dir
                        if len(output_dir)==0: output_dir = cel_file_dir
                        if len(cel_files)>0: assinged = 'yes' ### CEL files are present in this directory
                        else:
                            print_out = "No valid "+import_file+" files were found in the directory\n"+cel_file_dir+"\nPlease verify and try again."
                            IndicatorWindow(print_out,'Continue')
                    
                else:
                    print_out = "The directory containing "+import_file+" files has not\nbeen assigned! Select a directory before proceeding."
                    IndicatorWindow(print_out,'Continue')

            if array_type != 'RNASeq' and vendor != 'Agilent' and len(input_fastq_dir)==0 and '10X' not in vendor:
                ### Specific to Affymetrix CEL files
                cel_file_list_dir = exportCELFileList(cel_files_fn,cel_file_dir)
                """Determine if Library and Annotations for the array exist, if not, download or prompt for selection"""
                specific_array_types,specific_array_type = identifyArrayType(cel_files_fn); num_array_types = len(specific_array_types)
                #except Exception: pass; num_array_types=1; specific_array_type = None
                importSupportedArrayInfo()
                try:
                    sa = supproted_array_db[specific_array_type]; array_species = sa.Species(); cel_array_type = sa.ArrayType()
                except Exception: library_dir=''; array_species=''; annotation_dir=''; cel_array_type=''
                if backSelect == 'no':
                    ### Check for issues with arrays or user input options
                    if num_array_types>1: ### More than one array type found in the directory
                        print_out = 'Warning!!!!!!!\n\nMultiple array_types found ("'+specific_array_types[0]+'" and "'+specific_array_types[1]+'").\nIt is recommended you restart, otherwise, APT will try\n to process all different array types together as "'+specific_array_types[-1]+'".'
                        IndicatorWindow(print_out,'Continue with Existing')
                    if array_species != species and len(array_species)>0:
                        print_out = "The CEL files indicate that the proper\nspecies is "+array_species+", however, you\nindicated "+species+ ". The species indicated by the CEL\nfiles will be used instead."
                        IndicatorWindow(print_out,'Continue')
                        species = array_species
                        try: spdirs = read_directory('/AltDatabase/'+species)
                        except Exception: spdirs = []
                        if len(spdirs)==0:
                            print_out = 'Valid database directories were not found for this species.\nPlease re-install database.'
                            IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()

                    if cel_array_type != array_type and len(cel_array_type)>0:
                        print_out = "The CEL files indicate that the proper\narray type is "+cel_array_type+", however, you\nindicated "+array_type+ "." #The array type indicated by the CEL\nfiles will be used instead
                        #IndicatorWindow(print_out,'Continue')
                        fw = FeedbackWindow(print_out,'Use AltAnalyze Recommended',"Use Original Selected")
                        choice = fw.ButtonSelection()['button']
                        if choice == 'Use AltAnalyze Recommended': 
                            array_type = cel_array_type
                            option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type
                            option_db['array_type'].setArrayOptions(array_list)
                            #user_variables['array_type'] = array_type
                            ### See if the library and annotation files are on the server or are local
                        else: specific_array_type = ''; annotation_dir=''

                if specific_array_type == None:
                    if array_type == 'exon':
                        if species == 'Hs': specific_array_type = 'HuEx-1_0-st-v2'
                        if species == 'Mm': specific_array_type = 'MoEx-1_0-st-v2'
                        if species == 'Rn': specific_array_type = 'RaEx-1_0-st-v2'
                    elif array_type == 'gene':
                        if species == 'Hs': specific_array_type = 'HuGene-1_0-st-v1'
                        if species == 'Mm': specific_array_type = 'MoGene-1_0-st-v1'
                        if species == 'Rn': specific_array_type = 'RaGene-1_0-st-v1'
                    elif array_type == 'AltMouse': specific_array_type = 'altMouseA'
                    """ ### Comment this out to allow for different junction sub-types (likely do the above in the future)
                    elif array_type == 'junction':
                        if species == 'Hs': specific_array_type = 'HJAY_v2'
                        if species == 'Mm': specific_array_type = 'MJAY_v2'
                    """

                if specific_array_type in supproted_array_db:
                    input_cdf_file, annotation_dir, bgp_file, clf_file = getAffyFiles(specific_array_type,species)
                else: input_cdf_file=''; bgp_file = ''; clf_file = ''
                ### Remove the variable names for Library and Annotation file selection if these files are found
                option_list_library=[]
                if len(input_cdf_file)>0:
                    for i in option_list['InputLibraryFiles']:
                        if i != 'input_cdf_file': option_list_library.append(i)
                if len(annotation_dir)>0:
                    for i in option_list['InputLibraryFiles']:
                        if i != 'input_annotation_file': option_list_library.append(i)
                if len(option_list_library)==0:
                    option_list_library = option_list['InputLibraryFiles']

                """Identify and copy over any Libary or Annotation files on the computer"""                    
                if (len(input_cdf_file)==0 and len(annotation_dir) == 0) and backSelect == 'no':
                    ### Note: above line used to be "or" between the input_cdf_file and annotation_dir
                    ### this was discontinued in version 2.0.9 since the annotation file is no longer needed
                    ### unless the array type is not in the GO-elite database
                    assinged = 'no'
                    while assinged == 'no': ### Assigned indicates whether or not the CEL directory and CDF files are defined
                        if array_type == "3'array":
                            op = option_db['input_cdf_file']; input_cdf_file_label = op.Display()
                            op.setNotes('   note: the CDF file is apart of the standard library files for this array.   ')
                            input_cdf_file_label = string.replace(input_cdf_file_label,'PGF','CDF')
                            op.setDisplay(input_cdf_file_label)
                        if array_type == 'exon':
                            op = option_db['input_annotation_file']
                            new_notes = string.replace(op.Notes(),'this array','the Gene 1.0 array (NOT Exon)')
                            new_notes = string.replace(new_notes,'annotations','transcript cluster annotations')
                            new_display = string.replace(op.Display(),'your array','the Gene 1.0 array')
                            op.setDisplay(new_display)
                            op.setNotes(new_notes)
                        #if backSelect == 'no' or 'Library' == selected_parameters[-1]:
                        selected_parameters.append('Library')#; backSelect = 'no'
                        root = Tk(); root.title('AltAnalyze: Select Affymetrix Library and Annotation files')
                        gu = GUI(root,option_db,option_list_library,'')
                        #else: gu = PreviousResults(old_options)                    
                        if 'input_cdf_file' in option_list_library: ### Deals with Annotation Files
                            if 'input_cdf_file' in gu.Results():
                                input_cdf_file = gu.Results()['input_cdf_file']; input_cdf_file_lower = string.lower(input_cdf_file)
                                if array_type == "3'array":
                                    if '.cdf' in input_cdf_file_lower:
                                        clf_file='';bgp_file=''; assinged = 'yes'
                                        ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                        icf_list = string.split(input_cdf_file,'/'); cdf_short = icf_list[-1]
                                        destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                                        destination_parent = osfilepath(destination_parent+cdf_short)
                                        #print destination_parent
                                        #print input_cdf_file
                                        if destination_parent not in input_cdf_file:
                                            info_list = input_cdf_file,destination_parent; StatusWindow(info_list,'copy')
                                    else:
                                        print_out = "The file;\n"+input_cdf_file+"\ndoes not appear to be a valid Affymetix\nlibrary file. If you do not have library files, you must\ngo to the Affymetrix website to download."
                                        IndicatorWindow(print_out,'Continue')            
                                else:
                                    if '.pgf' in input_cdf_file_lower:
                                        ###Check to see if the clf and bgp files are present in this directory 
                                        icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
                                        clf_short = string.replace(cdf_short,'.pgf','.clf')
                                        kil_short = string.replace(cdf_short,'.pgf','.kil') ### Only applies to the Glue array
                                        if array_type == 'exon' or array_type == 'junction': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
                                        else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
                                        dir_list = read_directory(parent_dir)
                                        if clf_short in dir_list and bgp_short in dir_list:
                                            pgf_file = input_cdf_file
                                            clf_file = string.replace(pgf_file,'.pgf','.clf')
                                            kil_file = string.replace(pgf_file,'.pgf','.kil') ### Only applies to the Glue array
                                            if array_type == 'exon' or array_type == 'junction': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                                            else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                                            assinged = 'yes'
                                            ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                            destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                                            #print destination_parent
                                            #print input_cdf_file
                                            if destination_parent not in input_cdf_file:
                                                info_list = input_cdf_file,osfilepath(destination_parent+cdf_short); StatusWindow(info_list,'copy')
                                                info_list = clf_file,osfilepath(destination_parent+clf_short); StatusWindow(info_list,'copy')
                                                info_list = bgp_file,osfilepath(destination_parent+bgp_short); StatusWindow(info_list,'copy')
                                                if 'Glue' in pgf_file:
                                                    info_list = kil_file,osfilepath(destination_parent+kil_short); StatusWindow(info_list,'copy')
                                        else:
                                            print_out = "The directory;\n"+parent_dir+"\ndoes not contain either a .clf or antigenomic.bgp\nfile, required for probeset summarization."
                                            IndicatorWindow(print_out,'Continue')                                   
                                    else:
                                        print_out = "The file;\n"+input_cdf_file+"\ndoes not appear to be a valid Affymetix\nlibrary file. If you do not have library files, you must\ngo to the Affymetrix website to download."
                                        IndicatorWindow(print_out,'Continue')
                            else: 
                                print_out = "No library file has been assigned. Please\nselect a valid library file for this array."
                                IndicatorWindow(print_out,'Continue')                                
                        if 'input_annotation_file' in option_list_library: ### Deals with Annotation Files
                            assinged = 'yes'
                            if 'input_annotation_file' in gu.Results():
                                input_annotation_file = gu.Results()['input_annotation_file']; input_annotation_lower = string.lower(input_annotation_file)
                                if '.csv' in input_annotation_lower:
                                    ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                    icf_list = string.split(input_annotation_file,'/'); csv_short = icf_list[-1]
                                    destination_parent = 'AltDatabase/affymetrix/'+species+'/'
                                    #print destination_parent
                                    #print input_cdf_file
                                    if destination_parent not in input_cdf_file:
                                        info_list = input_annotation_file,filepath(destination_parent+csv_short); StatusWindow(info_list,'copy')
                                    sd = SupprotedArrays(specific_array_type,cdf_short,csv_short,species,array_type)
                                    supproted_array_db[specific_array_type] = sd
                                    try: exportSupportedArrayInfo()
                                    except Exception:
                                        print 'Cannot write Config/ArrayFileInfo.txt to the Config directory (like Permissions Error)'
                                        continue ### Occurs if the file is open... not critical to worry about       

        if run_from_scratch == 'Process Expression file':
            status = 'repeat'
            while status == 'repeat':
                if backSelect == 'no' or 'InputExpFiles' == selected_parameters[-1]:
                    root = Tk(); root.title('AltAnalyze: Select Expression File for Filtering')
                    selected_parameters.append('InputExpFiles'); backSelect = 'no'
                    gu = GUI(root,option_db,option_list['InputExpFiles'],'')
                else: gu = PreviousResults(old_options)
                try: input_exp_file = gu.Results()['input_exp_file']
                except KeyError: input_exp_file = '' ### Leave this blank so that the default directory is used
                try: input_stats_file = gu.Results()['input_stats_file']
                except KeyError: input_stats_file = '' ### Leave this blank so that the default directory is used
                #if array_type == 'exon':
                if 'steady-state' in input_exp_file or 'steady-state' in input_stats_file:
                    print_out = "Do not select steady-state expression files.."
                    IndicatorWindow(print_out,'Continue'); output_dir=''
                elif len(input_exp_file)>0:
                    try: output_dir = gu.Results()['output_dir']
                    except KeyError: output_dir = '' ### Leave this blank so that the default directory is used
                    try: cel_files, array_linker_db = ExpressionBuilder.getArrayHeaders(input_exp_file)
                    except Exception:
                        print_out = "Input Expression file does not have a valid format."
                        IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()  
                    if len(cel_files)>0: status = 'continue'
                    else:
                        if '.mtx' in input_exp_file:
                            array_type = '10XGenomics'
                            array_full == '10X Genomics sparse matrix'
                            vendor = '10x'
                            print_out = "The expression file:\n"+input_exp_file+"\nis a 10x Genomics matrix... change the Platform to 10x Genomics Aligned in the main menu."
                            IndicatorWindow(print_out,'Continue')
                        else:
                            print_out = "The expression file:\n"+input_exp_file+"\ndoes not appear to be a valid expression file. Check to see that\nthis is the correct tab-delimited text file."
                            IndicatorWindow(print_out,'Continue')
                else:
                    print_out = "No input expression file selected."
                    IndicatorWindow(print_out,'Continue')
                    
                if len(output_dir)<1:
                    ### Set output to the same directory or parent if none selected
                    if 'ExpressionInput' in input_exp_file: i = -2
                    else: i = -1
                    output_dir = string.join(string.split(input_exp_file,'/')[:i],'/')
                
            try: prior_platform = user_variables['prior_platform']
            except Exception: prior_platform = None
            if array_type == 'RNASeq' or prior_platform == 'RNASeq':
                steady_state = string.replace(input_exp_file,'.txt','-steady-state.txt')
                count = verifyFileLength(steady_state)
                if count == 0 or 'exp.' not in input_exp_file: #No counts file
                    systm = getGeneSystem(input_exp_file)
                    ### Wrong platform listed
                    array_type = "3'array"
                    prior_platform = 'RNASeq'
                    vendor = 'other:'+systm ### Ensembl linked system name
                    user_variables['manufacturer_selection'] = vendor
                    user_variables['prior_platform'] = prior_platform
                    if old_options==[] or 'marker_finder' not in old_options: ### If we haven't hit the back button
                        option_list,option_db = importUserOptions(array_type) ### will re-set the paramater values, so not good for back select
                        user_variables['array_type'] = array_type
                        
            if array_type == "3'array":
                ### This is the new option for expression filtering of non-RNASeq classified data
                try:
                    #print option_db['rpkm_threshold'].DefaultOption(),1
                    if 'rpkm_threshold' in option_db:
                        option_db['rpkm_threshold'].setArrayOptions('1')
                        if "other:Symbol" in vendor or "other:Ensembl" in vendor:
                            option_db['rpkm_threshold'].setDefaultOption('1')
                        if option_db['rpkm_threshold'].DefaultOption() == ['NA']:
                            option_db['rpkm_threshold'].setDefaultOption('1')
                        option_db['rpkm_threshold'].setDisplay('Remove genes expressed below (non-log)')
                    else:
                        option_db['rpkm_threshold'].setArrayOptions('0')
                        option_db['rpkm_threshold'].setDefaultOption('0')
                        option_db['rpkm_threshold'].setDisplay('Remove genes expressed below (non-log)')                     
                except Exception:
                    option_db['rpkm_threshold'].setArrayOptions('0')
                    option_db['rpkm_threshold'].setDefaultOption('0')
                    option_db['rpkm_threshold'].setDisplay('Remove genes expressed below (non-log)')
                    
            if "ExpressionInput" not in output_dir and len(input_exp_file)>1 and "ExpressionInput" not in input_exp_file:
                try:
                    ### If the user designates an output directory that doesn't contain ExpressionInput, move the exp-file there and rename
                    output_dir = output_dir + '/ExpressionInput' ### Store the result files here so that files don't get mixed up
                    try: os.mkdir(output_dir) ### Since this directory doesn't exist we have to make it
                    except OSError: null = [] ### Directory already exists
                    if 'exp.' not in input_exp_file: exp_prefix = 'exp.'
                    else: exp_prefix=''
                    moved_exp_dir = output_dir+'/'+exp_prefix+export.findFilename(input_exp_file)
                    alt_exp_dir = export.findParentDir(input_exp_file)+'/'+exp_prefix+export.findFilename(input_exp_file)
                    export.copyFile(input_exp_file, moved_exp_dir)
                    ### Do the same thing for a groups file
                    try: export.copyFile(string.replace(alt_exp_dir,'exp.','groups.'), string.replace(moved_exp_dir,'exp.','groups.'))
                    except: pass
                    ### Do the same thing for a comps file
                    try: export.copyFile(string.replace(alt_exp_dir,'exp.','comps.'), string.replace(moved_exp_dir,'exp.','comps.'))
                    except: pass
                    input_exp_file = moved_exp_dir
                    if len(input_stats_file)>1: ### Do the same for a stats file
                        if 'stats.' not in input_exp_file: stats_prefix = 'stats.'
                        else: stats_prefix=''
                        moved_stats_dir = output_dir+'/'+stats_prefix+export.findFilename(input_stats_file)
                        export.copyFile(input_stats_file, moved_stats_dir)
                        input_stats_file = moved_stats_dir
                except Exception: None

        if run_from_scratch != 'buildExonExportFiles': ### Update DBs is an option which has been removed from 1.1. Should be a separate menu item soon.
            expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = importDefaults(array_type,species)
            #print vendor
            if '10X' in vendor:
                option_db['rpkm_threshold'].setDefaultOption('1')
            if vendor == 'Affymetrix' or vendor == 'RNASeq':
                option_db['normalize_gene_data'].setArrayOptions(['NA']) ### Only use this option when processing Feature Extraction files or non-Affy non-RNA-Seq data
            if vendor == 'Agilent' and 'Feature Extraction' in run_from_scratch:
                option_db['normalize_gene_data'].setDefaultOption('quantile')
                option_db['normalize_gene_data'].setArrayOptions(['quantile']) ### Only set this as a default when performing Feature Extraction for Agilent data
            if run_from_scratch != 'Process AltAnalyze filtered' and run_from_scratch != 'Annotate External Results':
                proceed = 'no'
                option_db = check_moderated_support(option_db)
                while proceed == 'no':
                    if backSelect == 'no' or 'GeneExpression' == selected_parameters[-1]:
                        selected_parameters.append('GeneExpression'); backSelect = 'no'
                        root = Tk(); root.title('AltAnalyze: Expression Analysis Parameters')
                        gu = GUI(root,option_db,option_list['GeneExpression'],expr_defaults)
                    else: gu = PreviousResults(old_options)
                    try: rpkm_threshold = float(gu.Results()['rpkm_threshold'])
                    except Exception:
                        if array_type == 'RNASeq': rpkm_threshold = 1
                        else: rpkm_threshold = 'NA'

                    if array_type != "3'array":          
                        try: dabg_p = gu.Results()['dabg_p']
                        except Exception:
                            if array_type == 'RNASeq': dabg_p = 1
                            else: dabg_p = 'NA'
                        try: gene_exp_threshold = gu.Results()['gene_exp_threshold']
                        except Exception:
                            if array_type == 'RNASeq': gene_exp_threshold = 1
                            else: gene_exp_threshold = 'NA'
                        try: exon_rpkm_threshold = gu.Results()['exon_rpkm_threshold']
                        except Exception:
                            if array_type == 'RNASeq': exon_rpkm_threshold = 1
                            else: exon_rpkm_threshold = 'NA'
                        try: exon_exp_threshold = gu.Results()['exon_exp_threshold']
                        except Exception:
                            if array_type == 'RNASeq': exon_exp_threshold = 1
                            else: exon_exp_threshold = 'NA'
                        run_from_scratch = gu.Results()['run_from_scratch']
                        try: expression_threshold = gu.Results()['expression_threshold']
                        except Exception:
                            if array_type == 'RNASeq': expression_threshold = 0
                            else: expression_threshold = 'NA'
                        try: perform_alt_analysis = gu.Results()['perform_alt_analysis']
                        except Exception: perform_alt_analysis = 'just expression'
                        try: analyze_as_groups = gu.Results()['analyze_as_groups']
                        except Exception: analyze_as_groups = ''
                        if perform_alt_analysis == 'just expression': perform_alt_analysis = 'expression'
                        else: perform_alt_analysis = 'both'
                        try: avg_all_for_ss = gu.Results()['avg_all_for_ss']
                        except Exception: avg_all_for_ss = 'no'
                        excludeNonExpExons = True
                        if 'all exon aligning' in avg_all_for_ss or 'known' in avg_all_for_ss or 'expressed exons' in avg_all_for_ss:
                            if 'known exons' in avg_all_for_ss and array_type == 'RNASeq': excludeNonExpExons = False
                            if 'known junctions' in avg_all_for_ss and array_type == 'RNASeq':
                                fl.setUseJunctionsForGeneExpression(True)
                                excludeNonExpExons = False
                            avg_all_for_ss = 'yes'
                        else: avg_all_for_ss = 'no'
                    expression_data_format = gu.Results()['expression_data_format']
                    try: normalize_feature_exp = gu.Results()['normalize_feature_exp']
                    except Exception: normalize_feature_exp = 'NA'
                    try: normalize_gene_data = gu.Results()['normalize_gene_data']
                    except Exception: normalize_gene_data = 'NA'
                    include_raw_data = gu.Results()['include_raw_data']
                    run_goelite = gu.Results()['run_goelite']
                    visualize_results = gu.Results()['visualize_results']
                    run_lineage_profiler = gu.Results()['run_lineage_profiler']
                    probability_algorithm = gu.Results()['probability_algorithm']
                    try: FDR_statistic = gu.Results()['FDR_statistic']
                    except Exception: pass
                    try: batch_effects = gu.Results()['batch_effects']
                    except Exception: batch_effects = 'NA'
                    try: marker_finder = gu.Results()['marker_finder']
                    except Exception: marker_finder = 'NA'
                    if 'immediately' in run_goelite: run_goelite = 'yes'
                    else: run_goelite = 'no'
                    passed = 'yes'; print_out = 'Invalid threshold entered for '
                    if array_type != "3'array" and array_type !='RNASeq':
                        try:
                            dabg_p = float(dabg_p)
                            if dabg_p<=0 or dabg_p>1: passed = 'no'; print_out+= 'DABG p-value cutoff '
                        except Exception: passed = 'no'; print_out+= 'DABG p-value cutoff '
                    if array_type != "3'array":   
                        try:
                            try: rpkm_threshold = float(rpkm_threshold)
                            except Exception:
                                expression_threshold = float(expression_threshold)
                                if expression_threshold<1: passed = 'no'; print_out+= 'expression threshold '
                        except Exception: passed = 'no'; print_out+= 'expression threshold '   
                    if array_type == 'RNASeq':
                        try:
                            rpkm_threshold = float(rpkm_threshold)
                            if rpkm_threshold<0: passed = 'no'; print_out+= 'RPKM threshold '
                        except Exception: passed = 'no'; print_out+= 'RPKM threshold '
                        try:
                            exon_exp_threshold = float(exon_exp_threshold)
                            if exon_exp_threshold<0: passed = 'no'; print_out+= 'Exon expression threshold '
                        except Exception: passed = 'no'; print_out+= 'Exon expression threshold '
                        try:
                            exon_rpkm_threshold = float(exon_rpkm_threshold)
                            if exon_rpkm_threshold<0: passed = 'no'; print_out+= 'Exon RPKM threshold '
                        except Exception: passed = 'no'; print_out+= 'Exon RPKM threshold '
                        try:
                            gene_exp_threshold = float(gene_exp_threshold)
                            if gene_exp_threshold<0: passed = 'no'; print_out+= 'Gene expression threshold '
                        except Exception: passed = 'no'; print_out+= 'Gene expression threshold '
                    if visualize_results == 'yes':
                        try:
                            ### Tests to make sure these are installed - required for visualization
                            import matplotlib
                            from numpy import array
                            from scipy import rand
                        except Exception:
                            passed = 'no'; print_out = 'Support for matplotlib, numpy and scipy must specifically be installed to perform data visualization.\n'
                            print_out += traceback.format_exc() ### useful for seeing a warning window with the actuall error
                    if passed == 'no': IndicatorWindow(print_out,'Continue')
                    else: proceed = 'yes'
                if run_lineage_profiler == 'yes':
                    verifyLineageProfilerDatabases(species,'GUI')
            if (perform_alt_analysis == 'both') or (run_from_scratch == 'Process AltAnalyze filtered') or (run_from_scratch == 'Annotate External Results'):
                perform_alt_analysis = 'yes'

                if run_from_scratch == 'Process AltAnalyze filtered':
                    input_filtered_dir = ''
                    while len(input_filtered_dir)<1:
                        if backSelect == 'no' or 'InputFilteredFiles' == selected_parameters[-1]:
                            selected_parameters.append('InputFilteredFiles'); backSelect = 'no'
                            root = Tk(); root.title('AltAnalyze: Select AltAnalyze Filtered Probe set Files')
                            gu = GUI(root,option_db,option_list['InputFilteredFiles'],'')
                        else: gu = PreviousResults(old_options)
                        try:
                            input_filtered_dir = gu.Results()['input_filtered_dir']
                            if 'FullDataset' in input_filtered_dir: alt_exon_defaults[3] = 'all groups'
                        except Exception: input_filtered_dir = ''
                        if input_filtered_dir == '': 
                            print_out = "The directory containing filtered probe set text files has not\nbeen assigned! Select a valid directory before proceeding."
                            IndicatorWindow(print_out,'Continue')
                    fl = ExpressionFileLocationData('','','',''); dataset_name = 'filtered-exp_dir'
                    dirs = string.split(input_filtered_dir,'AltExpression'); parent_dir = dirs[0]
                    exp_file_location_db={}; exp_file_location_db[dataset_name]=fl

                if run_from_scratch == 'Annotate External Results':
                    input_filtered_dir = ''
                    while len(input_filtered_dir)<1:
                        if backSelect == 'no' or 'InputExternalFiles' == selected_parameters[-1]:
                            selected_parameters.append('InputExternalFiles'); backSelect = 'no'
                            root = Tk(); root.title('AltAnalyze: Select AltAnalyze Filtered Probe set Files')
                            gu = GUI(root,option_db,option_list['InputExternalFiles'],'')
                        else: gu = PreviousResults(old_options)
                        try: input_filtered_dir = gu.Results()['input_external_dir']
                        except Exception: input_filtered_dir = ''
                        if input_filtered_dir == '': 
                            print_out = "The directory containing external probe set text files has not\nbeen assigned! Select a valid directory before proceeding."
                            IndicatorWindow(print_out,'Continue')
                    fl = ExpressionFileLocationData('','','',''); dataset_name = 'external-results_dir'
                    dirs = string.split(input_filtered_dir,'AltExpression'); parent_dir = dirs[0]
                    exp_file_location_db={}; exp_file_location_db[dataset_name]=fl

                #print option_list[i:i+len(alt_exon_defaults)+len(functional_analysis_defaults)], alt_exon_defaults+functional_analysis_defaults;kill
                option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type

                if backSelect == 'yes':
                    for option in old_options: ### Set options to user selected
                        try: option_db[option].setDefaultOption(old_options[option])
                        except Exception: pass
                    
                if run_from_scratch == 'Process AltAnalyze filtered':
                    if array_type == 'RNASeq': cs_name = 'known exons'
                    else: cs_name = 'constitutive probesets'
                    functional_analysis_defaults.append(cs_name); option_list['AltAnalyze'].append('avg_all_for_ss')
                    if run_goelite == 'no': ### run_goelite will be set to no by default
                        functional_analysis_defaults.append('unpaired t-test'); option_list['AltAnalyze'].append('probability_algorithm')
                        functional_analysis_defaults.append('decide later'); option_list['AltAnalyze'].append('run_goelite')
                        
                if run_from_scratch == 'Annotate External Results':
                    ### Remove options relating to expression analysis when importing filtered probeset lists
                    options_to_exclude = ['analysis_method','p_threshold','gene_expression_cutoff','alt_exon_fold_cutoff','run_MiDAS']
                    options_to_exclude+= ['export_splice_index_values','probability_algorithm','run_goelite','analyze_all_conditions','calculate_splicing_index_p']
                    for option in options_to_exclude: del option_db[option]
                    
                proceed = 'no'
                while proceed == 'no':
                    if backSelect == 'no' or 'AltAnalyze' == selected_parameters[-1]:
                        selected_parameters.append('AltAnalyze'); backSelect = 'no'; proceed = 'no'

                        root = Tk(); root.title('AltAnalyze: Alternative Exon Analysis Parameters')
                        gu = GUI(root,option_db,option_list['AltAnalyze'],alt_exon_defaults+functional_analysis_defaults); #user_variables = {};
                        try: analyze_all_conditions = gu.Results()['analyze_all_conditions']
                        except KeyError: analyze_all_conditions = 'pairwise'
                        if analyze_all_conditions != 'pairwise':
                            print_out = 'Please note: When AltAnalyze compares all groups, the\nalternative exon fold to be filtered will be based on the\nlargest alternative exon fold for all possible comparisons.' 
                            IndicatorWindowSimple(print_out,'Continue')
                    else: gu = PreviousResults(old_options)
                    try: analysis_method = gu.Results()['analysis_method']
                    except Exception: analysis_method = analysis_method
                    try: p_threshold = gu.Results()['p_threshold']
                    except Exception: p_threshold = 0.05
                    try: gene_expression_cutoff = gu.Results()['gene_expression_cutoff']
                    except Exception: gene_expression_cutoff = 3
                    try: remove_intronic_junctions = gu.Results()['remove_intronic_junctions']
                    except Exception: remove_intronic_junctions = 'NA'
                    try: filter_probeset_types = gu.Results()['filter_probe_types']
                    except Exception: filter_probeset_types = 'core'
                    try: alt_exon_fold_cutoff = gu.Results()['alt_exon_fold_cutoff']
                    except KeyError: alt_exon_fold_cutoff = 2
                    try: permute_p_threshold = gu.Results()['permute_p_threshold']
                    except KeyError: permute_p_threshold = 0.05 ### Doesn't matter, not used
                    try:
                        additional_algorithms = gu.Results()['additional_algorithms']
                        additional_algorithms = AdditionalAlgorithms(additional_algorithms)
                    except KeyError: additionalAlgorithms = AdditionalAlgorithms('')
                    try:
                        additional_score = gu.Results()['additional_score']
                        additional_algorithms.setScore(additional_score)
                    except Exception:
                        try: additional_algorithms.setScore(2)
                        except Exception: pass
                    try: perform_permutation_analysis = gu.Results()['perform_permutation_analysis']
                    except KeyError: perform_permutation_analysis = perform_permutation_analysis
                    try: export_splice_index_values = gu.Results()['export_splice_index_values']
                    except KeyError: export_splice_index_values = export_splice_index_values
                    try: run_MiDAS = gu.Results()['run_MiDAS']
                    except KeyError: run_MiDAS = run_MiDAS
                    try: analyze_all_conditions = gu.Results()['analyze_all_conditions']
                    except KeyError: analyze_all_conditions = analyze_all_conditions
                    try: run_goelite = gu.Results()['run_goelite']
                    except KeyError: run_goelite = run_goelite
                    try: probability_algorithm = gu.Results()['probability_algorithm']
                    except KeyError: probability_algorithm = probability_algorithm
                    try:
                        avg_all_for_ss = gu.Results()['avg_all_for_ss']
                        if 'all exon aligning' in avg_all_for_ss or 'known' in avg_all_for_ss or 'core' in avg_all_for_ss or 'expressed exons' in avg_all_for_ss:
                            avg_all_for_ss = 'yes'
                        else: avg_all_for_ss = 'no'
                    except Exception:
                        try: avg_all_for_ss = avg_all_for_ss
                        except Exception: avg_all_for_ss = 'no'
                    if 'immediately' in run_goelite: run_goelite = 'yes'
                    else: run_goelite = 'no'
                    try: calculate_splicing_index_p = gu.Results()['calculate_splicing_index_p']
                    except KeyError: calculate_splicing_index_p = calculate_splicing_index_p
                    analyze_functional_attributes = gu.Results()['analyze_functional_attributes']
                    filter_for_AS = gu.Results()['filter_for_AS']
                    microRNA_prediction_method = gu.Results()['microRNA_prediction_method']
                    if analysis_method == 'splicing-index': p_threshold = float(p_threshold)
                    else:
                        try: p_threshold = float(permute_p_threshold)
                        except ValueError: permute_p_threshold = permute_p_threshold
                    if analysis_method == 'linearregres-rlm':
                        ### Test installation of rpy and/or R
                        x = [5.05, 6.75, 3.21, 2.66]; y = [1.65, 26.5, -5.93, 7.96]
                        try: s = statistics.LinearRegression(x,y,'no')
                        except Exception:
                            print_out = "The local installation of R and rpy is missing or\nis not properly configured. See the AltAnalyze ReadMe\nfor more information (may require loading AltAnalyze from source code)."
                            IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
                    passed = 'yes'; print_out = 'Invalid threshold entered for '
                    try: gene_expression_cutoff = float(gene_expression_cutoff)
                    except Exception: passed = 'no'; print_out+= 'gene expression cutoff'
                    try: alt_exon_fold_cutoff = float(alt_exon_fold_cutoff)
                    except Exception: passed = 'no'; print_out+= 'alternative exon fold change'
                    try: p_threshold = float(p_threshold)
                    except Exception: passed = 'no'; print_out+= 'alternative exon p-value'
                    if gene_expression_cutoff <= 1: passed = 'no'; print_out+= 'gene expression cutoff'
                    elif alt_exon_fold_cutoff < 1:
                        if analysis_method == 'splicing-index': passed = 'no'; print_out+= 'splicing-index fold change'
                        elif alt_exon_fold_cutoff < 0: passed = 'no'; print_out+= 'alternative exon fold change'
                    elif p_threshold <= 0: passed = 'no'; print_out+= 'alternative exon p-value'
                    if passed == 'no': IndicatorWindow(print_out,'Continue')
                    else: proceed = 'yes'
            if run_goelite == 'yes':
                option_db['get_additional'].setArrayOptions(['---']+importResourceList())
                option_db['get_additional'].setDefaultOption('---')
            
                ### Populate variables based on the existing imported data
                default_resources = option_db['resources_to_analyze'].ArrayOptions() ### Include alternative ontologies and gene-lists
                import_dir1 = '/AltDatabase/goelite/'+species+'/gene-mapp'
                import_dir2 = '/AltDatabase/goelite/'+species+'/gene-go'
                try:
                    gene_mapp_list = read_directory(import_dir1)
                    gene_mapp_list.sort()
                    for file in gene_mapp_list:
                        resource = string.split(file,'-')[-1][:-4]
                        if resource != 'MAPP' and resource not in default_resources and '.txt' in file:
                            default_resources.append(resource)
                except Exception: pass
                try:
                    gene_go_list = read_directory(import_dir2)
                    gene_go_list.sort()
                    for file in gene_go_list:
                        resource = string.split(file,'-')[-1][:-4]
                        if resource != 'GeneOntology' and resource not in default_resources and 'version' not in resource and '.txt' in file:
                            default_resources.append(resource)
                except Exception: pass
                option_db['resources_to_analyze'].setArrayOptions(default_resources)
        
                if run_from_scratch == 'Process AltAnalyze filtered':
                    ### Do not include gene expression analysis filters
                    option_list['GOElite'] = option_list['GOElite'][3:]; goelite_defaults = goelite_defaults[3:]
                if backSelect == 'no' or 'GOElite' == selected_parameters[-1]:
                    selected_parameters.append('GOElite'); backSelect = 'no'
                    root = Tk(); root.title('AltAnalyze: Pathway Analysis Parameters')
                    gu = GUI(root,option_db,option_list['GOElite'],goelite_defaults)
                else: gu = PreviousResults(old_options)
                if run_from_scratch != 'Process AltAnalyze filtered':
                    ge_fold_cutoffs = gu.Results()['ge_fold_cutoffs']
                    ge_pvalue_cutoffs = gu.Results()['ge_pvalue_cutoffs']
                    ge_ptype = gu.Results()['ge_ptype']
                filter_method = gu.Results()['filter_method']
                z_threshold = gu.Results()['z_threshold']
                returnPathways = gu.Results()['returnPathways']
                p_val_threshold = gu.Results()['p_val_threshold']
                change_threshold = gu.Results()['change_threshold']
                resources_to_analyze = gu.Results()['resources_to_analyze']
                pathway_permutations = gu.Results()['pathway_permutations']
                get_additional = gu.Results()['get_additional']
                ORA_algorithm = gu.Results()['ORA_algorithm']
                mod = gu.Results()['mod']
                ge_fold_cutoffs = float(ge_fold_cutoffs)
                change_threshold = float(change_threshold) - 1 ### This reflects the > statement in the GO-Elite filtering
                if ORA_algorithm == 'Fisher Exact Test':
                    pathway_permutations = 'FisherExactTest'
                if get_additional != '---':
                    analysis = 'getAdditionalOnlineResources'
                    values = species,get_additional
                    StatusWindow(values,analysis) ### display an window with download status
    except OSError:
        pass; sys.exit()
    """In this next section, create a set of GUI windows NOT defined by the options.txt file.
    These are the groups and comps files"""
    original_comp_group_list=[]; array_group_list=[]; group_name_list=[]
    if run_from_scratch != 'Process AltAnalyze filtered' and run_from_scratch != 'Annotate External Results': ### Groups and Comps already defined

        if run_from_scratch == 'Process CEL files' or run_from_scratch == 'Process RNA-seq reads' or 'Feature Extraction' in run_from_scratch or 'Process Chromium Matrix' in run_from_scratch:
            if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name+'.txt'
            
            groups_name = string.replace(dataset_name,'exp.','groups.')
            comps_name = string.replace(dataset_name,'exp.','comps.')
            batch_name = string.replace(groups_name,'groups.','batch.') ### may not apply
            if "ExpressionInput" not in output_dir:
                output_dir = output_dir + '/ExpressionInput' ### Store the result files here so that files don't get mixed up
                try: os.mkdir(output_dir) ### Since this directory doesn't exist we have to make it
                except OSError: null = [] ### Directory already exists
            exp_file_dir = output_dir+'/'+dataset_name
            ### store file locations (also use these later when running APT)
            stats_file_dir = string.replace(exp_file_dir,'exp.','stats.')
            groups_file_dir = string.replace(exp_file_dir,'exp.','groups.')
            comps_file_dir = string.replace(exp_file_dir,'exp.','comps.')
            batch_file_dir = string.replace(groups_file_dir, 'groups.','batch.')
            fl = ExpressionFileLocationData(exp_file_dir,stats_file_dir,groups_file_dir,comps_file_dir)
            exp_file_location_db={}; exp_file_location_db[dataset_name]=fl
            parent_dir = output_dir  ### interchangable terms (parent_dir used with expression file import)
            #print groups_file_dir
        if run_from_scratch == 'Process Expression file':
            if len(input_exp_file)>0:
                if len(input_stats_file)>1: ###Make sure the files have the same arrays and order first
                    try: cel_files2, array_linker_db2 = ExpressionBuilder.getArrayHeaders(input_stats_file)
                    except Exception:
                        print_out = "Input Expression file does not have a valid format."
                        IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()               
                    if cel_files2 != cel_files:
                        print_out = "The probe set p-value file:\n"+input_stats_file+"\ndoes not have the same array order as the\nexpression file. Correct before proceeding."
                        IndicatorWindow(print_out,'Continue')
                        
                ### Check to see if a groups/comps file already exists and add file locations to 'exp_file_location_db'
                ief_list = string.split(input_exp_file,'/'); parent_dir = string.join(ief_list[:-1],'/'); exp_name = ief_list[-1]
                dataset_name = string.replace(exp_name,'exp.','')
                groups_name = 'groups.'+dataset_name; comps_name = 'comps.'+dataset_name
                batch_name = string.replace(groups_name,'groups.','batch.') ### may not apply
                groups_file_dir = parent_dir+'/'+groups_name; comps_file_dir = parent_dir+'/'+comps_name
                batch_file_dir = string.replace(groups_file_dir, 'groups.','batch.')
                fl = ExpressionFileLocationData(input_exp_file,input_stats_file,groups_file_dir,comps_file_dir)
                dataset_name = exp_name
                exp_file_location_db={}; exp_file_location_db[exp_name]=fl
            else:
                ### This occurs if running files in the ExpressionInput folder. However, if so, we won't allow for GUI based creation of groups and comps files (too complicated and confusing for user).
                ### Grab all expression file locations, where the expression, groups and comps file exist for a dataset            
                exp_file_location_db = importExpressionFiles() ###Don't create 'array_group_list', but pass the 'exp_file_location_db' onto ExpressionBuilder                

        ### Import array-group and group comparisons. Only time relevant for probesetSummarization is when an error is encountered and re-running
        try: dir_files = read_directory(parent_dir)
        except Exception: dir_files=[]

        array_group_list=[]
        array_batch_list=[]
        if backSelect == 'yes':
            for cel_file in cel_files:
                if cel_file in user_variables:
                    group_name = user_variables[cel_file]; group = ''
                else:
                    group = ''; group_name = ''    
                agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd)
                if batch_effects == 'yes' or normalize_gene_data == 'group': ### Used during backselect (must include a 'batch' variable in the stored var name)
                    if (cel_file,'batch') in user_variables:
                        batch_name = user_variables[cel_file,'batch']; batch = ''
                    else:
                        batch = ''; batch_name = ''
                    agd = ArrayGroupData(cel_file,batch,batch_name); array_batch_list.append(agd); batch_db=[]
        elif run_from_scratch == 'buildExonExportFiles':
                fl = ExpressionFileLocationData('','','',''); fl.setExonBedBuildStatus('yes'); fl.setFeatureNormalization('none')
                fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type); fl.setOutputDir(output_dir); fl.setMultiThreading(multiThreading)
                exp_file_location_db={}; exp_file_location_db[dataset_name]=fl; parent_dir = output_dir
                perform_alt_analysis = 'expression'
        elif groups_name in dir_files:
            try:
                ### Try to import any current annotations and verify that the samples indicated in the input directory are in the corresponding groups file
                array_group_list,group_db = importArrayGroupsSimple(groups_file_dir,cel_files) #agd = ArrayGroupData(array_header,group,group_name)
            except Exception,e:
                ### Over-write these annotations if theres is a problem
                for cel_file in cel_files:
                    group = ''; group_name = ''    
                    agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd); group_db=[]
            if batch_effects == 'yes' or normalize_gene_data == 'group':
                if batch_name in dir_files: ### Almost identical format and output files (import existing if present here)
                    try:
                        array_batch_list,batch_db = importArrayGroupsSimple(batch_file_dir,cel_files) #agd = ArrayGroupData(array_header,group,group_name)
                    except Exception,e:
                        for cel_file in cel_files:
                            batch = ''; batch_name = ''
                            agd = ArrayGroupData(cel_file,batch,batch_name); array_batch_list.append(agd); batch_db=[]
                else:
                    for cel_file in cel_files:
                        batch = ''; batch_name = ''    
                        agd = ArrayGroupData(cel_file,batch,batch_name); array_batch_list.append(agd); batch_db=[]
            if comps_name in dir_files and len(group_db)>0:
                comp_group_list, null = ExpressionBuilder.importComparisonGroups(comps_file_dir)
                for group1,group2 in comp_group_list:
                    try:
                        group_name1 = group_db[int(group1)]; group_name2 = group_db[int(group2)]
                        original_comp_group_list.append((group_name1,group_name2)) ### If comparisons already exist, default to these
                    except KeyError:
                        print_out = 'The "comps." file for this dataset has group numbers\nnot listed in the "groups." file.'
                        #WarningWindow(print_out,'Exit'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
                        #print print_out
                        original_comp_group_list=[]
                        
        else:
            for cel_file in cel_files:
                group = ''; group_name = ''    
                agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd)
        
        if len(array_group_list)>0: ### Thus we are not analyzing the default (ExpressionInput) directory of expression, group and comp data.
            original_option_db,original_option_list = option_db,option_list
            if len(array_group_list)>200:
                ### Only display the top 200 and don't record edits
                option_db,option_list = formatArrayGroupsForGUI(array_group_list[:200])
            else:
                option_db,option_list = formatArrayGroupsForGUI(array_group_list)
            ###Force this GUI to repeat until the user fills in each entry, but record what they did add
            user_variables_long={}
            while len(user_variables_long) != len(option_db):
                if backSelect == 'no' or 'GroupArrays' == selected_parameters[-1]:
                    selected_parameters.append('GroupArrays'); backSelect = 'no'
                    root = Tk(); root.title('AltAnalyze: Assign files to a Group Annotation'); user_variables_long={}
                    #import copy; user_variables_original = copy.deepcopy(user_variables); user_variables={}
                    gu = GUI(root,option_db,option_list['GroupArrays'],'groups')
                else: gu = PreviousResults(old_options)
                try: predictGroups = gu.Results()['PredictGroups']
                except Exception: predictGroups = False
                for option in user_variables: ### By default, all arrays will be assigned a group of ''
                    try:
                        if len(user_variables[option])>0 and 'batch' not in option:
                            if option in option_db: user_variables_long[option]=[]
                    except Exception: pass
                
                ###Store the group names and assign group numbers
                group_name_db={}; group_name_list = []; group_number = 1
                if len(array_group_list)<=200:
                    for cel_file in option_list['GroupArrays']: ### start we these CEL files, since they are ordered according to their order in the expression dataset
                        group_name = gu.Results()[cel_file]
                        if group_name not in group_name_db:
                            if group_name != 'yes' and group_name !='no': ### Results for PredictGroups
                                group_name_db[group_name]=group_number; group_number+=1
                                group_name_list.append(group_name)
                else:
                    ### For very large datasets with hundreds of samples
                    for agd in array_group_list:
                        if agd.GroupName() not in group_name_list:
                            group_name_list.append(agd.GroupName())
                            group_name_db[agd.GroupName()]=agd.Group()

                if len(group_name_db)==2: analyze_all_conditions = 'pairwise' ### Don't allow multiple comparison analysis if only two conditions present
                
                ###Store the group names and numbers with each array_id in memory
                if len(array_group_list)<=200:
                    for agd in array_group_list:
                        cel_file = agd.Array()
                        group_name = gu.Results()[cel_file] ###Lookup the new group assignment entered by the user
                        group_number = group_name_db[group_name]
                        agd.setGroupName(group_name); agd.setGroup(group_number)
                if predictGroups == 'yes':
                    predictGroups = True; break
                elif predictGroups == 'no': predictGroups = False
                elif (len(user_variables_long) != len(option_db)) or len(group_name_db)<2:
                    if len(group_name_db)<2:
                        print_out = "At least two array groups must be established\nbefore proceeding."
                    else:
                        print_out = "Not all arrays have been assigned a group. Please\nassign to a group before proceeding (required)."
                    IndicatorWindow(print_out,'Continue')   
                    option_db,option_list = formatArrayGroupsForGUI(array_group_list) ### array_group_list at this point will be updated with any changes made in the GUI by the user

            if predictGroups == False:
                exported = 0 ### Export Groups file
                if len(array_group_list)>200:
                    print 'Not storing groups due to length'
                else:
                    while exported == 0:
                        try:
                            fl = exp_file_location_db[dataset_name]; groups_file = fl.GroupsFile()
                            exportGroups(exp_file_location_db,array_group_list)
                            exported = 1
                        except Exception:                 
                            print_out = "The file:\n"+groups_file+"\nis still open. This file must be closed before proceeding"
                            IndicatorWindow(print_out,'Continue')                 
                    exported = 0
                
                if batch_effects == 'yes' or normalize_gene_data == 'group':
                    option_db,option_list = formatArrayGroupsForGUI(array_batch_list, category = 'BatchArrays')
                    ###Force this GUI to repeat until the user fills in each entry, but record what they did add
                    user_variables_long={}
                    while len(user_variables_long) != len(option_db):
                        if backSelect == 'no' or 'BatchArrays' == selected_parameters[-1]:
                            selected_parameters.append('BatchArrays'); backSelect = 'no'
                            root = Tk(); root.title('AltAnalyze: Indicate Which Batch a File is From'); user_variables_long={}
                            #import copy; user_variables_original = copy.deepcopy(user_variables); user_variables={}
                            gu = GUI(root,option_db,option_list['BatchArrays'],'batch')
                        else: gu = PreviousResults(old_options)
                        for option in user_variables: ### By default, all arrays will be assigned a batch of ''
                            try:
                                if len(user_variables[option])>0 and 'batch' in option:
                                    if option[0] in option_db: user_variables_long[option]=[]
                            except Exception: pass
                        ###Store the batch names and assign batch numbers
                        batch_name_db={}; batch_name_list = []; batch_number = 1
                        #print option_list['BatchArrays']
                        for cel_file in option_list['BatchArrays']: ### start we these CEL files, since they are ordered according to their order in the expression dataset
                            batch_name = gu.Results()[cel_file,'batch']
                            if batch_name not in batch_name_db:
                                batch_name_db[batch_name]=batch_number; batch_number+=1
                                batch_name_list.append(batch_name)
                        if len(batch_name_db)==2: analyze_all_conditions = 'pairwise' ### Don't allow multiple comparison analysis if only two conditions present
                        
                        ###Store the batch names and numbers with each array_id in memory   
                        for agd in array_batch_list:
                            cel_file = agd.Array()
                            batch_name = gu.Results()[cel_file,'batch'] ###Lookup the new batch assignment entered by the user
                            batch_number = batch_name_db[batch_name]
                            agd.setGroupName(batch_name); agd.setGroup(batch_number)
                        if (len(user_variables_long) != len(option_db)) or len(batch_name_db)<2:
                            if len(batch_name_db)<2:
                                print_out = "At least two sample batchs must be established\nbefore proceeding."
                            else:
                                print_out = "Not all arrays have been assigned a batch. Please\nassign to a batch before proceeding (required)."
                            IndicatorWindow(print_out,'Continue')   
                            option_db,option_list = formatArrayGroupsForGUI(array_batch_list, category = 'BatchArrays') ### array_batch_list at this point will be updated with any changes made in the GUI by the user
    
                    exported = 0 ### Export Batch file
                    while exported == 0:
                        try:
                            fl = exp_file_location_db[dataset_name]
                            exportGroups(exp_file_location_db,array_batch_list,filetype='Batch')
                            exported = 1
                        except Exception:                 
                            print_out = "The file:\n"+batch_file_dir+"\nis still open. This file must be closed before proceeding"
                            IndicatorWindow(print_out,'Continue')                 
                    exported = 0
                i=2; px=0 ###Determine the number of possible comparisons based on the number of groups

                while i<=len(group_name_list): px = px + i - 1; i+=1
                group_name_list.reverse(); group_name_list.append(''); group_name_list.reverse() ### add a null entry first
                if px > 150: px = 150 ### With very large datasets, AltAnalyze stalls
                possible_comps = px
                ### Format input for GUI like the imported options.txt Config file, except allow for custom fields in the GUI class
                category = 'SetupComps'; option_db={}; option_list={}; cn = 0 #; user_variables={}
                while cn < px:
                    try: group1,group2 = original_comp_group_list[cn]
                    except IndexError: group1='';group2=''
                    cn+=1; option = 'comparison '+str(cn); array_options = group_name_list; displayed_title=option; display_object='pulldown_comps'; notes=[group1,group2]
                    od = OptionData(option,displayed_title,display_object,notes,array_options,'')
                    option_db[option] = od
                    try: option_list[category].append(option) ###group is the name of the GUI menu group
                    except KeyError: option_list[category] = [option]
     
                proceed = 'no'
                while proceed == 'no' and analyze_all_conditions != 'all groups':
                    identical_groups = 'no'; comp_groups_db={}; proceed = 'no'
                    if (backSelect == 'no' or 'SetupComps' == selected_parameters[-1]):
                        selected_parameters.append('SetupComps'); backSelect = 'no'
                        root = Tk(); root.title('AltAnalyze: Establish All Pairwise Comparisons')
                        gu = GUI(root,option_db,option_list['SetupComps'],'comps')
                    else: gu = PreviousResults(old_options)
                    ### Sort comparisons from user for export
                    for comparison in gu.Results():
                        try:
                            group_name = gu.Results()[comparison]
                            if len(group_name)>0 and 'comparison' in comparison: ### Group_names are by default blank
                                cn_main,cn_minor = string.split(comparison[11:],'-') ### e.g. 1-1 and 1-2
                                try:
                                    null = int(cn_main); null = int(cn_minor)
                                    try: comp_groups_db[cn_main].append([cn_minor,group_name])
                                    except KeyError: comp_groups_db[cn_main]=[[cn_minor,group_name]]
                                except Exception: pass
                        except Exception: pass
                    print_out = "You must pick at least one comparison group before proceeding."
                    if len(comp_groups_db)>0:
                        try:
                            comp_group_list=[]
                            for cn_main in comp_groups_db:
                                cg = comp_groups_db[cn_main]
                                cg.sort()
                                comp_group_list.append([cn_main,[group_name_db[cg[0][1]],group_name_db[cg[1][1]]]])
                                if cg[0][1] == cg[1][1]: identical_groups = 'yes' ### Thus the two groups in the comparisons are identical, flag
                            comp_group_list.sort()                    
                            proceed = 'yes'
                        except Exception:
                            print traceback.format_exc()
                            print_out = "You must pick at least two groups for each comparison."
                    if identical_groups == 'yes': proceed = 'no'; print_out = "The same group is listed as both the experimental and\ncontrol group in a comparison. Fix before proceeding."
                    if proceed == 'no': IndicatorWindow(print_out,'Continue')   
                    
                ### Export user modified comps files
                while exported == 0:
                    try:
                        fl = exp_file_location_db[dataset_name]; comps_file = fl.CompsFile()
                        if analyze_all_conditions != 'all groups': exportComps(exp_file_location_db,comp_group_list)
                        exported = 1
                    except Exception:
                        print_out = "The file:\n"+comps_file+"\nis still open. This file must be closed before proceeding"
                        IndicatorWindow(print_out,'Continue')

        ### See if there are any Affymetrix annotation files for this species
        import_dir = '/AltDatabase/affymetrix/'+species
        try: dir_list = read_directory(import_dir); fn_dir = filepath(import_dir[1:]); species_dir_found = 'yes'
        except Exception: fn_dir = filepath(import_dir); dir_list = []; species_dir_found = 'no'

        ### Used to check if the user has an Affymetrix CSV file around... no longer needed        
        """
        if (len(dir_list)<1 or species_dir_found == 'no') and array_type != 'exon':
            print_out = 'No Affymetrix annnotations file found in the directory:\n'+fn_dir
            print_out += '\n\nTo download, click on the below button, find your array and download the annotation CSV file'
            print_out += '\nlisted under "Current NetAffx Annotation Files". Extract the compressed zip archive to the'
            print_out += '\nabove listed directory and hit continue to include these annotations in your results file.'
    
            button_text = 'Download Annotations'; url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
            IndicatorLinkOutWindow(print_out,button_text,url)
            """
    """ ### Change made in version 2.1.2
    if len(input_fastq_dir)>0:
        array_type = "3'array"
        vendor = 'other:Ensembl' ### Ensembl linked system name
    """
    if microRNA_prediction_method == 'two or more': microRNA_prediction_method = 'multiple'
    else: microRNA_prediction_method = 'any'

    try: permute_p_threshold = float(permute_p_threshold)
    except ValueError: permute_p_threshold = permute_p_threshold
    
    try: dabg_p = float(dabg_p)
    except ValueError: dabg_p = dabg_p
    
    try: expression_threshold = float(expression_threshold)
    except ValueError: expression_threshold = expression_threshold

    try: alt_exon_fold_cutoff = float(alt_exon_fold_cutoff)
    except ValueError: alt_exon_fold_cutoff = alt_exon_fold_cutoff

    try: gene_expression_cutoff = float(gene_expression_cutoff)
    except ValueError: gene_expression_cutoff = gene_expression_cutoff    

    ### Find the current verison of APT (if user deletes location in Config file) and set APT file locations
    
    try: apt_location = getAPTLocations(file_location_defaults,run_from_scratch,run_MiDAS)
    except Exception: pass

    ### Set the primary parent directory for ExpressionBuilder and AltAnalyze (one level above the ExpressionInput directory, if present)
    for dataset in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]
        try: fl.setAPTLocation(apt_location)
        except Exception: pass
        if run_from_scratch == 'Process CEL files' or 'Feature Extraction' in run_from_scratch:
            fl.setInputCDFFile(input_cdf_file); fl.setCLFFile(clf_file); fl.setBGPFile(bgp_file); fl.setXHybRemoval(remove_xhyb)
            fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type); fl.setOutputDir(output_dir)
            fl.setChannelToExtract(channel_to_extract)
        elif 'Chromium' in run_from_scratch:
            fl.setChromiumSparseMatrix(sparse_matrix_file)
            #print fl.ChromiumSparseMatrix()
        elif run_from_scratch == 'Process RNA-seq reads':
            fl.setCELFileDir(cel_file_dir); fl.setOutputDir(output_dir); fl.setExonBedBuildStatus(build_exon_bedfile)
            fl.setRunKallisto(input_fastq_dir);
        if array_type != 'gene' and array_type != "3'array":
            compendiumPlatform = 'exon'
        fl = exp_file_location_db[dataset]; fl.setRootDir(parent_dir)
        fl.setFeatureNormalization(normalize_feature_exp)
        fl.setNormMatrix(normalize_gene_data)
        fl.setProbabilityStatistic(probability_algorithm)
        fl.setBatchEffectRemoval(batch_effects)
        fl.setMarkerFinder(marker_finder)
        fl.setProducePlots(visualize_results)
        fl.setPerformLineageProfiler(run_lineage_profiler)
        fl.setCompendiumType(compendiumType)
        fl.setCompendiumPlatform(compendiumPlatform)
        fl.setVendor(vendor)
        try: fl.setFDRStatistic(FDR_statistic)
        except Exception: pass
        try: fl.setExcludeLowExpressionExons(excludeNonExpExons)
        except Exception: fl.setExcludeLowExpressionExons(True)
        try: fl.setPredictGroups(predictGroups)
        except Exception: fl.setPredictGroups(False)
        try: fl.setPredictGroupsParams(gsp)
        except Exception: pass
        fl.setMultiThreading(multiThreading)
        if run_from_scratch == 'Process Expression file':
            fl.setRootDir(output_dir) ### When the data is not primary array data files, allow for option selection of the output directory
            fl.setOutputDir(output_dir)
        try: fl.setRPKMThreshold(rpkm_threshold)
        except Exception: pass
        try: fl.setGeneExpThreshold(gene_exp_threshold)
        except Exception: pass
    if array_type == 'RNASeq': ### Post version 2.0, add variables in fl rather than below
        fl.setRPKMThreshold(rpkm_threshold)
        fl.setExonExpThreshold(exon_exp_threshold)
        fl.setGeneExpThreshold(gene_exp_threshold)
        fl.setExonRPKMThreshold(exon_rpkm_threshold)
        fl.setJunctionExpThreshold(expression_threshold)
    try: fl.setMLP(mlp)
    except Exception: pass

    if predictGroups:
        ### Single-Cell Analysis Parameters
        try: option_db,option_list=original_option_db,original_option_list ### was re-set above... needed to get the propper data from the last loop
        except Exception: option_list,option_db = importUserOptions(array_type,vendor=vendor)
        selected_parameters.append('PredictGroups')
        supported_geneset_types = getSupportedGeneSetTypes(species,'gene-mapp')
        supported_geneset_types += getSupportedGeneSetTypes(species,'gene-go')
        option_db['GeneSetSelectionPredict'].setArrayOptions(['None Selected']+supported_geneset_types)
        option_db['PathwaySelectionPredict'].setArrayOptions(['None Selected'])
        #option_db['PathwaySelection'].setArrayOptions(supported_genesets)

        status = 'repeat'
        while status == 'repeat':
            root = Tk()
            root.title('AltAnalyze: Predict Cell Populations')
            ### Run in GUI and wait to be executed
            gu = GUI(root,option_db,option_list['PredictGroups'],'')
            ### Permission to run full analsyis is granted, proceed
            gsp = gu.Results()['gsp']
            status = 'continue'

        import RNASeq  
        expFile = fl.ExpFile()
        mlp_instance = fl.MLP()
        
        global logfile
        root_dir = export.findParentDir(expFile)
        root_dir = string.replace(root_dir,'/ExpressionInput','')
        time_stamp = AltAnalyze.timestamp()    
        logfile = filepath(root_dir+'AltAnalyze_report-'+time_stamp+'.log')
        
        count = verifyFileLength(expFile[:-4]+'-steady-state.txt')

        if count>1:
            expFile = expFile[:-4]+'-steady-state.txt'
        elif array_type=='RNASeq' or len(input_fastq_dir)>0 or len(sparse_matrix_file)>0:
            ### Indicates that the steady-state file doesn't exist. The exp. may exist, be could be junction only so need to re-build from bed files here
            values = species,exp_file_location_db,dataset,mlp_instance
            StatusWindow(values,'preProcessRNASeq') ### proceed to run the full discovery analysis here!!!
            if array_type=='RNASeq':
                expFile = expFile[:-4]+'-steady-state.txt'
        """
        else:
            print_out = 'WARNING... Prior to running ICGS, you must first run AltAnalyze\nusing assigned groups for this array type.'
            IndicatorWindow(print_out,'Continue')
            AltAnalyze.AltAnalyzeSetup((selected_parameters[:-1],user_variables)); sys.exit()"""
            
        values = expFile, mlp_instance, gsp, False
        StatusWindow(values,'predictGroups') ### proceed to run the full discovery analysis here!!!
        if len(graphic_links)>0:
            root = Tk()
            root.title('AltAnalyze: Evaluate ICGS Clustering Results')
            ### Review results in custom GUI for predicting groups
            gu = GUI(root,'PredictGroups',[],'')
            nextStep = gu.Results()['next']
            group_selected = gu.Results()['group_select']
            if nextStep == 'UseSelected':
                group_selected = group_selected[:-4]+'.txt'
                exp_file = fl.ExpFile()
                try:
                    exonExpFile,newExpFile,new_groups_dir = exportAdditionalICGSOutputs(expFile,group_selected,outputTSNE=False)
                    for exp_name in exp_file_location_db: break ### get name
                    fl.setExpFile(exonExpFile) ### Use the ICGS re-ordered and possibly OutlierFiltered for downstream analyses
                    comps_file = string.replace(newExpFile,'exp.','comps.')
                    fl.setGroupsFile(new_groups_dir)
                    fl.setCompsFile(string.replace(new_groups_dir,'groups.','comps.'))
                    del exp_file_location_db[exp_name]
                    exp_file_location_db[exp_name+'-ICGS'] = fl
                except Exception:
                    print traceback.format_exc()
                    pass ### Unknown error
                                    
                run_from_scratch = 'Process Expression file'
            else:
                #print 're-initializing window'
                AltAnalyze.AltAnalyzeSetup((selected_parameters,user_variables)); sys.exit()
        else:
            AltAnalyze.AltAnalyzeSetup((selected_parameters,user_variables)); sys.exit()
                           
    expr_var = species,array_type,vendor,constitutive_source,dabg_p,expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch, perform_alt_analysis
    alt_var = analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_cutoff,gene_expression_cutoff,remove_intronic_junctions,permute_p_threshold, perform_permutation_analysis, export_splice_index_values, analyze_all_conditions
    additional_var = calculate_splicing_index_p, run_MiDAS, analyze_functional_attributes, microRNA_prediction_method, filter_for_AS, additional_algorithms
    goelite_var = ge_fold_cutoffs,ge_pvalue_cutoffs,ge_ptype,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,pathway_permutations,mod,returnPathways

    return expr_var, alt_var, additional_var, goelite_var, exp_file_location_db

def getAPTLocations(file_location_defaults,run_from_scratch,run_MiDAS):
    from import_scripts import ResultsExport_module
    if 'APT' in file_location_defaults:
        fl = file_location_defaults['APT']
        apt_location = fl.Location() ###Only one entry for all species
        if len(apt_location)<1: ###If no APT version is designated, prompt the user to find the directory
            if run_from_scratch == 'CEL_summarize':
                print_out = 'To proceed with probeset summarization from CEL files,\nyou must select a valid Affymetrix Power Tools Directory.'
            elif run_MiDAS == 'yes': 
                print_out = "To proceed with running MiDAS, you must select\na valid Affymetrix Power Tools Directory."
            win_info = IndicatorChooseWindow(print_out,'Continue') ### Prompt the user to locate the APT directory
            apt_location = win_info.Folder()
            fl.SetLocation(apt_location)
            try: exportDefaultFileLocations(file_location_defaults)
            except Exception: pass
    return apt_location

def check_moderated_support(option_db):
    """ Excludes moderated t-test support when module import fails... shouldn't fail """
    try:
        from stats_scripts import mpmath
    except Exception,e:
        a = traceback.format_exc()
        GUIcriticalError(a)
        keep=[]
        od  = option_db['probability_algorithm']
        for i in od.ArrayOptions():
            if 'oderate' not in i: keep.append(i)
        od.setArrayOptions(keep) ### remove any moderated stats
        od.setDefaultOption(keep[0]) ### Change the default value to one of those in keep
    return option_db

def GUIcriticalError(log_report):
    log_file = filepath('GUIerror.log')
    data = open(log_file,'w')
    data.write(log_report); data.close()
    """
    if os.name == 'nt':
            try: os.startfile('"'+log_file+'"')
            except Exception:  os.system('open "'+log_file+'"')
    elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
    elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
    """
    
class GeneSelectionParameters:
    ### This class specifies parameters for filtering a large dataset for downstream gene or pathway analysis/visualization
    def __init__(self, species, platform, vendor): 
        self._species = species; self._platform = platform; self._vendor = vendor
        self._PathwaySelect = False; self._gene = False; self._GeneSet = False
        self._Normalize = False
    def Species(self): return self._species
    def Platform(self): return self._platform
    def Vendor(self): return self._vendor
    def setGeneSet(self,gene_set): self._GeneSet = gene_set
    def GeneSet(self):
        if isinstance(self._GeneSet, tuple) or isinstance(self._GeneSet, list):
            return tuple(self._GeneSet)
        else:
            return self._GeneSet
    def setPathwaySelect(self,pathway): self._PathwaySelect = pathway
    def setJustShowTheseIDs(self, justShowTheseIDs): self.justShowTheseIDs = justShowTheseIDs
    def setClusterGOElite(self, clusterGOElite): self.clusterGOElite = clusterGOElite
    def setStoreGeneSetName(self, geneSetName): self.geneSetName = geneSetName
    def StoreGeneSetName(self): return self.geneSetName
    def ClusterGOElite(self):
        if isinstance(self.clusterGOElite, tuple) or isinstance(self.clusterGOElite, list):
            self.clusterGOElite = list(self.clusterGOElite)
            if 'None Selected' in self.clusterGOElite: self.clusterGOElite.remove('None Selected')
            return self.clusterGOElite
        else:
            self.clusterGOElite = string.replace(self.clusterGOElite,'None Selected','')
            return [self.clusterGOElite]
    def PathwaySelect(self):
        if isinstance(self._PathwaySelect, tuple) or isinstance(self._PathwaySelect, list):
            return tuple(self._PathwaySelect)
        else:
            return self._PathwaySelect
    def setGeneSelection(self,gene): self._gene = gene
    def GeneSelection(self):
        try:
            genes = self._gene
            genes = string.replace(genes,'\r', ' ')
            genes = string.replace(genes,'\n', ' ') 
        except Exception:
            genes = self._gene
        return genes
    def JustShowTheseIDs(self):
        if 'None Selected' in self.justShowTheseIDs:
            return ''
        else:
            justShowTheseIDs = string.replace(self.justShowTheseIDs,'\r',' ')
            justShowTheseIDs = string.replace(justShowTheseIDs,'\n',' ')
            justShowTheseIDs = string.split(justShowTheseIDs,' ')
            try: justShowTheseIDs.remove('')
            except Exception: pass
            return justShowTheseIDs
    def GetGeneCorrelations(self):
        if len(self._gene)>0: return True
        else: return False
    def FilterByPathways(self):
        if self._GeneSet != 'None Selected': return True
        else: return False
    def setTranspose(self,transpose): self._transpose = transpose
    def Transpose(self):
        try: return self._transpose
        except Exception: return False
    def setOntologyID(self,OntologyID): self._OntologyID = OntologyID
    def OntologyID(self):
        try:
            return self._OntologyID
        except Exception: return ''
    def setIncludeExpIDs(self,IncludeExpIDs): self._IncludeExpIDs = IncludeExpIDs
    def IncludeExpIDs(self): return self._IncludeExpIDs
    def setNormalize(self,Normalize):
        if Normalize == 'NA': Normalize = False
        self._Normalize = Normalize
    def Normalize(self): return self._Normalize
    def setK(self,k): self.k = k
    def k(self):
        try: return self.k
        except: return None
    def K(self):
        try: return self.k
        except: return None
    def setExcludeGuides(self,excludeGuides): self.excludeGuides = excludeGuides
    def ExcludeGuides(self): return self.excludeGuides
    def setSampleDiscoveryParameters(self,ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering,dynamicCorrelation,
                removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff):
        ### For single-cell RNA-Seq data
        self.expressionCutoff = ExpressionCutoff
        self.countsCutoff = CountsCutoff
        self.rho_cutoff = rho_cutoff
        self.foldDiff = FoldDiff
        self.samplesDiffering = SamplesDiffering
        self.featurestoEvaluate = featurestoEvaluate
        self.restrictBy = restrictBy
        self.excludeCellCycle = excludeCellCycle
        self.column_metric = column_metric
        self.column_method = column_method
        self.removeOutliers = removeOutliers
        self.dynamicCorrelation = dynamicCorrelation
        if len(self._gene)>0:
            self._gene = self._gene + ' amplify' ### always amplify the selected genes if any
    def setExpressionCutoff(self,expressionCutoff):self.expressionCutoff = expressionCutoff
    def setCountsCutoff(self,countsCutoff):self.countsCutoff = countsCutoff
    def ExpressionCutoff(self):
        try: return float(self.expressionCutoff)
        except Exception: return False
    def setRhoCutoff(self,rho):
        self.rho_cutoff = rho
    def RhoCutoff(self):
        return float(self.rho_cutoff)
    def CountsCutoff(self):
        try: return int(float(self.countsCutoff))
        except Exception: return False
    def FoldDiff(self):
        try: return float(self.foldDiff)
        except Exception: return False
    def SamplesDiffering(self):
        try: return int(float(self.samplesDiffering))
        except Exception: return False
    def dynamicCorrelation(self):
        if self.dynamicCorrelation=='yes' or self.dynamicCorrelation==True:
            return True
        else:
            return False
    def amplifyGenes(self):
        if (self.FilterByPathways() != '' and self.FilterByPathways() !=False) or (self.GeneSelection() != '' and self.GeneSelection() != ' amplify'):
            return True
        else: return False
    def FeaturestoEvaluate(self): return self.featurestoEvaluate
    def RestrictBy(self):
        if self.restrictBy == True or self.restrictBy == 'yes' or self.restrictBy == 'protein_coding':
            return 'protein_coding'
        else:
            return None
    def RemoveOutliers(self):
        if self.removeOutliers == True or self.removeOutliers == 'yes':
            return True
        else:
            return False
    def ExcludeCellCycle(self):
        if self.excludeCellCycle == 'stringent' or self.excludeCellCycle == 'strict':
            return 'strict'  ### Also includes removing drivers correlated to any cell cycle genes, not just in the training set
        elif self.excludeCellCycle == False:
            return False
        elif self.excludeCellCycle == True or self.excludeCellCycle != 'no':
            return True
        else:
            return False
        
    def ColumnMetric(self): return self.column_metric
    def ColumnMethod(self): return self.column_method
    def MinEvents(self):
        return self.SamplesDiffering()-1
    def MedEvents(self):
        return (self.SamplesDiffering()-1)*2
    
def getSupportedGeneSetTypes(species,directory):
    try:
        geneset_types=[]
        current_geneset_dirs = unique.read_directory('/AltDatabase/goelite/'+species+'/'+directory)
        for geneset_dir in current_geneset_dirs:
            geneset_dir = string.join(string.split(geneset_dir,'-')[1:],'-')[:-4] ### remove the prefix gene system
            if geneset_dir == 'MAPP': geneset_dir = 'WikiPathways'
            if geneset_dir not in geneset_types:
                if len(geneset_dir)>1:
                    geneset_types.append(geneset_dir)
    except Exception:
        return []
    return geneset_types

def getSupportedGeneSystems(species,directory):
    system_names=[]
    
    current_system_dirs = unique.read_directory('/AltDatabase/goelite/'+species+'/'+directory)
    for system_dir in current_system_dirs:
        try:
            system_dir = string.split(system_dir,'-')[1][:-4] ### remove the prefix gene system
            if len(system_dir)>1:
                if system_dir not in system_names:
                    system_names.append(system_dir)
        except Exception: None
    system_names.append('Ensembl')
    system_names.append('HMDB')
    system_names = unique.unique(system_names)
    system_names.sort()
    return system_names

def listAllGeneSetCategories(species,geneset_type,directory):
    geneset_categories=[]
    if directory == 'gene-go':
        if geneset_type == 'GeneOntology': geneset_type = 'go'
        filename = 'AltDatabase/goelite/OBO/builds/'+geneset_type+'_annotations.txt'
        index = 1
    else:
        if geneset_type == 'WikiPathways': geneset_type = 'MAPP'
        filename = 'AltDatabase/goelite/'+species+'/'+directory+'/'+'Ensembl-'+geneset_type+'.txt'
        index = -1
        
    fn=filepath(filename)
    ### Imports a geneset category and stores pathway-level names
    i=0
    for line in open(fn,'rU').xreadlines():
        if i==0: i=1 ### Skip the header
        else:
            data = cleanUpLine(line)
            geneset_category = string.split(data,'\t')[index]
            if geneset_category not in geneset_categories:
                geneset_categories.append(geneset_category)
    geneset_categories.sort()
    return geneset_categories

def getValidExpFile(altanalyze_rawexp_dir):
    dir_files = read_directory(altanalyze_rawexp_dir)
    valid_file = ''
    for file in dir_files:
        if 'exp.' in file and 'state.txt' not in file and 'highExp' not in file:
            valid_file = altanalyze_rawexp_dir+'/'+file
            break
    return valid_file
        
def getValidSplicingScoreFile(altanalyze_rawsplice_dir):
    valid_dirs = ['splicing-index','FIRMA','ASPIRE','linearregres']
    dir_files = read_directory(altanalyze_rawsplice_dir)
    valid_folder = None
    for folder in valid_dirs:
        if folder in dir_files:
            valid_folder = folder
            break
    valid_file_dir = ''
    primary=''
    if valid_folder != None:
        child_dir = altanalyze_rawsplice_dir+'/'+valid_folder
        dir_files = read_directory(altanalyze_rawsplice_dir+'/'+valid_folder)
        for file in dir_files:
            if '.txt' in file:
                valid_file_dir = child_dir+'/'+file
                if '_vs_' not in file: ### You can have a folder with pairwise comps and all groups
                    primary = child_dir+'/'+file
    if len(primary)!=0: valid_file_dir = primary
    return valid_file_dir

def downloadInteractionDBs(species,windowType):
    analysis = 'getAdditionalOnlineResources' ### same option as updating gene-sets
    additional_resources=['Latest WikiPathways','KEGG','BioGRID','DrugBank','miRNA Targets','Transcription Factor Targets']
    get_additional = 'customSet',additional_resources
    values = species,get_additional
    StatusWindow(values,analysis,windowType=windowType) ### open in a TopLevel TK window (don't close current option selection menu)
    
if __name__ == '__main__':
    """
    dir = '/Users/saljh8/Desktop/dataAnalysis/FuKun/AltResults/Clustering/Combined-junction-exon-evidence.txt'
    expFile = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/ExpressionInput/exp.GGI-IG2.txt'
    group_selected = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/Kallisto/ICGS-NMF/FinalMarkerHeatmap.txt'
    array_type="3'array"; species='Mm'
    exportAdditionalICGSOutputs(expFile,group_selected,outputTSNE=True)
    sys.exit()"""
    #a = exportJunctionList(dir,limit=50)
    #print a;sys.exit()
    try:
        import multiprocessing as mlp
        mlp.freeze_support()
    except Exception:
        print 'Note: Multiprocessing not supported for this verison python.'
        mlp=None

    #getUpdatedParameters(array_type,species,run_from_scratch,file_dirs)    
    a = getUserParameters('yes')

    
