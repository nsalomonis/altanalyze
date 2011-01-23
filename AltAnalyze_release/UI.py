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
import statistics
import sys, string
import shutil
import os.path
import unique
import update; reload(update)
import export
import ExpressionBuilder
import time
import webbrowser
from sys import argv

try:
    import Tkinter 
    #import bwidget; from bwidget import *
    from Tkinter import *
    import PmwFreeze
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
    
def identifyCELfiles(dir):
    dir_list = read_directory(dir); dir_list2=[]; full_dir_list=[]
    for file in dir_list:
        file_lower = string.lower(file)
        if '.cel' in file_lower[-4:] and '.cel.' not in file_lower:
            dir_list2.append(file)
            file = dir+'/'+file
            full_dir_list.append(file)
    dir_list2.sort(); full_dir_list.sort()
    return dir_list2,full_dir_list

def identifyArrayType(full_dir_list):
    arrays={} ### Determine the type of unique arrays in each directory
    for filename in full_dir_list:
        fn=filepath(filename); ln=0
        for line in open(fn,'rU').xreadlines():
            if ln<150:
                data = cleanUpLine(line); ln+=1
                if 'sq' in data:
                    try:
                        array_info,null = string.split(data,'sq')
                        array_info = string.split(array_info,' ')
                        array_type = array_info[-1]
                        if '.' in array_type: array_type,null = string.split(array_type,'.')
                        arrays[array_type]=[]
                        break
                    except Exception: null=[]
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
    for array in arrays: array_ls.append(array)
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
    if backSelect == 'yes': warn = 'no'
    else: warn = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
        input_cdf_file = filename
        if '.pgf' in input_cdf_file:
            ###Check to see if the clf and bgp files are present in this directory 
            icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
            clf_short = string.replace(cdf_short,'.pgf','.clf')
            if array_type == 'exon' or array_type == 'junction': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
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
            dw = DownloadWindow(print_out,'Download by AltAnalyze','Select Local Files'); warn = 'no'
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
            else: filenames = [input_cdf_file]
            for filename in filenames:
                var_list = filename,'LibraryFiles'
                if debug_mode == 'no': StatusWindow(var_list,'download')
                else:
                    for filename in filenames:
                        continue_analysis = update.downloadCurrentVersion(filename,'LibraryFiles','')
                try: os.remove(filepath(filename)) ### Not sure why this works now and not before
                except Exception: null=[]
        else: library_dir = ''
    
    filename = 'AltDatabase/affymetrix/'+species+'/'+annot_file
    fn=filepath(filename); annotation_dir = filename
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception:
        if warn == 'yes':
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
                        except Exception: update.downloadCurrentVersion(filename,'AnnotationFiles',None)
                try: os.remove(filepath(filename)) ### Not sure why this works now and not before
                except Exception: null=[]
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
    def __init__(self,info_list,analysis_type):
        try:
            root = Tk()
            self._parent = root
            root.title('AltAnalyze 1.155')
            statusVar = StringVar() ### Class method for Tkinter. Description: "Value holder for strings variables."

            height = 250; width = 700
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = 'Download File Status Window',
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            
            group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Output')
            group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                
            Label(group.interior(),width=180,height=152,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)

            status = StringVarFile(statusVar,root) ### Likely captures the stdout
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
            file_location_defaults,db_version,new_species_codes = info_list
            try: sys.stdout = status; root.after(100,getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,self._parent))
            except Exception,e: getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,None)
        try: self._parent.mainloop()
        except Exception: None

    def deleteWindow(self): tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
    def quit(self):
        try: self._parent.quit(); self._parent.destroy(); sys.exit()
        except Exception: self._parent.quit(); sys.exit()
    
class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s):
        new = self.__stringvar.get()
        for c in s:
            #if c == '\n': self.__newline = 1
            if c == '\k': self.__newline = 1### This should not be found and thus results in a continous feed rather than replacing a single line
            else:
                if self.__newline: new = ""; self.__newline = 0
                new = new+c
        self.set(new)
    def set(self,s): self.__stringvar.set(s); self.__window.update()
    def get(self): return self.__stringvar.get()

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
	
class GUI:
    def __init__(self, parent, option_db, option_list, defaults): 
        self._parent = parent; self._option_list = option_list; self._option_db = option_db
        self._user_variables = user_variables; self.pathdb={}; i = -1
        enter_index=0; radio_index=0; dropdown_index=0; check_index=0 ### used to keep track of how many enter boxes we have
        self.default_dir = PathDir; self.default_file = PathFile
        
        filename = 'Config/icon.gif'; orient_type = 'left'
        if 'input_cel_dir' in option_list: filename = 'Config/aa_0.gif'
        if 'include_raw_data' in option_list: filename = 'Config/aa_1.gif'; orient_type = 'top'
        if 'filter_for_AS' in option_list: filename = 'Config/aa_2.gif'; orient_type = 'top'
        if 'pathway_permutations' in option_list: filename = 'Config/goelite.gif'
        
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        #except Exception: print filename; 'what?';kill
        self.pathdb={}; use_scroll = 'no'

        #if defaults == 'groups' or defaults == 'comps' or 'filter_for_AS' in option_list:
        if defaults != 'null':
            height = 400; width = 400
            if defaults == 'groups':
                notes = "For each CEL file, type in a name for the group it belongs to\n(e.g., 24hrs, 48hrs, 4days, etc.)."
                Label(self._parent,text=notes).pack(); label_text_str = 'AltAnalyze Group Names'
                if len(option_list)<15: height = 320; width = 320
            elif defaults == 'comps':
                notes = "Experimental Group\t\t\tControl Group     "
                label_text_str = 'AltAnalyze Pairwise Group Comparisons'
                if len(option_list)<5: height = 250; width = 400
            elif 'filter_for_AS' in option_list:
                label_text_str = 'AltAnalyze Alternative Exon Analysis Parameters'
                height = 350; width = 400; use_scroll = 'yes'
                if os.name != 'nt': width+=50
            elif 'pathway_permutations' in option_list:
                label_text_str = 'GO-Elite Parameters'
                height = 400; width = 425; use_scroll = 'yes'
            elif 'expression_data_format' in option_list:
                label_text_str = "AltAnalyze Expression Dataset Parameters"
                height = 350; width = 400; use_scroll = 'yes'
                if os.name != 'nt': width+=50
            else:
                label_text_str = "AltAnalyze Main Dataset Parameters"
                height = 300; width = 400; use_scroll = 'yes'
                
            if os.name != 'nt':height+=50; width+=100
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = label_text_str,
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            if defaults == 'comps':
                Label(self.frame,text=notes).pack()

        create_group = 'yes'
        if 'pathway_permutations' in option_list or 'expression_data_format' in option_list or 'filter_probeset_types' in option_list:
            if 'ge_ptype' in option_list:
                self.group_tag = 'GO-Elite Gene Expression Analysis Filters'; od = option_db['dabg_p']
                if od.ArrayOptions() == ['NA']: create_group = 'no'
            elif 'pathway_permutations' in option_list:
                self.group_tag = 'GO-Elite Over-Representation and Filtering Parameters'
            if 'dabg_p' in option_list:
                self.group_tag = 'Probesets Filtering Options'
                option_db['dabg_p']; od = option_db['dabg_p']
                if od.ArrayOptions() == ['NA']: create_group = 'no'
            elif 'expression_data_format' in option_list:
                self.group_tag = 'Gene Expression Analysis Options'
            if 'filter_probeset_types' in option_list:
                self.group_tag = 'Primary Alternative Exon Parameters'
            if create_group == 'yes': 
                custom_group = PmwFreeze.Group(self.sf.interior(),tag_text = self.group_tag)
                custom_group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
                insert_into_group = 'yes'
            else: insert_into_group = 'no'
        else: insert_into_group = 'no'
                
        object_directions = ['top','bottom','up','down']
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
                if 'pathway_permutations' in option_list or 'new_run' in option_list: orientation = 'vertical'
                elif 'run_from_scratch' in option_list: orientation = 'vertical'
                else: orientation = 'horizontal'
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
                
                #l = Label(group.interior(), text=self.title); l.pack(side=LEFT)        
                entry = Entry(group.interior(),textvariable=self.pathdb[option]);
                entry.pack(side='left',fill = 'both', expand = 1, padx = 10, pady = 2)
                button = Button(group.interior(), text="select "+od.DisplayObject(), width = 10, fg="red", command=filecallback); button.pack(side=LEFT, padx = 2,pady = 2)                    

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

                if 'species' in option:
                    if 'selected_species2' in option:
                        self.speciescomp2 = self.comp; self.speciescomp2.pack(anchor = 'w', padx = 10, pady = 0)
                    elif 'selected_species3' in option:
                        self.speciescomp3 = self.comp; self.speciescomp3.pack(anchor = 'w', padx = 10, pady = 0)
                    else: self.speciescomp = self.comp; self.speciescomp.pack(anchor = 'w', padx = 10, pady = 0)
                    self.speciescomp.invoke(selected_default) 
                elif 'array_type' in option:
                    self.arraycomp = self.comp; self.arraycomp.pack(anchor = 'w', padx = 10, pady = 0)
                    self.arraycomp.invoke(selected_default)
                elif 'manufacturer_selection' in option:
                    self.vendorcomp = self.comp; self.vendorcomp.pack(anchor = 'w', padx = 10, pady = 0)
                    self.vendorcomp.invoke(selected_default)
                else:
                    if insert_into_group == 'no':
                        if 'version' in option: pady_int = 0
                        else: pady_int = 1
                        self.comp.pack(anchor = 'w', padx = 10, pady = pady_int)
                    elif dropdown_index == 1: comp1 = self.comp
                    elif dropdown_index == 2: comp2 = self.comp
                    elif dropdown_index == 3: comp3 = self.comp
                    elif dropdown_index == 4: comp4 = self.comp
                    elif dropdown_index == 5: comp5 = self.comp
                    elif dropdown_index == 6: comp6 = self.comp
                    elif dropdown_index == 7: comp7 = self.comp
                    try: self.comp.invoke(selected_default)
                    except Exception: print selected_default, option;kill
                    if option == 'selected_version':
                        notes = 'Note: Available species may vary based on database selection    \n'
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
                except Exception: null=[]
                try: self.comp2.invoke(notes[1])
                except Exception: null=[]
                
            if 'simple_entry' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                if len(override_default)>0: self.default_option = override_default
                else: self.default_option = self.display_options[0]
                def enter_callback(tag,enter_callback=self.enter_callback,option=option):
                    enter_callback(tag,option)
                self.title = self.title + '\t '
                self.entry_field = PmwFreeze.EntryField(self.sf.interior(),
                        labelpos = 'w', label_text = self.title,
                        validate = enter_callback,
                        value = self.default_option
                ); self.entry_field.pack(padx = 10, pady = 1)

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
                try:
                    if float(self.default_option) <= 1: use_method = 'p'
                    else: use_method = 'i'
                except ValueError:
                    #self.default_option = 'CHANGE TO A NUMERIC VALUE'; use_method = 'i'
                    self.default_option = string.replace(self.default_option,'---','')
                    use_method = 'i'
		  
		"""
		if use_method == 'i' or use_method == 'p':
		    l = Label(parent_type, text=self.title); l.pack()         
		    self.entry_field = Entry(parent_type,textvariable=self.default_option);
		    self.entry_field.pack()"""
		    
                if use_method == 'p':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w', label_text = self.title, validate = custom_validate_p, 
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge')                   
                if use_method == 'i':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w', label_text = self.title, validate = custom_validate,
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge')
                if insert_into_group == 'no': self.entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)	
                elif enter_index == 1: self.entry_field1 = self.entry_field
                elif enter_index == 2: self.entry_field2 = self.entry_field
                elif enter_index == 3: self.entry_field3 = self.entry_field
                elif enter_index == 4: self.entry_field4 = self.entry_field  
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
                        label_text = self.title, hull_borderwidth = 2, hull_relief = 'ridge')
                self.checkbuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                ### Add some buttons to the checkbutton RadioSelect.
                for text in self.display_options:
                     if text != ['NA']: self.checkbuttons.add(text)
                self.checkbuttons.invoke(self.default_option)
                self.checkbuttons.invoke(self.default_option2)
                if len(notes)>0: Label(self._parent, text=notes).pack()
                
            if 'single-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                if defaults == 'comps': parent_type = self._parent; orient_type = 'top'
                if insert_into_group == 'yes': parent_type = custom_group.interior(); check_index+=1
                
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
                        if insert_into_group == 'no': self.checkbuttons.pack(side = orient_type, expand = 1, padx = 10, pady = 10)
                        elif check_index == 1:  checkbuttons1 = self.checkbuttons
                        elif check_index == 2:  checkbuttons2 = self.checkbuttons
                        elif check_index == 3:  checkbuttons3 = self.checkbuttons
                        elif check_index == 4:  checkbuttons4 = self.checkbuttons
                        ### Add some buttons to the checkbutton RadioSelect.
                        self.checkbuttons.add(self.title)
                        if self.default_option == 'yes': self.checkbuttons.invoke(self.title)
                        else: self._user_variables[option] = 'no'
                if len(notes)>0: Label(self._parent, text=notes).pack()

            custom_group_endpoints = ['ge_ptype', 'mod', 'expression_threshold', 'run_goelite', 'gene_expression_cutoff', 'microRNA_prediction_method']
            if option in custom_group_endpoints and insert_into_group == 'yes':
                ### This is employed when we want to place several items into a group frame together.
                ### Since this is a generic class, we need to setup special cases to do this, however,
                ### this same code could be used in other instances as well
                reorganize = 'no'
                self.group_tag = 'GO-Elite Over-Representation and Filtering Parameters'; pady_int = 5
                if 'run_goelite' in option_list: self.group_tag = 'Gene Expression Analysis Options'; pady_int = 1                
                if 'microRNA_prediction_method' in option_list: self.group_tag = 'Advanced Options'; pady_int = 1; reorganize = 'yes'

                try: checkbuttons1.pack(side = 'top', expand = 1, padx = 9, pady = 0)
                except Exception: null=[]
                try: checkbuttons2.pack(side = 'top', expand = 1, padx = 9, pady = 0)
                except Exception: null=[]
                try: checkbuttons3.pack(side = 'top', expand = 1, padx = 9, pady = 0)
                except Exception: null=[]
                try: checkbuttons4.pack(side = 'top', expand = 1, padx = 9, pady = 0)
                except Exception: null=[]
                try: radiobuttons2.pack(side = orient_type, expand = 1, padx = 10, pady = 5)
                except Exception: null=[]
                
                try: comp1.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                try: radiobuttons1.pack(side = orient_type, expand = 1, padx = 10, pady = 5)
                except Exception: null=[]
                if reorganize == 'yes':
                    try: comp2.pack(anchor = 'w', padx = 10, pady = pady_int)
                    except Exception: null=[]
                try: self.entry_field1.pack(fill = 'x', expand = 1, padx = 10, pady = 5)
                except Exception: null=[]
                try: self.entry_field2.pack(fill = 'x', expand = 1, padx = 10, pady = 5);
                except Exception: null=[]
                try: self.entry_field3.pack(fill = 'x', expand = 1, padx = 10, pady = 5)
                except Exception: null=[]
                try: self.entry_field4.pack(fill = 'x', expand = 1, padx = 8, pady = 5)
                except Exception: null=[]
                if reorganize == 'no':
                    try: comp2.pack(anchor = 'w', padx = 10, pady = pady_int)
                    except Exception: null=[]
                try: comp3.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                try: comp4.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                try: comp5.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                try: comp6.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                try: comp7.pack(anchor = 'w', padx = 10, pady = pady_int)
                except Exception: null=[]
                enter_index=0; radio_index=0; dropdown_index=0
                if 'ge_ptype' in option or 'expression_threshold' in option or 'gene_expression_cutoff' in option:
                    custom_group = PmwFreeze.Group(self.sf.interior(),tag_text = self.group_tag)
                    custom_group.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
                    insert_into_group = 'yes'
                
            #i+=1 ####Keep track of index
        #def quitcommand(): parent.destroy; sys.exit()
        #self.button = Button(text="   Quit  ", command=quitcommand)
        #self.button.pack(side = 'bottom', padx = 10, pady = 10)

        if 'input_cdf_file' in option_list: ### For the CEL file selection window, provide a link to get Library files
            button_text = 'Download Library Files'; d_url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
            self.d_url = d_url; text_button = Button(self._parent, text=button_text, command=self.Dlinkout); text_button.pack(side = 'left', padx = 5, pady = 5)
            
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

        if 'species' in option_list or 'selected_species1' in option_list:
            new_species_button = Button(self._parent, text='Add New Species', command=self.newSpecies)
            new_species_button.pack(side = 'left', padx = 5, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()
        
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
        if os.name == 'nt':
            try: os.startfile('"'+self.pdf_help_file+'"')
            except Exception:  os.system('open "'+self.pdf_help_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.pdf_help_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.pdf_help_file+'"')   
        self._tl.destroy()
	
    def openOnlineHelp(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
        self._tl.destroy()
    
    def linkout(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
        
    def Dlinkout(self):
        try: webbrowser.open(self.d_url)
        except Exception: null=[]
            
    def setvscrollmode(self, tag):
        self.sf.configure(vscrollmode = tag)

    def info(self):
        tkMessageBox.showinfo("title","message",parent=self._parent)

    def deleteWindow(self):
        tkMessageBox.showwarning("Quit","Use 'Quit' button to end program!",parent=self._parent)

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
        if 'dir' in option or 'folder' in option:
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
                        print_out += "Contact genmapp@gladstone.ucsf.edu if this error persists with your system information.\n"
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
                except Exception: null=[]
            exportDefaultFileLocations(file_location_defaults)

        if 'file' in option:
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
                        print_out += "Contact genmapp@gladstone.ucsf.edu if this error persists with your system information.\n"
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
            exportDefaultFileLocations(file_location_defaults)
        
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
        self._user_variables[option] = tag
        change_var = ''
        if option == 'dbase_version':
            ###Export new species info
            exportDBversion(tag); change_var = 'all'
            try: self.changeVendorSelection(); self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: null=[]
        elif option == 'species':
            try: self.changeArraySelection()
            except Exception: null=[]
        elif option == 'manufacturer_selection':
            try: self.changeSpeciesSelection(); self.changeArraySelection()
            except Exception: null=[]
        #elif option == 'array_type':
            #self.checkSpeciesArraySelection(array_type)
	elif option == 'analysis_method':
	    if tag == 'ASPIRE':
		 try: self.entry_field2.setentry('0.2')
		 except Exception: null=[]
		 self._user_variables['alt_exon_fold_cutoff'] = '0.2'
	    elif tag == 'linearregres':
		try: self.entry_field2.setentry('2')
		except Exception: null=[]
		self._user_variables['alt_exon_fold_cutoff'] = '2'
        elif option == 'selected_version':
            current_species_names = db_versions[tag]
            current_species_names.sort()
            try: self.speciescomp.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp2.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp3.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated

    def changeSpeciesSelection(self):
        vendor = self._user_variables['manufacturer_selection'] ### Get vendor (stored as global)
        current_species_names = getSpeciesList(vendor) ### Need to change species, manufacturers and array_type
        try: self.speciescomp.setitems(current_species_names)
        except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'species' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                last_selected_species = self._user_variables[i]
                if last_selected_species not in current_species_names:
                    try: self._user_variables[i] = current_species_names[0]
                    except Exception: null = []

    def checkSpeciesArraySelection(self,array_type):
        current_species_names = getSpeciesForArray(array_type)
        try: self.speciescomp.setitems(current_species_names)
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
        try: self.arraycomp.setitems(current_array_types)
        except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'array_type' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                if self._user_variables[i] not in current_array_types: ### If the current array type is supported by the new species selection, keep it the same
                    try: self._user_variables[i] = current_array_types[0]
                    except Exception: null = []

    def changeVendorSelection(self):
        species_name = self._user_variables['species'] ### Get species (stored as global)
        current_array_types, manufacturer_list = getArraysAndVendors(species,'')
        try: self.vendorcomp.setitems(manufacturer_list)
        except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        for i in self._option_list:
            if 'manufacturer_selection' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                try: self._user_variables[i] = manufacturer_list[0]
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

def getSpeciesList(vendor):
    try: current_species_dirs = unique.read_directory('/AltDatabase')
    except Exception: ### Occurs when the version file gets over-written with a bad directory name
        try:
            ### Remove the version file and wipe the species file
            os.remove(filepath('Config/version.txt'))
            #raw = export.ExportFile('Config/species.txt'); raw.close()
            os.mkdir(filepath('AltDatabase'))
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        except Exception: null = []
        try: elite_db_versions = returnDirectoriesNoReplace('/AltDatabase')
        except Exception:
            try: os.mkdir(filepath('AltDatabase'))
            except Exception: null=[]
            elite_db_versions = returnDirectoriesNoReplace('/AltDatabase')
        try: exportDBversion(elite_db_versions[0])
        except Exception: exportDBversion('')
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
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    exportVersionData(db_version,today,'Config/')

def exportVersionData(version,version_date,dir):   
    new_file = dir+'version.txt'
    data = export.ExportFile(new_file)
    data.write(str(version)+'\t'+str(version_date)+'\n'); data.close()

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
        data.write(config+'\t'+str(config_db[config])+'\n'); data.close()
    
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
    fln1,status1 = update.download(base_url+'Config/species_all.txt','Config/','')
    fln2,status2 = update.download(base_url+'Config/source_data.txt','Config/','')
    fln3,status3 = update.download(base_url+'Config/array_versions.txt','Config/','')
    try:
        if 'Internet' not in status3:
            print 'Finished downloading the latest configuration files.'; root.destroy()
        else:
            try: WarningWindow(status3,'Error Encountered!'); root.destroy(); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
            except Exception: print status3; root.destroy(); sys.exit()
    except Exception: null=[]

def getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,root):
    base_url = file_location_defaults['url'].Location()
    goelite_url = file_location_defaults['goelite'].Location()
    dbs_added = 0

    AltAnalyze_folders = read_directory(''); Cytoscape_found = 'no'
    for dir in AltAnalyze_folders:
        if 'Cytoscape_' in dir: Cytoscape_found='yes'
    if Cytoscape_found == 'no':
        fln,status = update.download(goelite_url+'Cytoscape/cytoscape.zip','','')
        if 'Internet' not in status: print "Cytoscape program folder downloaded."
  
    fln,status = update.download(goelite_url+'Databases/'+db_version+'Plus/OBO.zip','AltDatabase/goelite/','')
    if 'Internet' not in status: print "Gene Ontology structure files downloaded."
    
    for species_code in new_species_codes:
        #print [base_url+'AltDatabase/'+db_version+'/'+species_code+'.zip']
        fln,status = update.download(base_url+'AltDatabase/updated/'+db_version+'/'+species_code+'.zip','AltDatabaseNoVersion/','')
        if 'Internet' not in status:
            print 'Finished downloading the latest species database files.'
            dbs_added+=1
            #print goelite_url+'Databases/'+db_version+'Plus/'+species_code+'.zip'
        try: fln,status = update.download(goelite_url+'Databases/'+db_version+'Plus/'+species_code+'.zip','AltDatabase/goelite/','')
        except Exception: print "No species GO-Elite database found."
        if 'Internet' not in status: print "GO-Elite database installed." ; dbs_added+=1
        else: print "No species GO-Elite database found."
        try: os.mkdir(filepath('AltDatabase/'+species_code))
        except Exception: null=[]
    if dbs_added>0:
        print_out = "New species data successfully added to database."
        if root !='' and root !=None:
            try: InfoWindow(print_out,'Continue')
            except Exception: print print_out
        else: print print_out
        try: root.destroy()
        except Exception: null=[]
    else:
        if root !='' and root !=None: WarningWindow(status,'Error Encountered!'); root.destroy(); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        else: print status; root.destroy(); sys.exit()
        
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

def getSystemInfo():
    importSystemInfo()
    return system_codes

def importSystemInfo():
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

    return system_list,mod_list


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
    importSpeciesInfo(); species_names={}
    for species_full in species_codes:
        sc = species_codes[species_full]; abrev = sc.SpeciesCode()
        species_names[abrev] = species_full
    return species_names
    
def importSpeciesInfo():
    try:
        if integrate_online_species == 'yes': filename = 'Config/species_all.txt'
        else: filename = 'Config/species.txt'
    except Exception: filename = 'Config/species.txt'
    
    fn=filepath(filename); global species_list; species_list=[]; global species_codes; species_codes={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            try: abrev,species,algorithms = string.split(data,'\t')
            except Exception: abrev,species = string.split(data,'\t'); algorithms = ''
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

def importArrayGroupsSimple(expr_group_dir):
    array_group_list = []; group_db={}
    fn=filepath(expr_group_dir)
    for line in open(fn,'rU').xreadlines():             
        data = cleanUpLine(line)
        array_header,group,group_name = string.split(data,'\t')
        try: group = int(group); group_db[group]=group_name
        except ValueError: print group, group_name;kill
        agd = ArrayGroupData(array_header,group,group_name)
        array_group_list.append(agd)
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
        app,status,location,species = string.split(data,'\t')
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

def exportGroups(exp_file_location_db,array_group_list):
    ### If the user supplies new defaults, over-write the existing
    for dataset_name in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]; groups_file = fl.GroupsFile()
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

def annotateMetaProbesetGenes(summary_exp_file, expression_file, metaprobeset_file, species):
    metaprobeset_cv_file = string.replace(metaprobeset_file,species+'_',species+'_Conversion_')
    metaprobeset_cv_file = string.replace(metaprobeset_cv_file,'.mps','.txt')

    fn=filepath(metaprobeset_cv_file); uid_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        uid,ens_gene = string.split(data,'\t')
        uid_db[uid] = ens_gene

    export_data = export.ExportFile(expression_file)        
    fn=filepath(summary_exp_file); x=0
    for line in open(fn,'rU').xreadlines():
        if line[0] == '#': null=[]
        elif x == 0: export_data.write(line); x+=1
        else:
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            uid = t[0]; ens_gene = uid_db[uid]
            export_data.write(string.join([ens_gene]+t[1:],'\t')+'\n')
    export_data.close()

def reformatResidualFile(residual_exp_file,residual_destination_file):
    ### Re-write the residuals file so it has a single combined unique ID (arbitrary gene ID + probe ID)
    print 'Re-formatting and moving the calculated residuals file...'
    export_data = export.ExportFile(residual_destination_file)   
    fn=filepath(residual_exp_file); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0 and data[0]=='#': null=[]
        elif x == 0:
            x+=1; t = string.split(data,'\t')
            export_data.write(string.join(['UID']+t[5:],'\t')+'\n')
        else:
            t = string.split(data,'\t')
            uid = t[0]+'-'+t[2] ### arbitrary numeric gene ID + probes ID
            export_data.write(string.join([uid]+t[5:],'\t')+'\n')
    export_data.close()
    os.remove(residual_exp_file)
    
def probesetSummarize(exp_file_location_db,analyze_metaprobesets,probeset_type,species,root):
    for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
        fl = exp_file_location_db[dataset]
        apt_dir =fl.APTLocation()
        array_type=fl.ArrayType()  
        pgf_file=fl.InputCDFFile()
        clf_file=fl.CLFFile()
        bgp_file=fl.BGPFile()
        xhyb_remove = fl.XHybRemoval()
        cel_dir=fl.CELFileDir() + '/cel_files.txt'
        expression_file = fl.ExpFile()
        stats_file = fl.StatsFile()
        output_dir = fl.OutputDir() + '/APT-output'
        cache_dir = output_dir + '/apt-probeset-summarize-cache'
        architecture = fl.Architecture() ### May over-ride the real architecture if a failure occurs
        
        if xhyb_remove == 'yes' and (array_type == 'gene' or array_type == 'junction'): xhyb_remove = 'no' ### This is set when the user mistakenly selects exon array, initially
        if analyze_metaprobesets == 'yes':
            export_features = 'true'
            metaprobeset_file = filepath('AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'_'+probeset_type+'.mps')
        
        import subprocess; import platform
        print 'Processor architecture set =',architecture,platform.machine()
        if '/bin' in apt_dir: apt_file = apt_dir +'/apt-probeset-summarize' ### if the user selects an APT directory
        elif os.name == 'nt':
            if '32bit' in architecture: apt_file = apt_dir + '/PC/32bit/apt-probeset-summarize'; plat = 'Windows'
            elif '64bit' in architecture: apt_file = apt_dir + '/PC/64bit/apt-probeset-summarize'; plat = 'Windows'
        elif 'darwin' in sys.platform: apt_file = apt_dir + '/Mac/apt-probeset-summarize'; plat = 'MacOSX'
        elif 'linux' in sys.platform:
            if '32bit' in platform.architecture(): apt_file = apt_dir + '/Linux/32bit/apt-probeset-summarize'; plat = 'linux32bit'
            elif '64bit' in platform.architecture(): apt_file = apt_dir + '/Linux/64bit/apt-probeset-summarize'; plat = 'linux64bit'
        apt_file = filepath(apt_file)
        #print 'AltAnalyze has choosen APT for',plat
        print "Beginning probeset summarization of input CEL files with Affymetrix Power Tools (APT)..."
        if 'cdf' in pgf_file or 'CDF' in pgf_file:
            if xhyb_remove == 'yes' and array_type == 'AltMouse':
                kill_list_dir = osfilepath('AltDatabase/'+species+'/AltMouse/'+species+'_probes_to_remove.txt')
            else: kill_list_dir = osfilepath('AltDatabase/affymetrix/APT/probes_to_remove.txt')
            
            try:
                cdf_file = pgf_file; algorithm = 'rma'; pval = 'dabg'
                retcode = subprocess.call([
                apt_file, "-d", cdf_file, "--kill-list", kill_list_dir, "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir]) # "-a", pval,
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)
                    try: residual_dabg = output_dir+'/dabg.residuals.txt'; os.remove(residual_dabg)
                    except Exception: null = []
            except NameError: status = 'failed'
        else:
            if xhyb_remove == 'yes':
                kill_list_dir = osfilepath('AltDatabase/'+species+'/exon/'+species+'_probes_to_remove.txt')
            else: kill_list_dir = osfilepath('AltDatabase/affymetrix/APT/probes_to_remove.txt')
            try:
                algorithm = 'rma-sketch'; pval = 'dabg'
                if analyze_metaprobesets != 'yes':
                    retcode = subprocess.call([
                    apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir,
                    "-a", algorithm, "-a", pval, "-o", output_dir, "--cel-files", cel_dir])
                else:
                    retcode = subprocess.call([
                    apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir, "-m", metaprobeset_file,
                    "-a", algorithm, "-a", pval, "-o", output_dir, "--cel-files", cel_dir, "--feature-details", export_features])
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    #if analyze_metaprobesets == 'yes': annotateMetaProbesetGenes(summary_exp_file, expression_file, metaprobeset_file, species)
                    shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)

                    summary_exp_file = output_dir+'/'+pval+'.summary.txt'
                    #if analyze_metaprobesets == 'yes': annotateMetaProbesetGenes(summary_exp_file, stats_file, metaprobeset_file, species)
                    shutil.copyfile(summary_exp_file, stats_file)
                    os.remove(summary_exp_file)
                    
                    if analyze_metaprobesets == 'yes':
                        residual_destination_file = string.replace(expression_file,'exp.','residuals.')
                        residual_exp_file = output_dir+'/'+algorithm+'.residuals.txt' 
                        #shutil.copyfile(residual_exp_file, residual_destination_file);os.remove(residual_exp_file)
                        reformatResidualFile(residual_exp_file,residual_destination_file)
                        residual_dabg_file = output_dir+'/dabg.residuals.txt'; os.remove(residual_dabg_file)
            except NameError:  status = 'failed'
            
        cache_delete_status = export.deleteFolder(cache_dir)
        if status == 'failed':
            if architecture == '64bit' and platform.architecture()[0] == '64bit' and (os.name == 'nt' or 'linux' in sys.platform):
                print 'Warning! 64bit version of APT encountered an error, trying 32bit.'
                ### If the above doesn't work, try 32bit architecture instead of 64bit (assuming the problem is related to known transient 64bit build issues)
                for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
                    fl = exp_file_location_db[dataset]; fl.setArchitecture('32bit')
                probesetSummarize(exp_file_location_db,analyze_metaprobesets,probeset_type,species,root)            
            else:
                print_out = 'apt-probeset-summarize failed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.'
                try:
                    WarningWindow(print_out,'Exit')
                    root.destroy()
                except Exception: print print_out
                sys.exit()
        else:
            print 'CEL files successfully processed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.' 
            
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
            array_abrev, dabg_p, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, avg_all_for_ss, include_raw_data, run_goelite = string.split(data,'\t')
            if array_type == array_abrev:
                return dabg_p, expression_threshold, perform_alt_analysis, analyze_as_groups, expression_data_format, avg_all_for_ss, include_raw_data, run_goelite
            
        if '-alt' in filename:
            array_abrev, analysis_method, additional_algorithms, filter_probeset_types, analyze_all_conditions, p_threshold, alt_exon_fold_variable, additional_score, permute_p_threshold, gene_expression_cutoff, perform_permutation_analysis, export_splice_index_values, run_MiDAS, calculate_splicing_index_p, filter_for_AS = string.split(data,'\t')
            if array_type == array_abrev:
                return  [analysis_method, additional_algorithms, filter_probeset_types, analyze_all_conditions, p_threshold, alt_exon_fold_variable, additional_score, permute_p_threshold, gene_expression_cutoff, perform_permutation_analysis, export_splice_index_values, run_MiDAS, calculate_splicing_index_p, filter_for_AS]
            
        if '-funct' in filename:
            array_abrev, analyze_functional_attributes, microRNA_prediction_method = string.split(data,'\t')
            if array_type == array_abrev:
                return [analyze_functional_attributes,microRNA_prediction_method]

        if '-goelite' in filename:
            array_abrev, ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, pathway_permutations, mod = string.split(data,'\t')
            if array_type == array_abrev:
                return [ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, pathway_permutations, mod]
        
class OptionData:
    def __init__(self,option,displayed_title,display_object,notes,array_options,global_default):
        self._option = option; self._displayed_title = displayed_title; self._notes = notes
        self._array_options = array_options; self._display_object = display_object; self._global_default = global_default
    def Option(self): return self._option
    def Display(self): return self._displayed_title
    def setDisplay(self,display_title): self._displayed_title = display_title
    def DisplayObject(self): return self._display_object
    def Notes(self): return self._notes
    def setNotes(self,notes): self._notes = notes
    def DefaultOption(self): return self._default_option
    def setDefaultOption(self,default_option): self._default_option = default_option
    def setNotes(self,notes): self._notes = notes
    def ArrayOptions(self): return self._array_options
    def setArrayOptions(self,array_options): self._array_options = array_options
    def __repr__(self): return self.Report()

def importUserOptions(array_type):
    filename = 'Config/options.txt'; option_db={}; option_list_db={}
    fn=filepath(filename); x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        data = string.replace(data,'\k','\n') ###Used \k in the file instead of \n, since these are removed above
        t = string.split(data,'\t')
        option,mac_displayed_title,pc_displayed_title,pc_display2,linux_displayed_title,display_object,group,notes,description,global_default = t[:10]
    
        if os.name == 'nt':
            import platform
            if '64' in platform.machine(): displayed_title = pc_display2
            elif '32' in platform.machine(): displayed_title = pc_display2
            elif '64bit' in platform.architecture(): displayed_title = pc_display2
            else: displayed_title = pc_displayed_title
        elif 'darwin' in sys.platform: displayed_title = mac_displayed_title
        elif 'linux' in sys.platform: displayed_title = linux_displayed_title
        else: displayed_title = linux_displayed_title
        
        if x == 0:
            i = t.index(array_type) ### Index position of the name of the array_type selected by user (or arbitrary to begin with)
            x = 1
        else:
            array_options = t[i]
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
        tl = Toplevel(); tl.title('AltAnalyze')
        
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
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()
        text_button = Button(parent, text=self.button_text, command=self.button1); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button2 = Button(parent, text=self.button_text2, command=self.button2); text_button2.pack(side = 'bottom', padx = 5, pady = 5)
        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def button1(self): self._user_variables['button']=self.button_text; self._parent.destroy()
    def button2(self): self._user_variables['button']=self.button_text2; self._parent.destroy()
    def ButtonSelection(self): return self._user_variables
    def deleteWindow(self):
        tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)

class IndicatorWindowSimple:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
        
class IndicatorWindow:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
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
        import ResultsExport_module
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
        except Exception: print 'Warning encountered...exiting program';sys.exit()

class InfoWindow:
    def __init__(self,dialogue,header):
        try: tkMessageBox.showinfo(header, dialogue)
        except Exception:
            print dialogue
            #print 'Attempted to open a GUI that is not accessible...exiting program';sys.exit()
            #print "Analysis finished...exiting AltAnalyze."; sys.exit()
        
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
        about = 'AltAnalyze 1.155 beta.\n'
        about+= 'AltAnalyze is an open-source, freely available application covered under the\n'
        about+= 'Apache open-source license. Additional information can be found at:\n'
        about+= 'http://www.altanalyze.org\n'
        about+= '\nDeveloped by:\n\tNathan Salomonis\n\tBruce Conklin\nGladstone Institutes 2008'
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
        txt.insert(END, 'AltAnalyze 1.155.\n')
        txt.insert(END, 'AltAnalyze is an open-source, freely available application covered under the\n')
        txt.insert(END, 'Apache open-source license. Additional information can be found at:\n')
        txt.insert(END, "http://www.altanalyze.org\n", ('link', str(0)))
        txt.insert(END, '\nDeveloped by:\n\tNathan Salomonis\n\tBruce Conklin\nGladstone Institutes 2008')
        txt.tag_config('link', foreground="blue", underline = 1)
        txt.tag_bind('link', '<Button-1>', showLink)
        
    def deleteWindow(self):
        tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)

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
        except Exception: null=[] ### This error will more likely occur since no expression file has been created 
    return files_exported

def formatArrayGroupsForGUI(array_group_list):
        ### Format input for GUI like the imported options.txt Config file, except allow for custom fields in the GUI class
        category = 'GroupArrays'; option_db={}; option_list={}; 
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
        self._comps_file = comps_file
        import platform; self.architecture = platform.architecture()[0]
    def ExpFile(self): return self._exp_file
    def StatsFile(self): return self._stats_file
    def GroupsFile(self): return self._groups_file
    def CompsFile(self): return self._comps_file
    def StatsFile(self): return self._stats_file
    def setArchitecture(self,architecture): self.architecture = architecture
    def setAPTLocation(self,apt_location): self._apt_location = osfilepath(apt_location)
    def setInputCDFFile(self,cdf_file): self._cdf_file = osfilepath(cdf_file)
    def setCLFFile(self,clf_file): self._clf_file = osfilepath(clf_file)
    def setBGPFile(self,bgp_file): self._bgp_file = osfilepath(bgp_file)
    def setCELFileDir(self,cel_file_dir): self._cel_file_dir = osfilepath(cel_file_dir)
    def setArrayType(self,array_type): self._array_type = array_type
    def setOutputDir(self,output_dir): self._output_dir = output_dir
    def setRootDir(self,parent_dir):
        ### Get directory above ExpressionInput
        split_dirs = string.split(parent_dir,'ExpressionInput')
        root_dir = split_dirs[0]
        self._root_dir = root_dir + '/'
    def setXHybRemoval(self,xhyb): self._xhyb = xhyb
    def XHybRemoval(self): return self._xhyb
    def RootDir(self): return self._root_dir
    def APTLocation(self): return self._apt_location
    def InputCDFFile(self): return self._cdf_file
    def CLFFile(self): return self._clf_file
    def BGPFile(self): return self._bgp_file
    def CELFileDir(self): return self._cel_file_dir
    def ArrayType(self): return self._array_type
    def OutputDir(self): return self._output_dir
    def Architecture(self): return self.architecture
    def Report(self): return self.ExpFile()+'|'+str(len(self.StatsFile()))+'|'+str(len(self.GroupsFile()))+'|'+str(len(self.CompsFile()))
    def __repr__(self): return self.Report()

class AdditionalAlgorithms:
    def __init__(self, additional_algorithm):
        self._additional_algorithm = additional_algorithm
    def Algorithm(self): return self._additional_algorithm
    def setScore(self,score): self._score = score
    def Score(self): return self._score
    def __repr__(self): return self.Algorithm()
    
def getUpdatedParameters(array_type,species,run_from_scratch,file_dirs):
    ### Get default options for ExpressionBuilder and AltAnalyze
    na = 'NA'; log = 'log'; no = 'no'
    global user_variables; user_variables={}; global selected_parameters; selected_parameters = []

    run_goelite=no; change_threshold=na;pathway_permutations=na;mod=na; ge_ptype = 'rawp';resources_to_analyze = na
    ge_fold_cutoffs=2;ge_pvalue_cutoffs=0.05;filter_method=na;z_threshold=1.96;p_val_threshold=0.05
    
    option_list,option_db = importUserOptions(array_type)
    global root; root = Tk()
    root.title('AltAnalyze: Perform Additional Analyses')
    selected_parameters.append('AdditionalOptions'); backSelect = 'no'
    gu = GUI(root,option_db,option_list['AdditionalOptions'],'')
    new_run = gu.Results()['new_run']
    if new_run == 'Change Parameters and Re-Run': AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
    else:
        expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = importDefaults(array_type,species)
        proceed = 'no'
        while proceed == 'no':
            root = Tk(); root.title('AltAnalyze: Pathway Analysis Parameters')
            if run_from_scratch == 'Process AltAnalyze filtered':
                option_list['GOElite'] = option_list['GOElite'][3:]; goelite_defaults = goelite_defaults[3:]
            selected_parameters.append('GOElite'); backSelect = 'no'
            gu = GUI(root,option_db,option_list['GOElite'],goelite_defaults)
            if run_from_scratch != 'Process AltAnalyze filtered':
                ge_fold_cutoffs = gu.Results()['ge_fold_cutoffs']
                ge_pvalue_cutoffs = gu.Results()['ge_pvalue_cutoffs']
                ge_ptype = gu.Results()['ge_ptype']
            filter_method = gu.Results()['filter_method']
            z_threshold = gu.Results()['z_threshold']
            p_val_threshold = gu.Results()['p_val_threshold']
            change_threshold = gu.Results()['change_threshold']
            resources_to_analyze = gu.Results()['resources_to_analyze']
            pathway_permutations = gu.Results()['pathway_permutations']
            mod = gu.Results()['mod']
            try:
                z_threshold = float(z_threshold)
                change_threshold = float(change_threshold)
                p_val_threshold = float(p_val_threshold)
                pathway_permutations = int(pathway_permutations)
                if run_from_scratch != 'Process AltAnalyze filtered':
                    ge_fold_cutoffs = float(ge_fold_cutoffs)
                    ge_pvalue_cutoffs = float(ge_pvalue_cutoffs)
                    ge_fold_cutoffs = math.log(float(ge_fold_cutoffs),2)
                proceed = 'yes'
            except Exception:
                print_out = "Invalid numerical entry. Try again."
                IndicatorWindow(print_out,'Continue')
    try:
        criterion_input_folder, criterion_denom_folder, main_output_folder = file_dirs
        import GO_Elite
        ###Export dataset criterion using user-defined filters
        ExpressionBuilder.buildCriterion(ge_fold_cutoffs, ge_pvalue_cutoffs, ge_ptype, main_output_folder)
        #except Exception: null = []; # print 'No expression files to summarize'
        goelite_var = species,mod,pathway_permutations,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,file_dirs,''
        GO_Elite.remoteAnalysis(goelite_var,'UI')
        AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
    except OSError:
        print_out = "Unexpected error encountered. Please see log file."
        IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
   
def addOnlineSpeciesDatabases(backSelect):
    StatusWindow(file_location_defaults,'getOnlineDBConfig')
    #except Exception,e: print [e]; null = []

    importSystemInfo(); exportSystemInfo() ### By re-importing we incorporate new source data from the downloaded file    
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
        else: gu = PreviousResults(old_options)
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
        if species_full == 'Add Species':
            AltAnalyze.AltAnalyzeSetup(species_full); sys.exit()
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
                            if ac.Manufacturer() in ad.Manufacturer() and 'expression' in ac.ArrayName():
                                if sc not in compatible_species: compatible_species.append(sc)
                            ac.setSpeciesCodes(compatible_species)
                exportArrayInfo(array_codes)
        if len(new_species_codes) > 0:
            analysis = 'getOnlineEliteDatabase'
            values = file_location_defaults,db_version,new_species_codes
            StatusWindow(values,analysis)
            proceed = 'yes'
        else:
            print_out = "Please select a species before continuing."
            IndicatorWindow(print_out,'Try Again')
    
    #db_versions_vendors
    exportSpeciesInfo(existing_species_codes)
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

def getUserParameters(run_parameter):
    global AltAnalyze; import AltAnalyze
    if run_parameter == 'yes':
        try:
            try: MainMenu()
            except TclError: null=[]
        except NameError: print "Tcl is not availalbe on this machine. Please install Tkinter for Python before using AltAnalyze with the GUI. Running AltAnalyze via command line options are also availalble without the GUI."; sys.exit()
    global species; species=''; global user_variables; user_variables={}; global analysis_method; global array_type
    global PathDir; global PathFile; global file_location_defaults; global integrate_online_species; integrate_online_species = 'no'
    global option_db; global option_list; global analysis_status; analysis_status = 'continue'; global selected_parameters; selected_parameters=[]
    global backSelect
    
    ### Get default options for ExpressionBuilder and AltAnalyze

    na = 'NA'; log = 'log'; no = 'no'
    run_from_scratch=na; expression_threshold=na; perform_alt_analysis=na; expression_data_format=log
    include_raw_data=na; avg_all_for_ss=no; dabg_p=na;
    analysis_method=na; p_threshold=na; filter_probeset_types=na; alt_exon_fold_cutoff=na
    permute_p_threshold=na; perform_permutation_analysis=na; export_splice_index_values=no
    run_MiDAS=no; analyze_functional_attributes=no; microRNA_prediction_method=na
    gene_expression_cutoff=na; cel_file_dir=na; input_exp_file=na; input_stats_file=na; filter_for_AS=no
    calculate_splicing_index_p=no; run_goelite=no; ge_ptype = 'rawp'
    ge_fold_cutoffs=2;ge_pvalue_cutoffs=0.05;filter_method=na;z_threshold=1.96;p_val_threshold=0.05
    change_threshold=2;pathway_permutations=na;mod=na; analyze_all_conditions=no; resources_to_analyze=na
    additional_algorithms = na
                
    option_list,option_db = importUserOptions('exon')  ##Initially used to just get the info for species and array_type
    importSpeciesInfo()
    file_location_defaults = importDefaultFileLocations()
    importArrayInfo()

    try: elite_db_versions = returnDirectoriesNoReplace('/AltDatabase')
    except Exception:
        try: elite_db_versions=[]; os.mkdir(filepath('AltDatabase'))
        except Exception: null=[] ### directory already exists      
    try: gene_database_dir = unique.getCurrentGeneDatabaseVersion()
    except Exception: gene_database_dir=''
    if len(elite_db_versions)>0 and gene_database_dir == '':
        gene_database_dir = elite_db_versions[-1]; exportDBversion(elite_db_versions[-1])
    
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
                except Exception: null=[]
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
        if run_parameter == 'Add Species': species_full = 'Homo sapiens'; species = 'Hs'; vendor = 'Affymetrix'; specific_array = 'Exon ST array'
        if backSelect == 'yes' and 'array_type' in old_options: null=[]
        elif 'Homo sapiens' in current_species_names: species_full = 'Homo sapiens'; species = 'Hs'; vendor = 'Affymetrix'; specific_array = 'Exon ST array'
        elif 'Mus musculus' in current_species_names: species_full = 'Mus musculus'; species = 'Mm'; vendor = 'Affymetrix'; specific_array = 'Exon ST array'
        elif 'Rattus norvegicus' in current_species_names: species_full = 'Rattus norvegicus'; species = 'Rn'; vendor = 'Affymetrix'; specific_array = 'Exon ST array'
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
                root = Tk(); root.title('AltAnalyze: Select Species and Microarray Type')
                gu = GUI(root,option_db,option_list['ArrayType'],'')
            else: gu = PreviousResults(old_options)
            species_full = gu.Results()['species']

        if species_full == 'Add Species' or 'NewSpecies' == selected_parameters[-1]:
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
                    exportSpeciesInfo(species_codes)
                    try: os.mkdir(filepath('AltDatabase/'+new_species_code))
                    except Exception: null=[]
                    for array_system in array_codes:
                        ac = array_codes[array_system]
                        manufacturer=ac.Manufacturer()
                        compatible_species = ac.SpeciesCodes()
                        if manufacturer == new_manufacturer and 'expression array' in ac.ArrayName():
                            if new_species_code not in compatible_species: compatible_species.append(new_species_code)
                        ac.setSpeciesCodes(compatible_species)
                        exportArrayInfo(array_codes)
                    fn = filepath('AltDatabase/affymetrix/'+ new_species_code)
                    try: os.mkdir(fn)
                    except OSError: null = [] ### Directory already exists
                    AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
                else:
                    print_out = "Valid species data was not added. You must\nindicate a two letter species code and full species name."
                    IndicatorWindow(print_out,'Continue')  
        else: species = species_codes[species_full].SpeciesCode()    
        update_dbs = gu.Results()['update_dbs']
        array_full = gu.Results()['array_type']
        array_type = array_codes[array_full].ArrayCode()
        vendor = gu.Results()['manufacturer_selection']
        if update_dbs == 'yes':
            integrate_online_species = 'yes'
            addOnlineSpeciesDatabases(backSelect)
            AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        new_analysis_options=[]

        if array_type == 'gene':
            try: gene_database = unique.getCurrentGeneDatabaseVersion()
            except Exception: gene_database='00'
            if int(gene_database[-2:]) < 54:
                print_out = 'The AltAnalyze database indicated for Gene 1.0 ST\narray analysis is not supported for alternative exon\nanalysis. Please update to EnsMart54 or greater\nbefore proceeding.'
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
        if array_type == 'junction':
            try: gene_database = unique.getCurrentGeneDatabaseVersion()
            except Exception: gene_database='00'
            if int(gene_database[-2:]) < 55:
                print_out = 'The AltAnalyze database indicated for the JAY array\n is not supported for alternative exon analysis.\nPlease update to EnsMart55 or greater before\nproceeding.'
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
            downloaded_junction_db = 'no'; file_problem='no'
            while downloaded_junction_db == 'no': ### Used as validation in case internet connection is unavailable
                try: dirs = read_directory('/AltDatabase/'+species)
                except Exception: dirs=[]
                if 'junction' not in dirs or file_problem == 'yes':
                    if file_problem == 'yes':
                        print_out = 'Unknown junction installation error occured.\nPlease try again.'
                    else:
                        print_out = 'To perform a junction array analysis AltAnalyze must\nfirst download the appropriate junction array database.'
                    IndicatorWindow(print_out,'Download'); filename = 'AltDatabase/'+species+'/'+species+'_junction.zip'
                    dir = 'AltDatabase/updated/'+gene_database; var_list = filename,dir
                    if debug_mode == 'no': StatusWindow(var_list,'download')
                try: dirs = read_directory('/AltDatabase/'+species)
                except Exception: dirs=[]
                if 'junction' in dirs:
                    file_length = AltAnalyze.verifyFileLength('AltDatabase/'+species+'/junction/'+species+'_Ensembl_probesets.txt')
                    if file_length>0: downloaded_junction_db = 'yes'
                    else: file_problem = 'yes'
        if array_type == "3'array":
            for i in option_db['run_from_scratch'].ArrayOptions():
                if 'AltAnalyze' not in i:
                    if array_type == "3'array":
                        if 'CEL' in i and vendor != 'Affymetrix': proceed = 'no'
                        else: proceed = 'yes'
                    else: proceed = 'yes'
                    if proceed == 'yes': new_analysis_options.append(i)
            option_db['run_from_scratch'].setArrayOptions(new_analysis_options)

        if array_type == 'exon' or array_type == 'AltMouse' or array_type == 'gene':
            try: dirs = read_directory('/AltDatabase/'+species)
            except Exception: dirs=[]
            if len(dirs)==0:
                print_out = 'Valid database directories were not found for this array.\nPlease re-install database.'
                IndicatorWindow(print_out,'Continue'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()

        proceed = 'no'
        if len(new_analysis_options)!=1:
            if backSelect == 'no' or 'AnalysisType' == selected_parameters[-1]:
                selected_parameters.append('AnalysisType'); backSelect = 'no'
                root = Tk(); root.title('AltAnalyze: Select Analysis Method')
                gu = GUI(root,option_db,option_list['AnalysisType'],'')
            else: gu = PreviousResults(old_options)
            run_from_scratch = gu.Results()['run_from_scratch']
        else: run_from_scratch = 'Process Expression file'
        vendor = array_codes[array_full].Manufacturer()
        constitutive_source = array_codes[array_full].ConstitutiveSource()
        option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type

        if backSelect == 'yes':
            for option in old_options: ### Set options to user selected
                try: option_db[option].setDefaultOption(old_options[option])
                except Exception: null=[]
                
        if run_from_scratch == 'Process CEL files':
            """Designate CEL file directory, Dataset Name and Output Directory"""
            assinged = 'no'
            while assinged == 'no': ### Assigned indicates whether or not the CEL directory and CDF files are defined
                if species == 'Rn': del option_list['InputCELFiles'][-1] ### Don't examine xyb
                #print (((backSelect,selected_parameters)))
                if backSelect == 'no' or 'InputCELFiles' == selected_parameters[-1]:
                    selected_parameters.append('InputCELFiles'); backSelect = 'no'
                    root = Tk(); root.title('AltAnalyze: Select CEL files for APT')
                    gu = GUI(root,option_db,option_list['InputCELFiles'],'')
                else: gu = PreviousResults(old_options)
                dataset_name = gu.Results()['dataset_name']
                try: remove_xhyb = gu.Results()['remove_xhyb']
                except KeyError: remove_xhyb = 'no'
                if len(dataset_name)<1:
                    print_out = "Please provide a name for the dataset before proceeding."
                    IndicatorWindow(print_out,'Continue')
                elif 'input_cel_dir' in gu.Results():
                    cel_file_dir = gu.Results()['input_cel_dir']
                    cel_files,cel_files_fn=identifyCELfiles(cel_file_dir)
                    try: output_dir = gu.Results()['output_CEL_dir']
                    except KeyError: output_dir = cel_file_dir
                    if len(output_dir)==0: output_dir = cel_file_dir
                    if len(cel_files)>0: assinged = 'yes' ### CEL files are present in this directory
                    else:
                        print_out = "No valid .CEL files were found in the directory\n"+cel_file_dir+"\nPlease verify and try again."
                        IndicatorWindow(print_out,'Continue')
                else:
                    print_out = "The directory containing CEL files has not\nbeen assigned! Select a directory before proceeding."
                    IndicatorWindow(print_out,'Continue')
            cel_file_list_dir = exportCELFileList(cel_files_fn,cel_file_dir)

            """Determine if Library and Annotations for the array exist, if not, download or prompt for selection"""
            try: specific_array_types,specific_array_type = identifyArrayType(cel_files_fn); num_array_types = len(specific_array_types)
            except Exception: null=[]; num_array_types=1; specific_array_type = None
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
                elif array_type == 'junction':
                    if species == 'Hs': specific_array_type = 'HJAY_v2'
                    if species == 'Mm': specific_array_type = 'MJAY_v2'

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
            if (len(input_cdf_file)==0 or len(annotation_dir) == 0) and backSelect == 'no':
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
                                    print destination_parent
                                    print input_cdf_file
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
                                    if array_type == 'exon' or array_type == 'junction': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
                                    else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
                                    dir_list = read_directory(parent_dir)
                                    if clf_short in dir_list and bgp_short in dir_list:
                                        pgf_file = input_cdf_file
                                        clf_file = string.replace(pgf_file,'.pgf','.clf')
                                        if array_type == 'exon' or array_type == 'junction': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                                        else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                                        assinged = 'yes'
                                        ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                        destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                                        print destination_parent
                                        print input_cdf_file
                                        if destination_parent not in input_cdf_file:
                                            info_list = input_cdf_file,osfilepath(destination_parent+cdf_short); StatusWindow(info_list,'copy')
                                            info_list = clf_file,osfilepath(destination_parent+clf_short); StatusWindow(info_list,'copy')
                                            info_list = bgp_file,osfilepath(destination_parent+bgp_short); StatusWindow(info_list,'copy')
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
                                print destination_parent
                                print input_cdf_file
                                if destination_parent not in input_cdf_file:
                                    info_list = input_annotation_file,filepath(destination_parent+csv_short); StatusWindow(info_list,'copy')
                                sd = SupprotedArrays(specific_array_type,cdf_short,csv_short,species,array_type)
                                supproted_array_db[specific_array_type] = sd
                                try: exportSupportedArrayInfo()
                                except Exception: continue ### Occurs if the file is open... not critical to worry about       
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

        if run_from_scratch != 'update DBs': ### Update DBs is an option which has been removed from 1.1. Should be a separate menu item soon.
            expr_defaults, alt_exon_defaults, functional_analysis_defaults, goelite_defaults = importDefaults(array_type,species)
            
            if run_from_scratch != 'Process AltAnalyze filtered' and run_from_scratch != 'Annotate External Results':
                proceed = 'no'
                while proceed == 'no':
                    if backSelect == 'no' or 'GeneExpression' == selected_parameters[-1]:
                        selected_parameters.append('GeneExpression'); backSelect = 'no'
                        root = Tk(); root.title('AltAnalyze: Expression Analysis Parameters')
                        gu = GUI(root,option_db,option_list['GeneExpression'],expr_defaults)
                    else: gu = PreviousResults(old_options)
                    if array_type != "3'array":          
                        dabg_p = gu.Results()['dabg_p']
                        run_from_scratch = gu.Results()['run_from_scratch']
                        expression_threshold = gu.Results()['expression_threshold']
                        perform_alt_analysis = gu.Results()['perform_alt_analysis']
                        try: analyze_as_groups = gu.Results()['analyze_as_groups']
                        except Exception: analyze_as_groups = ''
                        if perform_alt_analysis == 'just expression': perform_alt_analysis = 'expression'
                        else: perform_alt_analysis = 'both'
                        try: avg_all_for_ss = gu.Results()['avg_all_for_ss']
                        except Exception: avg_all_for_ss = 'no'
                        if 'all exon aligning' in avg_all_for_ss: avg_all_for_ss = 'yes'
                        else: avg_all_for_ss = 'no'
                    expression_data_format = gu.Results()['expression_data_format']
                    include_raw_data = gu.Results()['include_raw_data']
                    run_goelite = gu.Results()['run_goelite']
                    if 'immediately' in run_goelite: run_goelite = 'yes'
                    else: run_goelite = 'no'
                    passed = 'yes'; print_out = 'Invalid threshold entered for '
                    if array_type != "3'array":    
                        try: expression_threshold = float(expression_threshold)
                        except Exception: passed = 'no'; print_out+= 'expression threshold'
                        try: dabg_p = float(dabg_p)
                        except Exception: passed = 'no'; print_out+= 'DABG p-value cutoff'
                        if expression_threshold<1: passed = 'no'; print_out+= 'expression threshold'
                        elif dabg_p<=0 or dabg_p>1: passed = 'no'; print_out+= 'DABG p-value cutoff'
                    if passed == 'no': IndicatorWindow(print_out,'Continue')
                    else: proceed = 'yes'
                    
            if (perform_alt_analysis == 'both') or (run_from_scratch == 'Process AltAnalyze filtered') or (run_from_scratch == 'Annotate External Results'):
                perform_alt_analysis = 'alt'

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
                        except Exception: null=[]
                    
                if run_from_scratch == 'Process AltAnalyze filtered':
                    functional_analysis_defaults.append('constitutive probesets'); option_list['AltAnalyze'].append('avg_all_for_ss')
                    if run_goelite == 'no':
                        functional_analysis_defaults.append('decide later'); option_list['AltAnalyze'].append('run_goelite')

                if run_from_scratch == 'Annotate External Results':
                    ### Remove options relating to expression analysis when importing filtered probeset lists
                    options_to_exclude = ['analysis_method','p_threshold','gene_expression_cutoff','alt_exon_fold_cutoff','run_MiDAS']
                    options_to_exclude+= ['export_splice_index_values','run_goelite','analyze_all_conditions','calculate_splicing_index_p']
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
                    try: filter_probeset_types = gu.Results()['filter_probeset_types']
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
                        except Exception: null=[]
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
                    try:
                        avg_all_for_ss = gu.Results()['avg_all_for_ss']
                        if 'all exon aligning' in avg_all_for_ss: avg_all_for_ss = 'yes'
                        else: avg_all_for_ss = 'no'
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
                if run_from_scratch == 'Process AltAnalyze filtered':
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
                p_val_threshold = gu.Results()['p_val_threshold']
                change_threshold = gu.Results()['change_threshold']
                resources_to_analyze = gu.Results()['resources_to_analyze']
                pathway_permutations = gu.Results()['pathway_permutations']
                mod = gu.Results()['mod']
                ge_fold_cutoffs = math.log(float(ge_fold_cutoffs),2)
    except OSError:
        null=[]; sys.exit()
    """In this next section, create a set of GUI windows NOT defined by the options.txt file.
    These are the groups and comps files"""
    original_comp_group_list=[]; array_group_list=[]; group_name_list=[]
    if run_from_scratch != 'Process AltAnalyze filtered' and run_from_scratch != 'Annotate External Results': ### Groups and Comps already defined

        if run_from_scratch == 'Process CEL files':
            if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name+'.txt'
            groups_name = string.replace(dataset_name,'exp.','groups.')
            comps_name = string.replace(dataset_name,'exp.','comps.')
            if "ExpressionInput" not in output_dir:
                output_dir = output_dir + '/ExpressionInput' ### Store the result files here so that files don't get mixed up
                try: os.mkdir(output_dir) ### Since this directory doesn't exist we have to make it
                except OSError: null = [] ### Directory already exists
            exp_file_dir = output_dir+'/'+dataset_name
            ### store file locations (also use these later when running APT)
            stats_file_dir = string.replace(exp_file_dir,'exp.','stats.')
            groups_file_dir = string.replace(exp_file_dir,'exp.','groups.')
            comps_file_dir = string.replace(exp_file_dir,'exp.','comps.')
            fl = ExpressionFileLocationData(exp_file_dir,stats_file_dir,groups_file_dir,comps_file_dir)
            exp_file_location_db={}; exp_file_location_db[dataset_name]=fl
            parent_dir = output_dir  ### interchangable terms (parent_dir used with expression file import)
            
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
                groups_file_dir = parent_dir+'/'+groups_name; comps_file_dir = parent_dir+'/'+comps_name
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
        if backSelect == 'yes':
            for cel_file in cel_files:
                if cel_file in user_variables: group_name = user_variables[cel_file]; group = ''
                else: group = ''; group_name = ''    
                agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd)      
        elif groups_name in dir_files:
            array_group_list,group_db = importArrayGroupsSimple(groups_file_dir) #agd = ArrayGroupData(array_header,group,group_name)
            if comps_name in dir_files:
                comp_group_list, null = ExpressionBuilder.importComparisonGroups(comps_file_dir)
                for group1,group2 in comp_group_list:
                    try: group_name1 = group_db[int(group1)]; group_name2 = group_db[int(group2)]
                    except KeyError:
                        print_out = 'The "comps." file for this dataset has group numbers\nnot listed in the "groups." file.'
                        WarningWindow(print_out,'Exit'); AltAnalyze.AltAnalyzeSetup('no'); sys.exit()
                    original_comp_group_list.append((group_name1,group_name2)) ### If comparisons already exist, default to these
        else:
            for cel_file in cel_files:
                group = ''; group_name = ''    
                agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd)
                
        if len(array_group_list)>0: ### Thus we are not analyzing the default (ExpressionInput) directory of expression, group and comp data.
            option_db,option_list = formatArrayGroupsForGUI(array_group_list)
            ###Force this GUI to repeat until the user fills in each entry, but record what they did add
            user_variables_long={}
            while len(user_variables_long) != len(option_db):
                if backSelect == 'no' or 'GroupArrays' == selected_parameters[-1]:
                    selected_parameters.append('GroupArrays'); backSelect = 'no'
                    root = Tk(); root.title('AltAnalyze: Assign CEL files to a Group Annotation'); user_variables_long={}
                    #import copy; user_variables_original = copy.deepcopy(user_variables); user_variables={}
                    gu = GUI(root,option_db,option_list['GroupArrays'],'groups')
                else: gu = PreviousResults(old_options)
                for option in user_variables: ### By default, all arrays will be assigned a group of ''
                    try:
                        if len(user_variables[option])>0:
                            if option in option_db: user_variables_long[option]=[]
                    except Exception: null=[]
                ###Store the group names and assign group numbers
                group_name_db={}; group_name_list = []; group_number = 1
                for cel_file in option_list['GroupArrays']: ### start we these CEL files, since they are ordered according to their order in the expression dataset
                    group_name = gu.Results()[cel_file]
                    if group_name not in group_name_db:
                        group_name_db[group_name]=group_number; group_number+=1
                        group_name_list.append(group_name)
                if len(group_name_db)==2: analyze_all_conditions = 'pairwise' ### Don't allow multiple comparison analysis if only two conditions present
                
                ###Store the group names and numbers with each array_id in memory   
                for agd in array_group_list:
                    cel_file = agd.Array()
                    group_name = gu.Results()[cel_file] ###Lookup the new group assignment entered by the user
                    group_number = group_name_db[group_name]
                    agd.setGroupName(group_name); agd.setGroup(group_number)
                if (len(user_variables_long) != len(option_db)) or len(group_name_db)<2:
                    if len(group_name_db)<2:
                        print_out = "At least two array groups must be established\nbefore proceeding."
                    else:
                        print_out = "Not all arrays have been assigned a group. Please\nassign to a group before proceeding (required)."
                    IndicatorWindow(print_out,'Continue')   
                    option_db,option_list = formatArrayGroupsForGUI(array_group_list) ### array_group_list at this point will be updated with any changes made in the GUI by the user

            i=2; px=0 ###Determine the number of possible comparisons based on the number of groups
            while i<=len(group_name_list): px = px + i - 1; i+=1
            group_name_list.reverse(); group_name_list.append(''); group_name_list.reverse() ### add a null entry first
            if px > 100: px = 100 ### With very large datasets, AltAnalyze stalls
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
                        if len(group_name)>0: ### Group_names are by default blank
                            cn_main,cn_minor = string.split(comparison[11:],'-') ### e.g. 1-1 and 1-2
                            try:
                                null = int(cn_main); null = int(cn_minor)
                                try: comp_groups_db[cn_main].append([cn_minor,group_name])
                                except KeyError: comp_groups_db[cn_main]=[[cn_minor,group_name]]
                            except Exception: null=[]
                    except Exception: null=[]
                print_out = "You must pick at least one comparison group before proceeding."
                if len(comp_groups_db)>0:
                    try:
                        comp_group_list=[]
                        for cn_main in comp_groups_db:
                            cg = comp_groups_db[cn_main]; cg.sort(); comp_group_list.append([cn_main,[group_name_db[cg[0][1]],group_name_db[cg[1][1]]]])
                            if cg[0][1] == cg[1][1]: identical_groups = 'yes' ### Thus the two groups in the comparisons are identical, flag
                        comp_group_list.sort()                    
                        proceed = 'yes'
                    except Exception: print_out = "You must pick at least two groups for each comparison."
                if identical_groups == 'yes': proceed = 'no'; print_out = "The same group is listed as both the experimental and\ncontrol group in a comparison. Fix before proceeding."
                if proceed == 'no': IndicatorWindow(print_out,'Continue')   
                
            ### Export user modified groups and comps files
            exported = 0
            while exported == 0:
                try:
                    fl = exp_file_location_db[dataset_name]; groups_file = fl.GroupsFile()
                    exportGroups(exp_file_location_db,array_group_list)
                    exported = 1
                except Exception:                 
                    print_out = "The file:\n"+groups_file+"\nis still open. This file must be closed before proceeding"
                    IndicatorWindow(print_out,'Continue')                 
            exported = 0
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
    apt_location = getAPTLocations(file_location_defaults,run_from_scratch,run_MiDAS)

    ### Set the primary parent directory for ExpressionBuilder and AltAnalyze (one level above the ExpressionInput directory, if present)
    for dataset in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]
        fl.setAPTLocation(apt_location)
        if run_from_scratch == 'Process CEL files':
            fl.setInputCDFFile(input_cdf_file); fl.setCLFFile(clf_file); fl.setBGPFile(bgp_file); fl.setXHybRemoval(remove_xhyb)
            fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type); fl.setOutputDir(output_dir)
        fl = exp_file_location_db[dataset]; fl.setRootDir(parent_dir)

    expr_var = species,array_type,vendor,constitutive_source,dabg_p,expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch, perform_alt_analysis
    alt_var = analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_cutoff,gene_expression_cutoff,permute_p_threshold, perform_permutation_analysis, export_splice_index_values, analyze_all_conditions
    additional_var = calculate_splicing_index_p, run_MiDAS, analyze_functional_attributes, microRNA_prediction_method, filter_for_AS, additional_algorithms
    goelite_var = ge_fold_cutoffs,ge_pvalue_cutoffs,ge_ptype,filter_method,z_threshold,p_val_threshold,change_threshold,resources_to_analyze,pathway_permutations,mod

    return expr_var, alt_var, additional_var, goelite_var, exp_file_location_db

def getAPTLocations(file_location_defaults,run_from_scratch,run_MiDAS):
    import ResultsExport_module
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
            exportDefaultFileLocations(file_location_defaults)
    return apt_location
    
if __name__ == '__main__':
    #getUpdatedParameters(array_type,species,run_from_scratch,file_dirs)
    
    a = getUserParameters('yes'); print a; sys.exit()
    
