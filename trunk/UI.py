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
import export
import ExpressionBuilder
import time
import webbrowser
from sys import argv
try:
    import Tkinter
    from Tkinter import *
    import PmwFreeze
    from Tkconstants import LEFT
    import tkMessageBox
    import tkFileDialog
except ImportError: print "\nPmw or Tkinter not found... proceeding with manual input"

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

def identifyCELfiles(dir):
    dir_list = read_directory(dir); dir_list2=[]; full_dir_list=[]
    for file in dir_list:
        file_lower = string.lower(file)
        if '.cel' in file_lower and '.cel.' not in file_lower:
            dir_list2.append(file)
            file = dir+'/'+file
            full_dir_list.append(file)
    return dir_list2,full_dir_list

def identifyArrayType(full_dir_list):
    arrays={} ### Determine the type of unique arrays in each directory
    for filename in full_dir_list:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():             
            data = cleanUpLine(line)
            if 'DatHeader' == data[:9]:                
                array_info,null = string.split(data,'sq')
                array_info = string.split(array_info,' ')
                array_type = array_info[-1]
                if '.' in array_type: array_type,null = string.split(array_type,'.')
                arrays[array_type]=[]
                break
    array_list = []
    for array in arrays: array_list.append(array)
    return array_list, array_type

def getAffyFiles(array_name,species):#('AltDatabase/affymetrix/LibraryFiles/'+library_file,species)
    sa = supproted_array_db[array_name]; library_file = sa.LibraryFile(); annot_file = sa.AnnotationFile(); original_library_file = library_file
    filename = 'AltDatabase/affymetrix/LibraryFiles/'+library_file
    fn=filepath(filename); library_dir=filename; bgp_file = ''; clf_file = ''
    import update; reload(update); warn = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
        input_cdf_file = filename
        if '.pgf' in input_cdf_file:
            ###Check to see if the clf and bgp files are present in this directory 
            icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
            clf_short = string.replace(cdf_short,'.pgf','.clf')
            if array_type == 'exon': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
            else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
            try: dir_list = read_directory(parent_dir)
            except Exception: dir_list = read_directory('/'+parent_dir)
            if clf_short in dir_list and bgp_short in dir_list:
                pgf_file = input_cdf_file; clf_file = string.replace(pgf_file,'.pgf','.clf')
                if array_type == 'exon': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
            else:
                print_out = "The directory;\n"+parent_dir+"\ndoes not contain either a .clf or antigenomic.bgp\nfile, required for probeset summarization."
                IndicatorWindow(print_out,'Continue')
    except Exception:
        print_out = "AltAnalyze was not able to find a library file\nfor your arrays. Would you like AltAnalyze to\nautomatically download these files?"
        dw = DownloadWindow(print_out,'Download','Continue'); warn = 'no'
        dw_results = dw.Results(); option = dw_results['selected_option']
        if option == 1:
            library_file = string.replace(library_file,'.cdf','.zip')
            filename = 'AltDatabase/affymetrix/LibraryFiles/'+library_file
            input_cdf_file = filename
            if '.pgf' in input_cdf_file:
                pgf_file = input_cdf_file; clf_file = string.replace(pgf_file,'.pgf','.clf')
                if array_type == 'exon': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                filenames = [pgf_file+'.gz',clf_file+'.gz',bgp_file+'.gz']
            else: filenames = [input_cdf_file]
            for filename in filenames:
                var_list = filename,'LibraryFiles'
                if debug_mode == 'no': StatusWindow(var_list,'download')
                else:
                    for filename in filenames:
                        update.downloadCurrentVersionUI(filename,'LibraryFiles','',Tk())
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
            dw = DownloadWindow(print_out,'Download','Continue'); warn = 'no'
            dw_results = dw.Results(); option = dw_results['selected_option']
        if option == 1:
            annot_file += '.zip'
            filenames = ['AltDatabase/affymetrix/'+species+'/'+annot_file]
            for filename in filenames:
                var_list = filename,'AnnotationFiles'
                if debug_mode == 'no': StatusWindow(var_list,'download')
                else:
                    for filename in filenames:
                        update.downloadCurrentVersionUI(filename,'AnnotationFiles','',Tk())
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
    print 'to...',file2
    shutil.copyfile(file1,file2)
    root.destroy()
    
class StatusWindow:
    def __init__(self,info_list,analysis_type):
            root = Tk()
            self._parent = root
            root.title('AltAnalyze 1.11 Beta')
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
            import update; reload(update)
            status = StringVarFile(statusVar,root) ### Likely captures the stdout
            if analysis_type == 'download':
                filename,dir = info_list
                sys.stdout = status; root.after(100,update.downloadCurrentVersionUI(filename,dir,'',self._parent))
            if analysis_type == 'copy':
                file1,file2 = info_list
                sys.stdout = status; root.after(100,copyFiles(file1,file2,self._parent))
            self._parent.mainloop()

    def deleteWindow(self): tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()
    
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

class GUI:
    def __init__(self, parent, option_db, option_list, defaults): 
        self._parent = parent; self._option_list = option_list; self._option_db = option_db
        self._user_variables = user_variables; i = 0
        
        filename = 'Config/icon.gif'
        if 'input_cel_dir' in option_list: filename = 'Config/aa_0.gif'
        if 'include_raw_data' in option_list: filename = 'Config/aa_1.gif'
        if 'filter_for_AS' in option_list: filename = 'Config/aa_2.gif'
        
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
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
                label_text_str = 'AltAnalyze Group Comparisons'
                if len(option_list)<5: height = 250; width = 400
            elif 'filter_for_AS' in option_list:
                label_text_str = 'AltAnalyze Alternative Exon Analysis Parameters'
                height = 420; width = 775; use_scroll = 'yes'
            else:
                label_text_str = "AltAnalyze Main Dataset Parameters"
                height = 300; width = 400; use_scroll = 'yes'
            self.sf = PmwFreeze.ScrolledFrame(self._parent,
                    labelpos = 'n', label_text = label_text_str,
                    usehullsize = 1, hull_width = width, hull_height = height)
            self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
            self.frame = self.sf.interior()
            if defaults == 'comps':
                Label(self.frame,text=notes).pack()
        for option in option_list:
            od = option_db[option]; self.title = od.Display(); notes = od.Notes()      
            self.display_options = od.ArrayOptions()
            #if len(defaults)>0: print i, od.DisplayObject(),self.display_options, option_list, defaults
            if option == 'array_type':
                valid_display_options=[]
                for array_name in self.display_options:
                    compatible_species = array_codes[array_name].SpeciesCodes()
                    if species in compatible_species: valid_display_options.append(array_name)
                self.display_options = valid_display_options
            if 'radio' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                ### Create and pack a RadioSelect widget, with radiobuttons.
                self._option = option
                def radiocallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                radiobuttons = PmwFreeze.RadioSelect(parent_type,                       
                        buttontype = 'radiobutton', orient = 'vertical',
                        labelpos = 'w', command = radiocallback, label_text = self.title,
                        hull_borderwidth = 2, hull_relief = 'ridge',
                ); radiobuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                ### print self.display_options
                ### Add some buttons to the radiobutton RadioSelect.
                for text in self.display_options:
                    if text != ['NA']: radiobuttons.add(text)
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                radiobuttons.invoke(self.default_option)
                if len(notes)>0: Label(self._parent, text=notes).pack()
            if 'button' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                
                self._option = option
                if mac_print_mode == 'yes': button_type = 'radiobutton'
                else: button_type = 'button'                
                ### Create and pack a horizontal RadioSelect widget.
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                def buttoncallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                horiz = PmwFreeze.RadioSelect(parent_type, buttontype = button_type, orient = 'horizontal',
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
                self.default_dir = ''                
                entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
                if len(od.GlobalDefault())>0: entrytxt.set(od.GlobalDefault())
                self.pathdb[option] = entrytxt ###This should normally occur, unless setGlobalDefault has been set
                #l = Label(group.interior(), text=self.title); l.pack(side=LEFT)        
                entry = Entry(group.interior(),textvariable=self.pathdb[option]); entry.pack(side=LEFT,fill = 'both', expand = 1, padx = 10, pady = 2)
                button = Button(group.interior(), text="select "+od.DisplayObject(), width = 10, fg="red", command=filecallback); button.pack(side=LEFT, padx = 2,pady = 2)                    
             
                if len(notes)>0: ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)

            if 'drop-down' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = self.display_options
                ### Pack these into a groups to maintain organization
                #group = PmwFreeze.Group(self._parent,tag_text = self.title)
                #group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                    
                def comp_callback1(tag,callback=self.callback,option=option):
                    callback(tag,option)

                self.comp = PmwFreeze.OptionMenu(self._parent,
                    labelpos = 'w', label_text = self.title,                                                 
                    items = self.default_option, command = comp_callback1,
                ); self.comp.pack(anchor = 'w', padx = 10, pady = 0)
                self.comp.invoke(self.default_option[0]) ###Just pick the first option
                
            if 'pulldown_comps' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = self.display_options
                ###From the option, create two new options, one for each group in the comparison
                option1 = option+'-1'; option2 = option+'-2'
                ### Pack these into a groups to maintain organization
                group = PmwFreeze.Group(self.sf.interior(),tag_text = self.title)
                group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                    
                def comp_callback1(tag,callback=self.callback,option1=option1):
                    callback(tag,option1)
                def comp_callback2(tag,callback=self.callback,option2=option2):
                    callback(tag,option2)

                #labelpos = 'w', label_text = self.title,  -inside of OptionMenu
                self.comp1 = PmwFreeze.OptionMenu(group.interior(),
                    items = self.default_option, menubutton_width = 20, command = comp_callback1,
                ); self.comp1.pack(side = LEFT, anchor = 'w', padx = 10, pady = 0)

                self.comp2 = PmwFreeze.OptionMenu (group.interior(), 
                    items = self.default_option, menubutton_width = 20, command = comp_callback2,
                ); self.comp2.pack(side = LEFT, anchor = 'w', padx = 10, pady = 0)

                self.comp1.invoke(notes[0]) ; self.comp2.invoke(notes[1])
                
            if 'simple_entry' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                self.default_option = self.display_options[0]
                def enter_callback(tag,enter_callback=self.enter_callback,option=option):
                    enter_callback(tag,option)
                self.title = self.title + '\t '
                entry_field = PmwFreeze.EntryField(self.sf.interior(),
                        labelpos = 'w', label_text = self.title,
                        validate = enter_callback,
                        value = self.default_option
                ); entry_field.pack(padx = 10, pady = 1)

            if 'enter' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                #print self.default_option, self.title; kill

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
                    self.default_option = ''; use_method = 'i'
                if use_method == 'p':
                    entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate_p, 
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                    ); entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                if use_method == 'i':
                    entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate,
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                    ); entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                if len(notes)>0: Label(self._parent, text=notes).pack()
            if 'multiple-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                self._option = option
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=option):
                    checkbuttoncallback(tag,state,option)                    
                ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                        buttontype = 'checkbutton', orient = 'vertical',
                        labelpos = 'w', command = self.checkbuttoncallback,
                        label_text = self.title, hull_borderwidth = 2, hull_relief = 'ridge',
                ); self.checkbuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                ### Add some buttons to the checkbutton RadioSelect.
                for text in self.display_options:
                     if text != ['NA']: self.checkbuttons.add(text)
                self.checkbuttons.invoke(self.default_option)
                self.checkbuttons.invoke(self.default_option2)
                if len(notes)>0: Label(self._parent, text=notes).pack()
                
            if 'single-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                if use_scroll == 'yes': parent_type = self.sf.interior()
                else: parent_type = self._parent
                self._option = option
                proceed = 'yes'
                """if option == 'export_splice_index_values':
                    if analysis_method != 'splicing-index': proceed = 'no' ### only export corrected constitutive ratios if splicing index method chosen"""
                if proceed == 'yes':
                    if len(defaults) <1: self.default_option = self.display_options[0]
                    else: self.default_option = defaults[i]
                    if self.default_option != 'NA':
                        def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=option):
                            checkbuttoncallback(tag,state,option)                
                        ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                        self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                                buttontype = 'checkbutton', command = checkbuttoncallback,
                                hull_borderwidth = 2, hull_relief = 'ridge',
                        ); self.checkbuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                        ### Add some buttons to the checkbutton RadioSelect.
                        self.checkbuttons.add(self.title)
                        if self.default_option == 'yes': self.checkbuttons.invoke(self.title)
                        else: self._user_variables[option] = 'no'
                if len(notes)>0: Label(self._parent, text=notes).pack()
            i+=1 ####Keep track of index
            
        #def quitcommand(): parent.destroy; sys.exit()
        #self.button = Button(text="   Quit  ", command=quitcommand)
        #self.button.pack(side = 'bottom', padx = 10, pady = 10)

        if 'input_cdf_file' in option_list: ### For the CEL file selection window, provide a link to get Library files
            button_text = 'Download Library Files'; url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
            self.url = url; text_button = Button(self._parent, text=button_text, command=self.linkout); text_button.pack(side = 'left', padx = 5, pady = 5)
            
        continue_to_next_win = Button(text = 'Continue', command = self._parent.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)

        back_button = Button(self._parent, text="Back", command=self.goBack) 
        back_button.pack(side = 'right', padx =10, pady = 5)
        
        quit_win = Button(self._parent, text="Quit", command=self.quit) 
        quit_win.pack(side = 'right', padx =10, pady = 5)

        button_text = 'Help'; url = 'ReadMe/help_main.htm'; url = filepath(url)
        self.url = url; help_button = Button(self._parent, text=button_text, command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()

    def goBack(self):
        self._parent.destroy()
        getUserParameters('no'); sys.exit()
        
    def linkout(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
            
    def setvscrollmode(self, tag):
        self.sf.configure(vscrollmode = tag)

    def info(self):
        tkMessageBox.showinfo("title","message",parent=self._parent)

    def deleteWindow(self):
        tkMessageBox.showwarning("Quit","Use 'Quit' button to end program!",parent=self._parent)

    def quit(self):
        self._parent.quit()
        self._parent.destroy()
        sys.exit()

    def continue_win(self):
        ### Currently not used - can be used to check data integrity before closing a window
        self._parent.quit()
        self._parent.destroy()
        sys.exit()
        
    def chooseDirectory(self,option):
        tag = tkFileDialog.askdirectory(parent=self._parent)
        self._user_variables[option] = tag

    def chooseFile(self,option):
        tag = tkFileDialog.askopenfile(parent=self._parent)
        self._user_variables[option] = tag.name

    def getPath(self,option):
        if 'dir' in option:
            dirPath = tkFileDialog.askdirectory(parent=self._parent) #initialdir=self.default_dir
        if 'file' in option:
            tag = tkFileDialog.askopenfile(parent=self._parent)
            try: dirPath = tag.name #initialdir=self.default_dir
            except AttributeError: dirPath = ''
        entrytxt = self.pathdb[option]
        entrytxt.set(dirPath)
        self._user_variables[option] = dirPath
        ### Allows the option_db to be updated for the next round (if an error is encountered)
        od = self._option_db[option]
        od.setGlobalDefault(dirPath)
        
    def Report(self,tag,option):
        output = tag
        return output
    def __repr__(self,tag,option): return self.Report(tag,option)
    
    def Results(self): return self._user_variables

    def custom_validate(self, text, option):
        self._user_variables[option] = text
        try: text = float(text);return 1
        except ValueError: return -1

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

 ################# Import Options and Defaults ##################

def getPrimaryUserParameters():
    print "\n******Analysis Options*******"
    species_code=''; array_code=''
    while species_code == '':
        x = 1; print "Select Species for Analyses"
        for species_name in species_list: print str(x)+')',species_name; x+=1
        try:
            inp = sys.stdin.readline(); inp = int(inp.strip())
            try: species = species_list[inp-1]; species_code = species_codes[species].SpeciesCode()
            except IndexError: print "!!!!!Requires a numeric value between 1-"+str(len(species_list)),'\n'
        except ValueError: print "!!!!!Requires a numeric value between 1-"+str(len(species_list)),'\n'
    
    while array_code == '':
        x = 1; array_list2=[]; print "\nSelect Array Type"
        for array_name in array_list:
            if species_code in array_codes[array_name].SpeciesCodes():
                print str(x)+')',array_name; x+=1; array_list2.append(array_name)
        try:
            inp = sys.stdin.readline(); inp = int(inp.strip())
            try: array = array_list2[inp-1]; array_code = array_codes[array].ArrayCode()
            except IndexError: print "!!!!!Requires a numeric value between 1-"+str(len(species_list)),'\n'
        except ValueError: print "!!!!!Requires a numeric value between 1-"+str(len(species_list)),'\n'
    if "3'" in array: ###Thus the user is analyzing a conventional 3' array
        print '\nMake sure that an Affymetrix .CSV annotation file for your microarray has been to "AltDatabases/Affymetrix/*species-code*/"'
        print "(Select enter/return to continue...)"; inp = sys.stdin.readline()

    manufacturer = array_codes[array].Manufacturer()
    constitutive_source = array_codes[array].ConstitutiveSource()
    
    print "\nSelect Analysis"
    print '1) Run analyses from scratch using expression file files ("ExpressionInput")'
    print '2) Directly run alternative exon analyses using AltAnalyze filtered files ("AltExpression")'
    print '3) Update existing AltAnalyze databases'
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": run_from_scratch = 'expression file'
    elif inp == "2": run_from_scratch = 'AltAnalyze filtered'
    elif inp == "3": run_from_scratch = 'update DBs'
    
    print "Proceed Using Default Parameters?"
    print '1) Use defaults (see "Config/default-**.txt")'
    print '2) Change defaults (see "Config/options.txt" and associated default files)'
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": proceed = 'yes'
    elif inp == "2": proceed = 'no'; sys.exit()
    
    return species_code, array_code, manufacturer, constitutive_source, run_from_scratch

class SupprotedArrays:
    def __init__(self, array_name, library_file, annotation_file, species, array_type):
        self.array_name = array_name; self.library_file = library_file; self.annotation_file = annotation_file
        self.species = species; self.array_type = array_type
    def ArrayName(self): return self.array_name
    def LibraryFile(self): return self.library_file
    def AnnotationFile(self): return self.annotation_file
    def Species(self): return self.species
    def ArrayType(self): return self.array_type
    def __repr__(self): return self.Report()
    
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
            
class SpeciesData:
    def __init__(self, abrev, species, algorithms):
        self._abrev = abrev; self._species = species; self._algorithms = algorithms
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def Algorithms(self): return self._algorithms
    def __repr__(self): return self.Report()
    
def importSpeciesInfo():
    filename = 'Config/species.txt'; x=0
    fn=filepath(filename); global species_list; species_list=[]; global species_codes; species_codes={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        abrev,species,algorithms = string.split(data,'\t')
        if x==0: x=1
        else:
            algorithms = string.split(algorithms,'|')
            species_list.append(species)
            sd = SpeciesData(abrev,species,algorithms)
            species_codes[species] = sd
    ### Add an other entry
    species_list.append('other')
    sd = SpeciesData('other','other',[])
    species_codes['other'] = sd

def exportSpeciesInfo(species_codes):
    fn=filepath('Config/species.txt'); data = open(fn,'w'); x=0
    header = string.join(['species_code','species_name','compatible_algorithms'],'\t')+'\n'
    data.write(header)
    for species in species_codes:
        if 'other' not in species:
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
    def __repr__(self): return self.Report()
    
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
            
class DefaultFileLocationData:
    def __init__(self, status, location, species):
        self._status = status; self._location = location; self._species = species
    def Status(self): return self._status
    def Location(self): return self._location
    def ResetLocation(self,location): self._location = location
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
        fl = DefaultFileLocationData(status, location, species)
        if x==0: x=1
        else:
            try: file_location_defaults[app].append(fl)
            except KeyError: file_location_defaults[app] = [fl]
    return file_location_defaults

def exportDefaultFileLocations(file_location_defaults):
    ### If the user supplies new defaults, over-write the existing
    fn=filepath('Config/default-files.csv'); data = open(fn,'w')
    data.write('Program/Download,Status,Location,Species\n')
    for app in file_location_defaults:
        fl_list = file_location_defaults[app]
        for fl in fl_list:
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

def probesetSummarize(exp_file_location_db,species,root):
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

        import subprocess; import platform
        if '/bin' in apt_dir: apt_file = apt_dir +'/apt-probeset-summarize' ### if the user selects an APT directory
        elif os.name == 'nt': apt_file = apt_dir + '/PC/apt-probeset-summarize'
        elif 'darwin' in sys.platform: apt_file = apt_dir + '/Mac/apt-probeset-summarize'
        elif 'linux' in sys.platform:
            if '32bit' in platform.architecture(): apt_file = apt_dir + '/Linux/32bit/apt-probeset-summarize'
            elif '64bit' in platform.architecture(): apt_file = apt_dir + '/Linux/64bit/apt-probeset-summarize'
        apt_file = filepath(apt_file)
        print "Begining probeset summarization of input CEL files with Affymetrix Power Tools (APT)..."
        if array_type == "3'array":
            try:
                cdf_file = pgf_file; algorithm = 'rma'
                retcode = subprocess.call([
                apt_file, "-d", cdf_file, "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir])
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)
            except NameError: status = 'failed'
            
        if array_type == 'gene':
            try:
                algorithm = 'rma-sketch'
                retcode = subprocess.call([
                apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file,
                "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir])
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)                    
            except NameError: status = 'failed'
            #apt-probeset-summarize -a rma-sketch -a plier-mm-sketch -p chip.pgf -c chip.clf -o output-dir *.cel

        if array_type == 'exon':
            if xhyb_remove == 'yes':
                kill_list_dir = osfilepath('AltDatabase/'+species+'/exon/'+species+'_probes_to_remove.txt')
            else: kill_list_dir = osfilepath('AltDatabase/affymetrix/APT/probes_to_remove.txt')
            try:
                algorithm = 'rma-sketch'; pval = 'dabg'
                retcode = subprocess.call([
                apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir,
                "-a", algorithm, "-a", pval, "-o", output_dir, "--cel-files", cel_dir])
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)

                    summary_exp_file = output_dir+'/'+pval+'.summary.txt'
                    shutil.copyfile(summary_exp_file, stats_file)
                    os.remove(summary_exp_file) 
            except NameError:  status = 'failed'
            
        cache_delete_status = export.deleteFolder(cache_dir)
        if status == 'failed':
            print_out = 'apt-probeset-summarize failed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.'
            WarningWindow(print_out,'Exit')
            root.destroy(); sys.exit()
        else:
            print 'CEL files successfully processed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.' 
            
def importDefaults(array_type,species):
    filename = 'Config/defaults-expr.txt'
    expr_defaults = importDefaultInfo(filename,array_type)
    #perform_alt_analysis, expression_data_format, dabg_p, expression_threshold, avg_all_for_ss, include_raw_data
    
    filename = 'Config/defaults-alt_exon.txt'
    alt_exon_defaults = importDefaultInfo(filename,array_type)
    #analysis_method, alt_exon_fold_variable, p_threshold, filter_probeset_types, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold,exportTransitResultsforAnalysis, export_splice_index_values = values
    
    filename = 'Config/defaults-funct.txt'
    functional_analysis_defaults = importDefaultInfo(filename,array_type)
    #analyze_functional_attributes,microRNA_prediction_method = functional_analysis_defaults
    return expr_defaults, alt_exon_defaults, functional_analysis_defaults

def importDefaultInfo(filename,array_type):
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if '-expr' in filename:
            array_abrev, dabg_p, expression_threshold, perform_alt_analysis, expression_data_format, avg_all_for_ss, include_raw_data = string.split(data,'\t')
            if array_type == array_abrev:
                return dabg_p, expression_threshold, perform_alt_analysis, expression_data_format, avg_all_for_ss, include_raw_data
            
        if '-alt' in filename:
            array_abrev, analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_variable, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, MiDAS_analysis, export_splice_index_values, calculate_splicing_index_p, filter_for_AS = string.split(data,'\t')
            if array_type == array_abrev:
                ### NOTE: p_threshold is used for MiDAS and t-test comparison p-values thresholds
                if MiDAS_analysis == 'yes': exportTransitResultsforAnalysis = 'yes'
                else: exportTransitResultsforAnalysis = 'no'
                return [analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_variable, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, exportTransitResultsforAnalysis, export_splice_index_values, calculate_splicing_index_p, filter_for_AS]
            
        if '-funct' in filename:
            array_abrev, analyze_functional_attributes, microRNA_prediction_method = string.split(data,'\t')
            if array_type == array_abrev:
                return [analyze_functional_attributes,microRNA_prediction_method]

class OptionData:
    def __init__(self,option,displayed_title,display_object,notes,array_options,global_default):
        self._option = option; self._displayed_title = displayed_title; self._notes = notes
        self._array_options = array_options; self._display_object = display_object; self._global_default = global_default
    def Option(self): return self._option
    def Display(self): return self._displayed_title
    def setDisplay(self,display_title): self._displayed_title = display_title
    def DisplayObject(self): return self._display_object
    def Notes(self): return self._notes
    def GlobalDefault(self): return self._global_default
    def setGlobalDefault(self,global_default): self._global_default = global_default
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
        option,displayed_title,display_object,group,notes,description,global_default = t[:7]
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
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()

class DownloadWindow:
    def __init__(self,message,option1,option2):
        self._user_variables = {}
        self.message = message; self.option1 = option1; self.option2 = option2
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        text_button = Button(parent, text=self.option1, command=self.selected1); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button2 = Button(parent, text=self.option2, command=self.selected2); text_button2.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
    def selected1(self):
        self._user_variables['selected_option']=1; self._parent.destroy()
    def selected2(self):
        self._user_variables['selected_option']=2; self._parent.destroy()
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
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()      
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
    def destroy_win(self): self._parent.quit(); self._parent.destroy()  
    def Folder(self): return self._tag
    
class WarningWindow:
    def __init__(self,warning,window_name):
        tkMessageBox.showerror(window_name, warning)

class InfoWindow:
    def __init__(self,dialogue,header):
        tkMessageBox.showinfo(header, dialogue)

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
                label_text = 'AltAnalyze version 1.11 Main', frame_borderwidth = 2,
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
        about = 'AltAnalyze 1.11 beta.\n'
        about+= 'AltAnalyze is an open-source, freely available application covered under the\n'
        about+= 'Apache open-source license. Additional information can be found at:\n'
        about+= 'http://www.genmapp.org/AltAnalyze\n'
        about+= '\nDeveloped by:\n\tNathan Salomonis\n\tBruce Conklin\nGladstone Institutes 2008'
        tkMessageBox.showinfo("About AltAnalyze",about,parent=self._parent)
        """
        
        def showLink(event):
            idx= int(event.widget.tag_names(CURRENT)[1])
            webbrowser.open(LINKS[idx])
        LINKS=('http://www.genmapp.org/AltAnalyze','')
        self.LINKS = LINKS
        tl = Toplevel() ### Create a top-level window separate than the parent        
        txt=Text(tl)

        #filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        #can = Canvas(tl); can.pack(side='left'); can.config(width=img.width(), height=img.height())        
        #can.create_image(2, 2, image=img, anchor=NW)
        
        txt.pack(expand=True, fill="both")
        txt.insert(END, 'AltAnalyze 1.11 beta.\n')
        txt.insert(END, 'AltAnalyze is an open-source, freely available application covered under the\n')
        txt.insert(END, 'Apache open-source license. Additional information can be found at:\n')
        txt.insert(END, "http://www.genmapp.org/AltAnalyze\n", ('link', str(0)))
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
            txt.insert(END, "http://www.genmapp.org/AltAnalyze\n", ('link', str(0)))
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

def formatArrayGroupsForGUI(array_group_list):
        ### Format input for GUI like the imported options.txt Config file, except allow for custom fields in the GUI class
        category = 'GroupArrays'; option_db={}; option_list={}; user_variables={}
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
    def ExpFile(self): return self._exp_file
    def StatsFile(self): return self._stats_file
    def GroupsFile(self): return self._groups_file
    def CompsFile(self): return self._comps_file
    def StatsFile(self): return self._stats_file
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
    def Report(self): return self.ExpFile()+len(self.StatsFile())+len(self.GroupsFile())+len(self.CompsFile())
    def __repr__(self): return self.Report()

def getUserParameters(skip_intro):
    if skip_intro == 'yes':
        try: MainMenu()
        except TclError: null=[]
    global species; species=''; global user_variables; user_variables={}; global analysis_method; global array_type
    global analysis_status; analysis_status = 'continue'
    ### Get default options for ExpressionBuilder and AltAnalyze

    na = 'NA'; log = 'log'; no = 'no'
    run_from_scratch=na; expression_threshold=na; perform_alt_analysis=na; expression_data_format=log
    include_raw_data=na; avg_all_for_ss=na; dabg_p=na;
    analysis_method=na; p_threshold=na; filter_probeset_types=na; alt_exon_fold_cutoff=na
    permute_p_threshold=na; perform_permutation_analysis=na; export_splice_index_values=no
    exportTransitResultsforAnalysis=no; analyze_functional_attributes=no; microRNA_prediction_method=na
    gene_expression_cutoff='any'; cel_file_dir=na; input_exp_file=na; input_stats_file=na; filter_for_AS=no
    calculate_splicing_index_p=no
            
    option_list,option_db = importUserOptions('exon')  ##Initially used to just get the info for species and array_type
    importSpeciesInfo()
    file_location_defaults = importDefaultFileLocations()
    array_list = importArrayInfo()
    try:   
        ###Update this informatin in option_db which will be over-written after the user selects a species and array_type
        option_db['species'].setArrayOptions(species_list)

        global root; root = Tk()
        root.title('AltAnalyze: Species Selection')
        gu = GUI(root,option_db,option_list['Species'],'')
        species_full = gu.Results()['species']
        species = species_codes[species_full].SpeciesCode()

        if species == 'other':
            species_added = 'no'
            while species_added == 'no':
                root = Tk()
                root.title('AltAnalyze: Add New Species Support')
                gu = GUI(root,option_db,option_list['NewSpecies'],'')
                new_species_code = gu.Results()['new_species_code']
                new_species_name = gu.Results()['new_species_name']
                allowed_array_systems = gu.Results()['allowed_array_systems']
                if len(new_species_code)==2 and len(new_species_name)>0 and len(allowed_array_systems)>0:
                    species_added = 'yes'
                    sd = SpeciesData(new_species_code,new_species_name,[''])
                    species_codes[new_species_name] = sd
                    exportSpeciesInfo(species_codes)
                    ac = array_codes[allowed_array_systems]; compatible_species = ac.SpeciesCodes()
                    if new_species_code not in compatible_species: compatible_species.append(new_species_code)
                    ac.setSpeciesCodes(compatible_species)
                    exportArrayInfo(array_codes)
                    fn = filepath('AltDatabase/affymetrix/'+ new_species_code)
                    try: os.mkdir(fn)
                    except OSError: null = [] ### Directory already exists
                    species = new_species_code; species_full = new_species_name
                else:
                    print_out = "Valid species data was not added. You must\nindicate a two letter species code and full species name."
                    IndicatorWindow(print_out,'Continue')  
        array_list2=[]
        for array_name in array_list:
            if species in array_codes[array_name].SpeciesCodes(): array_list2.append(array_name)

        array_list = array_list2 ### Filtered based on compatible species arrays
        option_db['array_type'].setArrayOptions(array_list)

        proceed = 'no'
        while proceed == 'no':        
            root = Tk()
            root.title('AltAnalyze: Select Microarray and Analysis Method')
            gu = GUI(root,option_db,option_list['AnalysisType'],'')
            array_full = gu.Results()['array_type']
            array_type = array_codes[array_full].ArrayCode()
            run_from_scratch = gu.Results()['run_from_scratch']
            if array_type == "3'array" or array_type == 'gene':
                if run_from_scratch == 'CEL files' or run_from_scratch == 'expression file':
                    proceed = 'yes' ### Only certain options are compatible with expression arrays
                else:
                    print_out = "The selected option is not compatible with gene expression only arrays."
                    IndicatorWindow(print_out,'Continue')   
            else: proceed = 'yes' ### All options are compatible with splicing arrays

        array_type_original = array_type
        array_type = string.replace(array_type,'gene',"3'array")
        manufacturer = array_codes[array_full].Manufacturer()
        constitutive_source = array_codes[array_full].ConstitutiveSource()
        option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type
        
        if run_from_scratch == 'CEL files':
            """Designate CEL file directory, Dataset Name and Output Directory"""
            assinged = 'no'
            while assinged == 'no': ### Assigned indicates whether or not the CEL directory and CDF files are defined
                root = Tk()
                root.title('AltAnalyze: Select CEL files for APT')
                if species == 'Rn': del option_list['InputCELFiles'][-1] ### Don't examine xyb
                gu = GUI(root,option_db,option_list['InputCELFiles'],'')
                dataset_name = gu.Results()['dataset_name']
                try: remove_xhyb = gu.Results()['remove_xhyb']
                except KeyError: remove_xhyb = 'no'
                if len(dataset_name)<1:
                    print_out = "Please provide a name for the dataset before proceeding."
                    IndicatorWindow(print_out,'Continue')
                elif 'input_cel_dir' in gu.Results():
                    cel_file_dir = gu.Results()['input_cel_dir']
                    cel_files,cel_files_fn=identifyCELfiles(cel_file_dir)
                    try: output_dir = gu.Results()['output_dir']
                    except KeyError: output_dir = cel_file_dir
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
            except Exception: null=[]
            importSupportedArrayInfo()
            try:
                sa = supproted_array_db[specific_array_type]; array_species = sa.Species(); cel_array_type = sa.ArrayType()
            except KeyError: library_dir=''; array_species=''; annotation_dir=''; cel_array_type=''
            ### Check for issues with arrays or user input options
            if num_array_types>1: ### More than one array type found in the directory
                print_out = 'Warning!!!!!!!\n\nMultiple array_types found ("'+specific_array_types[0]+'" and "'+specific_array_types[1]+'").\nIt is recommended you restart, otherwise, APT will try\n to process all different array types together as "'+specific_array_types[-1]+'".'
                IndicatorWindow(print_out,'Continue with Existing')
            if array_species != species and len(array_species)>0:
                print_out = "The CEL files indicate that the proper\nspecies is "+array_species+", however, you\nindicated "+species+ ". The species indicated by the CEL\nfiles will be used instead."
                IndicatorWindow(print_out,'Continue')
                species = array_species
            if cel_array_type != array_type and len(cel_array_type)>0:
                print_out = "The CEL files indicate that the proper\narray type is "+cel_array_type+", however, you\nindicated "+array_type+ ". The array type indicated by the CEL\nfiles will be used instead."
                IndicatorWindow(print_out,'Continue')
                array_type = cel_array_type
            ### See if the library and annotation files are on the server or are local
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
            if len(input_cdf_file)==0 or len(annotation_dir) == 0:
                assinged = 'no'
                while assinged == 'no': ### Assigned indicates whether or not the CEL directory and CDF files are defined
                    if array_type == "3'array":
                        op = option_db['input_cdf_file']; input_cdf_file_label = op.Display()
                        op.setNotes('   note: the CDF file is apart of the standard library files for this array.   ')
                        input_cdf_file_label = string.replace(input_cdf_file_label,'PGF','CDF')
                        op.setDisplay(input_cdf_file_label)      
                    root = Tk()
                    root.title('AltAnalyze: Select Affymetrix Library and Annotation files')
                    gu = GUI(root,option_db,option_list_library,'')
                    if 'input_cdf_file' in option_list_library: ### Deals with Annotation Files
                        if 'input_cdf_file' in gu.Results():
                            input_cdf_file = gu.Results()['input_cdf_file']; input_cdf_file_lower = string.lower(input_cdf_file)
                            if array_type == "3'array":
                                if '.cdf' in input_cdf_file_lower:
                                    clf_file='';bgp_file=''; assinged = 'yes'
                                    ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                    icf_list = string.split(input_cdf_file,'/'); cdf_short = icf_list[-1]
                                    destination_parent = 'AltDatabase/affymetrix/LibraryFiles/'
                                    info_list = input_cdf_file,filepath(destination_parent+cdf_short); StatusWindow(info_list,'copy')
                                else:
                                    print_out = "The file;\n"+input_cdf_file+"\ndoes not appear to be a valid Affymetix\nlibrary file. If you do not have library files, you must\ngo to the Affymetrix website to download."
                                    IndicatorWindow(print_out,'Continue')            
                            else:
                                if '.pgf' in input_cdf_file_lower:
                                    ###Check to see if the clf and bgp files are present in this directory 
                                    icf_list = string.split(input_cdf_file,'/'); parent_dir = string.join(icf_list[:-1],'/'); cdf_short = icf_list[-1]
                                    clf_short = string.replace(cdf_short,'.pgf','.clf')
                                    if array_type == 'exon': bgp_short = string.replace(cdf_short,'.pgf','.antigenomic.bgp')
                                    else: bgp_short = string.replace(cdf_short,'.pgf','.bgp')
                                    dir_list = read_directory(parent_dir)
                                    if clf_short in dir_list and bgp_short in dir_list:
                                        pgf_file = input_cdf_file
                                        clf_file = string.replace(pgf_file,'.pgf','.clf')
                                        if array_type == 'exon': bgp_file = string.replace(pgf_file,'.pgf','.antigenomic.bgp')
                                        else: bgp_file = string.replace(pgf_file,'.pgf','.bgp')
                                        assinged = 'yes'
                                        ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                        destination_parent = 'AltDatabase/affymetrix/Library/'
                                        info_list = input_cdf_file,filepath(destination_parent+cdf_short); StatusWindow(info_list,'copy')
                                        info_list = clf_file,filepath(destination_parent+clf_short); StatusWindow(info_list,'copy')
                                        info_list = bgp_file,filepath(destination_parent+bgp_short); StatusWindow(info_list,'copy')
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
                        if 'input_annotation_file' in gu.Results():
                            input_annotation_file = gu.Results()['input_annotation_file']; input_annotation_lower = string.lower(input_annotation_file)
                            if '.csv' in input_annotation_lower:
                                assinged = 'yes'
                                ###Thus the CDF or PDF file was confirmed, so copy it over to AltDatabase
                                icf_list = string.split(input_annotation_file,'/'); csv_short = icf_list[-1]
                                destination_parent = 'AltDatabase/affymetrix/'+species+'/'
                                info_list = input_annotation_file,filepath(destination_parent+csv_short); StatusWindow(info_list,'copy')
                                sd = SupprotedArrays(specific_array_type,cdf_short,csv_short,species,array_type)
                                supproted_array_db[specific_array_type] = sd
                                try: exportSupportedArrayInfo()
                                except Exception: continue ### Occurs if the file is open... not critical to worry about
                                
        if run_from_scratch == 'expression file':
            status = 'repeat'
            while status == 'repeat':
                root = Tk()
                root.title('AltAnalyze: Select Expression File for Filtering')
                gu = GUI(root,option_db,option_list['InputExpFiles'],'')
                try: input_exp_file = gu.Results()['input_exp_file']
                except KeyError: input_exp_file = '' ### Leave this blank so that the default directory is used
                try: input_stats_file = gu.Results()['input_stats_file']
                except KeyError: input_stats_file = '' ### Leave this blank so that the default directory is used
                #if array_type == 'exon':
                try: output_dir = gu.Results()['output_dir']
                except KeyError: output_dir = '' ### Leave this blank so that the default directory is used
                cel_files, array_linker_db = ExpressionBuilder.getArrayHeaders(input_exp_file)
                if len(cel_files)>0: status = 'continue'
                else:
                    print_out = "The expression file:\n"+input_exp_file+"\ndoes not appear to be a valid expression file. Check to see that\nthis is the correct tab-delimited text file."
                    IndicatorWindow(print_out,'Continue')
        if run_from_scratch != 'update DBs': ### Update DBs is an option which has been removed from 1.1. Should be a separate menu item soon.
            expr_defaults, alt_exon_defaults, functional_analysis_defaults = importDefaults(array_type,species)
            
            if run_from_scratch != 'AltAnalyze filtered':
                root = Tk(); root.title('AltAnalyze: Expression Analysis Parameters')
                gu = GUI(root,option_db,option_list['GeneExpression'],expr_defaults)
                if array_type != "3'array":          
                    dabg_p = gu.Results()['dabg_p']
                    run_from_scratch = gu.Results()['run_from_scratch']
                    expression_threshold = gu.Results()['expression_threshold']
                    perform_alt_analysis = gu.Results()['perform_alt_analysis']
                    avg_all_for_ss = gu.Results()['avg_all_for_ss']
                expression_data_format = gu.Results()['expression_data_format']
                include_raw_data = gu.Results()['include_raw_data']
                
            if (perform_alt_analysis == 'both') or (run_from_scratch == 'AltAnalyze filtered'):
                perform_alt_analysis = 'alt'

                if run_from_scratch == 'AltAnalyze filtered':
                    input_filtered_dir = ''
                    while len(input_filtered_dir)<1:
                        root = Tk(); root.title('AltAnalyze: Select AltAnalyze Filtered Probe set Files')
                        gu = GUI(root,option_db,option_list['InputFilteredFiles'],'')
                        if 'input_filtered_dir' in gu.Results():
                            input_filtered_dir = gu.Results()['input_filtered_dir']
                        else: 
                            print_out = "The directory containing filtered probe set text files has not\nbeen assigned! Select a valid directory before proceeding."
                            IndicatorWindow(print_out,'Continue')
                    fl = ExpressionFileLocationData('','','',''); dataset_name = 'filtered-exp_dir'
                    dirs = string.split(input_filtered_dir,'AltExpression'); parent_dir = dirs[0]
                    exp_file_location_db={}; exp_file_location_db[dataset_name]=fl
                #print option_list[i:i+len(alt_exon_defaults)+len(functional_analysis_defaults)], alt_exon_defaults+functional_analysis_defaults;kill
                option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type
                root = Tk(); root.title('AltAnalyze: Alternative Exon Analysis Parameters')
                gu = GUI(root,option_db,option_list['AltAnalyze'],alt_exon_defaults+functional_analysis_defaults); user_variables = {}

                analysis_method = gu.Results()['analysis_method']
                p_threshold = gu.Results()['p_threshold']
                gene_expression_cutoff = gu.Results()['gene_expression_cutoff']
                filter_probeset_types = gu.Results()['filter_probeset_types']
                alt_exon_fold_cutoff = gu.Results()['alt_exon_fold_cutoff']
                try: permute_p_threshold = gu.Results()['permute_p_threshold']
                except KeyError: permute_p_threshold = 0.05 ### Doesn't matter, not used
                try: perform_permutation_analysis = gu.Results()['perform_permutation_analysis']
                except KeyError: perform_permutation_analysis = perform_permutation_analysis
                try: export_splice_index_values = gu.Results()['export_splice_index_values']
                except KeyError: export_splice_index_values = export_splice_index_values
                try: exportTransitResultsforAnalysis = gu.Results()['exportTransitResultsforAnalysis']
                except KeyError: exportTransitResultsforAnalysis = exportTransitResultsforAnalysis
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
                        IndicatorWindow(print_out,'Continue'); getUserParameters('no'); sys.exit()

    except IndexError:  ###Occurs when PMW or Tkinter does not propperly load
        species, array_type, manufacturer, constitutive_source, run_from_scratch = getPrimaryUserParameters()
        expr_defaults, alt_exon_defaults, functional_analysis_defaults = importDefaults(array_type,species)

        dabg_p, expression_threshold, perform_alt_analysis, expression_data_format, avg_all_for_ss, include_raw_data = expr_defaults
        analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_cutoff, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, exportTransitResultsforAnalysis, export_splice_index_values = alt_exon_defaults
        analyze_functional_attributes,microRNA_prediction_method = functional_analysis_defaults
        
    """In this next section, create a set of GUI windows NOT defined by the options.txt file.
    These are the groups and comps files"""
    original_comp_group_list=[]; array_group_list=[]; group_name_list=[]
    if run_from_scratch != 'AltAnalyze filtered': ### Groups and Comps already defined

        if run_from_scratch == 'CEL files':
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
            
        if run_from_scratch == 'expression file':
            if len(input_exp_file)>0:
                if len(input_stats_file)>1: ###Make sure the files have the same arrays and order first
                    cel_files2, array_linker_db2 = ExpressionBuilder.getArrayHeaders(input_stats_file)
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
        if groups_name in dir_files and comps_name in dir_files:
            array_group_list,group_db = importArrayGroupsSimple(groups_file_dir) #agd = ArrayGroupData(array_header,group,group_name)
            comp_group_list, null = ExpressionBuilder.importComparisonGroups(comps_file_dir)
            for group1,group2 in comp_group_list:
                try: group_name1 = group_db[int(group1)]; group_name2 = group_db[int(group2)]
                except KeyError:
                    print_out = 'The "comps." file for this dataset has group numbers\nnot listed in the "groups." file.'
                    WarningWindow(print_out,'Exit'); getUserParameters('no'); sys.exit()
                original_comp_group_list.append((group_name1,group_name2)) ### If comparisons already exist, default to these
        else:
            array_group_list=[]
            for cel_file in cel_files:
                group = ''; group_name = '' ### temporarily assign group and group_name as blanks, until assigned
                agd = ArrayGroupData(cel_file,group,group_name); array_group_list.append(agd)
                        
        if len(array_group_list)>0: ### Thus we are not analyzing the default (ExpressionInput) directory of expression, group and comp data.
            option_db,option_list = formatArrayGroupsForGUI(array_group_list)
            ###Force this GUI to repeat until the user fills in each entry, but record what they did add
            user_variables_long={}
            while len(user_variables_long) != len(option_db):
                root = Tk(); root.title('AltAnalyze: Assign CEL files to a Group Annotation'); user_variables={}; user_variables_long={}
                gu = GUI(root,option_db,option_list['GroupArrays'],'groups')
                for option in user_variables: ### By default, all arrays will be assigned a group of ''
                    if len(user_variables[option])>0: user_variables_long[option]=[]
                ###Store the group names and assign group numbers
                group_name_db={}; group_name_list = []; group_number = 1
                for cel_file in option_list['GroupArrays']: ### start we these CEL files, since they are ordered according to their order in the expression dataset
                    group_name = gu.Results()[cel_file]
                    if group_name not in group_name_db:
                        group_name_db[group_name]=group_number; group_number+=1
                        group_name_list.append(group_name)
                
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
            category = 'SetupComps'; option_db={}; option_list={}; user_variables={}; cn = 0
            while cn < px:
                try: group1,group2 = original_comp_group_list[cn]
                except IndexError: group1='';group2=''
                cn+=1; option = 'comparison '+str(cn); array_options = group_name_list; displayed_title=option; display_object='pulldown_comps'; notes=[group1,group2]
                od = OptionData(option,displayed_title,display_object,notes,array_options,'')
                option_db[option] = od
                try: option_list[category].append(option) ###group is the name of the GUI menu group
                except KeyError: option_list[category] = [option]

            root = Tk(); root.title('AltAnalyze: Establish All Pairwise Comparisons')
            gu = GUI(root,option_db,option_list['SetupComps'],'comps'); comp_groups_db={}
            ### Sort comparisons from user for export
            for comparison in gu.Results():
                group_name = gu.Results()[comparison]
                if len(group_name)>0: ### Group_names are by default blank
                    cn_main,cn_minor = string.split(comparison[11:],'-') ### e.g. 1-1 and 1-2
                    try: comp_groups_db[cn_main].append([cn_minor,group_name])
                    except KeyError: comp_groups_db[cn_main]=[[cn_minor,group_name]]
            comp_group_list=[]
            for cn_main in comp_groups_db: cg = comp_groups_db[cn_main]; cg.sort(); comp_group_list.append([cn_main,[group_name_db[cg[0][1]],group_name_db[cg[1][1]]]]) 
            comp_group_list.sort()
            
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
                    exportComps(exp_file_location_db,comp_group_list)
                    exported = 1
                except Exception:
                    print_out = "The file:\n"+comps_file+"\nis still open. This file must be closed before proceeding"
                    IndicatorWindow(print_out,'Continue')
        ### See if there are any Affymetrix annotation files for this species
        import_dir = '/AltDatabase/affymetrix/'+species; dir_list = read_directory(import_dir)
        fn_dir = filepath(import_dir[1:])
        if len(dir_list)<1 and array_type != 'exon':
            print_out = 'No Affymetrix annnotations file found in the directory:\n'+fn_dir
            print_out += '\n\nTo download, click on the below button, find your array and download the annotation CSV file'
            print_out += '\nlisted under "Current NetAffx Annotation Files". Extract the compressed zip archive to the'
            print_out += '\nabove listed directory and hit continue to include these annotations in your results file.'
    
            button_text = 'Download Annotations'; url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
            IndicatorLinkOutWindow(print_out,button_text,url)
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
    apt_location = getAPTLocations(file_location_defaults,run_from_scratch,exportTransitResultsforAnalysis)

    ### Set the primary parent directory for ExpressionBuilder and AltAnalyze (one level above the ExpressionInput directory, if present)
    for dataset in exp_file_location_db:
        fl = exp_file_location_db[dataset_name]
        fl.setAPTLocation(apt_location)
        if run_from_scratch == 'CEL files':
            fl.setInputCDFFile(input_cdf_file); fl.setCLFFile(clf_file); fl.setBGPFile(bgp_file); fl.setXHybRemoval(remove_xhyb)
            fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type_original); fl.setOutputDir(output_dir)
        fl = exp_file_location_db[dataset]; fl.setRootDir(parent_dir)

    expr_var = species,array_type,manufacturer,constitutive_source,dabg_p,expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch, perform_alt_analysis
    alt_var = analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_cutoff,gene_expression_cutoff,permute_p_threshold, perform_permutation_analysis, export_splice_index_values
    additional_var = calculate_splicing_index_p, exportTransitResultsforAnalysis, analyze_functional_attributes, microRNA_prediction_method, filter_for_AS
    return expr_var, alt_var, additional_var, exp_file_location_db

def getAPTLocations(file_location_defaults,run_from_scratch,exportTransitResultsforAnalysis):
    import ResultsExport_module
    if 'APT' in file_location_defaults:
        for fl in file_location_defaults['APT']: 
          apt_location = fl.Location() ###Only one entry for all species
          if len(apt_location)<1: ###If no APT version is designated, prompt the user to find the directory
              if run_from_scratch == 'CEL_summarize':
                  print_out = 'To proceed with probeset summarization from CEL files,\nyou must select a valid Affymetrix Power Tools Directory.'
              elif exportTransitResultsforAnalysis == 'yes': 
                  print_out = "To proceed with running MiDAS, you must select\na valid Affymetrix Power Tools Directory."
              win_info = IndicatorChooseWindow(print_out,'Continue') ### Prompt the user to locate the APT directory
              apt_location = win_info.Folder()
              fl.ResetLocation(apt_location)
              exportDefaultFileLocations(file_location_defaults)
    return apt_location

if __name__ == '__main__':
    getUserParameters('yes'); sys.exit()
    ###Test probesetSummarize
    exp_file_location_db={}
    fl = ExpressionFileLocationData('','','','')
    exp_file_location_db['sk9_mutant'] = fl

    apt_location = 'AltDatabase/affymetrix/APT'
    input_cdf_file = 'C:/Documents and Settings/Nathan Salomonis/Desktop/zebrafish_libraryfile/CD_Zebrafish/Full/Zebrafish/LibFiles/Zebrafish.cdf'
    clf_file=''; bgp_file=''; array_type = "3'array"
    output_dir = 'C:/Documents and Settings/Nathan Salomonis/Desktop/test-AA'
    species = 'Hs'
    cel_file_dir = 'C:/Documents and Settings/Nathan Salomonis/My Documents/1-collaborations/Keerthi'
    
    fl.setAPTLocation(apt_location); fl.setInputCDFFile(input_cdf_file); fl.setCLFFile(clf_file); fl.setBGPFile(bgp_file); fl.setCELFileDir(cel_file_dir); fl.setArrayType(array_type); fl.setOutputDir(output_dir)
            
    probesetSummarize(exp_file_location_db,species,'')
    