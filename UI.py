import math
import statistics
import sys, string
import os.path
import unique
import time
from sys import argv
try:
    import Tkinter
    from Tkinter import *
    import PmwFreeze
    from Tkconstants import LEFT
    import tkMessageBox
    import tkFileDialog
except ImportError: print "\nPmw or Tkinter not found... proceeding with manual input"
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
    dir_list = os.listdir(dir + sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

################# GUI #################

class GUI:
    def __init__(self, parent, option_db, option_list, defaults): 
        self._parent = parent; self._option_list = option_list
        self._user_variables = user_variables
        i = 0
        for option in option_list:
            od = option_db[option]; self.title = od.Display()        
            self.display_options = od.ArrayOptions()
            #if len(defaults)>0: print i, od.DisplayObject(),self.display_options, option_list, defaults
            if option == 'array_type':
                valid_display_options=[]
                for array_name in self.display_options:
                    compatible_species = array_codes[array_name].SpeciesCodes()
                    if species in compatible_species: valid_display_options.append(array_name)
                self.display_options = valid_display_options
            if 'radio' in od.DisplayObject() and self.display_options != ['NA']:
                ### Create and pack a RadioSelect widget, with radiobuttons.
                self._option = option
                def radiocallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                radiobuttons = PmwFreeze.RadioSelect(self._parent,                       
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

            if 'button' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                def buttoncallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                horiz = PmwFreeze.RadioSelect(self._parent,
                        labelpos = 'w', command = buttoncallback,
                        label_text = self.title, frame_borderwidth = 2,
                        frame_relief = 'ridge'
                ); horiz.pack(fill = 'x', padx = 10, pady = 10)

                ### Add some buttons to the horizontal RadioSelect
                for text in self.display_options:
                    if text != ['NA']: horiz.add(text)
                horiz.invoke(self.default_option)
                
            if 'enter' in od.DisplayObject() and self.display_options != ['NA']:
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
                    self.default_option = 'CHANGE TO A NUMERIC VALUE'; use_method = 'i'
                if use_method == 'p':
                    entry_field = PmwFreeze.EntryField(self._parent,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate_p, 
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                    ); entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                if use_method == 'i':
                    entry_field = PmwFreeze.EntryField(self._parent,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate,
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                                       
                    ); entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                
            if 'multiple-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                if len(defaults) <1: self.default_option = self.display_options[0]
                else: self.default_option = defaults[i]
                ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                self.checkbuttons = PmwFreeze.RadioSelect(self._parent,
                        buttontype = 'checkbutton', orient = 'vertical',
                        labelpos = 'w', command = self.checkbuttoncallback,
                        label_text = self.title, hull_borderwidth = 2, hull_relief = 'ridge',
                ); self.checkbuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                ### Add some buttons to the checkbutton RadioSelect.
                for text in self.display_options:
                     if text != ['NA']: self.checkbuttons.add(text)
                self.checkbuttons.invoke(self.default_option)
                self.checkbuttons.invoke(self.default_option2)
                
            if 'single-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
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
                        self.checkbuttons = PmwFreeze.RadioSelect(self._parent,
                                buttontype = 'checkbutton', command = checkbuttoncallback,
                                hull_borderwidth = 2, hull_relief = 'ridge',
                        ); self.checkbuttons.pack(side = 'left', expand = 1, padx = 10, pady = 10)

                        ### Add some buttons to the checkbutton RadioSelect.
                        self.checkbuttons.add(self.title)
                        if self.default_option == 'yes': self.checkbuttons.invoke(self.title)
                        else: self._user_variables[option] = 'no'
            i+=1 ####Keep track of index

        #def quitcommand(): parent.destroy; sys.exit()
        
        #self.button = Button(text="   Quit  ", command=quitcommand)
        #self.button.pack(side = 'bottom', padx = 10, pady = 10)
        
        continue_to_next_win = Button(text = 'Continue', command = self._parent.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)

        quit_win = Button(self._parent, text="Quit", command=self.quit) 
        quit_win.pack(side = 'right', padx =10, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()

    def info(self):
        tkMessageBox.showinfo("title","message",parent=self._parent)

    def deleteWindow(self):
        tkMessageBox.showwarning("Quit","Use 'Quit' button to end program!",parent=self._parent)

    def quit(self):
        #print "quit starts"
        #print "cleaning up things..."
        self._parent.quit()
        self._parent.destroy()
        sys.exit()
        #print "quit ends"
        
    def Report(self,tag,option):
        output = tag
        return output
    def __repr__(self,tag,option): return self.Report(tag,option)
    
    def Results(self): return self._user_variables

    def custom_validate(self, text, option):
        #print [option],'text:', text
        self._user_variables[option] = text
        try:
            text = float(text);return 1
        except ValueError: return -1

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
    print '1) Run analyses from scratch using raw input files ("ExpressionInput")'
    print '2) Directly run alternative exon analyses using pre-processed files ("AltExpression")'
    print '3) Update existing AltAnalyze databases'
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": run_from_scratch = 'raw input'
    elif inp == "2": run_from_scratch = 'pre-processed'
    elif inp == "3": run_from_scratch = 'update DBs'
    
    print "Proceed Using Default Parameters?"
    print '1) Use defaults (see "Config/default-**.txt")'
    print '2) Change defaults (see "Config/options.txt" and associated default files)'
    inp = sys.stdin.readline(); inp = inp.strip()
    if inp == "1": proceed = 'yes'
    elif inp == "2": proceed = 'no'; sys.exit()
    
    return species_code, array_code, manufacturer, constitutive_source, run_from_scratch

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

class ArrayData:
    def __init__(self, abrev, array, manufacturer, constitutive_source, species):
        self._abrev = abrev; self._array = array; self._manufacturer = manufacturer; self._species = species
        self._constitutive_source = constitutive_source
    def ArrayCode(self): return self._abrev
    def ArrayName(self): return self._array
    def Manufacturer(self): return self._manufacturer
    def ConstitutiveSource(self): return self._constitutive_source
    def SpeciesCodes(self): return self._species
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

class FileLocationData:
    def __init__(self, status, location, species):
        self._status = status; self._location = location; self._species = species
    def Status(self): return self._status
    def Location(self): return self._location
    def Species(self): return self._species
    def __repr__(self): return self.Report()
    
def importDefaultFileLocations():
    filename = 'Config/default-files.csv'; x=0
    fn=filepath(filename); file_location_defaults={}
    for line in open(fn,'rU').readlines():
        line = string.replace(line,',','\t') ### Make tab-delimited (had to make CSV since Excel would impoperly parse otherwise)
        data = cleanUpLine(line)
        app,status,location,species = string.split(data,'\t')
        fl = FileLocationData(status, location, species)
        if x==0: x=1
        else:
            if species == 'all': file_location_defaults[app] = fl
            else:
                try: file_location_defaults[app].append(fl)
                except KeyError: file_location_defaults[app] = [fl]
    return file_location_defaults

class Defaults:
    def __init__(self, abrev, array, species):
        self._abrev = abrev; self._array = array; self._species = species
    def ArrayCode(self): return self._abrev
    def ArrayName(self): return self._array
    def Species(self): return self._species
    def __repr__(self): return self.Report()

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
            array_abrev, analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_variable, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, MiDAS_analysis, export_splice_index_values = string.split(data,'\t')
            if array_type == array_abrev:
                ### NOTE: p_threshold is used for MiDAS and t-test comparison p-values thresholds
                if MiDAS_analysis == 'yes': exportTransitResultsforAnalysis = 'yes'
                else: exportTransitResultsforAnalysis = 'no'
                if  analysis_method != 'splicing-index': exportTransitResultsforAnalysis='NA'; export_splice_index_values='NA'
                return [analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_variable, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, exportTransitResultsforAnalysis, export_splice_index_values]
            
        if '-funct' in filename:
            array_abrev, analyze_functional_attributes,microRNA_prediction_method = string.split(data,'\t')
            if array_type == array_abrev:
                return [analyze_functional_attributes,microRNA_prediction_method]

class OptionData:
    def __init__(self,option,displayed_title,display_object,description,array_options):
        self._option = option; self._displayed_title = displayed_title; self._description = description
        self._array_options = array_options; self._display_object = display_object
    def Option(self): return self._option
    def Display(self): return self._displayed_title
    def DisplayObject(self): return self._display_object
    def Description(self): return self._description
    def ArrayOptions(self): return self._array_options
    def setArrayOptions(self,array_options): self._array_options = array_options
    def __repr__(self): return self.Report()

def importUserOptions(array_type):
    filename = 'Config/options.txt'; option_db={}; option_list=[]
    fn=filepath(filename); x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        data = string.replace(data,'\k','\n') ###Used \k in the file instead of \n, since these are removed above
        t = string.split(data,'\t')
        option,displayed_title,display_object,description = t[:4]
        if x == 0:
            i = t.index(array_type) ### Index position of the name of the array_type selected by user (or arbitrary to begin with)
            x = 1
        else:
            array_options = t[i]
            array_options = string.split(array_options,'|')
            od = OptionData(option,displayed_title,display_object,description,array_options)
            option_db[option] = od
            option_list.append(option)
    return option_list,option_db

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
        filename = 'logo.gif'
        fn=filepath(filename)
        img = PhotoImage(file=fn)
        can = Canvas(parent)
        can.pack(fill=BOTH)
        can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)

        """      
        ### Create and pack a horizontal RadioSelect widget.
        def buttoncallback(tag,callback=self.callback):
            callback(tag)
        horiz = PmwFreeze.RadioSelect(parent,
                labelpos = 'w', command = buttoncallback,
                label_text = 'GO-Elite version 1.0 Main', frame_borderwidth = 2,
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
        about = 'AltAnalyze 1.01 beta.\n'
        about+= 'AltAnalyze is an open-source, freely available application covered under the\n'
        about+= 'Apache open-source license. Additional information can be found at:\n'
        about+= 'http://www.genmapp.org/AltAnalyze\n'
        about+= '\nDeveloped by:\n\tNathan Salomonis\n\tBruce Conklin\nGladstone Institutes 2008'
        tkMessageBox.showinfo("About AltAnalyze",about,parent=self._parent)

    def deleteWindow(self):
        tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)

    def callback(self, tag):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables['continue'] = tag        
                
def getUserParameters(skip_intro):
    if skip_intro == 'yes':
        try: MainMenu()
        except TclError: null=[]
    global species; species=''; global user_variables; user_variables={}; global analysis_method; global array_type 
    ### Get default options for ExpressionBuilder and AltAnalyze

    na = 'NA'; log = 'log'; no = 'no'
    run_from_scratch=na; expression_threshold=na; perform_alt_analysis=na; expression_data_format=log
    include_raw_data=na; avg_all_for_ss=na; dabg_p=na;
    analysis_method=na; p_threshold=na; filter_probeset_types=na; alt_exon_fold_cutoff=na
    permute_p_threshold=na; perform_permutation_analysis=na; export_splice_index_values=no
    exportTransitResultsforAnalysis=no; analyze_functional_attributes=no; microRNA_prediction_method=na
    gene_expression_cutoff='any'
            
    option_list,option_db = importUserOptions('exon')  ##Initially used to just get the info for species and array_type
    importSpeciesInfo()
    file_location_defaults = importDefaultFileLocations()
    array_list = importArrayInfo()
    try:   
        ###Update this informatin in option_db which will be over-written after the user selects a species and array_type
        option_db['species'].setArrayOptions(species_list)

        root = Tk()
        root.title('AltAnalyze: Main Dataset Parameters')
        gu = GUI(root,option_db,option_list[:1],'')
        species_full = gu.Results()['species']
        species = species_codes[species_full].SpeciesCode()

        x = 1; array_list2=[]
        for array_name in array_list:
            if species in array_codes[array_name].SpeciesCodes(): array_list2.append(array_name)

        array_list = array_list2 ### Filtered based on compatible species arrays
        option_db['array_type'].setArrayOptions(array_list)
        
        root = Tk(); i = 3
        root.title('AltAnalyze: Main Dataset Parameters')
        gu = GUI(root,option_db,option_list[1:i],'')
        array_full = gu.Results()['array_type']
        run_from_scratch = gu.Results()['run_from_scratch']    
        array_type = array_codes[array_full].ArrayCode()

        manufacturer = array_codes[array_full].Manufacturer()
        constitutive_source = array_codes[array_full].ConstitutiveSource()
    
        if run_from_scratch != 'update DBs':
            expr_defaults, alt_exon_defaults, functional_analysis_defaults = importDefaults(array_type,species)
            
            if run_from_scratch != 'pre-processed':

                option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type
                root = Tk(); root.title('AltAnalyze: Expression Analysis Parameters')
                gu = GUI(root,option_db,option_list[i:i+len(expr_defaults)],expr_defaults)
                if array_type != "3'array":          
                    dabg_p = gu.Results()['dabg_p']
                    run_from_scratch = gu.Results()['run_from_scratch']
                    expression_threshold = gu.Results()['expression_threshold']
                    perform_alt_analysis = gu.Results()['perform_alt_analysis']
                    avg_all_for_ss = gu.Results()['avg_all_for_ss']
                expression_data_format = gu.Results()['expression_data_format']
                include_raw_data = gu.Results()['include_raw_data']
                
            if (perform_alt_analysis == 'both') or (run_from_scratch == 'pre-processed'):
                perform_alt_analysis = 'alt'
                
                i = i+len(expr_defaults)
                #print option_list[i:i+len(alt_exon_defaults)+len(functional_analysis_defaults)], alt_exon_defaults+functional_analysis_defaults;kill
                option_list,option_db = importUserOptions(array_type)  ##Initially used to just get the info for species and array_type
                root = Tk(); root.title('AltAnalyze: Alternative Exon Analysis Parameters')
                gu = GUI(root,option_db,option_list[i:i+len(alt_exon_defaults)+len(functional_analysis_defaults)],alt_exon_defaults+functional_analysis_defaults); user_variables = []

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
                analyze_functional_attributes = gu.Results()['analyze_functional_attributes']
                microRNA_prediction_method = gu.Results()['microRNA_prediction_method']
                try: p_threshold = float(permute_p_threshold)
                except ValueError: permute_p_threshold = permute_p_threshold
    except NameError:  ###Occurs when PMW or Tkinter does not propperly load
        species, array_type, manufacturer, constitutive_source, run_from_scratch = getPrimaryUserParameters()
        expr_defaults, alt_exon_defaults, functional_analysis_defaults = importDefaults(array_type,species)

        dabg_p, expression_threshold, perform_alt_analysis, expression_data_format, avg_all_for_ss, include_raw_data = expr_defaults
        analysis_method, p_threshold, filter_probeset_types, alt_exon_fold_cutoff, gene_expression_cutoff, perform_permutation_analysis, permute_p_threshold, exportTransitResultsforAnalysis, export_splice_index_values = alt_exon_defaults
        analyze_functional_attributes,microRNA_prediction_method = functional_analysis_defaults
        
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

    #print [exportTransitResultsforAnalysis,analyze_functional_attributes,microRNA_prediction_method,avg_all_for_ss,include_raw_data];kill
    expr_var = species,array_type,manufacturer,constitutive_source,dabg_p,expression_threshold,avg_all_for_ss,expression_data_format,include_raw_data, run_from_scratch, perform_alt_analysis
    alt_var = analysis_method,p_threshold,filter_probeset_types,alt_exon_fold_cutoff,gene_expression_cutoff,permute_p_threshold, perform_permutation_analysis, export_splice_index_values
    additional_var = exportTransitResultsforAnalysis, analyze_functional_attributes, microRNA_prediction_method
    return expr_var, alt_var, additional_var, file_location_defaults

if __name__ == '__main__':
    getUserParameters()