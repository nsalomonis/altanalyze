### Original source from: https://github.com/icetime/pyinfor/blob/master/venn.py

try:
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        matplotlib.rcParams['backend'] = 'TkAgg'
        import matplotlib.pyplot as pylab
        matplotlib.rcParams['axes.linewidth'] = 0.5
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.sans-serif'] = 'Arial'
except Exception:
    None
    
import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

from matplotlib import pyplot as pylab
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable
import UI
import time, datetime
today = str(datetime.date.today()); today = string.split(today,'-'); today = today[0]+''+today[1]+''+today[2]
time_stamp = string.replace(time.ctime(),':',''); time_stamp = string.replace(time_stamp,'  ',' ')
time_stamp = string.split(time_stamp,' '); time_stamp = today+'-'+time_stamp[3]
venn_export = 'Standard_VennDiagram-'+time_stamp
venn_export_weighted = 'Weighted_VennDiagram-'+time_stamp

#--------------------------------------------------------------------
alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}

#--------------------------------------------------------------------
def venn(data, names=None, fill="number", show_names=True, show_plot=True, outputDir=False, **kwds):
    """
data: a list
names: names of groups in data
fill = ["number"|"logic"|"both"], fill with number, logic label, or both
show_names = [True|False]
show_plot = [True|False]
"""

    if data is None:
        raise Exception("No data!")
    if len(data) == 2:
        venn2(data, names, fill, show_names, show_plot, outputDir, **kwds)
    elif len(data) == 3:
        venn3(data, names, fill, show_names, show_plot, outputDir, **kwds)
    elif len(data) == 4:
        venn4(data, names, fill, show_names, show_plot, outputDir, **kwds)
    else:
        print len(data), 'files submitted, must be less than 4 and greater than 1...'
        #raise Exception("currently only 2-4 sets venn diagrams are supported")

#--------------------------------------------------------------------
def get_labels(data, fill="number"):
    """
to get a dict of labels for groups in data

input
data: data to get label for
fill = ["number"|"logic"|"both"], fill with number, logic label, or both

return
labels: a dict of labels for different sets

example:
In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
Out[12]:
{'001': '001: 0',
'010': '010: 5',
'011': '011: 0',
'100': '100: 3',
'101': '101: 2',
'110': '110: 2',
'111': '111: 3'}
"""

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)] # sets for separate groups
    s_all = set(chain(*data)) # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value
    """for i in set_collections:
        if len(set_collections[i])<100:
            print set_collections[i]"""
    labels={}
    if fill == "number":
        for k in set_collections:
            #labels[k] = len(set_collections[k])
            labels[k] = set_collections[k]
    elif fill == "logic":
        for k in set_collections: labels[k] = k
    elif fill == "both":
        for k in set_collections: labels[k] = ("%s: %d" % (k, len(set_collections[k])))
        #labels = {k: ("%s: %d" % (k, len(set_collections[k]))) for k in set_collections}
    else: # invalid value
        raise Exception("invalid value for fill")

    return labels

#--------------------------------------------------------------------
def venn2(data=None, names=None, fill="number", show_names=True, show_plot=True, outputDir=False, **kwds):
    global coordinates
    coordinates={}
    
    if (data is None) or len(data) != 2:
        raise Exception("length of data should be 2!")
    if (names is None) or (len(names) != 2):
        names = ("set 1", "set 2")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (11, 8)

    fig = pylab.figure(figsize=figsize)
    ax = fig.gca(); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_xlim(0, 11); ax.set_ylim(0, 8)

    # r: radius of the circles
    # (x1, y1), (x2, y2): center of circles
    r, x1, y1, x2, y2 = 2.0, 3.0, 4.0, 5.0, 4.0

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 2:
        colors = kwds['colors']
    else:
        colors = ['red', 'green']

    c1 = Circle((x1,y1), radius=r, alpha=0.5, color=colors[0])
    c2 = Circle((x2,y2), radius=r, alpha=0.5, color=colors[1])

    ax.add_patch(c1)
    ax.add_patch(c2)

    ## draw text
    #1
    pylab.text(round(x1-r/2), round(y1), len(labels['10']),fontsize=16, picker=True, **alignment); coordinates[round(x1-r/2), round(y1)]=labels['10']
    pylab.text(round(x2+r/2), round(y2), len(labels['01']),fontsize=16, picker=True, **alignment); coordinates[round(x2+r/2), round(y2)]=labels['01']
    # 2
    pylab.text(round((x1+x2)/2), round(y1), len(labels['11']),fontsize=16, picker=True, **alignment); coordinates[round((x1+x2)/2), round(y1)]=labels['11']
    # names of different groups
    if show_names:
        pylab.text(x1, y1-1.2*r, names[0], fontsize=16, **alignment)
        pylab.text(x2, y2-1.2*r, names[1], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    fig.canvas.mpl_connect('pick_event', onpick)

    try:
        if outputDir!=False:
            filename = outputDir+'/%s.pdf' % venn_export
            pylab.savefig(filename)
            filename = outputDir+'/%s.png' % venn_export
            pylab.savefig(filename, dpi=100) #,dpi=200
    except Exception:
        print 'Image file not saved...'
        
    if show_plot:
        pylab.show()

    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

#--------------------------------------------------------------------
def venn3(data=None, names=None, fill="number", show_names=True, show_plot=True, outputDir=False, **kwds):
    global coordinates
    coordinates={}
    if (data is None) or len(data) != 3:
        raise Exception("length of data should be 3!")
    if (names is None) or (len(names) != 3):
        names = ("set 1", "set 2", "set 3")
    
    labels = get_labels(data, fill=fill)
    
    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (11, 8.3)

    fig = pylab.figure(figsize=figsize) # set figure size
    ax = fig.gca()
    ax.set_aspect("equal") # set aspect ratio to 1
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_xlim(0, 11); ax.set_ylim(0, 8.3)

    # r: radius of the circles
    # (x1, y1), (x2, y2), (x3, y3): center of circles
    r, x1, y1, x2, y2 = 2.0, 3.0, 3.0, 5.0, 3.0
    x3, y3 = (x1+x2)/2.0, y1 + 3**0.5/2*r

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 3:
        colors = kwds['colors']
    else:
        colors = ['red', 'green', 'blue']

    c1 = Circle((x1,y1), radius=r, alpha=0.5, color=colors[0])
    c2 = Circle((x2,y2), radius=r, alpha=0.5, color=colors[1])
    c3 = Circle((x3,y3), radius=r, alpha=0.5, color=colors[2])
    for c in (c1, c2, c3):
        ax.add_patch(c)

    ## draw text
    # 1
    pylab.text(x1-r/2, round(y1-r/2), len(labels['100']),fontsize=16, picker=True, **alignment); coordinates[x1-r/2, round(y1-r/2)]=labels['100']
    pylab.text(x2+r/2, round(y2-r/2), len(labels['010']), fontsize=16, picker=True, **alignment); coordinates[x2+r/2, round(y2-r/2)]=labels['010']
    pylab.text((x1+x2)/2, round(y3+r/2), len(labels['001']), fontsize=16, picker=True, **alignment); coordinates[(x1+x2)/2, round(y3+r/2)]=labels['001']
    # 2
    pylab.text((x1+x2)/2, round(y1-r/2), len(labels['110']),fontsize=16,picker=True, **alignment); coordinates[(x1+x2)/2, round(y1-r/2)]=labels['110']
    pylab.text(x1, round(y1+2*r/3), len(labels['101']),fontsize=16,picker=True, **alignment); coordinates[x1, round(y1+2*r/3)]=labels['101']
    pylab.text(x2, round(y2+2*r/3), len(labels['011']),fontsize=16,picker=True, **alignment); coordinates[x2, round(y2+2*r/3)]=labels['011']
    # 3
    pylab.text((x1+x2)/2, round(y1+r/3),len(labels['111']), fontsize=16,picker=True, **alignment); coordinates[(x1+x2)/2, round(y1+r/3)]=labels['111']
    # names of different groups
    if show_names:
        pylab.text(x1-r, y1-r, names[0], fontsize=16, **alignment)
        pylab.text(x2+r, y2-r, names[1], fontsize=16, **alignment)
        pylab.text(x3, y3+1.2*r, names[2], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    fig.canvas.mpl_connect('pick_event', onpick)
    
    try:
        if outputDir!=False:
            filename = outputDir+'/%s.pdf' % venn_export
            pylab.savefig(filename)
            filename = outputDir+'/%s.png' % venn_export
            pylab.savefig(filename, dpi=100) #,dpi=200
    except Exception:
        print 'Image file not saved...'
        
    if show_plot:
        pylab.show()

    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

#--------------------------------------------------------------------
def venn4(data=None, names=None, fill="number", show_names=True, show_plot=True, outputDir=False, **kwds):
    global coordinates
    coordinates={}
    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (11, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize) # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110 # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5, label=names[0])) ### some disconnect with the colors and the labels (see above) so had to re-order here
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[2], alpha=0.5, label=names[2]))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[3], alpha=0.5, label=names[3]))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[1], alpha=0.5, label=names[1]))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, len(labels['1000']),fontsize=16, picker=True, **alignment);coordinates[120, 200]=labels['1000']
    pylab.text(280, 200, len(labels['0100']),fontsize=16, picker=True, **alignment);coordinates[280, 200]=labels['0100']
    pylab.text(155, 250, len(labels['0010']),fontsize=16, picker=True, **alignment);coordinates[155, 250]=labels['0010']
    pylab.text(245, 250, len(labels['0001']),fontsize=16, picker=True, **alignment);coordinates[245, 250]=labels['0001']
    # 2
    pylab.text(200, 115, len(labels['1100']),fontsize=16, picker=True, **alignment);coordinates[200, 115]=labels['1100']
    pylab.text(140, 225, len(labels['1010']),fontsize=16, picker=True, **alignment);coordinates[140, 225]=labels['1010']
    pylab.text(145, 155, len(labels['1001']),fontsize=16, picker=True, **alignment);coordinates[145, 155]=labels['1001']
    pylab.text(255, 155, len(labels['0110']),fontsize=16, picker=True, **alignment);coordinates[255, 155]=labels['0110']
    pylab.text(260, 225, len(labels['0101']),fontsize=16, picker=True, **alignment);coordinates[260, 225]=labels['0101']
    pylab.text(200, 240, len(labels['0011']),fontsize=16, picker=True, **alignment);coordinates[200, 240]=labels['0011']
    # 3
    pylab.text(235, 205, len(labels['0111']),fontsize=16, picker=True, **alignment);coordinates[235, 205]=labels['0111']
    pylab.text(165, 205, len(labels['1011']),fontsize=16, picker=True, **alignment);coordinates[165, 205]=labels['1011']
    pylab.text(225, 135, len(labels['1101']),fontsize=16, picker=True, **alignment);coordinates[225, 135]=labels['1101']
    pylab.text(175, 135, len(labels['1110']),fontsize=16, picker=True, **alignment);coordinates[175, 135]=labels['1110']
    # 4
    pylab.text(200, 175, len(labels['1111']),fontsize=16, picker=True, **alignment);coordinates[200, 175]=labels['1111']
    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

    leg = ax.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    fig.canvas.mpl_connect('pick_event', onpick)
    
    try:
        if outputDir!=False:
            filename = outputDir+'/%s.pdf' % venn_export
            pylab.savefig(filename)
            filename = outputDir+'/%s.png' % venn_export
            pylab.savefig(filename, dpi=100) #,dpi=200
    except Exception:
        print 'Image file not saved...'
        
    if show_plot:
        pylab.show()

    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass
#--------------------------------------------------------------------

def onpick(event):
        text = event.artist
        x = int(event.mouseevent.xdata)
        y = int(event.mouseevent.ydata)
        #print [text.get_text()]
  
        if (x,y) in coordinates:
            #print text.get_text()
            from visualization_scripts import TableViewer
            header = ['Associated Genes']
            tuple_list = []
            
            for gene in coordinates[(x,y)]:
                tuple_list.append([(gene)])
            TableViewer.viewTable(text.get_text(),header,tuple_list) #"""
            
def test():
    """ a test function to show basic usage of venn()"""

    # venn3()
    venn([range(10), range(5,15), range(3,8)], ["aaaa", "bbbb", "cccc"], fill="both", show_names=False)
    # venn2()
    venn([range(10), range(5,15)])
    venn([range(10), range(5,15)], ["aaaa", "bbbb"], fill="logic", show_names=False)
    # venn4()
    venn([range(10), range(5,15), range(3,8), range(4,9)], ["aaaa", "bbbb", "cccc", "dddd"], figsize=(12,12))

######### Added for AltAnalyze #########

def test2():
    """ a test function to show basic usage of venn()"""
    venn([['a','b','c','d','e'], ['a','b','c','d','e'], ['a','b','c','d','e']], ["a", "b", "c"], fill="number", show_names=True)

def compareDirectoryFiles(import_dir):
    file_list = UI.read_directory(import_dir)
    compareImportedTables(file_list,import_dir,importDir=import_dir)
    
def compareInputFiles(file_list,outputDir,display=True):
    compareImportedTables(file_list,outputDir,display=display)
    
def compareImportedTables(file_list,outputDir,importDir=False,considerNumericDirection=False,display=True):
    ### added for AltAnalyze
    print 'Creating Venn Diagram from input files...'
    import UI
    import export
    file_id_db={}
    file_list2=[]
    for file in file_list:
        x=0
        if '.txt' in file:
            if importDir !=False: ### When all files in a directory are analyzed
                fn=UI.filepath(import_dir+'/'+file)
            else:
                fn = file
                file = export.findFilename(fn) ### Only report the actual filename
            file_list2.append(file)
            for line in open(fn,'rU').xreadlines():
                if x == 0:
                    data_type = examineFields(line)
                    x+=1
                else:
                    data = UI.cleanUpLine(line)
                    t = string.split(data,'\t')
                    uid = t[0]
                    if string.lower(uid) == 'uid':
                        continue
                    valid = True
                    if data_type != 'first':
                        if data_type == 'comparison':
                            score = float(string.split(t[6],'|')[0])
                            if 'yes' not in t[5]:
                                valid = False ### not replicated independently
                        if data_type == 'reciprocal':
                            uid = t[8]+'-'+t[10]
                            score = float(t[1])
                        if data_type == 'single':
                            uid = t[6]
                            score = float(t[1])
                    else:
                        try:
                            score = float(t[1]) #t[2]
                        except Exception: score = None
                    if score != None and considerNumericDirection: ### change the UID so that it only matches if the same direction
                        if score>0:
                            uid+='+' ### encode the ID with a negative sign
                        else:
                            uid+='-' ### encode the ID with a negative sign
                    #if score>0:
                    if valid:
                        try: file_id_db[file].append(uid)
                        except Exception: file_id_db[file] = [uid]
                    
    id_lists=[]
    new_file_list=[]
    for file in file_list2: ### Use the sorted names
        if file in file_id_db:
            uids = file_id_db[file]
            id_lists.append(uids)
            new_file_list.append(file)
            #print file, len(new_file_list), len(uids)
            
    if len(file_id_db):
        if len(new_file_list)==2 or len(new_file_list)==3:
            SimpleMatplotVenn(new_file_list,id_lists,outputDir=outputDir,display=False) ### display both below
        venn(id_lists, new_file_list, fill="number", show_names=False, outputDir=outputDir, show_plot=display)
        
def examineFields(data):
    if 'independent confirmation' in data:
        data_type = 'comparison'
    elif 'norm-p2' in data:
        data_type = 'reciprocal'
    elif 'functional_prediction' in data:
        data_type = 'single'
    else:
        data_type = 'first' ### first column
    return data_type

def SimpleMatplotVenn(names,data,outputDir=False,display=True):
    """ Uses http://pypi.python.org/pypi/matplotlib-venn (code combined into one module) to export
    simple or complex, overlapp weighted venn diagrams as an alternative to the default methods in
    this module """

    import numpy as np
    pylab.figure(figsize=(11,7),facecolor='w')
    
    vd = get_labels(data, fill="number")
    set_labels=[]
    for i in names:
        set_labels.append(string.replace(i,'.txt',''))
        
    if len(set_labels)==2:
        from matplotlib_venn import venn2, venn2_circles
        set_colors = ('r', 'g')
        subsets = (vd['10'], vd['01'], vd['11'])
        v = venn2(subsets=subsets, set_labels = set_labels, set_colors=set_colors)
        c = venn2_circles(subsets=subsets, alpha=0.5, linewidth=1.5, linestyle='dashed')
        
    if len(set_labels)==3:
        from matplotlib_venn import venn3, venn3_circles
        set_colors = ('r', 'g', 'b')
        subsets = (vd['100'], vd['010'], vd['110'], vd['001'], vd['101'], vd['011'], vd['111'])
        v = venn3(subsets=subsets, set_labels = set_labels,set_colors=set_colors)
        c = venn3_circles(subsets=subsets, alpha=0.5, linewidth=1.5, linestyle='dashed')
        
    pylab.title("Overlap Weighted Venn Diagram",fontsize=24)

    try:
        if outputDir!=False:
            filename = outputDir+'/%s.pdf' % venn_export_weighted
            pylab.savefig(filename)
            filename = outputDir+'/%s.png' % venn_export_weighted
            pylab.savefig(filename, dpi=100) #,dpi=200
    except Exception:
        print 'Image file not saved...'
        
    if display:
        pylab.show()

    try:
        import gc
        fig.clf()
        pylab.close()
        gc.collect()
    except Exception:
        pass

 ######### End Added for AltAnalyze #########

if __name__ == '__main__':
    names = ['a','b','c']
    venn_data = {}
    venn_data['100']=1
    venn_data['010']=2
    venn_data['110']=3
    venn_data['001']=4
    venn_data['101']=5
    venn_data['011']=6
    venn_data['111']=7
    
    #SimpleMatplotVenn(names,venn_data); sys.exit()
    
    import_dir = '/Users/nsalomonis/Desktop/code/python/combine-lists/input'
    import_dir = '/Users/nsalomonis/Downloads/GSE31327_RAW/GO-Elite/input/compare'
    import_dir = '/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/AltResults/AlternativeOutput/ASPIRE-comp'
    import_dir = '/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/AltResults/AlternativeOutput/SI-comp'
    import_dir = '/Users/nsalomonis/Desktop/dataAnalysis/r4_Bruneau_TopHat/AltResults/AlternativeOutput/stringent-comp'
    import_dir = '/Users/nsalomonis/Desktop/dataAnalysis/collaborations/Faith/MarkerClustering/NC-type-markers'
    compareDirectoryFiles(import_dir)