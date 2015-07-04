from numpy import arange, sin, pi
import os.path, sys, shutil
import wx.lib.scrolledpanel
import wx.grid as gridlib
import os
import string, re
import wx
from wx import *

import RemoteViewer ### Import itself as a reference to it's location
dirfile = RemoteViewer
application_path = os.path.dirname(__file__)
import matplotlib
#matplotlib.use('WXAgg')
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

"""
os.chdir(os.path.dirname(application_path))
back_dir = ""
cwd_root = os.getcwd()
cwd_root = cwd_root.split("/")
cwd_root = cwd_root[1:-1]
for i in cwd_root:
    back_dir = back_dir + "/" + i

#print back_dir

os.chdir(back_dir)
print os.getcwd()
"""

#os.chdir(os.path.dirname(__file__))

class PageOne(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

class PageTwo(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour("white")
        myGrid = ""        
        
class PageThree(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour("white")

class PageFour(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

class PageFive(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)


class Main(wx.Frame):
    def __init__(self,parent,id):
        wx.Frame.__init__(self, parent, id,'AltAnalyze Results Viewer', size=(1500,810))
        
        self.SetBackgroundColour((230, 230, 230))
        
        #PANELS & WIDGETS              
        self.panel = wx.Panel(self, id=2, pos=(200,0), size=(900,45), style=wx.RAISED_BORDER)       
        self.panel.SetBackgroundColour((170, 220, 170))
        #self.panel2 = wx.lib.scrolledpanel.ScrolledPanel(self, id=3, size=(1400,605), pos=(200,50), style=wx.RAISED_BORDER)
        self.panel2 = wx.Panel(self, id=3, pos=(200,50), size=(1400,605), style=wx.RAISED_BORDER)
        #self.panel2.SetupScrolling()
        self.panel2.SetBackgroundColour((218, 218, 218))
        self.panel3 = wx.Panel(self, id=4, pos=(0,50), size=(200,625), style=wx.RAISED_BORDER)
        self.panel3.SetBackgroundColour("white")
        self.panel4 = wx.Panel(self, id=5, pos=(200,650), size=(1400,150), style=wx.RAISED_BORDER)
        self.panel4.SetBackgroundColour("black")

        self.panel_left = wx.Panel(self, id=12, pos=(0,0), size=(200,45), style=wx.RAISED_BORDER)
        self.panel_left.SetBackgroundColour((218, 218, 218))
        self.panel_right = wx.Panel(self, id=11, pos=(1100,0), size=(200,45), style=wx.RAISED_BORDER)
        self.panel_right.SetBackgroundColour((218, 218, 218))
        self.panel_right2 = wx.Panel(self, id=13, pos=(1300,0), size=(300,45), style=wx.RAISED_BORDER)
        self.panel_right2.SetBackgroundColour((218, 218, 218))

        self.sortbox = wx.TextCtrl(self.panel_right2, id=7, pos=(55,10), size=(40,25))
        wx.Button(self.panel_right2, id=8, label="Sort", pos=(5, 12), size=(40, 10))
        self.Bind(wx.EVT_BUTTON, self.SortTablefromButton, id=8)

        self.AscendingRadio = wx.RadioButton(self.panel_right2, id=17, label="Sort", pos=(100, 3), size=(12, 12))
        self.DescendingRadio = wx.RadioButton(self.panel_right2, id=18, label="Sort", pos=(100, 23), size=(12, 12))

        font = wx.Font(10, wx.DECORATIVE, wx.BOLD, wx.NORMAL)
        self.AscendingOpt = wx.StaticText(self.panel_right2, label="Ascending", pos=(115, 1))
        self.AscendingOpt.SetFont(font)
        self.DescendingOpt = wx.StaticText(self.panel_right2, label="Descending", pos=(115, 21))
        self.DescendingOpt.SetFont(font)

        self.filterbox = wx.TextCtrl(self.panel_right, id=9, pos=(60,10), size=(125,25))
        wx.Button(self.panel_right, id=10, label="Filter", pos=(0, 12), size=(50, 10))
        self.Bind(wx.EVT_BUTTON, self.FilterTablefromButton, id=10)        
                        
        self.control = wx.TextCtrl(self.panel4, id=6, pos=(1,1), size=(1400,150), style=wx.TE_MULTILINE)
        self.control.write("Welcome to AltAnalyze Results Viewer!" + "\n")
        self.Show(True)                

        self.input_stream = ""

        self.browser = wx.TreeCtrl(self.panel3, id=2000, pos=(0,0), size=(200,325)) 
        self.browser2 = wx.TreeCtrl(self.panel3, id=2001, pos=(0,325), size=(200,325))
        self.tree = self.browser                                           

        #DICTIONARY
        self.sortdict = {"A" : 0, "B" : 1, "C" : 2, "D" : 3, "E" : 4, "F" : 5, "G" : 6, "H" : 7, "I" : 8, "J" : 9, "K" : 10, "L" : 11, "M" : 12, "N" : 13, "O" : 14, "P" : 15, "Q" : 16, "R" : 17, "S" : 18, "T" : 19, "U" : 20, "V" : 21, "W" : 22, "X" : 23, "Y" : 24, "Z" : 25, "AA" : 26, "AB" : 27, "AC" : 28, "AD" : 29, "AE" : 30, "AF" : 31, "AG" : 32, "AH" : 33, "AI" : 34, "AJ" : 35, "AK" : 36, "AL" : 37, "AM" : 38, "AN" : 39, "AO" : 40, "AP" : 41, "AQ" : 42, "AR" : 43, "AS" : 44, "AT" : 45, "AU" : 46, "AV" : 47, "AW" : 48, "AX" : 49, "AY" : 50, "AZ" : 51}

        #SIZER
        ver = wx.BoxSizer(wx.VERTICAL)
        verpan2 = wx.BoxSizer(wx.VERTICAL)
        hpan1 = wx.BoxSizer(wx.HORIZONTAL)
        hpan2 = wx.BoxSizer(wx.HORIZONTAL)
        hpan3 = wx.BoxSizer(wx.HORIZONTAL)
        verpan2.Add(self.panel2, 8, wx.ALL|EXPAND, 2)  
        hpan1.Add(self.panel_left, 5, wx.ALL|EXPAND, 2)
        hpan1.Add(self.panel, 24, wx.ALL|EXPAND, 2)
        hpan1.Add(self.panel_right, 3, wx.ALL|EXPAND, 2)
        hpan1.Add(self.panel_right2, 3, wx.ALL|EXPAND, 2)
        hpan2.Add(self.panel3, 1, wx.ALL|EXPAND, 2)
        hpan2.Add(verpan2, 7, wx.ALL|EXPAND, 2)
        hpan3.Add(self.panel4, 1, wx.ALL|EXPAND, 2)
        ver.Add(hpan1, 1, wx.EXPAND)
        ver.Add(hpan2, 18, wx.EXPAND) 
        ver.Add(hpan3, 4, wx.EXPAND)
        self.browser.SetSize(self.panel3.GetSize())
        self.SetSizer(ver)  

        self.nb = wx.Notebook(self.panel2, style = wx.NB_BOTTOM)
        self.page1 = PageOne(self.nb)
        self.page2 = PageTwo(self.nb)
        self.page3 = PageThree(self.nb)
        self.nb.AddPage(self.page2, "Table")        
        self.nb.AddPage(self.page1, "PNG")
        self.nb.AddPage(self.page3, "Interactive")
        self.page3.SetBackgroundColour((218, 218, 218))
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer)
        self.page1.SetBackgroundColour("white")
        self.myGrid = gridlib.Grid(self.page2)
        self.myGrid.CreateGrid(100, 100)
        gridsizer = wx.BoxSizer(wx.VERTICAL)
        gridsizer.Add(self.myGrid)
        self.page2.SetSizer(gridsizer)
        self.page2.Layout()

	#INTERACTIVE PANEL LAYOUT

        wx.Button(self.page3, id=16, label="Demo", pos=(550, 15), size=(60, 10))
        self.Bind(wx.EVT_BUTTON, self.InteractiveDemo, id=16) 

        wx.Button(self.page3, id=43, label="Run", pos=(275, 150), size=(120, 10))
        self.Bind(wx.EVT_BUTTON, self.InteractiveDemo, id=43) 

        self.ln = wx.StaticLine(self.page3, pos=(5,100))
        self.ln.SetSize((415,10))
        
        IntTitleFont = wx.Font(15, wx.SWISS, wx.NORMAL, wx.BOLD)
        self.InteractiveTitle = wx.StaticText(self.page3, label="Main Dataset Parameters", pos=(10, 15))
        self.InteractiveTitle.SetFont(IntTitleFont)        

        self.IntFileTxt = wx.TextCtrl(self.page3, id=43, pos=(105,45), size=(375,20))       
                        
        self.InteractiveFileLabel = wx.StaticText(self.page3, label="Selected File:", pos=(10, 45))                        
                                                                                
        self.Yes1Label = wx.StaticText(self.page3, label="Yes", pos=(305, 80))
        self.No1Label = wx.StaticText(self.page3, label="No", pos=(375, 80))

        self.D_3DLabel = wx.StaticText(self.page3, label="3D", pos=(305, 120))
        self.D_2DLabel = wx.StaticText(self.page3, label="2D", pos=(375, 120))
                        
        self.Yes1Radio = wx.RadioButton(self.page3, id=40, pos=(285, 83), size=(12, 12), style=wx.RB_GROUP)
        self.No1Radio = wx.RadioButton(self.page3, id=41, pos=(355, 83), size=(12, 12))

        self.D_3DRadio = wx.RadioButton(self.page3, id=46, pos=(285, 123), size=(12, 12), style=wx.RB_GROUP)
        self.D_2DRadio = wx.RadioButton(self.page3, id=47, pos=(355, 123), size=(12, 12))

        Opt1Desc = wx.StaticText(self.page3, label="Display sample labels next to each object", pos=(10, 80))
        Opt2Desc = wx.StaticText(self.page3, label="Dimensions to display", pos=(10, 120))

        self.Bind(wx.EVT_BUTTON, self.InteractiveFileChoose, id=50)

        TxtBox = wx.BoxSizer(wx.VERTICAL)   
        TxtBox.Add(self.control, 1, wx.EXPAND)
        self.panel4.SetSizer(TxtBox)
        self.panel4.Layout()

        #SELECTION LIST
        self.TopSelectList = []
        self.SearchArray = []
        self.SearchArrayFiltered = []
        self.TopID = ""
        
        self.ColoredCellList = []

        #LOGO
        self.png = wx.Image("Config/logo.gif", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
        imgsizer_v = wx.BoxSizer(wx.VERTICAL)
        imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
        self.page1.SetSizer(imgsizer_v)
        self.page1.Layout()

        self.nb.SetSelection(1)

        browspan = wx.BoxSizer(wx.VERTICAL)
        browspan.Add(self.browser, 1, wx.EXPAND)
        browspan.Add(self.browser2, 1, wx.EXPAND)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
        self.panel3.SetSizer(browspan)
        
        self.PanelTitle = wx.StaticText(self.panel, label="", pos=(210, 15))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        #TOOLBAR BUTTONS        
        toolbar1 = wx.ToolBar(self.panel_left, -1, pos=(45,15), style= wx.TB_HORIZONTAL)
        toolbar1.AddLabelTool(1000, '', wx.Bitmap("Config/quit.png"))
        toolbar1.AddLabelTool(1001, '', wx.Bitmap("Config/folder_open.png"))
        toolbar1.AddLabelTool(1002, '', wx.Bitmap("Config/save.png"))
        toolbar1.AddLabelTool(1003, '', wx.Bitmap("Config/eraser.png"))
        toolbar1.AddLabelTool(1004, '', wx.Bitmap("Config/search.png"))
        toolbar1.Realize()               
        self.Bind(wx.EVT_TOOL, self.OnQuit, id=1000)        
        self.Bind(wx.EVT_TOOL, self.OnOpen, id=1001)
        self.Bind(wx.EVT_TOOL, self.OnSave, id=1002)
        self.Bind(wx.EVT_TOOL, self.ClearVisualPanel, id=1003)
        self.Bind(wx.EVT_TOOL, self.TreeSearch, id=1004)

        #STATUS BAR CREATE             
        status = self.CreateStatusBar()
        menubar = wx.MenuBar()
                
        file = wx.Menu()
        edit = wx.Menu()
        view = wx.Menu()
        search = wx.Menu()
        filter_table = wx.Menu()
        about = wx.Menu()
        
        open_menu = wx.Menu()
        open_menu.Append(120, 'Project')
        open_menu.Append(121, 'File')
        
        file.AppendMenu(101, '&Open\tCtrl+O', open_menu)
        
        file.Append(102, '&Save\tCtrl+S', 'Save the document')
        file.AppendSeparator()
        file.Append(103, 'Options', '')
        file.AppendSeparator()
        quit = wx.MenuItem(file, 105, '&Quit\tCtrl+Q', 'Quit the Application')
        file.AppendItem(quit)

        edit.Append(109, 'Undo', '')
        edit.Append(110, 'Redo', '')  
        edit.AppendSeparator()       
        edit.Append(106, '&Cut\tCtrl+X', '')
        edit.Append(107, '&Copy\tCtrl+C', '') 
        edit.Append(108, '&Paste\tCtrl+V', '')   
        edit.AppendSeparator()
        edit.Append(111, '&Select All\tCtrl+A', '')     

        view.Append(112, '&Clear Panel\tCtrl+.', '')

        search.Append(113, 'Tree', '')
        search.Append(114, 'Table', '')
        search.Append(115, 'Project', '')
        
        filter_table.Append(116, 'Filter', '')
        filter_table.Append(117, 'Sort', '')

        menubar.Append(file, "File")
        menubar.Append(edit, "Edit")
        menubar.Append(view, "View")
        menubar.Append(search, "Search")
        menubar.Append(filter_table, "Table")
        menubar.Append(about, "About")      
        self.SetMenuBar(menubar)
        
        #STATUS BAR BINDINGS        
        self.Bind(wx.EVT_MENU, self.OnOpen, id=120)
        self.Bind(wx.EVT_MENU, self.OnOpenSingleFile, id=121)
        self.Bind(wx.EVT_MENU, self.OnQuit, id=105)
        self.Bind(wx.EVT_MENU, self.ClearVisualPanel, id=112)
        self.Bind(wx.EVT_MENU, self.TreeSearch, id=113)
        self.Bind(wx.EVT_MENU, self.GridSearch, id=114)
        self.Bind(wx.EVT_MENU, self.FilterTable, id=116)
        self.Bind(wx.EVT_MENU, self.SortTable, id=117)
                     
    def OnQuit(self, event):
        popup = wx.MessageDialog(None, "Are you sure you want to quit?", "Warning", wx.YES_NO)
        popup_answer = popup.ShowModal()
        print popup_answer
        if(popup_answer == 5103):
            self.Close()
        else:
            return

    def InteractiveFileChoose(self, event):
        openFileDialog = wx.FileDialog(self, "", "", "", "", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     
            
        single_input_stream = openFileDialog.GetPath()
        self.ChooseFileLog.Clear()
        self.ChooseFileLog.write(str(single_input_stream))        

    def OnOpen(self, event):
        openFileDialog = wx.DirDialog(None, "Choose project", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)  
        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return 
        self.input_stream = openFileDialog.GetPath()
        if (len(self.input_stream) > 0):
            
            self.SearchArray = []
            self.SearchArrayFiltered = []

            self.control.write("Working..." + "\n")

            #FLAG COLLECT
            root = 'Data'
            for (dirpath, dirnames, filenames) in os.walk(root):
                for dirname in dirnames:
                    fullpath = os.path.join(dirpath, dirname)
                for filename in sorted(filenames):
                    if filename == "location.txt":
                        file_fullpath = os.path.join(dirpath, filename)
                        #self.control.write(str(file_fullpath) + "\n")
                        file_location = open(file_fullpath, "r")
                        fl_array = []
                        for line in file_location:
                            line = line.rstrip()
                            line = line.split("\r")
                            if len(line) > 1:
                                fl_array.append(line[0])
                                fl_array.append(line[1])
                            else:
                                fl_array.append(line[0])
                        #self.control.write(str(fl_array) + "\n")
                        file_location.close()
                        if(len(fl_array) == 3):
                            fl_array.append(dirpath)
                            self.SearchArray.append(fl_array)
                                   
                                                            
            self.control.write("Opening project at: " + self.input_stream + "\n")
            print self.input_stream  
            self.browser2.DeleteAllItems()

            #SEARCH USING FLAGS
            count = 0
            for FLAG in self.SearchArray:
                if(FLAG[0][-1] != "/"):
                    SearchingFlag = FLAG[0] + "/"
                SearchingFlag = FLAG[0]
                SearchingFlagPath = self.input_stream + "/" + SearchingFlag
                try:
                    SFP_Contents = os.listdir(SearchingFlagPath)
                    for filename in SFP_Contents:
                        Search_base = FLAG[1]
                        Search_base = Search_base.split(":")
                        Search_base = Search_base[1]
                        Split_Extension = str(FLAG[2])
                        Split_Extension = Split_Extension.split(":")
                        S_E = str(Split_Extension[1]).split(",")
                        GOOD_FLAG = 0
                        if(Search_base != "*"):
                            for i in S_E:
                                if(filename[-4:] == i):
                                    GOOD_FLAG = 1
                        
                        if(Search_base != "*"):
                            candidate =  re.findall(Search_base, filename) 
                        
                        if(Search_base == "*"):
                            candidate = "True"
                            GOOD_FLAG = 1
                        if (len(Search_base) == 0 or GOOD_FLAG == 0):
                            continue                       
                        if len(candidate) > 0:
                            self.SearchArrayFiltered.append(FLAG)
                except:
                    continue
                count = count + 1

            #AVAILABLE DATA SET
            
            #for i in self.SearchArrayFiltered:
            #    self.control.write(str(i) + "\n")
            
            try:
                shutil.rmtree("Config/AvailableData")
            except:
                pass
            
            for i in self.SearchArrayFiltered:
                    #self.control.write(str(i[3]) + "  " + str(i[0]) + "  " + str(i[1]) + str(i[4]) + "\n")
                    AvailablePath = "Config/Available" + i[3]
                    Path_List = AvailablePath.split("/")
                    Created_Directory = ""
                    for directorynum in range(len(Path_List)):
                        if directorynum == 0:
                            Created_Directory = Created_Directory + Path_List[directorynum]
                            try:
                                os.mkdir(Created_Directory)
                            except:
                                continue
                        else:
                            Created_Directory = Created_Directory + "/" + Path_List[directorynum]
                            try:
                                os.mkdir(Created_Directory)
                            except:
                                continue
            
            #TOP BROWSER SET
            root = 'Config/AvailableData'
            color_root = [253, 253, 253]
            self.tree.DeleteAllItems()
            self.ids = {root : self.tree.AddRoot(root)}
            for (dirpath, dirnames, filenames) in os.walk(root):
                for dirname in dirnames:
                    fullpath = os.path.join(dirpath, dirname)
                    self.ids[fullpath] = self.tree.AppendItem(self.ids[dirpath], dirname)
                    DisplayColor = [255, 255, 255]
                    DisplayColor[0] = color_root[0] - len(dirpath)
                    DisplayColor[1] = color_root[1] - len(dirpath)
                    DisplayColor[2] = color_root[2] - len(dirpath)
                    self.tree.SetItemBackgroundColour(self.ids[fullpath], DisplayColor)
                    for i in self.SearchArrayFiltered:
                        SearchRoot = "Config/Available" + i[3]
                        if(SearchRoot == fullpath):
                            SearchSplit = i[1].split(":")
                            SearchSplit = SearchSplit[1]
                            if(SearchSplit == "*"):
                                self.control.write(SearchSplit + "\n")
                            SearchSplit = SearchSplit + ";" + i[0]
                            SearchSplit = SearchSplit + ";" + i[2]
                            DisplayColor = [0, 255, 0]
                            DisplayColor[0] = color_root[0] - len(dirpath)
                            DisplayColor[2] = color_root[2] - len(dirpath)
                            self.tree.SetItemData(self.ids[fullpath],wx.TreeItemData(SearchSplit)) 
                            self.tree.SetItemBackgroundColour(self.ids[fullpath], DisplayColor)
            self.tree.SetItemBackgroundColour(self.ids[root], [240, 255, 240])
            self.tree.Expand(self.ids[root])
            self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.SelectedTopTreeID, self.tree)                        
                                    
            #OPENING DISPLAY
            try:
                self.LOGO.Destroy()
            except:
                pass        
            self.png = wx.Image("Config/no-image-available.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
            self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
            imgsizer_v = wx.BoxSizer(wx.VERTICAL)
            imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
            self.page1.SetSizer(imgsizer_v)
            self.page1.Layout() 
            self.control.write("Resetting grid..." + "\n")
            self.control.write("Currently displaying: " + "SUMMARY" + "\n")
            self.myGrid.ClearGrid()
            opening_display_folder = self.input_stream + "/ExpressionOutput"
            list_contents = os.listdir(opening_display_folder)
            target_file = ""
            for file in list_contents:
                candidate = re.findall("SUMMARY", file)
                if len(candidate) > 0:
                    target_file = file
                    break 
            target_file = opening_display_folder + "/" + target_file
            opened_target_file = open(target_file, "r")
            opened_target_file_contents = []
            for line in opened_target_file:
                line = line.rstrip()
                line = line.split("\t")
                opened_target_file_contents.append((line))

            self.table_length = len(opened_target_file_contents)
            for cell in self.ColoredCellList:
                self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
            self.ColoredCellList = []

            x_count = 0
            for item_list in opened_target_file_contents:
                y_count = 0
                for item in item_list:
                    self.myGrid.SetCellValue(x_count, y_count, item)
                    if(x_count == 0):
                        self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                    y_count = y_count + 1
                x_count = x_count + 1
            self.myGrid.AutoSize()   
            for i in range(50):
                colsize = self.myGrid.GetColSize(i)
                if(colsize > 200):
                    self.myGrid.SetColSize(i, 200)
            self.page2.Layout()     
            
            self.nb.SetSelection(0)                                                                                                                                                                                                                                                                                  

    def OnOpenSingleFile(self, event):
        openFileDialog = wx.FileDialog(self, "", "", "", "", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     
            
        single_input_stream = openFileDialog.GetPath()
        self.control.write(str(single_input_stream) + "\n")
        if single_input_stream[-4:] == ".txt":
            self.myGrid.ClearGrid()
            self.DirFileTxt = single_input_stream
            table_file = open(self.DirFileTxt, "r")
            table_file_contents = []
            for line in table_file:
                line = line.rstrip()
                line = line.split("\t")
                table_file_contents.append((line))
            for cell in self.ColoredCellList:
                self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
            self.ColoredCellList = []
            x_count = 0
            for item_list in table_file_contents:
                y_count = 0
                for item in item_list:
                    self.myGrid.SetCellValue(x_count, y_count, item)
                    if(x_count == 0):
                        self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                    y_count = y_count + 1
                x_count = x_count + 1 
            self.myGrid.AutoSize() 
            for i in range(50):
                colsize = self.myGrid.GetColSize(i)
                if(colsize > 200):
                    self.myGrid.SetColSize(i, 200)
            self.page2.Layout()                 
            
        if single_input_stream[-4:] == ".png":
            self.myGrid.ClearGrid()
            try:
                self.LOGO.Destroy()
            except:
                pass    
            self.png = wx.Image(single_input_stream, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
            self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
            imgsizer_v = wx.BoxSizer(wx.VERTICAL)
            imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
            self.page1.SetSizer(imgsizer_v)
            self.page1.Layout()        
        

    def OnSave(self, event):
        saveFileDialog = wx.FileDialog(self, "", "", "", "", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)        
        if saveFileDialog.ShowModal() == wx.ID_CANCEL:
            return
             
                
    def OnSearch(self, event):
        popup = wx.TextEntryDialog(None, "Enter filter for results.", "Search", "Enter search here.")
        if popup.ShowModal()==wx.ID_OK:
            answer=popup.GetValue()
            popup.Destroy()
        else:
            popup.Destroy()
            return
    
    def TreeSearch(self, event):
        popup = wx.TextEntryDialog(None, "Search the browser tree for directories and files.", "Search", "Enter search here.")
        if popup.ShowModal()==wx.ID_OK:
            answer=popup.GetValue()
            self.control.write("K" + str(answer) + "\n")
            os.chdir(os.path.dirname(__file__))
            rootman = "Config/AvailableData"
            search_box = []
            found = ""
            for (dirpath, dirnames, filenames) in os.walk(rootman):
                for dirname in dirnames:
                    fullpath = dirpath + "/" + dirname
                    search_box.append(fullpath)
            self.control.write("Searching..." + "\n")
            for path in search_box:
                path2 = path.split("/")
                search_candidate = path2[-1]
                self.control.write(search_candidate + " " + str(answer) + "\n")
                if(str(answer) == search_candidate):
                    found = path
                    break
            self.control.write(found + "\n")
            tree_recreate = found.split("/")
            treepath = ""
            self.control.write(str(range(len(tree_recreate))) + "\n")
            tree_length = len(tree_recreate)
            last_tree_value = len(tree_recreate) - 1
            for i in range(tree_length):                
                self.control.write(str(i) + "\n")
                if(i == 0):
                    self.tree.Expand(self.ids[tree_recreate[i]])
                    treepath = treepath + tree_recreate[i]
                    self.control.write(treepath + "\n")
                if(i > 0 and i < last_tree_value):
                    treepath = treepath + "/" + tree_recreate[i]
                    self.control.write(treepath + "\n")
                    self.tree.Expand(self.ids[treepath])                    
                if(i == last_tree_value):
                    treepath = treepath + "/" + tree_recreate[i]
                    self.control.write(treepath + "\n")
                    self.tree.SelectItem(self.ids[treepath])                    
            popup.Destroy()
        else:
            popup.Destroy()
            return

    def GridSearch(self, event):
        popup = wx.TextEntryDialog(None, "Search the table.", "Search", "Enter search here.")
        if popup.ShowModal()==wx.ID_OK:
            PageDownFound = "False"
            match_count = 0
            answer=popup.GetValue()
            for cell in self.ColoredCellList:
                self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
            self.ColoredCellList = []
            y_range = range(self.table_length)
            x_range = range(100)
            y_count = 0
            for number in y_range:                                
                x_count = 0
                for number in x_range:        
                    cellvalue = self.myGrid.GetCellValue(y_count, x_count)
                    gridmatch = re.findall(answer, cellvalue)
                    if(len(gridmatch) > 0):
                        if(PageDownFound == "False"):
                            PageScrollY = y_count
                            PageScrollX = x_count
                            PageDownFound = "True"
                        match_count = match_count + 1
                        self.ColoredCellList.append((y_count, x_count))
                        self.myGrid.SetCellBackgroundColour(y_count, x_count, (255, 255, 125))
                    x_count = x_count + 1
                y_count = y_count + 1
            self.myGrid.MakeCellVisible(PageScrollY, PageScrollX)
            terminal_list = []
            for cell in self.ColoredCellList:
                newrow = cell[0] + 1
                newcolumn = cell[1] + 1
                terminal_list.append((newrow, newcolumn))
            self.control.write(str(match_count) + " matches found for " + answer + "\n")
            self.control.write("At positions (row, column): " + str(terminal_list) + "\n")
            popup.Destroy()
            self.nb.SetSelection(0)
        else:
            popup.Destroy()
            return

    def FilterTable(self, event):
        popup = wx.TextEntryDialog(None, "Filter the table.", "Search", "Enter filter phrase.")
        if popup.ShowModal()==wx.ID_OK:
            self.myGrid.ClearGrid()
            answer=popup.GetValue()            
            try:
                table_file = open(self.DirFileTxt, "r")
                table_file_contents = []
                count = 0
                for line in table_file:
                    line = line.rstrip()
                    regex_test = re.findall(answer.upper(), line.upper())
                    line = line.split("\t")
                    if(len(regex_test) > 0 or count == 0):
                        table_file_contents.append((line))
                        count = count + 1

                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1
                self.myGrid.AutoSize()      
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()             
            except:
                self.control.write("Unable to open txt." + "\n") 
            self.nb.SetSelection(0)   
            
    def SortTable(self, event):
        popup = wx.TextEntryDialog(None, "Sort the table.", "Sort", "Which column to sort from?")
        if popup.ShowModal()==wx.ID_OK:
            self.myGrid.ClearGrid()
            answer=popup.GetValue()
            answer = answer.upper()
            try:
                table_file = open(self.DirFileTxt, "r")
                table_file_contents = []
                pre_sort2 = []
                header = []
                t_c = 0
                column_clusters_flat = 0
                for line in table_file:
                    try:
                        line = line.rstrip()
                        line = line.split("\t")
                        if(t_c == 0):
                            header.append((line))
                            t_c = t_c + 1
                            continue
                        if(line[0] == "column_clusters-flat"):
                            header.append((line))
                            column_clusters_flat = 1
                            continue                            
                        line_sort_select = line[self.sortdict[answer]]
                        pre_sort1 = []
                        count = 0
                        for i in line:
                            if(count == 0):
                                try:
                                    pre_sort1.append(float(line_sort_select))
                                except:
                                    pre_sort1.append(line_sort_select)
                                pre_sort1.append(i)
                            if(count == self.sortdict[answer]):
                                count = count + 1
                                continue
                            if(count != 0):
                                pre_sort1.append(i)
                            count = count + 1
                        pre_sort2.append((pre_sort1))
                    except:
                        continue
                
                

                table_file_contents.append(header[0])
                if(column_clusters_flat == 1):
                    table_file_contents.append(header[1])
                pre_sort2 = sorted(pre_sort2, reverse = True)
                for line in pre_sort2:
                    try:
                        final_count1 = 0
                        final_count2 = 1
                        send_list = []
                        for item in line:
                            if(final_count1 == 0):
                                send_list.append(line[final_count2])
                            if(final_count1 == self.sortdict[answer]):
                                send_list.append(str(line[0]))
                            if(final_count1 != 0 and final_count1 != self.sortdict[answer]):
                                if(final_count1 < self.sortdict[answer]):
                                    send_list.append(line[final_count2])
                                if(final_count1 > self.sortdict[answer]):
                                    send_list.append(line[final_count1])                                    
                            final_count1 = final_count1 + 1
                            final_count2 = final_count2 + 1 
                        table_file_contents.append((send_list))   
                    except:
                        continue           
                           
                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1
                self.myGrid.AutoSizeRows(True)      
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()             
            except:
                self.control.write("Unable to sort." + "\n")   
            self.nb.SetSelection(0) 

    def FilterTablefromButton(self, event):
            self.myGrid.ClearGrid()
            answer = self.filterbox.GetLineText(0)            
            try:
                try:
                    self.myGrid.DeleteRows(100, self.AppendTotal, True)        
                except:
                    pass

                table_file_contents = []
                count = 0
                for line in open(self.DirFileTxt,'rU').xreadlines():
                    line = line.rstrip()
                    regex_test = re.findall(answer.upper(), line.upper())
                    line = line.split("\t")
                    if(len(regex_test) > 0 or count == 0):
                        table_file_contents.append((line))
                        count = count + 1

                self.table_length = len(table_file_contents)
                self.control.write("Table Length: " + str(self.table_length) + "\n")
                if(self.table_length > 100):
                    self.AppendTotal = self.table_length - 100
                    self.myGrid.AppendRows(self.AppendTotal, True)

                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1
                self.myGrid.AutoSize()      
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()             
            except:
                self.control.write("Unable to open txt." + "\n") 
            self.nb.SetSelection(0)   
                    
                                                            
    def SortTablefromButton(self, event):
            answer = self.sortbox.GetLineText(0)
            self.control.write(str(answer) + "\n")
            self.myGrid.ClearGrid()
            answer = answer.upper()
            try:
                table_file = open(self.DirFileTxt, "r")
                table_file_contents = []
                pre_sort2 = []
                header = []
                t_c = 0
                column_clusters_flat = 0
                for line in table_file:
                    try:
                        line = line.rstrip()
                        line = line.split("\t")
                        if(t_c == 0):
                            header.append((line))
                            t_c = t_c + 1
                            continue
                        if(line[0] == "column_clusters-flat"):
                            header.append((line))
                            column_clusters_flat = 1
                            continue                            
                        line_sort_select = line[self.sortdict[answer]]
                        pre_sort1 = []
                        count = 0
                        for i in line:
                            if(count == 0):
                                try:
                                    pre_sort1.append(float(line_sort_select))
                                except:
                                    pre_sort1.append(line_sort_select)
                                pre_sort1.append(i)
                            if(count == self.sortdict[answer]):
                                count = count + 1
                                continue
                            if(count != 0):
                                pre_sort1.append(i)
                            count = count + 1
                        pre_sort2.append((pre_sort1))
                    except:
                        continue
                
                table_file_contents.append(header[0])
                if(column_clusters_flat == 1):
                    table_file_contents.append(header[1])
                if(self.DescendingRadio.GetValue() == True):
                    pre_sort2 = sorted(pre_sort2, reverse = True)                    
                if(self.AscendingRadio.GetValue() == True):
                    pre_sort2 = sorted(pre_sort2) 
                for line in pre_sort2:
                    try:
                        final_count1 = 0
                        final_count2 = 1
                        send_list = []
                        for item in line:
                            if(final_count1 == 0):
                                send_list.append(line[final_count2])
                            if(final_count1 == self.sortdict[answer]):
                                send_list.append(str(line[0]))
                            if(final_count1 != 0 and final_count1 != self.sortdict[answer]):
                                if(final_count1 < self.sortdict[answer]):
                                    send_list.append(line[final_count2])
                                if(final_count1 > self.sortdict[answer]):
                                    send_list.append(line[final_count1])                                    
                            final_count1 = final_count1 + 1
                            final_count2 = final_count2 + 1 
                        table_file_contents.append((send_list))   
                    except:
                        continue           
                           
                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1
                self.myGrid.AutoSizeRows(True)      
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()             
            except:
                self.control.write("Unable to sort." + "\n")   
            self.nb.SetSelection(0) 

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    def SelectedTopTreeID(self, event):
        item = event.GetItem()
        itemObject = self.tree.GetItemData(item).GetData()        
        SearchObject = itemObject.split(";")
        SearchSuffix = SearchObject[0]
        SearchPath = SearchObject[1]
        SearchExtension = SearchObject[2]
        SearchExtension = SearchExtension.split(":")
        SearchExtension = SearchExtension[1:]
        SearchExtension = SearchExtension[0]
        SearchExtension = SearchExtension.split(",")
        #SearchExtension = SearchExtension[1:-1]
        #SearchExtension = SearchExtension.replace("'", "")
        #SearchExtension = SearchExtension.split(",")
        #for i in SearchExtension:
            #self.control.write(str(i) + "\n")
        
        #SELECTION IMPLEMENT
        ID_Strings = []
        self.TopSelectList = []
        self.TopID = SearchSuffix
        root = self.input_stream + "/" + SearchPath
        root_display = self.input_stream + "/" + SearchPath
        self.control.write("Current selection: " + itemObject + "\n")
        #self.control.write(root + "\n")
        root_contents = os.listdir(root)
        root_contents_display = os.listdir(root)
        for obj in root_contents:
            #self.control.write(obj + "\n")
            if(SearchSuffix != "*"):
                FindList = re.findall(SearchSuffix, obj)
                if(len(FindList) > 0):
                    self.TopSelectList.append(obj)
        self.browser2.DeleteAllItems()
        for filename in root_contents:              
            if(SearchSuffix != "*"):
                FindList2 = re.findall(SearchSuffix, filename)
                if(len(FindList2) > 0):
                    #self.control.write(filename + "\n")
                    display_name = filename[0:-4]
                    ID_Strings.append(display_name)
            else:
                #self.control.write(filename + "\n")
                if(filename[-4] == "."):
                    display_name = filename[0:-4]
                    ID_Strings.append(display_name)                
        ID_Strings = list(set(ID_Strings))
        change_path = os.path.dirname(__file__) + "/UseDir"
        shutil.rmtree("UseDir")
        os.mkdir("UseDir")
        #self.control.write(ID_Strings[0] + "\n")
        os.chdir(change_path)            
        for marker in ID_Strings:
            try:
                os.mkdir(marker)
            except:
                pass
        os.chdir(os.path.dirname(__file__))
        root = "UseDir"                      
        color_root2 = [223, 250, 223]
        self.ids2 = {root : self.browser2.AddRoot(root)}
        for (dirpath, dirnames, filenames) in os.walk(root):
            color_root2[0] = color_root2[0] - 1
            color_root2[1] = color_root2[1] - 0
            color_root2[2] = color_root2[2] - 1
            for dirname in dirnames:
                #self.control.write(str(SearchExtension) + "\n")
                Extensions = dirname + "|" + str(SearchExtension) + "|" + str(SearchPath)
                fullpath = os.path.join(dirpath, dirname)
                self.ids2[fullpath] = self.browser2.AppendItem(self.ids2[dirpath], dirname)                
                self.browser2.SetItemData(self.ids2[fullpath],wx.TreeItemData(Extensions))
                self.browser2.SetItemBackgroundColour(self.ids2[fullpath], color_root2)
        self.browser2.SetItemBackgroundColour(self.ids2[root], [240, 255, 240])
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.SelectedBottomTreeID, self.browser2)
        self.browser2.ExpandAll()
        
        #OPENING DISPLAY
        display_file_selected = ""
        TXT_FLAG = 0
        PNG_FLAG = 0
        if(root_display[-1] != "/"):
            root_display = root_display + "/"
        for possible in root_contents_display:
            total_filepath = root_display + possible 
            if(possible[-4:] == ".txt"):
                self.control.write("Displaying file: " + str(total_filepath) + "\n")
                display_file_selected = total_filepath
                break

        TXT_FLAG = 0
        PNG_FLAG = 0
        self.control.write(str(os.listdir(root)) + "\n")  
        self.control.write(str(SearchExtension) + "\n")  
        for i in SearchExtension:
            if(i == ".txt"):
                TXT_FLAG = 1
                self.control.write(str(i) + "\n")
            if(i == ".png"):
                PNG_FLAG = 1
                self.control.write(str(i) + "\n")
        if(root_display[-1] != "/"):
            root_display = root_display + "/"
        Pitch = os.listdir(root)
        PitchSelect = Pitch[0]
        self.control.write(str(PitchSelect) + " " + root_display + "\n")                
        self.DirFileTxt = root_display + PitchSelect + ".txt"       
        DirFilePng = root_display + PitchSelect + ".png"        
        self.myGrid.ClearGrid()
        title_name = PitchSelect
        try:
            self.LOGO.Destroy()
        except:
            pass                    
        try:
            self.PanelTitle.Destroy()
        except:
            pass    
        font = wx.Font(18, wx.DECORATIVE, wx.BOLD, wx.NORMAL)
        self.PanelTitle = wx.StaticText(self.panel, label=title_name, pos=(5, 7))
        self.PanelTitle.SetFont(font)


        if(TXT_FLAG == 1):
            try:
                self.myGrid.DeleteRows(100, self.AppendTotal, True)        
            except:
                pass
            try:
                #font = wx.Font(16, wx.DECORATIVE, wx.BOLD, wx.NORMAL)
                #self.PanelTitle = wx.StaticText(self.panel, label=title_name, pos=(210, 15))
                #self.PanelTitle.SetFont(font)
                table_file = open(self.DirFileTxt, "r")
                table_file_contents = []
                for line in table_file:
                    line = line.rstrip()
                    line = line.split("\t")
                    table_file_contents.append((line))

                self.table_length = len(table_file_contents)
                self.control.write("Table Length: " + str(self.table_length) + "\n")
                if(self.table_length > 100):
                    self.AppendTotal = self.table_length - 100
                    self.myGrid.AppendRows(self.AppendTotal, True)                
                
                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1 
                self.myGrid.AutoSize() 
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()                 
            except:
                TXT_FLAG = 0
                self.control.write("Unable to open txt." + "\n")

        if(PNG_FLAG == 1):
            try:                
                self.png = wx.Image(DirFilePng, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
                self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
                imgsizer_v = wx.BoxSizer(wx.VERTICAL)
                imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
                self.page1.SetSizer(imgsizer_v)
                self.page1.Layout()        
            except:
                PNG_FLAG = 0
                self.control.write("Unable to open png." + "\n")
        
        if(PNG_FLAG == 1 and TXT_FLAG == 0):
            self.nb.SetSelection(1)
        if(PNG_FLAG == 0 and TXT_FLAG == 1):
            self.nb.SetSelection(0)
        if(PNG_FLAG == 1 and TXT_FLAG == 1):
            self.nb.SetSelection(1)

    def SelectedBottomTreeID(self, event):
        item = event.GetItem()
        itemObject = self.browser2.GetItemData(item).GetData()
        Parameters = itemObject.split("|")
        file_extension = Parameters[1][1:-1]
        file_extension.replace("'", "")
        file_extension = file_extension.split(",")
        file_exts = []
        TXT_FLAG = 0
        PNG_FLAG = 0
        for i in file_extension:
            i = i.replace("'", "")
            i = i.replace(" ", "")
            file_exts.append(i)
            #self.control.write(str(i) + "\n")
        for i in file_exts:
            if(i == ".txt"):
                TXT_FLAG = 1
            if(i == ".png"):
                PNG_FLAG = 1
        DirPath = self.input_stream + "/" + Parameters[2]
        if(DirPath[-1] != "/"):
            DirPath = DirPath + "/"
        DirFile = DirPath + Parameters[0]
        self.control.write("Displaying file: " + DirFile + "\n")
        title_name = DirFile.split("/")
        title_name = title_name[-1]
        self.control.write(title_name + "\n")
        
        self.DirFile = DirFile
        self.DirFileTxt = DirFile + ".txt"
        DirFilePng = DirFile + ".png"
        self.myGrid.ClearGrid()
        try:
            self.LOGO.Destroy()
        except:
            pass    
        try:
            self.PanelTitle.Destroy()
        except:
            pass    
        font = wx.Font(18, wx.DECORATIVE, wx.BOLD, wx.NORMAL)
        self.PanelTitle = wx.StaticText(self.panel, label=title_name, pos=(5, 7))
        self.PanelTitle.SetFont(font)

      
        if(PNG_FLAG == 1):
            try:
                self.png = wx.Image(DirFilePng, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
                self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
                imgsizer_v = wx.BoxSizer(wx.VERTICAL)
                imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
                self.page1.SetSizer(imgsizer_v)
                self.page1.Layout()        
            except:
                PNG_FLAG = 0
                self.control.write("Unable to open png." + "\n")
        
        if(TXT_FLAG == 1):
            try:
                self.myGrid.DeleteRows(100, self.AppendTotal, True)        
            except:
                pass
            try:
                table_file = open(self.DirFileTxt, "r")
                table_file_contents = []
                for line in table_file:
                    line = line.rstrip()
                    line = line.split("\t")
                    table_file_contents.append((line))
                     
                self.table_length = len(table_file_contents)
                self.control.write("Table Length: " + str(self.table_length) + "\n")
                if(self.table_length > 100):
                    self.AppendTotal = self.table_length - 100
                    self.myGrid.AppendRows(self.AppendTotal, True)

                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                for item_list in table_file_contents:
                    y_count = 0
                    for item in item_list:
                        self.myGrid.SetCellValue(x_count, y_count, item)
                        if(x_count == 0):
                            self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                        y_count = y_count + 1
                    x_count = x_count + 1
                self.myGrid.AutoSize()      
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()             
            except:
                TXT_FLAG = 0
                self.control.write("Unable to open txt." + "\n")
            DATASET_FIND_FLAG = re.findall("DATASET", self.DirFileTxt)
            if(len(DATASET_FIND_FLAG) > 0):
                try:
                    table_file = open(self.DirFileTxt, "r")
                    table_file_contents = []
                    pre_sort2 = []
                    header = []
                    t_c = 0
                    column_clusters_flat = 0
                    answer = "AC"
                    for line in table_file:
                        try:
                            line = line.rstrip()
                            line = line.split("\t")
                            if(t_c == 0):
                                header.append((line))
                                t_c = t_c + 1
                                continue
                            if(line[0] == "column_clusters-flat"):
                                header.append((line))
                                column_clusters_flat = 1
                                continue                            
                            line_sort_select = line[self.sortdict[answer]]
                            pre_sort1 = []
                            count = 0
                            for i in line:
                                if(count == 0):
                                    try:
                                        pre_sort1.append(float(line_sort_select))
                                    except:
                                        pre_sort1.append(line_sort_select)
                                    pre_sort1.append(i)
                                if(count == self.sortdict[answer]):
                                    count = count + 1
                                    continue
                                if(count != 0):
                                    pre_sort1.append(i)
                                count = count + 1
                            pre_sort2.append((pre_sort1))
                        except:
                            continue
                    
                    
    
                    table_file_contents.append(header[0])
                    if(column_clusters_flat == 1):
                        table_file_contents.append(header[1])
                    pre_sort2 = sorted(pre_sort2)
                    for line in pre_sort2:
                        try:
                            final_count1 = 0
                            final_count2 = 1
                            send_list = []
                            for item in line:
                                if(final_count1 == 0):
                                    send_list.append(line[final_count2])
                                if(final_count1 == self.sortdict[answer]):
                                    send_list.append(str(line[0]))
                                if(final_count1 != 0 and final_count1 != self.sortdict[answer]):
                                    if(final_count1 < self.sortdict[answer]):
                                        send_list.append(line[final_count2])
                                    if(final_count1 > self.sortdict[answer]):
                                        send_list.append(line[final_count1])                                    
                                final_count1 = final_count1 + 1
                                final_count2 = final_count2 + 1 
                            table_file_contents.append((send_list))   
                        except:
                            continue           
                            
                    for cell in self.ColoredCellList:
                        self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                    self.ColoredCellList = []
                    x_count = 0
                    for item_list in table_file_contents:
                        y_count = 0
                        for item in item_list:
                            self.myGrid.SetCellValue(x_count, y_count, item)
                            if(x_count == 0):
                                self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                            y_count = y_count + 1
                        x_count = x_count + 1
                    self.myGrid.AutoSizeRows(True)      
                    for i in range(50):
                        colsize = self.myGrid.GetColSize(i)
                        if(colsize > 200):
                            self.myGrid.SetColSize(i, 200)
                    self.page2.Layout()             
                except:
                    self.control.write("Unable to sort." + "\n")   
                self.nb.SetSelection(0) 


        if(PNG_FLAG == 1 and TXT_FLAG == 0):
            self.nb.SetSelection(1)
        if(PNG_FLAG == 0 and TXT_FLAG == 1):
            self.nb.SetSelection(0)
        if(PNG_FLAG == 1 and TXT_FLAG == 1):
            self.nb.SetSelection(1)

    def ClearVisualPanel(self, event):
        popup = wx.MessageDialog(None, "Are you sure you want to clear the visual panel?", "Warning", wx.YES_NO)
        popup_answer = popup.ShowModal()
        if(popup_answer == 5103):
            try:
                self.png.Destroy()
                self.panel2.Layout()
            except:
                pass  
            try:    
                self.myGrid.Destroy()
                self.panel2.Layout()
            except:
                pass
            popup.Destroy()
            self.control.write("Visual panel cleared." + "\n")
        else:
            return                            

    def InteractiveDemo(self, event):
        os.chdir(back_dir)
        self.IntFileTxt.Clear()
        self.IntFileTxt.write(self.DirFile)
        RegExMatch = re.findall("exp.", self.DirFile)
        if(len(RegExMatch) == 0):
            InputFile = self.DirFile.replace("-3D", "")
            InputFile = InputFile.replace("-PCA", "")
            InputFile = InputFile.replace("DataPlots/Clustering-", "ExpressionOutput/Clustering/")
            input_file_dir= InputFile + ".txt"            
        else:        
            InputFile = self.DirFile.replace("-3D", "")
            InputFile = InputFile.replace("-PCA", "")
            InputFile = InputFile.replace("DataPlots/Clustering-", "ExpressionInput/")
            input_file_dir= InputFile + ".txt"
        include_labels= 'yes'
        pca_algorithm = 'SVD'
        transpose = False
        plotType = '3D'
        display = True
        UI.performPCA(input_file_dir, include_labels, pca_algorithm, transpose, None, plotType=plotType, display=display)                                                                                                                                                                              

class ImageFrame(wx.Frame):
    title = "Image"

    def __init__(self):
        wx.Frame.__init__(self, None, title=self.title)        

def remoteViewer(projectPath=None):
    global app
    app = wx.PySimpleApp(False)
    fr = Main(parent=None,id=1)
    fr.Show()
    app.MainLoop()

if __name__ == "__main__":
    app = wx.PySimpleApp(False)
    fr = Main(parent=None,id=1)
    fr.Show()
    app.MainLoop()
