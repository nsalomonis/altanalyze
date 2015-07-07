from numpy import arange, sin, pi
import os.path, sys, shutil
import wx.lib.scrolledpanel
import wx.grid as gridlib
import os
import Tkinter
import string, re
import wx
import subprocess
from wx import *
import __init__
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
#matplotlib.use('WXAgg')
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

currentDirectory = str(os.getcwd()) ### NS-91615 alternative to __file__
parentDirectory = os.path.abspath(os.path.join(currentDirectory, os.pardir))  ### NS-91615 gives the parent AltAnalyze directory
sys.path.insert(1,parentDirectory)  ### NS-91615 adds the AltAnalyze modules to the system path to import clustering and others
import UI

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
 
        self.heatmap_translation = {}
        self.heatmap_run = {}
                      
        self.SetBackgroundColour((230, 230, 230))
        self.species=''
        #PANELS & WIDGETS              
        self.panel = wx.Panel(self, id=2, pos=(200,0), size=(600,45), style=wx.RAISED_BORDER)       
        self.panel.SetBackgroundColour((110, 150, 250))
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
        
        self.panel_right2.SetMaxSize([300, 45])        
        
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

        self.nb = wx.Notebook(self.panel2, id=7829, style = wx.NB_BOTTOM)
        self.page1 = wx.ScrolledWindow(self.nb, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
        self.page1.SetScrollRate( 5, 5 )
        self.page2 = PageTwo(self.nb)
        self.page3 = PageThree(self.nb)
        self.page4 = PageFour(self.nb)
        self.nb.AddPage(self.page2, "Table")        
        self.nb.AddPage(self.page1, "PNG")
        self.nb.AddPage(self.page3, "Interactive")
        self.page3.SetBackgroundColour((218, 218, 218))
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer)
        self.page1.SetBackgroundColour("white")
        self.myGrid = gridlib.Grid(self.page2, id=1002)
        self.myGrid.CreateGrid(100, 100)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_CLICK, self.GridRightClick, id=1002)
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.GridRowColor, id=1002)
        self.HighlightedCells = []
        gridsizer = wx.BoxSizer(wx.VERTICAL)
        gridsizer.Add(self.myGrid)
        self.page2.SetSizer(gridsizer)
        self.page2.Layout()

        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.InteractiveTabChoose, id=7829)

	#INTERACTIVE PANEL LAYOUT
        #Pca Setup
        
        self.RunButton1 = wx.Button(self.page3, id=43, label="Run", pos=(275, 150), size=(120, 10))
        self.Bind(wx.EVT_BUTTON, self.InteractiveRun, id=43) 

        self.Divider1 = self.ln = wx.StaticLine(self.page3, pos=(5,100))
        self.ln.SetSize((415,10))
        
        IntTitleFont = wx.Font(15, wx.SWISS, wx.NORMAL, wx.BOLD)
        self.InteractiveTitle = wx.StaticText(self.page3, label="Main Dataset Parameters", pos=(10, 15))
        self.InteractiveDefaultMessage = wx.StaticText(self.page3, label="No interactive options available.", pos=(10, 45))
        self.InteractiveTitle.SetFont(IntTitleFont)        

        self.IntFileTxt = wx.TextCtrl(self.page3, id=43, pos=(105,45), size=(375,20))       
                        
        self.InteractiveFileLabel = wx.StaticText(self.page3, label="Selected File:", pos=(10, 45))                        
                                                                                
        self.Yes1Label = wx.StaticText(self.page3, label="Yes", pos=(305, 80))
        self.No1Label = wx.StaticText(self.page3, label="No", pos=(375, 80))

        self.D_3DLabel = wx.StaticText(self.page3, label="3D", pos=(305, 120))
        self.D_2DLabel = wx.StaticText(self.page3, label="2D", pos=(375, 120))
                        
        self.Yes1Radio = wx.RadioButton(self.page3, id=40, pos=(285, 83), size=(12, 12), style=wx.RB_GROUP)
        self.No1Radio = wx.RadioButton(self.page3, id=41, pos=(355, 83), size=(12, 12))
        self.Yes1Radio.SetValue(True)        
        
        self.D_3DRadio = wx.RadioButton(self.page3, id=46, pos=(285, 123), size=(12, 12), style=wx.RB_GROUP)
        self.D_2DRadio = wx.RadioButton(self.page3, id=47, pos=(355, 123), size=(12, 12))
        self.D_3DRadio.SetValue(True)

        self.Opt1Desc = wx.StaticText(self.page3, label="Display sample labels next to each object", pos=(10, 80))
        self.Opt2Desc = wx.StaticText(self.page3, label="Dimensions to display", pos=(10, 120))

        self.IntFileTxt.Hide()
        self.InteractiveFileLabel.Hide()
        self.Yes1Label.Hide()
        self.No1Label.Hide()
        self.D_3DLabel.Hide()
        self.D_2DLabel.Hide()
        self.Yes1Radio.Hide()
        self.No1Radio.Hide()
        self.D_3DRadio.Hide()
        self.D_2DRadio.Hide()
        self.Opt1Desc.Hide()
        self.Opt2Desc.Hide()
        self.RunButton1.Hide()
        self.Divider1.Hide()
        
        #self.Heatmap_CCMetricLabel = wx.StaticText(self.page3, label="Select the column clustering metric", pos=(150, 75))
        #self.Heatmap_CCMethodLabel = wx.StaticText(self.page3, label="Select the column clustering method", pos=(150, 100))

        #self.Heatmap_CCMetricUI = wx.ComboBox(self.page3, 500, '', (10, 75), (120,25), phasesList, wx.CB_DROPDOWN)
        #self.Heatmap_CCMethodUI = wx.ComboBox(self.page3, 501, '', (10, 100), (120,25), phasesList, wx.CB_DROPDOWN)

        #TERMINAL SETUP

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
        self.png = wx.Image("logo.gif", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        #Open Button        
        ButtonMan = wx.Button(self.panel_left, id=1001, label="Open Project", pos=(0,0), size=(100,100)) 
        self.Bind(wx.EVT_BUTTON, self.OnOpen, id=1001)
        OpenSizer = wx.BoxSizer(wx.HORIZONTAL)
        OpenSizer.Add(ButtonMan, 1, wx.EXPAND)
        self.panel_left.SetSizer(OpenSizer)

        #STATUS BAR CREATE             
        status = self.CreateStatusBar()
        menubar = wx.MenuBar()
                
        file = wx.Menu()
        edit = wx.Menu()
        view = wx.Menu()
        search = wx.Menu()
        filter_table = wx.Menu()
        help_menu = wx.Menu()
        
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
        
        filter_table.Append(116, 'Filter', '')
        filter_table.Append(117, 'Sort', '')

        help_menu.AppendSeparator()
        help_menu.Append(139, 'Help', '')
        help_menu.Append(140, 'About', '')

        menubar.Append(file, "File")
        menubar.Append(edit, "Edit")
        menubar.Append(view, "View")
        menubar.Append(search, "Search")
        menubar.Append(filter_table, "Table")
        menubar.Append(help_menu, "Help")      
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
        self.Bind(wx.EVT_MENU, self.OnAbout, id=140)
        self.Bind(wx.EVT_MENU, self.OnHelp, id=139)      
        
        self.Layout()           
                                                       
    def OnQuit(self, event):
        popup = wx.MessageDialog(None, "Are you sure you want to quit?", "Warning", wx.YES_NO)
        popup_answer = popup.ShowModal()
        print popup_answer
        if(popup_answer == 5103):
            self.Close()
        else:
            return

    def GridRowColor(self, event):
        if len(self.HighlightedCells) > 0:
            for i in self.HighlightedCells:
                self.myGrid.SetCellBackgroundColour(i[0], i[1], (255, 255, 255)) 
            self.HighlightedCells = []
        self.GridRowEvent = event.GetRow() 
        for i in range(50):
            self.myGrid.SetCellBackgroundColour(self.GridRowEvent, i, (235, 255, 255))   
            self.HighlightedCells.append((self.GridRowEvent, i))  
          
            
    def GridRightClick(self, event):
        """
        Create and show a Context Menu
        """
        self.GridRowEvent = event.GetRow()       
        
        # only do this part the first time so the events are only bound once 
        if not hasattr(self, "popupID2"):
            self.popupID1 = wx.NewId()
            self.popupID2 = wx.NewId()
            self.Bind(wx.EVT_MENU, self.CellGraph, id=self.popupID1)
            self.Bind(wx.EVT_MENU, self.PrintGraphVariables, id=self.popupID2)
 
        # build the menu
        menu = wx.Menu()
        itemOne = menu.Append(self.popupID1, "Plot")
        itemTwo = menu.Append(self.popupID2, "Print Variables")
 
        # show the popup menu
        self.PopupMenu(menu)
        menu.Destroy()
                                    
    def CellGraph(self, event):
        Wikipathway_Flag = 0
        Protein_Flag = 0
        VarGridSet = []
        for i in range(100):
            p = self.myGrid.GetCellValue(0, i)
            VarGridSet.append(p)
        for i in VarGridSet:
            y = re.findall("WikiPathways", i)
            if len(y) > 0:
                Wikipathway_Flag = 1
                break
        if Wikipathway_Flag == 0:
            for i in VarGridSet:
                y = re.findall("Select Protein Classes", i)
                if len(y) > 0:
                    Protein_Flag = 1
                    break
        if Protein_Flag == 1:
            VariableBox = []
            for i in range(len(VarGridSet)):
                y = re.findall("avg", VarGridSet[i])
                if(len(y) > 0):
                    VariableBox.append(i)
        if Wikipathway_Flag == 1:
            VariableBox = []
            for i in range(len(VarGridSet)):
                y = re.findall("avg", VarGridSet[i])
                if(len(y) > 0):
                    VariableBox.append(i)
        q_barrel = []
        for i in VariableBox:
            q_box = []
            q = i
            for p in range(500):                
                if(q < 0):
                    break
                q = q - 1   
                FLAG_log_fold = re.findall("log_fold",VarGridSet[q])
                FLAG_adjp = re.findall("adjp",VarGridSet[q])
                FLAG_rawp = re.findall("rawp",VarGridSet[q])
                FLAG_wiki = re.findall("Wiki",VarGridSet[q])
                FLAG_pc = re.findall("Protein Classes",VarGridSet[q])
                FLAG_avg = re.findall("avg",VarGridSet[q])
                if(len(FLAG_log_fold) > 0 or len(FLAG_adjp) > 0 or len(FLAG_rawp) > 0 or len(FLAG_wiki) > 0 or len(FLAG_pc) > 0 or len(FLAG_avg) > 0):
                    break
                q_box.append(q)
            q_barrel.append((q_box))
        Values_List = []
        HeaderList = []
        TitleList = self.myGrid.GetCellValue(self.GridRowEvent, 0)
        for i in VariableBox:
            HeaderList.append(self.myGrid.GetCellValue(0, i))
        for box in q_barrel:
            output_box = []
            for value in box:
                output_var = self.myGrid.GetCellValue(self.GridRowEvent, value)
                output_box.append(float(output_var))
            Values_List.append((output_box))            
        self.control.write("Plotting values from: " + str(self.myGrid.GetCellValue(self.GridRowEvent, 0)) + "\n")
        Output_Values_List = []
        Output_std_err = []
        for box in Values_List:
            T = 0
            for item in box:
                T = T + item
            output_item = T / float(len(box))
            Output_Values_List.append(output_item)
        
        for box in Values_List:
            box_std = np.std(box)
            box_power = np.power((len(box)), 0.5)
            std_err = box_std / float(box_power)
            Output_std_err.append(std_err)                
            
        n_groups = len(Output_Values_List)
        
        means_men = Output_Values_List
        
        fig, ax = plt.subplots()
        
        index = np.arange(n_groups)
        bar_width = 0.35
        pos = bar_width / float(2)
        
        opacity = 0.4
        error_config = {'ecolor': '0.3'}
        
        rects1 = plt.bar((index + pos), Output_Values_List, bar_width,
                        alpha=opacity,
                        color='b',
                        yerr=Output_std_err,
                        label="")
        
        #plt.title(self.myGrid.GetCellValue(self.GridRowEvent, 2))
        plt.title(TitleList)
        plt.xticks(index + bar_width, HeaderList)
        plt.legend()
        
        plt.tight_layout()
        plt.show()

    def PrintGraphVariables(self, event):
        Wikipathway_Flag = 0
        Protein_Flag = 0
        VarGridSet = []
        for i in range(100):
            p = self.myGrid.GetCellValue(0, i)
            VarGridSet.append(p)
        for i in VarGridSet:
            y = re.findall("WikiPathways", i)
            if len(y) > 0:
                Wikipathway_Flag = 1
                break
        if Wikipathway_Flag == 0:
            for i in VarGridSet:
                y = re.findall("Select Protein Classes", i)
                if len(y) > 0:
                    Protein_Flag = 1
                    break
        if Protein_Flag == 1:
            VariableBox = []
            for i in range(len(VarGridSet)):
                y = re.findall("avg", VarGridSet[i])
                if(len(y) > 0):
                    VariableBox.append(i)
        if Wikipathway_Flag == 1:
            VariableBox = []
            for i in range(len(VarGridSet)):
                y = re.findall("avg", VarGridSet[i])
                if(len(y) > 0):
                    VariableBox.append(i)
        q_barrel = []
        for i in VariableBox:
            q_box = []
            q = i
            for p in range(500):                
                if(q < 0):
                    break
                q = q - 1   
                FLAG_log_fold = re.findall("log_fold",VarGridSet[q])
                FLAG_adjp = re.findall("adjp",VarGridSet[q])
                FLAG_rawp = re.findall("rawp",VarGridSet[q])
                FLAG_wiki = re.findall("Wiki",VarGridSet[q])
                FLAG_pc = re.findall("Protein Classes",VarGridSet[q])
                FLAG_avg = re.findall("avg",VarGridSet[q])
                if(len(FLAG_log_fold) > 0 or len(FLAG_adjp) > 0 or len(FLAG_rawp) > 0 or len(FLAG_wiki) > 0 or len(FLAG_pc) > 0 or len(FLAG_avg) > 0):
                    break
                q_box.append(q)
            q_barrel.append((q_box))
        self.control.write("Selected Row: " + str(self.myGrid.GetCellValue(self.GridRowEvent, 0)) + "\n")
        self.control.write("Selected Columns: " + str(q_barrel) + "\n")
        Values_List = []
        HeaderList = []
        for i in VariableBox:
            HeaderList.append(self.myGrid.GetCellValue(0, i))
        for box in q_barrel:
            output_box = []
            for value in box:
                output_var = self.myGrid.GetCellValue(self.GridRowEvent, value)
                output_box.append(float(output_var))
            Values_List.append((output_box))            
        self.control.write("Selected Values: " + str(Values_List) + "\n")
                
                                
    def InteractiveTabChoose(self, event):
        try:
            PCA_RegEx = re.findall("PCA", self.DirFile)
            Heatmap_RegEx = re.findall("hierarchical", self.DirFile)
            if(self.nb.GetSelection() == 2):
                if(len(PCA_RegEx) > 0 or len(Heatmap_RegEx) > 0):
                    self.InteractiveRun(event)
        except:
            pass          
        
    def getDatasetVariables(self):
        
        for file in os.listdir(self.input_stream):
            if 'AltAnalyze_report' in file and '.log' in file:
                log_file = self.input_stream+'/'+file
                log_contents = open(log_file, "rU")
                species = '	species: '
                platform = '	method: '
                for line in log_contents:
                    if species in line:
                        self.species = string.split(line,species)[1]
                    if platform in line:
                        self.platform = string.split(line,platform)[1]

    def OnOpen(self, event):
        openFileDialog = wx.DirDialog(None, "Choose project", "", wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)  
        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return 
        self.input_stream = openFileDialog.GetPath()
        if (len(self.input_stream) > 0):
            
            self.getDatasetVariables()
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
                            line = string.replace(line,'"','')
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
            
            try:
                shutil.rmtree("AvailableData")
            except:
                pass
            
            for i in self.SearchArrayFiltered:
                    AvailablePath = "Available" + i[3]
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
            root = 'AvailableData'
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
                        SearchRoot = "Available" + i[3]
                        if(SearchRoot == fullpath):
                            SearchSplit = i[1].split(":")
                            SearchSplit = SearchSplit[1]
                            SearchSplit = SearchSplit + ";" + i[0]
                            SearchSplit = SearchSplit + ";" + i[2]
                            DisplayColor = [130, 170, 250]
                            self.tree.SetItemData(self.ids[fullpath],wx.TreeItemData(SearchSplit)) 
                            self.tree.SetItemBackgroundColour(self.ids[fullpath], DisplayColor)
            self.tree.SetItemBackgroundColour(self.ids[root], [100, 140, 240])
            self.tree.Expand(self.ids[root])
            self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.SelectedTopTreeID, self.tree)                        
                                    
            #OPENING DISPLAY
            try:
                self.LOGO.Destroy()
            except:
                pass        
            self.png = wx.Image("no-image-available.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
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
                line = string.replace(line,'"','')
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
            self.DirFile = single_input_stream
            table_file = open(self.DirFileTxt, "r")
            table_file_contents = []
            for line in table_file:
                line = line.rstrip()
                line = string.replace(line,'"','')
                line = line.split("\t")
                if(len(table_file_contents) >= 5000):
                    break
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
            os.chdir(currentDirectory) ### NS-91615 alternative to __file__
            rootman = "AvailableData"
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
            if(self.table_length > 5100):
                y_range = range(5100)
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
                    line = string.replace(line,'"','')
                    regex_test = re.findall(answer.upper(), line.upper())
                    line = line.split("\t")
                    if(len(regex_test) > 0 or count == 0):
                        if(len(table_file_contents) >= 5100):
                            break
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
                        line = string.replace(line,'"','')
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
                        if(len(table_file_contents) >= 5100):
                            break
                        table_file_contents.append((send_list))   
                    except:
                        continue           
                
                n_table_file_contents = []
                if(answer.upper() == "A"):
                    for i in range(len(table_file_contents)):
                        if(i == 0):
                            max_length = len(table_file_contents[i])  
                        if(max_length < len(table_file_contents[i])):
                            n_l = table_file_contents[i][2:]                  
                        else:
                            n_l = table_file_contents[i]
                        n_table_file_contents.append((n_l)) 
                    
                    table_file_contents = n_table_file_contents                      
                                                                                 
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
                    line = string.replace(line,'"','')
                    regex_test = re.findall(answer.upper(), line.upper())
                    line = line.split("\t")
                    if(len(regex_test) > 0 or count == 0):
                        if(len(table_file_contents) >= 5100):
                            break
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
                        line = string.replace(line,'"','')
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
                        if(len(table_file_contents) >= 5100):
                            break
                        table_file_contents.append((send_list))   
                    except:
                        continue           

                n_table_file_contents = []
                if(answer.upper() == "A"):
                    for i in range(len(table_file_contents)):
                        if(i == 0):
                            max_length = len(table_file_contents[i])  
                        if(max_length < len(table_file_contents[i])):
                            n_l = table_file_contents[i][2:]                  
                        else:
                            n_l = table_file_contents[i]
                        n_table_file_contents.append((n_l)) 
                    
                    table_file_contents = n_table_file_contents                      
                           
                                                                                 
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
        
        #SELECTION IMPLEMENT
        ID_Strings = []
        self.TopSelectList = []
        self.TopID = SearchSuffix
        root = self.input_stream + "/" + SearchPath
        root_display = self.input_stream + "/" + SearchPath
        root_contents = os.listdir(root)
        root_contents_display = os.listdir(root)
        for obj in root_contents:
            if(SearchSuffix != "*"):
                FindList = re.findall(SearchSuffix, obj)
                if(len(FindList) > 0):
                    self.TopSelectList.append(obj)
        self.browser2.DeleteAllItems()
        for filename in root_contents:              
            if(SearchSuffix != "*"):
                FindList2 = re.findall(SearchSuffix, filename)
                if(len(FindList2) > 0):
                    display_name = filename[0:-4]
                    ID_Strings.append(display_name)
            else:
                if(filename[-4] == "."):
                    display_name = filename[0:-4]
                    ID_Strings.append(display_name)                
        ID_Strings = list(set(ID_Strings))
        change_path = currentDirectory + "/UseDir" ### NS-91615 alternative to __file__
        shutil.rmtree("UseDir")
        os.mkdir("UseDir")
        #self.control.write(ID_Strings[0] + "\n")
        os.chdir(change_path)            
        for marker in ID_Strings:
            try:
                os.mkdir(marker)
            except:
                pass
        os.chdir(currentDirectory) ### NS-91615 alternative to __file__
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
                T = re.findall("DATASET", fullpath)
                if(len(T) > 0):
                    self.browser2.SetItemBackgroundColour(self.ids2[fullpath], [250, 100, 100])                    
                else:
                    self.browser2.SetItemBackgroundColour(self.ids2[fullpath], [130, 170, 250])
        self.browser2.SetItemBackgroundColour(self.ids2[root], [110, 150, 250])
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
        #self.control.write(str(os.listdir(root)) + "\n")  
        #self.control.write(str(SearchExtension) + "\n")  
        for i in SearchExtension:
            if(i == ".txt"):
                TXT_FLAG = 1
                #self.control.write(str(i) + "\n")
            if(i == ".png"):
                PNG_FLAG = 1
                #self.control.write(str(i) + "\n")
        if(root_display[-1] != "/"):
            root_display = root_display + "/"
        Pitch = os.listdir(root)
        PitchSelect = Pitch[0]
        #self.control.write(str(PitchSelect) + " " + root_display + "\n")                
        self.DirFile = root_display + PitchSelect
        self.IntFileTxt.Clear()
        self.IntFileTxt.write(self.DirFile)
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
                column_lengths = []
                for line in table_file:
                    line = line.rstrip()
                    line = string.replace(line,'"','')
                    line = line.split("\t")
                    column_lengths.append(len(line))
                    table_file_contents.append((line))

                self.max_column_length = max(column_lengths)
                self.table_length = len(table_file_contents)
                if(self.table_length > 100 and self.table_length < 5000):
                    self.AppendTotal = self.table_length - 100
                    self.myGrid.AppendRows(self.AppendTotal, True)
                if(self.table_length >= 5000):
                    self.AppendTotal = 5000
                    self.myGrid.AppendRows(self.AppendTotal, True)
                
                for cell in self.ColoredCellList:
                    self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                self.ColoredCellList = []
                x_count = 0
                try:
                    for item_list in table_file_contents:
                        y_count = 0
                        for item in item_list:
                            self.myGrid.SetCellValue(x_count, y_count, item)
                            if(x_count == 0):
                                self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                            y_count = y_count + 1
                        x_count = x_count + 1 
                except:
                    pass
                self.myGrid.AutoSize() 
                for i in range(50):
                    colsize = self.myGrid.GetColSize(i)
                    if(colsize > 200):
                        self.myGrid.SetColSize(i, 200)
                self.page2.Layout()                 
            except:
                TXT_FLAG = 0
                self.control.write("Unable to open txt." + "\n")
        
        try:
            self.myGrid.AutoSize()
            for i in range(50):
                colsize = self.myGrid.GetColSize(i)
                if(colsize > 200):
                    self.myGrid.SetColSize(i, 200)
            self.page2.Layout()
        except:
            pass

        if(PNG_FLAG == 1):
            try:                
                open(DirFilePng, "r")
                self.png = wx.Image(DirFilePng, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
                self.LOGO = wx.StaticBitmap(self.page1, -1, self.png, (0,0), (self.png.GetWidth(), self.png.GetHeight()), style=wx.ALIGN_CENTER)
                imgsizer_v = wx.BoxSizer(wx.VERTICAL)
                imgsizer_v.Add(self.LOGO, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
                self.page1.SetSizer(imgsizer_v)
                self.page1.Layout()        
            except:
                PNG_FLAG = 0
                self.control.write("Unable to open png." + "\n")

        try:
            self.root_widget_id = 500
            self.root_widget_text = 550
            for i in range(self.root_widget_id, self.root_widget_end):
                self.heatmap_ids[i].Destroy()
            for i in range(self.root_widget_text, self.rwtend):
                self.heatmap_ids[i].Destroy()
            self.RunButton2.Destroy()
        except:
            pass
        
        self.InteractivePanelUpdate(event)                                                                
                                                                                                                                                                                                        
        if(PNG_FLAG == 1 and TXT_FLAG == 0):
            self.nb.SetSelection(1)
            self.Layout()
            self.page1.Layout()
        if(PNG_FLAG == 0 and TXT_FLAG == 1):
            self.nb.SetSelection(0)
        if(PNG_FLAG == 1 and TXT_FLAG == 1):
            self.nb.SetSelection(1)
            self.Layout()
            self.page1.Layout()

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
        
        self.DirFile = DirFile
        self.IntFileTxt.Clear()
        self.IntFileTxt.write(self.DirFile)
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
                open(DirFilePng, "r")
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
                column_lengths = []
                for line in table_file:
                    line = line.rstrip()
                    line = string.replace(line,'"','')
                    line = line.split("\t")
                    column_lengths.append(len(line))
                    table_file_contents.append((line))
                
                self.max_column_length = max(column_lengths)     
                self.table_length = len(table_file_contents)
                if(self.table_length > 100 and self.table_length < 5000):
                    self.AppendTotal = self.table_length - 100
                    self.myGrid.AppendRows(self.AppendTotal, True)
                if(self.table_length >= 5000):
                    self.AppendTotal = 5000
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
                    table_file = open(self.DirFileTxt, "rU")
                    table_file_contents = []
                    pre_sort2 = []
                    header = []
                    t_c = 0
                    column_clusters_flat = 0
                    answer = "AC"
                    for line in table_file:
                        try:
                            line = line.rstrip()
                            line = string.replace(line,'"','')
                            line = line.split("\t")
                            if(t_c == 0):
                                header.append((line))
                                t_c = t_c + 1
                                index=0
                                for i in line:
                                    if 'ANOVA-rawp' in i: answer = index
                                    index+=1
                                continue
                            if(line[0] == "column_clusters-flat"):
                                header.append((line))
                                column_clusters_flat = 1
                                continue                            
                            line_sort_select = line[answer]
                            pre_sort1 = []
                            count = 0
                            for i in line:
                                if(count == 0):
                                    try:
                                        pre_sort1.append(float(line_sort_select))
                                    except:
                                        pre_sort1.append(line_sort_select)
                                    pre_sort1.append(i)
                                if(count == answer):
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
                                if(final_count1 == answer):
                                    send_list.append(str(line[0]))
                                if(final_count1 != 0 and final_count1 != answer):
                                    if(final_count1 < answer):
                                        send_list.append(line[final_count2])
                                    if(final_count1 > answer):
                                        send_list.append(line[final_count1])                                    
                                final_count1 = final_count1 + 1
                                final_count2 = final_count2 + 1 
                            if(len(table_file_contents) >= 5100):
                                break
                            table_file_contents.append((send_list))   
                        except:
                            continue           
                            
                    for cell in self.ColoredCellList:
                        self.myGrid.SetCellBackgroundColour(cell[0], cell[1], wx.WHITE)
                    self.ColoredCellList = []
                    x_count = 0
                    try:
                        for item_list in table_file_contents:
                            y_count = 0
                            for item in item_list:
                                self.myGrid.SetCellValue(x_count, y_count, item)
                                if(x_count == 0):
                                    self.myGrid.SetCellFont(x_count, y_count, wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                                y_count = y_count + 1
                            x_count = x_count + 1
                    except:
                        pass
                    self.myGrid.AutoSizeRows(True)      
                    for i in range(50):
                        colsize = self.myGrid.GetColSize(i)
                        if(colsize > 200):
                            self.myGrid.SetColSize(i, 200)
                    self.page2.Layout()             
                except:
                    self.control.write("Unable to sort." + "\n")   
                self.nb.SetSelection(0) 
        
        self.InteractivePanelUpdate(event)
        try:
            self.myGrid.AutoSize()
            for i in range(50):
                colsize = self.myGrid.GetColSize(i)
                if(colsize > 200):
                    self.myGrid.SetColSize(i, 200)
            self.page2.Layout()
        except:
            pass
                            
        if(PNG_FLAG == 1 and TXT_FLAG == 0):
            self.nb.SetSelection(1)
            self.Layout()
            self.page1.Layout()
        if(PNG_FLAG == 0 and TXT_FLAG == 1):
            self.nb.SetSelection(0)
        if(PNG_FLAG == 1 and TXT_FLAG == 1):
            self.nb.SetSelection(1)
            self.Layout()
            self.page1.Layout()

    def InteractivePanelUpdate(self, event):
        self.IntFileTxt.Hide()
        self.InteractiveFileLabel.Hide()
        self.Yes1Label.Hide()
        self.No1Label.Hide()
        self.D_3DLabel.Hide()
        self.D_2DLabel.Hide()
        self.Yes1Radio.Hide()
        self.No1Radio.Hide()
        self.D_3DRadio.Hide()
        self.D_2DRadio.Hide()
        self.Opt1Desc.Hide()
        self.Opt2Desc.Hide()
        self.RunButton1.Hide()
        self.Divider1.Hide()        
        
        self.InteractiveDefaultMessage.Hide()        
        
        try:
            self.root_widget_id = 500
            self.root_widget_text = 550
            for i in range(self.root_widget_id, self.root_widget_end):
                self.heatmap_ids[i].Destroy()
            for i in range(self.root_widget_text, self.rwtend):
                self.heatmap_ids[i].Destroy()
            self.RunButton2.Destroy()
        except:
            pass

        PCA_RegEx = re.findall("PCA", self.DirFile)
        if(len(PCA_RegEx) > 0):            
            self.IntFileTxt.Show()
            self.InteractiveFileLabel.Show()
            self.Yes1Label.Show()
            self.No1Label.Show()
            self.D_3DLabel.Show()
            self.D_2DLabel.Show()
            self.Yes1Radio.Show()
            self.No1Radio.Show()
            self.D_3DRadio.Show()
            self.D_2DRadio.Show()
            self.Opt1Desc.Show()
            self.Opt2Desc.Show()
            self.RunButton1.Show()
            self.Divider1.Show()                                

        Heatmap_RegEx = re.findall("hierarchical", self.DirFile)
        if(len(Heatmap_RegEx) > 0):
            #Heatmap Setup
            os.chdir(parentDirectory)
            options_open = open(currentDirectory+"/options.txt", "rU")
            heatmap_array = []
            self.heatmap_ids = {}
            self.heatmap_translation = {}
            
            supported_geneset_types = UI.getSupportedGeneSetTypes(self.species,'gene-mapp')
            supported_geneset_types += UI.getSupportedGeneSetTypes(self.species,'gene-go')
            print supported_geneset_types
            for line in options_open:
                line = line.split("\t")
                variable_name,displayed_title,display_object,group,notes,description,global_default,options = line[:8]
                options = string.split(options,'|')
                if(group == "heatmap"):
                    if(display_object == "file"):
                        continue                    
                    od = UI.OptionData(variable_name,displayed_title,display_object,notes,options,global_default)
                    od.setDefaultOption(global_default)
                    if group == 'ClusterGOElite':
                        od.setArrayOptions(['None Selected','all']+supported_geneset_types)
                    if group == 'GeneSetSelection':
                        od.setArrayOptions(['None Selected']+supported_geneset_types)
                    if group == 'PathwaySelection':
                        od.setArrayOptions(['None Selected']+supported_geneset_types)
                    heatmap_array.append(od)
                    #heatmap_array.append((line[1], line[2], line[7], line[6]))           
            os.chdir(currentDirectory)
            
            root_widget_y_pos = 45 
            self.root_widget_id = 500
            self.root_widget_text = 550                 
            for od in heatmap_array:
                #od.VariableName()
                id = wx.NewId()
                self.heatmap_translation[od.VariableName()] = self.root_widget_id
                self.heatmap_ids[self.root_widget_text] = wx.StaticText(self.page3, self.root_widget_text, label=od.Display(), pos=(150, root_widget_y_pos))
                if(od.DisplayObject() == "comboBox" or od.DisplayObject() == "multiple-comboBox"):
                    self.heatmap_ids[self.root_widget_id] = wx.ComboBox(self.page3, self.root_widget_id, od.DefaultOption(), (10, root_widget_y_pos), (120,25), od.Options(), wx.CB_DROPDOWN)
                else:
                    self.heatmap_ids[self.root_widget_id] = wx.TextCtrl(self.page3, self.root_widget_id, "", (10, root_widget_y_pos), (120,25))
                self.root_widget_id = self.root_widget_id + 1
                self.root_widget_text = self.root_widget_text + 1
                root_widget_y_pos = root_widget_y_pos + 25
            
            self.rwtend = self.root_widget_text
            self.root_widget_end = self.root_widget_id       
    
            self.RunButton2 = wx.Button(self.page3, id=599, label="Run", pos=(175, (self.root_widget_end + 10)), size=(120, 10))
            self.Bind(wx.EVT_BUTTON, self.InteractiveRun, id=599) 
        
        if(len(PCA_RegEx) == 0 and len(Heatmap_RegEx) == 0):
            self.InteractiveDefaultMessage.Show() 

    def ClearVisualPanel(self, event):
        popup = wx.MessageDialog(None, "Are you sure you want to clear the visual panel?", "Warning", wx.YES_NO)
        popup_answer = popup.ShowModal()
        if(popup_answer == 5103):
            try:
                self.LOGO.Destroy()
                self.panel2.Layout()
            except:
                pass  
            try:    
                self.myGrid.ClearGrid()
                self.panel2.Layout()
            except:
                pass
            popup.Destroy()
            self.control.write("Visual panel cleared." + "\n")
        else:
            return                            

    def InteractiveRun(self, event):                        
            RegExHeat = re.findall("hierarchical", self.DirFile)
            
            if(len(RegExHeat) > 0):
                #print self.heatmap_translation
                for VariableName in self.heatmap_translation:
                    self.control.write(str(self.heatmap_ids[self.heatmap_translation[VariableName]].GetValue()) + " " + str(VariableName) + " " + str(self.heatmap_ids[self.heatmap_translation[VariableName]]) + "\n")
                    self.heatmap_translation[VariableName] = self.heatmap_ids[self.heatmap_translation[VariableName]].GetValue()
                self.control.write(self.DirFile + "\n")
                input_file_dir = self.DirFile + ".txt"
                column_metric = self.heatmap_translation['column_metric']; #self.control.write(column_metric + "\n")
                column_method = self.heatmap_translation['column_method']; #self.control.write(column_method + "\n")
                row_metric = self.heatmap_translation['row_metric']; #self.control.write(row_metric + "\n")
                row_method = self.heatmap_translation['row_method']; #self.control.write(row_method+ "\n") 
                color_gradient = self.heatmap_translation['color_selection']; #self.control.write(color_gradient + "\n")
                cluster_rows = self.heatmap_translation['cluster_rows']; #self.control.write(cluster_rows + "\n")
                cluster_columns = self.heatmap_translation['cluster_columns']; #self.control.write(cluster_columns + "\n")
                normalization = self.heatmap_translation['normalization']; #self.control.write(normalization + "\n")
                contrast = self.heatmap_translation['contrast']; #self.control.write(contrast + "\n")
                transpose = self.heatmap_translation['transpose']; #self.control.write(transpose + "\n")
                GeneSetSelection = self.heatmap_translation['GeneSetSelection']; #self.control.write(GeneSetSelection + "\n")
                PathwaySelection = self.heatmap_translation['PathwaySelection']; #self.control.write(PathwaySelection + "\n")
                OntologyID = self.heatmap_translation['OntologyID']; #self.control.write(OntologyID + "\n")
                GeneSelection = self.heatmap_translation['GeneSelection']; #self.control.write(GeneSelection + "\n")
                JustShowTheseIDs = self.heatmap_translation['JustShowTheseIDs']; #self.control.write(JustShowTheseIDs + "\n")
                CorrelationCutoff = self.heatmap_translation['CorrelationCutoff']; #self.control.write(CorrelationCutoff + "\n")
                HeatmapAdvanced = self.heatmap_translation['HeatmapAdvanced']; #self.control.write(HeatmapAdvanced + "\n")
                ClusterGOElite = self.heatmap_translation['ClusterGOElite']; #self.control.write(ClusterGOElite + "\n")
                heatmapGeneSets = self.heatmap_translation['heatmapGeneSets']; #self.control.write(heatmapGeneSets + "\n")
                if cluster_rows == 'no': row_method = None
                if cluster_columns == 'no': column_method = None
                
                if self.DirFile not in self.heatmap_run:
                    self.heatmap_run[self.DirFile]=None
                    ### occurs when automatically running the heatmap
                    column_method = None
                    row_method = None
            
                translate={'None Selected':'','Exclude Cell Cycle Effects':'excludeCellCycle','Top Correlated Only':'top','Positive Correlations Only':'positive','Perform Iterative Discovery':'driver', 'Intra-Correlated Only':'IntraCorrelatedOnly', 'Perform Monocle':'monocle'}
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
                        if 'top' in GeneSelection or 'driver' in GeneSelection or 'excludeCellCycle' in GeneSelection or 'positive' in GeneSelection or 'IntraCorrelatedOnly' in GeneSelection:
                            GeneSelection+=' amplify'
                    except Exception: pass
                GeneSetSelection = string.replace(GeneSetSelection,'\n',' ')
                GeneSetSelection = string.replace(GeneSetSelection,'\r',' ')
                        
                try: CorrelationCutoff = float(gu.Results()['CorrelationCutoff'])
                except Exception: CorrelationCutoff=None
                if transpose == 'yes': transpose = True
                else: transpose = False
        
                vendor = 'RNASeq'
                color_gradient = string.replace(color_gradient,'-','_')
                if len(GeneSetSelection)>0 or GeneSelection != '':
                    gsp = UI.GeneSelectionParameters(self.species,self.platform,vendor)
                    gsp.setGeneSet(GeneSetSelection)
                    gsp.setPathwaySelect(PathwaySelection)
                    gsp.setGeneSelection(GeneSelection)
                    gsp.setOntologyID(OntologyID)
                    gsp.setTranspose(transpose)
                    gsp.setNormalize(normalization)
                    gsp.setJustShowTheseIDs(justShowTheseIDs)
                    try: gsp.setClusterGOElite(clusterGOElite)
                    except Exception: pass
                    if rho!=None:
                        try:
                            float(rho)
                            gsp.setRhoCutoff(rho)
                        except Exception: print 'Must enter a valid Pearson correlation cutoff (float)'
                    transpose = gsp ### this allows methods that don't transmit this object to also work
                if row_method == 'no': row_method = None
                if column_method == 'no': column_method = None
    
                UI.createHeatMap(input_file_dir, row_method, row_metric, column_method, column_metric, color_gradient, transpose, contrast, None, display=True)
                    
            else:
                os.chdir(parentDirectory)
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
                if(self.Yes1Radio.GetValue() == True):
                    include_labels= 'yes'
                else:
                    include_labels= 'no'    
                pca_algorithm = 'SVD'
                transpose = False
                if(self.D_3DRadio.GetValue() == True):
                    plotType = '3D'
                else:
                    plotType = '2D'
                display = True
                UI.performPCA(input_file_dir, include_labels, pca_algorithm, transpose, None, plotType=plotType, display=display)
                os.chdir(currentDirectory)
            self.InteractivePanelUpdate(event)                                                                                                                                                                     
        

    def OnAbout(self, event):
        dial = wx.MessageDialog(None, 'AltAnalyze Results Viewer\nVersion 0.5\n2015', 'About', wx.OK)
        dial.ShowModal()

    def OnHelp(self, event):
        os.chdir(parentDirectory)
        ManualPath = str(os.getcwd()) + "/Documentation/Source_files/ViewerManual.doc"
        subprocess.Popen(['open', ManualPath])
        os.chdir(currentDirectory)        
                        
class ImageFrame(wx.Frame):
    title = "Image"

    def __init__(self):
        wx.Frame.__init__(self, None, title=self.title)        

def remoteViewer(app):
    currentDirectory = str(os.getcwd())
    os.chdir(currentDirectory+'/Config') 
    #app = wx.App(False)
    fr = Main(parent=None,id=1)
    fr.Show()
    app.MainLoop()
    os.chdir(currentDirectory)
    
if __name__ == "__main__":
    app = wx.App(False)
    fr = Main(parent=None,id=1)
    fr.Show()
    app.MainLoop()
