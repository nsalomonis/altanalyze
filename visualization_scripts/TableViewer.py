"""
Code obtained from: https://github.com/eiinfuva/dynamic-typesetting-ac/blob/master/multilistbox.py

"""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

from Tkinter import *
import webbrowser

MOVE_LINES = 0
MOVE_PAGES = 1
MOVE_TOEND = 2

class TableViewer(Frame):
    def __init__(self, master, lists, command=None, **options):
        defaults = {
            'background': None,
            'borderwidth': 2,
            'font': None,
            'foreground': None,
            'height': 10,
            'highlightcolor': None,
            'highlightthickness': 1,
            'relief': SUNKEN,
            'takefocus': 1,
            }

        aliases = {'bg':'background', 'fg':'foreground', 'bd':'borderwidth'}

        for k in aliases.keys ():
            if options.has_key (k):
                options [aliases[k]] = options [k]
            
        for key in defaults.keys():
            if not options.has_key (key):
                options [key] = defaults [key]

        apply (Frame.__init__, (self, master), options)
        self.lists = []

        # MH (05/20049
        # These are needed for sorting
        self.colmapping={}
        self.origData = None

        #  Keyboard navigation.
        
        self.bind ('<Up>',    lambda e, s=self: s._move (-1, MOVE_LINES))
        self.bind ('<Down>',  lambda e, s=self: s._move (+1, MOVE_LINES))
        self.bind ('<Prior>', lambda e, s=self: s._move (-1, MOVE_PAGES))
        self.bind ('<Next>',  lambda e, s=self: s._move (+1, MOVE_PAGES))
        self.bind ('<Home>',  lambda e, s=self: s._move (-1, MOVE_TOEND))
        self.bind ('<End>',   lambda e, s=self: s._move (+1, MOVE_TOEND))
        if command:
            self.bind ('<Return>', command)

        # Columns are a frame with listbox and label in it.
        
        # MH (05/2004):
        # Introduced a PanedWindow to make the columns resizable
        
        m = PanedWindow(self, orient=HORIZONTAL, bd=0, 
            background=options['background'], showhandle=0, sashpad=1)
        m.pack(side=LEFT, fill=BOTH, expand=1)
        for label, width in lists:
            lbframe = Frame(m)
            m.add(lbframe, width=width)
            # MH (05/2004)
            # modified this, click to sort
            b = Label(lbframe, text=label, borderwidth=1, relief=RAISED)
            b.pack(fill=X)
            b.bind('<Button-1>', self._sort)

            self.colmapping[b]=(len(self.lists),1)
            
            lb = Listbox (lbframe,
                          width=width,
                          height=options ['height'],
                          borderwidth=0,
                          font=options ['font'],
                          background=options ['background'],
                          selectborderwidth=0,
                          relief=SUNKEN,
                          takefocus=FALSE,
                          exportselection=FALSE)
            lb.pack (expand=YES, fill=BOTH)
            self.lists.append (lb)

            # Mouse features
            
            lb.bind ('<B1-Motion>', lambda e, s=self: s._select (e.y))
            lb.bind ('<Button-1>',  lambda e, s=self: s._select (e.y))
            lb.bind ('<Leave>',     lambda e: 'break')
            lb.bind ('<B2-Motion>', lambda e, s=self: s._b2motion (e.x, e.y))
            lb.bind ('<Button-2>',  lambda e, s=self: s._button2 (e.x, e.y))
            if command:
                lb.bind ('<Double-Button-1>', command)

        sbframe = Frame (self)
        sbframe.pack (side=LEFT, fill=Y)
        l = Label (sbframe, borderwidth=1, relief=RAISED)
        l.bind ('<Button-1>', lambda e, s=self: s.focus_set ())
        l.pack(fill=X)
        sb = Scrollbar (sbframe,
                        takefocus=FALSE,
                        orient=VERTICAL,
                        command=self._scroll)
        sb.pack (expand=YES, fill=Y)
        self.lists[0]['yscrollcommand']=sb.set

        return


    # MH (05/2004)
    # Sort function, adopted from Rick Lawson 
    # http://tkinter.unpythonic.net/wiki/SortableTable
    
    def _sort(self, e):
        # get the listbox to sort by (mapped by the header button)
        b=e.widget
        col, direction = self.colmapping[b]

        # get the entire table data into mem
        tableData = self.get(0,END)
        if self.origData == None:
            import copy
            self.origData = copy.deepcopy(tableData)

        rowcount = len(tableData)

        #remove old sort indicators if it exists
        for btn in self.colmapping:
            lab = btn.cget('text')
            if lab[0]=='[': btn.config(text=lab[4:])

        btnLabel = b.cget('text')
        #sort data based on direction
        if direction==0:
            tableData = self.origData
        else:
            if direction==1: b.config(text='[+] ' + btnLabel)
            else: b.config(text='[-] ' + btnLabel)
            # sort by col
            tableData.sort(key=lambda x: x[col], reverse=direction<0)

        #clear widget
        self.delete(0,END)

        # refill widget
        for row in range(rowcount):
            self.insert(END, tableData[row])

        # toggle direction flag 
        if direction==1: direction=-1
        else: direction += 1
        self.colmapping[b] = (col, direction) 
        

    def _move (self, lines, relative=0):
        """
        Move the selection a specified number of lines or pages up or
        down the list.  Used by keyboard navigation.
        """
        selected = self.lists [0].curselection ()
        try:
            selected = map (int, selected)
        except ValueError:
            pass

        try:
            sel = selected [0]
        except IndexError:
            sel = 0

        old  = sel
        size = self.lists [0].size ()
        
        if relative == MOVE_LINES:
            sel = sel + lines
        elif relative == MOVE_PAGES:
            sel = sel + (lines * int (self.lists [0]['height']))
        elif relative == MOVE_TOEND:
            if lines < 0:
                sel = 0
            elif lines > 0:
                sel = size - 1
        else:
            print "TableViewer._move: Unknown move type!"

        if sel < 0:
            sel = 0
        elif sel >= size:
            sel = size - 1
        
        self.selection_clear (old, old)
        self.see (sel)
        self.selection_set (sel)
        return 'break'


    def _select (self, y):
        """
        User clicked an item to select it.
        """
        row = self.lists[0].nearest (y)
        label = labels[row]
        try: webbrowser.open('http://www.genecards.org/cgi-bin/carddisp.pl?gene='+label)
        except Exception: null=[]
        self.selection_clear (0, END)
        self.selection_set (row)
        self.focus_set ()
        return 'break'


    def _button2 (self, x, y):
        """
        User selected with button 2 to start a drag.
        """
        for l in self.lists:
            l.scan_mark (x, y)
        return 'break'


    def _b2motion (self, x, y):
        """
        User is dragging with button 2.
        """
        for l in self.lists:
            l.scan_dragto (x, y)
        return 'break'


    def _scroll (self, *args):
        """
        Scrolling with the scrollbar.
        """
        for l in self.lists:
            apply(l.yview, args)

    def curselection (self):
        """
        Return index of current selection.
        """
        return self.lists[0].curselection()


    def delete (self, first, last=None):
        """
        Delete one or more items from the list.
        """
        for l in self.lists:
            l.delete(first, last)


    def get (self, first, last=None):
        """
        Get items between two indexes, or one item if second index
        is not specified.
        """
        result = []
        for l in self.lists:
            result.append (l.get (first,last))
        if last:
            return apply (map, [None] + result)
        return result

            
    def index (self, index):
        """
        Adjust the view so that the given index is at the top.
        """
        for l in self.lists:
            l.index (index)


    def insert (self, index, *elements):
        """
        Insert list or tuple of items.
        """
        for e in elements:
            i = 0
            for l in self.lists:
                l.insert (index, e[i])
                i = i + 1
        if self.size () == 1:
            self.selection_set (0)
            

    def size (self):
        """
        Return the total number of items.
        """
        return self.lists[0].size ()


    def see (self, index):
        """
        Make sure given index is visible.
        """
        for l in self.lists:
            l.see (index)


    def selection_anchor (self, index):
        """
        Set selection anchor to index.
        """
        for l in self.lists:
            l.selection_anchor (index)


    def selection_clear (self, first, last=None):
        """
        Clear selections between two indexes.
        """
        for l in self.lists:
            l.selection_clear (first, last)


    def selection_includes (self, index):
        """
        Determine if given index is selected.
        """
        return self.lists[0].selection_includes (index)


    def selection_set (self, first, last=None):
        """
        Select a range of indexes.
        """
        for l in self.lists:
            l.selection_set (first, last)

def viewTable(title,header,list_values):
    global labels
    tk = Toplevel()
    Label(tk, text=title).pack()
    column_values = []
    labels=[]
    import string
    for i in header:
        width = 120
        column_values.append((i,width))
        
    mlb = TableViewer (tk,
                        tuple(column_values),
                        height=20,
                        bg='white')
    for i in list_values:
        mlb.insert (END, i)
        labels.append(i[0])
    mlb.pack (expand=YES,fill=BOTH)
    def deleteWindow():
        tk.quit()
        tk.destroy() ### just quit instead
    try: tk.clipboard_clear()
    except Exception: pass
    try: tk.clipboard_append(string.join(labels,'\n'))
    except Exception: pass

    Button (tk, text="Close", command = deleteWindow).pack ()
    tk.protocol("WM_DELETE_WINDOW", deleteWindow)
    tk.mainloop()

if __name__ == '__main__':
    title = 'test'
    list_values = [(1,2,3)]
    header = ['1','2','3']
    header = ['Associated Genes']
    list_values = [[('CAMK2D')],[('TTN')]]
    viewTable(title,header,list_values)