'''Example 09

Shows how to embed visvis figures in a wxPython GUI.
Visvis can also be embedded in Qt4 (PyQt), GTK, and FLTK.
Based on http://code.google.com/p/visvis/wiki/example_embeddingInWx
'''

import wx
from pycalfem import *
from pycalfem_utils import *
from pycalfem_mesh import *
import pycalfem_vis as pcv
import visvis as vv

# Create a visvis app instance, which wraps a wx application object.
# This needs to be done *before* instantiating the main window. 

app = vv.use('qt')


class MainWindow(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Embedding in WX", size=(560, 420))
        
        # Make a panel with buttons
        self.panel = wx.Panel(self)
        but1 = wx.Button(self.panel, -1, 'Calc')
        but2 = wx.Button(self.panel, -1, 'Plot')
        but3 = wx.Button(self.panel, -1, 'Clear')
        
        #Make panel sizer and embed stuff
        self.panelsizer = wx.BoxSizer(wx.VERTICAL)
        self.panelsizer.Add(but1)
        self.panelsizer.Add(but2)
        self.panelsizer.Add(but3)
        self.panel.SetSizer(self.panelsizer)
        
        # Make figure using "self" as a parent
        Figure = app.GetFigureClass()
        self.fig = Figure(self)
        
        # Make window sizer and embed stuff
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add(self.panel, 1, wx.EXPAND)
        self.sizer.Add(self.fig._widget, 2, wx.EXPAND)
        
        # Make callback
        but1.Bind(wx.EVT_BUTTON, self._Calc)
        but2.Bind(wx.EVT_BUTTON, self._Plot)
        but3.Bind(wx.EVT_BUTTON, self._Clear)
        
        # Apply window sizers        
        self.SetSizer(self.sizer)
        self.SetAutoLayout(True)
        self.Layout()   
        
        # Finish
        self.Show()
    
  
    def _Calc(self, event):
        #Calculations are taken from example 08.
        # Constants:
        t = 0.2
        v = 0.35
        E = 2.1e9
        ptype = 1
        ep = [ptype,t]
        D=hooke(ptype, E, v)
        
        # Create Geometry:
        g = "examplegeo\ex8.geo"
        self.elType = 3 #3 Quads
        self.dofsPerNode = 2
        mesher = GmshMesher(geoData = g,
                            gmshExecPath = None, #"gmsh\gmsh.exe"
                            elType = self.elType, 
                            dofsPerNode= self.dofsPerNode)
        self.coords, self.edof, self.dofs, self.bdofs, _ = mesher.create()
        
        # Assem systems matrix:
        nDofs = size(self.dofs)
        ex, ey = coordxtr(self.edof, self.coords, self.dofs)
        K = zeros([nDofs,nDofs])
        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = planqe(elx, ely, ep, D)
            assem(eltopo, K, Ke)
        
        # Solve:
        f = zeros([nDofs,1])
        bc = array([],'i')
        bcVal = array([],'i')
        bc, bcVal = applybc(self.bdofs, bc, bcVal, 5, 0.0, 0)
        applyforce(self.bdofs, f, 7, 10e5, 1)
        self.a, _ = solveq(K,f,bc,bcVal)
    
      
    def _Plot(self, event):
        # Make sure our figure is the active one
        # If only one figure, this is not necessary.
        #vv.figure(self.fig.nr)
        
        # Clear it:
        vv.clf()
        
        # Plot:
        pcv.drawDisplacements(self.a, self.coords, self.edof, self.dofsPerNode, self.elType, doDrawUndisplacedMesh=True, title="Example 09")
    
    def _Clear(self, event):
        vv.clf() #Clear current figure

# Two ways to create the application and start the main loop
if True:
    # The visvis way. Will run in interactive mode when used in IEP or IPython.
    app.Create()
    m = MainWindow()
    app.Run()

else:
    # The native way.
    wxApp = wx.App()    
    m = MainWindow()
    wxApp.MainLoop()