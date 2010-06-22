"""DynamicLattice: supports for animation of two-dimensional arrays
and rectangular lattices, by mapping values of an input array to a
grayscale colormap.

---

Typical usage:

>>> dynlat = DynamicLattice((Nx, Ny), zmin, zmax)

initializes a DynamicLattice instance (named, e.g., dynlat)
to animate an array with shape (Nx, Ny), with arrays values
mapped uniformly on a grayscale colormap:

grayscale = (array[i,j]-zmin)/(zmax-zmin),
with clipping at grayscale=0 for array[i,j] <= zmin
and grayscale=1 for array[i,j] >= zmax.

>>> dynlat.display(a)

updates the dynlat display with the current values of the array a

---

DynamicLattice also supports limited mouse interaction to delineate
selected boxes in the display window.

>>> dynlat.IsBoxSelected()

returns True if a rectangular box has been selected with the mouse

>>> dynlat.GetMouseBox()

returns (x0,y0,x1,y1) describing a box with lower and upper coordinates
(x0, y0) and (x1, y1), respectively, if a box has been selected

---

>>> dynlat.setTitle(s)

sets the DynamicLattice window title to the string s

"""

import math
import Tkinter 
import Image, ImageTk

DefaultImageSize = (400,400)

root = Tkinter.Tk()
root.withdraw()

class DynamicLattice (Tkinter.Label):

    def __init__(self, shape, size=DefaultImageSize, mode='P', 
                 zmin=0.0, zmax=1.0):
        top = Tkinter.Toplevel()
        top.title('DynamicLattice')
        Tkinter.Label.__init__(self, top)
	self.shape = shape
	self.size = size
	self.mode = mode
	self.zmin = zmin
	self.zmax = zmax
	self.zrange = zmax-zmin
        self.canvas = Tkinter.Canvas(top, 
				     width=self.size[0], height=self.size[1])
        Tkinter.Widget.bind(self.canvas, "<Button-1>",
                           self.mouseDown)
        Tkinter.Widget.bind(self.canvas, "<Button1-ButtonRelease>",
                           self.mouseUp)
        self.mousecoords = []
	self.im = Image.new(self.mode, self.shape)
	self.displayIm = self.im.resize(self.size)
	if self.mode == '1':
	    self.tkimage = \
			 ImageTk.BitmapImage(self.displayIm,
					     foreground="white")
	else:
	    self.tkimage = \
			 ImageTk.PhotoImage(self.displayIm)
	self.canvas.create_image(0, 0, anchor=Tkinter.NW, image=self.tkimage)
        self.canvas.pack()

    def setTitle(self, title):
        self.master.title(title)

    def mouseDown(self, event):
        x0 = self.canvas.canvasx(event.x)
        y0 = self.canvas.canvasx(event.y)
        sx0 = int(self.shape[0] * float(x0)/self.size[0])
        sy0 = int(self.shape[1] * float(y0)/self.size[1])
        self.mousecoords = [sx0, sy0]

    def mouseUp(self, event):
        sx0, sy0 = self.mousecoords
        x1 = self.canvas.canvasx(event.x)
        y1 = self.canvas.canvasx(event.y)
        sx1 = int(self.shape[0] * float(x1)/self.size[0])
        sy1 = int(self.shape[1] * float(y1)/self.size[1])
        X0, X1 = min(sx0, sx1), max(sx0, sx1)
        Y0, Y1 = min(sy0, sy1), max(sy0, sy1)
        self.mousecoords = [X0, Y0, X1, Y1]
    
    def IsBoxSelected(self):
        return len(self.mousecoords)==4

    def GetMouseBox(self):
        mc = self.mousecoords[:]
        self.mousecoords = []
        return mc

    def display(self, array, site=None):
        if site is not None:
            self.im.putpixel(site, self.grayscale(array[site]))
	else:
	    for i in range(array.shape[0]):
                for j in range(array.shape[1]):
                    self.im.putpixel((i,j), self.grayscale(array[i,j]))
	self.displayIm = self.im.resize(self.size)
        self.tkimage.paste(self.displayIm)
        self.canvas.update()

    def grayscale(self, value):
        sval = (value-self.zmin)/(self.zrange)
        if sval < 0.: sval = 0.
        if sval > 1.: sval = 1.
        return int(255.*(sval))

    def set_zmin(self, zmin):
        self.zmin = zmin
        self.zrange = self.zmax - self.zmin
        
    def set_zmax(self, zmax):
        self.zmax = zmax
        self.zrange = self.zmax - self.zmin
        
	    
