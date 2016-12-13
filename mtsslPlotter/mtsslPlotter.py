import os
import sys
import threading
import warnings
warnings.filterwarnings('ignore')
import string
import math
import numpy
from PIL import Image
import pymol
from pymol import stored, cmd
import matplotlib
matplotlib.use('WXAgg', warn=False, force=False)
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
import wx
import cPickle as pickle
try:
	import colormaps as cmaps
except:
	pass
from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub
from ObjectListView import ObjectListView, ColumnDefn, OLVEvent

##########################################################################################
##  Main Window
##########################################################################################
class MainWindow(wx.Frame):
	def __init__(self):
		wx.Frame.__init__(self, None, -1, "mtsslPlotter",  wx.DefaultPosition, wx.Size(800,700))
		self.Bind(wx.EVT_CLOSE, self.closeWindow)
		self.splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
		self.leftP = LeftPanel(self.splitter)
		self.rightP = Plotter(self.splitter)
		self.splitter.SplitVertically(self.leftP, self.rightP)
		self.splitter.SetMinimumPaneSize(1)
		self.splitter.SetSashGravity(0)
		self.splitter.SetSashPosition(100, redraw=True)
		self.Bind(wx.EVT_SPLITTER_DCLICK, self.toggleSplitter, None)
		self.sizer = wx.BoxSizer(wx.VERTICAL)
		self.sizer.Add(self.splitter, 1, wx.EXPAND)
		self.SetSizer(self.sizer)
		self.statusBar = wx.StatusBar(self, -1)
		self.statusBar.SetFieldsCount(1)
		self.SetStatusBar(self.statusBar)
		self.initMenu()
		pub.subscribe(self.OnUpdateStatusBar, 'update.statusBar')

	def toggleSplitter(self, Sender):
		if self.splitter.GetSashPosition() > 1:
			self.splitter.SetSashPosition(1)

	def initMenu(self):
		menubar = wx.MenuBar()
		# File menu
		fileMenu = wx.Menu()
		openItem = wx.MenuItem(fileMenu, 11, '&Load Session...\tCtrl+O')
		self.Bind(wx.EVT_MENU, self.OnOpen, id=11)
		saveItem = wx.MenuItem(fileMenu, 12, '&Save Session as...\tCtrl+S')
		self.Bind(wx.EVT_MENU, self.OnSave, id=12)
		savePlotDataItem = wx.MenuItem(fileMenu, 41, '&Save Selected Plot Data as...\tCtrl+Shift+S')
		self.Bind(wx.EVT_MENU, self.OnSavePlotData, id=41)
		refreshItem = wx.MenuItem(fileMenu, 13, '&Refresh\tCtrl+R')
		self.Bind(wx.EVT_MENU, self.OnRefresh, id=13)
		exitItem = wx.MenuItem(fileMenu, 14, '&Exit\tAlt+F4')
		self.Bind(wx.EVT_MENU, self.OnClose, id=14)
		deerAnalysisItem = wx.MenuItem(fileMenu, 15, '&Load DeerAnalysis Result...\tShift+Ctrl+O')
		self.Bind(wx.EVT_MENU, self.OnDeerAnalysis, id=15)
		fileMenu.AppendItem(openItem)
		fileMenu.AppendItem(saveItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(savePlotDataItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(refreshItem)
		fileMenu.AppendItem(deerAnalysisItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(exitItem)
		
		#Edit menu
		editMenu = wx.Menu()
		deletePlotItem = wx.MenuItem(editMenu, 42, '&Delete Selected Plots\tBACK')
		self.Bind(wx.EVT_MENU, self.OnDelete, id=42)
		editMenu.AppendItem(deletePlotItem)
		
		# Help menu
		helpMenu = wx.Menu()
		aboutItem = wx.MenuItem(helpMenu, 31, '&About\tCtrl+A')
		self.Bind(wx.EVT_MENU, self.OnAbout, id=31)
		helpMenu.AppendItem(aboutItem)
		
		menubar.Append(fileMenu, '&File')
		menubar.Append(editMenu, '&Edit')
		menubar.Append(helpMenu, '&Help')
		self.SetMenuBar(menubar)

	def OnSavePlotData(self, event):
		pub.sendMessage('file.savePlot')

	def OnDelete(self, event):
		pub.sendMessage('edit.delete')

	def OnDeerAnalysis(self, event):
		dlg = wx.FileDialog(self, message="Choose a DeerAnalysis result...", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style=wx.OPEN | wx.CHANGE_DIR)
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
		dlg.Destroy()
		try:
			result = numpy.loadtxt(path)
		except:
			pub.sendMessage('update.statusBar', data = "Cannot open file!")
			return
		graphData = {}
		graphData['Title'] = os.path.basename(path)
		graphData['Type'] = "DistanceDistribution"
		graphData['xTitle'] = "Distance (Angstrom)"
		graphData['yTitle'] = "P(r)"
		graphData['xData'] = result[:,0]*10
		graphData['yData'] = (result[:,1]-result[:,1].min())/(result[:,1].max()-result[:,1].min())#result[:,1]
		graphData['zData'] = 0
		graphData["xlim"] = [0, 100]
		graphData["ylim"] = [0, 1]
		stored.plots.append(graphData)
		pub.sendMessage('list.updated')
		pub.sendMessage('update.statusBar', data = "DeerAnalysis Result loaded!")

	def OnUpdateStatusBar(self, data):
		self.statusBar.SetStatusText(data)

	def OnRefresh(self, event):
		pub.sendMessage('list.updated')
		pub.sendMessage('update.statusBar', data = "Synchronized with PyMOL.")

	def OnClear(self, event):
		self.graphs = []
		self.graphsList.Clear()
		event.Skip()
	
	def closeWindow(self, event):
		dlg = wx.MessageDialog(self, 
			"Do you really want to close this application?",
			"Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
		result = dlg.ShowModal()
		dlg.Destroy()
		if result == wx.ID_OK:
			stored.mtsslplot = 0
			print stored.mtsslplot
			self.Destroy()
			

	def OnOpen(self, event):
		dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style=wx.OPEN | wx.CHANGE_DIR)
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
			print "Unpickling plots...",
			stored.plots = pickle.load(open(path, 'rb'))
			print "Done!\n"
		dlg.Destroy()
		pub.sendMessage('list.updated')
		pub.sendMessage('update.statusBar', data = "Plots loaded!")

	def OnSave(self, event):
		dialog = wx.FileDialog(None, message="Save data as ...", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style = wx.SAVE | wx.OVERWRITE_PROMPT)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			dataFile = open(path, 'w')
			print "\nPickling plots...",
			pickle.dump(stored.plots, open(path, 'wb')) 
			print "Done!\n"
			dataFile.close()
		dialog.Destroy()
		pub.sendMessage('update.statusBar', data = "Plots saved!")

	def OnClose(self, event):
		stored.mtsslplot = 0
		self.Destroy()
	
	def OnAbout(self, event):
		message = "mtsslPlotter 2.0\nG. Hagelueken & Dinar Abdullin\n2015"
		caption = "About"
		dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.ICON_INFORMATION)
		dlg.ShowModal()
		dlg.Destroy()
		

##########################################################################################
##  PlotListEntry (Objects that are listed in the Left Panel)
##########################################################################################

class PlotListEntry(object):
	def __init__(self, olv, type, title, color, idx):
		self.type = type
		self.title = title
		self.index = idx
		self.icon = self.createIcon(olv, color)
		self.isChecked = False
	
	def toggleState(self):
		if self.isChecked == False:
			self.isChecked = True
		else:
			self.isChecked = False
	
	def setCheckedState(self, state):
		if state == False:
			self.isChecked = False
		elif state == True:
			self.isChecked = True
		
	
	def recolorIcon(self, filename, color):
		im = Image.open(filename)
		data = numpy.array(im)
		#convert color from hex to rgb
		value = color.lstrip('#')
		lv = len(value)
		rgb = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
		r1, g1, b1 = 0, 0, 0 # Original value
		r2, g2, b2 = rgb[0], rgb[1],rgb[2] # Value that we want to replace it with
		red, green, blue = data[:,:,0], data[:,:,1], data[:,:,2]
		mask = (red == r1) & (green == g1) & (blue == b1)
		data[:,:,:3][mask] = [r2, g2, b2]
		im = Image.fromarray(data)
		im.save(os.path.dirname(__file__)+"/icons/tmp.png")
	
	def createIcon(self, olv, color=None):
		#print "in createIcon"
		if self.type == "DistanceDistribution":
			filename = os.path.dirname(__file__)+"/icons/dist.png"
			self.recolorIcon(filename, color)	
			filename = os.path.dirname(__file__)+"/icons/tmp.png"
		elif self.type == "DistanceMap":
			filename = os.path.dirname(__file__)+"/icons/dm.png"
		elif self.type == "DistancePlot":
			filename = os.path.dirname(__file__)+"/icons/diag.png"
			self.recolorIcon(filename, color)	
			filename = os.path.dirname(__file__)+"/icons/tmp.png"
		elif self.type == "AccessibilityPlot":
			filename = os.path.dirname(__file__)+"/icons/acc.png"
			self.recolorIcon(filename, color)	
			filename = os.path.dirname(__file__)+"/icons/tmp.png"
		elif self.type == "DistanceDistributions":
			filename = os.path.dirname(__file__)+"/icons/dist.png"
		elif self.type == "DifferenceMap":
			filename = os.path.dirname(__file__)+"/icons/ddm.png"
		img = wx.Image(filename, wx.BITMAP_TYPE_ANY)
		img = img.Scale(16, 16, wx.IMAGE_QUALITY_HIGH)
		img = img.ConvertToBitmap()
		indexInOlv = olv.AddImages(img)
		return indexInOlv

##########################################################################################
##  LeftPanel (List View)
##########################################################################################

class LeftPanel(wx.Panel):	
	def __init__(self, parent):
		wx.Panel.__init__(self, parent=parent)
		self.plots = []
		self.pickedColors = []
		self.colors = self.setupColors()
		self.plotsOlv = ObjectListView(self, wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER, useAlternateBackColors=False)
		self.plotsOlv.Bind(OLVEvent.EVT_ITEM_CHECKED, self.checkboxClicked)
		self.plotsOlv.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.onDoubleClick)
		#self.plotsOlv.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.onRightClick)
		self.getPlots()
		self.setColumns()
		self.plotsOlv.CreateCheckStateColumn()
		self.plotsOlv.SetObjects(self.plots)
		bsizer = wx.BoxSizer()
		#bsizer.Add(lb, 1, wx.EXPAND)
		bsizer.Add(self.plotsOlv, 1, wx.EXPAND)
		self.SetSizer(bsizer)
		self.checkedPlots = {}
		pub.subscribe(self.refreshList, 'list.updated')
		pub.subscribe(self.checkboxClicked, 'checkbox.clicked')
		pub.subscribe(self.onSavePlot, 'file.savePlot')
		pub.subscribe(self.onDelete, 'edit.delete')
	
	def onDoubleClick(self, event):
		#print event.m_itemIndex
		dlg = wx.ColourDialog(self)
		# Ensure the full colour dialog is displayed, 
		# not the abbreviated version.
		dlg.GetColourData().SetChooseFull(True)
		if dlg.ShowModal() == wx.ID_OK:
			data = dlg.GetColourData()
		color = data.GetColour().Get()
		htmlColor = "#%02x%02x%02x" % color
		dlg.Destroy()
		index = event.m_itemIndex
		stored.plots[index]['Color'] = htmlColor
		self.plots[index].icon = self.plots[index].createIcon(self.plotsOlv, htmlColor)
		self.plotsOlv.RefreshObject(self.plots[index])
		pub.sendMessage('item.checked', data = self.checkedPlots)

	def onDelete(self):
		if self.plotsOlv.GetSelectedObjects():
			objectList = self.plotsOlv.GetSelectedObjects()
			for object in objectList:
				index = object.index
				stored.plots.pop(index)
				self.plots.pop(index)
			self.refreshList()
			pub.sendMessage('update.statusBar', data = "Plots deleted!")
		else:
			pub.sendMessage('update.statusBar', data = "Cannot delete plots! Please select one ore more plots!")
		
	def onSavePlot(self):
		if self.plotsOlv.GetSelectedObject():
			index = self.plotsOlv.GetSelectedObject().index
			data = stored.plots[index]
			outputString =  "%s\n" %data['Title']
			outputString += "%s\n" %data['Type']
			outputString += "%s,%s\n" %(data['xTitle'], data['xTitle'])
			for idx, row in enumerate(data['xData'].tolist()):
				z = 0
				outputString += "%1.3f,%1.3f\n" %(data['xData'][idx], data['yData'][idx])
			if stored.plots[index]['Type'] == "DifferenceMap" or stored.plots[index]['Type'] == "DistanceMap":
				outputString += "z-Data:\n"
				numpy.set_printoptions(threshold=numpy.nan)
				zData = numpy.array_str(data['zData'])
				outputString += zData
			dialog = wx.FileDialog(None, message="Save plot data as ...", defaultDir=os.getcwd(), defaultFile=stored.plots[index]['Title'], wildcard="*.dat", style = wx.SAVE | wx.OVERWRITE_PROMPT)
			if dialog.ShowModal() == wx.ID_OK:
				path = dialog.GetPath()
				dataFile = open(path, 'w')
				print "\nSaving data...",
				dataFile.write(outputString)
				print "Done!\n"
				dataFile.close()
			dialog.Destroy()
			pub.sendMessage('update.statusBar', data = "Plot data saved!")
		else:
			pub.sendMessage('update.statusBar', data = "Cannot save data! Please select one plot!")
		
		
			
	def setupColors(self):
		cnames = {
					'aliceblue':            '#F0F8FF',
					'antiquewhite':         '#FAEBD7',
					'aqua':                 '#00FFFF',
					'aquamarine':           '#7FFFD4',
					'azure':                '#F0FFFF',
					'beige':                '#F5F5DC',
					'bisque':               '#FFE4C4',
					'black':                '#000000',
					'blanchedalmond':       '#FFEBCD',
					'blue':                 '#0000FF',
					'blueviolet':           '#8A2BE2',
					'brown':                '#A52A2A',
					'burlywood':            '#DEB887',
					'cadetblue':            '#5F9EA0',
					'chartreuse':           '#7FFF00',
					'chocolate':            '#D2691E',
					'coral':                '#FF7F50',
					'cornflowerblue':       '#6495ED',
					'cornsilk':             '#FFF8DC',
					'crimson':              '#DC143C',
					'cyan':                 '#00FFFF',
					'darkblue':             '#00008B',
					'darkcyan':             '#008B8B',
					'darkgoldenrod':        '#B8860B',
					'darkgray':             '#A9A9A9',
					'darkgreen':            '#006400',
					'darkkhaki':            '#BDB76B',
					'darkmagenta':          '#8B008B',
					'darkolivegreen':       '#556B2F',
					'darkorange':           '#FF8C00',
					'darkorchid':           '#9932CC',
					'darkred':              '#8B0000',
					'darksalmon':           '#E9967A',
					'darkseagreen':         '#8FBC8F',
					'darkslateblue':        '#483D8B',
					'darkslategray':        '#2F4F4F',
					'darkturquoise':        '#00CED1',
					'darkviolet':           '#9400D3',
					'deeppink':             '#FF1493',
					'deepskyblue':          '#00BFFF',
					'dimgray':              '#696969',
					'dodgerblue':           '#1E90FF',
					'firebrick':            '#B22222',
					'floralwhite':          '#FFFAF0',
					'forestgreen':          '#228B22',
					'fuchsia':              '#FF00FF',
					'gainsboro':            '#DCDCDC',
					'ghostwhite':           '#F8F8FF',
					'gold':                 '#FFD700',
					'goldenrod':            '#DAA520',
					'gray':                 '#808080',
					'green':                '#008000',
					'greenyellow':          '#ADFF2F',
					'honeydew':             '#F0FFF0',
					'hotpink':              '#FF69B4',
					'indianred':            '#CD5C5C',
					'indigo':               '#4B0082',
					'ivory':                '#FFFFF0',
					'khaki':                '#F0E68C',
					'lavender':             '#E6E6FA',
					'lavenderblush':        '#FFF0F5',
					'lawngreen':            '#7CFC00',
					'lemonchiffon':         '#FFFACD',
					'lightblue':            '#ADD8E6',
					'lightcoral':           '#F08080',
					'lightcyan':            '#E0FFFF',
					'lightgoldenrodyellow': '#FAFAD2',
					'lightgreen':           '#90EE90',
					'lightgray':            '#D3D3D3',
					'lightpink':            '#FFB6C1',
					'lightsalmon':          '#FFA07A',
					'lightseagreen':        '#20B2AA',
					'lightskyblue':         '#87CEFA',
					'lightslategray':       '#778899',
					'lightsteelblue':       '#B0C4DE',
					'lightyellow':          '#FFFFE0',
					'lime':                 '#00FF00',
					'limegreen':            '#32CD32',
					'linen':                '#FAF0E6',
					'magenta':              '#FF00FF',
					'maroon':               '#800000',
					'mediumaquamarine':     '#66CDAA',
					'mediumblue':           '#0000CD',
					'mediumorchid':         '#BA55D3',
					'mediumpurple':         '#9370DB',
					'mediumseagreen':       '#3CB371',
					'mediumslateblue':      '#7B68EE',
					'mediumspringgreen':    '#00FA9A',
					'mediumturquoise':      '#48D1CC',
					'mediumvioletred':      '#C71585',
					'midnightblue':         '#191970',
					'mintcream':            '#F5FFFA',
					'mistyrose':            '#FFE4E1',
					'moccasin':             '#FFE4B5',
					'navajowhite':          '#FFDEAD',
					'navy':                 '#000080',
					'oldlace':              '#FDF5E6',
					'olive':                '#808000',
					'olivedrab':            '#6B8E23',
					'orange':               '#FFA500',
					'orangered':            '#FF4500',
					'orchid':               '#DA70D6',
					'palegoldenrod':        '#EEE8AA',
					'palegreen':            '#98FB98',
					'paleturquoise':        '#AFEEEE',
					'palevioletred':        '#DB7093',
					'papayawhip':           '#FFEFD5',
					'peachpuff':            '#FFDAB9',
					'peru':                 '#CD853F',
					'pink':                 '#FFC0CB',
					'plum':                 '#DDA0DD',
					'powderblue':           '#B0E0E6',
					'purple':               '#800080',
					'red':                  '#FF0000',
					'rosybrown':            '#BC8F8F',
					'royalblue':            '#4169E1',
					'saddlebrown':          '#8B4513',
					'salmon':               '#FA8072',
					'sandybrown':           '#FAA460',
					'seagreen':             '#2E8B57',
					'seashell':             '#FFF5EE',
					'sienna':               '#A0522D',
					'silver':               '#C0C0C0',
					'skyblue':              '#87CEEB',
					'slateblue':            '#6A5ACD',
					'slategray':            '#708090',
					'snow':                 '#FFFAFA',
					'springgreen':          '#00FF7F',
					'steelblue':            '#4682B4',
					'tan':                  '#D2B48C',
					'teal':                 '#008080',
					'thistle':              '#D8BFD8',
					'tomato':               '#FF6347',
					'turquoise':            '#40E0D0',
					'violet':               '#EE82EE',
					'wheat':                '#F5DEB3',
					'white':                '#FFFFFF',
					'whitesmoke':           '#F5F5F5',
					'yellow':               '#FFFF00',
					'yellowgreen':          '#9ACD32'}
		return cnames
	
	def setColumns(self, data=None):
		self.plotsOlv.SetColumns([
			#imageGetter: "icon" is a property of plotListEntry objects
			ColumnDefn("Title", "left", 500, "title", imageGetter="icon")
			])
		#print "hallo"

	def refreshList(self):
		self.plotsOlv.ClearAll()
		self.plots = []
		self.pickedColors = []
		self.colors = self.setupColors()
		self.getPlots()
		self.setColumns()
		self.plotsOlv.CreateCheckStateColumn()
		self.plotsOlv.SetObjects(self.plots)
		self.checkedPlots = {}

	def getPlots(self):
		for idx, plot in enumerate(stored.plots):
			#print self.colors
			if not "Color" in plot:
				allowedColor = False
				preferredCnames = {
					'blue':		'#0000FF',
					'green':	'#008000',
					'red':		'#FF00FF',
					'orange':	'#FFA500',
					'violet':	'#EE82EE',
					'magenta':	'#FF0000'}
				while not allowedColor:
					if len(preferredCnames) > 0:
						color = preferredCnames.popitem()[1]
					else:
						color = self.colors.popitem()[1]
					if not color in self.pickedColors:
						allowedColor = True
					self.pickedColors.append(color)
					#print color
				plot['Color'] = color
			else:
				pass
			olv = self.plotsOlv
			self.plots.append(PlotListEntry(olv, plot['Type'],
											plot['Title'], plot['Color'], idx))
		#print self.plots

	def checkboxClicked(self, data):
		data.rowModel.toggleState()
		index = data.rowModel.index
		isChecked = data.rowModel.isChecked
		currentPlot = stored.plots[index]
		currentId = index
		currentType = currentPlot["Type"]
		allowed = True
		nonSuperposeableTypes = ["DifferenceMap","DistanceMap"]
		if len(self.checkedPlots) == 0 and isChecked:
			#this is the first checked plot, just add it
			self.checkedPlots[currentId] = currentPlot
			#print "Added first plot"
		elif len(self.checkedPlots) > 0 and isChecked:
			#check if of same type as other checked plots
			for key, value in self.checkedPlots.iteritems():
				#value will be a graph dict from mtsslWizard
				if currentType != value["Type"] or currentType in nonSuperposeableTypes:
					allowed = False
					#print "False, cannot add plot"
			if allowed:
				self.checkedPlots[currentId] = currentPlot
				#print "True, plot added"
			else:
				#self.lb.Check(index, check=False)
				print "New:"
				print self.plotsOlv.GetCheckState(data.rowModel)
				self.plotsOlv.SetCheckState(data.rowModel, False)
				#data.rowModel.setCheckedState(False)
				pub.sendMessage('update.statusBar', data = "Plots are not of same type or not superposeable!")
				print self.plotsOlv.GetCheckState(data.rowModel)
				self.plotsOlv.RefreshObject(data.rowModel)
		elif len(self.checkedPlots) > 0 and not isChecked:
			self.checkedPlots.pop(currentId, None)
			#print "removed plot from plot list"
		pub.sendMessage('item.checked', data = self.checkedPlots)

##########################################################################################
## Plotter (right Panel)
##########################################################################################

class Plotter(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent=parent)
		self.initGui()
		self.Centre()
		self.Show()
		self.graphsToPlot = []
		self.currentGraph = {}
		#listen to messages from checklistbox
		pub.subscribe(self.__onItemChecked, 'item.checked')

	#the list of plots changed
	def __onItemChecked(self, data):
		self.graphsToPlot = []
		for key, value in data.iteritems():
			self.graphsToPlot.append(value)
		self.plotGraphs()

	def initGui(self):
		# Create main panel
		panel = self
		panel.SetFont( wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
		# Set sizers
		sizer0 = wx.BoxSizer(wx.VERTICAL)
		panel.SetSizer(sizer0)
		# Matplotlib item
		panelcolour = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE).Get()
		facecolor = [channel/float(255) for channel in panelcolour]
		self.figure = plt.figure(figsize=(-1, -1), dpi=80, facecolor=facecolor, linewidth=1.0)#Figure(figsize=(-1, -1), dpi=80, facecolor=facecolor, linewidth=1.0)
		self.figure.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.98)
		self.canvas = FigCanvas(panel, -1, self.figure)
		self.mainAxes = plt.subplot2grid((1,1),(0, 0))
		self.mainAxes.set_xlabel('X-Axis')
		self.mainAxes.set_ylabel('Y-Axis')
		self.mainAxes.set_xlim(0.0, 1.0)
		self.mainAxes.set_ylim(0.0, 1.0)
		matplotlib.rcParams.update({'font.size': 12})
		self.mainAxes.grid(False)
		sizer0.Add(self.canvas, 1, flag=wx.GROW|wx.LEFT|wx.RIGHT, border=20)
		self.toolbar = NavigationToolbar(self.canvas)
		sizer0.Add(self.toolbar, 0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.LEFT|wx.RIGHT, border=20)
		sizer0.Add((-1, 10))
		self.mouseMove = self.canvas.mpl_connect('motion_notify_event', self.updateStatusBar)
		self.onKey = self.canvas.mpl_connect('key_press_event', self.onKey)
		self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
		self.plotLogo()
		# Redefine close event
		#self.Bind(wx.EVT_CLOSE, self.OnClose)

	def ChangeCursor(self,event):
		self.canvas.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

	def onKey(self, event):
		if event.inaxes and event.key == "s":
			x, y = event.xdata, event.ydata
			try:
				cmd.delete("mtsslPlot")
				cmd.select("mtsslPlot", "resi %i or resi %i" %(x,y))
			except:
				print "Could not select: resi %i or resi %i" %(x, y)

	def updateStatusBar(self, event):
		#could be made faster using animation blits:
		#http://scipy-cookbook.readthedocs.org/items/Matplotlib_Animations.html
		if event.inaxes: 
			x, y = event.xdata, event.ydata
			if event.inaxes == self.mainAxes:
				if self.currentGraph["Type"] == "DistanceMap" or self.currentGraph["Type"] == "DifferenceMap":
					z = self.currentGraph["zData"][x,y]
					statusString = "x = %i; y = %i; z = %i" %(x, y, z)
					pub.sendMessage('update.statusBar', data = statusString)
				
				#plot crosshair
				try:
					self.mainCrosshairX.remove()
					self.mainCrosshairY.remove()
					self.upperCrosshairY.remove()
					self.rightCrosshairX.remove()
				except:
					pass
				try:
					self.mainCrosshairX = self.mainAxes.axhline(y, color = "red")
					self.mainCrosshairY = self.mainAxes.axvline(x, color = "red")
					self.upperCrosshairY = self.axes1.axvline(x, color = "red")
					self.rightCrosshairX = self.axes3.axhline(y, color = "red")
					self.figure.tight_layout()		
					self.canvas.draw()
				except:
					pass
			else:
				statusString = "x = %i; y = %i" %(x, y)
				pub.sendMessage('update.statusBar', data = statusString)

	def plotGraphs(self):
		# Which type of graph is it
		if len(self.graphsToPlot) > 0:
			graphType = self.graphsToPlot[0]["Type"]
			if (graphType == 'DistanceMap'):
				self.plotDistanceMap(self.graphsToPlot)
			if (graphType == 'DifferenceMap'):
				self.plotDifferenceMap(self.graphsToPlot)
			if (graphType == 'DistancePlot'):
				self.plotDistancePlot(self.graphsToPlot)
			if (graphType == 'AccessibilityPlot'):
				self.plotAccessibilityPlot(self.graphsToPlot)
			if (graphType == 'DistanceDistribution'):
				self.plotDistanceDistribution(self.graphsToPlot)
			self.currentGraph = self.graphsToPlot[0]
			pub.sendMessage('update.statusBar', data = "Data plotted.")
		else:
			self.figure.clf()
			self.plotLogo()
			self.canvas.draw()
			pub.sendMessage('update.statusBar', data = "Nothing to plot. Please select a plot.")

	def refreshFigure(self):
		self.figure.clear()
		self.mainAxes.set_xlabel('X-Axis')
		self.mainAxes.set_ylabel('Y-Axis')
		self.mainAxes.set_xlim(0.0, 1.0)
		self.mainAxes.set_ylim(0.0, 1.0)
		matplotlib.rcParams.update({'font.size': 12})

	def plotDistanceMap(self, graphs):
		#there should only be one graph in graphs
		graph = graphs[0]
		# Read on x,y,z
		x = graph['xData'] - 0.5 * numpy.ones(len(graph['xData']))
		y = graph['yData'] - 0.5 * numpy.ones(len(graph['yData']))
		id = graph["id"]
		accessibilityPlots = []
		for plot in stored.plots:
			try:
				if id == plot["id"] and plot["Type"] == "AccessibilityPlot":
					accessibilityPlots.append(plot)
			except:
				pass
		#print accessibilityPlots
		acc1x = accessibilityPlots[0]["xData"]
		acc1y = accessibilityPlots[0]["yData"]		
		acc2x = accessibilityPlots[1]["xData"]
		acc2y = accessibilityPlots[1]["yData"]		
				
		xlim = graph['xlim']
		ylim = graph['ylim']
		X, Y = numpy.meshgrid(x, y)
		Z = graph['zData']
		
		# Define colormap
		cmap = colors.ListedColormap(['blue', 'green', 'yellow', 'orange', 'red'])
		cmap.set_under('white')
		cmap.set_over('red')
		cmap.set_bad(color = "grey")
		bounds = [1,15,40,60,80,100]
		norm = colors.BoundaryNorm(bounds, cmap.N)
		
		# Draw plot
		self.mainAxes = plt.subplot2grid((5,5),(1,0), colspan=4, rowspan=4)
		self.axes1 = plt.subplot2grid((5,5),(0,0), colspan=4, sharex=self.mainAxes)
		self.axes3 = plt.subplot2grid((5,5),(1,4), rowspan=4, sharey=self.mainAxes)
		self.axes4 =  plt.subplot2grid((5,5),(0,4))
		self.axes1.plot(acc1x, acc1y, 
				linestyle='solid',
				linewidth=0.5, 
				color='dimgrey',  
				markersize=6, 
				markeredgewidth=1)
		img = self.mainAxes.pcolormesh(x, y, Z, cmap=cmap, norm=norm)
		self.mainAxes.set_xlim(xlim)
		self.mainAxes.set_ylim(ylim)
		self.mainAxes.set_xlabel(graph['xTitle'])
		self.mainAxes.set_ylabel(graph['yTitle'])
		self.axes3.plot(acc2y, acc2x, 
				linestyle='solid',
				linewidth=0.5, 
				color='dimgrey', 
				markersize=6, 
				markeredgewidth=1)
		
		# Cosmetics
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes1.xaxis.set_ticks_position('top')
		self.axes1.set_yticks(self.axes1.get_ylim())
		self.axes1.set_ylabel("Access.")
		self.axes3.yaxis.set_ticks_position('right')
		self.axes3.set_xticks(self.axes1.get_ylim())
		self.axes3.set_xlabel("Access.")
		# Draw colorbar
		colorbar = plt.colorbar(img, cax=self.axes4)
		colorbar.ax.set_xlabel('Angstrom')
		colorbar.ax.xaxis.set_label_position('top')
		#colorbar.ax.xaxis.labelpad = 20
		self.figure.tight_layout()	
		self.canvas.draw()

	def plotDifferenceMap(self, graphs):
		#there should only be one graph in graphs
		graph = graphs[0]
		# Read on x,y,z
		x = graph['xData'] - 0.5 * numpy.ones(len(graph['xData']))
		y = graph['yData'] - 0.5 * numpy.ones(len(graph['yData']))
		id = graph["id"]
		accessibilityPlots = []
		for plot in stored.plots:
			if id == plot["id"] and plot["Type"] == "AccessibilityPlot":
				accessibilityPlots.append(plot)
		#print accessibilityPlots
		acc1x = accessibilityPlots[0]["xData"]
		acc1y = accessibilityPlots[0]["yData"]		
		acc2x = accessibilityPlots[1]["xData"]
		acc2y = accessibilityPlots[1]["yData"]

		X, Y = numpy.meshgrid(x, y)
		Z = graph['zData']	   
		xlim = graph['xlim']
		ylim = graph['ylim']
		# Draw surface plot
		self.mainAxes = plt.subplot2grid((5,5),(1,0), colspan=4, rowspan=4)
		self.axes1 = plt.subplot2grid((5,5),(0,0), colspan=4, sharex=self.mainAxes)
		self.axes3 = plt.subplot2grid((5,5),(1,4), rowspan=4, sharey=self.mainAxes)
		self.axes4 =  plt.subplot2grid((5,5),(0,4))
		self.axes1.plot(acc1x, acc1y, 
				linestyle='solid',
				linewidth=0.5, 
				color='dimgrey',  
				markersize=6, 
				markeredgewidth=1)
		try:
			plt.register_cmap(name='plasma', cmap=cmaps.plasma)
			plt.set_cmap(cmaps.plasma)
			cmap = plt.get_cmap(cmaps.plasma)
		except:
			cmap = matplotlib.cm.jet
		cmap.set_bad(color = "grey")
		img = self.mainAxes.pcolormesh(X, Y, Z, cmap=cmap, vmin=0.0, vmax=numpy.amax(Z))
		self.mainAxes.set_xlim(xlim)
		self.mainAxes.set_ylim(ylim)
		self.mainAxes.set_xlabel(graph['xTitle'])
		self.mainAxes.set_ylabel(graph['yTitle'])
		self.axes3.plot(acc2y, acc2x, 
				linestyle='solid',
				linewidth=0.5, 
				color='dimgrey', 
				markersize=6, 
				markeredgewidth=1)
		# Cosmetics
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.mainAxes.xaxis.set_minor_locator(xminorLocator)
		self.mainAxes.yaxis.set_minor_locator(yminorLocator)
		self.mainAxes.tick_params(direction='out', length=6, width=1)
		self.mainAxes.tick_params(which='minor', direction='out', length=3, width=1)
		self.mainAxes.xaxis.labelpad = 15
		self.mainAxes.yaxis.labelpad = 15
		self.axes1.xaxis.set_ticks_position('top')
		self.axes1.set_yticks(self.axes1.get_ylim())
		self.axes1.set_ylabel("Access.")
		self.axes3.yaxis.set_ticks_position('right')
		self.axes3.set_xticks(self.axes1.get_ylim())
		self.axes3.set_xlabel("Access.")
		# Draw colorbar
		colorbar = plt.colorbar(img, cax=self.axes4)
		colorbar.ax.set_xlabel('Angstrom')
		colorbar.ax.xaxis.set_label_position('top')
		#colorbar.ax.xaxis.labelpad = 20
		self.figure.tight_layout()		
		self.canvas.draw()

	def plotDistancePlot(self, graphs):	
		# Plot graph
		axes = plt.subplot2grid((1,1),(0,0))
		for graph in graphs:
			axes.plot(graph['xData'], graph['yData'], 
						   linestyle='solid',
						   linewidth=3, 
						   color=graph['Color'], 
						   #marker='o',
						   #markerfacecolor='white', 
						   #markeredgecolor='black', 
						   #markersize=3.5, 
						   #markeredgewidth=1
						   )
		Xmin = 0
		Xmax = graph['xData'].max()+1
		Ymin = 0
		Ymax = 100
		axes.set_xlim(graph['xlim'])
		axes.set_ylim(Ymin, Ymax)
		axes.set_xlabel(graph['xTitle'])
		axes.set_ylabel(graph['yTitle'])
		# Color the distance ranges
		axes.fill_between([Xmin, Xmax], 15, Ymin, color='blue', alpha=0.15)
		axes.fill_between([Xmin, Xmax], 40, 15, color='green', alpha=0.15)
		axes.fill_between([Xmin, Xmax], 60, 40, color='yellow', alpha=0.15)
		axes.fill_between([Xmin, Xmax], 80, 60, color='orange', alpha=0.15)
		axes.fill_between([Xmin, Xmax], Ymax, 80, color='red', alpha=0.15)
		# Cosmetics
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		axes.xaxis.set_minor_locator(xminorLocator)
		axes.yaxis.set_minor_locator(yminorLocator)
		axes.tick_params(direction='out', length=6, width=1)
		axes.tick_params(which='minor', direction='out', length=3, width=1)
		axes.xaxis.labelpad = 15
		axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()

	def plotDistanceDistribution(self, graphs):
		axes = plt.subplot2grid((1,1),(0,0))
		for graph in graphs:
			index=numpy.zeros(1)
			for i in range (0,graph['yData'].shape[0]):
				if (graph['yData'][i] <= 0):
					index = numpy.append(index,i)
					index = numpy.delete(index, 0)
			Y = numpy.delete(graph['yData'], index)
			X = numpy.delete(graph['xData'], index)
			axes.plot(X, Y, 
						   linestyle='solid',
						   linewidth=3, 
						   color=graph['Color'], 
						   #markerfacecolor='white', 
						   #markeredgecolor='black', 
						   markersize=3.5, 
						   #markeredgewidth=1
						   )
		Xmin = 0
		Xmax = graph['xData'].max()+1
		Ymin = 0
		Ymax = graph['yData'].max() * 1.05
		axes.set_xlim(graph['xlim'])
		axes.set_ylim(Ymin, Ymax)
		axes.set_xlabel(graph['xTitle'])
		axes.set_ylabel(graph['yTitle'])
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		axes.xaxis.set_minor_locator(xminorLocator)
		axes.yaxis.set_minor_locator(yminorLocator)
		axes.tick_params(direction='out', length=6, width=1)
		axes.tick_params(which='minor', direction='out', length=3, width=1)
		axes.xaxis.labelpad = 15
		axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()

	def plotAccessibilityPlot(self, graphs):
		axes = plt.subplot2grid((1,1),(0,0))
		for graph in graphs:
			axes.plot(graph['xData'], graph['yData'], 
						   linestyle='solid',
						   linewidth=3, 
						   color=graph['Color'], 
						   marker='.',
						   #markerfacecolor='white', 
						   #markeredgecolor='black', 
						   markersize=3.5, 
						   #markeredgewidth=1
						   )
		Xmin = 0
		Xmax = graph['xData'].max()+1
		Ymin = 0
		Ymax = graph['yData'].max() * 1.05
		axes.set_xlim(graph['xlim'])
		axes.set_ylim(Ymin, Ymax)
		axes.set_xlabel(graph['xTitle'])
		axes.set_ylabel(graph['yTitle'])
		# Cosmetics
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		axes.xaxis.set_minor_locator(xminorLocator)
		axes.yaxis.set_minor_locator(yminorLocator)
		axes.tick_params(direction='out', length=6, width=1)
		axes.tick_params(which='minor', direction='out', length=3, width=1)
		axes.xaxis.labelpad = 15
		axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()
	
	def plotLogo(self):
		axes = plt.subplot2grid((1,1),(0,0))
		img = mpimg.imread(os.path.dirname(__file__)+"/logo.png")
		axes.imshow(img)