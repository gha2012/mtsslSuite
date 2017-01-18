# -*- coding: utf-8 -*-
import os
import sys
import threading
from threading import Thread
import string
import math
import numpy
import copy
from operator import attrgetter
import pymol
from protein import Protein
from pymol import cmd, stored
from docker import Docker
import wx
import wx.lib.scrolledpanel
import cPickle as pickle
import cProfile
import EnhancedStatusBar as ESB
from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub
from ObjectListView import ObjectListView, ColumnDefn, OLVEvent
import wx.lib.agw.flatnotebook as fnb
from random import randint, uniform

##########################################################################################

class MainWindow(wx.Frame):
	""" Main GUI of mtsslDock"""
	#-------------------------------------------------------------------------------------

	def __init__(self, *args, **kwds):
		#setup GUI
		wx.Frame.__init__(self, None, -1, "mtsslDock",	 wx.DefaultPosition, wx.Size(950,768))
		self.panel = wx.Panel(self)
		self.notebook_1 = fnb.FlatNotebook(self.panel, wx.ID_ANY, agwStyle=fnb.FNB_NODRAG|fnb.FNB_FF2)
		self.notebook_1.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onNotebookPageChanged)
		self.notebook_1.Bind(fnb.EVT_FLATNOTEBOOK_PAGE_CLOSING, self.onNotebookPageClosing)
		self.notebook_1_pane_1 = wx.Panel(self.notebook_1, wx.ID_ANY)
		self.bounding_box_1 = wx.StaticBox(self.notebook_1_pane_1, label='Settings')
		
		currentPymolObjects = cmd.get_object_list('(all)')
		self.label_1 = wx.StaticText(self.notebook_1_pane_1, wx.ID_ANY, "Fixed molecule:")
		self.combo_box_1 = wx.ComboBox(self.notebook_1_pane_1, 0, style=wx.CB_READONLY, choices=currentPymolObjects)
		self.label_2 = wx.StaticText(self.notebook_1_pane_1, wx.ID_ANY, "Moving molecule:")
		self.combo_box_2 = wx.ComboBox(self.notebook_1_pane_1, 0, style=wx.CB_READONLY, choices=currentPymolObjects)		
		self.label_3 = wx.StaticText(self.notebook_1_pane_1, wx.ID_ANY, "Score clashes:")
		self.combo_box_3 = wx.ComboBox(self.notebook_1_pane_1, 0, style=wx.CB_READONLY, choices=["Yes", "No"])
		self.label_4 = wx.StaticText(self.notebook_1_pane_1, wx.ID_ANY, "Symmetry:")
		self.combo_box_4 = wx.ComboBox(self.notebook_1_pane_1, 0, style=wx.CB_READONLY, choices=["None", "C2", "C3", "C4", "C5", "C6", "C7", "C8"])

		self.scrollingPanel = wx.lib.scrolledpanel.ScrolledPanel(self.notebook_1_pane_1)
		self.scrollingPanel.SetupScrolling()
		self.scrollingPanel.SetAutoLayout(1)
		self.bounding_box_2 = wx.StaticBox(self.scrollingPanel, label='Restraints')
		self.statusBar = ESB.EnhancedStatusBar(self, -1)
		self.statusBar.SetFieldsCount(3)
		self.statusBar.SetStatusWidths([-1, 188, 1])
		self.progress = wx.Gauge(self.statusBar, range=100)
		self.statusBar.AddWidget(self.progress, ESB.ESB_EXACT_FIT, ESB.ESB_EXACT_FIT, pos=1)
		self.SetStatusBar(self.statusBar)
		
		self.__set_properties()
		self.__do_layout()
		self.SetMinSize(wx.Size(600,600))
		self.initMenu()
		self._init_toolbar()

		pub.subscribe(self.onDockingResultReady, "docking.ready")
		pub.subscribe(self.onDockingUpdate, "docking.update")
		pub.subscribe(self.onRestraintSelected, "restraint.selected")
		pub.subscribe(self.onRestraintDeselected, "restraint.deselected")
		
		self.restraintViews = []
		self.selectedRestraintViews = []
		self.resultsPages = []
		self.runNumber = 1
		self.anchorNumber = 1
		#self.onAddRestraint(None)
		self.statusBar.SetStatusText("Idle.")	
		
	#-------------------------------------------------------------------------------------

	def _init_toolbar(self):
		iconsDir = os.path.dirname(__file__) + "/icons"
		self.toolbar = self.CreateToolBar(style=(wx.TB_TEXT))

		dockIcon = wx.Bitmap(os.path.join(iconsDir, "start.png"))
		dockButton = self.toolbar.AddLabelTool(2, "Start docking run", dockIcon, longHelp="Start docking run using the current set of restraints.")
		stopDockIcon = wx.Bitmap(os.path.join(iconsDir, "stop.png"))
		stopDockButton = self.toolbar.AddLabelTool(5, "Abort docking run", stopDockIcon, longHelp="Abort the currently running docking job.")
		addRestraintIcon = wx.Bitmap(os.path.join(iconsDir, "addRestraint.png"))
		addRestraintButton = self.toolbar.AddLabelTool(1, "Add restraint", addRestraintIcon, longHelp="Add a new restraint.")
		deleteRestraintIcon = wx.Bitmap(os.path.join(iconsDir, "deleteRestraint.png"))
		deleteRestraintButton = self.toolbar.AddLabelTool(7, "Delete restraints", deleteRestraintIcon, longHelp="Delete selected restraints.")
		addAnchorIcon = wx.Bitmap(os.path.join(iconsDir, "addAnchor.png"))
		addAnchorButton = self.toolbar.AddLabelTool(8, "Add anchor", addAnchorIcon, longHelp="Add a new anchor point at current center of view in PyMOL.")
		createInPymolIcon = wx.Bitmap(os.path.join(iconsDir, "pymol.png"))
		createInPymolButton = self.toolbar.AddLabelTool(6, "Create in PyMOL", createInPymolIcon, longHelp="Create solutions of this run in Pymol.")
		synchroniseIcon = wx.Bitmap(os.path.join(iconsDir, "synchronise.png"))
		synchroniseButton = self.toolbar.AddLabelTool(9, "Sync", synchroniseIcon, longHelp="Synchronise GUI with PyMOL.")
	
		self.Bind(wx.EVT_TOOL, self.onAddRestraint, addRestraintButton)
		self.Bind(wx.EVT_TOOL, self.onDeleteRestraint, deleteRestraintButton)
		self.Bind(wx.EVT_TOOL, self.onAddAnchor, addAnchorButton)
		self.Bind(wx.EVT_TOOL, self.startDockingThread, dockButton)
		self.Bind(wx.EVT_TOOL, self.stopDockingThread, stopDockButton)
		self.Bind(wx.EVT_TOOL, self.createInPymolButtonClicked, createInPymolButton)
		self.Bind(wx.EVT_TOOL, self.onRefresh, synchroniseButton)
		self.toolbar.Realize()
		self.toolbar.EnableTool(6, False)

	#-------------------------------------------------------------------------------------
		
	def initMenu(self):
		menubar = wx.MenuBar()
		#File menu
		fileMenu = wx.Menu()
		loadRestraintsItem = wx.MenuItem(fileMenu, 11, '&Load restraints...\tCtrl+O')
		self.Bind(wx.EVT_MENU, self.onLoadRestraints, id=11)
		saveRestraintsItem = wx.MenuItem(fileMenu, 12, '&Save restraints...\tCtrl+S')
		self.Bind(wx.EVT_MENU, self.onSaveRestraints, id=12)
		refreshItem = wx.MenuItem(fileMenu, 13, '&Refresh\tCtrl+R')
		self.Bind(wx.EVT_MENU, self.onRefresh, id=13)
		quitItem = wx.MenuItem(fileMenu, 14, '&Quit\tCtrl+Q')
		self.Bind(wx.EVT_MENU, self.onQuit, id=14)
		loadSessionItem = wx.MenuItem(fileMenu, 15, '&Load session...\tCtrl+Shft+O')
		self.Bind(wx.EVT_MENU, self.onLoadSession, id=15)
		saveSessionItem = wx.MenuItem(fileMenu, 16, '&Save session...\tCtrl+Shft+S')
		self.Bind(wx.EVT_MENU, self.onSaveSession, id=16)
		fileMenu.AppendItem(loadRestraintsItem)
		fileMenu.AppendItem(saveRestraintsItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(loadSessionItem)
		fileMenu.AppendItem(saveSessionItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(refreshItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(quitItem)
		#dock menu
		dockMenu = wx.Menu()
		addRestraintItem = wx.MenuItem(dockMenu, 21, '&Add restraint\tCtrl+A')
		self.Bind(wx.EVT_MENU, self.onAddRestraint, id=21)
		deleteRestraintItem = wx.MenuItem(dockMenu, 25, '&Delete selected restraints\tCtrl+A')
		self.Bind(wx.EVT_MENU, self.onDeleteRestraint, id=25)
		startDockingItem = wx.MenuItem(dockMenu, 22, '&Start docking run\tCtrl+D')
		self.Bind(wx.EVT_MENU, self.startDockingThread, id=22)
		stopDockingItem = wx.MenuItem(dockMenu, 23, '&Abort docking run\tCtrl+Z')
		self.Bind(wx.EVT_MENU, self.stopDockingThread, id=23)
		showSettingsItem = wx.MenuItem(dockMenu, 24, '&Advanced settings...')
		self.Bind(wx.EVT_MENU, self.onShowSettingsDialog, id=24)
		dockMenu.AppendItem(addRestraintItem)
		dockMenu.AppendItem(deleteRestraintItem)
		dockMenu.AppendSeparator()
		dockMenu.AppendItem(startDockingItem)
		dockMenu.AppendItem(stopDockingItem)
		dockMenu.AppendSeparator()
		dockMenu.AppendItem(showSettingsItem)
		# Help menu
		helpMenu = wx.Menu()
		aboutItem = wx.MenuItem(helpMenu, 31, '&About\tCtrl+A')
		self.Bind(wx.EVT_MENU, self.onAbout, id=31)
		helpMenu.AppendItem(aboutItem)
		
		menubar.Append(fileMenu, '&File')
		menubar.Append(dockMenu, '&Dock')
		menubar.Append(helpMenu, "&Help")
		self.SetMenuBar(menubar)

	#-------------------------------------------------------------------------------------

	def __set_properties(self):
		self.SetTitle("mtsslDock")
		self.settings = {}
		self.settings['numberOfPopulations'] = 10
		self.settings['numberOfGenerations'] = 1000
		self.settings['numberOfChromosomes'] = 400
		self.settings['numberOfRigidBodyCycles'] = 50
		self.settings['scoreClashes'] = True
		self.settings['scoreCOGdiff'] = True
		self.settings['symmetry'] = "None"
		if self.settings["scoreClashes"]:
			self.combo_box_3.SetValue("Yes")
		else:
			self.combo_box_3.SetValue("No")
		self.combo_box_4.SetValue(self.settings["symmetry"])


	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		self.mainSizer = wx.BoxSizer(wx.VERTICAL)
		notebookSizer = wx.BoxSizer(wx.HORIZONTAL)
		settingsSizer = wx.StaticBoxSizer(self.bounding_box_1, wx.VERTICAL)
		
		self.scrollingPanelSizer = wx.StaticBoxSizer(self.bounding_box_2, wx.VERTICAL)
		self.scrollingPanel.SetSizer(self.scrollingPanelSizer)
		
		settingsSizer.Add(self.label_1, 0,flag = wx.ALL, border=7)
		settingsSizer.Add(self.combo_box_1, 0, flag = wx.ALL, border=3)
		settingsSizer.Add(self.label_2, 0, flag = wx.ALL, border=7)
		settingsSizer.Add(self.combo_box_2, flag = wx.ALL, border=3)
		settingsSizer.Add(self.label_3, 0,flag = wx.ALL, border=7)
		settingsSizer.Add(self.combo_box_3, 0, flag = wx.ALL, border=3)
		settingsSizer.Add(self.label_4, 0, flag = wx.ALL, border=7)
		settingsSizer.Add(self.combo_box_4, flag = wx.ALL, border=3)

		notebookSizer.Add(settingsSizer, 0, flag = wx.EXPAND | wx.ALL, border=7)
		notebookSizer.Add(self.scrollingPanel, 1, flag = wx.EXPAND | wx.ALL, border=7)
		self.notebook_1_pane_1.SetSizer(notebookSizer)
		self.notebook_1.AddPage(self.notebook_1_pane_1, "Docking")
		self.mainSizer.Add(self.notebook_1, 1, flag = wx.EXPAND | wx.ALL, border=7)
		self.panel.SetSizer(self.mainSizer)
		self.Layout()

	#-------------------------------------------------------------------------------------

	def onAbout(self, event):
		message = "mtsslDock 2.0\nG. Hagelueken\n2016"
		caption = "About"
		dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.ICON_INFORMATION)
		dlg.ShowModal()
		dlg.Destroy()

	#-------------------------------------------------------------------------------------

	def onRestraintSelected(self, index):
		restraintView = self.restraintViews[index-1]
		self.selectedRestraintViews.append(restraintView)
		#print self.selectedRestraintViews

	#-------------------------------------------------------------------------------------
	
	def onRestraintDeselected(self, index):
		for idx, restraintView in enumerate(self.selectedRestraintViews):
			if restraintView.id == index:
				self.selectedRestraintViews.pop(idx)
		#print self.selectedRestraintViews

	#-------------------------------------------------------------------------------------

	def onAddAnchor(self, event):
		objectName = "mD_anchor_%i" %self.anchorNumber
		cmd.pseudoatom(object=objectName, name="CA")
		cmd.show("spheres", "mD_anchor_%i" %self.anchorNumber)
		self.anchorNumber += 1

	#-------------------------------------------------------------------------------------

	def onDeleteRestraint(self, id):
		selectedIds = []
		for restraintView in self.selectedRestraintViews:
			selectedIds.append(restraintView.id)
		for id in selectedIds:
			for idx, restraintView in enumerate(self.restraintViews):
				if restraintView.id == id:
					restraintViewToDelete = self.restraintViews.pop(idx)
					restraintViewToDelete.Destroy()
				else:
					pass
		self.selectedRestraintViews = []
		#renumber remaining restraintViews
		for idx, restraintView in enumerate(self.restraintViews):
			restraintView.id = idx + 1
			restraintView.label_0.SetLabel("Restraint: %i" %(idx+1))
		self.scrollingPanel.Layout()
		self.scrollingPanel.FitInside()

	#-------------------------------------------------------------------------------------

	def onLoadSession(self, event):
		results = {}
		dlg = wx.FileDialog(self, message="Choose a session file", defaultDir=os.getcwd(), defaultFile="", wildcard="*.session", style=wx.OPEN | wx.CHANGE_DIR)
		#step 1: restore restraints
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
			print "Unpickling restraints...",
			sessionDict = pickle.load(open(path, 'rb'))
			print "Done!\n"
		dlg.Destroy()
		self.onLoadRestraints(event, sessionDict["restraints"])
		self.combo_box_1.SetValue(sessionDict["restraints"][0]["proteinAname"])
		self.combo_box_2.SetValue(sessionDict["restraints"][0]["proteinBname"])
		#step 2: restore settings
		self.settings = sessionDict["settings"]
		#step 3: restore resultsPages
		pickledResultsPages = sessionDict["resultsPages"]
		for pickledResultsPage in pickledResultsPages:
			self.onDockingResultReady({"dockingResult":pickledResultsPage["dockingResults"],
										"settings":pickledResultsPage["settings"]})
		self.statusBar.SetStatusText("Session loaded!")

	#-------------------------------------------------------------------------------------

	def onSaveSession(self, event):
		sessionDict = {}
		resultsPages = []
		for resultsPage in self.resultsPages:
			resultsPages.append({"dockingResults":resultsPage.dockingResults,
										"settings":resultsPage.settings})
			
		sessionDict["restraints"] = self.prepareRestraintsForDocking()
		sessionDict["resultsPages"] = resultsPages
		sessionDict["settings"] = self.settings
		 
		dialog = wx.FileDialog(None, message="Save session as ...", defaultDir=os.getcwd(), defaultFile=".session", wildcard="*.session", style = wx.SAVE | wx.OVERWRITE_PROMPT)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			dataFile = open(path, 'w')
			print "\nPickling session...",
			pickle.dump(sessionDict, open(path, 'wb')) 
			print "Done!\n"
			dataFile.close()
		dialog.Destroy()
		message = wx.MessageDialog(None, "Don't forget to save the PyMOL session!", 'mtsslDock', wx.OK)
		message.ShowModal()
		message.Destroy()
		self.statusBar.SetStatusText("Session saved!")

	#-------------------------------------------------------------------------------------

	def onRefresh(self, event):
		#refresh currentPymolObjects in Dropdowns
		currentPymolObjects = cmd.get_object_list('(all)')
		self.combo_box_1.SetItems(currentPymolObjects)
		self.combo_box_2.SetItems(currentPymolObjects)
		currentPymolObjects = cmd.get_object_list('*_label or mD_anchor_*')
		for restraintView in self.restraintViews:
			restraintView.combo_box_1.SetItems(currentPymolObjects)
			restraintView.combo_box_2.SetItems(currentPymolObjects)
		
		#sync checkstate in resultslists
		for resultsPage in self.resultsPages:
			olv = resultsPage.topPanel.resultsOlv
			for solution in olv.GetObjects():
				if solution.isChecked:
					cmd.enable(solution.name)
				else:
					cmd.disable(solution.name)
		#if objName in cmd.get_names(enabled_only=1):
		#print "objName is enabled"

	#-------------------------------------------------------------------------------------
	
	def onShowSettingsDialog(self, event):
		dlg = SettingsDialog(self.settings)
		res = dlg.ShowModal()
		if res == wx.ID_OK:
			self.settings['numberOfPopulations'] = int(dlg.text_ctrl_4.GetValue())
			self.settings['numberOfGenerations'] = int(dlg.text_ctrl_2.GetValue())
			self.settings['numberOfChromosomes'] = int(dlg.text_ctrl_1.GetValue())
			self.settings['numberOfRigidBodyCycles'] = int(dlg.text_ctrl_3.GetValue())
			#self.settings['scoreClashes'] = dlg.combo_box_2.GetValue()
			#self.settings['scoreCOGdiff'] = dlg.combo_box_3.GetValue()
			#self.settings['symmetry'] = dlg.combo_box_1.GetValue()
		dlg.Destroy()

	#-------------------------------------------------------------------------------------
	
	def onQuit(self, event):
		dlg = wx.MessageDialog(self, 
			"Do you really want to close this mtsslDock?",
			"Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
		result = dlg.ShowModal()
		dlg.Destroy()
		if result == wx.ID_OK:
			self.Destroy()

	#-------------------------------------------------------------------------------------

	def onNotebookPageClosing(self, event):
		print "Page about to close..."
		tabid = event.GetSelection()
		if tabid == 0:
			event.Veto()
		else:
			self.resultsPages.pop(tabid-1)

	#-------------------------------------------------------------------------------------

	def onNotebookPageChanged(self, event):
		self.statusBar.SetStatusText("Idle.")
		tabid = self.notebook_1.GetSelection()
		if tabid > 0:
			self.toolbar.EnableTool(6, True)
			self.toolbar.EnableTool(1, False)
			self.toolbar.EnableTool(7, False)
			self.toolbar.EnableTool(8, False)
		else:
			try:
				self.toolbar.EnableTool(6, False)
				self.toolbar.EnableTool(1, True)
				self.toolbar.EnableTool(7, True)
				self.toolbar.EnableTool(8, True)
			except:
				pass

	#-------------------------------------------------------------------------------------
		
	def createInPymolButtonClicked(self, event):
		idx = self.notebook_1.GetSelection() - 1
		dockingResult = self.resultsPages[idx].dockingResults
		id = self.resultsPages[idx].id
		self.createSolutionsInPymol(dockingResult, id)

	#-------------------------------------------------------------------------------------

	def hsv_to_rgb(self, h, s, v):
		"""Converts HSV value to RGB values
		Hue is in range 0-359 (degrees), value/saturation are in range 0-1 (float)

		http://stackoverflow.com/questions/1586147/how-to-generate-random-greenish-colors

		Direct implementation of:
		http://en.wikipedia.org/wiki/HSL_and_HSV#Conversion_from_HSV_to_RGB
		"""
		h, s, v = [float(x) for x in (h, s, v)]

		hi = (h / 60) % 6
		hi = int(round(hi))

		f = (h / 60) - (h / 60)
		p = v * (1 - s)
		q = v * (1 - f * s)
		t = v * (1 - (1 - f) * s)

		if hi == 0:
			return v, t, p
		elif hi == 1:
			return q, v, p
		elif hi == 2:
			return p, v, t
		elif hi == 3:
			return p, q, v
		elif hi == 4:
			return t, p, v
		elif hi == 5:
			return v, p, q

	def colorObject(self, selection):
		h = randint(90, 140) # Select random green'ish hue from hue wheel
		s = uniform(0.2, 1)
		v = uniform(0.3, 1)
		r, g, b = self.hsv_to_rgb(h, s, v)
		cmd.set_color("random_color", [r, g, b])
		cmd.color("random_color", selection)

	def createSolutionsInPymol(self, dockingResults, id):
		try:
			objectPrefix = "mD"
			dockingRunNumber = id
			nonClashingSolution = 1
			clashingSolution = 1
			proteinA = dockingResults.proteinA
			proteinB = dockingResults.proteinB
			for population in dockingResults.populations:
				solution = population.chromosomes[0]
				if solution.clashes <= 5:
					nameOfSolution = "%s-%i_sol-%i" % (objectPrefix, dockingRunNumber, nonClashingSolution)
					solution.name = nameOfSolution
					proteinB.moveInPymol(nameOfSolution, solution, 1)
					cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 0, 0, None)
					if population.symmetry != "None":
						cmd.align(nameOfSolution, proteinA.pymolString, target_state=1, mobile_state=1)
					nonClashingSolution += 1

				elif solution.clashes > 5:
					nameOfSolution = "%s-%i_clash-%i" % (objectPrefix, dockingRunNumber, clashingSolution)
					solution.name = nameOfSolution
					proteinB.moveInPymol(nameOfSolution, solution, 1)
					cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 0, 0, None)
					if population.symmetry != "None":
						cmd.align(nameOfSolution, proteinA.pymolString, target_state=1, mobile_state=1)
					clashingSolution += 1
			cmd.group("%s-%i" % (objectPrefix, dockingRunNumber), "%s-%i*" % (objectPrefix, dockingRunNumber))
			self.colorObject("%s-%i" % (objectPrefix, dockingRunNumber))
			pub.sendMessage("add.pymol")
		except Exception,e:
			self.statusBar.SetStatusText("Cannot create solutions in PyMOL.")
			print str(e)

	#-------------------------------------------------------------------------------------

	def onDockingUpdate(self, progress):
		self.progress.SetValue(int(progress*100))

	#-------------------------------------------------------------------------------------

	def onLoadRestraints(self, event, restraints=[]):
		if len(restraints) == 0:
			dlg = wx.FileDialog(self, message="Choose a restraints file", defaultDir=os.getcwd(), defaultFile="", wildcard="*.restraints", style=wx.OPEN | wx.CHANGE_DIR)
			if dlg.ShowModal() == wx.ID_OK:
				path = dlg.GetPath()
				print "Unpickling restraints...",
				restraints = pickle.load(open(path, 'rb'))
				print "Done!\n"
			dlg.Destroy()
		currentNumberOfRestraints = 0
		for view in self.restraintViews:
			view.Destroy()
		self.restraintViews = []
		for restraint in restraints:
			newRestraintView = RestraintView(self.scrollingPanel, currentNumberOfRestraints + 1)
			newRestraintView.id = currentNumberOfRestraints + 1
			try:
				newRestraintView.label_0.SetValue(restraint["Name"])
			except:
				pass
			newRestraintView.combo_box_1.SetValue(restraint["anchorAname"])
			newRestraintView.combo_box_2.SetValue(restraint["anchorBname"])
			newRestraintView.text_ctrl_1.SetValue(str(restraint["distance"]))
			newRestraintView.text_ctrl_2.SetValue(str(restraint["width"]))
			try:
				coord = cmd.get_model(restraint["anchorAname"], 1).get_coord_list()
				newRestraintView.anchorAcoord = coord[0]
				coord = cmd.get_model(restraint["anchorBname"], 1).get_coord_list()
				newRestraintView.anchorBcoord = coord[0]
			except:
				self.statusBar.SetStatusText("Cannot get coordinates of anchors. Structure loaded?")
			self.restraintViews.append(newRestraintView)
			self.scrollingPanelSizer.Add(newRestraintView, 0, flag = wx.EXPAND | wx.ALL, border=10)
			currentNumberOfRestraints += 1
		self.scrollingPanel.Layout()
		self.scrollingPanel.FitInside()
		self.statusBar.SetStatusText("Restraints loaded!")

	#-------------------------------------------------------------------------------------

	def onSaveRestraints(self, event):
		restraints = self.prepareRestraintsForDocking()
		dialog = wx.FileDialog(None, message="Save restraints as ...", defaultDir=os.getcwd(), defaultFile=".restraints", wildcard="*.restraints", style = wx.SAVE | wx.OVERWRITE_PROMPT)
		if dialog.ShowModal() == wx.ID_OK:
			path = dialog.GetPath()
			dataFile = open(path, 'w')
			print "\nPickling restraints...",
			pickle.dump(restraints, open(path, 'wb')) 
			print "Done!\n"
			dataFile.close()
		dialog.Destroy()
		self.statusBar.SetStatusText("Restraints saved!")

	#-------------------------------------------------------------------------------------

	def onAddRestraint(self, event):
		currentNumberOfRestraints = len(self.restraintViews)
		newRestraintView = RestraintView(self.scrollingPanel, currentNumberOfRestraints + 1)
		self.restraintViews.append(newRestraintView)
		self.scrollingPanelSizer.Add(newRestraintView, 0, flag = wx.EXPAND | wx.ALL, border=10)
		#The two commands below refresh the scrolling panel so that scroll bars are drawn,
		#if necessary
		self.scrollingPanel.Layout()
		self.scrollingPanel.FitInside()

	#-------------------------------------------------------------------------------------

	def stopDockingThread(self, event):
		try:
			self.dockingThread.stop()
		except:
			pass
		self.progress.SetValue(0)
		self.toolbar.EnableTool(2, True)
		self.statusBar.SetStatusText("Run aborted!")

	#-------------------------------------------------------------------------------------

	def startDockingThread(self, event):
		self.statusBar.SetStatusText("Docking... please be patient!")
		self.toolbar.EnableTool(2, False)
		try:
			restraints = self.prepareRestraintsForDocking()
		except:
			self.statusBar.SetStatusText("Cannot prepare restraints. Something wrong with input.")
			self.stopDockingThread(None)
			return
		
		if len(self.restraintViews) > 0 and len(self.combo_box_1.GetValue()) > 0 and len(self.combo_box_2.GetValue()) > 0:
			if self.combo_box_3.GetValue() == "Yes":
				self.settings['scoreClashes'] = True
			elif self.combo_box_3.GetValue() == "No":
				self.settings['scoreClashes'] = False
			self.settings['symmetry'] = self.combo_box_4.GetValue()
			self.dockingThread = DockingThread(self.runNumber, restraints, self.settings)
		else:
			self.statusBar.SetStatusText("Cannot start docking run. Something wrong with input.")
			self.stopDockingThread(None)

	#-------------------------------------------------------------------------------------

	def onDockingResultReady(self, result):
		dockingResult = result["dockingResult"]
		settings = copy.deepcopy(result["settings"])
		restraints = result["restraints"]
		#print result["settings"]
		resultsPageContents = ResultsPageContents(self.notebook_1, self.runNumber, dockingResult, settings, restraints)
		self.notebook_1.AddPage(resultsPageContents, "Result-%i"%self.runNumber)
		self.resultsPages.append(resultsPageContents)
		self.notebook_1.Layout()
		self.runNumber += 1
		self.statusBar.SetStatusText("Docking run finished.")
		self.progress.SetValue(0)
		self.toolbar.EnableTool( 2, True )
		#for resultPage in self.resultsPages:
		#	print resultPage.settings
		#	print 

	#-------------------------------------------------------------------------------------

	def prepareRestraintsForDocking(self):
		restraints = []
		for restraintView in self.restraintViews:
			if restraintView.active:
				restraint = {}
				restraint["id"] = restraintView.id
				restraint["Name"] = restraintView.getRestraintName()
				restraint["proteinAname"] = self.combo_box_1.GetValue()
				restraint["proteinBname"] = self.combo_box_2.GetValue()
				restraint["anchorAname"] = restraintView.combo_box_1.GetValue()
				restraint["anchorBname"] = restraintView.combo_box_2.GetValue()
				restraint["anchorAcoord"] = restraintView.anchorAcoord
				restraint["anchorBcoord"] = restraintView.anchorBcoord
				restraint["distance"] = float(restraintView.text_ctrl_1.GetValue())
				restraint["width"] = float(restraintView.text_ctrl_2.GetValue())
				restraint["weight"] = 1.0
				stored.allCaAtomsA = []
				stored.allCaAtomsB = []
				cmd.iterate_state(1, "%s & name CA" % self.combo_box_1.GetValue(), 'stored.allCaAtomsA.append((x,y,z))')
				cmd.iterate_state(1, "%s & name CA" % self.combo_box_2.GetValue(), 'stored.allCaAtomsB.append((x,y,z))')
				restraint["proteinAcalpha"] = numpy.array(stored.allCaAtomsA)
				restraint["proteinBcalpha"] = numpy.array(stored.allCaAtomsB)
				restraints.append(restraint)
		return restraints

	#-------------------------------------------------------------------------------------

##########################################################################################

class ResultListEntry(object):

	"""On object of this class represents an entry in the list of docking results on a results
	page"""
	
	#-------------------------------------------------------------------------------------

	def __init__(self, olv, id, name, chi2, rmsd, clashes, isChecked):
		self.index = id
		self.name = name
		self.chi2 = chi2 
		self.rmsd = rmsd 
		self.clashes = clashes 
		self.isChecked = isChecked

	#-------------------------------------------------------------------------------------
	
	def toggleState(self):
		if self.isChecked == False:
			self.isChecked = True
		else:
			self.isChecked = False

	#-------------------------------------------------------------------------------------
	
	def setCheckedState(self, state):
		if state == False:
			self.isChecked = False
		elif state == True:
			self.isChecked = True

	#-------------------------------------------------------------------------------------

##########################################################################################

class ResultsTable(wx.Panel):

	"""Class that shows a table with detailed information for each docking result"""

	#-------------------------------------------------------------------------------------

	def __init__(self, parent, dockingResults, settings, restraints):
		wx.Panel.__init__(self, parent)
		self.restraints = restraints
		self.settings = settings
		self.dockingResults = dockingResults
		self.bounding_box = wx.StaticBox(self, label='Details for selected solution')
		self.grid_1 = MyGrid(self, wx.ID_ANY, wx.WANTS_CHARS)
		self.grid_1.SetRowLabelSize(0)
		self.radio_box_1 = wx.RadioBox(self, wx.ID_ANY, "Available result tables", choices=["Restraints vs docked distances", "Optimisation cycle vs parameters", "Settings"], majorDimension=1,)
		self.radio_box_1.Bind(wx.EVT_RADIOBOX, self.onRadioBoxChanged)
		self.__set_properties()
		self.__do_layout()
		pub.subscribe(self.onSelectionChanged, "update.selection")

	#-------------------------------------------------------------------------------------

	def __set_properties(self):
		try:
			index = self.selectedResult.index
		except:
			index = 0
		rows = len(self.dockingResults.populations[index].log.split("\n"))
		self.grid_1.CreateGrid(rows, 11)
		self.grid_1.SetColLabelValue(0, "Cycle")
		self.grid_1.SetColLabelValue(1, "rotX (°)")
		self.grid_1.SetColLabelValue(2, "rotY (°)")
		self.grid_1.SetColLabelValue(3, "rotZ (°)")
		self.grid_1.SetColLabelValue(4, u"tX (\u212B)")
		self.grid_1.SetColLabelValue(5, u"tY (\u212B)")
		self.grid_1.SetColLabelValue(6, u"tZ (\u212B)")
		self.grid_1.SetColLabelValue(7, u"R.m.s.d. (\u212B)")
		self.grid_1.SetColLabelValue(8, "Chi2")
		self.grid_1.SetColLabelValue(9, "Clashes")
		#self.grid_1.SetColLabelValue(10, u"COG (\u212B)")
		self.grid_1.SetColLabelValue(10, "Overall Score")
		self.grid_1.EnableEditing(False)

	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		sizer_1 = wx.StaticBoxSizer(self.bounding_box, wx.VERTICAL)
		sizer_1.Add(self.radio_box_1, 0, wx.EXPAND | wx.ALL, border=7)
		sizer_1.Add(self.grid_1, 1, wx.EXPAND | wx.ALL, border=7)
		self.SetSizer(sizer_1)

	#-------------------------------------------------------------------------------------

	def onSelectionChanged(self, selectedResult):
		self.selectedResult = selectedResult
		self.onRadioBoxChanged(None)
		

	#-------------------------------------------------------------------------------------

	def onRadioBoxChanged(self, event):
		selection = self.radio_box_1.GetSelection()
		index = self.selectedResult.index
		self.grid_1.ClearGrid()
		if selection == 2:
			self.grid_1.SetColLabelValue(0, "Parameter")
			self.grid_1.SetColLabelValue(1, "Value")
			self.grid_1.SetColLabelValue(2, "")
			self.grid_1.SetColLabelValue(3, "")
			self.grid_1.SetColLabelValue(4, "")
			self.grid_1.SetColLabelValue(5, "")
			self.grid_1.SetColLabelValue(6, "")
			self.grid_1.SetColLabelValue(7, "")
			self.grid_1.SetColLabelValue(8, "")
			self.grid_1.SetColLabelValue(9, "")
			self.grid_1.SetColLabelValue(10, "")
			self.grid_1.SetColLabelValue(11, "")
			
			self.grid_1.SetCellValue(0, 0, "# Optimisation cycles")
			self.grid_1.SetCellValue(0, 1, str(self.settings["numberOfGenerations"]))
			self.grid_1.SetCellValue(1, 0, "# Runs")
			self.grid_1.SetCellValue(1, 1, str(self.settings["numberOfPopulations"]))
			self.grid_1.SetCellValue(2, 0, "# Starting Structures")
			self.grid_1.SetCellValue(2, 1, str(self.settings["numberOfChromosomes"]))
			self.grid_1.SetCellValue(3, 0, "# Refinement cycles")
			self.grid_1.SetCellValue(3, 1, str(self.settings["numberOfRigidBodyCycles"]))
			self.grid_1.SetCellValue(4, 0, "Score clashes")
			self.grid_1.SetCellValue(4, 1, str(self.settings["scoreClashes"]))
			self.grid_1.SetCellValue(5, 0, "Symmetry")
			self.grid_1.SetCellValue(5, 1, str(self.settings["symmetry"]))
			#self.grid_1.SetCellValue(6, 0, "Score COGdiff")
			#self.grid_1.SetCellValue(6, 1, str(self.settings["scoreCOGdiff"]))
		
		if selection == 1:
			self.grid_1.SetColLabelValue(0, "Cycle")
			self.grid_1.SetColLabelValue(1, "rotX (°)")
			self.grid_1.SetColLabelValue(2, "rotY (°)")
			self.grid_1.SetColLabelValue(3, "rotZ (°)")
			self.grid_1.SetColLabelValue(4, u"tX (\u212B)")
			self.grid_1.SetColLabelValue(5, u"tY (\u212B)")
			self.grid_1.SetColLabelValue(6, u"tZ (\u212B)")
			self.grid_1.SetColLabelValue(7, u"R.m.s.d. (\u212B)")
			self.grid_1.SetColLabelValue(8, "Chi2")
			self.grid_1.SetColLabelValue(9, "Clashes")
			#self.grid_1.SetColLabelValue(10, u"COG (\u212B)")
			self.grid_1.SetColLabelValue(10, "Overall Score")
			
			log = self.dockingResults.populations[index].log
			log = log.split("\n")
			for idx, row in enumerate(log):
				row = row.split(" ")
				for idy, col in enumerate(row):
					self.grid_1.SetCellValue(idx, idy, col)
		if selection == 0:
			self.grid_1.SetColLabelValue(0, "Restraint")
			self.grid_1.SetColLabelValue(1, u"Docked distance(\u212B)")
			self.grid_1.SetColLabelValue(2, u"Experimental distance(Width)(\u212B)")
			self.grid_1.SetColLabelValue(3, u"Difference (\u212B)")
			self.grid_1.SetColLabelValue(4, "")
			self.grid_1.SetColLabelValue(5, "")
			self.grid_1.SetColLabelValue(6, "")
			self.grid_1.SetColLabelValue(7, "")
			self.grid_1.SetColLabelValue(8, "")
			self.grid_1.SetColLabelValue(9, "")
			self.grid_1.SetColLabelValue(10, "")
			#self.grid_1.SetColLabelValue(11, "")
			restraintIds = []
			for restraint in self.restraints:
				restraintIds.append(restraint["Name"])
			expDistances = self.dockingResults.expDistances
			expErrors = self.dockingResults.expErrors
			dockedDistances = self.dockingResults.populations[index].chromosomes[0].trialDistances
			for idx, values in enumerate (zip(dockedDistances, expDistances, expErrors, restraintIds)):
				self.grid_1.SetCellValue(idx, 0, "%s" %(values[3]))
				self.grid_1.SetCellValue(idx, 1, "%1.1f" %values[0])
				self.grid_1.SetCellValue(idx, 2, "%1.1f (%1.1f)" %(values[1], values[2]))
				self.grid_1.SetCellValue(idx, 3, "%1.1f" %(numpy.abs(values[0]-values[1])))
		self.grid_1.AutoSizeColumns()

	#-------------------------------------------------------------------------------------

##########################################################################################

class ResultsList(wx.Panel):

	"""Class that shows a list of docking results."""

	#-------------------------------------------------------------------------------------

	def __init__(self, parent, dockingResults):
		wx.Panel.__init__(self, parent)
		self.dockingResults = dockingResults
		self.resultsOlv = ObjectListView(self, wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER, useAlternateBackColors=False)
		self.resultsOlv.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onListBoxSelectionChanged)
		self.resultsOlv.Bind(OLVEvent.EVT_ITEM_CHECKED, self.checkboxClicked)
		#self.resultsOlv.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.onDoubleClick)
		self.prepareResultsForOlv()
		self.setOlvColumns()
		#self.resultsOlv.InstallCheckStateColumn(self.pymolColumn)
		self.resultsOlv.SetObjects(self.resultListEntries)
		self.bounding_box = wx.StaticBox(self, label='Solutions')
		self.__set_properties()
		self.__do_layout()
		self.allChecked = True
		pub.subscribe(self.onAddPymolColumn, "add.pymol")

	#-------------------------------------------------------------------------------------

	def __set_properties(self):
		pass

	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		sizer_1 = wx.StaticBoxSizer(self.bounding_box, wx.VERTICAL)
		sizer_1.Add(self.resultsOlv, 1, wx.EXPAND | wx.ALL, border=7)
		self.SetSizer(sizer_1)

	#-------------------------------------------------------------------------------------
	
	def onAddPymolColumn(self):
		self.setOlvColumnsWithPymol()
		self.resultsOlv.InstallCheckStateColumn(self.pymolColumn)

	#-------------------------------------------------------------------------------------
	
	def prepareResultsForOlv(self):
		self.resultListEntries = []
		for idx, population in enumerate(self.dockingResults.populations):
			#chromosomes[0] gets the best individual
			chromosome = population.chromosomes[0]
			resultListEntry = ResultListEntry(self.resultsOlv, idx, chromosome.name, chromosome.chi2, chromosome.rmsd, chromosome.clashes, True)
			self.resultListEntries.append(resultListEntry)

	#-------------------------------------------------------------------------------------

	def setOlvColumnsWithPymol(self, data=None):
		self.pymolColumn = ColumnDefn("PyMOL", fixedWidth=70, checkStateGetter="isChecked")
		self.resultsOlv.SetColumns([
			self.pymolColumn,
			ColumnDefn("Name", "left", 300, "name", isSpaceFilling=False),
			ColumnDefn("Chi2", "left", 100, "chi2", stringConverter="%.2f"),
			ColumnDefn(u"Rmsd (\u212B)", "left", 100, "rmsd", stringConverter="%.2f"),
			ColumnDefn("Clashes", "left", 100, "clashes", stringConverter="%i", isSpaceFilling=True),
			])

	#-------------------------------------------------------------------------------------
	
	def setOlvColumns(self, data=None):
		#self.pymolColumn = ColumnDefn("PyMOL", fixedWidth=70, checkStateGetter="isChecked")
		self.resultsOlv.SetColumns([
			ColumnDefn("Name", "left", 300, "name", isSpaceFilling=False),
			ColumnDefn("Chi2", "left", 100, "chi2", stringConverter="%.2f"),
			ColumnDefn(u"Rmsd (\u212B)", "left", 100, "rmsd", stringConverter="%.2f"),
			ColumnDefn("Clashes", "left", 100, "clashes", stringConverter="%i", isSpaceFilling=True),
			])
			

	#-------------------------------------------------------------------------------------
		
	def checkboxClicked(self, data):
		name = self.resultsOlv.GetSelectedObject().name
		if self.resultsOlv.GetSelectedObject().isChecked:
			#print "enable: %s" %name
			cmd.enable(name)
		else:
			#print "disable: %s" %name
			cmd.disable(name)

	#-------------------------------------------------------------------------------------

	def onListBoxSelectionChanged(self, event):
		name = self.resultsOlv.GetSelectedObject().name
		pub.sendMessage("update.selection", selectedResult = self.resultsOlv.GetSelectedObject())	

	#-------------------------------------------------------------------------------------

##########################################################################################

class ResultsPageContents(wx.Panel):

	"""Class that shows the resuls of a docking run"""

	#-------------------------------------------------------------------------------------

	def __init__(self, parent, id, dockingResults, settings, restraints):
		wx.Panel.__init__(self, parent)
		self.settings = settings
		self.id = id
		self.dockingResults = dockingResults
		self.splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE)
		self.topPanel = ResultsList(self.splitter, dockingResults)
		self.bottomPanel = ResultsTable(self.splitter, dockingResults, settings, restraints)
		self.splitter.SplitVertically(self.topPanel, self.bottomPanel)
		self.__set_properties()
		self.__do_layout()

	#-------------------------------------------------------------------------------------

	def __set_properties(self):
		self.splitter.SetMinimumPaneSize(1)
		self.splitter.SetSashGravity(0.5)

	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		sizer_1 = wx.BoxSizer(wx.VERTICAL)
		sizer_1.Add(self.splitter, 1, flag = wx.EXPAND | wx.ALL, border=7)
		self.SetSizer(sizer_1)
		self.Layout()

	#-------------------------------------------------------------------------------------

##########################################################################################

class RestraintView(wx.Panel):

	"""Class that shows a restraint view"""

	#-------------------------------------------------------------------------------------

	def __repr__(self):
		return "id: %i; selected: %r" %(self.id, self.selected)

	#-------------------------------------------------------------------------------------

	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
		self.Bind(wx.EVT_LEFT_UP, self.onClick)
		self.SetBackgroundColour('#dddddd')
		self.active = True
		self.selected = False
		self.id = id
		self.label_0 = wx.TextCtrl(self, wx.ID_ANY, 'Restraint: %i' %self.id)
		self.activeCheckBox = wx.CheckBox(self, wx.ID_ANY, label = "Include in docking run?")
		self.activeCheckBox.Bind(wx.EVT_CHECKBOX, self.onCheckStateChanged)
		
		#Comboboxes
		currentPymolObjects = cmd.get_object_list('*_label or mD_anchor_*')
		self.label_1 = wx.StaticText(self, wx.ID_ANY, "Anchor point on fixed molecule:")
		self.combo_box_1 = wx.ComboBox(self, 0, style=wx.CB_READONLY, choices=currentPymolObjects)
		self.combo_box_1.Bind(wx.EVT_COMBOBOX, self.onSelectionChanged)
		
		self.label_2 = wx.StaticText(self, wx.ID_ANY, "Anchor point on moving molecule:")
		self.combo_box_2 = wx.ComboBox(self, 1, style=wx.CB_READONLY, choices=currentPymolObjects)
		self.combo_box_2.Bind(wx.EVT_COMBOBOX, self.onSelectionChanged)
		
		#TextCtrls
		self.label_3 = wx.StaticText(self, wx.ID_ANY, u"Distance (\u212B)")
		self.text_ctrl_1 = wx.TextCtrl(self, wx.ID_ANY, "0.0")
		self.label_4 = wx.StaticText(self, wx.ID_ANY, u"Distribution width (\u212B)")
		self.text_ctrl_2 = wx.TextCtrl(self, wx.ID_ANY, "3.0")
		self.__set_properties()
		self.__do_layout()

	#-------------------------------------------------------------------------------------

	def __set_properties(self):
		self.activeCheckBox.SetValue(True)
		self.label_0.SetMinSize((200, 25))
		self.label_1.SetMinSize((200, 25))
		self.label_2.SetMinSize((200, 25))
		self.label_3.SetMinSize((150, 25))
		self.label_4.SetMinSize((150, 25))

	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		sizer_0 = wx.BoxSizer(wx.HORIZONTAL)
		sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
		sizer_3 = wx.BoxSizer(wx.VERTICAL)
		sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
		sizer_5 = wx.BoxSizer(wx.VERTICAL)
		sizer_6 = wx.BoxSizer(wx.VERTICAL)
		sizer_7 = wx.BoxSizer(wx.VERTICAL)
		
		sizer_2.Add(sizer_3, 1, flag=wx.EXPAND | wx.ALL, border=3)
		sizer_2.Add(sizer_4, 0, flag=wx.EXPAND | wx.ALL, border=3)
		sizer_4.Add(sizer_5, 0, flag=wx.EXPAND | wx.TOP, border=23)
		sizer_4.Add(sizer_6, 0, flag=wx.EXPAND | wx.TOP, border=20)
		sizer_4.Add(sizer_7, 0, flag=wx.EXPAND | wx.ALL, border=3)
		
		sizer_3.Add(self.label_0, 0, 0, 0)
		sizer_3.Add(self.activeCheckBox, 0, 0, 0)
		sizer_3.Add(self.label_1, 0, flag=wx.EXPAND | wx.ALL, border=3)
		sizer_3.Add(self.combo_box_1, 1, wx.EXPAND, 0)
		sizer_3.Add(self.label_2, 0, flag=wx.EXPAND | wx.TOP, border=3)
		sizer_3.Add(self.combo_box_2, 1, wx.EXPAND, 0)
		sizer_5.AddSpacer((1,15))
		sizer_5.Add(self.label_3, 0, flag=wx.EXPAND | wx.LEFT, border=7)

		sizer_5.Add(self.text_ctrl_1, 0, flag=wx.EXPAND | wx.LEFT, border=7)
		sizer_5.AddSpacer((1,5))
		sizer_5.Add(self.label_4, 0, flag=wx.EXPAND | wx.LEFT, border=7)
		sizer_5.Add(self.text_ctrl_2, 0, flag=wx.EXPAND | wx.LEFT, border=7)
		sizer_0.Add(sizer_2, 1, flag=wx.EXPAND | wx.ALL, border=7)
		self.SetSizer(sizer_0)
		sizer_0.Fit(self)
		self.Layout()

	#-------------------------------------------------------------------------------------

	def getRestraintName(self):
		return self.label_0.GetValue()

	#-------------------------------------------------------------------------------------

	def onClick(self, event):
		if self.selected:
			self.selected = False
			self.BackgroundColour = '#dddddd'
			pub.sendMessage("restraint.deselected", index=self.id)
			self.Refresh()
		else:
			self.selected = True
			self.BackgroundColour = "#3ea100"
			pub.sendMessage("restraint.selected", index=self.id)
			self.Refresh()
		self.Layout()

	#-------------------------------------------------------------------------------------

	def onCheckStateChanged(self, event):
		self.active = self.activeCheckBox.GetValue()

	#-------------------------------------------------------------------------------------

	def onSelectionChanged(self, event):
		id = event.GetId()
		if id == 0:
			coord = cmd.get_model(self.combo_box_1.GetValue(), 1).get_coord_list()
			self.anchorAcoord = coord[0]
			print self.anchorAcoord
		if id == 1:
			coord = cmd.get_model(self.combo_box_2.GetValue(), 1).get_coord_list()
			self.anchorBcoord = coord[0]
			print self.anchorBcoord

	#-------------------------------------------------------------------------------------

	def updateComboBoxes(self):
		currentPymolObjects = cmd.get_object_list('(all)')
		self.combo_box_1.SetItems(currentPymolObjects)
		self.combo_box_2.SetItems(currentPymolObjects)

	#-------------------------------------------------------------------------------------

##########################################################################################

class DockingThread(Thread):
	"""Class that creates a new thread for a docking run."""
	
	#----------------------------------------------------------------------
	
	def __init__(self, runNumber, restraints, settings):
		"""Init Worker Thread Class."""
		Thread.__init__(self)
		self._stop = threading.Event()
		self.runNumber = runNumber
		self.restraints = restraints
		self.settings = settings
		self.start()
 
	#----------------------------------------------------------------------
	
	def run(self):
		"""Run Worker Thread."""
		while not self._stop.isSet():
			dockingRun = Docker(self.runNumber, self.restraints, self.settings)
			dockingResult, settings = dockingRun.dock()
			result = {}
			result["dockingResult"] = dockingResult
			result["settings"] = settings
			result["restraints"] = self.restraints
			wx.CallAfter(pub.sendMessage, "docking.ready", result=result)
			self.stop()
		
	#-------------------------------------------------------------------------------------

	def stop(self):
		wx.CallAfter(pub.sendMessage, "docking.abort")
		self._stop.set()

	#-------------------------------------------------------------------------------------

	def stopped(self):
		return self._stop.isSet()

	#-------------------------------------------------------------------------------------


##########################################################################################

class SettingsDialog(wx.Dialog):
	
	""" Class that shows a settings dialog."""
	
	#-------------------------------------------------------------------------------------

	def __init__(self, currentSettings):
		wx.Dialog.__init__(self, None, title="Advanced settings")
		self.text_ctrl_2 = wx.TextCtrl(self, wx.ID_ANY, "%i" %currentSettings["numberOfGenerations"])
		self.label_1 = wx.StaticText(self, wx.ID_ANY, "#Optimisation cycles")
		
		self.text_ctrl_1 = wx.TextCtrl(self, wx.ID_ANY, "%i" %currentSettings["numberOfChromosomes"])
		self.label_2 = wx.StaticText(self, wx.ID_ANY, "#Starting structures")
		
		self.text_ctrl_3 = wx.TextCtrl(self, wx.ID_ANY, "%i" %currentSettings["numberOfRigidBodyCycles"])
		self.label_3 = wx.StaticText(self, wx.ID_ANY, "#Refinement cycles")
		
		#self.combo_box_1 = wx.ComboBox(self, 0, style=wx.CB_READONLY, choices=["None", "C2", "C3", "C4", "C5", "C6", "C7", "C8"])
		#self.label_4 = wx.StaticText(self, wx.ID_ANY, "Symmetry")
		#self.combo_box_2 = wx.ComboBox(self, 0, style=wx.CB_READONLY, choices=["True", "False"])
		#self.label_5 = wx.StaticText(self, wx.ID_ANY, "Score clashes")
		#self.combo_box_3 = wx.ComboBox(self, 0, style=wx.CB_READONLY, choices=["True", "False"])
		#self.label_7 = wx.StaticText(self, wx.ID_ANY, "Score COGdiff")
		self.text_ctrl_4 = wx.TextCtrl(self, wx.ID_ANY, "%i" %currentSettings["numberOfPopulations"])
		self.label_6 = wx.StaticText(self, wx.ID_ANY, "#Runs")
		
		self.okBtn = wx.Button(self, wx.ID_OK)
		#self.__set_properties(currentSettings)
		self.__do_layout()

	#-------------------------------------------------------------------------------------

	def __set_properties(self, currentSettings):
		pass
		#self.combo_box_1.SetValue(currentSettings["symmetry"])
		#self.combo_box_2.SetValue("%s" %currentSettings["scoreClashes"])
		#self.combo_box_3.SetValue("%s" %currentSettings["scoreCOGdiff"])

	#-------------------------------------------------------------------------------------

	def __do_layout(self):
		sizer_1 = wx.BoxSizer(wx.VERTICAL)
		grid_sizer_1 = wx.GridSizer(6, 2, 0, 0)
		grid_sizer_1.Add(self.text_ctrl_2, 0, wx.ALL, border=7)
		grid_sizer_1.Add(self.label_1, 0,wx.ALL, border=7)
		grid_sizer_1.Add(self.text_ctrl_1, 0, wx.ALL, border=7)
		grid_sizer_1.Add(self.label_2, 0, wx.ALL, border=7)
		grid_sizer_1.Add(self.text_ctrl_3, 0, wx.ALL, border=7)
		grid_sizer_1.Add(self.label_3, 0, wx.ALL, border=7)
		#grid_sizer_1.Add(self.combo_box_1, 0,wx.ALL, border=7)
		#grid_sizer_1.Add(self.label_4, 0, wx.ALL, border=7)
		#grid_sizer_1.Add(self.combo_box_2, 0, wx.ALL, border=7)
		#grid_sizer_1.Add(self.label_5, 0, wx.ALL, border=7)
		#grid_sizer_1.Add(self.combo_box_3, 0, wx.ALL, border=7)
		#grid_sizer_1.Add(self.label_7, 0, wx.ALL, border=7)
		grid_sizer_1.Add(self.text_ctrl_4, 0,wx.ALL, border=7)
		grid_sizer_1.Add(self.label_6, 0, wx.ALL, border=7)
		sizer_1.Add(grid_sizer_1, 1, wx.ALL, border=7)
		sizer_1.Add(self.okBtn, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, border=7)
		self.SetSizer(sizer_1)
		sizer_1.Fit(self)
		self.Layout()

	#-------------------------------------------------------------------------------------

##########################################################################################

class MyGrid(wx.grid.Grid):
	""" A Copy enabled grid class
	based on:
	http://stackoverflow.com/questions/28509629/work-with-ctrl-c-and-ctrl-v-to-copy-and-paste-into-a-wx-grid-in-wxpython
	"""
	def __init__(self, parent, id, style):
		wx.grid.Grid.__init__(self, parent, id, wx.DefaultPosition, wx.DefaultSize, style)
		wx.EVT_KEY_DOWN(self, self.OnKey)

	def OnKey(self, event):
		# If Ctrl+C is pressed...
		if event.CmdDown() and event.GetKeyCode() == 67:
			self.copy()
		# Skip other Key events
		if event.GetKeyCode():
			event.Skip()
			return

	def copy(self):
		# Number of rows and cols
		#print self.GetSelectionBlockBottomRight()
		#print self.GetGridCursorRow()
		#print self.GetGridCursorCol()
		if self.GetSelectionBlockTopLeft() == []:
			rows = 1
			cols = 1
			iscell = True
		else:
			rows = self.GetSelectionBlockBottomRight()[0][0] - self.GetSelectionBlockTopLeft()[0][0] + 1
			cols = self.GetSelectionBlockBottomRight()[0][1] - self.GetSelectionBlockTopLeft()[0][1] + 1
			iscell = False
		# data variable contain text that must be set in the clipboard
		data = ''
		# For each cell in selected range append the cell value in the data variable
		# Tabs '\t' for cols and '\r' for rows
		for r in range(rows):
			for c in range(cols):
				if iscell:
					data += str(self.GetCellValue(self.GetGridCursorRow() + r, self.GetGridCursorCol() + c))
				else:
					data += str(self.GetCellValue(self.GetSelectionBlockTopLeft()[0][0] + r, self.GetSelectionBlockTopLeft()[0][1] + c))
				if c < cols - 1:
					data += '\t'
			data += '\n'
		# Create text data object
		clipboard = wx.TextDataObject()
		# Set data object value
		clipboard.SetText(data)
		# Put the data in the clipboard
		if wx.TheClipboard.Open():
			wx.TheClipboard.SetData(clipboard)
			wx.TheClipboard.Close()
		else:
			wx.MessageBox("Can't open the clipboard", "Error")
