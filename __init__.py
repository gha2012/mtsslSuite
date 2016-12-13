import pymol
import os
from pymol import cmd
import wx
from pymol import stored
import threading

try:
	from mtsslWizard import mtsslWizard
except Exception, e:
	print "Cannot import mtsslWizard: %s" %e
try:
	from mtsslPlotter import mtsslPlotter
except Exception, e:
	print "Cannot import mtsslPlotter: %s" %e
try:
	from mtsslTrilaterate import mtsslTrilaterate
except Exception, e:
	print "Cannot import mtsslTrilaterate: %s" %e

try:
	from mtsslDock import mtsslDock
except Exception, e:
	print "Cannot import mtsslDock: %s" %e

def __init__(self):
	#add menu
	self.menuBar.addmenu('mtsslSuite', 'mtsslSuite')
	#mtsslWizard
	self.menuBar.addmenuitem('mtsslSuite', 'command',
							 'MtsslWizard',
							 label = 'mtsslWizard',
							 command = lambda s=self : open_mtsslWizard())
	#mtsslPlotter
	self.menuBar.addmenuitem('mtsslSuite', 'command',
							 'mtsslPlotter',
							 label = 'mtsslPlotter',
							 command = lambda s=self : open_mtsslPlotter())
							 
	#mtsslTrilaterate
	self.menuBar.addmenuitem('mtsslSuite', 'command', 'mtsslTrilaterate',
							 label='mtsslTrilaterate',
							 command=lambda s=self: open_mtsslTrilaterate())
							 
	#mtsslDock
	self.menuBar.addmenuitem('mtsslSuite', 'command', 'mtsslDock',
							 label = 'mtsslDock',
							 command = lambda s=self : open_mtsslDock())

##########################################################################################
#start mtsslWizard																		 #
##########################################################################################
def open_mtsslWizard():
	#print mtsslWizard
	wiz = mtsslWizard.MtsslWizard()
	cmd.set_wizard(wiz)
	
##########################################################################################
#start mtsslPlotter																		 #
##########################################################################################
class mtsslPlotterApp(wx.App):
	def OnInit(self):
		frame = mtsslPlotter.MainWindow()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True

def run_mtsslPlotter():
	app = mtsslPlotterApp(0)
	app.MainLoop()

def open_mtsslPlotter():
	if hasattr(stored, 'mtsslplot'):
		if (stored.mtsslplot == 1):
			print "mtsslPlotter is already opened!"
		else:
			t = threading.Thread(target=run_mtsslPlotter,args=())
			t.setDaemon(0)
			t.start()
			stored.mtsslplot = 1
	else:
		t = threading.Thread(target=run_mtsslPlotter,args=())
		t.setDaemon(0)
		t.start()
		stored.mtsslplot = 1
	try:
		print "Available plots:"
		print len(stored.plots)
	except:
		print "Creating plots list."
		stored.plots = []
		
##########################################################################################
#start mtsslTrilaterate																	 #
##########################################################################################
class mtsslTrilaterateApp(wx.App):
	def OnInit(self):
		frame = mtsslTrilaterate.mtsslTrilaterate(None, -1, "mtsslTrilaterate")
		frame.Show(True)
		self.SetTopWindow(frame)
		return True

def run_mtsslTrilaterate():
	app = mtsslTrilaterateApp(0)
	app.MainLoop()


def open_mtsslTrilaterate():
	t = threading.Thread(target=run_mtsslTrilaterate, args=())
	t.setDaemon(1),
	t.start()
	
##########################################################################################
#start mtsslDock																		 #
##########################################################################################
class mtsslDockApp(wx.App):
	def OnInit(self):
		frame = mtsslDock.MainWindow()
		frame.Show(True)
		self.SetTopWindow(frame)
		return True

def run_mtsslDock():
	app = mtsslDockApp(0)
	app.MainLoop()

def open_mtsslDock():
	t = threading.Thread(target=run_mtsslDock,args=())
	t.setDaemon(0)
	t.start()
