from pymol.wizard import Wizard
from pymol import cmd, stored, util
import pymol
import math
import sys, os
import label
import multiprocessing
import numpy
import ensemble
import rotamer
import atom
import distanceDistribution
import distanceMap
import uuid

class MtsslWizard(Wizard):
	def __init__(self):
		"""main mtsslWizard class"""
		Wizard.__init__(self)
		self.printStartupInfo()
		self.setupMenus()
		self.reset()
		self.path = os.path.dirname(__file__)
	
	def reset(self):
		default_thoroughness = "thorough"
		default_label = "R1"
		default_cutoff = 3.4
		default_clashes = 0
		default_mode = "Search"
		default_vdwRestraints = "tight"
		internalClash_cutoff = 2.5
		self.homoOligomerMode = False
		self.running = False
		self.toggleStatesCaption="Toggle states: ON"
		self.object_prefix = "mW-"
		self.pickCount = 0
		self.object_count = 0
		self.vdwRestraints="tight"
		self.thoroughness = default_thoroughness
		self.cutoff = default_cutoff
		self.clashes = default_clashes
		self.set_currentLabel(default_label)
		self.residue1Name = None
		self.residue2Name = None
		self.pickedObject1 = None
		self.pickedObject2 = None
		self.numberOfLabel = 0
		self.mode = default_mode
		self.writeToFile = False
		self.label = None
		#create array for plots
		try:
			print stored.plots
		except:
			print "Created plots list."
			stored.plots = []
		cmd.set("mouse_selection_mode",1) # set selection mode to residue
		cmd.deselect()
		cmd.unpick()
		cmd.delete("sele*")
		cmd.delete("_indicate*")
		cmd.delete("pk*")
		cmd.delete("*_tmp*")
		cmd.refresh_wizard()	

	def set_currentLabel(self, name):
		switcher = {
		"R1": label.Label.fromfile("labels/r1.txt"),
		"K1": label.Label.fromfile("labels/k1.txt"),
		"PROXYL": label.Label.fromfile("labels/prox.txt"),
		"DOTA": label.Label.fromfile("labels/dta.txt"),
		"URIP": label.Label.fromfile("labels/urip.txt"),
		"C": label.Label.fromfile("labels/c.txt"),
		"TPA": label.Label.fromfile("labels/tpa.txt"),
		"TRITYL": label.Label.fromfile("labels/trt.txt"),
		"TRITYL-1": label.Label.fromfile("labels/trt-1.txt"),
		"TRITYL-2": label.Label.fromfile("labels/trt-2.txt"),
		"TRITYL-3": label.Label.fromfile("labels/trt-3.txt"),
		"TRITYL-4": label.Label.fromfile("labels/trt-4.txt"),
		"TRITYL-5": label.Label.fromfile("labels/trt-5.txt"),
		"dU": label.Label.fromfile("labels/dU.txt"),
		"EGFR": label.Label.fromfile("labels/egfr.txt"),
		"Rx": label.Label.fromfile("labels/rx.txt"),
		"TTRITYL": label.Label.fromfile("labels/trtdna.txt")
		}
		selectedLabel = switcher.get(name)
		#self.set_vdwRestraints(options[currentLabel]().defaultVdw)
		print selectedLabel.info
		self.currentLabel = selectedLabel
		self.set_vdwRestraints(selectedLabel.defaultVdw)
		self.cmd.refresh_wizard()

	def setupMenus(self):
		self.menu["currentLabel"] = [
						[2, "\\559Protein", ""],
						[1, "R1","cmd.get_wizard().set_currentLabel('R1')"],
						[1, "PROXYL","cmd.get_wizard().set_currentLabel('PROXYL')"],
						[1, "DOTA","cmd.get_wizard().set_currentLabel('DOTA')"],
						[1, "TRITYL","cmd.get_wizard().set_currentLabel('TRITYL')"],
						[1, "TRITYL-1","cmd.get_wizard().set_currentLabel('TRITYL-1')"],
						[1, "TRITYL-2","cmd.get_wizard().set_currentLabel('TRITYL-2')"],
						[1, "TRITYL-3","cmd.get_wizard().set_currentLabel('TRITYL-3')"],
						[1, "TRITYL-4","cmd.get_wizard().set_currentLabel('TRITYL-4')"],
						[1, "TRITYL-5","cmd.get_wizard().set_currentLabel('TRITYL-5')"],
						[1, "K1","cmd.get_wizard().set_currentLabel('K1')"],
						[1, "Rx","cmd.get_wizard().set_currentLabel('Rx')"],
						[1, "EGFR","cmd.get_wizard().set_currentLabel('EGFR')"],
						[2, "\\559DNA", ""],
						[1, "URIP","cmd.get_wizard().set_currentLabel('URIP')"],
						[1, "C","cmd.get_wizard().set_currentLabel('C')"],
						[1, "TPA","cmd.get_wizard().set_currentLabel('TPA')"],
						[1, "dU","cmd.get_wizard().set_currentLabel('dU')"],
						[1, "TTRITYL","cmd.get_wizard().set_currentLabel('TTRITYL')"]
						]

		self.menu["mode"] = [
							[1, "Search","cmd.get_wizard().set_mode('Search')"],
							[1, "Measure","cmd.get_wizard().set_mode('Measure')"],
							[1, "Distance Map","cmd.get_wizard().set_mode('Distance Map')"],
							]

		self.menu["currentLabel1"] = [
						[ 2, "\\559Protein", ""],
						[1, "R1","cmd.get_wizard().set_currentLabel('R1')"],
						[1, "PROXYL","cmd.get_wizard().set_currentLabel('PROXYL')"],
						]
						
		self.menu["thoroughness"] = [
						[1, "painstaking","cmd.get_wizard().set_thoroughness('painstaking')"],
						[1, "thorough search","cmd.get_wizard().set_thoroughness('thorough')"],
						[1, "quick search","cmd.get_wizard().set_thoroughness('quick')"],
						]
						
		self.menu["vdwRestraints"] = [
						[ 2, "\\559NOTE:\\559 ", ""],
						[ 2, "\\559'loose': vdW cutoff 2.5 A, 5 clashes allowed", ""],
						[ 2, "\\559'tight': vdW cutoff 3.4 A, 0 clashes allowed", ""],
						[1, "loose","cmd.get_wizard().set_vdwRestraints('loose')"],
						[1, "tight","cmd.get_wizard().set_vdwRestraints('tight')"],
						]
	
		self.menu["writeToFile"] = [
						[1, "Yes","cmd.get_wizard().set_writeToFile(True)"],
						[1, "No","cmd.get_wizard().set_writeToFile(False)"]
						]

		self.menu["homoOligomerMode"] = [
						[ 2, "\\955NOTE:\\559Use this to speed up the ", ""],
						[ 2, "\\559calculation for homooligomers. ", ""],
						[ 2, "\\559Avoids calculation of redundant maps. ", ""],
						[1, "Yes","cmd.get_wizard().set_homoOligomerMode(True)"],
						[1, "No","cmd.get_wizard().set_homoOligomerMode(False)"]
						]

	def printStartupInfo(self):
		print ""
		print "******************************************************************************"
		print "* MtsslWizard by gha.                                                        *"
		print "* Version 2.0                                                                *"
		print "* March 2016                                                                 *"
		print "******************************************************************************"
		print ""
	
	def get_panel(self):
		if self.mode == "Search":
			return [
					[ 1, "Mtssl Wizard",""],
					[ 3, "Mode: %s"%self.mode,"mode"],
					[ 3, "Label: %s"%self.currentLabel.name,"currentLabel"],
					[ 3, "Speed: %s"%self.thoroughness,"thoroughness"],
					[ 3, "vdW restraints: %s"%self.vdwRestraints,"vdwRestraints"],
					[ 2, "Search conformers!","cmd.get_wizard().run()"],
					[ 2, self.toggleStatesCaption,"cmd.get_wizard().toggle_states()"],
					[ 2, "Delete last label","cmd.get_wizard().delete_last()"],
					[ 2, "Remove solvent","cmd.get_wizard().remove_solvent()"],
					[ 2, "Reset","cmd.get_wizard().reset()"],
					[ 2, "Done","cmd.set_wizard()"],
					]
		elif self.mode == "Measure":
			return [
					[ 1, "Mtssl Wizard",""],
					[ 3, "Mode: %s"%self.mode,"mode"],
					[ 3, "Results to file?: %s"%self.writeToFile,"writeToFile"],
					[ 2, "Reset","cmd.get_wizard().reset()"],
					[ 2, "Done","cmd.set_wizard()"]
					]
		elif self.mode == "Distance Map":
			return [
					[ 1, "MtsslWizard", ""],
					[ 3, "Mode: %s"%self.mode,"mode"],
					[ 3, "Label: %s"%self.currentLabel.name,"currentLabel1"],
					[ 3, "Homooligomer mode?: %s"%self.homoOligomerMode,"homoOligomerMode"],
					[ 3, "Results to file?: %s"%self.writeToFile,"writeToFile"],
					[ 2, "Reset","cmd.get_wizard().reset()"],
					[ 2, "Done","cmd.set_wizard()"]
			]
	def set_thoroughness(self,thoroughness):
		self.thoroughness = thoroughness
		self.cmd.refresh_wizard()

	def remove_solvent(self):
		cmd.remove("solvent")
		cmd.refresh_wizard()
		
	def set_writeToFile(self,writeToFile):
		self.writeToFile = writeToFile
		self.cmd.refresh_wizard()
		
	def set_homoOligomerMode(self,homoOligomerMode):
		self.homoOligomerMode = homoOligomerMode
		self.cmd.refresh_wizard()
		
	def set_cutoff(self,cutoff):
		self.cutoff = cutoff
		self.cmd.refresh_wizard()
		
	def set_clashes(self,clashes):
		self.clashes = clashes
		self.cmd.refresh_wizard()
	
	def set_mode(self, mode):
		self.mode = mode
		#reset pickCount and deselect everything
		self.pickCount = 0
		cmd.deselect()
		cmd.delete("sele*")
		self.cmd.refresh_wizard()
	
	def set_vdwRestraints(self, vdwRestraints):
		self.vdwRestraints = vdwRestraints
		if vdwRestraints == "loose":
			self.set_cutoff(2.5)
			self.set_clashes(5)
		elif vdwRestraints == "tight":
			self.set_cutoff(3.4)
			self.set_clashes(0)
		self.cmd.refresh_wizard()
	
	def delete_all(self):
		print "Deleting everything..."
		cmd.delete("%s*" %self.object_prefix)
 
	def delete_last(self):
		try:
			print self.numberOfLabel
			if self.numberOfLabel >= 1:
				cmd.delete("%s%s*" %(self.object_prefix, str(self.numberOfLabel)))
				self.numberOfLabel -= 1
		except pymol.CmdException, pmce:
			print pmce

	def cleanupAfterRun(self):
		cmd.deselect()
		self.pickCount = 0
		cmd.delete("pk*")
		cmd.delete("sele*")
		cmd.delete("*_tmp*")
		cmd.delete("*tmp*")
		cmd.delete("_indicate*")
		cmd.delete("labelEnvironment*")
		cmd.delete("currentLabel")
		self.cmd.refresh_wizard()
	
	def toggle_states(self):
		if cmd.get("all_states")=='on':
			self.toggleStatesCaption='Toggle states: OFF'
			cmd.set("all_states",0)
		elif cmd.get("all_states")=='off':
			self.toggleStatesCaption='Toggle states: ON'
			cmd.set("all_states",1)
		self.cmd.refresh_wizard()

	def get_prompt(self):
		if self.mode == 'Search':
			if self.pickCount == 0 and self.currentLabel.uid != 'Rx':
				self.prompt = [ 'Select a residue to label...']
			elif  self.pickCount == 0 and self.currentLabel.uid == 'Rx':
				self.prompt = [ 'Select first anchor point...']
			elif  self.pickCount == 1 and self.currentLabel.uid == 'Rx':
				self.prompt = [ 'Select second anchor point...']
		if self.pickCount == 0 and self.mode == 'Measure':
			self.prompt = [ 'Select first label...']
		if self.pickCount == 1 and self.mode == 'Measure':
			self.prompt = [ 'Select second label...']
		if self.pickCount == 0 and self.mode == 'Distance Map':
			self.prompt = [ 'Select first object...']
		if self.pickCount == 1 and self.mode == 'Distance Map':
			self.prompt = [ 'Select second object...']
		return self.prompt

	def createSelectionMacro(self, selection):
		obj_list = cmd.get_names('objects')
		selectedObject=""
		for object in obj_list:
			if cmd.get_type(object)=="object:molecule":
				if cmd.count_atoms("(%s and (sele))" %(object)):
					selectedObject = object
					break
		my_dict = {"my_list" : []}
		cmd.iterate(selection, "my_list.append((segi,chain,resi,resn))", space=my_dict)
		my_list = my_dict['my_list']
		macro = "%s-%s-%s-%s-%s" %(selectedObject, my_list[0][0], my_list[0][1], my_list[0][2] ,my_list[0][3])
		return macro

	def do_select(self, name):
		try:
			self.do_pick(0)
		except pymol.CmdException, pmce:
			print pmce

	def do_pick(self, picked_bond):
		if self.mode == "Search" or self.mode == "Measure" or self.mode == "Distance Map":
			#Catch repeated selections in "Search" mode
			if self.pickCount > 0 and (self.mode == "Search" and not self.currentLabel.uid == 'Rx'):
				self.pickCount = 0
				cmd.delete("pk*")
				cmd.delete("sele*")
				cmd.delete("*_tmp*")
				cmd.delete("*tmp*")
				cmd.delete("_indicate*")
				cmd.delete("labelEnvironment*")
				cmd.deselect()
				#refresh wizard, so that the prompt is updated
				self.cmd.refresh_wizard()
				return
			#first click
			if self.pickCount == 0:
				self.residue1Name, self.pickedObject1 = self.click()
				if self.mode == "Measure" or self.mode == "Distance Map" or (self.mode == "Search" and self.currentLabel.uid == "Rx"):
					cmd.delete("sele*")
					cmd.deselect()
				self.pickCount += 1

			#second click
			elif self.pickCount == 1 and (self.mode == "Measure" \
											or self.mode == "Distance Map" \
											or (self.mode == "Search" and self.currentLabel.uid == "Rx")):
				self.residue2Name, self.pickedObject2 = self.click() 
				cmd.deselect()
				self.run()
				self.pickCount = 0
				cmd.deselect()
		#refresh wizard, so that the prompt is updated
		self.cmd.refresh_wizard()

	def click(self):
		residueName = self.createSelectionMacro("(sele)") 
		# transfer the click selection to a named selection
		cmd.select(residueName, "(sele)")
		# find the name of the object which contains the selection
		new_name = None
		obj_list = cmd.get_names("objects")
		for object in obj_list:
			if cmd.get_type(object)=="object:molecule":
				if cmd.count_atoms("(%s and %s)"%(object, residueName)):
					pickedObject = object
					break
		if pickedObject == None:
			print "MtsslWizard: object not found."
 		self.cmd.refresh_wizard()
		return residueName, pickedObject

	def run(self):
		my_view = cmd.get_view()
		#mark Search
		if self.mode == 'Search':
			#show message
			cmd.wizard("message", "Searching conformers...")
			cmd.refresh()
			if self.currentLabel.uid != "Rx":
				self.search()
				print "Creating Rotamers in PyMOL...",
				for aRotamer in self.currentLabel.ensembles["mW"].rotamers:
					self.createRotamerInPymol(aRotamer, "mW")
# 				if self.thoroughness == "painstaking" and self.currentLabel.uid == "R1":
# 					scoringWeights = {"totalContactWeight": 0.0, "typeOfContactWeight": 1.0}
# 					for aRotamer in self.currentLabel.ensembles["mW"].rotamers:
# 						aRotamer.score(scoringWeights)
# 					self.currentLabel.ensembles["mW"].sortRotamers("chi2")
# 					numberOfRotamers = len(self.currentLabel.ensembles["mW"].rotamers)
# 					newEnsemble = ensemble.Ensemble()
# 					newEnsemble.name = "contactFit"
# 					newEnsemble.rotamers = self.currentLabel.ensembles["mW"].rotamers[0:20]
# 					self.currentLabel.ensembles["contactFit"] = newEnsemble
# 					for aRotamer in self.currentLabel.ensembles["contactFit"].rotamers:
# 						self.createRotamerInPymol(aRotamer, "contactFit")
				print "done!"

			elif self.currentLabel.uid == "Rx":
				ca1 = numpy.array(cmd.get_model(self.residue1Name + " & name CA", 1).get_coord_list()[0])
				ca2 = numpy.array(cmd.get_model(self.residue2Name + " & name CA", 1).get_coord_list()[0])
				try:
					cb1 = numpy.array(cmd.get_model(self.residue1Name + " & name CB", 1).get_coord_list()[0])
				except:
					cb1 = ca1
				try:
					cb2 = numpy.array(cmd.get_model(self.residue2Name + " & name CB", 1).get_coord_list()[0])
				except:
					cb2 = ca2
					
				environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s or %s) and not (%s or %s)" \
									%(self.pickedObject1, \
									 self.residue1Name, \
									 self.residue2Name, \
									 self.residue1Name, \
									 self.residue2Name), 1).get_coord_list())
				anchor1rotamers = self.currentLabel.calculateCone(ca1, cb1, environmentatoms, numberOfAtoms=8000)
				anchor2rotamers = self.currentLabel.calculateCone(ca2, cb2, environmentatoms, numberOfAtoms=8000)
				solutions1 = []
				solutions2 = []
				for anchor1rotamer in anchor1rotamers:
					anAtom = anchor1rotamer.atoms["N1"]
					solutions1.append(anAtom.coordinate)
				for anchor2rotamer in anchor2rotamers:
					anAtom = anchor2rotamer.atoms["N1"]
					solutions2.append(anAtom.coordinate)
				
				#determine common accessible volume
				distances1 = self.currentLabel.quick_map(solutions1, cb2)
				distances2 = self.currentLabel.quick_map(solutions2, cb1)
				indices1 = numpy.where(numpy.any(distances1 > 6, axis=1))
				indices2 = numpy.where(numpy.any(distances2 > 6, axis=1))
				solutions1 = numpy.delete(solutions1, indices1, 0)
				solutions2 = numpy.delete(solutions2, indices2, 0)
				solutions = numpy.concatenate((solutions1, solutions2))
				
				#create resulting ensemble				
				newEnsemble = ensemble.Ensemble()
				newEnsemble.name = "mW"
				id = 0
				newRotamers = []
				if len(solutions) > 0:
					for solution in solutions:
						newRotamer = rotamer.Rotamer()
						thisAtom = atom.Atom()
						thisAtom.coordinate = solution
						thisAtom.name = "N1"
						thisAtom.element = "N"
						newRotamer.id = id
						newRotamer.atoms["N1"] = thisAtom
						id += 1
						newRotamers.append(newRotamer)
				else:
					print "Did not find any possible N1 locations. Are the two anchorpoints too far apart?"
				newEnsemble.rotamers = newRotamers
				self.currentLabel.ensembles[newEnsemble.name] = newEnsemble
				cmd.load("%s/labels/%s" %(self.path, self.currentLabel.pdbFile), "currentLabel")
				for aRotamer in self.currentLabel.ensembles["mW"].rotamers:
					self.createRotamerInPymol(aRotamer, "mW")
			self.numberOfLabel += 1
			self.finalCosmetics()
			#dismiss message
			cmd.wizard()
		
		#mark Measure
		elif self.mode == "Measure":
			cmd.wizard("message", "Calculating distances...")
			cmd.refresh()
			print "\n\n\nDistance calculation:\n"
			print "The dashed lines are the c-beta distance (green),\nand the distance between the geometric averages\nof the two ensembles (yellow).\n"
			print "The following statistics refer to the distribution\nof the individual distances between all conformers (may take a while):\n"
			
			#find out what the selections are
			stored.label1 = []
			stored.label2 = []
			stored.label1Coordinates = []
			stored.label2Coordinates = []
			stored.atomNames1 = []
			stored.atomNames2 = []
			
			#extract label info
			cmd.iterate(self.residue1Name, 'stored.label1.append(segi)')
			cmd.iterate(self.residue2Name, 'stored.label2.append(segi)')
			cmd.iterate(self.residue1Name, 'stored.atomNames1.append(name)')
			cmd.iterate(self.residue2Name, 'stored.atomNames2.append(name)')
			try:
				label1 = label.Label.fromfile("labels/%s.txt" %stored.label1[0])
				cmd.iterate_state(0, "%s & name %s" %(self.residue1Name, label1.spinLocation), 'stored.label1Coordinates.append((x,y,z))')
			except:
				cmd.iterate_state(0, "%s" %(self.residue1Name), 'stored.label1Coordinates.append((x,y,z))')
			try:	
				label2 = label.Label.fromfile("labels/%s.txt" %stored.label2[0])
				cmd.iterate_state(0, "%s & name %s" %(self.residue2Name, label2.spinLocation), 'stored.label2Coordinates.append((x,y,z))')
			except:
				cmd.iterate_state(0, "%s" %(self.residue2Name), 'stored.label2Coordinates.append((x,y,z))')
			#calculate distances
			distances = distanceDistribution.DistanceDistribution() 
			dist = distances.calculateDistanceDistribution(stored.label1Coordinates, stored.label2Coordinates)
			
			#create pseudoatom at average coordinate of each ensemble and display the distance between them
			atoms1=numpy.array(stored.label1Coordinates)
			atoms2=numpy.array(stored.label2Coordinates)
			avgAtoms1=numpy.average(atoms1, axis=0)
			avgAtoms2=numpy.average(atoms2, axis=0)
			self.createPseudoatom (avgAtoms1, "tmp_average1", 1)
			self.createPseudoatom (avgAtoms2, "tmp_average2", 1)
			cmd.distance(self.object_prefix+"avg","tmp_average1 & name PS1","tmp_average2 & name PS1")
			cmd.delete("tmp_average1")
			cmd.delete("tmp_average2")
			
			#cbeta distance if cbeta is present in both selections
			#cBetaDistance = 0.0
			if any("CB" in atom for atom in stored.atomNames1) and any("CB" in atom for atom in stored.atomNames2):
				cmd.distance(self.object_prefix+"cBeta", self.residue1Name+" & name CB",self.residue2Name+" & name CB")
				#for some reason, cmd.distance does not return the correct distance. Although it is shown in the viewer...
				#get_distance gives the correct distance, but does not create the object in the viewer.
				cBetaDistance = cmd.get_distance(self.residue1Name+" & name CB",self.residue2Name+" & name CB")
				cmd.set("dash_color", "green", self.object_prefix+"cBeta")

			
			histogram=numpy.histogram(dist, numpy.arange(100))
			envelopePlot = numpy.zeros((100,2))
			envelopePlot[0:99] = numpy.column_stack((histogram[1][0:len(histogram[1])-1], histogram[0]))
			
			#put point in mid of bin
			envelopePlot[:,0] += 0.5 
			normEnvelopePlot = numpy.copy(envelopePlot)
			normEnvelopePlot[:,1] = normEnvelopePlot[:,1]/numpy.amax(histogram[0])
			
			#combine dist and histogram to single array before output
			output = numpy.column_stack((envelopePlot, normEnvelopePlot[:,1]))
			averageDistance = numpy.average(dist)
			
			#make graph dictionary for mtsslPlotter
			graphtitle = "%s-%s" %(self.residue1Name, self.residue2Name)
			xlim = [0, 100]
			ylim = [0, 1]
			plotDictionary = self.makeGraphDataDictionary (graphtitle, "DistanceDistribution", "Distance (Angstrom)", "Relative Probability", output[:,0], output[:,2], 0, xlim, ylim)
			stored.plots.append(plotDictionary)
			print "Distribution plot added to memory. Inspect it with mtsslPlotter."
			
			#Copy to clipboard
			header = "Dist.   Count   Norm.Count\n"
			outputStr = header + numpy.array_str(output)
			outputStr = outputStr.replace("[", "")
			outputStr = outputStr.replace("]", "")
			self.copyStringToClipboard(outputStr)
			
			#Write to file
			if self.writeToFile:
				try:
					filename = "%s-%s" %(self.residue1Name, self.residue2Name)
					numpy.savetxt(filename+".txt", output, delimiter='\t')
					print "Written to file:"
					print "%s/%s" %(os.getcwd(), filename)
				except:
					print "Writing to file failed!"
			print self.calculateStatistics2(dist)
			try:
				if cBetaDistance > 0.0:
					print "Cbeta distance: %3.1f" %cBetaDistance
			except:
				print "No Cbeta distance."
			cmd.wizard()	
			
		#mark Distance Map	
		elif self.mode == "Distance Map":
			#show message
			cmd.wizard("message", "Calculating distance maps. Please be patient...")
			cmd.refresh()
			dm = distanceMap.DistanceMap(self.writeToFile, self.currentLabel, self.homoOligomerMode)
			if self.pickedObject1 == self.pickedObject2:
				dm.intraDistanceMap(self.pickedObject1)
			else:
				dm.interDistanceMap(self.pickedObject1, self.pickedObject2)
			print "Done!"
			cmd.wizard()
			
		
		self.cleanupAfterRun()
		cmd.set_view(my_view)
	
	def makeGraphDataDictionary (self, title, type, xTitle, yTitle, xData, yData, zData, xlim, ylim):
		graphData = {}
		graphData['Title'] = title
		graphData['Type'] = type
		graphData['xTitle'] = xTitle
		graphData['yTitle'] = yTitle
		graphData['xData'] = xData
		graphData['yData'] = yData
		graphData['zData'] = zData
		graphData["xlim"] = xlim
		graphData["ylim"] = ylim
		return graphData	
	
	def calculateStatistics2(self, distances):
		statisticsResult=""
		average = numpy.average(distances)
		median = numpy.median(distances)
		stddev = numpy.std(distances)
		longest = numpy.amax(distances)
		shortest = numpy.amin(distances)
		statisticsResult+= "Average of distribution: %3.2f\n" %average
		statisticsResult+= "Median of distribution: %3.2f\n" %median
		statisticsResult+= "Std. dev. of distribution: %3.2f\n" %stddev
		statisticsResult+= "Shortest distance: %3.2f\n" % shortest
		statisticsResult+= "Longest distance: %3.2f" %longest
		return statisticsResult
		
	
	def search(self):
		print "\n\n\nNew run:\n"
		#load the label and superpose onto selected position
		cmd.load("%s/labels/%s" %(self.path, self.currentLabel.pdbFile), "currentLabel")
		print "Attempting superposition..."
		if not self.superpose():
			print "Superposition does not work."
			print "Possible reasons:"
			print "1) Glycine? Mutate to Ala first."
			print "2) Trying to attach DNA label to Protein or vice versa?"
			if len(self.currentLabel.errorMessage) > 0:
				print "3) %s" %self.currentLabel.errorMessage
			self.cleanupAfterRun(my_view)
			return
		else:
			print "Superposition worked!"
		#if self.currentLabel.rotate == False:
		#	return
		#prepare movingAtoms array of label, put into correct order...
		stored.movingAtoms = []
		for i in range (0, len(self.currentLabel.atomNames)):
			xyz = cmd.get_model("%s & name %s" %("currentLabel", self.currentLabel.atomNames[i] ), 1).get_coord_list()
			stored.movingAtoms.extend(xyz)
		self.currentLabel.movingAtoms = numpy.array(stored.movingAtoms)
		
		#create object with only the atoms around the label to speed everything up 
		protein ="%s &! %s within %f of %s" %(self.pickedObject1, \
											  self.residue1Name, \
											  self.currentLabel.radius, \
											  "currentLabel")
		cmd.create ("labelEnvironment", "byres %s" %protein)
		stored.environmentAtomCoordinates = []
		stored.environmentAtomNames = []
		stored.environmentAtomResidueNames = []
		cmd.iterate_state(1, protein, 'stored.environmentAtomCoordinates.append((x,y,z))')
		cmd.iterate(protein, 'stored.environmentAtomNames.append(name)')
		cmd.iterate(protein, 'stored.environmentAtomResidueNames.append(resn)')
		environmentAtomCoordinates = numpy.array(stored.environmentAtomCoordinates)
		environmentAtomNames = numpy.array(stored.environmentAtomNames)
		environmentAtomResidueNames = numpy.array(stored.environmentAtomResidueNames)
		environmentAtomInfo = [environmentAtomCoordinates, environmentAtomNames, environmentAtomResidueNames]
		numberOfCPUs = multiprocessing.cpu_count()
		numberOfTries = self.currentLabel.numberOfTries[self.thoroughness]
		numberOfRotamers = self.currentLabel.numberToFind[self.thoroughness]
		processes = []
		chunkSize = int(math.ceil(numberOfRotamers/float(numberOfCPUs)))
		chunks = []
		sum = 0
		i = 0
		while sum + chunkSize < numberOfRotamers:
			chunks.append(chunkSize)
			sum += chunkSize
			i += 1
		chunks.append(numberOfRotamers - i * chunkSize)
		
		numberOfProcesses = len(chunks)
		
		queue = multiprocessing.Queue()
		newEnsemble = ensemble.Ensemble()
		newEnsemble.name = "mW"
		self.currentLabel.ensembles[newEnsemble.name] = newEnsemble
		#only use multiprocessing on mac or linux
		if os.name != "nt":
			#print chunks
			print "Trying to find %s rotamers. Using %i cores." %(numberOfRotamers, numberOfProcesses),
			for i in range (numberOfProcesses):
				p = multiprocessing.Process(target = self.currentLabel.generateEnsembleMulti,
											args = (self.currentLabel.movingAtoms,\
											environmentAtomInfo, chunks[i], numberOfTries,\
											newEnsemble.name, \
											False, self.cutoff, self.clashes,\
											queue))
				p.start()
				processes.append(p)
		else:
			print "Trying to find %s rotamers. Using 1 core." %(numberOfRotamers),
			numberOfProcesses = 1
			self.currentLabel.generateEnsembleMulti(self.currentLabel.movingAtoms,\
											environmentAtomInfo, numberOfRotamers, numberOfTries,\
											newEnsemble.name, \
											False, self.cutoff, self.clashes,\
											queue)
		resultsDictionary = {}
		for i in range (numberOfProcesses):
			resultsDictionary.update(queue.get())
		for p in processes:
			p.join()
		print "Done."
		print "Collecting results...",
		newRotamers = []
		for resultList in resultsDictionary.values():
			for result in resultList:
				newRotamers.append(result)
		print "Done! Found %i rotamers" %len(newRotamers)
		#reassign ids. They are not unique with multiprocessing
		for idx, rotamer in enumerate(newRotamers):
			rotamer.id = idx
		newEnsemble.rotamers = newRotamers
		
	def createRotamerInPymol(self, rotamer, ensembleName):
		atomNames = self.currentLabel.atomNames
		for atomName in atomNames:
			thisAtom = rotamer.atoms[atomName]
			stored.xyz = []
			stored.xyz = thisAtom.coordinate
			try:
				cmd.alter_state(1,"%s & name %s " %("currentLabel", atomName), "(x,y,z)=stored.xyz")
			except:
				print "Could not alter coordinate"
			cmd.alter("%s & name %s " %("currentLabel", atomName), "segi='%s'" %self.currentLabel.uid)
		cmd.create("%s_%s_%s" %(self.residue1Name, self.currentLabel.identifier, ensembleName), "currentLabel", 1, rotamer.id+1)

	def finalCosmetics(self): #make everything look nice
		for ensemble in self.currentLabel.ensembles:
			if len(self.currentLabel.ensembles[ensemble].rotamers) > 0:
				#show all conformations and do some coloring
				cmd.set("all_states",1)
				self.toggleStatesCaption='Toggle states: ON'
				cmd.color("blue","%s_%s_%s" %(self.residue1Name, self.currentLabel.identifier, ensemble))
				cmd.color("red", "%s_%s_%s & name %s" %(self.residue1Name, self.currentLabel.identifier, ensemble, self.currentLabel.highlight))
				util.cnc("%s_%s_%s" %(self.residue1Name, self.currentLabel.identifier, ensemble))
				identifierLabel = "%s|%s|%s|%s|%s" %(self.residue1Name, self.currentLabel.identifier, self.vdwRestraints, self.thoroughness, ensemble)
				#pseudoatom at average N1 position
				stored.label = []
				cmd.iterate_state(0, "%s_%s_%s & name %s" %(self.residue1Name, self.currentLabel.identifier, ensemble, self.currentLabel.spinLocation), 'stored.label.append((x,y,z))')
				atoms1 = numpy.array(stored.label)
				#create pseudoatom at average coordinate of each ensemble
				avgAtoms = numpy.average(atoms1,axis = 0)
				self.createPseudoatom (avgAtoms, "%s_%s_label" %(self.residue1Name, ensemble), 1)
				cmd.set("sphere_scale", "0.5", "%s_%s_label" %(self.residue1Name, ensemble))
				cmd.label("%s_%s_label" %(self.residue1Name, ensemble), `identifierLabel`)
				cmd.show("label")
				cmd.show("spheres", "name PS1")
				cmd.group("%s%s" %(self.object_prefix, str(self.numberOfLabel)), "%s*, labelEnvironment_%s,%s*" %("currentLabel", "currentLabel", self.residue1Name))

	def createPseudoatom (self, coordinates, objectName, state):
		x=float(coordinates[0])
		y=float(coordinates[1])
		z=float(coordinates[2])
		posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
		cmd.pseudoatom(pos=posString, object=objectName, state=state)	
	
	def copyStringToClipboard(self, string):
		try:
			import pyperclip
			pyperclip.copy(string)
			print "Copied to clipboard."
			return
		except:
			pass
		try:
			import xerox
			xerox.copy(string)
			print "Copied to clipboard."
			return
		except:
			pass
		print "Copy to clipboard failed. Try to install either the 'pyperclip' or 'xerox' module for Python."
		
	def superpose(self):
		#get the position of the selected residue's O atom
		stored.xyz = []
		if self.currentLabel.modifiedAA:
			cmd.iterate_state(1,"%s & name O" %self.residue1Name, "stored.xyz.append([x,y,z])")
		args=[]
		i = 0
		while i < len(self.currentLabel.atomsForSuperposition):
			args.append("%s & name %s" %("currentLabel", self.currentLabel.atomsForSuperposition[i]))
			args.append("%s & name %s" %(self.residue1Name, self.currentLabel.atomsForSuperposition[i]))
			i += 1
		#print args
		if apply(cmd.pair_fit, args):
			#set the label's O atom to the stored position
			if self.currentLabel.modifiedAA:
				cmd.alter_state(1,"%s & name O" %"currentLabel", "(x,y,z)=stored.xyz.pop(0)")
			return True
		else:
			return False