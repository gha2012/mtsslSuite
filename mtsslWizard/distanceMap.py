import numpy
import numpy.ma as ma
import pymol
from pymol import cmd, stored
import scipy.spatial.distance
import random
import string
import uuid

class DistanceMap:

	def __init__(self, writeToFile, label, homoOligomer):
		self.writeToFile = writeToFile
		self.label = label
		self.homoOligomerMode = homoOligomer

	def makeGraphDataDictionary (self, id, title, type, xTitle, yTitle, xData, yData, zData, xlim, ylim):
		graphData = {}
		graphData["id"] = id
		graphData["Title"] = title
		graphData["Type"] = type
		graphData["xTitle"] = xTitle
		graphData["yTitle"] = yTitle
		graphData["xData"] = xData
		graphData["yData"] = yData
		graphData["zData"] = zData
		graphData["xlim"] = xlim
		graphData["ylim"] = ylim
		return graphData	

	def calculateCone(self, ca, cb, environmentAtoms, numberOfAtoms = 1000):
		# generate umbrella of trial atoms
		#print self.label
		p = 2*self.label.trialAtomSphereRadius * numpy.random.rand(3, numberOfAtoms)- self.label.trialAtomSphereRadius
		p = p[:, sum(p* p, 0)** .5 <= self.label.trialAtomSphereRadius]
		p = p.transpose()
		p = p + cb
		distances = self.quick_map(ca, p)
		indices = numpy.where(numpy.any(distances < self.label.exclusionSphereRadius, axis=1))
		p = numpy.delete(p, indices, 0)

		#compute distances between trial sphere and environment
		distances = self.quick_map(p, environmentAtoms)
		#check for clashes and discard clashing atoms from sphere
		indices = numpy.where(numpy.any(distances < 3.5, axis=1))
		solutions = numpy.delete(p, indices,0)
		#print solutions
		numberOfConeAtoms = numpy.shape(solutions)[0]

		#plotting the spheres in Pymol
		ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
		if numberOfConeAtoms < 1:
			return ca, numberOfConeAtoms, solutions
		else:
			#calculate the average cone coordinate
			averageConeCoordinate = numpy.average(solutions,axis=0)
			return averageConeCoordinate, numberOfConeAtoms, solutions

	def interDistanceMap(self, pymolObject1, pymolObject2):
		chains_object1 = cmd.get_chains(pymolObject1)
		chains_object2 = cmd.get_chains(pymolObject2)
		print chains_object1, chains_object2
		#iterate over chains in object 1
		for chain_1 in chains_object1:
			intraDistObj1, obj1x, obj1y, diagonalObj1, accObj1_1, accObj1_1, limits1 = self.calculateDistanceMap(pymolObject1, chain_1, pymolObject1, chain_1)
			cBetaDistObj1, obj1x, obj1y, diagonalObj1, accObj1_1, accObj1_1, limits2 = self.calculateDistanceMap(pymolObject1, chain_1, pymolObject1, chain_1, cBetaDM=True)
			
			#iterate over chains in object 2
			for chain_2 in chains_object2:
				id = uuid.uuid1()
				intraDistObj2, obj2x, obj2y, diagonalObj2, accObj2_2, accObj2_2, limits3 = self.calculateDistanceMap(pymolObject2, chain_2, pymolObject2, chain_2)
				cBetaDistObj2, obj2x, obj2y, diagonalObj2, accObj2_2, accObj2_2, limits4 = self.calculateDistanceMap(pymolObject2, chain_2, pymolObject2, chain_2, cBetaDM=True)
			
				#calculate inter chain distances for object 1-2
				interDistObj1_2, junk, junk, diagonal1_2, accObj1_1, accObj2_2, limits5 = self.calculateDistanceMap(pymolObject1, chain_1, pymolObject2, chain_2)
				interDistObj2_1, junk, junk, diagonal2_1, accObj2_2, accObj1_2, limits6 = self.calculateDistanceMap(pymolObject2, chain_2, pymolObject1, chain_1)
			
				#calculate difference distance matrix based on cbetas
				if intraDistObj1.shape == intraDistObj2.shape:
					xlim = [limits2["minChain_1"], limits2["maxChain_1"]]
					ylim = [limits4["minChain_2"], limits4["maxChain_2"]]
					differences = numpy.abs(intraDistObj1-intraDistObj2)
					x = numpy.linspace(1, differences.shape[0], num = differences.shape[0])
					y = numpy.linspace(1, differences.shape[1], num = differences.shape[1])
					cbetaDifferences = numpy.abs(cBetaDistObj1-cBetaDistObj2)
					z = numpy.transpose(cbetaDifferences)
					zm = ma.masked_where(numpy.isnan(z),z)
					fileName = "%s-%s_%s-%s_ddm_cbeta" %(pymolObject1, chain_1, pymolObject2, chain_2)
					self.plotMap(id, fileName, "DifferenceMap", "Number of Residue", "Number of Residue", x, y, zm, xlim, ylim)
					self.plotMap(id, "%s_accessibility" %pymolObject1, "AccessibilityPlot", "Number of Residue", "Acc.", x, accObj1_1, 0, xlim, 0)
					self.plotMap(id, "%s_accessibility" %pymolObject2, "AccessibilityPlot", "Number of Residue", "Acc.", y, accObj2_2, 0, ylim, 0)
				else:
					fileName = "%s-%s_%s-%s_ddm_cbeta" %(pymolObject1, chain_1, pymolObject2, chain_2)
					print "%s: Matrices have different shapes. Cannot calculate DDM." %(fileName)
					print intraDistObj1.shape, intraDistObj2.shape
					print "Please trim the molecules in PyMOL and try again."
					print intraDistObj1.shape, intraDistObj2.shape
			
				
				#plot inter distance map
				z = numpy.transpose(interDistObj1_2)
				x = numpy.linspace(1, z.shape[0], num = z.shape[0])
				y = numpy.linspace(1, z.shape[1], num = z.shape[1])#x
				xlim = [limits5["minChain_1"], limits5["maxChain_1"]]
				ylim = [limits5["minChain_2"], limits5["maxChain_2"]]
				zm = ma.masked_where(numpy.isnan(z),z)
				fileName = "%s_%s_%s-%s" %(pymolObject1, pymolObject2, chain_1, chain_2)
				self.plotMap(id, fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, zm, xlim, ylim)
				if cmd.count_atoms(pymolObject1) > 1 and cmd.count_atoms(pymolObject2) > 1:
					self.plotMap(id, fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, zm, xlim, ylim)
					self.plotMap(id, "%s_accessibility" %pymolObject1, "AccessibilityPlot", "Number of Residue", "Acc.", x, accObj1_1, 0, xlim, 0)
					self.plotMap(id, "%s_accessibility" %pymolObject2, "AccessibilityPlot", "Number of Residue", "Acc.", x, accObj2_2, 0, xlim, 0)
				elif cmd.count_atoms(pymolObject1) == 1 and cmd.count_atoms(pymolObject2) > 1:
					print xlim, ylim, z.shape
					self.plotMap(id, fileName, "DistancePlot", "Number of Residue", "Distance", numpy.linspace(1, z.shape[1], num = z.shape[1]), zm[:,0], 0, ylim, 100)
				elif cmd.count_atoms(pymolObject1) > 1 and cmd.count_atoms(pymolObject2) == 1:
					self.plotMap(id, fileName, "DistancePlot", "Number of Residue", "Distance", numpy.linspace(1, z.shape[1], num = z.shape[1]), zm[0,:], 0, xlim, 100)
				else:
					pass	

	def plotMap(self, id, fileName, type, xTitle, yTitle, x, y, z, xlim, ylim):
		plotDictionary = self.makeGraphDataDictionary(id, fileName, type, xTitle, yTitle, x, y, z, xlim, ylim)
		stored.plots.append(plotDictionary)
		if self.writeToFile:
			numpy.savetxt(fileName, z)			

	def calculateApproximateLabelCoordinate(self, residues):
		#iterate over residues in chain 1 and calculate approximate label position
		labelPositions = []
		cBetaPositions = []
		numberOfConeAtomsList = []
		for residue in residues:
			selectionString = "%s & chain %s & resi %s" %(residue[0], residue[1], residue[2])
			ca = numpy.array(cmd.get_model(selectionString + " & name CA", 1).get_coord_list()[0])
			try:
				cb = numpy.array(cmd.get_model(selectionString + " & name CB", 1).get_coord_list()[0])
			except:
				cb = ca
			environmentatoms = numpy.array(cmd.get_model("(%s within 10 of %s) &! (%s)" %(residue[0], selectionString, selectionString), 1).get_coord_list())
			averageConeCoordinate, numberOfConeAtoms, solutions = self.calculateCone(ca, cb, environmentatoms)
			labelPositions.append(averageConeCoordinate)
			cBetaPositions.append(cb)
			numberOfConeAtomsList.append(numberOfConeAtoms)
		return labelPositions, cBetaPositions, numberOfConeAtomsList

	def quick_map(self, atoms1, atoms2):
		# if there is only one atom it has to be duplicated for quick_dist2 to work
		duplicated = False
		if len(numpy.shape(atoms1)) == 1:
			duplicated = True
			atoms1 = numpy.tile(atoms1, (2,1))
		if len(numpy.shape(atoms2)) == 1:
			duplicated = True
			atoms2 = numpy.tile(atoms2, (2,1))
		dist = scipy.spatial.distance.cdist(atoms1, atoms2)
		#remove the duplication depending on which selection contained the single atom
		if duplicated and dist.shape[0] == 2 and not dist.shape[1] == 2:
			dist=numpy.reshape(dist[0,:], (-1, 1))

		elif duplicated and dist.shape[1] == 2 and not dist.shape[0] == 2:
			dist=numpy.reshape(dist[:,0], (-1, 1))

		elif duplicated and dist.shape[0] == 2 and dist.shape[1] == 2:
			dist=numpy.reshape(dist[:1,0], (-1, 1))
		return dist

	def calculateDistanceMap(self, pymolObject1, chain_1, pymolObject2, chain_2, cBetaDM=False):
		print "Working on: %s_%s-%s_%s" %(pymolObject1, chain_1, pymolObject2, chain_2)
		stored.residueList = []
		chain_1_residues = []
		chain_1_labelPositions = []
		chain_1_cBetaPositions = []
		chain_1_numberOfConeAtoms = []
		if cmd.count_atoms(pymolObject1) > 1:
			cmd.iterate("%s & chain %s & name CA & polymer" %(pymolObject1, chain_1), "stored.residueList.append((model, chain, resi,resn))")
			chain_1_residues = stored.residueList
			stored.residueList = []
			chain_1_labelPositions, chain_1_cBetaPositions, chain_1_numberOfConeAtoms = self.calculateApproximateLabelCoordinate(chain_1_residues)
		else:
			cmd.iterate(pymolObject1, "stored.residueList.append((model, chain, resi,resn))")
			chain_1_residues = stored.residueList
			chain_1_residues[0] =(pymolObject1, chain_1, "1","CA")
			chain_1_labelPositions.append(numpy.array(cmd.get_model(pymolObject1, 1).get_coord_list()[0]))
			chain_1_cBetaPositions.append(numpy.array(cmd.get_model(pymolObject1, 1).get_coord_list()[0]))
			chain_1_numberOfConeAtoms.append(0)
		stored.residueList = []
		chain_2_labelPositions = []
		chain_2_cBetaPositions = []
		chain_2_numberOfConeAtoms = []
		chain_2_residues = []
		if cmd.count_atoms(pymolObject2) > 1:
			cmd.iterate("%s & chain %s & name CA & polymer" %(pymolObject2, chain_2), "stored.residueList.append((model, chain, resi,resn))")
			chain_2_residues = stored.residueList
			chain_2_labelPositions, chain_2_cBetaPositions, chain_2_numberOfConeAtoms = self.calculateApproximateLabelCoordinate(chain_2_residues)
		else:
			cmd.iterate(pymolObject2, "stored.residueList.append((model, chain, resi,resn))")
			chain_2_residues = stored.residueList
			chain_2_residues[0] =(pymolObject2, chain_2, "1","CA")
			chain_2_labelPositions.append(numpy.array(cmd.get_model(pymolObject2, 1).get_coord_list()[0]))
			chain_2_cBetaPositions.append(numpy.array(cmd.get_model(pymolObject2, 1).get_coord_list()[0]))
			chain_2_numberOfConeAtoms.append(0)
		
		#print chain_1_residues, chain_2_residues
		#account for truncations, gaps, ...
		#find highest residue number
		minChain_1 = 100000
		minChain_2 = 100000
		maxChain_1 = 0
		maxChain_2 = 0
		for residue in chain_1_residues:
			if int(residue[2]) > maxChain_1:
				maxChain_1 = int(residue[2])
			if int(residue[2]) < minChain_1:
				minChain_1 = int(residue[2])
		for residue in chain_2_residues:
			if int(residue[2]) > maxChain_2:
				maxChain_2 = int(residue[2])
			if int(residue[2]) < minChain_2:
				minChain_2 = int(residue[2])
		if maxChain_1 >= maxChain_2:
			maxLength = maxChain_1
		else:
			maxLength = maxChain_2

		plotLimits = {"minChain_1": minChain_1, "maxChain_1": maxChain_1, "minChain_2": minChain_2, "maxChain_2": maxChain_2}
	
		#create nan arrays with maxLength rows
		new_chain_1_labelPositions = numpy.empty([maxLength, 3])
		new_chain_1_labelPositions[:] = numpy.NAN
		new_chain_2_labelPositions = numpy.empty([maxLength, 3])
		new_chain_2_labelPositions[:] = numpy.NAN
		new_chain_1_cBetaPositions = numpy.empty([maxLength, 3])
		new_chain_1_cBetaPositions[:] = numpy.NAN
		new_chain_2_cBetaPositions = numpy.empty([maxLength, 3])
		new_chain_2_cBetaPositions[:] = numpy.NAN
		new_chain_1_numberOfConeAtoms = numpy.empty([maxLength, 1])
		new_chain_1_numberOfConeAtoms[:] = numpy.NAN
		new_chain_2_numberOfConeAtoms = numpy.empty([maxLength, 1])
		new_chain_2_numberOfConeAtoms[:] = numpy.NAN

		#put the data in the nan arrays
		for idx, residue in enumerate(chain_1_residues):
			new_chain_1_labelPositions[int(residue[2])-1,:] = numpy.array(chain_1_labelPositions[idx])
			new_chain_1_cBetaPositions[int(residue[2])-1,:] = numpy.array(chain_1_cBetaPositions[idx])
			new_chain_1_numberOfConeAtoms[int(residue[2])-1,:] = numpy.array(chain_1_numberOfConeAtoms[idx])
		for idx, residue in enumerate(chain_2_residues):
			new_chain_2_labelPositions[int(residue[2])-1,:] = numpy.array(chain_2_labelPositions[idx])	
			new_chain_2_cBetaPositions[int(residue[2])-1,:] = numpy.array(chain_2_cBetaPositions[idx])
			new_chain_2_numberOfConeAtoms[int(residue[2])-1,:] = numpy.array(chain_2_numberOfConeAtoms[idx])
		
		#calculate the distance matrix
		if cBetaDM:
			distances = self.quick_map(new_chain_1_cBetaPositions, new_chain_2_cBetaPositions)
		else:
			distances = self.quick_map(new_chain_1_labelPositions, new_chain_2_labelPositions)
		x = numpy.linspace(1, len(new_chain_1_labelPositions)+1, num = len(new_chain_1_labelPositions)+1)
		y = numpy.linspace(1, len(new_chain_2_labelPositions)+1, num = len(new_chain_2_labelPositions)+1)
		diagonal = numpy.diagonal(distances)
		diagonal = ma.masked_where(numpy.isnan(diagonal), diagonal)
		distances =  ma.masked_where(numpy.isnan(distances), distances)
		new_chain_1_numberOfConeAtoms = ma.masked_where(numpy.isnan(new_chain_1_numberOfConeAtoms),new_chain_1_numberOfConeAtoms)
		new_chain_2_numberOfConeAtoms = ma.masked_where(numpy.isnan(new_chain_2_numberOfConeAtoms),new_chain_2_numberOfConeAtoms)
		return distances, x, y, diagonal, new_chain_1_numberOfConeAtoms, new_chain_2_numberOfConeAtoms, plotLimits

	def intraDistanceMap(self, pymolObject1):
		labelPositions = []
		stored.residueList = []
		chains = cmd.get_chains(pymolObject1)
		#for homoOligomer mode
		chains1 = chains
		#iterate over chains
		if self.homoOligomerMode == True:
			chains1 = chains[0]
		for chain_1 in chains1:
			for chain_2 in chains:
				id = uuid.uuid1()
				try:
					distances, x, y, diagonal, accObj1, accObj2, limits = self.calculateDistanceMap(pymolObject1, chain_1, pymolObject1, chain_2)
				except:
					print "Could not calculate distance map."
					break
				#plot distance map
				xlim = [limits["minChain_1"], limits["maxChain_1"]]
				ylim = [limits["minChain_2"], limits["maxChain_2"]]
				z = numpy.transpose(distances)
				zm = ma.masked_where(numpy.isnan(z),z)
				fileName = "%s_%s-%s_distanceMatrix" %(pymolObject1, chain_1, chain_2)
				plotDictionary = self.makeGraphDataDictionary(id, fileName, "DistanceMap", "Number of Residue", "Number of Residue", x, y, zm, xlim, ylim)
				stored.plots.append(plotDictionary)
				if self.writeToFile:
					numpy.savetxt(fileName, distances)
		
				#Plot accessibility for chain 1
				fileName = "%s_%s_accessibilty" %(pymolObject1, chain_1)
				x = numpy.linspace(1, accObj1.shape[0], num = accObj1.shape[0])
				if self.writeToFile:
					numpy.savetxt("%s_%s_accessibility.txt" %(pymolObject1, chain_1), numpy.column_stack((x,accObj1)))
				plotDictionary = self.makeGraphDataDictionary(id, fileName, "AccessibilityPlot", "Number of Residue", "Relative Accessibility", x, accObj1, 0, xlim, 0)
				stored.plots.append(plotDictionary)
				
				#Plot accessibility for chain 2
				fileName = "%s_%s_accessibilty" %(pymolObject1, chain_2)
				x = numpy.linspace(1, accObj2.shape[0], num = accObj2.shape[0])
				if self.writeToFile:
					numpy.savetxt("%s_%s_accessibility.txt" %(pymolObject1, chain_2), numpy.column_stack((x, accObj2)))
				plotDictionary = self.makeGraphDataDictionary(id, fileName, "AccessibilityPlot", "Number of Residue", "Relative Accessibility", x, accObj2, 0, ylim, 0)
				stored.plots.append(plotDictionary)

				#diagonal, only if both chains have the same number of residues
				if distances.shape[0] == distances.shape[1]:
					fileName = "%s_%s-%s_diagonal" %(pymolObject1,chain_1,chain_2)
					x = numpy.linspace(1, diagonal.shape[0], num = diagonal.shape[0])
					if self.writeToFile:
						numpy.savetxt(fileName+".txt", numpy.column_stack((x, diagonal)))
					plotDictionary = self.makeGraphDataDictionary(id, fileName, "DistancePlot", "Number of Residue", "Distance (Angstrom)", x, diagonal, 0, xlim, ylim)
					stored.plots.append(plotDictionary)