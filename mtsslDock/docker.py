import numpy
from protein import Protein
from chromosome import Chromosome
from population import Population
from environment import Environment
import copy
import os
import multiprocessing
import wx
from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub
from pymol import cmd
import cProfile

class Docker:
	def __init__(self,
			dockingRunNumber,
			restraints,
			settings):
		self.dockingRunNumber = dockingRunNumber
		self.settings = settings
		self.numberOfPopulations = settings['numberOfPopulations']
		self.numberOfGenerations = settings['numberOfGenerations']
		self.numberOfChromosomes = settings['numberOfChromosomes']
		self.numberOfRigidBodyCycles = settings['numberOfRigidBodyCycles']
		self.scoreClashes = settings['scoreClashes']
		self.scoreCOGdiff = settings['scoreCOGdiff']
		self.symmetry = settings['symmetry']
		self.scoreVdw = True
		self.abort = False
		self.restraints = restraints
		self.objectPrefix = "mD"
		pub.subscribe(self.onAbortDocking, "docking.abort")

	#-------------------------------------------------------------------------------------	

	#for testing purposes
	def worker1(self, environment, idx, resultQueue, progressQueue):
		cProfile.runctx("self.worker(environment, idx, resultQueue, progressQueue)", globals(), locals())

	#-------------------------------------------------------------------------------------

	def worker(self, environment, idx, resultQueue, progressQueue):
		#print "clashes: ", self.scoreClashes
		#print "cogDiff: ", self.scoreCOGdiff
		#environment.scoreClashes = self.scoreClashes
		#print environment.scoreClashes
		#environment.scoreCOGdiff = self.scoreCOGdiff
		if idx >= 0:
			population = environment.populations[idx]
			environment.applySelectionPressure(idx)
			population.rank()
			population.log += population.chromosomes[0].printChromosomeWithoutClashes()
			# step one
			generation = 0
			for i in range(0, self.numberOfGenerations):
				progressQueue.put(1.0)
				population.sacrifice(0.05)
				population.produceOffspring(rigidBody=False)
				environment.applySelectionPressure(idx, generation=generation)
				generation += 1
			
			# step two: rigid body, only small changes are made
			for i in range(0, self.numberOfRigidBodyCycles):
				progressQueue.put(1.0)
				population.sacrifice(0.95)
				population.produceOffspring(rigidBody=True)
				environment.applySelectionPressure(idx, generation=generation)
				generation += 1
			population.rank()
			resultQueue.put(population)
		#Windows...
		elif idx == -1:
			results = []
			cycles = 0.0
			factor = len(environment.populations)*self.numberOfGenerations+self.numberOfRigidBodyCycles
			for idy, population in enumerate(environment.populations):
				environment.applySelectionPressure(idy)
				population.rank()
				population.log += population.chromosomes[0].printChromosomeWithoutClashes()
				# step one
				generation = 0
				for i in range(0, self.numberOfGenerations):
					if not self.abort:
						wx.CallAfter(pub.sendMessage, "docking.update", progress=cycles/factor)
						population.sacrifice(0.05)
						population.produceOffspring(rigidBody=False)
						environment.applySelectionPressure(idy, generation=generation)
						generation += 1
						cycles += 1.0
				population.rank()
				# step two: rigid body
				for i in range(0, self.numberOfRigidBodyCycles):
					if not self.abort:
						wx.CallAfter(pub.sendMessage, "docking.update", progress=cycles/factor)
						population.sacrifice(0.95)
						population.produceOffspring(rigidBody=True)
						environment.applySelectionPressure(idy, generation=generation)
						generation += 1
						cycles += 1.0
				population.rank()
						
	#-------------------------------------------------------------------------------------
		
	def onAbortDocking(self):
		self.abort = True
		for p in self.processes:
			p.terminate()

	#-------------------------------------------------------------------------------------

	def dock(self):
		#setup tables with distances, errors and weights
		expDistances = []
		expErrors = []
		weights = []
		constraintNames = []
		
		for restraint in self.restraints:
			constraintNames.append("%s-%s"%(restraint["anchorAname"], restraint["anchorBname"]))
			expDistances.append(restraint["distance"])
			expErrors.append(restraint["width"])
			weights.append(restraint["weight"])
		
		expDistances = numpy.array(expDistances)
		expErrors = numpy.array(expErrors)
		weights = numpy.array(weights)
		
		# setup Proteins
		anchorAcoords = []
		for restraint in self.restraints:
			anchorAcoords.append(restraint["anchorAcoord"])
		
		proteinAname = self.restraints[0]["proteinAname"] 
		anchorAcalphas = self.restraints[0]["proteinAcalpha"]
		proteinA = Protein(anchorAcoords, anchorAcalphas, proteinAname)

		anchorBcoords = []
		for restraint in self.restraints:
			anchorBcoords.append(restraint["anchorBcoord"])
		proteinBname = self.restraints[0]["proteinBname"] 
		anchorBcalphas = self.restraints[0]["proteinBcalpha"]
		proteinB = Protein(anchorBcoords, anchorBcalphas, proteinBname)

		# move both to origin
		proteinA.moveToOrigin(proteinA.labelAtomsCog)
		proteinB.moveToOrigin(proteinB.labelAtomsCog)

		# setup populations
		print "Starting..."
		
		populations = []
		for i in range(0, self.numberOfPopulations):
			if self.symmetry != "None":
				population = Population(self.numberOfChromosomes, self.symmetry)
			elif self.symmetry == "None":
				population = Population(self.numberOfChromosomes, "None")
			population.name = "%i" % (i + 1)
			populations.append(population)

		# put them into an environment and evolve
		environment1 = Environment(populations, proteinA, proteinB, expDistances, expErrors, weights, self.scoreClashes)
		environment1.constraintNames = constraintNames

		self.processes = []
		resultQueue = multiprocessing.Queue()
		progressQueue = multiprocessing.Queue()
		numberOfProcesses = len(environment1.populations)
		if os.name != "nt":
			for idx, population in enumerate(environment1.populations):
					environment = copy.deepcopy(environment1)
					p = multiprocessing.Process(target = self.worker, args = (environment, idx, resultQueue, progressQueue))
					p.start()
					self.processes.append(p)
					cycles = 0
					maxCycles = numberOfProcesses * (self.numberOfGenerations + self.numberOfRigidBodyCycles)
			while True:
					cycles += progressQueue.get()
					progress = cycles/(numberOfProcesses*(self.numberOfGenerations+self.numberOfRigidBodyCycles))
					#send message to main thread
					wx.CallAfter(pub.sendMessage, "docking.update", progress=progress)
					if cycles >= maxCycles:
							break
			resultsList = [resultQueue.get() for p in self.processes]
			for p in self.processes:
					p.join()
			environment1.populations = resultsList
		else:
			print "Windows... Using 1 core."
			self.worker(environment1, -1, resultQueue, progressQueue)
			self.abort = False
		
		# name solutions
		nonClashingSolution = 1
		clashingSolution = 1
		for population in environment1.populations:
			solution = population.chromosomes[0]
			if solution.clashes <= 5:
				nameOfSolution = "%s-%i_sol-%i" % (self.objectPrefix, self.dockingRunNumber, nonClashingSolution)
				solution.name = nameOfSolution
				nonClashingSolution += 1
			
			elif solution.clashes > 5:
				nameOfSolution = "%s-%i_clash-%i" % (self.objectPrefix, self.dockingRunNumber, clashingSolution)
				solution.name = nameOfSolution
				clashingSolution += 1
		return environment1, self.settings


	#-----------------------------------------------------------------------------------------

	def createPseudoatom(coordinates, objectName, state):
		for coordinate in coordinates:
			x = float(coordinate[0])
			y = float(coordinate[1])
			z = float(coordinate[2])
			posString = "[%3.2f,%3.2f,%3.2f]" % (x, y, z)
			cmd.pseudoatom(pos=posString, object=objectName, state=state)

	#-----------------------------------------------------------------------------------------
