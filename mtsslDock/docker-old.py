import numpy
from protein import Protein
from chromosome import Chromosome
from population import Population
from environment import Environment
import copy
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
		self.symmetry = settings['symmetry']
		self.farAwayPenalty = True
		self.scoreVdw = True
		#self.symmetry = symmetry
		self.restraints = restraints
		self.objectPrefix = "mD"
		pub.subscribe(self.onAbortDocking, "docking.abort")
	
	def worker1(self, environment, idx, resultQueue, progressQueue):
		cProfile.runctx("self.worker(environment, idx, resultQueue, progressQueue)", globals(), locals())

	def worker(self, environment, idx, resultQueue, progressQueue):
		population = environment.populations[idx]
		environment.applySelectionPressure(idx)
		population.log += population.chromosomes[0].printChromosomeWithoutClashes()
		# first evolution
		for i in range(0, self.numberOfGenerations):
			#progressQueue.put(1.0)
			population.sacrifice(0.05)
			population.produceOffspring(rigidBody=False)
			environment.applySelectionPressure(idx)
			population.log += population.chromosomes[0].printChromosomeWithoutClashes()
		population.sacrifice(0.95)
		population.resetFitness()
		environment.scoreClashes = self.scoreClashes
		environment.farAwayPenalty = self.farAwayPenalty
		#cProfile.runctx("environment.applySelectionPressure(idx)", globals(), locals())
		environment.applySelectionPressure(idx)
		population.log += population.chromosomes[0].printChromosomeWithoutClashes()
		environment.scoreVdw = self.scoreVdw
		# rigid body
		for i in range(0, self.numberOfRigidBodyCycles):
			#progressQueue.put(1.0)
			population.sacrifice(0.95)
			population.produceOffspring(rigidBody=True)
			environment.applySelectionPressure(idx)
			population.log += population.chromosomes[0].printChromosomeWithClashes()
		resultQueue.put(population)
		
	def onAbortDocking(self):
		for p in self.processes:
			p.terminate()

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

		#######
		# evolve
		#######

		zeroChromosome = Chromosome("None")
		zeroChromosome.genes = numpy.array([0, 0, 0, 0, 0, 0])

		# setup populations
		print "setting up populations..."
		
		
		populations = []
		for i in range(0, self.numberOfPopulations):
			if self.symmetry == "C2":
				population = Population(self.numberOfChromosomes, "C2")
			elif self.symmetry == "None":
				population = Population(self.numberOfChromosomes, "None")
			population.name = "%i" % (i + 1)
			populations.append(population)

		# put them into an environment and evolve
		environment1 = Environment(populations, proteinA, proteinB, expDistances, expErrors, weights, False, False, False)
		environment1.constraintNames = constraintNames
		#environment1.applySelectionPressure()
		#for population in environment1.populations:
		#	population.log += population.chromosomes[0].printChromosomeWithoutClashes()


		self.processes = []
		resultQueue = multiprocessing.Queue()
		progressQueue = multiprocessing.Queue()
		numberOfProcesses = len(environment1.populations)
		for idx, population in enumerate(environment1.populations):
				environment = copy.deepcopy(environment1)
				p = multiprocessing.Process(target = self.worker1, args = (environment, idx, resultQueue, progressQueue))
				p.start()
				self.processes.append(p)
		
		cycles = 0
		maxCycles = numberOfProcesses * (self.numberOfGenerations + self.numberOfRigidBodyCycles)
		#while True:
		#	cycles += progressQueue.get()
		#	progress = cycles/(numberOfProcesses*(self.numberOfGenerations+self.numberOfRigidBodyCycles))
		#	#send message to main thread
		#	wx.CallAfter(pub.sendMessage, "docking.update", progress=progress)
		#	if cycles >= maxCycles:
		#		break
		resultsList = [resultQueue.get() for p in self.processes]
		for p in self.processes:
			p.join()
		environment1.populations = resultsList

		# create solutions
		print ""
		print "Solutions:"
		nonClashingSolution = 1
		clashingSolution = 1
		for population in environment1.populations:
			#createPseudoatom(self.labelPositionsProteinB, "tmpSolution-labels", 1)
			tmpProtein = Protein(proteinB.originalLabelAtoms, proteinB.originalLabelAtoms, "tmpSolution-labels")
			solution = population.chromosomes[0]
			# print solution.printChromosomeWithClashes()
			if solution.clashes <= 5:
				nameOfSolution = "%s-%i_sol-%i" % (self.objectPrefix, self.dockingRunNumber, nonClashingSolution)
				solution.name = nameOfSolution

				proteinB.moveInPymol(nameOfSolution, solution, 1)
				#tmpProtein.moveInPymol("%s-labels" % nameOfSolution, solution, 1)
				cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 1, 0, None)
				#cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), "%s-labels" % nameOfSolution, 1, 0, None)
				nonClashingSolution += 1

			elif solution.clashes > 5:
				nameOfSolution = "%s-%i_clash-%i" % (self.objectPrefix, self.dockingRunNumber, clashingSolution)
				solution.name = nameOfSolution
				proteinB.moveInPymol(nameOfSolution, solution, 1)
				#tmpProtein.moveInPymol("%s-labels" % nameOfSolution, solution, 1)
				cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 1, 0, None)
				#cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), "%s-labels" % nameOfSolution, 1, 0, None)
				clashingSolution += 1
		cmd.group("%s-%i" % (self.objectPrefix, self.dockingRunNumber), "%s-%i*" % (self.objectPrefix, self.dockingRunNumber))
		#cmd.set_view(myView)
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