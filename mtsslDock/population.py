from chromosome import Chromosome
import random
from operator import attrgetter
import copy

class Population:

	#-------------------------------------------------------------------------------------

	def __init__(self, size, symmetry):
		self.size = size
		self.symmetry = symmetry
		self.chromosomes = self.randomSpawn(symmetry)
		self.sacrificed = 0
		self.rank()
		self.bestChromosome = self.chromosomes[0]
		self.name = ""
		self.log = ""

	#-------------------------------------------------------------------------------------

	def resetFitness(self):
		for chromosome in self.chromosomes:
			chromosome.fitness = 10000000000
			self.bestChromosome.fitness = 10000000000

	#-------------------------------------------------------------------------------------

	def randomSpawn(self, symmetry):
		chromosomes = []
		for i in range(0, self.size):
			chromosome = Chromosome(symmetry)
			chromosomes.append(chromosome)
		return chromosomes

	#-------------------------------------------------------------------------------------

	def rank(self):
		self.chromosomes = sorted(self.chromosomes, key=attrgetter('fitness'))

	#-------------------------------------------------------------------------------------

	def sacrifice(self, mutationFrequency):
		numberToSacrifice = int(len(self.chromosomes) * mutationFrequency)
		for i in range(numberToSacrifice):
			delinquent = self.tournament()["looser"]
			self.chromosomes.pop(self.tournament()["looserIndex"])
		self.sacrificed = numberToSacrifice

	#-------------------------------------------------------------------------------------

	def printChromosomes(self):
		for chromosome in self.chromosomes:
			print chromosome.printChromosomeWithoutClashes()

	#-------------------------------------------------------------------------------------

	def getGeneticOperator(self):
		if random.choice([True, False]):
			if random.choice([True, False]):
				return "smallCreepMutation"
			else:
				return "randomMutation"
		else:
			if random.choice([True, False]):
				return "singlePointCrossover"
			else:
				return "exchangeCrossover"

	#-------------------------------------------------------------------------------------

	def smallCreepMutation(self, parent):
		childs = []
		child = Chromosome(self.symmetry)
		genes = parent.getGenes()
		child.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		position = random.randrange(6)
		oldValue = parent.genes[position]
		newValue = 0
		# rotation
		if position < 3:
			creep = 0
			while creep == 0:
				creep = random.uniform(-5.0, 5.0)
			newValue = oldValue + creep
			if newValue > 360:
				newValue -= 360
			elif newValue < 0:
				newValue += 360
		# translation
		else:
			creep = 0
			while creep == 0:
				creep = random.uniform(-2.0, 2.0)
			newValue = oldValue + creep
		child.genes[position] = newValue
		child.fitness = 10000000000
		childs.append(child)
		return childs

	#-------------------------------------------------------------------------------------

	def randomMutation(self, parent):
		childs = []
		child = Chromosome(self.symmetry)
		genes = parent.getGenes()
		child.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		position = random.randrange(6)
		oldValue = parent.genes[position]
		newValue = 0
		# rotation
		if position < 3:
			newValue = random.randrange(360)
		# translation
		else:
			newValue = random.randrange(-50, 50)
		child.genes[position] = newValue
		child.fitness = 10000000000
		childs.append(child)
		return childs

	#-------------------------------------------------------------------------------------

	def singlePointCrossover(self, parent1, parent2):
		childs = []
		position = random.randrange(6)
		child1 = Chromosome(self.symmetry)
		genes = parent1.getGenes()
		child1.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		child2 = Chromosome(self.symmetry)
		genes = parent2.getGenes()
		child2.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		for i in range(0, 5):
			if i <= position:
				child1.genes[i] = parent1.genes[i]
				child2.genes[i] = parent2.genes[i]
			else:
				child1.genes[i] = parent2.genes[i]
				child2.genes[i] = parent1.genes[i]
		child1.fitness = 10000000000
		child2.fitness = 10000000000
		childs.append(child1)
		childs.append(child2)
		return childs

	#-------------------------------------------------------------------------------------

	def exchangeCrossover(self, parent1, parent2):
		childs = []
		position = random.randrange(6)
		child1 = Chromosome(self.symmetry)
		genes = parent1.getGenes()
		child1.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		child2 = Chromosome(self.symmetry)
		genes = parent2.getGenes()
		child2.setGenes(genes[0], genes[1], genes[2], genes[3], genes[4], genes[5])
		tmpValue = child1.genes[position]
		child1.genes[position] = child2.genes[position]
		child2.genes[position] = tmpValue
		child1.fitness = 10000000000
		#child1.log = ''
		child2.fitness = 10000000000
		#child2.log = ''
		childs.append(child1)
		childs.append(child2)
		return childs

	def tournament(self, size=5):
		#print "New tournament:"
		participants = []
		worstFitness = 0
		bestFitness = 10000000000
		resultDict = {}
		#select participants of this tournament
		for i in range(size):
			index = random.randrange(len(self.chromosomes))
			participant = (self.chromosomes[index])
			participants.append(participant)
		#set initial looser and winner
		currentWinner = participants[0]
		currentLooser = participants[1]
		for participant in participants:
			if participant.fitness < bestFitness:
				currentWinner = participant
				currentWinnerIndex = index
				bestFitness = participant.fitness
			if participant.fitness > worstFitness:
				currentLooser = participant
				currentLooserIndex = index
				worstFitness = participant.fitness
		try:
			resultDict["winner"] = currentWinner
			resultDict["looser"] = currentLooser
			resultDict["winnerIndex"] = currentWinnerIndex
			resultDict["looserIndex"] = currentLooserIndex
		except:
			print "Could not determine tournament winner and/or looser!"
			
			
		return resultDict

	#-------------------------------------------------------------------------------------

	def produceOffspring(self, rigidBody):
		j = 0
		while j < self.sacrificed:
			geneticOperator = ""

			# make sure smallCreepMutation is selected for rigid body refinements
			if not rigidBody:
				geneticOperator = self.getGeneticOperator()
			else:
				geneticOperator = "smallCreepMutation"
			childs = []
			if geneticOperator == "smallCreepMutation":
				if rigidBody:
					# take fittest chromosome as parent for rigid body
					parent = self.bestChromosome
				else:
					parent = self.tournament()["winner"]
				childs = self.smallCreepMutation(parent)
				j += 1

			elif geneticOperator == "randomMutation":
				parent = self.tournament()["winner"]
				childs = self.randomMutation(parent)
				j += 1

			elif geneticOperator == "singlePointCrossover" and (self.sacrificed - j <= 2):
				parent1 = self.tournament()["winner"]
				parent2 = self.tournament()["winner"]
				childs = self.singlePointCrossover(parent1, parent2)
				j += 2

			elif geneticOperator == "exchangeCrossover" and (self.sacrificed - j <= 2):
				parent1 = self.tournament()["winner"]
				parent2 = self.tournament()["winner"]
				childs = self.exchangeCrossover(parent1, parent2)
				j += 2

			# force translation along twofold (global x-Axis) to 0 for Cn symmetry to avoid screw axis
			if self.symmetry != "None":
				for child in childs:
					child.genes[3] = 0

			self.chromosomes += childs

#-----------------------------------------------------------------------------------------
