import copy
import numpy
from chromosome import Chromosome
from population import Population
import scipy.spatial.distance

class Environment:

	def __init__(self, populations, proteinA, proteinB, expDistances, expErrors, weights, scoreClashes):
		self.populations = populations
		self.proteinA = proteinA
		self.proteinB = proteinB
		self.expDistances = expDistances
		self.expErrors = expErrors
		self.weights = weights
		self.scoreClashes = scoreClashes
		#self.scoreCOGdiff = False
		self.constraintNames = []
		print "Score clashes: ", self.scoreClashes

	def quickClash(self, atomsProteinA, atomsProteinB, cutoff):
		clashes = 0
		#print len(atomsProteinA), len(atomsProteinB)
		try:
			dist = scipy.spatial.distance.cdist(atomsProteinA, atomsProteinB)
			#print atomsProteinA, atomsProteinB
			#print dist
			#print dist < cutoff
			#print numpy.nonzero(dist < cutoff)[0]
			#clashes = len(numpy.nonzero(dist < cutoff)[0]
			clashes = dist[(dist < cutoff) & (dist > 0)]
			#print clashes
		except:
			print "Error while determining clashes... Just one atom in protein?"
		return len(clashes)
	
	def quickDist2(self, atoms1, atoms2):
		# if there is only one atom it has to be duplicated for quick_dist2 to work
		duplicatedA = False
		duplicatedB = False
		if len(atoms1) == 1:
			duplicatedA = True
			atoms1 = numpy.tile(atoms1, (2, 1))
		if len(atoms2) == 1:
			duplicatedB = True
			atoms2 = numpy.tile(atoms2, (2, 1))
		dist = scipy.spatial.distance.cdist(atoms1, atoms2)
		dist = numpy.diag(dist)
		return dist

	#-------------------------------------------------------------------------------------

	def getFitness(self, trialDistances, trialAtomsProteinA, trialAtomsProteinB):
		#print len(trialAtomsProteinA),len(trialAtomsProteinB)
		rmsd = 0
		chi2 = 0
		clashes = 0
		cogA = numpy.average(trialAtomsProteinA, axis=0)
		cogB = numpy.average(trialAtomsProteinB, axis=0)
		cogDiff = numpy.linalg.norm(cogA - cogB)
		#print cogDiff
		rmsd = numpy.sqrt(numpy.mean(numpy.square(self.expDistances - trialDistances)))
		chi2 = numpy.sum(numpy.square((self.expDistances - numpy.sqrt(numpy.square(trialDistances))) / self.expErrors))
		if self.scoreClashes:
			clashes = self.quickClash(trialAtomsProteinA, trialAtomsProteinB, 3.5)
			#print len(trialAtomsProteinA),len(trialAtomsProteinB)
			#print "clashes: ", clashes

		fitnessResults = [rmsd, chi2, clashes, cogDiff]
		return fitnessResults

	#-------------------------------------------------------------------------------------

	def applySelectionPressure(self, idx, generation=0):
		population = self.populations[idx]
		# Determine fitness of childs
		counter = 0
		newBest = False
		for chromosome in population.chromosomes:
			if chromosome.fitness == 10000000000:
				tmpProteinB = copy.copy(self.proteinB)
				if chromosome.symmetry == "None":
					originalPosition = tmpProteinB.calculatePositionFromChromosome(chromosome, onlyLabelAtoms=True)
					trialDistances = self.quickDist2(self.proteinA.labelAtoms, tmpProteinB.labelAtoms)
				else:
					originalPosition = tmpProteinB.calculatePositionFromChromosome(chromosome, onlyLabelAtoms=True)
					trialDistances = self.quickDist2(originalPosition, tmpProteinB.labelAtoms)


				if self.scoreClashes:
					originalPosition = tmpProteinB.calculatePositionFromChromosome(chromosome, onlyLabelAtoms=False)
					if chromosome.symmetry != "None":
						fitness = self.getFitness(trialDistances, originalPosition, tmpProteinB.allAtoms)
					else:
						fitness = self.getFitness(trialDistances, self.proteinA.allAtoms, tmpProteinB.allAtoms)
				else:
					fitness = self.getFitness(trialDistances, originalPosition, tmpProteinB.labelAtoms)
				chromosome.rmsd = fitness[0]
				chromosome.chi2 = fitness[1]
				chromosome.clashes = fitness[2]
				chromosome.cogDiff = fitness[3]
				
				# add clash penalty if more than 5 clashes
				if self.scoreClashes == True and chromosome.clashes > 5:
					chromosome.fitness = chromosome.chi2 + chromosome.clashes * chromosome.chi2
				else:
					chromosome.fitness = chromosome.chi2
				
				#if self.scoreClashes == False:
				#	chromosome.log += chromosome.printChromosomeWithClashes()
				#else:
				#	chromosome.log += chromosome.printChromosomeWithoutClashes()
				
				if chromosome.fitness < population.bestChromosome.fitness:
					population.bestChromosome = chromosome
					population.bestChromosome.generationNumber = generation
					newBest = True
				chromosome.trialDistances = trialDistances
			counter += 1
		if newBest:
			population.log += population.bestChromosome.printChromosomeWithClashes()
			#if self.scoreClashes:
			#	chromosome.generationNumber = generation
			#	population.log += population.bestChromosome.printChromosomeWithClashes()
			#else:
			#	chromosome.generationNumber = generation
			#	population.log += population.bestChromosome.printChromosomeWithoutClashes()

#-----------------------------------------------------------------------------------------
