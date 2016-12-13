import numpy
class Chromosome:

	#-------------------------------------------------------------------------------------

	def __init__(self, symmetry):
		self.alive = True
		self.maxTranslation = 100
		self.symmetry = symmetry
		if self.symmetry == "None":
			self.genes = self.generateRandomGenes()
		else:
			self.genes = self.generateCnGenes()
		self.generationNumber = 0
		self.fitness = 10000000000
		self.chi2 = 10000000000
		self.rmsd = 10000000000
		self.clashes = 10000000000
		self.chi2red = 10000000000
		self.cogDiff = 10000000000
		self.log = ""
		self.trialDistances = numpy.zeros((3, 3))
		self.name = ""

	#-------------------------------------------------------------------------------------

	def getGenes(self):
		return self.genes

	#-------------------------------------------------------------------------------------

	def setGenes(self, rotx, roty, rotz, tx, ty, tz):
		self.genes = numpy.array([rotx, roty, rotz, tx, ty, tz])

	#-------------------------------------------------------------------------------------

	def generateRandomGenes(self):
		# translation between -maxTranslation and maxTranslation
		translation = self.maxTranslation * numpy.random.uniform(-1, 1, size=3)
		rotation = 360 * numpy.random.uniform(0, 1, size=3)
		return numpy.concatenate((rotation, translation), axis=0)

	#-------------------------------------------------------------------------------------

	def generateCnGenes(self):
		# translation between -maxTranslation and maxTranslation
		translation = self.maxTranslation * numpy.random.uniform(-1, 1, size=3)
		rotation = 360 * numpy.random.uniform(0, 1, size=3)
		translation[0] = 0
		return numpy.concatenate((rotation, translation), axis=0)

	#-------------------------------------------------------------------------------------
	def printChromosomeWithClashes(self):
		#string = "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t rmsd:%1.2f\t chi2:%1.2f\t clashes:%i\t cogDist:%1.2f\n" %(self.genes[0],self.genes[1],self.genes[2],self.genes[3],self.genes[4],self.genes[5], self.rmsd, self.chi2, self.clashes, self.cogDiff)
		string = "%i %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.4f %i %1.2f %1.2f\n" % (self.generationNumber, self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.genes[4], self.genes[5], self.rmsd, self.chi2, self.clashes, self.cogDiff, self.fitness)
		# print string
		return string

	def printChromosomeWithoutClashes(self):
		#string = "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t rmsd:%1.2f\t chi2:%1.2f\t clashes:%i\n" %(self.genes[0],self.genes[1],self.genes[2],self.genes[3],self.genes[4],self.genes[5], self.rmsd, self.chi2, self.clashes)
		string = "%i %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.4f - %1.2f %1.2f\n" % (self.generationNumber, self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.genes[4], self.genes[5], self.rmsd, self.chi2, self.cogDiff, self.fitness)
		# print string
		return string