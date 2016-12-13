import numpy
from pymol import cmd
import math
import os
import cPickle as pickle

class Protein:

	#-------------------------------------------------------------------------------------

	def __init__(self, labelAtoms, allAtoms, pymolString):
		self.originalAllAtoms = numpy.copy(allAtoms)
		self.originalLabelAtoms = numpy.copy(labelAtoms)
		self.labelAtoms = numpy.copy(labelAtoms)
		self.allAtoms = numpy.copy(allAtoms)
		self.labelAtomsCog = self.getCog(self.labelAtoms)
		self.pymolString = pymolString
		path = os.path.dirname(__file__)
		self.mat100 = pickle.load(open(path + "/mat100.pickle", 'rb'))
		self.mat010 = pickle.load(open(path + "/mat010.pickle", 'rb'))
		self.mat001 = pickle.load(open(path + "/mat001.pickle", 'rb'))
		self.approximateRadius = self.getApproximateRadius(allAtoms)

	#-------------------------------------------------------------------------------------

	def getApproximateRadius(self, atoms):
		cog = self.getCog(atoms)
		maxRadius = 0.0
		for atom in atoms:
			radius = numpy.linalg.norm(cog - atom)
			if radius > maxRadius:
				maxRadius = radius
		print maxRadius
		return maxRadius

	#-------------------------------------------------------------------------------------

	def getCog(self, atoms):
		return numpy.average(atoms, axis=0)

	#-------------------------------------------------------------------------------------

	def moveInPymol(self, solutionPymolString, chromosome, state):
		# recreate protein B at solution coordinates
		cmd.create(solutionPymolString, self.pymolString, 1, state)
		# translate to origin with proteinBcog. IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system!
		cmd.translate(list(-1 * self.labelAtomsCog.reshape(-1,)), solutionPymolString, state, 0, None)

		# rotate and translate according to solution
		translation = chromosome.genes[3:6]
		rotX = chromosome.genes[0]
		# IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system! Also set rotation origin to 0,0,0!
		cmd.rotate([1, 0, 0], rotX, solutionPymolString, state, 0, None, [0, 0, 0])
		rotY = chromosome.genes[1]
		cmd.rotate([0, 1, 0], rotY, solutionPymolString, state, 0, None, [0, 0, 0])
		rotZ = chromosome.genes[2]
		cmd.rotate([0, 0, 1], rotZ, solutionPymolString, state, 0, None, [0, 0, 0])
		cmd.translate(list(translation.reshape(-1,)), solutionPymolString, state, 0, None)

		if chromosome.symmetry != "None":
			if chromosome.symmetry == "C2":
				monomers = 2
			elif chromosome.symmetry == "C3":
				monomers = 3
			elif chromosome.symmetry == "C4":
				monomers = 4
			elif chromosome.symmetry == "C5":
				monomers = 5
			elif chromosome.symmetry == "C6":
				monomers = 6
			elif chromosome.symmetry == "C7":
				monomers = 7
			elif chromosome.symmetry == "C8":
				monomers = 8
			angle = 360.0/monomers
			for i in range(1, monomers):
				cmd.create(solutionPymolString, solutionPymolString, 1, state + i)
				cmd.rotate([1, 0, 0], i*angle, solutionPymolString, state+i, 0, None, [0, 0, 0])

	#-------------------------------------------------------------------------------------

	def moveToOrigin(self, cog):
		self.allAtoms = self.allAtoms - cog
		self.labelAtoms = self.labelAtoms - cog

	#-------------------------------------------------------------------------------------

	def calculatePositionFromChromosome(self, chromosome, onlyLabelAtoms):

		if onlyLabelAtoms:
			newPosition = numpy.copy(self.labelAtoms)
		else:
			newPosition = numpy.copy(self.allAtoms)

		orientationVector = numpy.array([[1, 0, 0]])
		zVector = numpy.array([[0, 0, 1]])
		# rotate around x
		rotationMatrix = self.mat100["%1.1f" %int(chromosome.genes[0])]
		newPosition = self.rotatePoints(newPosition, rotationMatrix)
		orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
		zVector = self.rotatePoints(zVector, rotationMatrix)

		# rotate around y
		rotationMatrix = self.mat010["%1.1f" %int(chromosome.genes[1])]
		newPosition = self.rotatePoints(newPosition, rotationMatrix)
		orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
		zVector = self.rotatePoints(zVector, rotationMatrix)

		# rotate around z
		rotationMatrix = self.mat001["%1.1f" %int(chromosome.genes[2])]
		newPosition = self.rotatePoints(newPosition, rotationMatrix)
		orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
		zVector = self.rotatePoints(zVector, rotationMatrix)

		# translate along x, y, z
		newPosition = newPosition + numpy.array([chromosome.genes[3], 0, 0])
		newPosition = newPosition + numpy.array([0, chromosome.genes[4], 0])
		newPosition = newPosition + numpy.array([0, 0, chromosome.genes[5]])
		originalPosition = newPosition
		
		# return position now, if no symmetry
		if chromosome.symmetry == "None":
			if onlyLabelAtoms:
				self.labelAtoms = newPosition
			else:
				self.allAtoms = newPosition
			return newPosition

		# create Cn symmetric dimer
		elif chromosome.symmetry != "None":
			if chromosome.symmetry == "C2":
				angle = 180
			elif chromosome.symmetry == "C3":
				angle = 120
			elif chromosome.symmetry == "C3":
				angle = 120
			elif chromosome.symmetry == "C4":
				angle = 90
			elif chromosome.symmetry == "C5":
				angle = 72
			elif chromosome.symmetry == "C6":
				angle = 60
			elif chromosome.symmetry == "C7":
				angle = 51.4
			elif chromosome.symmetry == "C8":
				angle = 45
			
			cnAxis = numpy.array([1, 0, 0])
			# rotate around Cn axis
			rotationMatrix = self.setupRotationMatrix(angle, numpy.array([1, 0, 0]))
			newPosition = self.rotatePoints(newPosition, rotationMatrix)
			if onlyLabelAtoms:
				self.labelAtoms = newPosition
			else:
				self.allAtoms = newPosition
			#the original Position is before the Cn symmetry is applied
			return originalPosition

	#-------------------------------------------------------------------------------------

	def getAngle(self, a, b):
		return numpy.rad2deg(numpy.arccos(numpy.dot(a / numpy.linalg.norm(a), b / numpy.linalg.norm(b))))
	#-------------------------------------------------------------------------------------

	def setupRotationMatrix(self, theta, axis):
		"""
		Return the rotation matrix associated with counterclockwise rotation about
		the given axis by theta radians.
		http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
		"""
		axis = numpy.asarray(axis)
		theta = numpy.asarray(numpy.deg2rad(theta))
		axis = axis/math.sqrt(numpy.dot(axis, axis))
		a = math.cos(theta/2.0)
		b, c, d = -axis*math.sin(theta/2.0)
		aa, bb, cc, dd = a*a, b*b, c*c, d*d
		bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
		return numpy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
						 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
						 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

	def rotatePoints(self, points, rotationMatrix):
		points = numpy.array(points)
		rotatedPoints = numpy.dot(points, rotationMatrix.T)
		return numpy.array(rotatedPoints)

