import multiprocessing
import math
import numpy
import scipy.spatial.distance
import sys, os

class DistanceDistribution:
	def __init__(self):
		self.name = ""
		self.labels = {}
		self.spinCenters = {}
	
	def getChunks(self, numberOfItems, numberOfProcesses):
		chunkSize = int(math.ceil(numberOfItems/float(numberOfProcesses)))
		chunks = []
		sum = 0
		i = 0
		while sum + chunkSize < numberOfItems:
			chunks.append(chunkSize)
			sum += chunkSize
			i += 1
		chunks.append(numberOfItems - i * chunkSize)
		#print chunks
		return chunks

	
	def calculateDistanceDistribution(self, spinCenter1Coordinates, spinCenter2Coordinates):
		#Setup parallel computation of distances	
		processes = []
		numberOfProcesses = multiprocessing.cpu_count()
		distanceQueue = multiprocessing.Queue()
		chunks = []
		reverse = False
		if len(spinCenter1Coordinates) >= len(spinCenter2Coordinates):
			chunks = self.getChunks(len(spinCenter1Coordinates), numberOfProcesses)
		else:
			chunks = self.getChunks(len(spinCenter2Coordinates), numberOfProcesses)
			reverse = True
			
		startPosition = 0
		stopPosition = chunks[0]
		outOfWork = False
		#only use multiprocessing on mac or linux
		if os.name != "nt":
			for i in range (numberOfProcesses):
				#get chunk of data
				#print i, startPosition, stopPosition
				if not reverse:
					chunk = spinCenter1Coordinates[startPosition:stopPosition]
				else:
					chunk = spinCenter2Coordinates[startPosition:stopPosition]
				if i < numberOfProcesses - 1:
					try:
						startPosition += chunks[i]
						stopPosition += chunks[i+1]
					except:
						outOfWork = True
			
				#calculate distances
				if not reverse and not outOfWork:
					p = multiprocessing.Process(target = self.quickDistMulti, args = (chunk, spinCenter2Coordinates, distanceQueue))
					p.start()
					processes.append(p)
				elif not outOfWork:
					p = multiprocessing.Process(target = self.quickDistMulti, args = (chunk, spinCenter1Coordinates, distanceQueue))
					p.start()
					processes.append(p)
				elif outOfWork:
					numberOfProcesses -= 1
		else:
			numberOfProcesses = 1
			self.quickDistMulti(spinCenter1Coordinates, spinCenter2Coordinates, distanceQueue)
			
		distanceDictionary = {}
		for i in range(numberOfProcesses):
			distanceDictionary.update(distanceQueue.get())
		
		for p in processes:
			p.join()
		
		print "Done."
		print "Collecting results...", 
		dist = numpy.vstack(distanceDictionary.values())
		print "Done!"
		return dist

	def quickDistMulti(self, atoms1, atoms2, distanceQueue):
		try:
			currentProcess = multiprocessing.current_process()
			uniqueId = currentProcess._identity
		except:
			pass
		distanceDictionary = {}
		#print len(atoms1), len(atoms2)
		
		# if there is only one atom it has to be duplicated for quick_dist to work
		duplicated = False
		atoms1 = numpy.array(atoms1)
		atoms2 = numpy.array(atoms2)
		
		if len(numpy.shape(atoms1)) == 1:
			duplicated = True
			atoms1=numpy.tile(atoms1, (2,1))
		if len(numpy.shape(atoms2)) == 1:
			duplicated = True
			atoms2=numpy.tile(atoms2, (2,1))
		
		dist = scipy.spatial.distance.cdist(atoms1, atoms2)
	
		#remove the duplication depending on which selection contained the single atom
		if duplicated and dist.shape[0] == 2 and not dist.shape[1] == 2:
			dist=numpy.reshape(dist[0,:], (-1, 1))
	
		elif duplicated and dist.shape[1] == 2 and not dist.shape[0] == 2:
			dist=numpy.reshape(dist[:,0], (-1, 1))
	
		elif duplicated and dist.shape[0] == 2 and dist.shape[1] == 2:
			dist=numpy.reshape(dist[:1,0], (-1, 1))
		else:
			dist=numpy.reshape(dist, (-1, 1))
		
		distanceDictionary[uniqueId] = dist
		distanceQueue.put(distanceDictionary)