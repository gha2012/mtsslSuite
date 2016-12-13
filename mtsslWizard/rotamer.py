import atom
import numpy

class Rotamer:	
	def __init__(self):
		self.id = 0
		self.chi2 = 0
		self.scaledChi2 = 0
		self.atoms = {}
		self.atomNames = ["SG","SD","CE","C2","C3","C4","C5","C6","C7","C8","C9","O1","N1"]
		#self.weights =   [5.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
		self.weights =   [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
		self.totalContacts = 0
		self.atomCoordinates = []
		self.scaledMeanValues = {"SG":[0.0,0.1,0.0,0.2,0.7,0.0,0.6,0.4,0.0],
						   "SD":[0.0,0.0,0.0,0.5,0.4,0.0,0.6,0.4,0.0],
						   "CE":[0.0,0.0,0.0,0.4,0.4,0.0,0.6,0.3,0.0],
						   "C2":[0.0,0.0,0.0,0.0,0.2,0.0,0.4,0.4,0.0],
						   "C3":[0.0,0.0,0.0,0.2,0.2,0.0,0.6,0.3,0.0],
						   "C4":[0.0,0.0,0.0,0.3,0.4,0.1,0.5,0.4,0.0],
						   "C5":[0.0,0.0,0.0,0.0,0.2,0.0,0.5,0.4,0.1],
						   "C6":[0.0,0.0,0.0,0.3,0.4,0.1,0.5,0.4,0.1],
						   "C7":[0.0,0.0,0.0,0.3,0.5,0.0,0.5,0.3,0.1],
						   "C8":[0.0,0.0,0.0,0.3,0.4,0.0,0.5,0.4,0.0],
						   "C9":[0.0,0.0,0.0,0.4,0.5,0.0,0.6,0.3,0.1],
						   "O1":[0.0,0.0,0.0,0.0,0.3,0.1,0.4,0.4,0.1],
						   "N1":[0.0,0.2,0.0,0.4,0.2,0.1,0.5,0.4,0.1]}
		self.scaledStdDevs = {"SG":[0.0,0.3,0.0,0.3,0.3,0.1,0.2,0.1,0.0],
						"SD":[0.0,0.2,0.0,0.3,0.3,0.1,0.1,0.2,0.1],
						"CE":[0.2,0.2,0.0,0.3,0.3,0.0,0.2,0.2,0.0],
						"C2":[0.0,0.0,0.0,0.1,0.4,0.1,0.3,0.3,0.2],
						"C3":[0.0,0.0,0.0,0.3,0.3,0.2,0.3,0.3,0.1],
						"C4":[0.0,0.0,0.0,0.4,0.4,0.2,0.3,0.3,0.1],
						"C5":[0.0,0.0,0.0,0.0,0.4,0.0,0.3,0.3,0.2],
						"C6":[0.1,0.1,0.0,0.3,0.4,0.1,0.3,0.2,0.2],
						"C7":[0.1,0.1,0.0,0.3,0.4,0.1,0.3,0.2,0.2],
						"C8":[0.0,0.0,0.0,0.3,0.4,0.1,0.3,0.3,0.1],
						"C9":[0.0,0.2,0.0,0.4,0.4,0.0,0.3,0.3,0.2],
						"O1":[0.0,0.0,0.0,0.2,0.5,0.2,0.3,0.3,0.2],
						"N1":[0.0,0.4,0.2,0.4,0.4,0.2,0.2,0.2,0.2]}
		self.totalMeanValues = {"SG":[0.0,0.1,0.0,0.7,2.4,0.1,6.4,4.9,0.2],
							    "SD":[0.0,0.0,0.0,1.8,1.9,0.1,4.9,3.0,0.1],
							    "CE":[0.0,0.0,0.0,1.3,0.7,0.0,4.2,2.6,0.1],
							    "C2":[0.0,0.0,0.0,0.0,0.2,0.0,3.1,2.3,0.2],
							    "C3":[0.0,0.0,0.0,0.3,0.2,0.0,3.0,1.7,0.1],
							    "C4":[0.0,0.0,0.0,0.6,0.6,0.0,3.1,2.0,0.2],
							    "C5":[0.0,0.0,0.0,0.0,0.2,0.0,3.2,2.3,0.5],
							    "C6":[0.0,0.0,0.0,1.3,1.2,0.4,3.3,1.9,0.3],
							    "C7":[0.0,0.0,0.0,0.8,1.2,0.2,3.6,2.0,0.5],
							    "C8":[0.0,0.0,0.0,1.1,1.4,0.1,3.7,2.2,0.2],
							    "C9":[0.0,0.0,0.0,1.6,1.6,0.1,4.3,2.5,0.3],
							    "O1":[0.0,0.2,0.0,1.0,0.4,0.2,2.8,2.5,0.5],
							    "N1":[0.0,0.0,0.0,0.1,0.4,0.1,2.2,1.9,0.3]}
		self.totalStdDevs = {"SG":[0.0,0.3,0.0,0.9,1.6,0.3,3.4,2.4,0.4],
							 "SD":[0.0,0.2,0.0,1.5,1.5,0.3,2.5,2.2,0.3],
							 "CE":[0.2,0.2,0.2,1.4,0.7,0.2,2.6,2.0,0.4],
							 "C2":[0.0,0.0,0.0,0.2,0.5,0.2,2.7,1.9,0.5],
							 "C3":[0.0,0.0,0.0,0.5,0.5,0.2,2.2,1.4,0.3],
							 "C4":[0.0,0.0,0.0,0.8,0.8,0.2,2.7,1.5,0.4],
							 "C5":[0.0,0.0,0.0,0.0,0.5,0.0,2.6,1.9,0.8],
							 "C6":[0.0,0.2,0.0,1.5,1.1,0.7,3.0,1.5,0.5],
							 "C7":[0.2,0.2,0.0,1.1,1.2,0.4,2.8,1.8,0.8],
							 "C8":[0.0,0.0,0.0,1.5,2.0,0.3,3.1,1.9,0.4],
							 "C9":[0.2,0.2,0.0,1.7,2.0,0.3,3.4,2.7,0.5],
							 "O1":[0.0,0.5,0.2,1.5,0.8,0.5,2.4,3.2,0.8],
							 "N1":[0.0,0.0,0.0,0.3,0.6,0.3,2.3,3.3,0.5]}
		for atomName in self.atomNames:
			newAtom = atom.Atom()
			self.atoms[atomName] = newAtom

	def getHistogram(self, atomName, scaled=False):
		thisAtom = self.atoms[atomName]
		histogram = [thisAtom.counts["nearHydrophobic"],
					 thisAtom.counts["nearPolar"],
					 thisAtom.counts["nearCharged"],
					 thisAtom.counts["medHydrophobic"],
					 thisAtom.counts["medPolar"],
					 thisAtom.counts["medCharged"],
					 thisAtom.counts["farHydrophobic"],
					 thisAtom.counts["farPolar"],
					 thisAtom.counts["farCharged"]]
		if scaled:
			histogram = []
			for shell in ["near","med","far"]:
				countsInShell = thisAtom.counts["%sHydrophobic"%shell] + \
								thisAtom.counts["%sPolar"%shell] + \
								thisAtom.counts["%sCharged"%shell]
				
				for key in ["nearHydrophobic",
							"nearPolar",
							"nearCharged",
							"medHydrophobic",
							"medPolar",
							"medCharged",
							"farHydrophobic",
							"farPolar",
							"farCharged"]:
					if shell in key:
						result = 10000
						try:
							result = float(thisAtom.counts[key])/float(countsInShell)
						except Exception as e: 
							#print(e)
							if countsInShell == 0:
								result = 0
						histogram.append(result)
				
			return histogram
		return histogram
	
	def printTable(self):
		table = ""
		scaledTable = ""
		table += str(self.id)
		scaledTable += str(self.id)
		for atomName in self.atomNames:
			thisAtom = self.atoms[atomName]
			table+= "\n%s: \t\t %i,%i,%i \t\t %i,%i,%i \t\t %i,%i,%i \t\t\t %i" %(thisAtom.name,thisAtom.counts["nearHydrophobic"],thisAtom.counts["nearPolar"],
																				 thisAtom.counts["nearCharged"],thisAtom.counts["medHydrophobic"],thisAtom.counts["medPolar"],
																				 thisAtom.counts["medCharged"],thisAtom.counts["farHydrophobic"],thisAtom.counts["farPolar"],
																				 thisAtom.counts["farCharged"], thisAtom.totalNeighbours())
			#scaledTable+= "\n%s: \t\t %1.1f,%1.1f,%1.1f \t\t %1.1f,%1.1f,%1.1f \t\t %1.1f,%1.1f,%1.1f \t\t\t %i" %(thisAtom.name,thisAtom.nearHydrophobic/thisAtom.totalNearNeighbours(),thisAtom.nearPolar/thisAtom.totalNearNeighbours(),
			#																 thisAtom.nearCharged/thisAtom.totalNearNeighbours(),thisAtom.medHydrophobic/thisAtom.totalMedNeighbours(),thisAtom.medPolar/thisAtom.totalMedNeighbours(),
			#																 thisAtom.medCharged/thisAtom.totalMedNeighbours(),thisAtom.farHydrophobic/thisAtom.totalFarNeighbours(),thisAtom.farPolar/thisAtom.totalFarNeighbours(),
			#																 thisAtom.farCharged/thisAtom.totalFarNeighbours(), thisAtom.totalNeighbours())																	 
		print table
	
	def score(self, scoringWeights):
		#sumOverAtoms(sumOverShell((n-<n>)2/sigma)))
		chi2 = 0
		totalContactWeight = scoringWeights["totalContactWeight"]
		typeOfContactWeight = scoringWeights["typeOfContactWeight"]
		for idx, atom in enumerate(self.atomNames):
			self.totalContacts += self.atoms[atom].totalNeighbours()
			weight = self.weights[idx]
			scaledHistogram = self.getHistogram(atom, True)
			totalHistogram = self.getHistogram(atom, False)
			#overallNumberOfContacts = self.atoms[atom].totalNeighbours()
			contactChi2 = 0
			typeChi2 = 0
			sum = 0
			for i in range (0, 9):
				if self.scaledStdDevs[atom][i] > 0:
					typeChi2 += typeOfContactWeight * weight * (scaledHistogram[i]-self.scaledMeanValues[atom][i])**2/(self.scaledStdDevs[atom][i])**2
				else:
					#stdDev might be zero. Set stddev to 0.1 for such cases.
					typeChi2 += typeOfContactWeight * (scaledHistogram[i]-self.scaledMeanValues[atom][i])**2/0.1
				if self.totalStdDevs[atom][i] > 0:
					contactChi2 += totalContactWeight * (totalHistogram[i]-self.totalMeanValues[atom][i])**2/(self.totalStdDevs[atom][i])**2
				else:
					#stdDev might be zero. Set stddev to 0.1 for such cases.
					contactChi2 += totalContactWeight * (totalHistogram[i]-self.totalMeanValues[atom][i])**2/0.1
					#print "StdDev = 0!"
				
			chi2 += contactChi2 + typeChi2
			#print contactChi2, typeChi2
		self.chi2 = chi2