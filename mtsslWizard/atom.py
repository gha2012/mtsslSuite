import numpy

class Atom:
	def __init__(self):
		self.name = ""
		self.resn = ""
		self.resi = ""
		self.element = ""
		self.coordinate = numpy.zeros(3)
		self.isSpinCenter = False
		self.indices = {}
		self.counts = {"nearHydrophobic":0, "nearPolar":0,"nearCharged":0,
					"medHydrophobic":0,"medPolar":0,"medCharged":0,
					"farHydrophobic":0,"farPolar":0,"farCharged":0}

	def totalNeighbours(self):
		total = self.counts["nearHydrophobic"] + \
				self.counts["nearPolar"] + \
				self.counts["nearCharged"] + \
				self.counts["medHydrophobic"] + \
				self.counts["medPolar"] + \
				self.counts["medCharged"] + \
				self.counts["farHydrophobic"] + \
				self.counts["farPolar"] + \
				self.counts["farCharged"]
		return total
		
	def totalNearNeighbours(self):
		total = self.counts["nearHydrophobic"] + \
				self.counts["nearPolar"] + \
				self.counts["nearCharged"]
		return float(total)

	def totalMedNeighbours(self):
		total = self.counts["medHydrophobic"] + \
				self.counts["medPolar"] + \
				self.counts["medCharged"]
		return float(total)

	def totalFarNeighbours(self):
		total = self.counts["medHydrophobic"] + \
				self.counts["medPolar"] + \
				self.counts["medCharged"]
		return float(total)
		
	def isCharged(self, atomName, resn):
		if atomName in ["OD1","OD2","OE1","OE2","NH1","NH2","NE", "ND1","NE2","NZ"] and resn in ["ASP","GLU","HIS","LYS","ARG"] or atomName == "OXT":
			return True
		else:
			return False

	def isPolar(self, atomName, resn):
		if atomName in ["O", "N", "OG","OG1", "OD", "OE", "ND", "NE", "NZ", "NH","NH1","NH2","OD1","OD2","ND1","ND2","OE1","OE2","NE1","NE2","SG","OG2", "O1", "N1", "O","OH"] and not self.isCharged(atomName, resn):
			return True
		else:
			return False

	def ishydrophobic(self, atomName, resn):
		if atomName in ["CA", "CB", "CG", "CD", "CE", "C","CG1","CG2","CD1","CD2", "SD","CE1","CE2","CE3","CZ","CZ1","CZ2","CZ3","CH1","CH2", "C3", "C2", "C4", "C5", "C6", "C7", "C8", "C9"]:
			return True
		else:
			return False
