from operator import attrgetter

class Ensemble:
	def __init__(self):
		self.name = ""
		self.rotamers = []
	
	def calculateScaledChi2(self):
		sortedRotamers = self.sortRotamers("chi2")
		lowestChi2 = sortedRotamers[0].chi2
		highestChi2 = sortedRotamers[len(sortedRotamers)-1].chi2
		for rotamer in self.rotamers:
			rotamer.scaledChi2 = (rotamer.chi2 - lowestChi2)/(highestChi2 - lowestChi2)
	
	def sortRotamers(self, key, reverse = False):
		#sort the rotamers according to key
		sortedRotamers = []
		if not reverse:
			sortedRotamers = sorted(self.rotamers, key=attrgetter(key))
		else:
			sortedRotamers = sorted(self.rotamers, key=attrgetter(key), reverse = True)
		#self.rotamers = sortedRotamers
		first = getattr(sortedRotamers[0], key)
		last = getattr(sortedRotamers[len(sortedRotamers)-1], key)
		avg = 0
		for rotamer in sortedRotamers:
			avg += getattr(rotamer, key)
		avg/=len(sortedRotamers)
		#print "First %s: %1.2f" %(key, first)
		#print "Last  %s: %1.2f" %(key, last)
		#print "Avg.  %s: %1.2f" %(key, avg)
		return sortedRotamers

	def writePDB(self, filename = "output.pdb"):
		structure = BPStructure(1)
		for idx, rotamer in enumerate(self.rotamers):
			model = BPModel(idx + 1)
			chain = BPChain("A")
			residue = BPResidue((" ", 1, " "), "R1A", "1")
			for key in rotamer.atoms:
				atomName = key
				#print key
				element = rotamer.atoms[key].element
				x = float(rotamer.atoms[key].coordinate[0])
				y = float(rotamer.atoms[key].coordinate[1])
				z = float(rotamer.atoms[key].coordinate[2])
				atom = BPAtom(atomName, (x,y,z), 10.0, 1.0, " ", " %s "%atomName, 1, "%s"%element)
				residue.add(atom)
			chain.add(residue)
			model.add(chain)
			structure.add(model)
		io = PDBIO()
		io.set_structure(structure)
		io.save("output.tmp")
		
		#tweak output file so that it can be read as multiple states by Pymol
		bad_words = ["END","TER"]
		good_words = ["ENDMDL"]
		with open("output.tmp") as oldfile, open(filename, 'w') as newfile:
			for line in oldfile:
				if not any(bad_word in line for bad_word in bad_words) or any (good_word in line for good_word in good_words):
					newfile.write(line)
		if len(self.rotamers) > 1:
			with open(filename, 'r+') as f:
				content = f.read()
				f.seek(0, 0)
				f.write("NUMMDL    %i\n" %len(self.rotamers) + content)
		os.remove("output.tmp")
		print "Written to file: %s/%s" %(os.getcwd(), filename)