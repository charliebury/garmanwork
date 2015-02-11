# -*- coding: utf-8 -*-

# function here is designed to read in a list of objects merged over multiple datasets
# (PDBmulti) and write a txt file containing the density metrics for each atom in the 
# structure, over the different datasets
def objlist2txt(PDBmulti,where):

	# make a new txt file in location where
	txtfile = open(where+'PDBmulti.txt','w')
	txtfile.write('atomnum\t')
	txtfile.write('atominfo\t')
	for i in range(0,len(PDBmulti[0].mindensity)):
		txtfile.write('mindens_%s\t' %(str(i+1)))
	for i in range(0,len(PDBmulti[0].mindensity)):
		txtfile.write('maxdens_%s\t' %(str(i+1)))
	for i in range(0,len(PDBmulti[0].mindensity)):
		txtfile.write('meandens_%s\t' %(str(i+1)))
	for i in range(0,len(PDBmulti[0].mindensity)):
		txtfile.write('mediandens_%s\t' %(str(i+1)))
	txtfile.write('\n')

	# sort by atomnumber
	PDBmulti.sort(key=lambda x: x.atomnum)

	# for each atom in PDBmulti
	for atom in PDBmulti:
		txtfile.write('%s\t' %(str(atom.atomnum)))
		txtfile.write('%s%s%s%s\t' %(str(atom.atomtype),str(atom.residuenum),str(atom.basetype),str(atom.chaintype)))
		for i in range(0,len(PDBmulti[0].mindensity)):
			txtfile.write('%s\t' %(str(atom.mindensity[i])))
		for i in range(0,len(PDBmulti[0].mindensity)):
			txtfile.write('%s\t' %(str(atom.maxdensity[i])))
		for i in range(0,len(PDBmulti[0].mindensity)):
			txtfile.write('%s\t' %(str(atom.meandensity[i])))
		for i in range(0,len(PDBmulti[0].mindensity)):
			txtfile.write('%s\t' %(str(atom.mediandensity[i])))
		txtfile.write('\n')

	txtfile.close()


