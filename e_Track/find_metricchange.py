# -*- coding: utf-8 -*-

def find_Bfactorchange(initialPDB,laterPDB):
	# function to determine the bfactor change between the initial and
	# and later dataset --> becomes an object attribute for the later 
	# dataset

	print '------------------------------------------------------------'
	print 'Determining Bfactor change between initial and later dataset'
	for atom in initialPDB:
		for otheratom in laterPDB:
			if (atom.atomtype == otheratom.atomtype and
			   atom.basetype == otheratom.basetype and
			   atom.chaintype == otheratom.chaintype and
			   atom.residuenum == otheratom.residuenum):

				otheratom.Bfactorchange = float(otheratom.Bfactor) - float(atom.Bfactor)
	print '---> success...'


def find_Bdamagechange(initialPDB,laterPDB):
	# function to determine the bdamage change between the initial and
	# and later dataset --> becomes an object attribute for the later 
	# dataset

	print '------------------------------------------------------------'
	print 'Determining Bdamage change between initial and later dataset'
	for atom in initialPDB:
		for otheratom in laterPDB:
			if (atom.atomtype == otheratom.atomtype and
			   atom.basetype == otheratom.basetype and
			   atom.chaintype == otheratom.chaintype and
			   atom.residuenum == otheratom.residuenum):

				otheratom.bdamchange = float(otheratom.bdam) - float(atom.bdam)
	print '---> success...'


