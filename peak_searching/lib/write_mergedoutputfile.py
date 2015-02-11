# -*- coding: utf-8 -*-

# write .txt output file for merged damage sites from multiple
# damage datasets once filtering has been performed


def mergedoutput_writer(where,damage_sites,num_damsets):

	outputfile = open(where+'mergedoutputfile.txt','w')

	outputfile.write('atomtype'+'\t'+'resinum'+'\t'+'basetype'+'\t'+'chaintype' +'\t')
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('X_%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('Y_%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('Z_%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('atomnum_%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('peaknum_%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('mindist__%s\t' %(str(i+1)))
	for i in range(0,len(damage_sites[0].peaknum)):
		outputfile.write('peakvalue__%s\t' %(str(i+1)))
	outputfile.write('\n')

	for dsite in damage_sites:
		outputfile.write(str(dsite.atomtype)+'\t'+str(dsite.residuenum)+\
				         '\t'+str(dsite.basetype)+'\t'+str(dsite.chaintype)+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.X_coord[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.Y_coord[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.Z_coord[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.atomnum[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.peaknum[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.mindist[i])+'\t')
		for i in range(0,len(damage_sites[0].peaknum)):
			outputfile.write(str(dsite.peakvalue[i])+'\t')
		outputfile.write('\n')

 