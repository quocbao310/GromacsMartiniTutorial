import os, glob
import numpy as np
from pymol import cmd
from pymol import stored
from InterfaceResidues import interfaceResidues
import argparse

d3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
	'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 
	'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL':'V', 'TRP': 'W', 'TYR': 'Y'}

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f

def DetectingInterfaceResidues(chainA, chainB, PDBFile):

	
	InterfaceResA = []
	InterfaceResB = []
	
	cmd.load(PDBFile)

	
	interface = interfaceResidues(os.path.basename(PDBFile).split('.')[0], 
					cA='c. ' + chainA, 
					cB='c. ' + chainB, 
					cutoff=1.0,		
					selName="interface")
	
			
	for residue in interface:
		stored.list = []
		if residue[0] == 'chA':
			cmd.iterate('resi ' + str(residue[1]) + ' and chain ' + chainA,
				    'stored.list.append(resn)')
				    
			InterfaceResA.append([residue[1] , np.unique(stored.list)[0], 
				      		chainA])
			    			  
				    
		elif residue[0] == 'chB':		    
			cmd.iterate('resi ' + str(residue[1]) + ' and chain ' + chainB,
				    'stored.list.append(resn)')				    
			InterfaceResB.append([residue[1] , np.unique(stored.list)[0], 
				      		chainB])					
	
	InterfaceResA = [d3to1[element[1]] + element[2] + element[0] for element in InterfaceResA]	  
	InterfaceResB = [d3to1[element[1]] + element[2] + element[0] for element in InterfaceResB]
	return InterfaceResA, InterfaceResB

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument("-chain1", help="Name of the first chain", required=True)
	parser.add_argument("-chain2", help="Name of the second chain", required=True)
	parser.add_argument("-PDBFile", type=validate_file, help="Complex PDB file", required=True)
	args = parser.parse_args()
	
	
	InterfaceA, InterfaceB = DetectingInterfaceResidues(args.chain1, args.chain2, args.PDBFile)


	with open('InterfaceRes.txt', 'w') as f:
		f.write(','.join(InterfaceA) + '\n')
		f.write(','.join(InterfaceB) + '\n')
		
