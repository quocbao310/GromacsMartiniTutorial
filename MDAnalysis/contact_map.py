from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis
from scipy.spatial import distance
from matplotlib import cm
from PIL import Image
import os, glob
import re
import argparse
import cv2

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f

def plotting_interface_matrix(matrix, i, LR_res_list, Leptin_res_list):
	
	fig, ax= plt.subplots(figsize=(10,10))
	im = ax.imshow(matrix, interpolation='nearest', cmap=cm.gray, 
			vmin=0.0, vmax=100.0)

	ax.set_title('Average pair distance [nm]')
	ax.set_ylabel('Leptin Receptor')
	ax.set_xlabel('Leptin')
	ax.invert_yaxis()
	# Move left and bottom spines outward by 10 points
	ax.spines['left'].set_position(('outward', 10))
	ax.spines['bottom'].set_position(('outward', 10))
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	ax.set_xticks(range(len(Leptin_res_list)))
	ax.set_yticks(range(len(LR_res_list)))
	ax.set_xticklabels(Leptin_res_list)
	ax.set_yticklabels(LR_res_list)	
	
	plt.setp(ax.get_xticklabels(), rotation=60, ha='right', rotation_mode="anchor")
	plt.colorbar(im, fraction=0.046, pad=0.04)
	fig.tight_layout()	
	plt.savefig('frame_{}.png'.format(i)) # Save frame


def plotting_matrix(matrix, i):
	
	fig, ax= plt.subplots(figsize=(10,10))
	im = ax.imshow(matrix, interpolation='nearest', cmap=cm.gray, 
			vmin=0.0, vmax=100.0)

	ax.set_title('Average pair distance [nm]')
	ax.set_ylabel('Leptin Receptor')
	ax.set_xlabel('Leptin')
	ax.invert_yaxis()
	# Move left and bottom spines outward by 10 points
	ax.spines['left'].set_position(('outward', 10))
	ax.spines['bottom'].set_position(('outward', 10))
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	
	plt.setp(ax.get_xticklabels(), rotation=60, ha='right', rotation_mode="anchor")
	plt.colorbar(im)
	fig.tight_layout()	
	plt.savefig('frame_{}.png'.format(i)) # Save frame

def plotting_matrix_chain(matrix, i, chainID):
	
	fig, ax= plt.subplots(figsize=(10,10))
	im = ax.imshow(matrix, interpolation='nearest', cmap=cm.gray, 
			vmin=0.0, vmax=100.0)

	ax.set_title('Average pair distance [nm]')
	ax.set_ylabel(chainID)
	ax.set_xlabel(chainID)
	ax.invert_yaxis()
	# Move left and bottom spines outward by 10 points
	ax.spines['left'].set_position(('outward', 10))
	ax.spines['bottom'].set_position(('outward', 10))
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	
	
	plt.setp(ax.get_xticklabels(), rotation=60, ha='right', rotation_mode="anchor")
	plt.colorbar(im)
	fig.tight_layout()	
	plt.savefig('frame_{}.png'.format(i)) # Save frame

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)
	parser.add_argument("-interface", type=validate_file, help="File containing interface residues", required=False)
	parser.add_argument("-chain", help="To create protein contact map of specified chain", required=False)
	args = parser.parse_args()
	

	################################################################
	################## Trajectory and topology #####################
	################################################################	
	
	XTC = args.traj
	PDB = args.tpr
	u = MDAnalysis.Universe(PDB, XTC)
	
	
	if bool(args.interface):
		f = open(args.interface, "r").readlines()	
		
		LR_residues = [res.strip() for res in f if res[1] == 'A']
		LR_residues = [s for s in re.findall(r'\d+', LR_residues[0])]
		
		Leptin_residues = [res.strip() for res in f if res[1] == 'B']
		Leptin_residues = [s for s in re.findall(r'\d+', Leptin_residues[0])]

		
		for ts in u.trajectory:
			if ts.frame == 0 or (ts.frame % 100 == 0):
				
				LR_res_pattern = ''
				for resid in LR_residues:
					LR_res_pattern += 'resid ' + resid + ' or '
				LR_res_pattern  = LR_res_pattern[:-3]
				LR_BB = u.select_atoms('name BB and ('+ LR_res_pattern +') and segid A')
				
				
				Leptin_res_pattern = ''
				for resid in Leptin_residues:
					Leptin_res_pattern += 'resid ' + resid + ' or '
				Leptin_res_pattern  = Leptin_res_pattern[:-3]
				Leptin_BB = u.select_atoms('name BB and ('+ Leptin_res_pattern +') and segid B')				
				
				matrix = np.zeros((len(LR_BB.positions), len(Leptin_BB.positions)))
				
				##############################################
				####### Calculating pairwise distance ########
				##############################################

				for i in range(len(LR_BB.positions)):
					for j in range(len(Leptin_BB.positions)):
						distances = distance.euclidean(LR_BB.positions[i],
										Leptin_BB.positions[j])

						matrix[i][j] = distances

				flip_matrix = np.flip(matrix,0)
				plotting_interface_matrix(flip_matrix, ts.frame, LR_residues, Leptin_residues)

	elif bool(args.chain):
		
		for ts in u.trajectory:
	
			if (ts.frame == 0) or (ts.frame % 100 == 0):

				row = u.select_atoms('name BB and segid ' + str(args.chain))
				column =  u.select_atoms('name BB and segid ' + str(args.chain))
				
				matrix = np.zeros((len(row.positions), len(column.positions)))
				
				##############################################
				####### Calculating pairwise distance ########
				##############################################
				for i in range(len(row.positions)):
					for j in range(len(column.positions)):
						distances = distance.euclidean(row.positions[i],
										column.positions[j])
										
						matrix[i][j] = distances	
						
				flip_matrix = np.flip(matrix,0)
				plotting_matrix_chain(flip_matrix, ts.frame, 'Chain ' + str(args.chain))			
		
	else:
		for ts in u.trajectory:
	
			if (ts.frame == 0) or (ts.frame % 100 == 0):
		
				LR_BB = u.select_atoms('name BB and segid A')
				
				Leptin_BB = u.select_atoms('name BB and segid B')
				
				matrix = np.zeros((len(LR_BB.positions), len(Leptin_BB.positions)))
			
				##############################################
				####### Calculating pairwise distance ########
				##############################################
				for i in range(len(LR_BB.positions)):
					for j in range(len(Leptin_BB.positions)):
						distances = distance.euclidean(LR_BB.positions[i],
										Leptin_BB.positions[j])

						matrix[i][j] = distances
			
				flip_matrix = np.flip(matrix,0)
				plotting_matrix(flip_matrix, ts.frame)				


	# Create the frames
	frames = []
	imgs = glob.glob('frame_*.png')
	imgs.sort(key=lambda f: int(re.sub('\D', '', f)))
	for i in imgs:
		new_frame = Image.open(i)
		frames.append(new_frame)

	
	# Create a movie
	img_array = []
	for filename in imgs:
		img = cv2.imread(filename)
		height, width, layers = img.shape
		size = (width,height)
		img_array.append(img)
	
	out = cv2.VideoWriter('project.avi',cv2.VideoWriter_fourcc(*'DIVX'), 1, size)

	for i in range(len(img_array)):
		out.write(img_array[i])
	out.release()
