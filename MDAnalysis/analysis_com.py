from __future__ import print_function
import numpy as np
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import MDAnalysis
from scipy.spatial import distance
from matplotlib import cm
from PIL import Image
import os, glob
import re
import argparse

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): tpr gro pdb", required=True)
	args = parser.parse_args()	

	u = MDAnalysis.Universe(args.tpr, args.traj)
	
	distance_com = []
	frame = []
	
	for ts in u.trajectory:
		
		A = u.select_atoms("segid A and name BB").center_of_geometry()
		B = u.select_atoms("segid B and name BB").center_of_geometry()
		distances = distance.euclidean(A, B)
		distance_com.append(distances)
		frame.append(ts.frame)		

	plt.plot(frame, distance_com)  # Plot the chart
	plt.savefig('COM_MoleculeDistance.png')

