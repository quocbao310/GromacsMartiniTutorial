import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance
import subprocess
import argparse
import os

coords = []

def onclick(event):

	global coords
	
	print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata))
	
	coords.append((float(format(event.xdata, ".5f")), float(format(event.ydata, ".5f"))))

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = 'Interacting with the principle component map to extract favorable states')
	parser.add_argument("--input", type=validate_file, help='Input file (2dproj.xvg)')
	parser.add_argument('--traj', help='Trajectory (traj.xtc)')
	parser.add_argument('--tpr', help='Topology (topol.tpr)')
	parser.add_argument('--dt', type=int, help='Time step (ps)')
	

	args = parser.parse_args()
	
	PCA = []
	
	with open(args.input, "r") as f:
		for line in f.readlines():
			if not (line.startswith("#") or line.startswith("@")):
				PCA.append(line)	
				
	PCA = [i.strip().split() for i in PCA]
	
	PC1 = [round(float(i[0]), 5) for i in PCA]
	PC2 = [round(float(i[1]), 5) for i in PCA]
	
	# Plotting PCA
	fig = plt.figure()
	plt.plot(PC1,PC2,c='k')
	plt.xlabel('PC1')
	plt.ylabel('PC2')
	cid = fig.canvas.mpl_connect('button_press_event', onclick)

	plt.show()
	

	# Start extracting frame
	commands = []
	
	for i in coords:
		distances = [distance.euclidean(i, [pc1, pc2]) for pc1, pc2 in zip(PC1, PC2)]	
		index = distances.index(min(distances))
		print(index)
		print(PC1[index], PC2[index])
		
		p = subprocess.Popen(['gmx_2021.2', 'trjconv', 
					'-f', args.traj, '-s', args.tpr,
					'-dump', str(int(index) * args.dt), '-o', 'frame_' + str(PC1[index]) + '_' + str(PC2[index]) + '_' + str(int(index) * args.dt) + '.pdb'],
					stdin=subprocess.PIPE)		
		p.communicate(b'1\n')
		p.wait()

		commands.append('Command just used: ' + 'gmx_2021.2 trjconv -f ' + args.traj + ' -s '+ args.tpr + ' -dump ' + str(int(index) * 10) + ' -o frame_' + str(PC1[index]) + '_' + str(PC2[index]) + '_' + str(int(index) * 10) + '.pdb')		
		

	for i in commands:
		print(i)

