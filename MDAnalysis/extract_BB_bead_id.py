import sys, os, re

if __name__=="__main__":

	file_input = sys.argv[1]
	
	with open(file_input, 'r') as f:
		file_input_contents = f.readlines()
	
	list_of_BB_index = []
		
	for line in file_input_contents:
		if (line[13:15] == 'BB'):
			list_of_BB_index.append(line[7:11])
	
	print(len(list_of_BB_index))		
	# Write to an index file	
	out_file = open('index.ndx', 'a')
	out_file.write('[ BB ]\n')
		
	for i in range(0, len(list_of_BB_index), 15):
		for idex in list_of_BB_index[i:i+15]:
			out_file.write(idex)
		out_file.write('\n')
			
