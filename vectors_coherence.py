# Mauro Alberti  2010-08-24 - http://www.malg.eu; http://www.gistrutturale.it
# GNU General Public License v. 3

import os, sys
from math import *
import numpy as np


kernel_side = 3 # size of the kernel

# returns the magnitude of the vector field as a new numpy array 	
def magnitude(x_grid, y_grid):

	return np.sqrt(x_grid**2 + y_grid**2)

	
def input_error():
	print \
""" 
vectors_coherence_07.py - script for determining vector field coherence
				from Imcorr output
				
Usage: python vectors_coherence_07.py parameter_file     
Example: python vectors_coherence_07.py param.txt

Parameters in parameter file
  input file: vector file
  output_file: output vector file with flow coherence parameters
"""
	sys.exit(1)

	
# extraction of script arguments
try:
  par_fn = sys.argv[1]
except:
	input_error()
	
	
# open input parameter file 
try:
	parameter_file = open(par_fn, 'r')
except:
  print 'Error in opening parameter file'
  sys.exit(1)

  
# reads parameter values
infile_name =  parameter_file.readline().split()[0]
outfile_name =  parameter_file.readline().split()[0]


# open input data file 
try:
	data_file = open(infile_name, 'r')
except:
  print 'Error in opening input data file'
  sys.exit(1)
 

# open output data file 
try:
	output_file = open(outfile_name, 'w')
except:
  print 'Error in creating output file'
  sys.exit(1)

  

# reads header line of data file
header_line = data_file.readline()


# reads data and stores in lists
px_j_list = []; px_i_list = []; x0_list = []; y0_list = []; x1_list = []; y1_list = []
dx_list = []; dy_list = []; x_err_list = []; y_err_list = [] 
magn_displ_list = []; displ_dir_list = []; corrstr_list = []; resflg_list = []

for line in data_file:
	(px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg) = line.split(',')

	px_j_list.append(px_j); px_i_list.append(px_i); x0_list.append(x0); y0_list.append(y0); x1_list.append(x1); y1_list.append(y1)
	dx_list.append(dx); dy_list.append(dy); x_err_list.append(x_err); y_err_list.append(y_err)
	magn_displ_list.append(magn_displ); displ_dir_list.append(displ_dir); corrstr_list.append(corrstr); resflg_list.append(resflg)


# converts lists into numpy arrays
px_j_array = np.array(px_j_list, int)
px_i_array = np.array(px_i_list, int) 
x0_array = np.array(x0_list, float) 
y0_array = np.array(y0_list, float) 
x1_array = np.array(x1_list, float) 
y1_array = np.array(y1_list, float)
dx_array = np.array(dx_list, float) 
dy_array = np.array(dy_list, float) 
x_err_array = np.array(x_err_list, float) 
y_err_array = np.array(y_err_list, float) 
magn_displ_array = np.array(magn_displ_list, float) 
displ_dir_array = np.array(displ_dir_list, float) 
corrstr_array = np.array(corrstr_list, float) 
resflg_array = np.array(resflg_list, int)

	
# min and max of px_j and px_i
min_px_j, max_px_j = min(px_j_array), max(px_j_array)
min_px_i, max_px_i = min(px_i_array), max(px_i_array)

# min and max of x0 and y0
min_x0, max_x0 = min(x0_array), max(x0_array)
min_y0, max_y0 = min(y0_array), max(y0_array)


delta_px_j_array = px_j_array[1:] - px_j_array[:-1]
delta_px_j_array = np.where(delta_px_j_array <= 0, 999999, delta_px_j_array)
incr_px_j = min(delta_px_j_array)
if incr_px_j == 999999:
  print 'Error in grid geometry processing (px_j values)'
  sys.exit(1)	


delta_px_i_array = px_i_array[1:] - px_i_array[:-1]
delta_px_i_array = np.where(delta_px_i_array <= 0, 999999, delta_px_i_array)
incr_px_i = min(delta_px_i_array)
if incr_px_i == 999999:
  print 'Error in grid geometry processing (px_i values)'
  sys.exit(1)
  
  
# definition of the output array
out_grid_colums = ((max_px_j - min_px_j)/incr_px_j) + 1
out_grid_lines = ((max_px_i - min_px_i)/incr_px_i) + 1
  
  
# reshape arrays
px_j_array = px_j_array.reshape(out_grid_lines,out_grid_colums)
px_i_array = px_i_array.reshape(out_grid_lines,out_grid_colums) 
x0_array = x0_array.reshape(out_grid_lines,out_grid_colums) 
y0_array = y0_array.reshape(out_grid_lines,out_grid_colums) 
x1_array = x1_array.reshape(out_grid_lines,out_grid_colums) 
y1_array = y1_array.reshape(out_grid_lines,out_grid_colums)
dx_array = dx_array.reshape(out_grid_lines,out_grid_colums) 
dy_array = dy_array.reshape(out_grid_lines,out_grid_colums) 
x_err_array = x_err_array.reshape(out_grid_lines,out_grid_colums) 
y_err_array = y_err_array.reshape(out_grid_lines,out_grid_colums) 
magn_displ_array = magn_displ_array.reshape(out_grid_lines,out_grid_colums) 
displ_dir_array = displ_dir_array.reshape(out_grid_lines,out_grid_colums) 
corrstr_array = corrstr_array.reshape(out_grid_lines,out_grid_colums) 
resflg_array = resflg_array.reshape(out_grid_lines,out_grid_colums)


# vectorizing
dx_00 = dx_array[:-2,:-2]*corrstr_array[:-2,:-2]
dx_10 = dx_array[1:-1,:-2]*corrstr_array[1:-1,:-2]
dx_20 = dx_array[2:,:-2]*corrstr_array[2:,:-2]

dx_01 = dx_array[:-2,1:-1]*corrstr_array[:-2,1:-1]
#dx_11 = dx_array[1:-1,1:-1]*corrstr_array[1:-1,1:-1]
dx_21 = dx_array[2:,1:-1]*corrstr_array[2:,1:-1]

dx_02 = dx_array[:-2,2:]*corrstr_array[:-2,2:]
dx_12 = dx_array[1:-1,2:]*corrstr_array[1:-1,2:]
dx_22 = dx_array[2:,2:]*corrstr_array[2:,2:]

dy_00 = dy_array[:-2,:-2]*corrstr_array[:-2,:-2]
dy_10 = dy_array[1:-1,:-2]*corrstr_array[1:-1,:-2]
dy_20 = dy_array[2:,:-2]*corrstr_array[2:,:-2]

dy_01 = dy_array[:-2,1:-1]*corrstr_array[:-2,1:-1]
#dy_11 = dy_array[1:-1,1:-1]*corrstr_array[1:-1,1:-1]
dy_21 = dy_array[2:,1:-1]*corrstr_array[2:,1:-1]

dy_02 = dy_array[:-2,2:]*corrstr_array[:-2,2:]
dy_12 = dy_array[1:-1,2:]*corrstr_array[1:-1,2:]
dy_22 = dy_array[2:,2:]*corrstr_array[2:,2:]

dx_sum = dx_00 + dx_10 + dx_20 + dx_01 +dx_21 + dx_02 + dx_12 + dx_22
dy_sum = dy_00 + dy_10 + dy_20 + dy_01 +dy_21 + dy_02 + dy_12 + dy_22


# sum of weights for border cells
corrstr_correct = np.where(magn_displ_array > 0.0, corrstr_array, 0.0)
corrstr_sum = corrstr_correct[:-2,:-2] + corrstr_correct[1:-1,:-2] + corrstr_correct[2:,:-2] + \
	corrstr_correct[:-2,1:-1] + corrstr_correct[2:,1:-1] + \
	corrstr_correct[:-2,2:] + corrstr_correct[1:-1,2:] + corrstr_correct[2:,2:]
corrstr_sum = np.where(corrstr_sum > 0.0, corrstr_sum, np.NaN)	


# components of mean values of weighted sum of border displacements	
dx_wmean = dx_sum/corrstr_sum
dy_wmean = dy_sum/corrstr_sum
magnitude_wmean = np.sqrt(dx_wmean**2 + dy_wmean**2)

# magnitude difference
kernbord_magnitude = (np.sqrt(dx_00**2 + dy_00**2)+np.sqrt(dx_10**2 + dy_10**2)+np.sqrt(dx_20**2 + dy_20**2)+ \
np.sqrt(dx_01**2 + dy_01**2)+np.sqrt(dx_21**2 + dy_21**2)+ \
np.sqrt(dx_02**2 + dy_02**2)+np.sqrt(dx_12**2 + dy_12**2)+np.sqrt(dx_22**2 + dy_22**2))/corrstr_sum

kerncent_magnitude = np.where(magn_displ_array > 0.0, magn_displ_array, np.NaN)

magnitude_diff = np.zeros((out_grid_lines,out_grid_colums))*np.NaN
magnitude_diff[1:-1, 1:-1] = abs(kerncent_magnitude[1:-1,1:-1] - kernbord_magnitude)


# orientation difference
scalar_product = dx_array[1:-1,1:-1]*dx_wmean + dy_array[1:-1,1:-1]*dy_wmean
acos_angdiff = scalar_product/(kerncent_magnitude[1:-1,1:-1]*magnitude_wmean)

acos_angdiff = np.where(acos_angdiff < -1.0, -1.0, acos_angdiff)
acos_angdiff = np.where(acos_angdiff > 1.0, 1.0, acos_angdiff)

ang_diff = np.zeros((out_grid_lines,out_grid_colums))*np.NaN
ang_diff[1:-1, 1:-1] = np.arccos(acos_angdiff)*180.0/pi


# creates output string
output_string = 'px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg,magn_diff,ang_diff\n' # header line

for b in range(out_grid_colums): 
	for a in range(out_grid_lines):
		if not(isnan(kerncent_magnitude[a,b])):
			output_string = output_string + (str(px_j_array[a,b])+','+str(px_i_array[a,b])+','+str(x0_array[a,b])+','+str(y0_array[a,b])+','+ \
				str(x1_array[a,b])+','+str(y1_array[a,b])+','+str(dx_array[a,b])+','+str(dy_array[a,b])+','+str(x_err_array[a,b])+','+ \
				str(y_err_array[a,b])+','+str(magn_displ_array[a,b])+','+str(displ_dir_array[a,b])+','+str(corrstr_array[a,b])+','+ \
				str(resflg_array[a,b])+','+str(magnitude_diff[a,b])+','+str(ang_diff[a,b])+'\n')

output_string = output_string.replace(",nan",",-999")
		
		
# writes results in output file
output_file.write(output_string)






