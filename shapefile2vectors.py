# Mauro Alberti  2010-08-24 - http://www.malg.eu; http://www.gistrutturale.it
# GNU General Public License v. 3

import os, sys
from math import *
import numpy as np

try:
	from osgeo import ogr
except:
	import ogr
	
		
def input_error():

	print \
""" 
shapefile2vectors_05.py - script for converting line shapefile 
						to vector grid 

Usage: python shapefile2vectors_05.py input_shapefile out_vectorgrid     
Example: python shapefile2vectors_05.py vectors.shp vectors.txt

"""
	sys.exit(1)


	
# extraction of script arguments
try:
  input_shape_name = sys.argv[1]
  output_vectorgrid_name = sys.argv[2]
except:
	input_error()


# open input shapefile

driver = ogr.GetDriverByName('ESRI Shapefile') 
inshape = driver.Open(input_shape_name)
if inshape is None:
  print 'error in reading input shapefile'
  sys.exit(1)
  
	
# open output vector grid file 
try:
	outfile = open(output_vectorgrid_name, 'w')
except:
  print 'error in creating output vector grid file'
  sys.exit(1)  


# reads lines from input shapefile

px_j_list = []; px_i_list = []
x0_list = []; y0_list = []
x1_list = []; y1_list = []
xerr_list = []; yerr_list = []
corrstr_list = []; resflg_list = []


layer = inshape.GetLayer()
feature = layer.GetNextFeature()
while feature:

	line = feature.GetGeometryRef()	

	px_j_list.append(int(feature.GetField('px_j')))
	px_i_list.append(int(feature.GetField('px_i')))
	
	x0_list.append(float(line.GetX(0)))
	y0_list.append(float(line.GetY(0)))
	x1_list.append(float(line.GetX(1)))
	y1_list.append(float(line.GetY(1)))

	xerr_list.append(float(feature.GetField('x_err')))
	yerr_list.append(float(feature.GetField('y_err')))	

	corrstr_list.append(float(feature.GetField('corrstr')))
	resflg_list.append(int(feature.GetField('resflg')))
			
	feature.Destroy()
	feature = layer.GetNextFeature()

inshape.Destroy()


# number of input records
num_input_recs = len(px_j_list)	

# conversion from lists to arrays
px_j_array = np.asarray(px_j_list, 'i')
px_i_array = np.asarray(px_i_list, 'i')

x0_array = np.asarray(x0_list, 'f')
y0_array = np.asarray(y0_list, 'f')
x1_array = np.asarray(x1_list, 'f')
y1_array = np.asarray(y1_list, 'f')

dx_array = x1_array - x0_array
dy_array = y1_array - y0_array

magn_displ_array = np.sqrt(dx_array**2 + dy_array**2)

displ_dir_array =  np.arctan2(dx_array,dy_array)*180.0/pi
displ_dir_array = np.where((displ_dir_array < 0.0), displ_dir_array + 360.0, displ_dir_array)
displ_dir_array = np.where((displ_dir_array > 360.0), displ_dir_array - 360.0, displ_dir_array)

xerr_array = np.asarray(xerr_list, 'f')
yerr_array = np.asarray(yerr_list, 'f')

corrstr_array = np.asarray(corrstr_list, 'f')
resflg_array = np.asarray(resflg_list, 'i')

	
# min and max of px_j and px_i
min_px_j, max_px_j = min(px_j_array), max(px_j_array)
min_px_i, max_px_i = min(px_i_array), max(px_i_array)

# min and max of x0 and y0
min_x0, max_x0 = min(x0_array), max(x0_array)
min_y0, max_y0 = min(y0_array), max(y0_array)

pixel_size_x = (max_x0-min_x0)/(max_px_j-min_px_j)
pixel_size_y = (max_y0-min_y0)/(max_px_i-min_px_i)


delta_px_j_array = abs(px_j_array[1:] - px_j_array[:-1])
delta_px_j_array = np.where(delta_px_j_array == 0, 999999, delta_px_j_array)
incr_px_j = min(delta_px_j_array)

delta_px_i_array = abs(px_i_array[1:] - px_i_array[:-1])
delta_px_i_array = np.where(delta_px_i_array == 0, 999999, delta_px_i_array)
incr_px_i = min(delta_px_i_array)


# definition of the output array
out_grid_colums = ((max_px_j - min_px_j)/incr_px_j) + 1
out_grid_lines = ((max_px_i - min_px_i)/incr_px_i) + 1


# header of output file
out_string_array = []

# writes null values for all grid records
for a in range(out_grid_lines):
	for b in range(out_grid_colums):
		px_j = min_px_j + incr_px_j*b
		px_i = min_px_i + incr_px_i*a
		x0 = min_x0 + b*incr_px_j*pixel_size_x
		y0 = min_y0 + a*incr_px_i*pixel_size_y
		x1 = x0
		y1 = y0
		out_string = '%d,%d,%f,%f,%f,%f,0.0,0.0,0.0,0.0,0.0,-9999,0.0,3\n' % (px_j,px_i,x0,y0,x1,y1)
		out_string_array.append(out_string)

# substitutes valid values 
for n in range(num_input_recs):

	out_string = '%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n' % (px_j_array[n],px_i_array[n],x0_array[n],y0_array[n],x1_array[n],y1_array[n],dx_array[n],dy_array[n],xerr_array[n],yerr_array[n],magn_displ_array[n],displ_dir_array[n],corrstr_array[n],resflg_array[n])
	
	position_in_outmatrix = ((px_i_array[n] - min_px_i)/incr_px_i)*out_grid_colums + ((px_j_array[n] - min_px_j)/incr_px_j) 
	
	out_string_array[position_in_outmatrix] = out_string
	
	
# writes and closes output file
outfile.write('px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg\n')
for out_string in out_string_array:
	outfile.write(out_string)
outfile.close()


	
	
	
		
		
	
	





