# Mauro Alberti  2010-08-24 - http://www.malg.eu; http://www.gistrutturale.it
# GNU General Public License v. 3


import os, sys
from math import *
import numpy as np

		
def input_error():
	print \
""" 
imcorr2vectors_01.py - script for converting imcorr output
	           to a linear shapefile

Usage: python imcorr2vectors_01.py parameter_file     
Example: python imcorr2vectors_01.py param.txt

Parameters in parameter file
  imcorr_output: text file created from Imcorr
  cell_size: cell size of analysed raster images
  x_min: geographic x value of bottom-left
  y_max: geographic y value of top-right
  out_file: output displacement vectors 
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
  print 'error in reading parameter file'
  sys.exit(1)

  
# reads parameter values
imcorr_fn =  parameter_file.readline().split()[0]
cell_size =  float(parameter_file.readline().split()[0])
x_min =  float(parameter_file.readline().split()[0])
y_max =  float(parameter_file.readline().split()[0])
outgeo_fnm =  parameter_file.readline().split()[0]

  
# open input imcorr file 
try:
	imcorr_file = open(imcorr_fn, 'r')
except:
  print 'error in reading input Imcorr file'
  sys.exit(1)
  
  
# open output file 
try:
	outfile = open(outgeo_fnm, 'w')
except:
  print 'error in creating output file'
  sys.exit(1)
  
  
# reads data from imcorr file
px_j_list = []; px_i_list = []; totdispl_px_list = []; corrstr_list = []; resflg_list = []
xdispl_px_list = []; ydispl_px_list = []; xerr_px_list = []; yerr_px_list = []
for line in imcorr_file:
	(px_j, px_i, totdispl_px, corrstr, resflg, xdispl_px, ydispl_px, xerr_px, yerr_px) = line.split()
	
	px_j_list.append(px_j); px_i_list.append(px_i); totdispl_px_list.append(totdispl_px)
	corrstr_list.append(corrstr); resflg_list.append(resflg); xdispl_px_list.append(xdispl_px)
	ydispl_px_list.append(ydispl_px); xerr_px_list.append(xerr_px); yerr_px_list.append(yerr_px)

num_recs = len(px_j_list)


# converts lists into numpy arrays
px_j_array = np.array(px_j_list, int) 
px_i_array = np.array(px_i_list, int)  
totdispl_px_array = np.array(totdispl_px_list, float)  
corrstr_array = np.array(corrstr_list, float)  
resflg_array = np.array(resflg_list, int) 
xdispl_px_array = np.array(xdispl_px_list, float)  
ydispl_px_array = np.array(ydispl_px_list, float)  
xerr_px_array = np.array(xerr_px_list, float)  
yerr_px_array = np.array(yerr_px_list, float) 


# conversions from graphic to geographic coordinates  
   
x0_array = x_min + px_j_array*cell_size
y0_array = y_max - px_i_array*cell_size
dx_array = xdispl_px_array*cell_size
dy_array = -ydispl_px_array*cell_size

x1_array = x0_array + dx_array
y1_array = y0_array + dy_array	

x_err_array = xerr_px_array*cell_size
y_err_array = yerr_px_array*cell_size

magn_displ_array = np.sqrt(dx_array**2+dy_array**2)

dir_array = (np.arctan2(dx_array,dy_array)*180.0/pi)
dir_array = np.where(dir_array < 0.0, dir_array + 360.0, dir_array)
dir_array = np.where(dir_array > 360.0, dir_array - 360.0, dir_array)
displ_dir_array = np.where(magn_displ_array > 1.0e-01, dir_array, -999)


# writes results in output file
outfile.write('px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg\n')
for n in range(num_recs):
	outfile.write(str(px_j_array[n])+','+str(px_i_array[n])+','+str(x0_array[n])+','+str(y0_array[n])+','+ \
		str(x1_array[n])+','+str(y1_array[n])+','+str(dx_array[n])+','+str(dy_array[n])+','+str(x_err_array[n])+','+ \
		str(y_err_array[n])+','+str(magn_displ_array[n])+','+str(displ_dir_array[n])+','+str(corrstr_array[n])+','+ \
		str(resflg_array[n])+'\n')
	

