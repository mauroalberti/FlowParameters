# Mauro Alberti  2010-08-24 - http://www.malg.eu; http://www.gistrutturale.it
# GNU General Public License v. 3


import sys, os
from math import *

try:
	from osgeo import ogr
except:
	import ogr
	
	
# Imcorr output record
class Geo_rec:
		
	def __init__(self,px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg,magn_diff=-999,ang_diff=-999):	
		self.px_j = int(px_j)	
		self.px_i = int(px_i)
		self.x0 = float(x0)			
		self.y0 = float(y0)
		self.x1 = float(x1)			
		self.y1 = float(y1)
		self.dx = float(dx)			
		self.dy = float(dy)
		self.x_err = float(x_err)			
		self.y_err = float(y_err)
		self.magn_displ = float(magn_displ)			
		self.displ_dir = float(displ_dir)		
		self.corrstr = float(corrstr)
		self.resflg = int(resflg)		
		self.magn_diff = float(magn_diff)			
		self.ang_diff = float(ang_diff)
		


def read_line_file(ingeofile):

	# reads header line
	header_line = ingeofile.readline()	
	num_input_fields = len(header_line.split(','))

	georec_list = []

	if num_input_fields == 16:
		for line in ingeofile:
			(px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg,magn_diff,ang_diff) = line.split(',')
			geo_record = Geo_rec(px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg,magn_diff,ang_diff)
			georec_list.append(geo_record)
	elif num_input_fields == 14:
		for line in ingeofile:
			(px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg) = line.split(',')
			geo_record = Geo_rec(px_j,px_i,x0,y0,x1,y1,dx,dy,x_err,y_err,magn_displ,displ_dir,corrstr,resflg)
			georec_list.append(geo_record)
	else:
	  print 'Field number error in input file'
	  sys.exit(1)	

	return (georec_list, num_input_fields)
	

	
def input_error():
	print \
""" 
vectors2shapefile_02.py - script for converting imcorr output
						ta linear shapefile

Usage: vectors2shapefile_02.py parameter_file     
Example: vectors2shapefile_02.py param.txt

Parameters in parameter file
  input_geovectors.out: text file (output of imcorr2geog_xx.py/vectors_coherence_xx.py)
  out_shapefile: output displacement vectors (line shapefile)
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
infilename =  parameter_file.readline().split()[0]
out_shapefnm =  parameter_file.readline().split()[0]

  
# open input file 
try:
	ingeofile = open(infilename, 'r')
except:
  print 'error in reading input file'
  sys.exit(1)
  

# reads data from file
georec_list, num_input_fields = read_line_file(ingeofile) 


# Creation of output shapefile

driver = ogr.GetDriverByName('ESRI Shapefile')

if os.path.exists(out_shapefnm):
    driver.DeleteDataSource(out_shapefnm)

out_shape = driver.CreateDataSource(out_shapefnm)
if out_shape is None:
    print 'Could not create output shapefile'
    sys.exit(1)
out_layer = out_shape.CreateLayer('out_lines', geom_type=ogr.wkbLineString)
	
	
# add fields to the output shapefile
px_i_fieldDef = ogr.FieldDefn('px_i', ogr.OFTInteger)
out_layer.CreateField(px_i_fieldDef)
px_j_fieldDef = ogr.FieldDefn('px_j', ogr.OFTInteger)
out_layer.CreateField(px_j_fieldDef)


x0_fieldDef = ogr.FieldDefn('x0', ogr.OFTReal)
out_layer.CreateField(x0_fieldDef)
y0_fieldDef = ogr.FieldDefn('y0', ogr.OFTReal)
out_layer.CreateField(y0_fieldDef)
x1_fieldDef = ogr.FieldDefn('x1', ogr.OFTReal)
out_layer.CreateField(x1_fieldDef)
y1_fieldDef = ogr.FieldDefn('y1', ogr.OFTReal)
out_layer.CreateField(y1_fieldDef)

dx_fieldDef = ogr.FieldDefn('dx', ogr.OFTReal)
out_layer.CreateField(dx_fieldDef)
dy_fieldDef = ogr.FieldDefn('dy', ogr.OFTReal)
out_layer.CreateField(dy_fieldDef)

magn_displ_fieldDef = ogr.FieldDefn('magn_displ', ogr.OFTReal)
out_layer.CreateField(magn_displ_fieldDef)
displ_dir_fieldDef = ogr.FieldDefn('displ_dir', ogr.OFTReal)
out_layer.CreateField(displ_dir_fieldDef)


x_err_fieldDef = ogr.FieldDefn('x_err', ogr.OFTReal)
out_layer.CreateField(x_err_fieldDef)
y_err_fieldDef = ogr.FieldDefn('y_err', ogr.OFTReal)
out_layer.CreateField(y_err_fieldDef)


corrstr_fieldDef = ogr.FieldDefn('corrstr', ogr.OFTReal)
out_layer.CreateField(corrstr_fieldDef)
resflg_fieldDef = ogr.FieldDefn('resflg', ogr.OFTInteger)
out_layer.CreateField(resflg_fieldDef)

if num_input_fields == 16:	
	magn_diff_fieldDef = ogr.FieldDefn('magn_diff', ogr.OFTReal)
	out_layer.CreateField(magn_diff_fieldDef)
	ang_diff_fieldDef = ogr.FieldDefn('ang_diff', ogr.OFTReal)
	out_layer.CreateField(ang_diff_fieldDef)



# get the layer definition of the output shapefile
outshape_featdef = out_layer.GetLayerDefn()


# processes imcorr records
for rec in georec_list:
	
	# adds line
	line = ogr.Geometry(ogr.wkbLineString)
	line.AddPoint(rec.x0,rec.y0)
	line.AddPoint(rec.x1,rec.y1)	

	result_feature = ogr.Feature(outshape_featdef)
	result_feature.SetGeometry(line)
	result_feature.SetField('px_i', rec.px_i)
	result_feature.SetField('px_j', rec.px_j)

	result_feature.SetField('x0', rec.x0)
	result_feature.SetField('y0', rec.y0)

	result_feature.SetField('x1', rec.x1)
	result_feature.SetField('y1', rec.y1)

	result_feature.SetField('dx', rec.dx)
	result_feature.SetField('dy', rec.dy)

	result_feature.SetField('magn_displ', rec.magn_displ)
	result_feature.SetField('displ_dir', rec.displ_dir)

	result_feature.SetField('x_err', rec.x_err)
	result_feature.SetField('y_err', rec.y_err)

	result_feature.SetField('corrstr', rec.corrstr)
	result_feature.SetField('resflg', rec.resflg)

	if num_input_fields == 16:
		result_feature.SetField('magn_diff', rec.magn_diff)
		result_feature.SetField('ang_diff', rec.ang_diff)

	out_layer.CreateFeature(result_feature)

	# destroy not longer used objects
	line.Destroy()
	result_feature.Destroy() 

out_shape.Destroy()
 

