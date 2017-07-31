# Mauro Alberti http://malg.eu, http://www.gistrutturale.it
# GNU General Public License v. 3

# 2010-10-12: modified output esri grid function to deal with null values INSIDE input grids
# 2010-10-28: modified calculation of dvx_dy 
# 2010-10-29: set vtk null value to -99.0


import os, sys, string
from math import *
import numpy as np



# 3D point class
class Point:
	def __init__(self, x, y, z=0.0):
		self.x = x
		self.y = y
		self.z = z

	# calculates euclidean distance between two points
	def distance(self, pt):
		return sqrt((self.x - pt.x)**2 + (self.y - pt.y)**2 + (self.z - pt.z)**2)


# 3D rectangle class - aligned along cartesian axes
class Rectangle:
	def __init__(self, pt_init, pt_end):
		self.pt_init = pt_init
		self.pt_end = pt_end

	def xrange(self):
		return abs(self.pt_end.x-self.pt_init.x)

	def yrange(self):
		return abs(self.pt_end.y-self.pt_init.y)
		
	def horiz_area(self):
		return self.xrange()*self.yrange()

		
# rectangular spatial domain class
class SpatialDomain:

	def __init__(self):
		self.extent = None	
		self.unit = None			
		
	def set_extent(self, pt_init, pt_end):
		self.extent = Rectangle(pt_init, pt_end)

	def set_unit(self, string_unit):		
		self.unit = string_unit		
		
	def start_point(self):
		return self.extent.pt_init

	def end_point(self):
		return self.extent.pt_end

		
# gridded dataset expressed by an array, the measure unit and the no-data value		
class GriddedData:
	def __init__(self, array_values, sg_unit=''):
		self.unit = sg_unit	
		self.values = array_values		

		
# field class
class Field:

	def __init__(self):
		self.spatialdomain = None	
		self.griddeddata = None			

	def set_spatialdomain(self, pt_init, pt_end):
		self.spatialdomain = SpatialDomain()
		self.spatialdomain.set_extent(pt_init, pt_end)	

	# returns the lower-left corner of the grid domain 		
	def get_init_point(self):
		return self.spatialdomain.start_point()		
		
	# returns the upper-right corner of the grid domain 
	def get_end_point(self):
		return self.spatialdomain.end_point()	

	# returns the width (x range) of the spatial domain
	def field_width(self):
		return self.spatialdomain.extent.xrange()

	# returns the height (y range) of the spatial domain
	def field_height(self):
		return self.spatialdomain.extent.yrange()		
	
	# row number of the grid domain 		
	def nrows(self):
		return np.shape(self.griddeddata.values)[0]		

	# column number of the grid domain 		
	def ncols(self):
		return np.shape(self.griddeddata.values)[1]		
			
	# returns the cell size in the x direction of the gridded dataset
	def cellsize_x(self):
		return self.field_width()/float(self.ncols())

	# returns the cell size in the y direction of the gridded dataset
	def cellsize_y(self):
		return self.field_height()/float(self.nrows())

	# returns the mean cell size 
	def cellsize_mean(self):
		return (self.cellsize_x()+self.cellsize_y())/2.0		
		
	# checks equivalence between the parameters of two fields
	def check_geoequiv(self, sf, tolerance): 	
		distance_init_pts = self.get_init_point().distance(sf.get_init_point()) # distance between the start points (i.e. lower-left) of the two grids
		distance_end_pts = self.get_end_point().distance(sf.get_end_point()) # distance between the end points (i.e. upper-right) of the two grids	
		# checks if the two grid geographical parameters can be considered equivalent
		if self.ncols() != sf.ncols() or self.nrows() != sf.nrows() or \
			abs((self.cellsize_mean() - sf.cellsize_mean())/self.cellsize_mean()) > tolerance or \
			distance_init_pts > abs(tolerance*self.cellsize_mean()) or \
			distance_end_pts > abs(tolerance*self.cellsize_mean()):    
			return False
		else:
			return True

			
# 2D scalar field class		
class ScalarField(Field):

	"""
	2D scalar field 
	
    """
	

	def set_scalargrid(self, array_values):		
		self.griddeddata = GriddedData(array_values)


	# sets the field values from ESRI ascii grid parameters 		
	def fromesri(self, esri_params, grid_values):
		ncols = esri_params[0]
		nrows = esri_params[1]
		xstart = esri_params[2]
		ystart = esri_params[3]
		cellsize = esri_params[4]
		esri_nullvalue = esri_params[5]

		pt_init = Point(xstart, ystart)
		pt_end = Point(xstart+cellsize*ncols, ystart+cellsize*nrows)		
		self.set_spatialdomain(pt_init, pt_end)	
			
		grid_values = np.asarray(grid_values).reshape(nrows, ncols)
		values = np.where(abs(grid_values-esri_nullvalue)> 0.1, grid_values, np.NaN)		
		self.set_scalargrid(values)

		
	def write_esrigrid(self, outgrid_fn, esri_nullvalue):

		# checking existance of output slope grid
		if os.path.exists(outgrid_fn):
		  os.remove(outgrid_fn)

		outputgrid = open(outgrid_fn, 'w') #create the output ascii file

		if outputgrid is None:
		  print 'could not create output grid file'
		  sys.exit(1)

		# writes header of grid ascii file
		outputgrid.write('NCOLS %d\n' % self.ncols())
		outputgrid.write('NROWS %d\n' % self.nrows())
		outputgrid.write('XLLCORNER %.6f\n' % self.get_init_point().x)
		outputgrid.write('YLLCORNER %.6f\n' % self.get_init_point().y)
		outputgrid.write('CELLSIZE %.6f\n' % self.cellsize_mean())
		outputgrid.write('NODATA_VALUE %.6f\n' % esri_nullvalue)

		esrigrid_outvalues = np.where(np.isnan(self.griddeddata.values), esri_nullvalue, self.griddeddata.values)
		
		#output of results
		for i in range(0, self.nrows()):
			for j in range(0, self.ncols()):
				outputgrid.write('%.6f ' % (esrigrid_outvalues[i,j]))
			outputgrid.write('\n')

		outputgrid.close()

		
	def write2vtk(self, outvtk, name, vtk_nullvalue):
		outvtk.write('SCALARS %s float\n' % name)
		outvtk.write('LOOKUP_TABLE default\n')
		vtk_outvalues = np.where(np.isnan(self.griddeddata.values), vtk_nullvalue, self.griddeddata.values)
		vtk_outvalues = np.flipud(vtk_outvalues)
		for i in range(self.nrows()):
			for j in range(self.ncols()):
				outvtk.write('%f ' % (vtk_outvalues[i,j]))
		outvtk.write('\n')


# gridded dataset expressed by an array, the measure unit and the no-data value		
class VectorGrid:
	def __init__(self, array_values_x, array_values_y, vg_unit=''):
		self.unit = vg_unit	
		self.xvalues = array_values_x		
		self.yvalues = array_values_y
		
		
class VectorField(Field):

	def set_vectorgrid(self, array_values_x, array_values_y):		
		self.griddeddata = VectorGrid(array_values_x, array_values_y)	

	# row number of the grid domain 		
	def nrows(self):
		return np.shape(self.griddeddata.xvalues)[0]		

	# column number of the grid domain 		
	def ncols(self):
		return np.shape(self.griddeddata.xvalues)[1]	
		
	# returns the cell size in the x direction of the gridded dataset
	def cellsize_x(self):
		return self.field_width()/float(self.ncols())

	# returns the cell size in the y direction of the gridded dataset
	def cellsize_y(self):
		return self.field_height()/float(self.nrows())
		
	# calculates divergence 
	def divergence(self):
		
		dvx_dx = np.gradient(self.griddeddata.xvalues)[1]
		dvy_dy = -(np.gradient(self.griddeddata.yvalues)[0])		

		divergence_fld = ScalarField()
		divergence_fld.set_spatialdomain(self.get_init_point(), self.get_end_point())
		divergence_fld.set_scalargrid((dvx_dx + dvy_dy)/self.cellsize_mean())
		
		return divergence_fld
		

	# calculates curl 
	def curl(self):
	
		dvy_dx = np.gradient(self.griddeddata.yvalues)[1]
		dvx_dy = -(np.gradient(self.griddeddata.xvalues)[0])
		
		curl_fld = ScalarField()
		curl_fld.set_spatialdomain(self.get_init_point(), self.get_end_point())
		curl_fld.set_scalargrid((dvy_dx - dvx_dy)/self.cellsize_mean())
		
		return curl_fld


	# calculates speed variations along flow lines
	def vect_magn_grad(self):	

		vx, vy = self.griddeddata.xvalues, self.griddeddata.yvalues
		
		dir_array = np.arctan2(vx, vy)
		
		vect_magn = np.sqrt(vx**2 + vy**2)
		dm_dy, dm_dx = np.gradient(vect_magn)
		dm_dy = - dm_dy
		
		acceleration = dm_dx * np.sin(dir_array) + dm_dy * np.cos(dir_array)
		
		vect_magn_grad_fld = ScalarField()
		vect_magn_grad_fld.set_spatialdomain(self.get_init_point(), self.get_end_point())
		vect_magn_grad_fld.set_scalargrid(acceleration/self.cellsize_mean())
		
		return vect_magn_grad_fld
		
		
	def write2vtk(self, outvtk, vtk_nullvalue):

		outvtk.write('# vtk DataFile Version 3.0\n')
		outvtk.write('Vector field\n\n')
		outvtk.write('ASCII\n')
		outvtk.write('DATASET STRUCTURED_POINTS\n')
		outvtk.write('DIMENSIONS %d %d 1\n' % (self.ncols(), self.nrows()))
		outvtk.write('ORIGIN %f %f 0.0\n' % (self.get_init_point().x, self.get_init_point().y))
		outvtk.write('SPACING %f %f 1\n\n' % (self.cellsize_mean(), self.cellsize_mean()))
		 
		outvtk.write('POINT_DATA %d\n\n' % (self.ncols()*self.nrows()))
			
		# writing vector field values
		outvtk.write('VECTORS vector float\n')
		vx_values = np.where(np.isnan(self.griddeddata.xvalues), vtk_nullvalue, self.griddeddata.xvalues)
		vy_values = np.where(np.isnan(self.griddeddata.yvalues), vtk_nullvalue, self.griddeddata.yvalues)
		vx_values = np.flipud(vx_values); vy_values = np.flipud(vy_values)
		for i in range(vx_field.nrows()):
			for j in range(vx_field.ncols()):
				outvtk.write('%f %f 0.0\n' % (vx_values[i,j], vy_values[i,j]))

				
# reads grid geographical parameters from an esri ascii grid 		
def read_esri_params(asciigrid):
	try:

		# reads header lines content and sets variables
		ncols = int(asciigrid.readline().split()[1])
		nrows = int(asciigrid.readline().split()[1])
		cellx_line = asciigrid.readline().split()
		celly_line = asciigrid.readline().split()
		cellsize = float(asciigrid.readline().split()[1])
		nodata_value = float(asciigrid.readline().split()[1])
		
		# checks and extracts lower-left corner/center
		xstart_type = cellx_line[0].lower()
		ystart_type = celly_line[0].lower()
		if xstart_type[:3] != 'xll' or ystart_type[:3] != 'yll' or \
		(xstart_type[3:] != ystart_type[3:]): 
			return None				
		start_type = xstart_type[3:] # extracts 'corner' or 'center'
		xstart_value = float(cellx_line[1])	# lower-left x-value
		ystart_value = float(celly_line[1])	# lower-left y-value
		
		# defines lower-left corner coordinates 
		if start_type == 'corner':
			xstart_value = xstart_value
			ystart_value = ystart_value				
		elif start_type == 'center':
			xstart_value = xstart_value - cellsize/2.0
			ystart_value = ystart_value - cellsize/2.0
		else:
			return None


		# validates geographical parameters of grid		

		if ncols > 0 and \
		   nrows > 0 and \
		   xstart_value != None and \
		   ystart_value != None and \
		   cellsize > 0.0 and \
		   nodata_value != None:
			return (ncols, nrows, xstart_value, ystart_value, cellsize, nodata_value)
		else: 
			return None
		
	except:
	
		return None
		
		
# reads grid data (geographical parameters and values) from an esri ascii file		
def read_from_esrigrid(in_ascii_grid):	
	# open grid input file 
	try:
		ingrid = open(in_ascii_grid, 'r')
	except:
		print 'unable to open %s grid' % in_ascii_grid
		sys.exit(1)

	# reads ESRI parameters
	esri_params = read_esri_params(ingrid)	
	if esri_params == None:
		print 'unable to read geographical parameters of %s grid' % in_ascii_grid
		sys.exit(1)
	
	# reads input grid data 
	grid_values = [float(z) for z in ingrid.read().split()]
	if grid_values == None:
		print 'unable to read data of %s grid' % in_ascii_grid
		sys.exit(1)

	# closes input grid file
	ingrid.close()		
				
	# results
	scalar_fld = ScalarField()
	scalar_fld.fromesri(esri_params, grid_values)
	return scalar_fld


def input_error():
	print \
""" 
vector_field_par.py - script for determining vector field parameters
				
Usage: python vector_field_par.py parameter_file     
Example: python vector_field_par.py par.txt

Parameters in parameter file
  in_vx_grid: ascii grid with x cartesian components of vector field (input)
  in_vy_grid: ascii grid with y cartesian components of vector field (input)
  out_divergence_grid: ascii grid with divergence values (output)
  out_curl_grid_fn: ascii grid with curl magnitudes (output)
  out_acceleration: ascii grid with gradients of vector magnitudes alogn flow lines (output)
  out_vtk_file: vtk file storing all input and output data (output)
  
  vx.asc vy.asc div.asc curl.asc acc.asc output.vtk
  
  
  
"""
	sys.exit(1)
	
	
	
############
# main

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
in_vx_grid_fn =  parameter_file.readline().split()[0]
in_vy_grid_fn =  parameter_file.readline().split()[0]
out_div_grid_fn =  parameter_file.readline().split()[0]
out_curl_grid_fn =  parameter_file.readline().split()[0]
out_acc_grid_fn =  parameter_file.readline().split()[0]
out_vtk_fn =  parameter_file.readline().split()[0]

    
# read input data (from vx and vy grids)	
vx_field = read_from_esrigrid(in_vx_grid_fn)
if vx_field == None: 
	print 'Unable to read vx data'
	sys.exit(1)	
vy_field = read_from_esrigrid(in_vy_grid_fn)
if vy_field == None: 
	print 'Unable to read vy data'
	sys.exit(1)
	
	
# check equivalence of gegraphical parameters between the two grids			
tolerance = 1.0/1000.0  # max allowed relative difference between the two grid parameters
if not vx_field.check_geoequiv(vy_field, tolerance):
	print 'The geographical parameters of the two grids are not equivalent'
	sys.exit(1)

	
# creation of vector field
vect_fld = VectorField()
vect_fld.set_spatialdomain(vx_field.get_init_point(), vx_field.get_end_point())
vect_fld.set_vectorgrid(vx_field.griddeddata.values, vy_field.griddeddata.values)
	
	
# divergence, curl and speed variations calculations
divergence_fld = vect_fld.divergence()	
curl_fld = vect_fld.curl()
vect_magn_grad_fld = vect_fld.vect_magn_grad()


# output to esri-format ascii grid
esri_nullvalue = -9999
divergence_fld.write_esrigrid(out_div_grid_fn, esri_nullvalue)					
curl_fld.write_esrigrid(out_curl_grid_fn, esri_nullvalue)
vect_magn_grad_fld.write_esrigrid(out_acc_grid_fn, esri_nullvalue)
 
 
# output to vtk file 
vtk_nullvalue = -99.0 
outvtk = open(out_vtk_fn, 'w') 
if outvtk is None:
  print 'could not create output vtk file'
  sys.exit(1)
vect_fld.write2vtk(outvtk, vtk_nullvalue) # writes header of vtk file

divergence_fld.write2vtk(outvtk, 'divergence', vtk_nullvalue)
curl_fld.write2vtk(outvtk, 'curl', vtk_nullvalue)
vect_magn_grad_fld.write2vtk(outvtk, 'vect_magn_grad', vtk_nullvalue)

	
# closes output vtk file		
outvtk.close()






