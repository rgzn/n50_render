#!/usr/bin/python2.7

from osgeo import gdal
from osgeo import osr

import sys
import os
import math
import struct
import urllib
import subprocess
import time
from io import BytesIO
from PIL import Image

palette_method = "pngquant"

s3 = {}
try:
	from boto.s3.connection import S3Connection
	from boto.s3.key import Key
	file = open(os.path.expanduser("~/.awssecret"), "r")	
	lines = file.readlines()
	file.close()
	s3['connection'] = S3Connection(lines[0].strip(), lines[1].strip())
except:
	pass

resampling_list = ('average','near','bilinear','cubic','cubicspline','lanczos','aspect','max','min')

import math

class GlobalMercator(object):

	def __init__(self, tileSize=256):
		"Initialize the TMS Global Mercator pyramid"
		self.tileSize = tileSize
		self.initialResolution = 2 * math.pi * 6378137 / self.tileSize
		# 156543.03392804062 for tileSize 256 pixels
		self.originShift = 2 * math.pi * 6378137 / 2.0
		# 20037508.342789244

	def LatLonToMeters(self, lat, lon ):
		"Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"

		mx = lon * self.originShift / 180.0
		my = math.log( math.tan((90 + lat) * math.pi / 360.0 )) / (math.pi / 180.0)

		my = my * self.originShift / 180.0
		return mx, my

	def MetersToLatLon(self, mx, my ):
		"Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum"

		lon = (mx / self.originShift) * 180.0
		lat = (my / self.originShift) * 180.0

		lat = 180 / math.pi * (2 * math.atan( math.exp( lat * math.pi / 180.0)) - math.pi / 2.0)
		return lat, lon

	def PixelsToMeters(self, px, py, zoom):
		"Converts pixel coordinates in given zoom level of pyramid to EPSG:900913"

		res = self.Resolution( zoom )
		mx = px * res - self.originShift
		my = py * res - self.originShift
		return mx, my
		
	def MetersToPixels(self, mx, my, zoom):
		"Converts EPSG:900913 to pyramid pixel coordinates in given zoom level"
				
		res = self.Resolution( zoom )
		px = (mx + self.originShift) / res
		py = (my + self.originShift) / res
		return px, py
	
	def PixelsToTile(self, px, py):
		"Returns a tile covering region in given pixel coordinates"

		tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
		ty = int( math.ceil( py / float(self.tileSize) ) - 1 )
		return tx, ty

	def PixelsToRaster(self, px, py, zoom):
		"Move the origin of pixel coordinates to top-left corner"
		
		mapSize = self.tileSize << zoom
		return px, mapSize - py
		
	def MetersToTile(self, mx, my, zoom):
		"Returns tile for given mercator coordinates"
		
		px, py = self.MetersToPixels( mx, my, zoom)
		return self.PixelsToTile( px, py)

	def TileBounds(self, tx, ty, zoom):
		"Returns bounds of the given tile in EPSG:900913 coordinates"
		
		minx, miny = self.PixelsToMeters( tx*self.tileSize, ty*self.tileSize, zoom )
		maxx, maxy = self.PixelsToMeters( (tx+1)*self.tileSize, (ty+1)*self.tileSize, zoom )
		return ( minx, miny, maxx, maxy )

	def TileLatLonBounds(self, tx, ty, zoom ):
		"Returns bounds of the given tile in latutude/longitude using WGS84 datum"

		bounds = self.TileBounds( tx, ty, zoom)
		minLat, minLon = self.MetersToLatLon(bounds[0], bounds[1])
		maxLat, maxLon = self.MetersToLatLon(bounds[2], bounds[3])
		 
		return ( minLat, minLon, maxLat, maxLon )
		
	def Resolution(self, zoom ):
		"Resolution (meters/pixel) for given zoom level (measured at Equator)"
		
		# return (2 * math.pi * 6378137) / (self.tileSize * 2**zoom)
		return self.initialResolution / (2**zoom)
		
	def ZoomForPixelSize(self, pixelSize ):
		"Maximal scaledown zoom of the pyramid closest to the pixelSize."
		
		for i in range(32):
			if pixelSize > self.Resolution(i):
				if i!=0:
					return i-1
				else:
					return 0 # We don't want to scale up
		
	def GoogleTile(self, tx, ty, zoom):
		"Converts TMS tile coordinates to Google Tile coordinates"
		
		# coordinate origin is moved from bottom-left to top-left corner of the extent
		return tx, (2**zoom - 1) - ty

	def QuadTree(self, tx, ty, zoom ):
		"Converts TMS tile coordinates to Microsoft QuadTree"
		
		quadKey = ""
		ty = (2**zoom - 1) - ty
		for i in range(zoom, 0, -1):
			digit = 0
			mask = 1 << (i-1)
			if (tx & mask) != 0:
				digit += 1
			if (ty & mask) != 0:
				digit += 2
			quadKey += str(digit)
			
		return quadKey

#---------------------

class GDAL2Tiles(object):

	# -------------------------------------------------------------------------
	def process(self):
		"""The main processing function, runs all the main steps of processing"""
		
		# Opening and preprocessing of the input file
		self.open_input()

		for tz in range(self.tmaxz, self.tminz-1, -1):
			if self.tmaxz == tz and not(self.options.overview):
				if self.options.wms:
					print "fetching tiles from wms"
					self.wms_fetch_tiles(tz)
				else:
					print "Writing tiles for zoom level", tz
					self.write_tiles(tz)
			else:
				print "Compositing tiles for zoom level", tz
				self.composite_tiles(tz)
		
	# -------------------------------------------------------------------------
	def error(self, msg, details = "" ):
		"""Print an error message and stop the processing"""

		if details:
			self.parser.error(msg + "\n\n" + details)
		else:	
			self.parser.error(msg)
		
	# -------------------------------------------------------------------------
	def progressbar(self, complete = 0.0):
		"""Print progressbar for float value 0..1"""
		
		gdal.TermProgress_nocb(complete)

	# -------------------------------------------------------------------------
	def stop(self):
		"""Stop the rendering immediately"""
		self.stopped = True

	# -------------------------------------------------------------------------
	def __init__(self, arguments ):
		"""Constructor function - initialization"""
		
		self.stopped = False
		self.input = None
		self.output = None

		# Tile format
		self.tilesize = 256
		self.tiledriver = 'PNG'
		self.tileext = 'png'
		
		# Should we read bigger window of the input raster and scale it down?
		# Note: Modified leter by open_input()
		# Not for 'near' resampling
		# Not for Wavelet based drivers (JPEG2000, ECW, MrSID)
		# Not for 'raster' profile
		self.scaledquery = True
		# How big should be query window be for scaling down
		# Later on reset according the chosen resampling algorightm
		self.querysize = 4 * self.tilesize

		# Should we use Read on the input file for generating overview tiles?
		# Note: Modified later by open_input()
		# Otherwise the overview tiles are generated from existing underlying tiles
		self.overviewquery = False
		
		# RUN THE ARGUMENT PARSER:
		
		self.optparse_init()
		self.options, self.args = self.parser.parse_args(args=arguments)
		if not self.args:
			self.error("No input file specified")

		if not(self.options.bucket==None):
			s3['bucket'] = s3['connection'].get_bucket(self.options.bucket)
		# POSTPROCESSING OF PARSED ARGUMENTS:

		# Workaround for old versions of GDAL
		try:
			if (self.options.verbose and self.options.resampling == 'near') or gdal.TermProgress_nocb:
				pass
		except:
			self.error("This version of GDAL is not supported. Please upgrade to 1.6+.")
			#,"You can try run crippled version of gdal2tiles with parameters: -v -r 'near'")
		
		# Is output directory the last argument?

		# Test output directory, if it doesn't exist
		if os.path.isdir(self.args[-1]) or ( len(self.args) > 1 and not os.path.exists(self.args[-1])):
			self.output = self.args[-1]
			self.args = self.args[:-1]

		# More files on the input not directly supported yet
		
		if (len(self.args) > 1):
			self.error("Processing of several input files is not supported.",
			"""Please first use a tool like gdal_vrtmerge.py or gdal_merge.py on the files:
gdal_vrtmerge.py -o merged.vrt %s""" % " ".join(self.args))
			# TODO: Call functions from gdal_vrtmerge.py directly
			
		self.input = self.args[0]
		
		# Default values for not given options
		
		if not self.output:
			# Directory with input filename without extension in actual directory
			self.output = os.path.splitext(os.path.basename( self.input ))[0]
				
		if self.options.url and not self.options.url.endswith('/'):
			self.options.url += '/'
		if self.options.url:
			self.options.url += os.path.basename( self.output ) + '/'

		# Supported options
		
		self.colors = 256
		if(self.options.colors > 0):
			self.colors = int(self.options.colors)

		self.resampling = None
		
		if self.options.resampling == 'average':
			try:
				if gdal.RegenerateOverview:
					pass
			except:
				self.error("'average' resampling algorithm is not available.", "Please use -r 'near' argument or upgrade to newer version of GDAL.")
		
		elif self.options.resampling == 'near':
			self.resampling = gdal.GRA_NearestNeighbour
			self.querysize = self.tilesize
			
		elif self.options.resampling == 'bilinear':
			self.resampling = gdal.GRA_Bilinear
			self.querysize = self.tilesize * 2

		elif self.options.resampling == 'cubic':
			self.resampling = gdal.GRA_Cubic

		elif self.options.resampling == 'cubicspline':
			self.resampling = gdal.GRA_CubicSpline

		elif self.options.resampling == 'lanczos':
			self.resampling = gdal.GRA_Lanczos
		
		# User specified zoom levels
		self.tminz = None
		self.tmaxz = None
		if self.options.zoom:
			minmax = self.options.zoom.split('-',1)
			minmax.extend([''])
			min, max = minmax[:2]
			self.tminz = int(min)
			if max:
				self.tmaxz = int(max)
			else:
				self.tmaxz = int(min) 
		
		# Output the results

		if self.options.verbose:
			print("Options:", self.options)
			print("Input:", self.input)
			print("Output:", self.output)
			print("Cache: %s MB" % (gdal.GetCacheMax() / 1024 / 1024))
			print('')

	# -------------------------------------------------------------------------
	def optparse_init(self):
		"""Prepare the option parser for input (argv)"""
		
		from optparse import OptionParser, OptionGroup
		usage = "Usage: %prog [options] input_file(s) [output]"
		p = OptionParser(usage, version="%prog")
		p.add_option("-r", "--resampling", dest="resampling", type='choice', choices=resampling_list,
						help="Resampling method (%s) - default 'average'" % ",".join(resampling_list))
		p.add_option('-s', '--s_srs', dest="s_srs", metavar="SRS",
						  help="The spatial reference system used for the source input data")
		p.add_option('-z', '--zoom', dest="zoom",
						  help="Zoom levels to render (format:'2-5' or '10').")
		p.add_option('-e', '--resume', dest="resume", action="store_true",
						  help="Resume mode. Generate only missing files.")
		p.add_option('-a', '--srcnodata', dest="srcnodata", metavar="NODATA",
			  			  help="NODATA transparency value to assign to the input data")
		p.add_option("-v", "--verbose",
						  action="store_true", dest="verbose",
						  help="Print status messages to stdout")
		p.add_option("-o", "--overview", dest="overview", action="store_true",
					help="Skip base layer; generate overview files only.")
		p.add_option("-w", "--wms", dest="wms",
					help="Fetch base level tiles from WMS server")
		p.add_option("-c", '--crop', dest="crop", action="store_true", help="assume that image is 100% alpha if cropped near borders")
		p.add_option("-n", '--colors', dest="colors", help="# of colortable entries (default is 256).")
		p.add_option("-f", '--fetch', dest="fetch", help="Base URL to fetch tiles from.")
		p.add_option("-b", '--bucket', dest="bucket", help="S3 bucket URL for upload (requires timkay's aws)")
		p.add_option("-g", '--grayscale', dest="grayscale", action="store_true", help="Source file is a 1-layer grayscale image")
		p.add_option("-2", '--sixteen', dest="long", action="store_true", help="Write destination as a single layer 16 bit image")
		p.add_option("-1", '--eight', dest="short", action="store_true", help="Write destination as a single layer 8 bit image")

		p.set_defaults(verbose=False, profile="mercator", kml=False, url='',
		webviewer='all', copyright='', resampling='average', resume=False,
		googlekey='INSERT_YOUR_KEY_HERE', yahookey='INSERT_YOUR_YAHOO_APP_ID_HERE')

		self.parser = p
		
	# -------------------------------------------------------------------------
	def open_input(self):
		"""Initialization of the input raster, reprojection if necessary"""
		
		gdal.AllRegister()

		# Initialize necessary GDAL drivers
		
		self.out_drv = gdal.GetDriverByName( self.tiledriver )
		self.mem_drv = gdal.GetDriverByName( 'MEM' )
		self.in_drv = gdal.GetDriverByName( self.tiledriver )
		gdal.SetConfigOption("GDAL_FORCE_CACHING", "true");
		
		if not self.out_drv:
			raise Exception("The '%s' driver was not found, is it available in this GDAL build?", self.tiledriver)
		if not self.mem_drv:
			raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?")
		
		# Open the input file
		
		if self.input:
			self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
		else:
			raise Exception("No input file was specified")

		if self.options.verbose:
			print("Input file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

		if not self.in_ds:
			# Note: GDAL prints the ERROR message too
			self.error("It is not possible to open the input file '%s'." % self.input )
			
		# Read metadata from the input file
		if self.in_ds.RasterCount == 0:
			self.error( "Input file '%s' has no raster band" % self.input )
			
		if self.in_ds.GetRasterBand(1).GetRasterColorTable():
			# TODO: Process directly paletted dataset by generating VRT in memory
			print "WARNING - Color table found, will be ignored.  You probably want to convert this file to RGB/RGBA.  Poceeding Anyway"

		# Get NODATA value
		self.in_nodata = []
		for i in range(1, self.in_ds.RasterCount+1):
			if self.in_ds.GetRasterBand(i).GetNoDataValue() != None:
				self.in_nodata.append( self.in_ds.GetRasterBand(i).GetNoDataValue() )
		if self.options.srcnodata:
			nds = list(map( float, self.options.srcnodata.split(',')))
			if len(nds) < self.in_ds.RasterCount:
				self.in_nodata = (nds * self.in_ds.RasterCount)[:self.in_ds.RasterCount]
			else:
				self.in_nodata = nds

		if self.options.verbose:
			print("NODATA: %s" % self.in_nodata)

		#
		# Here we should have RGBA input dataset opened in self.in_ds
		#

		if self.options.verbose:
			print("Preprocessed file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

		# Spatial Reference System of the input raster


		self.in_srs = None
		
		if self.options.s_srs:
			self.in_srs = osr.SpatialReference()
			self.in_srs.SetFromUserInput(self.options.s_srs)
			self.in_srs_wkt = self.in_srs.ExportToWkt()
		else:
			self.in_srs_wkt = self.in_ds.GetProjection()
			if not self.in_srs_wkt and self.in_ds.GetGCPCount() != 0:
				self.in_srs_wkt = self.in_ds.GetGCPProjection()
			if self.in_srs_wkt:
				self.in_srs = osr.SpatialReference()
				self.in_srs.ImportFromWkt(self.in_srs_wkt)
			#elif self.options.profile != 'raster':
			#	self.error("There is no spatial reference system info included in the input file.","You should run gdal2tiles with --s_srs EPSG:XXXX or similar.")

		# Spatial Reference System of tiles
		
		self.out_srs = osr.SpatialReference()

		self.out_srs.ImportFromEPSG(3857)
		
		# Are the reference systems the same? Reproject if necessary.

		self.out_ds = None
		
		if (self.in_ds.GetGeoTransform() == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)) and (self.in_ds.GetGCPCount() == 0):
			self.error("There is no georeference - neither affine transformation (worldfile) nor GCPs. You can generate only 'raster' profile tiles.",
			"Either gdal2tiles with parameter -p 'raster' or use another GIS software for georeference e.g. gdal_transform -gcp / -a_ullr / -a_srs")			
		
		
		if not self.out_ds:
			self.out_ds = self.in_ds

		#
		# Here we should have a raster (out_ds) in the correct Spatial Reference system
		#

		# Get alpha band (either directly or from NODATA value)
		self.alphaband = self.out_ds.GetRasterBand(1).GetMaskBand()

		self.out_gt = self.out_ds.GetGeoTransform()
			
		#originX, originY = self.out_gt[0], self.out_gt[3]
		#pixelSize = self.out_gt[1] # = self.out_gt[5]
		
		# Test the size of the pixel
		
		# MAPTILER - COMMENTED
		#if self.out_gt[1] != (-1 * self.out_gt[5]) and self.options.profile != 'raster':
			# TODO: Process corectly coordinates with are have swichted Y axis (display in OpenLayers too)
			#self.error("Size of the pixel in the output differ for X and Y axes.")
			
		# Report error in case rotation/skew is in geotransform (possible only in 'raster' profile)
		if (self.out_gt[2], self.out_gt[4]) != (0,0):
			self.error("Georeference of the raster contains rotation or skew. Such raster is not supported. Please use gdalwarp first.")
			# TODO: Do the warping in this case automaticaly

		#
		# Here we expect: pixel is square, no rotation on the raster
		#

		# Output Bounds - coordinates in the output SRS
		self.ominx = self.out_gt[0]
		self.omaxx = self.out_gt[0]+self.out_ds.RasterXSize*self.out_gt[1]
		self.omaxy = self.out_gt[3]
		self.ominy = self.out_gt[3]-self.out_ds.RasterYSize*self.out_gt[1]
		# Note: maybe round(x, 14) to avoid the gdal_translate behaviour, when 0 becomes -1e-15

		if self.options.verbose:
			print("Bounds (output srs):", round(self.ominx, 13), self.ominy, self.omaxx, self.omaxy)

		#
		# Calculating ranges for tiles in different zoom levels
		#

		self.mercator = GlobalMercator() # from globalmaptiles.py
		
		# Function which generates SWNE in LatLong for given tile
		self.tileswne = self.mercator.TileLatLonBounds
		
		# Generate table with min max tile coordinates for all zoomlevels
		self.tminmax = list(range(0,32))
		for tz in range(0, 32):
			tminx, tminy = self.mercator.MetersToTile( self.ominx, self.ominy, tz )
			tmaxx, tmaxy = self.mercator.MetersToTile( self.omaxx, self.omaxy, tz )
			# crop tiles extending world limits (+-180,+-90)
			tminx, tminy = max(0, tminx), max(0, tminy)
			tmaxx, tmaxy = min(2**tz-1, tmaxx), min(2**tz-1, tmaxy)
			self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

		# TODO: Maps crossing 180E (Alaska?)

		# Get the minimal zoom level (map covers area equivalent to one tile) 
		if self.tminz == None:
			self.tminz = self.mercator.ZoomForPixelSize( self.out_gt[1] * max( self.out_ds.RasterXSize, self.out_ds.RasterYSize) / float(self.tilesize) )

		# Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
		if self.tmaxz == None:
			self.tmaxz = self.mercator.ZoomForPixelSize( self.out_gt[1] )
		
		if self.options.verbose:
			print("Bounds (latlong):", self.mercator.MetersToLatLon( self.ominx, self.ominy), self.mercator.MetersToLatLon( self.omaxx, self.omaxy))
			print('MinZoomLevel:', self.tminz)
			print("MaxZoomLevel:", self.tmaxz, "(", self.mercator.Resolution( self.tmaxz ),")")



	# -------------------------------------------------------------------------
	def wms_fetch_tiles(self, tz):

		# Set the bounds
		tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
		tilebands = 4
		ds = self.out_ds		
		tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
		ti = 0

		for ty in range(tmaxy, tminy-1, -1):
			for tx in range(tminx, tmaxx+1):

				if self.stopped:
					break
				ti += 1

				tilefilename = os.path.join(self.output, str(tz), str(tx), "%s.%s" % ((2**tz - 1) - ty, self.tileext))
				if self.options.verbose:
					print(ti,'/',tcount, tilefilename)

				# Create directories for the tile
				if not os.path.exists(os.path.dirname(tilefilename)):
					os.makedirs(os.path.dirname(tilefilename))

				if os.path.exists(tilefilename) and self.options.resume:
					if self.options.verbose:
						print("Tile generation skiped because of --resume")
					else:
						self.progressbar( ti / float(tcount) )
					continue

				(miny, minx, maxy, maxx) = self.mercator.TileLatLonBounds(tx, ty, tz)

				dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, 3)
				
				url = self.options.wms + "&version=1.1.1&request=GetMap&srs=EPSG:4326&width=256&height=256&format=image/png&styles=&bbox=" + str(minx) + "," + str(miny) + "," + str(maxx) + "," + str(maxy)
				res = urllib.urlretrieve(url, tilefilename)
				if not(res[1].type=="image/png"):
					print res[1].type, open(tilefilename).read()
					os.remove(tilefilename)

				if os.path.exists(tilefilename):
					tmp_in = gdal.Open(tilefilename, gdal.GA_ReadOnly)
					tmp_out = self.mem_drv.Create('', self.tilesize, self.tilesize, 3)
					self.pct2rgb(tmp_in, tmp_out)
					data = tmp_out.ReadRaster(0, 0, self.tilesize, self.tilesize, self.tilesize, self.tilesize, band_list=list(range(1,4)))
					dstile.WriteRaster(0, 0, self.tilesize, self.tilesize, data, band_list=list(range(1,4)))
					tmp_in = None
					tmp_out = None

					self.write_tile_png_palette(dstile, tilefilename)				
						
				if(not(s3['bucket']==None)):
					k = Key(s3['bucket'])
					k.key = tilefilename
					k.set_metadata("Content-Type", "image/png")
					k.set_metadata("Cache-Control", "max-age=2592000")
					k.set_contents_from_filename(tilefilename)
					
				del dstile
					
				if not self.options.verbose:
					self.progressbar( ti / float(tcount) )



	# -------------------------------------------------------------------------
	def write_tiles(self, tz):
		"""Generation of the base tiles (the lowest in the pyramid) directly from the input raster"""
		
		# Set the bounds
		tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]

		ds = self.out_ds
		if self.options.long or self.options.short:
			tilebands = 1
		else:
			tilebands = 4
		querysize = self.querysize
		
		
		#print tminx, tminy, tmaxx, tmaxy
		tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
		#print tcount
		ti = 0
		
		#SARSOFT
		for ty in range(tmaxy, tminy-1, -1):
			for tx in range(tminx, tmaxx+1):
				if self.stopped:
					break
				ti += 1

				tilefilename = os.path.join(self.output, str(tz), str(tx), "%s.%s" % ((2**tz - 1) - ty, self.tileext))
				if self.options.verbose:
					print(ti,'/',tcount, tilefilename) #, "( TileMapService: z / x / y )"

				# Create directories for the tile
				if not os.path.exists(os.path.dirname(tilefilename)):
					os.makedirs(os.path.dirname(tilefilename))
					
				if os.path.exists(tilefilename):
					if(self.options.resume):
						if self.options.verbose:
							print("Tile generation skiped because of --resume")
						else:
							self.progressbar( ti / float(tcount) )
						continue

				b = self.mercator.TileBounds(tx, ty, tz)

				#print "\tgdalwarp -ts 256 256 -te %s %s %s %s %s %s_%s_%s.tif" % ( b[0], b[1], b[2], b[3], "tiles.vrt", tz, tx, ty)

				# Don't scale up by nearest neighbour, better change the querysize
				# to the native resolution (and return smaller query tile) for scaling
				
				rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1])
				nativesize = wb[0]+wb[2] # Pixel size in the raster covering query geo extent

				# TODO find a way to read nodata values, scale raster down appropriately					
				# Tile bounds in raster coordinates for ReadRaster query
				rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1], querysize=querysize)

				rx, ry, rxsize, rysize = rb
				wx, wy, wxsize, wysize = wb
				
				mb = ds.GetRasterBand(1).GetMaskBand()
#				print "read:",rx,",",ry,"+",rxsize,",",rysize
#				print "write:",wx,",",wy,"+",wxsize,",",wysize

						
#					b = (max(b[0], self.ominx), max(b[1], self.ominy), min(b[2], self.omaxx), min(b[3], self.omaxy))
#					(rx2, ry2, rxsize2, rysize2), (wx2, wy2, wxsize2, wysize2) = self.geo_query(ds, b[0], b[3], b[2], b[1], querysize=querysize)
					

				if wxsize == 0 or wysize == 0:
					continue	
				# Query is in 'nearest neighbour' but can be bigger in then the tilesize
				# We scale down the query to the tilesize by supplied algorithm.

				bands = 3
				if self.options.grayscale:
					bands=1
				
				# Tile dataset in memory
				data = ds.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize, band_list=list(range(1,bands+1)))
				alpha = self.alphaband.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize, gdal.GDT_Int32)
				if self.options.grayscale:
					alpha = struct.pack('i'*wxsize*wysize, *([255] * (wxsize*wysize)))
				if self.options.long:
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_UInt16)
					dsquery = self.mem_drv.Create('', querysize, querysize, tilebands, gdal.GDT_UInt16)
				elif self.options.short:
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_Byte)
					dsquery = self.mem_drv.Create('', querysize, querysize, tilebands, gdal.GDT_Byte)
					# enable below line for canopy shading, to prefill unwritten areas with 255 rather than 0
					dsquery.GetRasterBand(1).WriteRaster(0, 0, querysize, querysize, struct.pack('i' * 1 * querysize * querysize, *([255] * querysize * querysize)), buf_type=gdal.GDT_Int32)
				else:
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)
					dsquery = self.mem_drv.Create('', querysize, querysize, tilebands)
					
				total = sum(struct.unpack('i' * 1 * wxsize * wysize, alpha))
				
				allalpha = (total == 255*querysize*querysize)
				if (not(allalpha) and not(self.options.crop)) or (self.options.crop and (tx == tminx or tx == tmaxx or ty == tminy or ty == tmaxy)):
					if self.options.fetch and not(os.path.exists(tilefilename)):
						res = urllib.urlretrieve(self.options.fetch + str(tz) + "/" + str(tx) + "/%s.%s" % ((2**tz - 1) - ty, self.tileext), tilefilename)
						if not(res[1].type=="image/png"):
							os.remove(tilefilename)
					if os.path.exists(tilefilename):
						tmp_in = gdal.Open(tilefilename, gdal.GA_ReadOnly)
						if self.options.long:
							tmp_out = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_UInt16)
						elif self.options.short:
							tmp_out = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_Byte)
						else:
							tmp_out = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)
						self.pct2rgb(tmp_in, tmp_out)
						tdata = tmp_out.ReadRaster(0, 0, self.tilesize, self.tilesize, querysize, querysize, band_list=list(range(1,tmp_out.RasterCount+1)))
						dsquery.WriteRaster(0, 0, querysize, querysize, tdata, band_list=list(range(1,tmp_out.RasterCount+1)))
						tmp_in = None
						tmp_out = None
				if allalpha or self.options.crop:
					if self.options.verbose:
						print "Tile at ", tx, ", ", ((2**tz - 1) - ty), " is 100% alpha"
					if self.options.long:
						dsquery.GetRasterBand(1).WriteRaster(wx, wy, wxsize, wysize, data, buf_type=gdal.GDT_Float32)
					elif self.options.short:
						dsquery.GetRasterBand(1).WriteRaster(wx, wy, wxsize, wysize, data, buf_type=gdal.GDT_Byte)
					else:
						if self.options.grayscale:
							dsquery.GetRasterBand(1).WriteRaster(wx, wy, wxsize, wysize, data)
							dsquery.GetRasterBand(2).WriteRaster(wx, wy, wxsize, wysize, data)
							dsquery.GetRasterBand(3).WriteRaster(wx, wy, wxsize, wysize, data)
						else:
							dsquery.WriteRaster(wx, wy, wxsize, wysize, data, band_list=list(range(1,4)))
						dsquery.WriteRaster(wx, wy, wxsize, wysize, alpha, band_list=[tilebands], buf_type=gdal.GDT_Int32)
				else:
					if self.options.verbose:
						print "Tile at ", tx, ", ", ((2**tz - 1) - ty), " is incomplete -- copying pixel by pixel"
					xratio = float(wxsize)/float(rxsize)
					yratio = float(wysize)/float(rysize)
					my = 0
					for ypixel in range(ry, ry+rysize):
						alpha = struct.unpack('i'*rxsize, mb.ReadRaster(rx, ypixel, rxsize, 1, rxsize, 1, gdal.GDT_Int32))						
						for xpixel in range(rx, rx+rxsize):
							pxalpha = alpha[xpixel-rx]
							if pxalpha > 0:
								origx = wx + int(round((xpixel-rx)*xratio))
								origy = wy + int(round((ypixel-ry)*yratio))
								runx = int(math.ceil(xratio))
								runy = int(math.ceil(yratio))
								if origx+runx > querysize:
									runx = querysize-origx
								if origy+runy > querysize:
									runy = querysize-origy
								if self.options.long or self.options.short:
									pxdata = ds.ReadRaster(xpixel, ypixel, 1, 1, 1, 1, band_list=list(range(1,2)))
									dsquery.GetRasterBand(1).WriteRaster(origx, origy, runx, runy, pxdata, 1, 1, buf_type=gdal.GDT_Byte)
								elif self.options.grayscale:
									pxdata = ds.ReadRaster(xpixel, ypixel, 1, 1, 1, 1, band_list=list(range(1,2)))
									dsquery.GetRasterBand(1).WriteRaster(origx, origy, runx, runy, pxdata, 1, 1)
									dsquery.GetRasterBand(2).WriteRaster(origx, origy, runx, runy, pxdata, 1, 1)
									dsquery.GetRasterBand(3).WriteRaster(origx, origy, runx, runy, pxdata, 1, 1)
									dsquery.GetRasterBand(4).WriteRaster(origx, origy, runx, runy, struct.pack('i', 255), 1, 1)
								else:
									pxdata = ds.ReadRaster(xpixel, ypixel, 1, 1, 1, 1, band_list=list(range(1,5)))
									dsquery.WriteRaster(origx, origy, runx, runy, pxdata, 1, 1, band_list=list(range(1,5)))
				self.scale_query_to_tile(dsquery, dstile, tilefilename)
				del dsquery
				
				del data
				if self.options.long or self.options.short:
					self.out_drv.CreateCopy(tilefilename, dstile, 0, ["ZLEVEL=9"])
				else:
					self.write_tile_png_palette(dstile, tilefilename)
						
				if(not(self.options.bucket==None)):
					k = Key(s3['bucket'])
					k.key = tilefilename
					k.set_metadata("Content-Type", "image/png")
					k.set_metadata("Cache-Control", "max-age=2592000")
					k.set_contents_from_filename(tilefilename)

				del dstile
					
				if not self.options.verbose:
					self.progressbar( ti / float(tcount) )
		
	# -------------------------------------------------------------------------
	def composite_tiles(self, tz):
		"""Generation of the overview tiles (higher in the pyramid) based on existing tiles"""
		
		tilebands = 4
		if self.options.long or self.options.short:
			tilebands = 1
		
		ti = 0
		
		tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
		tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
		for ty in range(tmaxy, tminy-1, -1): #range(tminy, tmaxy+1):
			for tx in range(tminx, tmaxx+1):
				if self.stopped:
					break
					
				ti += 1
				tilefilename = os.path.join( self.output, str(tz), str(tx), "%s.%s" % ((2**tz - 1) - ty, self.tileext) )

				if self.options.verbose:
					print(ti,'/',tcount, tilefilename) #, "( TileMapService: z / x / y )"
				
				if self.options.resume and os.path.exists(tilefilename):
					if self.options.verbose:
						print("Tile generation skiped because of --resume")
					else:
						self.progressbar( ti / float(tcount) )
					continue

				# Create directories for the tile
				if not os.path.exists(os.path.dirname(tilefilename)):
					os.makedirs(os.path.dirname(tilefilename))

				if self.options.long:
					dsquery = self.mem_drv.Create('', 2*self.tilesize, 2*self.tilesize, tilebands, gdal.GDT_UInt16)
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_UInt16)
				elif self.options.short:
					dsquery = self.mem_drv.Create('', 2*self.tilesize, 2*self.tilesize, tilebands, gdal.GDT_Byte)
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_Byte)
					# enable below line for canopy shading, to prefill unwritten areas with 255 rather than 0
					dstile.GetRasterBand(1).WriteRaster(0, 0, self.tilesize, self.tilesize, struct.pack('i' * 1 * self.tilesize * self.tilesize, *([255] * self.tilesize * self.tilesize)), buf_type=gdal.GDT_Int32)
					dsquery.GetRasterBand(1).WriteRaster(0, 0, 2*self.tilesize, 2*self.tilesize, struct.pack('i' * 4 * self.tilesize * self.tilesize, *([255] * 4 * self.tilesize * self.tilesize)), buf_type=gdal.GDT_Int32)
				else:
					dsquery = self.mem_drv.Create('', 2*self.tilesize, 2*self.tilesize, tilebands)
					dstile = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)

				# TODO: Implement more clever walking on the tiles with cache functionality
				# probably walk should start with reading of four tiles from top left corner
				# Hilbert curve...

				childrenexist = False
				for y in range(2*ty,2*ty+2):
					for x in range(2*tx,2*tx+2):
						child_tile_filename = os.path.join( self.output, str(tz+1), str(x), "%s.%s" % ((2**(tz+1) - 1) - y, self.tileext))
						if os.path.exists(child_tile_filename):
							childrenexist = True

				# Read the tiles and write them to query window
				children = []
				colors = 0
				if childrenexist or self.options.overview:
					for y in range(2*ty,2*ty+2):
						for x in range(2*tx,2*tx+2):
							cy = (2**(tz+1) - 1) - y
							cx = x
							cz = tz + 1
							child_tile_filename = os.path.join( self.output, str(cz), str(cx), "%s.%s" % (cy, self.tileext))
							if self.options.fetch:
								allAlpha = False
								ds = self.out_ds
								b = self.mercator.TileBounds(x, y, cz)
								rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1])
								nativesize = wb[0]+wb[2] # Pixel size in the raster covering query geo extent
								rb, wb = self.geo_query( ds, b[0], b[3], b[2], b[1], querysize=self.querysize)					
								rx, ry, rxsize, rysize = rb
								mb = ds.GetRasterBand(1).GetMaskBand()
								alpha = self.alphaband.ReadRaster(rx, ry, rxsize, rysize, 256, 256, gdal.GDT_Int32)					
								total = sum(struct.unpack('i' * 1 * 256 * 256, alpha))
								if total == 255*256*256:
									allAlpha = True
								if not os.path.exists(os.path.dirname(child_tile_filename)):
									os.makedirs(os.path.dirname(child_tile_filename))
								if((not(self.options.bucket==None) and allAlpha==False) or not(os.path.exists(child_tile_filename))):
									res = urllib.urlretrieve(self.options.fetch + str(cz) + "/" + str(cx) + "/%s.%s" % (cy, self.tileext), child_tile_filename)
									if not(res[1].type=="image/png"):
										os.remove(child_tile_filename)
							if(os.path.exists(child_tile_filename)):
								dsquerytile = gdal.Open(child_tile_filename)
								if (ty==0 and y==1) or (ty!=0 and (y % (2*ty)) != 0):
									tileposy = 0
								else:
									tileposy = self.tilesize
								if tx:
									tileposx = x % (2*tx) * self.tilesize
								elif tx==0 and x==1:
									tileposx = self.tilesize
								else:
									tileposx = 0
								if self.options.long:
									dstmp = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_UInt16)
								elif self.options.short:
									dstmp = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands, gdal.GDT_Byte)
								else:
									dstmp = self.mem_drv.Create('', self.tilesize, self.tilesize, tilebands)
								colors = max(colors, self.pct2rgb(dsquerytile, dstmp))
								dsquery.WriteRaster( tileposx, tileposy, self.tilesize, self.tilesize,
									dstmp.ReadRaster(0,0,self.tilesize,self.tilesize),
									band_list=list(range(1,tilebands+1)))
								children.append( [x, y, tz+1] )
								dstmp = None

				if len(children) == 0:
					continue

				self.scale_query_to_tile(dsquery, dstile, tilefilename)

				if dstile.RasterCount == 4 and False:
					alpha = struct.unpack('i' * dstile.RasterXSize * dstile.RasterYSize, dstile.GetRasterBand(4).ReadRaster(0, 0, dstile.RasterXSize, dstile.RasterYSize, dstile.RasterXSize, dstile.RasterYSize, gdal.GDT_Int32))
					for x in range(0, dstile.RasterXSize):
						for y in range(0, dstile.RasterYSize):
							pxvalue = alpha[y*dstile.RasterYSize + x]
							if pxvalue > 0 and pxvalue < 255:
								dstile.GetRasterBand(4).WriteRaster(x, y, 1, 1, struct.pack('i', 0))

				if self.options.long or self.options.short:
					self.out_drv.CreateCopy(tilefilename, dstile, 0, ["ZLEVEL=9"])
				else:
					self.write_tile_png_palette(dstile, tilefilename)				

				if(not(self.options.bucket==None)):
					k = Key(s3['bucket'])
					k.key = tilefilename
					k.set_metadata("Content-Type", "image/png")
					k.set_metadata("Cache-Control", "max-age=2592000")
					k.set_contents_from_filename(tilefilename)

	# -------------------------------------------------------------------------
	def geo_query(self, ds, ulx, uly, lrx, lry, querysize = 0):
		"""For given dataset and query in cartographic coordinates
		returns parameters for ReadRaster() in raster coordinates and
		x/y shifts (for border tiles). If the querysize is not given, the
		extent is returned in the native resolution of dataset ds."""

		geotran = ds.GetGeoTransform()
		rx= int((ulx - geotran[0]) / geotran[1] + 0.001)
		ry= int((uly - geotran[3]) / geotran[5] + 0.001)
		rxsize= int((lrx - ulx) / geotran[1] + 0.5)
		rysize= int((lry - uly) / geotran[5] + 0.5)

		if not querysize:
			wxsize, wysize = rxsize, rysize
		else:
			wxsize, wysize = querysize, querysize

		# Coordinates should not go out of the bounds of the raster
		wx = 0
		if rx < 0:
			rxshift = abs(rx)
			wx = int( wxsize * (float(rxshift) / rxsize) )
			wxsize = wxsize - wx
			rxsize = rxsize - int( rxsize * (float(rxshift) / rxsize) )
			rx = 0
		if rx+rxsize > ds.RasterXSize:
			wxsize = int( wxsize * (float(ds.RasterXSize - rx) / rxsize) )
			rxsize = ds.RasterXSize - rx

		wy = 0
		if ry < 0:
			ryshift = abs(ry)
			wy = int( wysize * (float(ryshift) / rysize) )
			wysize = wysize - wy
			rysize = rysize - int( rysize * (float(ryshift) / rysize) )
			ry = 0
		if ry+rysize > ds.RasterYSize:
			wysize = int( wysize * (float(ds.RasterYSize - ry) / rysize) )
			rysize = ds.RasterYSize - ry

		return (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)

	# -------------------------------------------------------------------------
	def pct2rgb(self, source, dest):
		try:
			idband = source.GetRasterBand(1)
			ct = idband.GetRasterColorTable()
			if source.RasterCount == 1 and ct is not None:
				tmp_ds = self.mem_drv.Create('', source.RasterXSize, source.RasterYSize, 4)
				for ypixel in range(0, source.RasterYSize):
					indexes = struct.unpack('i' * source.RasterXSize, idband.ReadRaster(0, ypixel, source.RasterXSize, 1, source.RasterXSize, 1, gdal.GDT_Int32))
					r = [0 for x in range(0, source.RasterXSize)]
					g = [0 for x in range(0, source.RasterXSize)]
					b = [0 for x in range(0, source.RasterXSize)]
					a = [0 for x in range(0, source.RasterXSize)]
					for xpixel in range(0, source.RasterXSize):
						index = indexes[xpixel]
						try:
							r[xpixel] = ct.GetColorEntry(index)[0]
							g[xpixel] = ct.GetColorEntry(index)[1]
							b[xpixel] = ct.GetColorEntry(index)[2]
							a[xpixel] = ct.GetColorEntry(index)[3]
						except:
							pass
					dest.GetRasterBand(1).WriteRaster(0, ypixel, source.RasterXSize, 1, struct.pack('i' * source.RasterXSize, *r), source.RasterXSize, 1, gdal.GDT_Int32)
					dest.GetRasterBand(2).WriteRaster(0, ypixel, source.RasterXSize, 1, struct.pack('i' * source.RasterXSize, *g), source.RasterXSize, 1, gdal.GDT_Int32)
					dest.GetRasterBand(3).WriteRaster(0, ypixel, source.RasterXSize, 1, struct.pack('i' * source.RasterXSize, *b), source.RasterXSize, 1, gdal.GDT_Int32)
					dest.GetRasterBand(4).WriteRaster(0, ypixel, source.RasterXSize, 1, struct.pack('i' * source.RasterXSize, *a), source.RasterXSize, 1, gdal.GDT_Int32)
				return ct.GetCount()
			else:
				data = source.ReadRaster(0, 0, source.RasterXSize, source.RasterYSize)
				dest.WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, data)
				return 0
		except:
			print sys.exc_info()
			import traceback
			traceback.print_tb(sys.exc_info()[2])
			print "Error copying tile data"

	# -------------------------------------------------------------------------
 	def rgb2pct(self, source, dest, colors=0):
 		if colors==0:
 			colors=self.colors
		ct = gdal.ColorTable()
		s1 = source.GetRasterBand(1)
		s2 = source.GetRasterBand(1)
		s3 = source.GetRasterBand(1)
		if source.RasterCount >= 3:
			s2 = source.GetRasterBand(2)
			s3 = source.GetRasterBand(3)
		err = gdal.ComputeMedianCutPCT(s1, s2, s3, colors, ct)
		dest.GetRasterBand(1).SetRasterColorTable(ct)
		gdal.DitherRGB2PCT(s1, s2, s3, dest.GetRasterBand(1), ct)
		
	def rgb2rgb(self, source, dest):
		if dest.RasterCount == 4 and source.RasterCount == 3:
			dest.GetRasterBand(1).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(1).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
			dest.GetRasterBand(2).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(2).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
			dest.GetRasterBand(3).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(3).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))

			a = [255 for x in range(0, source.RasterXSize*source.RasterYSize)]
			dest.GetRasterBand(4).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, struct.pack('i' * source.RasterXSize * source.RasterYSize, *a), source.RasterXSize, source.RasterYSize, gdal.GDT_Int32)
		elif source.RasterCount == 2 and dest.RasterCount == 4:
			dest.GetRasterBand(1).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(1).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
			dest.GetRasterBand(2).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(1).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
			dest.GetRasterBand(3).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(1).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
			dest.GetRasterBand(4).WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, source.GetRasterBand(2).ReadRaster(0, 0, source.RasterXSize, source.RasterYSize))
		else:
			data = source.ReadRaster(0, 0, source.RasterXSize, source.RasterYSize)
			dest.WriteRaster(0, 0, source.RasterXSize, source.RasterYSize, data)
	
	def rgb2pil(self, source):
		data = []
		for ypixel in range(0, source.RasterYSize):
			r = struct.unpack('i' * source.RasterXSize, source.GetRasterBand(1).ReadRaster(0, ypixel, source.RasterXSize, 1, source.RasterXSize, 1, gdal.GDT_Int32))
			g = struct.unpack('i' * source.RasterXSize, source.GetRasterBand(2).ReadRaster(0, ypixel, source.RasterXSize, 1, source.RasterXSize, 1, gdal.GDT_Int32))
			b = struct.unpack('i' * source.RasterXSize, source.GetRasterBand(3).ReadRaster(0, ypixel, source.RasterXSize, 1, source.RasterXSize, 1, gdal.GDT_Int32))
			a = struct.unpack('i' * source.RasterXSize, source.GetRasterBand(4).ReadRaster(0, ypixel, source.RasterXSize, 1, source.RasterXSize, 1, gdal.GDT_Int32))
			for xpixel in range(0, source.RasterXSize):
				data.append((r[xpixel], g[xpixel], b[xpixel], a[xpixel]))
		im = Image.new("RGBA", (256,256))
		im.putdata(data)
		return im
	
	def write_tile_png_palette(self, dstile, tilefilename):
		if palette_method == "pillow":
			# Pillow's RGB->palette conversion has problems where large white areas
			# (like islands) get colored with the blue of surrounding water.
			self.rgb2pil(dstile).convert('P', palette=Image.ADAPTIVE).save(tilefilename, 'PNG')
			return

		# first write tile to disk, uncompressed
		self.out_drv.CreateCopy(tilefilename, dstile, 0, ["ZLEVEL=1"])
		if palette_method == "pngquant":
	        # Call pngquant https://pngquant.org/
	        # pngquant's quality=70 seems to match image size of Pillow's conversion, above
			subprocess.call([
				"pngquant",
				"--quality", "70",
				"--force",                   # allow overwrite
				"--output", tilefilename,    # overwrite source file
				tilefilename,
			])

		tilesidecarxmlfilename = tilefilename + ".aux.xml"
		if os.path.exists(tilesidecarxmlfilename):
			os.remove(tilesidecarxmlfilename)


	def resample_aspect(self, source, dest):
		src_data = struct.unpack('B' * 512*512, source.GetRasterBand(1).ReadRaster(0, 0, 512, 512, 512, 512, gdal.GDT_Byte))
		dst_data = [0 for i in range(0,256*256)]
		for xpixel in range(0, 256):
			for ypixel in range(0, 256):
				samples = []
				samples.append(src_data[xpixel*2 + ypixel*2*512]*1.5)
				samples.append(src_data[xpixel*2+1 + ypixel*2*512]*1.5)
				samples.append(src_data[xpixel*2 + (ypixel*2+1)*512]*1.5)
				samples.append(src_data[xpixel*2+1 + (ypixel*2+1)*512]*1.5)
				x = y = 0
				for angle in samples:
					angle = angle * math.pi / 180
					x += math.cos(angle)
					y += math.sin(angle)
				average_angle = math.atan2(y, x) * 180 / math.pi
				if average_angle < 0:
					average_angle = average_angle + 360
				if average_angle > 360:
					average_angle = average_angle - 360
				dst_data[xpixel + ypixel*256] = int(max(min(average_angle, 359), 0) / 1.5)
		dest.GetRasterBand(1).WriteRaster(0, 0, 256, 256, struct.pack('B' * 256 * 256, *dst_data), 256, 256, gdal.GDT_Byte)

	def resample_max(self, source, dest):
		src_data = struct.unpack('H' * 512*512, source.GetRasterBand(1).ReadRaster(0, 0, 512, 512, 512, 512, gdal.GDT_UInt16))
		dst_data = [0 for i in range(0,256*256)]
		for xpixel in range(0, 256):
			for ypixel in range(0, 256):
				val = max(src_data[xpixel*2 + ypixel*2*512], src_data[xpixel*2+1 + ypixel*2*512], src_data[xpixel*2 + (ypixel*2+1)*512], src_data[xpixel*2+1 + (ypixel*2+1)*512])
				dst_data[xpixel + ypixel*256] = val
		dest.GetRasterBand(1).WriteRaster(0, 0, 256, 256, struct.pack('H' * 256 * 256, *dst_data), 256, 256, gdal.GDT_UInt16)

	def resample_min(self, source, dest):
		src_data = struct.unpack('H' * 512*512, source.GetRasterBand(1).ReadRaster(0, 0, 512, 512, 512, 512, gdal.GDT_UInt16))
		dst_data = [0 for i in range(0,256*256)]
		for xpixel in range(0, 256):
			for ypixel in range(0, 256):
				val = min(src_data[xpixel*2 + ypixel*2*512], src_data[xpixel*2+1 + ypixel*2*512], src_data[xpixel*2 + (ypixel*2+1)*512], src_data[xpixel*2+1 + (ypixel*2+1)*512])
				dst_data[xpixel + ypixel*256] = val
		dest.GetRasterBand(1).WriteRaster(0, 0, 256, 256, struct.pack('H' * 256 * 256, *dst_data), 256, 256, gdal.GDT_UInt16)

	# -------------------------------------------------------------------------
	def scale_query_to_tile(self, dsquery, dstile, tilefilename='', resampling=None):
		"""Scales down query dataset to the tile dataset"""

		querysize = dsquery.RasterXSize
		tilesize = dstile.RasterXSize
		tilebands = dstile.RasterCount
		
		if resampling == None:
			resampling = self.resampling;

		if self.options.resampling == 'average':

			# Function: gdal.RegenerateOverview()
			for i in range(1,tilebands+1):
				res = gdal.RegenerateOverview( dsquery.GetRasterBand(i),
					dstile.GetRasterBand(i), 'average' )
				if res != 0:
				    self.error("RegenerateOverview() failed on %s, error %d" % (tilefilename, res))
		
		elif self.options.resampling == 'aspect':
		    self.resample_aspect(dsquery, dstile)

		elif self.options.resampling == 'max':
		    self.resample_max(dsquery, dstile)

		elif self.options.resampling == 'min':
		    self.resample_min(dsquery, dstile)

		elif tilesize == querysize:
			self.rgb2rgb(dsquery, dstile)
		else:

			# Other algorithms are implemented by gdal.ReprojectImage().
			dsquery.SetGeoTransform( (0.0, tilesize / float(querysize), 0.0, 0.0, 0.0, tilesize / float(querysize)) )
			dstile.SetGeoTransform( (0.0, 1.0, 0.0, 0.0, 0.0, 1.0) )

			res = gdal.ReprojectImage(dsquery, dstile, None, None, resampling)
			if res != 0:
			    self.error("ReprojectImage() failed on %s, error %d" % (tilefilename, res))
			

# =============================================================================
# =============================================================================
# =============================================================================

if __name__=='__main__':
	argv = gdal.GeneralCmdLineProcessor( sys.argv )
	if argv:
		gdal2tiles = GDAL2Tiles( argv[1:] )
		gdal2tiles.process()

