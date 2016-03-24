#!/usr/bin/python

""" 
Prepare MAIA data for Gaia Science Alerts:
-- run astrometry.net on the *REDUCED* MAIA images
-- run SExtractor
-- create new event
-- upload to calibration servers

To run astrometry.net you need to tell backend.conf where your index files are (see line 55).

By using this code you assume any and all responsibility for stress, hair loss or any other effects of usage!

Things I am aware of:
-- interference pattern in R band images
-- U flux is really low -- but that's apparently due to problems with coating on the CCD window
-- MAIA has a special set of filters U, G, R: U~SDSS u, G~SDSS g, R~SDSS r+i; therefore use SDSS u, g, r combo! In case the object is outside SDSS footprint, APASS filters seem to be best choice (APASS g and APASS r, there is no good u match)
-- I don't really trust the automatic filter matching...USNO filters...scanned plates...really???

Some guidelines:
-- at least EXPTIME=200s for a ~17mag object. There's no exp. time calc for MAIA. Otherwise SN is really poor.
-- when choosing new targets to follow-up, (preferentially) choose those that are in SDSS footprint, then force the filters.
-- when preparing for the run make an ASCII tab separated table with following header format:
# OBJNAME IVO RA DEC URL
Make damn sure that OBJNAME is exactly the same as the OBJECT FITS header keyword. Will make your life a lot easier later.

-- corrected bug in single file/dir mode where all observations would have the same MJD
-- some readability improvements

17/01/2015
________________________________________________________________________
________________________________________________________________________

-- files from the laptop were accidentaly erased. Current incarnation is a copy of the server code

06/02/2015
 """

import os, sys, glob, subprocess
from astropy.io import fits
import astropy.time as aptime
import numpy as np

def REGfile(DIRNAME, RAEV, DECEV):
	regfile = open(wd+DIRNAME+'/'+DIRNAME+'.reg', 'w')
	regfile.write('# Region file format: DS9 version 4.1\nglobal color=red dashlist=8 3 width=1 font=\"helvetica 12 bold roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\ncircle('+RAEV+','+DECEV+',2\") # text={'+DIRNAME+'}')
	#~ regfile.write('RAEV+','DECEV+',2 \# text={'+DIRNAME+'}')
	regfile.close()

def mffMAIA(wd):
	""" Split MEFs and and fix'em. """
	filters={1:'ri', 2:'g', 3:'u'}
	
	flist = glob.glob(wd+"c*OBJ.fits")
	print "FLIST:", [n.split("/")[-1] for n in flist], "\n" #flist.split["/"][-1]
	print "WD:", wd
	print os.getcwd()
	print flist
	
	for f in flist:
		for EXT in xrange(1,4):
			hdu=fits.PrimaryHDU(fits.getdata(f, EXT), header=fits.getheader(f,0))
			"""conservative value for gain=4e-"""
			hdu.header.set('RON', 4, '[e-] Readout noise')
			hdu.header.set('FILTER', filters[EXT], 'MAIA filter')
			try:
				exthdu=fits.getheader(f,EXT)
				MJD = aptime.Time(exthdu['DATE-AVG'], format='isot', scale='ut1').mjd # Better to have the middle of the observation. Note that BJD is not calculated at the moment but for objects that don't vary rapidly this is not that important...
			except KeyError:
				dateavg = fits.getheader(f,0)['DATE-AVG']
				MJD = aptime.Time(dateavg, format='isot', scale='ut1').mjd # When DATE-AVG is missing from the extension(s)
				print "--------------DATE-AVG-----From the main header---------"
				print MJD
			
			hdu.header.set('MJD', MJD, 'Modified Julian Date')
			EXPTIME = hdu.header['EXPTIME']
			
			print "Written:", f.split('.')[-2]+'_EXT'+repr(EXT)+'.fits'
			print "This is f:",f
			print "This is f.split[-2]:",f.split('.')[-2]
			
			hdu.writeto(f.split('.')[-2]+'_EXT'+repr(EXT)+'.fits', clobber=True)

	bashcode=open("backend.cfg", "w")
	bashcode.write("inparallel\nminwidth 0.1\nmaxwidth 1\ndepths 10 20 30 40 50 60 70 80 90 100\ncpulimit 10\nadd_path /home/lovro/Programi/astrometry.net-0.50/index_files\nautoindex")
	#~ bashcode.write("inparallel\nminwidth 0.1\nmaxwidth 1\ndepths 10 20 30 40 50 60 70 80 90 100\ncpulimit 30\nadd_path /home/lovro/Desktop/astrometry.net\nautoindex")
	bashcode.close()

	bashcode=open("default.conv", "w")
	bashcode.write("CONV NORM\n1 2 1\n2 4 2\n1 2 1")
	bashcode.close()

def astrometry(wd):
	""" Run astrometry.net code. """

	print "------------Running astrometry.net-----------------"
	
	flist=glob.glob(wd+"c*OBJ_EXT*.fits")
	print flist
	
	RA=fits.getheader(flist[0])['RA']
	DEC=fits.getheader(flist[0])['DEC']
	
	bashcode=open("bashcode.sh", "w")
	bashcode.write("#! /bin/bash\ncat file.list | \\\nwhile read list\ndo\nsolve-field ${list} -b backend.cfg -p -O -L 8 -H 15 -u aw -3 "+RA+" -4 "+DEC+" -5 0.3 -z 4 -S none -M none -R none -B none -U none\ndone")
	bashcode.close()
	
	fout = open("file.list", "w")
	
	for i in flist:
		fout.write("%s \n" % i)
	fout.close()

	subprocess.Popen(['sh', 'bashcode.sh']).wait()
		
def sex(wd):

	""" Run SExtractor (only after all files have been split and prepared!). """
	print "---------------Now starting SExtractor----------"
	
	""" TODO: Check exactly which coordinates are exported by SExtractor. """
	
	#~ default_param=open("default.param", "w")
	#~ default_param.write("MAG_AUTO\nMAGERR_AUTO\nALPHA_J2000\nDELTA_J2000")
	#~ default_param.close()
	#~ 
	#~ p=subprocess.Popen('sextractor -dd > default.sex', shell=True)
	#~ p.terminate()
	#~ 
	#~ flist=glob.glob(wd+"c*OBJ_EXT*.new")
	#~ 
	#~ for f in flist:
		#~ fout=f.split('.')[0]+'.cat'
		#~ subprocess.Popen(['sextractor',f,'-catalog_name',fout]).wait()
		#~ 
		#~ 
		#~ 
	
	
	
	""" Run SExtractor (only after all files have been split and prepared!). """
	print "---------------Now starting SExtractor----------"
		
	default_param=open("default.param", "w")
	default_param.write("MAG_AUTO\nMAGERR_AUTO\nALPHA_J2000\nDELTA_J2000")
	default_param.close()
	
	#~ p=subprocess.Popen('sextractor -dd > default.sex', shell=True)
	#~ p.terminate()
	
	flist=glob.glob(wd+"c*OBJ_EXT*.new")
	
	for f in flist:
		#~ fout=f.split('.')[-1]+'.cat'
		fout=f.replace('.new', '.cat')
		check=f.replace('.new', '.ap')

		#~ subprocess.Popen(['sextractor',f,'-catalog_name',fout, '-CHECKIMAGE_TYPE','APERTURES', '-CHECKIMAGE_NAME', check]).wait()
		
		#~ subprocess.Popen(['sextractor',f,'-catalog_name',fout, '-CHECKIMAGE_TYPE','APERTURES', '-CHECKIMAGE_NAME', check, '-DETECT_MINAREA', '4', '-DETECT_THRESH', '3', '-DEBLEND_NTHRESH','16', '-DEBLEND_MINCONT', '0.00005']).wait() # SUBARU version

		subprocess.Popen(['sextractor',f,'-catalog_name',fout, '-CHECKIMAGE_TYPE','APERTURES', '-CHECKIMAGE_NAME', check, '-PIXEL_SCALE', '0.276', '-GAIN_KEY', 'HIERARCH OGE_DET_OUT1_GAIN', '-DETECT_MINAREA', '10', '-DETECT_THRESH', '3', '-DEBLEND_NTHRESH','16', '-DEBLEND_MINCONT', '0.00005']).wait() # Testin version
		
		#~ subprocess.Popen(['sextractor',f,'-catalog_name',fout, '-CHECKIMAGE_TYPE','APERTURES', '-CHECKIMAGE_NAME', check, '-PIXEL_SCALE', '0.215', '-GAIN_KEY', 'HIERARCH OGE_DET_OUT1_GAIN', '-DETECT_MINAREA', '100', '-DETECT_THRESH', '4', '-DEBLEND_NTHRESH','8', '-DEBLEND_MINCONT', '0.5']).wait() # Testin version

		print check
		
def make_new_event(hashtag, ivo, RAEV, DECEV, URLEV):
	print "Inserting New Event"
	subprocess.Popen(['curl', '-F', hashtag, '-F', 'EventID='+ivo, '-F', 'ra='+RAEV, '-F', 'dec='+DECEV, '-F', 'url='+URLEV,  "http://gsaweb.ast.cam.ac.uk/followup/cgi/newevent"]).wait()
	
def upload(ivo, hashtag, wd, DIRNAME, sdss):
	filtersSDSS={'ri':'SDSS/r', 'g':'SDSS/g', 'u':'SDSS/u'}
	#~ filtersSDSS={1:'SDSS/r', 2:'SDSS/g', 3:'SDSS/u'}
	#~ filtersAPASS={'ri':'APASS/r', 'g':'APASS/g', 'u':'APASS/u'}
	
	print "DIRNAME:", DIRNAME
	print "wd:", wd
	
	flist=glob.glob(wd+"c*OBJ_EXT*.cat")
	#~ print "**********************************************************\n"
	#~ print "FLIST:", [n.split("/")[-1] for n in flist]
	
	print flist
	
	for f in flist:
		MJD = fits.getheader(f.replace(".cat",".new"))['MJD']
		EXPTIME = fits.getheader(f.replace(".cat",".new"))['EXPTIME']
		print "\nNow uploading:", f.split("/")[-1]
		print "F:",f
		print ivo, MJD, EXPTIME, DIRNAME
		if sdss == '1':
			filt=fits.getheader(f.replace('.cat', '.new'))['FILTER']
			filt=filtersSDSS[filt]
		else:
			filt = 'no'
		subprocess.Popen(['curl', '-F', 'sexCat=@'+f+';filename=test.cat', '-F','EventID='+ivo, '-F', 'MJD='+str(MJD), '-F', 'expTime='+str(EXPTIME), '-F', 'noPlot=1', '-F', 'forceFilter='+str(filt), '-F', 'outputFormat=json', '-F', 'matchDist=2', '-F', hashtag, "http://gsaweb.ast.cam.ac.uk/followup/cgi/upload"]).wait()

		""" Below are some options, if you want to do a dry run and (not) force the filters. Generally it seems that SDSS works the best, APASS is good if m~<17. There is no MAIA u band equivalent filter in APASS. """
		
		# '-F', 'dryRun=1',
		# 'forceFilter=no'
		# 'forceFilter='+str(filters[filt])
	
def cleanstuff():
	try:
		os.remove('backend.cfg')
		os.remove('default.conv')
		os.remove('file.list')
		os.remove('default.param')
		os.remove('bashcode.sh')
		os.remove('*axy*')
		os.remove(wd+'/*fits*')
	except OSError:
		pass

hashtag ="SUPERSECRET"	

""" Next seven lines only needed when running on one object/dir at a time, i.e. "diagnostic mode". """
#~ wd=sys.argv[1] # Directory to be processed
#~ ivo=sys.argv[2]
#~ DIRNAME=wd.split("/")[0]
#~ mff(wd)
#~ astrometry(wd)
#~ fout=sex(wd)
#~ upload(ivo, hashtag, wd, DIRNAME)
#~ cleanstuff()

""" Use stuff below if you're running from a list. """

wd=sys.argv[1] # Directory containing subdirs to be processed
objlist=np.genfromtxt(wd+"GSA_targets.list", names=True, dtype=None)

print "objlist", objlist

for event in objlist:
	DIRNAME=event['evname']
	ivo="ivo://"+DIRNAME
	RAEV=str(event['ra'])
	DECEV=str(event['dec'])
	sdss = str(event['sdss'])
	#~ sdss = '1'
	# URLEV=event['url']
	URLEV='http://gaia.ac.uk/selected-gaia-science-alerts#alerts'
	print "\n----------------> ivo:", ivo
	
	#~ make_new_event(hashtag, ivo, RAEV, DECEV, URLEV)
	for path in glob.glob(wd+DIRNAME+'/'):
		##~ MJD, EXPTIME=mff(path)
		mffMAIA(path)
		astrometry(path)
		fout=sex(path)
		REGfile(DIRNAME, RAEV, DECEV)
		upload(ivo, hashtag, path, DIRNAME, sdss)
		cleanstuff()
