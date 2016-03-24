#!/usr/bin/python
import numpy as np
from astropy.io import fits
import glob, sys, os
from select import select
from ds9 import *
import scipy.stats as stats

DIR = sys.argv[1] # Top level dir i.e., /blah/yyyymmdd/raw/

def defstuff():
	global RAW, CD, extensions, d, ds9x
	RAW = DIR +'/raw/'			# Top level directory with calibrated images
	CD = DIR + '/CALIBRATED/'	# Dir with raw files
	extensions=(1,2,3)			# Convenient way to call fits extensions
	#~ ds9x=ds9_xpans()
	#~ d=ds9()

def det_confs(DIR): # Retrieve all possible detector configurations
	dc = [] # List of all detector configurations ('BINX', 'BINY', 'WINDOWED', 'NAXIS1', 'NAXIS2')
	for i in glob.glob(RAW+"*fits"):
		tmp = []
		for KW in ['BINX', 'BINY']:
			tmp.append(fits.getval(i, KW))
		for KW in ['NAXIS1', 'NAXIS2']:
			tmp.append(fits.getval(i, KW, ext=1))
		dc.append(tmp)
	print "Determining all detector configurations."
	dc = set(map(tuple, dc)) # Unique detector configurations
	return dc 

def data_mgmt(dc):
	""" Creates neccessary dirs. """
	for i in dc:
		BIN=CD+'/'+str(i[0])+'x'+str(i[1])+'/'
		WD=BIN+str(i[-2])+'x'+str(i[-1])+'/' # Bottom level dir with calibrated and master frames
		if not os.path.exists(WD):
			os.makedirs(WD) # os.makedirs creates all neccessary dirs

	filelist = glob.glob(RAW + '*_OBJ.fits')
	
	objlist=[]
	for i in filelist:
		objlist.append(fits.getval(i, 'OBJECT'))
	
	objlist=np.unique(objlist)
	for i in objlist:
		if not os.path.exists(CD+'/'+i):
			os.makedirs(CD+'/'+i)

def CCD_sections(BIN):
	""" TRIM, OS, PS """

	""" 1x1 bining """
	if BIN==(1,1):
		TRIM1=(50,2095) ## CCD1 TRIM
		TRIM=(55,2098) ## CCD2+3 TRIM
		VR=(3,3072) ## Virtual rows TRIM
		PS1=(2100, 2125) ## CCD1 prescan
		PS=(25,50) ## ...
		OS1=(0,46)
		OS=(2103,2150)
		#~ (2150x3100)

	#~ """ 2x2 bining """
	elif BIN==(2,2):
		TRIM1=(23,1045) ## CCD1 TRIM
		TRIM=(28, 1051) ## CCD2+3 TRIM
		VR=(2,1536) ## Virtual rows TRIM
		PS1=(1047, 1056) ## CCD1 prescan
		PS=(16,27) ## ...
		OS1=(0,22)
		OS=(1052,1074)
		#~ CCDw = 1500 # CCD width
			
	#~ """ 3x3 bining review weird and ... too conservative?"""
	elif BIN==(3,3):
		TRIM1=(52,651) ## CCD1 TRIM ==> Warning nonlinear behavior in 3x3!!!
		TRIM=(52, 651) ## CCD2+3 TRIM
		VR=(2,1024) ## Virtual rows TRIM
		PS1=(698, 702) ## CCD1 prescan
		PS=(11,17) ## ...
		OS1=(0,11)
		OS=(704,715)
		#~ CCDw = 1500 # CCD width
	
	#~ """ 4x4 """
	elif BIN==(4,4):
		TRIM1=(35,515) ## CCD1 TRIM ==> weird things hapening with zero level again...
		TRIM=(35, 515) ## CCD2+3 TRIM
		VR=(2,768) ## Virtual rows TRIM
		PS1=(522, 530) ## CCD1 prescan
		PS=(8,13) ## ...
		OS1=(0,7)
		OS=(528,536)
		#~ CCDw = 1500 # CCD width

	return TRIM, TRIM1, VR, PS, PS1, OS, OS1

def check_exist(WD, ftype,i):
	timeout = 5
	if os.path.exists(WD+ftype):
		print "Master frame for", ftype,  i[0], 'x', i[1], 'binning and', i[-2], 'x', i[-1], 'windowing already exists. \nWould you like to create a new one? [y/n]'
		rlist, _, _ = select([sys.stdin], [], [], timeout)
		if rlist:
			B = raw_input()
		else:
			B='y'
	else:
		B='y'
	return B

def plotnp(im,fr):
	#~ d.set('frame delete all')
	d.set('tile yes')
	d.set('zscale')
	d.set('frame frameno ' +`fr`)
	d.set_np2arr(im)
	d.set('zoom to fit')
	d.set('cmap invert yes')

def slidingavg(mscrow, window_size):
    window = numpy.ones(int(window_size))/float(window_size)
    return np.convolve(mscrow, window, 'same')

def median_row(OSC, PSC, TR, im):
	mscrow = np.median(np.hstack((im[:,PSC[0]:PSC[1]], im[:,OSC[0]:OSC[1]])), axis=1)# median scan rows
	#~ sigmarow = np.std(np.hstack((im[:,PSC[0]:PSC[1]], im[:,OSC[0]:OSC[1]])), axis=1)# median scan rows
	mscrow = slidingavg(mscrow, 15)
	sigmarow = np.std(mscrow)
	return mscrow, sigmarow

def make_master_bias(dc):

	print "Making master BIAS"

	for i in dc:
		TRIM, TRIM1, VR, PS, PS1, OS, OS1 = CCD_sections((i[0], i[1]))
		filelist = []
		for f in glob.glob(RAW+"*BIAS*fits"):
			ccd_conf = []
			for KW in ['BINX', 'BINY']:
				ccd_conf.append(fits.getval(f, KW))
			for KW in ['NAXIS1', 'NAXIS2']:
				ccd_conf.append(fits.getval(f, KW, ext=1))
				if tuple(ccd_conf)==i:
					filelist.append(f)
		lfl = len(filelist)
		if lfl > 0:
			BIN=CD+'/'+str(i[0])+'x'+str(i[1])+'/'
			WD=BIN+str(i[-2])+'x'+str(i[-1])+'/' # Bottom level dir with calibrated and master frames
			B=check_exist(WD, 'MB.fits', i)
			if B=='n':
				pass
			else:
				hdul = fits.HDUList()
				hdul.append(fits.ImageHDU())
				x = np.array(range(0,i[-1]))
				for EXT in (extensions):
					print "##################################################"
					print "Stacking "+`lfl`+' '+`i[-2]`+'x'+`i[-1]`+' channel '+`EXT`+' bias frames!'
					if EXT==1:
						OSC=OS1
						PSC=PS1
					else:
						OSC=OS
						PSC=PS
					stack_arr = np.zeros( (len(filelist), i[-1], i[-2]) )
					for n, fn in enumerate(filelist):
						print "Files left:",`lfl-n`+'/'+`lfl`
						im = fits.getdata(fn, ext=EXT)
						row, col = [ range(0,d) for d in np.shape(im)]
						mscrow = np.median(np.hstack((im[:,PSC[0]:PSC[1]], im[:,OSC[0]:OSC[1]])), axis=1)# median scan rows
						crow = np.polyfit(row,mscrow, 11)
						frow = np.poly1d(crow) # fit rows
						stack_arr[n] = (im.T - frow(row)).T
						#~ stack_arr[n] = B
					median_stack = np.median(stack_arr, axis=0)
					hdul.append(fits.ImageHDU(median_stack))
				hdul[0].header.set("CALIBR", "T")
				hdul[0].header.set("INSTRUME", "MAIA")
				hdul[0].header.set("CALMODE", "MASTER BIAS")
				hdul[0].header.set("BINX", i[0])
				hdul[0].header.set("BINY", i[1])
				hdul.writeto(WD+"MB.fits", clobber=True)
				print "############################################################"
				#~ for EXT in (extensions):
					#~ print "Median BIAS EXT("+repr(EXT)+"):"+repr(np.round(np.median(hdul[EXT].data),2))
					#~ print "STD(EXT"+repr(EXT)+"):"+repr( np.round(np.std(hdul[EXT].data),2) )
	print "Completed master bias!"

def make_master_flats(dc):
	""" Make master FLAT for each MAIA CCD."""

	## Make EXTcheck: is there always the same number of extensions in each file
	print "Making master flats"
	
	## Choose extensions you are using
	
	for flat_type in ['FFS']: #  Currently FFD is unsupported. If you have FFDs, add them to the list but you must have ONLY FFDs or ONLY FFSs in the dir. Otherwise the first element in the list will get overwritten!
		#~ print "\n", flat_type, "\n"
		for i in dc:
			TRIM, TRIM1, VR, PS, PS1, OS, OS1 = CCD_sections((i[0], i[1]))
			filelist = []
			for f in glob.glob(RAW+'*'+flat_type+'*fits'):
				ccd_conf = []
				header0 = fits.getheader(f)
				header1 = fits.getheader(f, ext=1)
				if header0['OBSMODE']==flat_type:
					for KW in ['BINX', 'BINY']:
						ccd_conf.append(header0[KW])
					for KW in ['NAXIS1', 'NAXIS2']:
						ccd_conf.append(header1[KW])
						if tuple(ccd_conf)==i:
							filelist.append(f)
			lfl = len(filelist)
			if lfl > 0:
				BIN=CD+'/'+str(i[0])+'x'+str(i[1])+'/'
				WD=BIN+str(i[-2])+'x'+str(i[-1])+'/' # Bottom level dir with calibrated and master frames
				B=check_exist(WD, 'MF.fits', i)
				if B=='n':
					pass
				else:
					hdul = fits.HDUList()
					hdul.append(fits.ImageHDU())
					#~ MB = fits.open(WD+'MB.fits')
					x = np.array(range(0,i[-1]))
					for EXT in (extensions):
						print "##################################################"
						print "Stacking "+`lfl`+' '+`i[-2]`+'x'+`i[-1]`+' channel '+`EXT`+' flat frames!'
						if EXT==1:
							PSC=PS1
							OSC=OS1
							TR=TRIM1
						else:
							PSC=PS1
							OSC=OS
							TR=TRIM
						sc = -1 # counts how many flats have mean>limit
						for n, fn in enumerate(filelist):
							print "Files left:",`lfl-n`+'/'+`lfl`
							im = fits.getdata(fn, ext=EXT)
							meanval = np.mean(im[VR[0]:VR[1], TR[0]:TR[1]])
							#~ maxval = np.max(im[VR[0]:VR[1], TR[0]:TR[1]])
							maxval = stats.scoreatpercentile(im[VR[0]:VR[1], TR[0]:TR[1]], 90)
							exptime = fits.getheader(fn)['EXPTIME']
							#~ if meanval > 15000. and meanval < 40000. and maxval < 50000. and exptime>5.:
							if meanval > 16000. and meanval < 40000. and exptime>=5.:
								sc+=1
								#~ im[im<1]=1
								mscrow, sigmarow = median_row(OSC, PSC, TR, im)
								sh = np.shape(im)
								for y in range(0, sh[0]):
									im[y] = im[y]-mscrow[y]
								F=im
								norm = np.median(F[VR[0]:VR[1], TR[0]:TR[1]])
								F = F/norm #+np.min(F)+0.0001
								if sc==0:
									stack_arr = F
								else:
									stack_arr = np.dstack((stack_arr, F))
							else:
								print "****************************************"
								print "Rejected", fn, "AVG =", meanval, "EXPTIME =", exptime
								print "****************************************"
						print 'Will stack a total of', np.shape(stack_arr)[2], 'flats'
						MF = np.median(stack_arr, axis=2)
						hdul.append(fits.ImageHDU(MF))
						hdul[EXT].header.set("NAXIS1", np.shape(MF)[1])
						hdul[EXT].header.set("NAXIS2", np.shape(MF)[0])
					hdul[0].header.set("CALIBR", "T")
					hdul[0].header.set("INSTRUME", "MAIA")
					hdul[0].header.set("BINX", i[0])
					hdul[0].header.set("BINY", i[1])
					hdul[0].header.set("CALMODE", "MASTER FLAT")
					hdul.writeto(WD+"MF.fits", clobber=True)
					print "############################################################"
	print "Completed master flats"

def make_fake_masters(hdul, Mtype, ccd_shape, TD):
	""" In case master bias and/or master flat for a given configuration 
	are missing, try to reconstruct them from full frame masters. """
	
	if Mtype == "MB.fits":
		fn = Mtype
		print "NO MASTER BIAS AVAILABLE -- CROPPING A REPLACEMENT from FF BIAS"
	else:
		fn = Mtype
		print "NO MASTER FLAT AVAILABLE -- CROPPING A REPLACEMENT from FF FLAT"

	TRIM, TRIM1, VR, PS, PS1, OS, OS1 = CCD_sections((ccd_shape[1], ccd_shape[1]))
	BIN = `ccd_shape[1]` + 'x' + `ccd_shape[2]`
	FFDIR = np.sort(glob.glob(CD + '/' + BIN + '/*'))[-1]	 # Dir with the FF for given bining
	
	M = fits.HDUList()
	M.append(fits.ImageHDU())
	FFM = fits.open(FFDIR+'/MF.fits')
	
	for EXT in extensions:
		YSTART = int((hdul[0].header['YSTART']).strip('[,]')) #  Window y start pixel, virtual rows included
		YEND = YSTART + int(hdul[0].header['WINY']) # End of window y direction, virtual rows included
		data = FFM[EXT].data[YSTART:YEND]
		M.append(fits.ImageHDU(data))
		M[EXT].header.set("NAXIS1", np.shape(M[EXT].data)[1])
		M[EXT].header.set("NAXIS2", np.shape(M[EXT].data)[0])
		M[0].header.set("RECONSTR", 'Y')
		M.writeto(TD+'/'+fn, clobber=True)

	return M
	
def make_calibrated(dc):
	""" Calibrate science frames. """
	## Make EXTcheck: is there always the same number of extensions in each file
	
	filelist = glob.glob(RAW + '*_OBJ.fits')
	lfl = len(filelist)
	print "################################################################################"
	print "A total of", lfl, "science frames will be processed!"
	for n, fn in enumerate(filelist):
		#~ hdul = fits.HDUList()
		#~ hdul.append(fits.ImageHDU())
		print "################################################################################"
		print 'SCI frame' + fn.split('/')[-1] + '(' + `n` + '/' + `lfl` + ')'
		RS = fits.open(fn) # Raw SCI frame
		RS.info()
		c = []
		for KW in ['OBJECT', 'BINX', 'BINY']:
			c.append(RS[0].header[KW])
		for KW in ['NAXIS1', 'NAXIS2']:
			c.append(RS[1].header[KW])

		TRIM, TRIM1, VR, PS, PS1, OS, OS1 = CCD_sections((c[1], c[2]))
		TD = CD + `c[1]` + 'x' + `c[2]` + '/' + `c[3]`+ 'x' + `c[4]` + '/' # Target directory
		
		#~ if os.path.exists(TD+'MB.fits'):
			#~ MB = fits.open(TD+'MB.fits')
		#~ else:
			#~ MB = make_fake_masters(RS, 'MB.fits', c, TD)
		
		if os.path.exists(TD+'MF.fits'):
			MF = fits.open(TD+'MF.fits')
		else:
			MF = make_fake_masters(RS, 'MF.fits', c, TD)

		print "Files left:",`lfl-n`+'/'+`lfl`
		for EXT in extensions:
			if EXT==1:
				TR=TRIM1
				OSC=OS1
				PSC=PS1
			else:
				TR=TRIM
				OSC=OS
				PSC=PS1
			im = RS[EXT].data
			mscrow, sigmarow = median_row(OSC, PSC, TR, im)
			sh = np.shape(im)
			for y in range(0, sh[0]):
				im[y] = im[y]-mscrow[y]
			S = im/MF[EXT].data
			# Testirati radi li stvarno ono sto mislim da radi!!!
			if fits.getheader(fn)['WINDOWED'] == 'TRUE': # Trim final images, if windowing is enforced special treatment is required
				ystart = int(fits.getheader(fn)['YSTART'].replace('[','').replace(']','')) # YSTART
				winy = int(fits.getheader(fn)['WINY']) # Window width y-dim
				ystop = ystart+winy
				if ystart>=35:
					ystart=0
				elif ystart<35:
					ystart = 35-ystart
				if ystop<=3070:
					ystop=np.shape(S)[0]
				elif ystop>3070:
					ystop = 3070-ystop
				print ystart, ystop
				S = S[ystart:ystop,55:2095]
			else:
				#~ S = S[35:3070,55:2095]
				S = S[VR[0]:VR[1], TR[0]:TR[1]]
			#~ hdul.append(fits.ImageHDU(S))
			RS[EXT].data=S
			RS[EXT].header.set("NAXIS1", np.shape(RS[EXT].data)[1])
			RS[EXT].header.set("NAXIS2", np.shape(RS[EXT].data)[0])
			#~ plotnp(S, 2*EXT)
		#~ print RS.info()
		#~ raw_input()
		RS[0].header.set("CALIBR", "T")
		fn=CD+'/'+c[0]+'/c'+fn.split('/')[-1]
		RS.writeto(fn, clobber=True)

defstuff()
dc = det_confs(DIR)
data_mgmt(dc)
#~ make_master_bias(dc)
make_master_flats(dc)
make_calibrated(dc)
