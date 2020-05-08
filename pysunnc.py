def world2pix(ra, dec, ra0, dec0):
	from astropy import wcs
	w=wcs.WCS(naxis=2)
	w.wcs.crpix=[0.0, 0.0]
	w.wcs.cdelt=[1.0, 1.0]
	w.wcs.crval=[ra0, dec0]
	w.wcs.ctype=["RA---TAN", "DEC--TAN"]
	world=np.transpose(np.vstack([tab['ra'], tab['dec']]))
	pix=w.all_world2pix(world, 1)
	return (pix[:, 0], pix[:, 1])

def pix2world(xx, yy, ra0, dec0):
	from astropy import wcs
	w=wcs.WCS(naxis=2)
	w.wcs.crpix=[0.0, 0.0]
	w.wcs.cdelt=[1.0, 1.0]
	w.wcs.crval=[ra0, dec0]
	w.wcs.ctype=["RA---TAN", "DEC--TAN"]
	pix=np.transpose(np.vstack([xx, yy]))
	world=w.all_pix2world(pix, 1)
	return (world[:, 0], world[:, 1])

def savefig(fig, filename, openfig=True):
	import os
	import matplotlib.pyplot as plt
	fig.tight_layout()
	fig.savefig(filename)
	plt.close(fig)
	if openfig: os.system('open '+filename)

def pload(filename):
	import pickle
	obj=pickle.load(open(filename, 'rb'))
	return obj

def pdump(obj, filename):
	import pickle
	pickle.dump(obj, open(filename, 'wb'))

def ang2pc(arcsec, distance):
	import numpy as np
	return np.deg2rad(arcsec/3600.0)*distance

def d2m(d):
	import numpy as np
	return 5.0*np.log10(d)-5.0

def m2d(m):
	return 10**(m/5.0+1.0)

def drizit(flcinput, drcoutput, skysub=True, grow=3):
	from drizzlepac import astrodrizzle#{{{#
	astrodrizzle.AstroDrizzle(flcinput, output=drcoutput, \
		clean=True, \
		runfile='', \
		coeffs=True, \
		context=True, \
		group='', \
		build=True, \
		crbit=4096, \
		static=True, \
		static_sig=4.0, \
		skysub=skysub, \
		skywidth=0.1, \
		skystat='median', \
		skylower=-100.0, \
		skyupper=None, \
		skyclip=5, \
		skylsigma=4.0, \
		skyusigma=4.0, \
		driz_separate=True, \
		driz_sep_kernel='turbo', \
		driz_sep_pixfrac=1.0, \
		driz_sep_fillval=None, \
		driz_sep_bits=96, \
		driz_sep_rot=None, \
		driz_sep_scale=None, \
		median=True, \
		median_newmasks=True, \
		combine_type='minmed', \
		combine_nsigma='4 3', \
		combine_nlow=0, \
		combine_nhigh=1, \
		combine_lthresh=None, \
		combine_hthresh=None, \
		combine_grow=1, \
		blot=True, \
		blot_interp='poly5', \
		blot_sinscl=1.0, \
		driz_cr=True, \
		driz_cr_corr=False, \
		driz_cr_snr='3.5 3.0', \
		driz_cr_grow=grow, \
		driz_cr_scale='1.2 0.7', \
		driz_combine=True, \
		final_kernel='square', \
		final_pixfrac=1.0, \
		final_fillval=None, \
		final_bits=96, \
		final_rot=None, \
		final_scale=None)#}}}#

def gotods9(regfilename, xpos, ypos, color='green'):
	line1='# Region file format: DS9 version 4.1'#{{{#
	line2='global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" '
	line2=line2+'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
	line3='image'
	print(line1, file=open(regfilename, 'w'))
	print(line2, file=open(regfilename, 'a'))
	print(line3, file=open(regfilename, 'a'))
	for ii in range(len(xpos)):
		line='point('+"%8.3f"%(xpos[ii]+0.5)+','+"%8.3f"%(ypos[ii]+0.5)+') # point=circle'
		print(line, file=open(regfilename, 'a'))#}}}#

def matchstars(file1, file2):
	import numpy as np#{{{#
	with open(file1) as f: content1=f.readlines()
	idx1s=[]
	idstr1s=[]
	for ii in range(len(content1)):
		if 'text={' in content1[ii]:
			idstr=content1[ii].split('text={')[1].split('}')[0]
			idx1s.append(ii)
			idstr1s.append(idstr)
	idxsort1s=np.argsort(idstr1s)
	with open(file2) as f: content2=f.readlines()
	idx2s=[]
	idstr2s=[]
	for ii in range(len(content2)):
		if 'text={' in content2[ii]:
			idstr=content2[ii].split('text={')[1].split('}')[0]
			idx2s.append(ii)
			idstr2s.append(idstr)
	idxsort2s=np.argsort(idstr2s)
	if len(idx1s)!=len(idx2s):
		print('error')
		return 0
	else:
		nmatch=len(idx1s)
		x1s=np.zeros(nmatch, dtype=float)
		y1s=np.zeros(nmatch, dtype=float)
		x2s=np.zeros(nmatch, dtype=float)
		y2s=np.zeros(nmatch, dtype=float)
		for ii in range(nmatch):
			str1=content1[idx1s[idxsort1s[ii]]]
			str2=content2[idx2s[idxsort2s[ii]]]
			str1=str1.split('point(')[1].split(') #')[0]
			str2=str2.split('point(')[1].split(') #')[0]
			x1=float(str1.split(',')[0])
			y1=float(str1.split(',')[1])
			x2=float(str2.split(',')[0])
			y2=float(str2.split(',')[1])
			x1s[ii]=x1
			y1s[ii]=y1
			x2s[ii]=x2
			y2s[ii]=y2#}}}#
		return (nmatch, x1s, y1s, x2s, y2s)

def linearfitmap(nmatch, x1s, y1s, x2s, y2s, \
	xyorig=None, uvorig=None,  mode='general', \
	nclip=3, sigma=3.0, minobj=3, center=[0.0, 0.0], verbose=False):
	import numpy as np#{{{#
	from drizzlepac import linearfit
	xy1=np.vstack([x1s, y1s]).T
	xy2=np.vstack([x2s, y2s]).T
	result=linearfit.iter_fit_all(np.vstack([x1s, y1s]).T, np.vstack([x2s, y2s]).T, \
		np.arange(nmatch), np.arange(nmatch), \
		xyorig=xyorig, uvorig=uvorig,  mode=mode, \
		nclip=nclip, sigma=sigma, minobj=minobj, center=center, verbose=verbose)#}}}#
	return (center, result)

def linearfittran(center, result, x2s, y2s):
	import numpy as np#{{{#
	xy2=np.vstack([x2s, y2s]).T
	xy2[:, 0]=xy2[:, 0]+center[0]
	xy2[:, 1]=xy2[:, 1]+center[1]
	xy2to1=np.dot(xy2, result['fit_matrix'])+result['offset']
	xy2to1[:, 0]=xy2to1[:, 0]-center[0]
	xy2to1[:, 1]=xy2to1[:, 1]-center[1]
	x2to1s=xy2to1[:, 0]
	y2to1s=xy2to1[:, 1]#}}}#
	return (x2to1s, y2to1s)

def geomapscript(wd, x1s, y1s, x2s, y2s, nx, ny, \
	scriptfile='a.cl', coofile='coo', dbfile='db', \
	fitgeometry="general", function="polynomial", maxiter=10, reject=3.0, \
	xxorder=2, xyorder=2, xxterms="half", yxorder=2, yyorder=2, yxterms="half"):
	with open(wd+coofile,'w+') as f:#{{{#
		for ii in range(len(x1s)):
			f.write("%8.2f %8.2f %8.2f %8.2f\n"%(x1s[ii], y1s[ii], x2s[ii], y2s[ii]))
	with open(wd+scriptfile, 'w+') as f:
		f.write('geomap ("'+coofile+'", "'+dbfile+'", 1, '+str(nx)+', 1, '+str(ny)+',\n')
		st='fitgeometry="'+fitgeometry+'", function="'+function+'", '
		st=st+'maxiter='+str(maxiter)+', reject='+str(reject)+',\n'
		f.write(st)
		st='xxorder='+str(xxorder)+', xyorder='+str(xyorder)+', xxterms="'+xxterms+'", '
		st=st+'yxorder='+str(yxorder)+', yyorder='+str(yyorder)+', yxterms="'+yxterms+'",\n'
		f.write(st)
		st='transforms="", results="", calctype="real", '
		st=st+'verbose=yes, interactive=no, graphics="stdgraph", cursor="")'#}}}#

def geoxytranscript(wd, u1s, v1s, nx, ny, \
	scriptfile='b.cl', xyfile='xy', newxyfile='newxy', coofile='coo', dbfile='db'):
	with open(wd+xyfile,'w+') as f:#{{{#
		if np.size(u1s)==1:
			f.write("%8.2f %8.2f\n"%(u1s, v1s))
		else:
			for ii in range(np.size(u1s)):
				f.write("%8.2f %8.2f\n"%(u1s[ii], v1s[ii]))
	with open(wd+scriptfile, 'w+') as f:
		f.write('geoxytran ("'+xyfile+'", "'+newxyfile+'", "'+dbfile+'", "'+coofile+'",\n')
		f.write('geometry="geometric", direction="forward", xref=INDEF, yref=INDEF, xmag=INDEF,\n')
		f.write('ymag=INDEF, xrotation=INDEF, yrotation=INDEF, xout=INDEF, yout=INDEF,\n')
		f.write('xshift=INDEF, yshift=INDEF, xcolumn=1, ycolumn=2,\n')
		f.write('calctype="real", xformat="", yformat="", min_sigdigit=7)')#}}}#

def dpread(zzdp, out):
	import os#{{{#
	import numpy as np
	from astropy.table import Table

	#get filters
	with open(zzdp+'.columns') as f: content=f.readlines()
	filters=[]
	for jj in range(len(content)):
		if 'Total counts,' in content[jj]:
			filters.append(content[jj].split(', ')[-1].replace('\n', ''))
	#get rows and names
	rh=np.array([3, 4, 6, 7, 8, 10, 11])-1
	r0=np.array([16, 18, 20, 21, 22, 23, 24])-1
	usecols=[rh]
	dtype=[('xpos', float), ('ypos', float), \
		('snr', float), ('shp', float), \
		('rnd', float), ('crd', float), \
		('obj', int)]
	for jj in range(len(filters)):
		usecols.append(r0+13*jj)
		dtype=dtype+[('mag'+filters[jj], float), ('err'+filters[jj], float), \
			('snr'+filters[jj], float), ('shp'+filters[jj], float), \
			('rnd'+filters[jj], float), ('crd'+filters[jj], float), \
			('qfg'+filters[jj], int)]
	usecols=list(np.hstack(usecols))
	#read
	srcTab=Table(np.loadtxt(zzdp, usecols=usecols, dtype=dtype))
	#write
	srcTab.write(out+'.dat', format='ascii', overwrite=True)
	gotods9(out+'.reg', srcTab['xpos'], srcTab['ypos'])#}}}#

def readsynspec(filename):
	import numpy as np#{{{#
	from astropy.table import Table
	usecols=np.arange(9)
	dtype=[('energy', float), ('incident', float), ('trans', float), \
		('diffout', float), ('netrans', float), ('reflc', float), \
		('total', float), ('reflin', float), ('outlin', float)]
	tab=Table(np.loadtxt(filename, usecols=usecols, dtype=dtype))
	#photon energy is in units of Ry (1 Ry = 13.606 eV), get frequency
	tab['freq']=tab['energy']*13.605693122994/4.135667696e-15 # in Hz
	#the flux is given in nuLnu, convert to Lnu
	for col in tab.colnames[1:9]: tab[col+'_freq']=tab[col]/tab['freq'] # in erg/s/Hz
	#get wavelength
	tab['wave']=2.99792458e18/tab['freq'] # in Angstrom
	#convert flux into Llambda
	for col in tab.colnames[1:9]: tab[col+'_wave']=tab[col]/tab['wave'] # in erg/s/Angstrom
	return tab#}}}#

def readspecbpass(filename):
	import numpy as np#{{{#
	from astropy.table import Table
	spec=Table(np.loadtxt(filename))
	spec.rename_column('col0', 'wave')
	for ii in range(1, len(spec.colnames)): \
		spec.rename_column(spec.colnames[ii], str(round(6+0.1*(ii-1), 1)))
	return spec#}}}#

def addmag(mag1, mag2):
	import numpy as np#{{{#
	mag=mag1-2.5*np.log10(1+10**((mag1-mag2)/2.5))
	return mag#}}}#

def submag(mag, mag1):
	import numpy as np#{{{#
	mag2=mag-2.5*np.log10(1-10**((mag-mag1)/2.5))
	return mag2#}}}#

def reddenspec(sp, EBV, ext):
	import astropy.units as u#{{{#
	import numpy as np
	import pysynphot as S
	idx=np.logical_and(sp.wave>(1.0/10.0)*1e4, sp.wave<(1.0/0.3)*1e4)
	wave=sp.wave[idx]*u.AA
	flux=sp.flux[idx]
	ee=ext.extinguish(wave, Ebv=EBV)
	flux_ext=flux*ee
	sp_ext=S.ArraySpectrum(sp.wave[idx], flux_ext, waveunits='angstrom', fluxunits='flam')
	return sp_ext#}}}#

def readbpstarmodel(filename):
	import numpy as np#{{{#
	from astropy.table import Table
	usecols=np.array([1, 2, \
		3, 4, 5, 6, \
		7, 8, 9, \
		11, 12, 13, 14, 15, 16, \
		17, 18, 19, 20, 21, 22, 23, 24, 25, \
		26, 27, \
		28, 29, 30, 31, 32, 33, \
		34, 35, 36, \
		38, 39, \
#       40, 41, 42, 43, 44, 45, 46, \
		42, 43, 44, 45, \
		47, 48, 49, \
		50, 51, 52])-1
	dtype=[('step', float), ('age', float), \
		('LogR1', float), ('LogT1', float), ('LogL1', float), ('Mact1', float), \
		('massHeCore', float), ('massCOCore', float), ('massONeCore', float), \
		('Xsurf1', float), ('Ysurf1', float), ('Csurf1', float), ('Nsurf1', float), ('Osurf1', float), ('Nesurf1', float), \
		('massH1', float), ('massHe1', float), ('massC1', float), ('massN1', float), ('massO1', float), ('massNe1', float), ('massMg1', float), ('massSi1', float), ('massFe1', float), \
		('EnvBinEnergy', float), ('TotStarBinEnergy', float), \
		('massRemWeakSN', float), ('massEjWeakSN', float), ('massRemSN', float), ('massEjSN', float), ('massRemSuperSN', float), ('massEjSuperSN', float), \
		('J', float), ('period', float), ('LogSep', float), \
		('Mact2', float), ('Mtot', float), \
#       ('DM1wind', float), ('DM2wind', float), ('DM1acc', float), ('DM2acc', float), ('DM1roche', float), ('DM2roche', float), ('DJ', float), \
		('DM1acc', float), ('DM2acc', float), ('DM1roche', float), ('DM2roche', float), \
		('LogR2', float), ('LogT2', float), ('LogL2', float), \
		('rocheOverFlux2', float), ('modelimf', float), ('mixedimf', float)]
	evTab=Table(np.loadtxt(filename, usecols=usecols, dtype=dtype))
	return evTab#}}}#
