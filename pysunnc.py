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
	newxy2=np.dot(xy2, result['fit_matrix'])+result['offset']
	newxy2[:, 0]=newxy2[:, 0]-center[0]
	newxy2[:, 1]=newxy2[:, 1]-center[1]
	newx2s=newxy2[:, 0]
	newy2s=newxy2[:, 1]#}}}#
	return (newx2s, newy2s)

def geomapscript(wd, x1s, y1s, x2s, y2s, nx, ny, \
	scriptfile='a.cl', coofile='coo', dbfile='db', \#{{{#
	fitgeometry="general", function="polynomial", maxiter=10, reject=3.0, \
	xxorder=2, xyorder=2, xxterms="half", yxorder=2, yyorder=2, yxterms="half"):
	with open(wd+coofile,'w+') as f:
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
	scriptfile='b.cl', xyfile='xy', newxyfile='newxy', coofile='coo', dbfile='db'):#{{{#
	with open(wd+xyfile,'w+') as f:
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

