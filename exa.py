#!/usr/bin/env python

from numpy import mean
from numpy import std
from numpy import array
from numpy import sort
from numpy import arange
from numpy import argsort

import os,sys,string,re,shutil,getopt
from math import log10
from math import sqrt
from snoopy2 import *
import snoopy2
import math
import sys


help ="################################################################ \n\
help = Usage:   exa.py filename             \n\
################################################################"

def readdoc(file_in):
    f=file(file_in,'r')
    vecname,vecU,vecB,vecV,vecR=[],[],[],[],[]
    vecI,vecUerr,vecBerr,vecVerr,vecRerr,vecIerr=[],[],[],[],[],[]
    night=f.readline()
    JD=night.split()[1]
    telescope=night.split()[2]
    photometric=night.split()[3]
    s=f.readlines()
    for ii in s:
        if ii[0]!='#':
            c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11=ii.split() 
            vecname.append(c1)
            vecU.append(float(c2))
            vecUerr.append(float(c3))
            vecB.append(float(c4))
            vecBerr.append(float(c5))
            vecV.append(float(c6))
            vecVerr.append(float(c7))
            vecR.append(float(c8))
            vecRerr.append(float(c9))
            vecI.append(float(c10))
            vecIerr.append(float(c11))
    return  vecname,vecU,vecB,vecV,vecR,vecI,vecUerr,vecBerr,vecVerr,vecRerr,vecIerr,telescope,JD,photometric

def readdoc_i(file_in):
    f=file(file_in,'r')
    vecname,vecU,vecB,vecV=[],[],[],[]
    vecUerr,vecBerr,vecVerr=[],[],[]
    night=f.readline()
    JD=night.split()[1]
    telescope=night.split()[2]
    photometric=night.split()[3]
    s=f.readlines()
    for ii in s:
        if ii[0]!='#':
            c1,c2,c3,c4,c5,c6,c7=ii.split() 
            vecname.append(c1)
            vecU.append(float(c2))
            vecUerr.append(float(c3))
            vecB.append(float(c4))
            vecBerr.append(float(c5))
            vecV.append(float(c6))
            vecVerr.append(float(c7))
    return  vecname,vecU,vecB,vecV,vecUerr,vecBerr,vecVerr,telescope,JD,photometric

def pval0(_xx, p):
    _y=+p[0]+0*_xx
    return _y

##############################################################################

if len(sys.argv)==1:
    print help
    sys.exit()

p0=[0,0]


interactive=False
userPath=os.getcwd()

options, args = getopt.getopt(sys.argv[2:],'i',['interactive'])
for opt,arg in options:
    if opt in ('-i', '--interactive'): interactive = True


if sys.argv[1][0]=='@':
   ff = open(sys.argv[1][1:])
   files = ff.readlines()
   imglist = []
   for ff in files: 
          if not ff=='\n' and ff[0]!='#':
                 ff=re.sub('\n','',ff)
                 ff=re.sub(' ','',ff)
                 imglist.append(ff)
elif ',' in sys.argv[1]: imglist = string.split(sys.argv[1],sep=',')
else:
    imglist = [sys.argv[1]]

system=raw_input('Which photometric system? UBVRI [0] (or press enter), JHK [1], ugriz [2]')
if not system: system=0
try: system=int(system)
except: system=0


stars={}
sequence_star=[]
telescop,JD,fff={},[],[]
photometric={}
U,B,V,R,I={},{},{},{},{}
u,g,r,i,z={},{},{},{},{}
J,H,K={},{},{}
for img in imglist:
  if img and img[0]!='#':
    print '#####'+str(img)
    try:
        if system==0:
          vecname,vecU,vecB,vecV,vecR,vecI,vecUerr,vecBerr,vecVerr,vecRerr,vecIerr,_telescope,_JD,_photometric=readdoc(img)
	elif system==1:
	    vecname,vecU,vecB,vecV,vecUerr,vecBerr,vecVerr,_telescope,_JD,_photometric=readdoc_i(img)
	elif system==2:
	    vecname,vecU,vecB,vecV,vecR,vecI,vecUerr,vecBerr,vecVerr,vecRerr,vecIerr,_telescope,_JD,_photometric=readdoc(img)
        JD.append(_JD)
	fff.append(_JD)
        print string.count(JD,_JD)
        if string.count(JD,_JD)>1:
            print 'WARNING: Two nights have same JD'
            print 'second night labeled with (a)'
            _JD=str(_JD)+'(a)'
            JD[-1]=_JD
	    fff[-1]=_JD
        stars[_JD]=vecname
        telescop[_JD]=_telescope
	photometric[_JD]=_photometric
	if system==0:
           U[_JD]=vecU
           B[_JD]=vecB
           V[_JD]=vecV
           R[_JD]=vecR
           I[_JD]=vecI
	elif system==1:
           J[_JD]=vecU
           H[_JD]=vecB
           K[_JD]=vecV
	elif system==2:
           u[_JD]=vecU
           g[_JD]=vecB
           r[_JD]=vecV
           i[_JD]=vecR
           z[_JD]=vecI	   	
        #telescop.append(_telescope)
        for ii in vecname:
            if string.count(sequence_star,ii)==0:
                sequence_star.append(ii)
    except:
        print '####'
        print '#### WARNING: file format incorrect !'
	print '####'
  else:
    	print '####'
  	print '#### WARNING: Empty space in the list !!'
	print '####'

fff.sort()
JDsort=[]
for ii in fff:
    ddd=string.count(str(JD)[:string.find(str(JD),ii)],',')
    JDsort.append(ddd)

bands={}
if system==0:
  bands['U']=U
  bands['B']=B
  bands['V']=V
  bands['R']=R
  bands['I']=I
elif system==2:
  bands['u']=u
  bands['g']=g
  bands['r']=r
  bands['i']=i
  bands['z']=z
elif system==1:  
  bands['J']=J
  bands['H']=H
  bands['K']=K

bandelista=['UBVRI','JHK','ugriz']
bandaref=['V','J','r']
refbanda=bandaref[system]
listabande=bandelista[system]
mag_stella={}
mag_stella_err={}
numero=[]
mediav=[]
for j in sequence_star:
     numero.append(string.count(str(stars),"'"+j+"'"))
     if system==0:
          mag_stella[j]=[0,0,0,0,0]
          mag_stella_err[j]=[0,0,0,0,0]
     elif system==1:
           mag_stella[j]=[0,0,0]
           mag_stella_err[j]=[0,0,0]    
     elif system==2:
          mag_stella[j]=[0,0,0,0,0]
          mag_stella_err[j]=[0,0,0,0,0]     
     mag=[]
     for ii in bands[refbanda]:
         if string.count(str(stars[ii]),"'"+j+"'"):
             if abs(bands[refbanda][ii][stars[ii].index(j)])<=30:
                  mag.append(bands[refbanda][ii][stars[ii].index(j)])
	          print bands[refbanda][ii][stars[ii].index(j)],ii
     print mag
     if mag:
        mediav.append(mean(array(mag)))
     else:
        mediav.append(9999)
        	

order=argsort(mediav)
stellaorder=[]
for j in range(len(sequence_star)):
	stellaorder.append(sequence_star[order[j]])

bandanum=0
stellanum=0
_stella=stellaorder[stellanum]
_banda=listabande[bandanum]
answ0='y'
while answ0=='y':
    _stella=stellaorder[stellanum]
    _banda=listabande[bandanum]
    answ='n'
    while answ=='n':
        print 'SEQUENCE STARS'
        print sequence_star
	print 'Number of nights the star was measured'
	print numero
        stella=raw_input('Which sequence star ['+_stella+'] ? ')
	if not stella: stella=_stella
	_stella=stella
        try:
            inde=sequence_star.index(stella)
            jj=sequence_star[inde]
            answ='y'
        except:
            print 'ERROR: star "'+stella+'" not in the list ! '
            print 'try again '
    
    banda=raw_input('Which BAND ['+_banda+'] ? ')
    if not banda: banda=_banda
    _banda=banda
    print '### star= '+str(jj)
    print '### BAND= '+str(banda)
    ll=0
    xx=[]
    yy=[]
    zz=[]
    _label=[]
    position=[]
    for kk in range(len(bands[banda])):
        ii=JD[JDsort[kk]]
        if string.count(stars[ii],jj)==1 and abs(bands[banda][ii][stars[ii].index(jj)])<=30:
            ll=ll+1
            print ll,ii,bands[banda][ii][stars[ii].index(jj)]
            xx.append(ll)
            yy.append(bands[banda][ii][stars[ii].index(jj)])
	    if photometric[ii]=='ph':
		    zz.append(1)
	    elif photometric[ii]=='fit':
		    zz.append(2)
	    else:
	            zz.append(0)
            position.append([ll,ll])
	    titolo2='%6.6s   %15.15s' % (str(telescop[ii]),str(ii))
	    _label.append(titolo2)
####################################################################
    if xx:
      if len(xx)>=2:
        media=mean(array(yy))
        error=std(array(yy))
        phmin=min(xx)
        phmax=max(xx)
        mm=min(yy)-(max(yy)-min(yy))/10.
        mma=max(yy)+(max(yy)-min(yy))/10.
        phr=phmax-phmin
        if len(str(phr))==3:
           phr=int((phmax-phmin)/100)*15
        else:
          phr=int((phmax-phmin)/10)   
        if phr==0: 
		phr=1
        mmr=abs(int((mma-mm)/10))
        if mmr==0:
           mmr=abs((mma-mm)/10)	
        p=[media,error]
        xfit=arange(phmin,phmax,.1)
        yfit=pval0(xfit, p)
        pnew=p

        pnew=src.grafico_exa(xx,yy,position,zz,_label,p,banda,xfit,yfit,phmin,phmax,mm,mma,interactive,pnew)
        print '####'+str(pnew[0])
        mag_stella[jj][listabande.index(banda)]=pnew[0]
        mag_stella_err[jj][listabande.index(banda)]=pnew[1]
      else:
          print 'Warning: only one value !!!!'
          mag_stella[jj][listabande.index(banda)]=yy[0]
          mag_stella_err[jj][listabande.index(banda)]=0
    else:
        print 'Warning no measurements for star '+str(jj)+' in band '+str(banda)+' !'

    bandanum=bandanum+1
    if bandanum==len(listabande):
           stellanum=stellanum+1
  	   bandanum=0
    if stellanum==len(sequence_star):
          bandanum=0
  	  stellanum=0
	  print '########################################################## '
	  print '#### '
	  print "####            All the stars should have been measured  !!!"
	  print '#### '
	  print '########################################################## '
    answ0=raw_input('again ? [y/n]  [y]')
    if not answ0:
          answ0='y'

print ''
print '#### SEQUENCE STARS: MAGNITUDE'
f=open('sequence_star_mag.list','w')
for j in stellaorder:
   if system==0 or system==2:
      aaa = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t' % (str(j),str(mag_stella[j][0]),str(mag_stella_err[j][0]),str(mag_stella[j][1]),str(mag_stella_err[j][1]),str(mag_stella[j][2]))
      bbb = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s' % (str(mag_stella_err[j][2]),str(mag_stella[j][3]),str(mag_stella_err[j][3]),str(mag_stella[j][4]),str(mag_stella_err[j][4]))
   elif system==1: 
      aaa = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t' % (str(j),str(mag_stella[j][0]),str(mag_stella_err[j][0]),str(mag_stella[j][1]),str(mag_stella_err[j][1]),str(mag_stella[j][2]))
      bbb = '%6.6s\t' % (str(mag_stella_err[j][2]))   
   print aaa+bbb
   f.write(aaa+bbb+'\n')
f.close()

print ''
print '#### SEQUENCE STARS: LIST FORM'

f=open('sequence_star.list','w')
for j in stellaorder:
   if system==0:
          aaa = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s' % (str(j),str(mag_stella[j][2]),str(mag_stella[j][1]-mag_stella[j][2]),str(mag_stella[j][0]-mag_stella[j][1]),str(mag_stella[j][2]-mag_stella[j][3]),str(mag_stella[j][3]-mag_stella[j][4]))
   elif system==1: 
          aaa = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t' % (str(j),str(mag_stella[j][0]),str(mag_stella[j][0]-mag_stella[j][1]),str(mag_stella[j][1]-mag_stella[j][2]))   
   elif system==2:
          aaa = '%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s' % (str(j),str(mag_stella[j][2]),str(mag_stella[j][1]-mag_stella[j][2]),str(mag_stella[j][0]-mag_stella[j][1]),str(mag_stella[j][2]-mag_stella[j][3]),str(mag_stella[j][3]-mag_stella[j][4]))
   print aaa
   f.write(aaa+'\n')
f.close()

