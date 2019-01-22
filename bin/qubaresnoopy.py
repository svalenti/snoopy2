#!/usr/bin/env python

from numpy import mean
import os,sys,string,re,shutil
from math import log10
from math import sqrt
from snoopy2 import *
import getopt
import snoopy2
from astropy.io import fits as pyfits


help ="################################################################ \n\
help = Usage:   ph_new.py filename             \n\
                input format: The output file of svstandard.py (eg lp_ru149_20080212)    \n\
                [-r  list] list format: the input file of svstandard.py (eg lista_lp_20080212)   \n\
                -z          no error in the zeropoint             \n\
################################################################"

#######################################################################
if len(sys.argv)==1:
    print help
    sys.exit()

parent_dir = os.getcwd()+'/'

imglist=src.readlist(sys.argv[1])

def erroremag(z0,z1,m0,m1,c0,c1,position): #  z=zeropoint,m=magnitude,colorterm  
    if position==0:   #    error for first band in the color: (e.g.  V in VR) 
        dm0=1+(c0/(1-(c0-c1)))
        dz0=1+(c0/(1-(c0-c1)))
        dm1=(-1)*(c0/(1-(c0-c1)))
        dz1=(-1)*(c0/(1-(c0-c1)))
        dc0=(z0+m0-z1-m1)*(1+c1)*(1/(1-(c0-c1))**2)
        dc1=(-1)*(z0+m0-z1-m1)*(c0)*(1/(1-(c0-c1))**2)
    elif position==1:   #    error for second band in the color: (e.g.  R in VR) 
        dm0=1-(c1/(1-(c0-c1)))
        dz0=1-(c1/(1-(c0-c1)))
        dm1=(-1)*(c1/(1-(c0-c1)))
        dz1=(-1)*(c1/(1-(c0-c1)))
        dc0=(z0+m0-z1-m1)*(1-c0)*(1/(1-(c0-c1))**2)
        dc1=(z0+m0-z1-m1)*(c1)*(1/(1-(c0-c1))**2) 
    else:
        # added to make the pipeline working, but error not good
        dm0=1
        dz0=0
        dm1=1
        dz1=0
        dc0=0
        dc1=0
    return dc0,dc1,dz0,dz1,dm0,dm1

####################################### check ##########################
_telescope=''
_system=''
check=0
_zp=False
_telescope=src.telescope(imglist[0])
print '### TELESCOPE= '+_telescope

knew=src.atmospheric_site2()
kk={}

options, args = getopt.getopt(sys.argv[2:],'s:,z',['system=','zeroerror'])
for opt,arg in options:
    if opt in ('-s', '--system'): _system = int(arg)
    if opt in ('-z', '--zeroerror'): _zp = True

bandsvec='UBVRIugrizJHK'
mag,magerr,magstr,zeropoint,zeroerr,colore,colorevalue,colorerr,magerr_new,mag_new={},{},{},{},{},{},{},{},{},{}
for i in bandsvec:
    mag[i],magerr[i],magstr[i],zeropoint[i],zeroerr[i],colorevalue[i],colorerr[i],mag_new[i],magerr_new[i]=0,0,0,0,0,0,0,0,0
    colore[i]='none'

filtri=[]
airmass={}
exptime={}
JD={}
listsystem=[]
for img in imglist:
  if img:
    if img[-5:]=='.fits': 
        imgs=img[:-5]
    else: 
        imgs=img
    if not _system:
        system=src.check_system(_telescope,img,Stdout=False)
    else:
        system=_system
    if system not in listsystem: listsystem.append(system)
    _header=src.read_parameter(_telescope,system)
    _filter=src.filter(img,_header,_telescope)
    filter=src.filtername(_telescope,_filter,system)
    if filter=='unknown':
        filter=_filter
    filtri.append(filter)
    print '### FILTER = '+str(filter)
    try:
       _airmass=float(src.airmass(img,_header,_telescope))
       airmass[filter]=_airmass
       print "### AIRMASS header check ......ok "
       print "### AIRMASS = "+str(_airmass)
    except:
        print 'ERROR: NO AIRMASS HEADER FOUND'
        sys.exit()

    try:
       _exptime=src.exptime(img,_header,_telescope)
       exptime[filter]=_exptime
       print "### exptime header check ......ok "
       print "### EXPTIME = "+str(_exptime)
    except:
       print 'ERROR: NO EXPTIME HEADER FOUND'
       sys.exit()

    try:
       _constant=pyfits.getheader(img)['QUBACONS'] 
       _zeropoint=float(string.split(pyfits.getheader(img)['QUBACONS'])[1]) 
       _colore=string.split(pyfits.getheader(img)['QUBACONS'])[2] 
       _colorevalue=float(string.split(pyfits.getheader(img)['QUBACONS'])[3]) 
       zeropoint[filter]=raw_input('zero point (with filter = '+filter+')?['+str(_zeropoint)+'] ')
       if not zeropoint[filter]: zeropoint[filter]=_zeropoint
       zeropoint[filter]=float(zeropoint[filter])
       colore[filter]=raw_input('color term [BV,VR,RI] (with filter '+filter+')?['+str(_colore)+']  ')
       if not colore[filter]: colore[filter]=_colore
       colorevalue[filter]=raw_input('color term value (with filter '+filter+')?['+str(_colorevalue)+']  ')
       if not colorevalue[filter]: colorevalue[filter]=_colorevalue
       colorevalue[filter]=float(colorevalue[filter])
       try:
           _zeroerr=float(string.split(pyfits.getheader(img)['QUBACONS'])[5]) 
           zeroerr[filter]=float(_zeroerr)
           _colorerr=float(string.split(pyfits.getheader(img)['QUBACONS'])[6]) 
           colorerr[filter]=float(_colorerr)
       except:
           zeroerr[filter]=float(0.0)
           colorerr[filter]=float(0.0)

       print "### CONSTANTS header check ......ok "
       print "### CONSTANT = "+str(_constant)
    except:
       print 'WARNING: NO CONSTANT FOR THE FRAME !!!'
       answ=raw_input('do you want to enter the constant interactively  [y/n] ? [y] ')
       if not answ: answ='y'
       if answ=='y':
          zeropoint[filter]=raw_input('Zero point (with filter = '+filter+')? ')
          zeropoint[filter]=float(zeropoint[filter])
          colore[filter]=raw_input('Colour term [BV,VR,RI] (with filter '+filter+')? ')
          colorevalue[filter]=raw_input('Colour term value (with filter '+filter+')? ')
          colorevalue[filter]=float(colorevalue[filter])
          colorerr[filter]=float(0.0)
          zeroerr[filter]=float(0.0)
       else:
          print 'RUN SVSTD.PY ON THE LIST  !!!'
          sys.exit()

       
    try:
       magsn=float(string.split(pyfits.getheader(img)['QUBASN1'])[-3]) 
       magsnerr1=float(string.split(pyfits.getheader(img)['QUBASN1'])[-2]) 
       magsnerr2=float(string.split(pyfits.getheader(img)['QUBASN1'])[-1])
       magerr[filter]=max(magsnerr1,magsnerr2)
       mag[filter]=magsn
       print "### instrumental magnitude header check ......ok "
       print "### MAG = "+str(magsn)
    except:
       print 'ERROR: NO MAGNITUDE FOUND !!!'
       print 'RUN SVSN.PY  !!!'
#       sys.exit()

    try:
       _JD=src.JD(img,_header,_telescope)
       print "### JD header check ......ok "
       print "### JD = "+str(_JD)
       JD[filter]=_JD
    except:
        print 'ERROR: NO JD HEADER FOUND'
        sys.exit()

    magstr[filter]=mag[filter]+2.5*log10(1)-knew[filter]*airmass[filter]
    
for jj in filtri:
    colore_new=(zeropoint[colore[jj][0]]-zeropoint[colore[jj][1]]+(magstr[colore[jj][0]]-magstr[colore[jj][1]]))/(1-(colorevalue[colore[jj][0]]-colorevalue[colore[jj][1]]))
    if abs(colore_new)>=7:
        #print jj,colore[jj],colorevalue[jj]
        #print 'WARNING: colour term quite high '+str(colore_new)+'  !! '
        answ=raw_input('Do you want to use a different colour [y/n] ? [y] ')
        if not answ: answ='y'
        if answ in ['y','Y','yes','Yes','YES']:
            _colore_new=colore_new
            colore_new=raw_input('which is the colour term '+colore[jj]+' of the SN at JD '+str(JD[jj])+' ? ['+str(_colore_new)+'] ? ')
            if not colore_new: colore_new=float(_colore_new)
            else: colore_new=float(colore_new)

    dc0,dc1,dz0,dz1,dm0,dm1=erroremag(zeropoint[colore[jj][0]],zeropoint[colore[jj][1]],magstr[colore[jj][0]],magstr[colore[jj][1]],colorevalue[colore[jj][0]],colorevalue[colore[jj][1]],string.find(colore[jj],jj))

    #####################
    #   color error not considered
    dc0=0.0    #    
    dc1=0.0    #
    #####################
    #    option no zeropoint error
    if _zp:
        dz0,dz1=0.0,0.0
    #####################
    magerr_new[jj]=sqrt((dm0*magerr[colore[jj][0]])**2+(dz0*zeroerr[colore[jj][0]])**2+(dm1*magerr[colore[jj][1]])**2+(dz1*zeroerr[colore[jj][1]])**2+(dc0*colorerr[colore[jj][0]])**2+(dc1*colorerr[colore[jj][1]])**2)
    mag_new[jj]=zeropoint[jj]+colorevalue[jj]*colore_new+magstr[jj]


print listsystem

if len(listsystem)==1:
    system=listsystem[0]
else:
    system='mix'

print '###############################################################################'
JDvec=[]
for i in range(len(filtri)):
   JDvec.append(JD[filtri[i]])
   print '#',str(filtri[i]),' ',str(colore[filtri[i]]),' ',str(zeropoint[filtri[i]]),'  ',str(colorevalue[filtri[i]])
if system==0:
  print '#############'
  print '#       JD           U              B              V              R              I'
  print '# %8.8s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(mean(JDvec)),str(mag_new['U']+0.0005),str(magerr_new['U']+0.0005),str(mag_new['B']+0.0005),str(magerr_new['B']+0.0005),str(mag_new['V']+0.0005),str(magerr_new['V']+0.0005),str(mag_new['R']+0.0005),str(magerr_new['R']+0.0005),str(mag_new['I']+0.0005),str(magerr_new['I']+0.0005))
  print '###############################################################################'
elif system==1:
  print '#############'
  print '#       JD           J              H              K '
  print '# %8.8s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(mean(JDvec)),str(mag_new['J']+0.0005),str(magerr_new['J']+0.0005),str(mag_new['H']+0.0005),str(magerr_new['H']+0.0005),str(mag_new['K']+0.0005),str(magerr_new['K']+0.0005))
  print '###############################################################################'
elif system==2:
  print '#############'
  print '#       JD           u              g              r              i              z'
  print '# %8.8s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(mean(JDvec)),str(mag_new['u']+0.0005),str(magerr_new['u']+0.0005),str(mag_new['g']+0.0005),str(magerr_new['g']+0.0005),str(mag_new['r']+0.0005),str(magerr_new['r']+0.0005),str(mag_new['i']+0.0005),str(magerr_new['i']+0.0005),str(mag_new['z']+0.0005),str(magerr_new['z']+0.0005))
  print '###############################################################################'
elif system=='mix':
    print '###############################################################################'
    for i in range(len(filtri)):
        print '# %8.8s\t%8.8s\t%6.6s %5.5s\t' % (str(filtri[i]),str(JDvec[i]),str(mag_new[filtri[i]]+0.0005),str(magerr[filtri[i]]+0.0005))
    print '###############################################################################'

