#!/usr/bin/env python

from numpy import array
from numpy import mean
from numpy import zeros
import os,sys,string,re,shutil,getopt
from math import log10
from math import sqrt
from snoopy2 import *
import snoopy2
import pyfits

help ="################################################################ \n\
help = Usage:   ph_new.py  filename             \n\
              -s system   landolt (0) infrared(1), sloan(2)     \n\
              -z          no error in the zeropoint             \n\
################################################################"

def readec(file_in):
    f=file(file_in,'r')
    non_buona=f.readline()
    non_buona=non_buona + f.readline()
    non_buona=non_buona + f.readline()
    non_buona=non_buona + f.readline()
    non_buona=non_buona + f.readline()
    vec1,vec2,vec3,vec4,vec5=[],[],[],[],[]
    s=f.readline()
    while s[0:1]!='#':
        #print 's = '+ s
        c1,c2,c3,c4,c5 = s.split()
        vec1.append(c1)
        vec2.append(float(c2))
        vec3.append(float(c3))
        vec4.append(float(c4))
        vec5.append(float(c5))
        s=f.readline()
	if not s:
	      s='#'
    rest=f.readlines()
    if not rest:
        snname,sn_ap1,sn_ap2,sn_fit,sn_fit_err,sn_art_err=[],[],[],[],[],[]
	non_buona2=''
    else:
	non_buona2=''
    	snname,sn_ap1,sn_ap2,sn_fit,sn_fit_err,sn_art_err=[],[],[],[],[],[]
        for i in rest:
            if i[0]=='#':
                non_buona2=non_buona2+i
            else:
               c1,c2,c3,c4,c5,c6 = i.split()
               snname.append(c1)
               sn_ap1.append(c2)
               sn_ap2.append(c3)
               sn_fit.append(c4)
               sn_fit_err.append(c6)
               sn_art_err.append(c5)
    return vec1,vec2,vec3,vec4,vec5,snname,sn_ap1,sn_ap2,sn_fit,sn_fit_err,sn_art_err,non_buona,non_buona2


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


##########################################################################################################################

if len(sys.argv)==1:
    print help
    sys.exit()


logincl=src.open_program()
#from pyraf import iraf
userPath=os.getcwd()
_telescope=''
system=''
_system=''
_zp=False
check=0


################################################################
options, args = getopt.getopt(sys.argv[2:],'s:,z',['system=','zeroerror'])
for opt,arg in options:
    if opt in ('-s', '--system'): _system = int(arg)
    if opt in ('-z', '--zeroerror'): _zp = True

imglist=src.readlist(sys.argv[1])

_telescope=src.telescope(imglist[0])
print '### TELESCOPE= '+_telescope
kk=src.atmospheric_site2()

bandsvec='UBVRIugrizJHK'
magstdstr,magsnstr,magstdstrerr,magsnstrerr,magstr,zeropoint,zeroerr,colore,colorerr,colorevalue={},{},{},{},{},{},{},{},{},{}
magerr_newsn={}
for i in bandsvec:
    magerr_newsn[i],magstr[i],zeropoint[i],zeroerr[i],colorevalue[i],colorerr[i],magstdstr[i],magsnstr[i],magstdstrerr[i],magsnstrerr[i]=0,0,0,0,0,0,0,0,0,0
    colore[i]='none'

##########################################################################

filtri=[]
airmass={}
exptime={}
JD={}

answ=raw_input('phot or fit (p/f) [f]? ' )
if not answ: answ='f'
if answ=='f':
    colonna=2
    colonnaerr=colonna+1
else:
    colonna=1
    colonnaerr=colonna+2

tmp=[]
for img in imglist:
	if img[-5:]=='.fits': img=img[:-5]
	tmp.append(int(string.split(os.popen('wc -l '+img+'.ec').read())[0]))

for i in tmp:
   if i!=tmp[0]:
       print '###### Error: files .ec have different lengths !!!'
       src.close_program(logincl)
    
listsystem=[]
for imglong in imglist:
  if imglong:
    if imglong[-5:]=='.fits': img=imglong[:-5]
    if _system not in [1,2,0]:
        system=src.check_system(_telescope,imglong,Stdout=False)
    else:
        system=_system
    if system not in listsystem: listsystem.append(system)
    _header=src.read_parameter(_telescope,system)
    _object=src.objects(imglong,_header,_telescope)
    _date=src.date(imglong,_header,_telescope)
    _filter=src.filter(imglong,_header,_telescope)
    filter=src.filtername(_telescope,_filter,system)
    if filter=='unknown':
        filter=_filter
    filtri.append(filter)
    print '### FILTER = '+str(filter)
    try:
       _airmass=float(src.airmass(imglong,_header,_telescope))
       airmass[filter]=_airmass
       print "### AIRMASS header check ......ok "
       print "### AIRMASS = "+str(_airmass)
    except:
        print 'ERROR: NO AIRMASS HEADER FOUND'
        src.close_program(logincl)

    try:
       _exptime=src.exptime(imglong,_header,_telescope)
       exptime[filter]=_exptime
       print "### exptime header check ......ok "
       print "### EXPTIME = "+str(_exptime)
    except:
       print 'ERROR: NO EXPTIME HEADER FOUND'
       src.close_program(logincl)

    try:
       _constant=pyfits.getheader(imglong)['QUBACONS']
       _zeropoint=float(string.split(pyfits.getheader(imglong)['QUBACONS'])[1]) 
       _colore=string.split(pyfits.getheader(imglong)['QUBACONS'])[2] 
       _colorevalue=float(string.split(pyfits.getheader(imglong)['QUBACONS'])[3]) 
       zeropoint[filter]=raw_input('zero point (with filter = '+filter+')?['+str(_zeropoint)+'] ')
       if not zeropoint[filter]: zeropoint[filter]=_zeropoint
       zeropoint[filter]=float(zeropoint[filter])
       colore[filter]=raw_input('colour term [BV,VR,RI] (with filter '+filter+')?['+str(_colore)+']  ')
       if not colore[filter]: colore[filter]=_colore
       colorevalue[filter]=raw_input('colour term value (with filter '+filter+')?['+str(_colorevalue)+']  ')
       if not colorevalue[filter]: colorevalue[filter]=_colorevalue
       colorevalue[filter]=float(colorevalue[filter])
       try:
           _zeroerr=float(string.split(pyfits.getheader(imglong)['QUBACONS'])[5]) 
           zeroerr[filter]=float(_zeroerr)
           _colorerr=float(string.split(pyfits.getheader(imglong)['QUBACONS'])[6]) 
           colorerr[filter]=float(_colorerr)
       except:
           zeroerr[filter]=float(0.0)
           colorerr[filter]=float(0.0)

       print "### CONSTANTS header check ......ok "
       print "### CONSTANT = "+str(_constant)
    except:
       print 'WARNING: NO CONSTANT FOR THE FRAME !!!'
       answ=raw_input('Do you want to enter the constant interactively  [y/n] ? [y] ')
       if not answ: answ='y'
       if answ=='y':
          zeropoint[filter]=raw_input('Zero point (with filter = '+filter+')? ')
          zeropoint[filter]=float(zeropoint[filter])
          colore[filter]=raw_input('Colour term [BV,VR,RI] (with filter '+filter+')? ')
          colorevalue[filter]=raw_input('Colour term value (with filter '+filter+')? ')
          colorevalue[filter]=float(colorevalue[filter])
          zeroerr[filter]=float(0.0)
          colorerr[filter]=float(0.0)
       else:
          print 'RUN SVSTD.PY ON THE LIST  !!!'
          sys.exit()
#          src.close_program(logincl)
       
    try:
       _JD=src.JD(imglong,_header,_telescope)
       print "### JD header check ......ok "
       print "### JD = "+str(_JD)
       JD[filter]=_JD
    except:
        print 'ERROR: NO JD HEADER FOUND'
        sys.exit()
#        src.close_program(logincl)

    try:
        identstars={}
        identsn={}
        num,fwhm,ph,fitmag,fiterr,snname,sn_ap1,sn_ap2,sn_fit,sn_fit_err,sn_art_err,ppp,ccc=readec(img+'.ec')
        for i in range(len(num)):
            identstars[num[i]]=[fwhm[i],ph[i],fitmag[i],fiterr[i]]
	if snname:    
           for i in range(len(snname)):
               identsn[snname[i]]=[sn_ap1[i],sn_ap2[i],sn_fit[i],sn_fit_err[i],sn_art_err[i]]
    except:
        print '### ERROR: problem reading .ec file (filter '+str(filter)+') !'

    magstdstr[filter]=[]
    magsnstr[filter]=[]
    magstdstrerr[filter]=[]
    magsnstrerr[filter]=[]
#    for i in identstars:
    for i in num:    
        magstdstr[filter].append(float(identstars[i][colonna])+2.5*log10(1)-kk[filter]*airmass[filter])
        magstdstrerr[filter].append(float(identstars[i][colonnaerr]))
    if 	identsn:
      for i in range(len(identsn)):
          try:
              magsnstr[filter].append(float(identsn[snname[i]][colonna])+2.5*log10(1)-kk[filter]*airmass[filter])
              magsnstrerr[filter].append(max(float(identsn[snname[i]][colonnaerr]),float(identsn[snname[i]][colonnaerr+1])))
          except:
              magsnstr[filter].append(9999)
              magsnstrerr[filter].append(0.0)
  else:
  	print '####'
  	print '#### WARNING: empty space in the list !!'
	print '####'
#################     sequence stars
#
mag_new={}
magerr_new={}
for jj in filtri:
    mag_new[jj]=[]
    magerr_new[jj]=[]
    colore_new=[]
    if string.count(str(filtri),jj)==1:
        try: 
            if  magstdstr[colore[jj][0]] and magstdstr[colore[jj][1]]:
                for ii in range(len(magstdstr[jj])):
                    colore_new.append((zeropoint[colore[jj][0]]-zeropoint[colore[jj][1]]+(magstdstr[colore[jj][0]][ii]-magstdstr[colore[jj][1]][ii]))/(1-(colorevalue[colore[jj][0]]-colorevalue[colore[jj][1]])))
                    dc0,dc1,dz0,dz1,dm0,dm1=erroremag(zeropoint[colore[jj][0]],zeropoint[colore[jj][1]],magstdstr[colore[jj][0]][ii],magstdstr[colore[jj][1]][ii],colorevalue[colore[jj][0]],colorevalue[colore[jj][1]],string.find(colore[jj],jj))
                    mag_new[jj].append(zeropoint[jj]+colorevalue[jj]*colore_new[ii]+magstdstr[jj][ii])
                    #####################
                    #   color error not considered
                    dc0=0.0    #    
                    dc1=0.0    #
                    #####################  no zero error
                    if _zp:
                        dz0,dz1=0.0,0.0
                    #####################
                    magerr_new[jj].append(sqrt((dm0*magstdstrerr[colore[jj][0]][ii])**2+(dz0*zeroerr[colore[jj][0]])**2+(dm1*magstdstrerr[colore[jj][1]][ii])**2+(dz1*zeroerr[colore[jj][1]])**2+(dc0*colorerr[colore[jj][0]])**2+(dc1*colorerr[colore[jj][1]])**2))
            else:
                print ''
                print 'WARNING: COLOUR CORRECTION NOT APPLIED for '+str(jj)+'  '+str(colore[jj])+' !!!!!'
                print 'this is ok if you have only one band or you insert wrong color term !!!!!'
                print ''
                for ii in range(len(magstdstr[jj])):
                    mag_new[jj].append(zeropoint[jj]+magstdstr[jj][ii])
                    magerr_new[jj].append(sqrt((dm0*magstdstrerr[colore[jj][0]][ii])**2+(dz0*zeroerr[colore[jj][0]])**2+(dm1*magstdstrerr[colore[jj][1]][ii])**2+(dz1*zeroerr[colore[jj][1]])**2))
            for ii in range(len(magstdstr[jj])):
                if mag_new[jj][ii]>=99:
                    mag_new[jj][ii]=9999
                    magstdstrerr[jj][ii]=0.0
                    magerr_new[jj][ii]=0.0
        except:
            mag_new[jj]=list(array(magstdstr[filtri[0]])-array(magstdstr[filtri[0]])+9999)
            magstdstrerr[jj]=list(array(magstdstr[filtri[0]])-array(magstdstr[filtri[0]])-0.0005)
            magerr_new[jj]=list(array(magstdstr[filtri[0]])-array(magstdstr[filtri[0]])-0.0005)
            print '#  '+jj

#
##########     Supernova
#
if identsn:
    mag_newsn={}
    magerr_newsn={}
    for jj in filtri:
      if string.count(str(filtri),jj)==1:
          mag_newsn[jj]=[]
          magerr_newsn[jj]=[]
          colore_newsn=[]
          try:
            if  magsnstr[colore[jj][0]] and magsnstr[colore[jj][1]]:
                for ii in range(len(magsnstr[jj])):
                    colore_newsn.append((zeropoint[colore[jj][0]]-zeropoint[colore[jj][1]]+(magsnstr[colore[jj][0]][ii]-magsnstr[colore[jj][1]][ii]))/(1-(colorevalue[colore[jj][0]]-colorevalue[colore[jj][1]])))
                    dc0,dc1,dz0,dz1,dm0,dm1=erroremag(zeropoint[colore[jj][0]],zeropoint[colore[jj][1]],magsnstr[colore[jj][0]][ii],magsnstr[colore[jj][1]][ii],colorevalue[colore[jj][0]],colorevalue[colore[jj][1]],string.find(colore[jj],jj))
                    mag_newsn[jj].append(zeropoint[jj]+colorevalue[jj]*colore_newsn[ii]+magsnstr[jj][ii])
                    #####################
                    #   color error not considered
                    dc0=0.0    #    
                    dc1=0.0    #
                    #####################
                    #     no zero error  
                    if _zp:
                        dz0,dz1=0.0,0.0
                    #####################
                    magerr_newsn[jj].append(sqrt((dm0*magsnstrerr[colore[jj][0]][ii])**2+(dz0*zeroerr[colore[jj][0]])**2+(dm1*magsnstrerr[colore[jj][1]][ii])**2+(dz1*zeroerr[colore[jj][1]])**2+(dc0*colorerr[colore[jj][0]])**2+(dc1*colorerr[colore[jj][1]])**2))
            else:
                print ''
                print 'WARNING: COLOUR CORRECTION NOT APPLIED for '+str(jj)+'  '+str(colore[jj])+' !!!!!'
                print 'this is ok if you have only one band or you insert wrong color term !!!!!'
                print ''
                for ii in range(len(magsnstr[jj])):
                    mag_newsn[jj].append(zeropoint[jj]+magsnstr[jj][ii])
            if mag_newsn[jj][ii]>=99:
                mag_newsn[jj][ii]=9999
          except:
              mag_newsn[jj]=list(array(magsnstr[filtri[0]])-array(magsnstr[filtri[0]])+9999)
              magsnstrerr[jj]=list(array(magsnstr[filtri[0]])-array(magsnstr[filtri[0]])-0.0005)
              magerr_newsn[jj]=list(array(magsnstr[filtri[0]])-array(magsnstr[filtri[0]])-0.0005)

JDvec=[]
for i in range(len(filtri)):
   JDvec.append(JD[filtri[i]])
   print '#',str(filtri[i]),' ',str(colore[filtri[i]])
   print '# ',str(zeropoint[filtri[i]]),'  ',str(colorevalue[filtri[i]])
print '# '

filtrisort='UuBgVrRiIzJHK'
if len(listsystem)==1:
    system=listsystem[0]
else:
    system='mix'
    filtri2=''
    for i in filtrisort: 
        if i in filtri:
            filtri2=filtri2+i
    filtri=filtri2

for band in filtrisort:
    if band not in filtri:
        mag_new[band]=list(zeros(len(num)))
        magstdstrerr[band]=(zeros(len(num)))
        magerr_new[band]=(zeros(len(num)))
        if snname:
            mag_newsn[band]=(zeros(len(snname)))
            magerr_newsn[band]=(zeros(len(snname)))

for band in filtrisort:
    for i in range(0,len(mag_new[band])): 
        if mag_new[band][i]<=4 or mag_new[band][i]>=99:
            mag_new[band][i]=9999.
            magstdstrerr[band][i]=0.0
            magerr_new[band][i]=0.0
    if identsn:
        for i in range(0,len(mag_newsn[band])): 
            if mag_newsn[band][i]<=4 or mag_newsn[band][i]>=99:
                mag_newsn[band][i]=9999.
                magerr_newsn[band][i]=0.0

if system==0:
   print '###############################################################################'
   print '# '
   print '#       name          U              B              V              R              I'
   for i in range(len(num)):
#       print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t'    % (str(num[i]),str(mag_new['U'][i]+0.0005),str(magstdstrerr['U'][i]+0.0005),str(mag_new['B'][i]+0.0005),str(magstdstrerr['B'][i]+0.0005),str(mag_new['V'][i]+0.0005),str(magstdstrerr['V'][i]+0.0005), \
#       str(mag_new['R'][i]+0.0005),str(magstdstrerr['R'][i]+0.0005),str(mag_new['I'][i]+0.0005),str(magstdstrerr['I'][i]+0.0005))
       print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t'    % (str(num[i]),str(mag_new['U'][i]+0.0005),str(magerr_new['U'][i]+0.0005),str(mag_new['B'][i]+0.0005),str(magerr_new['B'][i]+0.0005),str(mag_new['V'][i]+0.0005),str(magerr_new['V'][i]+0.0005), \
       str(mag_new['R'][i]+0.0005),str(magerr_new['R'][i]+0.0005),str(mag_new['I'][i]+0.0005),str(magerr_new['I'][i]+0.0005))

   if snname:
      for i in range(len(snname)):
          print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(snname[i]),str(mag_newsn['U'][i]+0.0005),str(magerr_newsn['U'][i]+0.0005),str(mag_newsn['B'][i]+0.0005),str(magerr_newsn['B'][i]+0.0005),str(mag_newsn['V'][i]+0.0005),str(magerr_newsn['V'][i]+0.0005), \
	  str(mag_newsn['R'][i]+0.0005),str(magerr_newsn['R'][i]+0.0005),str(mag_newsn['I'][i]+0.0005),str(magerr_newsn['I'][i]+0.0005))    

elif system==1:
   print '###############################################################################'
   print '# '
   print '#       name          J              H              K'
   for i in range(len(num)):
       print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(num[i]),str(mag_new['J'][i]+0.0005),str(magerr_new['J'][i]+0.0005),str(mag_new['H'][i]+0.0005),str(magerr_new['H'][i]+0.0005), \
        str(mag_new['K'][i]+0.0005),str(magerr_new['K'][i]+0.0005))
   if snname:
      for i in range(len(snname)):
          print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(snname[i]),str(mag_newsn['J'][i]+0.0005),str(magerr_newsn['J'][i]+0.0005),str(mag_newsn['H'][i]+0.0005), \
	  str(magerr_newsn['H'][i]+0.0005),str(mag_newsn['K'][i]+0.0005),str(magerr_newsn['K'][i]+0.0005))    
elif system==2:
   print '###############################################################################'
   print '# '
   print '#       name          u              g              r              i              z'
   for i in range(len(num)):
       print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(num[i]),str(mag_new['u'][i]+0.0005),str(magerr_new['u'][i]+0.0005), str(mag_new['g'][i]+0.0005),str(magerr_new['g'][i]+0.0005), str(mag_new['r'][i]+0.0005),str(magerr_new['r'][i]+0.0005),str(mag_new['i'][i]+0.0005), str(magerr_new['i'][i]+0.0005),str(mag_new['z'][i]+0.0005),str(magerr_new['z'][i]+0.0005))
   if snname:
      for i in range(len(snname)):
          print ' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t' % (str(snname[i]),str(mag_newsn['u'][i]+0.0005),str(magerr_newsn['u'][i]+0.0005),str(mag_newsn['g'][i]+0.0005),str(magerr_newsn['g'][i]+0.0005),str(mag_newsn['r'][i]+0.0005), str(magerr_newsn['r'][i]+0.0005),str(mag_newsn['i'][i]+0.0005),str(magerr_newsn['i'][i]+0.0005),str(mag_newsn['z'][i]+0.0005),str(magerr_newsn['z'][i]+0.0005))

elif system=='mix':
   print '###############################################################################'
   print '# '
   title='#     name   '
   for band in filtri:
       title=title+'   '+str(band)+'   '
   print title
   for i in range(len(num)):
       stringa='%6.6s\t' % (str(num[i]))
       for band in filtri:
           stringa=stringa+'%6.6s %5.5s\t' % (str(mag_new[band][i]+0.0005),str(magstdstrerr[band][i]+0.0005)) 
       print stringa
   for i in range(len(snname)):
       stringa2='%6.6s\t' % (str(snname[i]))
       for band in filtri:
           stringa2=stringa2+'%6.6s %5.5s\t' % (str(mag_newsn[band][i]+0.0005),str(magsnstrerr[band][i]+0.0005))
       print stringa2
print '# '
print '###############################################################################'


photometric=raw_input('Do you know if the night was photometric ? (p) photometric or (n) non-photometric [n]')
if not photometric: photometric = 'n'
if photometric == 'n': photometric = 'noSTD'
else: photometric = 'ph'

if colonna==2:
    _type='fit.doc'
else:
    _type='ph.doc'

output=_object+'_'+_date+'_'+_telescope+'_'+_type

ff=open(output,'w')
ff.write('***  %7.7s  %6.6s  %5.5s\n' % (str(mean(JDvec)),str(_telescope),str(photometric)))
if system==0:
  ff.write('#       name          U              B              V              R              I\n')
  ff.write('#\n')
  for i in range(len(num)):#  magerr_newsn
        ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(num[i]),str(mag_new['U'][i]+0.0005),str(magerr_new['U'][i]+0.0005),str(mag_new['B'][i]+0.0005),str(magerr_new['B'][i]+0.0005),str(mag_new['V'][i]+0.0005),str(magerr_new['V'][i]+0.0005),str(mag_new['R'][i]+0.0005),str(magerr_new['R'][i]+0.0005),str(mag_new['I'][i]+0.0005),str(magerr_new['I'][i]+0.0005)))
  if snname:
     for i in range(len(snname)):
         ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(snname[i]),str(mag_newsn['U'][i]+0.0005),str(magerr_newsn['U'][i]+0.0005),str(mag_newsn['B'][i]+0.0005),str(magerr_newsn['B'][i]+0.0005),str(mag_newsn['V'][i]+0.0005),str(magerr_newsn['V'][i]+0.0005),str(mag_newsn['R'][i]+0.0005),str(magerr_newsn['R'][i]+0.0005),str(mag_newsn['I'][i]+0.0005),str(magerr_newsn['I'][i]+0.0005)))
elif system==1:
  ff.write('#       name          J              H             K\n')
  ff.write('#\n')
  for i in range(len(num)):
        ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(num[i]),str(mag_new['J'][i]+0.0005),str(magerr_new['J'][i]+0.0005),str(mag_new['H'][i]+0.0005),str(magerr_new['H'][i]+0.0005),str(mag_new['K'][i]+0.0005),str(magerr_new['K'][i]+0.0005)))
  if snname:
     for i in range(len(snname)):
         ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(snname[i]),str(mag_newsn['J'][i]+0.0005),str(magerr_newsn['J'][i]+0.0005),str(mag_newsn['H'][i]+0.0005),str(magerr_newsn['H'][i]+0.0005),str(mag_newsn['K'][i]+0.0005),str(magerr_newsn['K'][i]+0.0005)))
elif system==2:
  ff.write('#       name          u              g              r              i              z\n')
  ff.write('#\n')
  for i in range(len(num)):
        ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(num[i]),str(mag_new['u'][i]+0.0005),str(magerr_new['u'][i]+0.0005),str(mag_new['g'][i]+0.0005), str(magerr_new['g'][i]+0.0005),str(mag_new['r'][i]+0.0005),str(magerr_new['r'][i]+0.0005),str(mag_new['i'][i]+0.0005),str(magerr_new['i'][i]+0.0005),str(mag_new['z'][i]+0.0005),str(magerr_new['z'][i]+0.0005)))
  if snname:
     for i in range(len(snname)):
         ff.write(' %6.6s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t%6.6s %5.5s\t\n' % (str(snname[i]),str(mag_newsn['u'][i]+0.0005),str(magerr_newsn['u'][i]+0.0005),str(mag_newsn['g'][i]+0.0005),str(magerr_newsn['g'][i]+0.0005), \
	 str(mag_newsn['r'][i]+0.0005),str(magerr_newsn['r'][i]+0.0005),str(mag_newsn['i'][i]+0.0005),str(magerr_newsn['i'][i]+0.0005),str(mag_newsn['z'][i]+0.0005),str(magerr_newsn['z'][i]+0.0005)))
elif system=='mix':
   title='#     name   '
   for band in filtri:
       title=title+'   '+str(band)+'   '
   ff.write(title)
   ff.write('#\n')
   for i in range(len(num)):
       stringa='%6.6s\t' % (str(num[i]))
       for band in filtri:
           stringa=stringa+'%6.6s %5.5s\t' % (str(mag_new[band][i]+0.0005),str(magerr_new[band][i]+0.0005)) 
       stringa=stringa+'\n'
       ff.write(stringa)
   for i in range(len(snname)):
       stringa2='%6.6s\t' % (str(snname[i]))
       for band in filtri:
           stringa2=stringa2+'%6.6s %5.5s\t' % (str(mag_newsn[band][i]+0.0005),str(magerr_newsn[band][i]+0.0005))
       stringa2=stringa2+'\n'
       ff.write(stringa2)
ff.close()

