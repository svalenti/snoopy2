#!/usr/bin/env python
help ="################################################################ \n\
help = Usage: ecastro.py filename                                      \n \
                   [ -m,--magsel = 16.0]  catalogue selection mag        \n \
                   [-i,--interactive]     mode for coordinate map      \n \
################################################################"
import getopt
import os,sys,shutil,string
from numpy import *

pypath = os.path.expandvars('$HOME')
if not os.path.isfile('login.cl'):
    shutil.copyfile(pypath+'/iraf/login.cl','login.cl')
from pyraf import iraf

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)

magsel = '20.0'
interactive = False   # default interative mode 

if len(sys.argv)==1:
    print help
    sys.exit()
img = sys.argv[1]

options, args = getopt.getopt(sys.argv[2:],'m:,i',['magsel=','interactive'])
for opt,arg in options:
    if opt in ('-m', '--magsel'):      magsel = arg
    if opt in ('-i', '--interactive'): interactive = True

iraf.display(img,1)
iraf.delete('tmp*.*',verify='no')

##########    GET FIELD COORDINATES
RA,DEC=0,0
try:
    coord = iraf.hselect(img,'RA,DEC','yes',Stdout=1)
    RA = float(string.split(coord[0],'\t')[0]) 
    DEC = float(string.split(coord[0],'\t')[1]) 
except:
    print "RA,DEC not found"

ra = iraf.clDms(RA)
dec = iraf.clDms(DEC)
answ = raw_input('>> RA,Dec ['+str(ra)+' '+str(dec)+'] ? ')
if len(answ)>0: 
    ra,dec = string.split(answ)
    RA = iraf.real(ra)
    DEC = iraf.real(dec) 
iraf.printf("RA=%h DEC=%m\n ",RA,DEC)

ff = open('tmp.reg','w')
ff.write(ra+' '+dec+' '+'10 10')
ff.close()
iraf.agetim('tmp.reg','tmpdss',imsurve='dss1@cadc',wcsedit='yes')
iraf.display('tmpdss',2)

answ = raw_input('>> Change cuts (y/n) ? [n]')
if len(answ)>0:
    while string.lower(answ)=='y': 
        z1,z2 = string.split(raw_input('z1 z2 ? '))
        iraf.display('tmpdss',2,z1=z1,z2=z2,zscale='no',zrange='no')
        answ = raw_input('>> Change cuts (y/n) ? [n]')
        if len(answ)==0: answ='n' 

iraf.agetcat('tmpdss','tmp.cat',catalog='usno2@cadc')
iraf.afiltcat('tmp.cat','tmp.catsel',fexpr="mag1<="+str(magsel))
iraf.fields('tmp.catsel','2,3', Stdout='tmp.tv')
iraf.wcsctran('tmp.tv','tmp.pix','tmpdss',inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
iraf.tvmark(2,'tmp.pix',mark="circle",number='yes',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=2)

print '>> Identify (min. 2, preferably 3) reference stars (mark with "a")' 
iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',radii=10,nxoffse=10,nyoffse=10,color=214,txsize=2)
xycoo = iraf.fields('tmp.coo','1,2,13',Stdout=1)
print '>> Identify reference stars'
idcat = []
for i in range(len(xycoo)):
    idcat.append(int(raw_input('Star '+str(i+1)+'= ? '))) 

ff = open('tmp.tv','r')
rr = ff.readlines()
ff.close()
gg = open('tmp.ccmap','w')
fw = []
for i in range(len(idcat)):
    _rr = string.split(rr[idcat[i]-1])
    _x,_y,_fw = string.split(xycoo[i])
    gg.write(_x+' '+_y+' '+_rr[0]+' '+_rr[1]+' \n')
    fw.append(float(_fw))
gg.close()
iraf.hedit(img,'Seeing',mean(fw),add='yes',update='yes',verify='yes')

iraf.ccmap('tmp.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)

iraf.ccfind('tmp.cat','tmp.star',img,lngcolu=2,latcolu=3,lngunit='degrees',usewcs='yes')
iraf.ccmap('tmp.star','tmp.ccdb',images=img,fitgeome='general',xcolum=10,ycolum=11,lngcolum=2,latcolumn=3,lngunit='degrees',update='yes',interact=interactive)

iraf.display(img,1)  

iraf.wcsctran('tmp.tv','tmp.pix3',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
iraf.tvmark(1,'tmp.pix3',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=2)
#iraf.delete("tmp*",verify='no')
