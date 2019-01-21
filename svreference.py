#!/usr/bin/env python

help ="################################################################ \n\
help = Usage: svreference.py filename                                   \n\
                   [ -m,--magsel = 16.0]  catalogue selection mag       \n\
                   [-i,--interactive]     mode for coordinate map       \n\
                   [-b,--boxsize = '10 10'] boxsize for ref image in arcmins default 10 by 10  \n\
################################################################"

import getopt,subprocess
import os,sys,shutil,string,re,time
import snoopy2
from snoopy2 import *
from numpy import abs
from numpy import mean
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
    
if len(sys.argv)==1:
    print help
    sys.exit()

logincl=src.open_program()

from pyraf import iraf

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)

magsel = '20.0'
boxsize='10 10'
interactive = False   # default interative mode 
userPath=os.getcwd()

img = sys.argv[1]

options, args = getopt.getopt(sys.argv[2:],'m:,i,b:,p:,',['magsel=','interactive','userPath='])
for opt,arg in options:
    if opt in ('-m', '--magsel'):      magsel = arg
    if opt in ('-i', '--interactive'): interactive = True
    if opt in ('-p', '--path'): userPath = arg
    if opt in ('-b', '--boxsize'):      boxsize = arg

_z1,_z2,goon=src.display_image(img,1,'','',True)
if not goon: src.close_program(logincl)    

src.delete('tmp*.*')

def xpa(arg):
    subproc = subprocess.Popen('xpaset -p ds9 '+arg,shell=True)
    subproc.communicate()

ds9 = subprocess.Popen("ps -U"+str(os.getuid())+"|grep -v grep | grep ds9",shell=True,stdout=subprocess.PIPE).stdout.readlines()
if len(ds9)== 0 :
    subproc = subprocess.Popen('ds9',shell=True)
    time.sleep(3)

##########    GET FIELD COORDINATES
check=0
_telescope=src.telescope(img)
if _telescope:
   print '### TELESCOPE= '+_telescope
   _system=src.check_system(_telescope,img,Stdout=True)
   if _system in [1,2,0]:
       check=src.check_tel(img,_telescope,_system)

if check:
    _header=src.read_parameter(_telescope,_system)
    _RA=src.RA(img,_header,_telescope)
    _DEC=src.DEC(img,_header,_telescope)
    print _RA,_DEC
else:
    try:
        _RA=pyfits.getheader(img)['RA']
        _DEC=pyfits.getheader(img)['DEC']
        print _RA,_DEC
    except:
        _RA=''
        _DEC=''

print ''
print '### note: if Dec is positive, do not write the "+"'
print' (e.g. 10:10:13 13:34:30 or 10.2345  13.543)'
print ''
_ra = raw_input('>> Ra ['+str(_RA)+'] ? ')
if not _ra: _ra=_RA
_dec = raw_input('>> Dec [ '+str(_DEC)+'] ? ')
if not _dec: _dec=_DEC

if string.count(str(_ra),':')>=1:
    ra1,ra2,ra3=string.split(_ra,':')
    ran=(float(ra1)+(float(ra2)+float(ra3)/60.)/60.)
else:
    ran=_ra
#    ra11=string.replace(_ra,':',' ')
#    dec11=string.replace(dec,':',' ')

if string.count(str(_dec),':')>=1:
    dec1,dec2,dec3=string.split(_dec,':')
#    if float(dec1)>=0:
    if string.count(str(dec1),'-')==0:
        decn=float(dec1)+(float(dec2)+float(dec3)/60.)/60.
    else:    
        decn=(-1)*(abs(float(dec1))+(float(dec2)+float(dec3)/60.)/60.)
else:
    decn=_dec

try:
    sistema=os.popen('uname -s').read()
    if 'Darwin' in sistema:
        environ='mac'
    elif 'Linux' in sistema:
        environ='linux'
    elif 'altro' in sistema:
        environ='altro'
    else:
        environ=''
except:
        environ=''

if not environ:
    print 'Warning: operative system not recognised !!!!'
    print 'you can try with one of the following but it is not sure it will work !!!!'
    environ=raw_input('which OS [linux/mac/exit] ? [mac]')
    if not environ: environ='mac'
    if environ=='exit':
        sys.exit()


#################################
if environ=='mac':
    survey, mag = 'dss','$\Rmag'
elif environ=='linux':
    survey, mag = 'dss','$\Rmag'

# survey = 'dss' working on ubuntu erkki

xpa('frame 2')
xpa('zscale contrast 0.1')
xpa(survey+' size '+boxsize)
xpa(survey+' coord '+str(ran)+' '+str(decn))
xpa('zoom to fit')
xpa('catalog ua2')
xpa('catalog server cds')
xpa('catalog size 10 10 arcmin')
xpa('catalog save "'+userPath+'/tmp.cat"')
xpa('catalog filter "'+mag+'<'+str(magsel)+'"')
xpa('catalog symbol shape text')

if environ=='mac':
    xpa('catalog save sb "'+userPath+'/tmp.catsel"')
elif environ=='linux':
    xpa('catalog save sb "'+userPath+'/tmp.catsel"')

xpa('raise')

iraf.fields('tmp.catsel','1,2', Stdout='tmp.tv')

answ='no'

while answ!='yes':
    print '>> Identify (minimum 2, preferably 3) reference stars (mark with "a") in your image' 
    iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
    iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',radii=10,nxoffse=10,nyoffse=10,color=214,txsize=2)
    xycoo = iraf.fields('tmp.coo','1,2,13',Stdout=1)
    print '>> Identify reference stars'
    idcat = []
    for i in range(len(xycoo)):
        idcat.append(int(raw_input('Star '+str(i+1)+'= ? '))) 

    ff = open('tmp.catsel')
    rr = ff.readlines()
    ff.close()
    for i in range(len(rr)):
        if  rr[i].split('\t')[0] =='_RAJ2000': istart = i
    _racat,_decat = [],[]
    for r in rr[istart+2:]:
        _racat.append(float(r.split('\t')[0]))
        _decat.append(float(r.split('\t')[1]))

    gg = open('tmp.ccmap','w')
    fw = []
    for i in range(len(idcat)):
        _x,_y,_fw = xycoo[i].split()
        gg.write(_x+' '+_y+' '+str(_racat[idcat[i]-1])+' '+\
                 str(_decat[idcat[i]-1])+' \n')
        fw.append(float(_fw))
    gg.close()

    src.updateheader(img,0,'Seeing',mean(fw))
    gg = open('tmp.coo','w')
    for i in range(len(_racat)):
        gg.write(str(_racat[i])+' '+str(_decat[i])+' \n')
    gg.close()

    iraf.ccmap('tmp.ccmap','tmp.ccdb',images=img,fitgeome='rscale',\
                   xcolum=1,ycolum=2,lngcolum=3,latcolumn=4,lngunit='degrees',\
                   update='yes',interact=False,verbose='yes')
    tmp = iraf.ccfind('tmp.coo','tmp.star',img,lngcolu=1,latcolu=2,\
                          lngunit='degrees',usewcs='yes',verbose='yes',Stderr=1)

    print '>>>', tmp[-2]
    if int(tmp[-2].split()[1])<=3: 
        print "!!! Error: stars are not matched !!!"
        sys.exit()

    iraf.ccmap('tmp.star','tmp.ccdb',solution='tmp.sol',images=img,\
                   fitgeome='rscale',xcolum=3,ycolum=4,lngcolum=1,\
                   latcolumn=2,lngunit='degrees',\
                   update='yes',interact=False)

    iraf.ccmap('tmp.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)

    _z1,_z2,goon=src.display_image(img,1,_z1,_z2,False)
    if not goon: src.close_program(logincl)    

    src.delete("tmp.pix3")
    iraf.wcsctran('tmp.tv','tmp.pix3',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    iraf.tvmark(1,'tmp.pix3',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=2)

    answ=raw_input(' ok ? [y/n] [y]')
    if not answ: answ='yes'
    if answ in ['YES','Yes','Y','y','yes']:
        answ='yes'
    else:
        src.delete("tmp.coo")
        src.delete("tmp.star")

src.delete("tmp*")

_z1,_z2,goon=src.display_image(img,1,_z1,_z2,False)
if not goon: src.close_program(logincl)    

print '>> Identify sequence stars (mark with "a")' 
iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
iraf.wcsctran('tmp.coo','tmp.tv2',img,inwcs='logical',units='degrees degrees',outwcs='world',columns='1 2',formats='%10.8f %10.8f')


sss=iraf.fields('tmp.coo','1,2', Stdout=1)
lll=iraf.fields('tmp.tv2','1,2', Stdout=1)

nn=len(sss)

testo="# BEGIN CATALOG HEADER \n\
# nfields "+str(nn)+" \n\
#     id 1 0 c INDEF %15s \n\
#     ra 2 0 d degrees %10.5f \n\
#     dec 3 0 d degrees %10.5f \n\
#     magV  4 0 r INDEF %6.2f \n\
#     magBV 5 0 r INDEF %6.2f \n\
#     magUB 6 0 r INDEF %6.2f \n\
#     magVR 7 0 r INDEF %6.2f \n\
#     magRI 8 0 r INDEF %6.2f \n\
# END CATALOG HEADER \n\
#\n"

testo1="  \tV\tBV\tUB\tVR\tRI\t"

for i in range(len(sss)):
    testo= testo+'  '+str(i+1)+'\t'+lll[i+2]+testo1+'\n'

src.delete('file_coordinate.list')
gg = open('file_coordinate.list','w')
gg.write(testo)
gg.close()

print testo
print '#################################################'
print '### EDIT the file file_coordinate.list and include the V absolute magnitude and the color terms !!!!'
print '### This should be done using output from qubstd'
print '### See Manual step 1 for further details.'
src.delete("tmp*")

#src.close_program(logincl)
