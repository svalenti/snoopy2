#!/usr/bin/env python

from numpy import mean
from numpy import array
from numpy import compress
from numpy import std
from numpy import average
from numpy import median
from numpy import abs
from snoopy2 import *
import getopt
import os,sys,shutil,string,re
import snoopy2
import glob
from astropy.io import fits as pyfits


####################  READ OPTION ###########################################
bb=[]
aa=glob.glob(snoopy2.__path__[0]+'/coordinate_std/optical/*list')
for i in aa:
    bb.append(i[len(snoopy2.__path__[0]+'/coordinate_std/optical/'):-5])

help ="################################################################       \n\
help = Usage:   svsn.py filename                                              \n\
                input            filelist (iraf format)                       \n\
            [-p --psf psffile or psflist]  (iraf format)                      \n\
            [-i,--interactive]     interactive mode (check the measurement)   \n\
            [-r,--recenter]     interactive mode (check the measurement)   \n\
   	    [-s,--system value]     Specific photometric system  0 opt 1 inf 2 sloan   \n\
	    [-z,--size value] half size of stamp around the object      \n\
	    [-c,--coordinate]      position of the object from list           \n\
	    [-b,--background_region value]  region of background interpolation     \n\
	                                 (until the border of the stamp)      \n\
	    [-x --xorder  value]   xorder for background                   \n\
	    [-y --yorder  value]   yorder for background	              \n\
	    [-n,--iteration value]     number of background iterations            \n\
	    [-R,-- ra]     (12:12:23.3) or 12.25      \n\
	    [-D,-- dec]    (10:12:14.4) or 10.345       \n\
################################################################"


_interactive = False   # default interative mode
coordinate = False 
arterr=0
_fb0=3
_size0=7
_xbgord0=2
_ybgord0=2
_numiter0=3
_size=''
_fb=''
_xord=''
_yord=''
_numiter=''
_recenter= False
out = ''
system = ''
_ra=''
_dec=''
subdirectory=['optical/','infrared/','sloan/']

if len(sys.argv)<=1:
    print help
    sys.exit() 


logincl=src.open_program()
from pyraf import iraf

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)

imglist=src.readlist(sys.argv[1])

psflista=''

options, args = getopt.getopt(sys.argv[2:],'p:,i,s:,c,z:,b:,x:,y:,n:,r,R:,D:',['psflist=','interactive','system=','coordinate','size=','background_region=','xorder=','yorder=','iteraction=','recenter','ra=','dec='])
for opt,arg in options:
    if opt in ('-p', '--psf'): psflista = arg
    if opt in ('-i', '--interactive'): _interactive = True
    if opt in ('-s', '--system'): 
        system = int(arg)
        if system not in [0,1,2]:
            print 'ERROR: system not recognised !!'
            src.close_program(logincl)

    if opt in ('-R', '--ra'): 
        _ra = arg
        try:
            _ra=float(_ra)
            print 'ra= '+str(_ra)
        except:
            try:
                _ra=float(string.split(_ra,':')[0])+float(string.split(_ra,':')[1])/60.+float(string.split(_ra,':')[2])/3600.
                print 'ra= '+str(_ra)
            except:
                print 'ERROR: ra format not recognised !!'
                src.close_program(logincl)

    if opt in ('-D', '--dec'): 
        _dec = arg    
        try:
            _dec=float(_dec)
            print 'dec= '+str(_dec)
        except:
            try:
                if '-' not in str(_dec):
                    _dec=float(string.split(_dec,':')[0])+float(string.split(_dec,':')[1])/60.+float(string.split(_dec,':')[2])/3600.
                else:
                    _dec=(-1)*(abs(float(string.split(_dec,':')[0]))+float(string.split(_dec,':')[1])/60.+float(string.split(_dec,':')[2])/3600.)
                print 'dec= '+str(_dec)
            except:
                print 'ERROR: dec format not recognised !!'
                src.close_program(logincl)

    if opt in ('-c', '--coordinate'): coordinate = True
    if opt in ('-r', '--recenter'): _recenter= True
    if opt in ('-z', '--size'): 
                           try:  
			      _size = float(arg)
			      print ' size = '+str(_size)
			   except: 
			         print 'ERROR: size format not recognised !!'
			         src.close_program(logincl)
    if opt in ('-b', '--background_region'): 
                           try: 
			      _fb=float(arg) 
			      print ' length of bg region = '+str(_fb)
                           except:  
			         print 'ERROR: length of bg region format not recognised !!'
			         src.close_program(logincl)
    if opt in ('-x', '--xorder'): 
                           try:
			       _xord = int(arg)
			       print ' xorder = '+str(_xord)
                           except:  
			         print 'ERROR: xorder format not recognised !!'
			         src.close_program(logincl)			   
    if opt in ('-y', '--yorder'): 
                           try:
			       _yord = int(arg)
			       print ' xorder = '+str(_xord)
                           except:  
			         print 'ERROR: yorder format not recognised !!'
			         src.close_program(logincl)
    if opt in ('-n', '--iteration'): 
                           try:  
			      _numiter = int(arg)
			      print 'number of iterations = '+str(_numiter)
                           except:  
			         print 'ERROR: iteration format not recognised !!'
			         src.close_program(logincl)
####################################### check ##########################
check=0
_telescope=src.telescope(imglist[0])
print '### TELESCOPE= '+_telescope
if not system:
     system=src.check_system(_telescope,imglist[0],Stdout=True)
else:
    if system==0:
        print '###  OPTICAL '
        print '##### UBVRI  #######'
        print '###'
    elif system==1:
        print '###   INFRARED '
        print '##### JHK  #######'
        print '###'
    elif system==2:
        print '###   SLOAN '
        print '##### ugriz  #######'
        print '###'
check=src.check_tel(imglist[0],_telescope,system)
if check==0:
        print '####################### '
        print "#### Error with the header !!!!"
        print '### telescope not in the list '
        print '### if you want to continueanyway,'
        print '### run "svother.py" to correct the header '
        print '####################### '

##########################################################################

if psflista!='':
    if psflista[0]=='@':
        ff = open(psflista[1:])
        files = ff.readlines()
        psflist = []
        for ff in files: psflist.append(ff[:-1])
    elif ',' in psflista:
        psflist = string.split(psflista,sep=',')
    else:
        psflist = [psflista]
######################## SET   IRAF   PARAMETERS  #######################
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.imcoords(_doprint=0)
from iraf import digiphot
from iraf import daophot
from iraf import ptools

#datamin = -200
#datamax0 = 60000
zmag = 0 

#iraf.digiphot.daophot.datapars.datamin = datamin
#iraf.digiphot.daophot.datapars.datamax = datamax0
iraf.digiphot.daophot.photpars.zmag = zmag

##########################        DEFINE FITSN  #########################################
def fitsn(img,coordlist,_recenter,fwhm0,a4,original,sn,residual,displ):
    src.delete("apori")
    src.delete(img+".sn.mag")
###############################
    from pyraf import iraf
    import string,re,os,sys
    from numpy import log10
    iraf.imcoords(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)
    from iraf import digiphot
    from iraf import daophot
    from iraf import ptools
    a1 = int(fwhm0)
    a2 = int(2.*fwhm0+.5)
    a3 = int(3.*fwhm0+.5)
    a4 = int(4.*fwhm0+.5)
    ap = str(a1)+","+str(a2)+","+str(a3)
######################
    if _recenter:        answ='no'
    else:                answ='yes'
################################
    iraf.noao.digiphot.daophot.photpars.zmag = 0
    iraf.noao.digiphot.daophot.datapars.readnoi = 1.4
    iraf.noao.digiphot.daophot.datapars.epadu = 13
    iraf.noao.digiphot.daophot.datapars.datamin = -100
    iraf.noao.digiphot.daophot.datapars.datamax = 51000

    iraf.noao.daophot.fitskypars.annulus=a3
    iraf.noao.daophot.photpars.apertures = ap
    iraf.noao.digiphot.daophot.datapars.exposure = _header['hed_exptime'] #'exptime'
    iraf.noao.digiphot.daophot.datapars.airmass = _header['hed_airmass']  #'airmass'
    iraf.noao.digiphot.daophot.datapars.filter = 'filter2'
    iraf.noao.digiphot.daophot.daopars.psfrad = a4
    iraf.noao.digiphot.daophot.daopars.fitrad = fwhm0  
    iraf.noao.digiphot.daophot.daopars.sannulus = int(a4)
    iraf.noao.digiphot.daophot.daopars.recenter =  answ
    iraf.noao.digiphot.daophot.daopars.fitsky = 'yes'
    #  fitskypars.salgorithm = "constant"
    #  fitskypars.skyvalue = 0
    print '\n### recentering: '+str(answ)
    if displ:
        iraf.noao.digiphot.daophot.phot(original,coordlist,"apori",veri='no')   
        iraf.noao.digiphot.daophot.phot(sn,coordlist,img+".sn.mag",veri='no')   
    else:
        iraf.noao.digiphot.daophot.phot(original,coordlist,"apori",veri='no',verb='no')   
        iraf.noao.digiphot.daophot.phot(sn,coordlist,img+".sn.mag",veri='no',verb='no')   

    src.delete(img+".sn.als")
    iraf.allstar(sn,img+".sn.mag",img+".psf",img+".sn.als","",residual,veri='no',verb='no')
    src.delete("snfit.fits")
    iraf.imarith(sn+'.fits',"-",residual+'.fits',"snfit.fits")
#    src.pyimarith(sn+'.fits',"-",residual+'.fits',"snfit.fits")
    src.delete("skyfit.fits")
    iraf.imarith(original+'.fits',"-","snfit.fits","skyfit.fits")
#    src.pyimarith(original+'.fits',"-","snfit.fits","skyfit.fits")
    iraf.txsort(img+".sn.als","ID")
    tmptbl = iraf.txdump(img+".sn.als","mag,merr,xcenter,ycenter",expr='yes', Stdout=1)
    magerr,fitmag,centx,centy=[],[],[],[]
    for i in tmptbl:
        try:
            fitmag.append(float(string.split(i)[0]))
        except:
            fitmag.append(string.split(i)[0])
        try:
            magerr.append(float(string.split(i)[1]))
        except:
            magerr.append(string.split(i)[1])
        centx.append(float(string.split(i)[2]))
        centy.append(float(string.split(i)[3]))
    tmptbl=iraf.txdump("apori","mag",expr='yes', Stdout=1)
    apori1,apori2,apori3=[],[],[]
    for i in tmptbl:
        try:
            apori1.append(float(string.split(i)[0]))
        except:
            apori1.append(string.split(i)[0])
        try:            
            apori2.append(float(string.split(i)[1]))
        except:
            apori2.append(string.split(i)[1])
        try:
            apori3.append(float(string.split(i)[2]))
        except:
            apori3.append(string.split(i)[2])
            
    iraf.txsort(img+".sn.mag","YCENTER")
    tmptbl=iraf.txdump(img+".sn.mag","mag,magerr",expr='yes', Stdout=1) 

    if displ:
        print "********************************************************************"
        print "ID <apmag on original>  <apmag on bgsubt> fitmag truemag err_fit"         
        print "     ",a1,"       ",a2,"      ",a3,"        ",a1,"     ",a2,"     ",a3 

    apmag1,apmag2,apmag3,truemag=[],[],[],[]
    for i in range(len(tmptbl)):
        apmag1.append(string.split(tmptbl[i])[0])
        apmag2.append(string.split(tmptbl[i])[1])
        apmag3.append(string.split(tmptbl[i])[2])
        try:
            truemag.append(fitmag[i]+float(apco0))
        except:
            truemag.append('INDEF')
        if displ:
            print i,apori1[i],apori2[i],apori3[i],apmag1[i],apmag2[i],apmag3[i],fitmag[i],truemag[i],magerr[i]
    if displ:
        print "********************************************************************"

    print fitmag,tmptbl,str(apco0)
    if displ:
        _tmp1,_tmp2,goon=src.display_image(original+'.fits',1, z11, z22, False, _xcen=.25, _ycen=.25, _xsize=.3, _ysize=.3)
        #iraf.display(original,1,fill='yes',xcen=.25,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11) 
        z01 = z11-midpt
        z02 = z22-midpt 
        s1 = 1
        s2 = -int(fwhm0)
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" ORIGINAL")
        ff.close()    
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
        _tmp1,_tmp2,goon=src.display_image('snfit.fits',1, z01, z02, False, _xcen=.25, _ycen=.75, _xsize=.3, _ysize=.3, _erase='no')
        #iraf.display("snfit",1,erase='no',fill='yes',xcen=.25,ycen=.75,xsize=.3,ysize=.3, zscal='no',zrang='no', z2=z02, z1=z01)
        src.delete("tmptbl")
        tmptbl0=iraf.txdump(img+".sn.als","xcen,ycen",expr='yes',Stdout=1)
        ff=open('tmptbl','w')
        for i in tmptbl0:
            ff.write(i+'\n')
        ff.close()    
        lra = int((2*float(size)*float(fwhm0))*2)
        iraf.tvmark(1,"tmptbl",autol='no',mark="circle", number='yes',nyoffset=lra,radi=a2,txsize=2,inter='no')
        s1 = 1
        s2 = -1*int(fwhm0)
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" FITTED")
        ff.close()    
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
        _tmp1,_tmp2,goon=src.display_image('skyfit.fits',1, z11, z22, False, _xcen=.75, _ycen=.25, _xsize=.3, _ysize=.3, _erase='no')
        #iraf.display("skyfit",1,erase='no',fill='yes',xcen=.75,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11)
        s1 = 1
        s2 = -1*int(fwhm0)
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" RESIDUAL")
        ff.close()    
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
    return apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy

###########################  DEFINE MANUAL CHANGING  #########################################

def manusn(dmag0,apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy):
   from numpy import zeros
   fdmag = 10**(-0.4*float(dmag0))
   src.delete(",_snfit.fit?,skyfit.fit?,_snfit.ar?")
   if truemag[0]=='INDEF':
       print 'ACTUNG'
       magerr[0]=0.0
       src.delete('test.fits')
       os.system('echo '+str(iraf.field(img+'.sn.coo', field='1,2', Stdout=1)[0])+' '+dmag0+' > dddd')
       iraf.addstar("snfit","dddd",img+".psf","_snfit",nstar=1,veri='no',simple='yes',verb='yes')
       src.delete('dddd')
   else:
       src.pyimarith("snfit.fits","*",fdmag,"_snfit.fits")   
   src.pyimarith("original.fits","-","_snfit.fits","skyfit.fits")
   _tmp1,_tmp2,goon=src.display_image('original.fits',1, z11, z22, False, _xcen=.25, _ycen=.25, _xsize=.3, _ysize=.3)
   #iraf.display("original",1,fill='yes',xcen=.25,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1= z11) 
   z01 = z11-midpt
   z02 = z22-midpt 
   s1 = 1
   s2 = -int(fwhm0)
   src.delete("tmptbl")
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" ORIGINAL") 
   ff.close()
   iraf.tvmark(1,"tmptbl",autol='yes',mark="none",inter='no',label='yes',txsize=2)

   #iraf.display("_snfit",1,erase='no',fill='yes',xcen=.25,ycen=.75,xsize=.3,ysize=.3, zscal='no', zrang='no', z2=z02, z1=z01 )
   _tmp1,_tmp2,goon=src.display_image('_snfit.fits',1, z01, z02, False, _xcen=.25, _ycen=.75, _xsize=.3, _ysize=.3, _erase='no')
   src.delete("tmptbl")
   tmptbl=iraf.txdump(img+".sn.als","xcen,ycen",expr='yes', Stdout=1)
   ff=open('tmptbl','w')
   for i in tmptbl:
       ff.write(i) 
   ff.close()  
   lra = int((2*size*fwhm0)*2)
   iraf.tvmark(1,"tmptbl",autol='no',mark="circle", number='yes',nyoffset=lra,radi=a2,txsize=2,inter='no')
   s1 = 1
   s2 = -int(fwhm0)
   src.delete("tmptbl")
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" FITTED") 
   ff.close()  
   iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
   _tmp1,_tmp2,goon=src.display_image('skyfit.fits',1, z11, z22, False, _xcen=.75, _ycen=.25, _xsize=.3, _ysize=.3, _erase='no')
   #iraf.display("skyfit",1,erase='no',fill='yes',xcen=.75,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11) 
   s1 = 1
   s2 = -int(fwhm0)
   src.delete("tmptbl")
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" RESIDUAL") 
   ff.close()
   iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)

   newmag=list(zeros(len(truemag)))
   for i in range(0,len(truemag)):
       try:
           newmag[i]=float(truemag[i])+float(dmag0)
       except:
           newmag[i]=float(dmag0)
           magerr[i]=0.0
   print truemag
   print newmag
   #newmag = list(array(truemag)+float(dmag0))   

   print "***************************************************************************" 
   print "#id  x_ori   y_ori     x     y    ap_ori ap_bgsub  fit_mag  err_art  err_fit" 
   for i in range(len(fitmag)):
       print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(newmag[i]),"  ",str(arterr),"  ",str(magerr[i])
   print "**************************************************************************"

   return apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,newmag


###################   DEFINE ERROR     ######################################

def errore(img,coordlist,size,truemag,fwhm0,leng0,_interactive,_numiter):
    if not _numiter: _numiter=3
    dartf=100
    while dartf>=size-1:
        if _interactive:
            artfac0 =raw_input('>>> Dispersion of artificial star positions (in units of FWHM) [1] ')
            if not artfac0: artfac0=1
        else:
            artfac0=1
        try:  
            artfac0=float(artfac0)
            if float(artfac0)>=size-1:
                print '!!! WARNING: '+str(artfac0)+' too large (max '+str(size)+'- 1)'
                print 'try again....'
            else:
                dartf = artfac0
        except:
            print '#### WARNING: '+str(artfac0)+' should be a number !!!!'
            print 'try again....'
    src.delete("tmpar?")

    src.delete('artskyfit.fits')
    os.system('cp skyfit.fits artskyfit.fits')
    i=0
    tmpart=[]
    while i<=8:
        src.delete("reserr.fit?,artbg.fit?,artstar.fit?,artres.fit?,artfit.fit?")
        artrad = fwhm0/2.
        #artseed = artseed+1234
        artx =  int(i/3.)-1
        if i<=2: arty = artx+i
        if 3<=i<=5: arty = artx-1+i-3
        if i>=6: arty = artx-2+i-6

        ff=open(img+".sn.coo",'r')
        ss=ff.readline()
        ff.close()
        xbb=float(string.split(ss)[0])
        ybb=float(string.split(ss)[1])

        xbb=xbb+artx*fwhm0*artfac0
        ybb=ybb+arty*fwhm0*artfac0
        #print xbb,ybb

        src.delete(coordlist)
        ff=open(coordlist,'w')
        ff.write(str(xbb)+'  '+str(ybb)+'  '+str(truemag[0])+"  1")
        ff.close()
                        
        xb1 = int(float(xbb)-fwhm0*float(leng0)/2)
        xb2 = int(float(xbb)+fwhm0*float(leng0)/2) 
        yb1 = int(float(ybb)-fwhm0*float(leng0)/2)
        yb2 = int(float(ybb)+fwhm0*float(leng0)/2)
        sec="1 "+str(xb1)+" 1 "+str(nay)+'\n'
        sec=sec+str(xb2)+' '+str(nax)+" 1 "+str(nay)+'\n'
        sesc=sec+str(xb1)+' '+str(xb2)+" 1 "+str(yb1)+'\n'
        sec=sec+str(xb1)+' '+str(xb2)+' '+str(yb2)+' '+str(nay)+'\n'
        ff = open('sec','w')
        ff.write(sec)
        ff.close()

        src.delete("reserr.ar?")
        src.delete("artlist.ma?")
        src.delete("artsky.fit?")
        src.delete("artbg.fit?")
        src.delete("artbgs.fit?")
        src.delete("artsn.fit?")
        src.delete("artres.fit?")
        src.delete("artlist.al?")

        iraf.addstar("artskyfit",coordlist,img+".psf","reserr",nstar=1,veri='no',simple='yes',verb='no')  # reserr = skyfit + artificial star
########
        inp = "artbg.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]"
        out = "artsky.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]" 
        iraf.imsurfit("reserr","artbg",xorder=xbgord0,yorder=ybgord0,regions="section",section="sec")  
        midpt=float(iraf.imstat("artbg",field="mean", Stdout=1)[1])
        iraf.imcopy('reserr.fits','artsky.fits')
        iraf.imcopy(inp,'artbgs.fits')
        iraf.imcopy("artbgs.fits",out)
        #iraf.imarith("reserr","-","artsky","artsn",calctype="r",pixtype="r",verb='no') 
        src.pyimarith("reserr.fits","-","artsky.fits","artsn.fits") 
        src.pyimarith("artsn.fits","+",midpt,"artsn.fits")

        artap1,artap2,artap3,artmag1,artmag2,artmag3,artfitmag,arttruemag,artmagerr,artcentx,artcenty=fitsn(img,coordlist,_recenter,fwhm0,a4,'reserr','artsn','artres','')

        for ii in range(0,_numiter):
            src.delete("reserr.ar?")
            src.delete("artlist.ma?")
            src.delete("artsky.fit?")
            src.delete("artbg.fit?")
            src.delete("artbgs.fit?")
            src.delete("artsn.fit?")
            src.delete("artres.fit?")
            src.delete("artlist.al?")

            iraf.imsurfit("skyfit","artbg",xorder=xbgord0,yorder=ybgord0,regions="section",section="sec")  
            midpt=float(iraf.imstat("artbg",field="mean", Stdout=1)[1])
            iraf.imcopy("reserr.fits","artsky.fits")
            iraf.imcopy(inp,"artbgs.fits")
            iraf.imcopy("artbgs.fits",out)

            #iraf.imarith("reserr","-","artsky","artsn",calctype="r",pixtype="r",verb='no') 
            src.pyimarith("reserr.fits","-","artsky.fits","artsn.fits")
            src.pyimarith("artsn.fits","+",midpt,"artsn.fits")
            artap1,artap2,artap3,artmag1,artmag2,artmag3,artfitmag,arttruemag,artmagerr,artcentx,artcenty=fitsn(img,coordlist,_recenter,fwhm0,a4,'reserr','artsn','artres','')
####### 
        if i==0: era='yes'
        else:    era='no'
        artx = .5+.25*artx
        arty = .5+.25*arty
        _tmp1,_tmp2,goon=src.display_image('skyfit.fits',1,'', '', False, _xcen=artx, _ycen=arty, _xsize=.25, _ysize=.25, _erase=era)
        # aaa=iraf.display("skyfit",1,fill='yes',erase=era,xcen=artx,ycen=arty,xsize=.25,ysize=.25,Stdout=1)
	try:
            tmpart.append(float(arttruemag[0]))
#	   tmpart.append(float(iraf.txdump("artlist.als","mag",expr='yes', Stdout=1)[0]))
	except: pass
        i=i+1

    for i in tmpart:  print i 
    
    print " ########## "
    try:
        media=mean(array(tmpart))
        arterr=std(array(tmpart))
        arterr2=std(compress((average(tmpart)-std(tmpart)<array(tmpart))&(array(tmpart)<average(tmpart)+std(tmpart)),array(tmpart)))
    except:
        media=0
        arterr=0
        arterr2=0
    print '### average = %6.6s \t arterr= %6.6s ' % (str(media),str(arterr))
    print '###  %6.6s \t (error at 1 sigma rejection) ' % (str(arterr2))
    src.delete("reserr.fit?,artbg.fit?,artstar.fit?,artres.fit?,artfit.fit?,artskyfit.fit?")
    src.delete("reserr.ar?")
    src.delete("artlist.co?")
    return arterr2,arterr

#################################    CHECK HEADER    ######################################


warn='##################################\n'

for imglong in imglist:
  if imglong:
    if imglong[-5:]=='.fits': img=imglong[:-5]
    else: img=imglong
    _header=src.read_parameter(_telescope,system)
    _gain=src.gain(imglong,_header,_telescope)
    _ron=src.ron(imglong,_header,_telescope)
    _datamin=src.ccdmin(imglong,_header,_telescope)
    _datamax=src.ccdmax(imglong,_header,_telescope)
    iraf.noao.digiphot.daophot.datapars.readnoi = _ron
    iraf.noao.digiphot.daophot.datapars.epadu = _gain
    iraf.noao.digiphot.daophot.datapars.datamin = _datamin
    iraf.noao.digiphot.daophot.datapars.datamax = _datamax

    print _gain,_ron,_datamax,_datamin    

    _filter=src.filter(imglong,_header,_telescope)
    filter=src.filtername(_telescope,_filter,system)
    if filter=='unknown':
        print '### WARNING: filter not recognised in this photometric system '
        filter=_filter
    print '##########  '+str(filter)+'  ##############'

    check=0
    if psflista=='':
        check=1
        if not os.path.isfile(img+'.psf.fits'):
                check=0
        else:
                psfimage=img+'.psf.fits'
    else:
        psfimage=psflist[imglist.index(img+'.fits')]
        print '##############################################'
        print psfimage+' <--> '+img+'.fits'
        print '##############################################'
        answ=raw_input('### Is this correct ? [y/n] ? [y] ')
        if not answ: answ='y'
        if answ=='n' or answ=='no':
            print '#### Choose the right psf files and try again !' 
            src.close_program(logincl)
        else:
            check=1
            

    if check==0:
        print '### ERROR: NO PSF FILE  !!!'
        print '###  psf image not found or   '
        print '### the filter is not the same !!'
        print '### run svpsf.py and try again or '
        print '### use the "-p option" with the right psf file  '
    else:
        print "### psf check ......ok "

    try:
        fwhm0=float(pyfits.getheader(imglong)['QUBAFWHM'])
        print "### FWHM header check ......ok "
        print "### FWHM = "+str(fwhm0)
    except:
        print 'WARNING: NO SEEING HEADER KEYWORD FOUND'
        answ=raw_input('### write the FWHM or press enter  ?  ')
        if not answ:
            print '#### run qubapsf and try again !!'
            src.close_program(logincl)
        try:
            fwhm0=float(answ)
            print "### FWHM header check ......ok "
            print "### FWHM = "+str(fwhm0)
        except:
            print '#### WARNING: FWHM not valid !!'
	    src.close_program(logincl)
    try:
        apco0=float(pyfits.getheader(imglong)['QUBAAPCO'])
        print "### aperture correction header check ......ok "
        print "### APCO = "+str(apco0)
    except:
        print 'WARNING: NO APERTURE CORRECTION HEADER KEYWORD FOUND'
        answ=raw_input('### Write the aperture Correction or press enter  ?  ')
        if not answ:
            print '#### run qubapsf and try again !!'
	    src.close_program(logincl)
        try:
            apco0=float(answ)
            print "### Aperture Correction header check ......ok "
            print "### Aperture Correction = "+str(apco0)
        except:
            print '#### WARNING: FWHM not valid !!'
	    src.close_program(logincl)

        
    a1 = int(fwhm0)
    a2 = int(2.*fwhm0+.5)
    a3 = int(3.*fwhm0+.5)
    a4 = int(4.*fwhm0+.5)
    ap = str(a1)+","+str(a2)+","+str(a3)
    
#  fitskypars.salgorithm = "constant"
#  fitskypars.skyvalue = 0

    iraf.daophot.fitskypars.annulus=a3
    iraf.daophot.photpars.apertures = ap
    iraf.digiphot.daophot.datapars.exposure = _header['hed_exptime']
    iraf.digiphot.daophot.datapars.airmass = _header['hed_airmass']
    iraf.digiphot.daophot.datapars.filter = _header['hed_filter1']

######################    POINT TO SN #######################################

    iraf.set(stdimage='imt2048')
    _z1,_z2='',''
    if _interactive:
        _z1,_z2,goon=src.display_image(img+'.fits',1,'','',True,_xsize=1,_ysize=1)
    else:
        _z1,_z2,goon=src.display_image(img+'.fits',1,'','',False)

    if _ra and _dec:
        os.system('rm -rf tmp.*')
        ff=open('tmp.tv','w')
        ff.write(str(_ra*15.)+' '+str(_dec))
        ff.close()
        iraf.wcsctran('tmp.tv','tmp.pix',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
        iraf.tvmark(1,'tmp.pix',mark="circle",number='yes',radii=10,nxoffse=5,nyoffse=5,color=214,txsize=2)
        xx0,yy0=string.split(iraf.fields('tmp.pix','1,2',Stdout=1)[2])
        print xx0,yy0
        os.system('rm -rf tmp.*')
    else:
        if coordinate:
            xx0,yy0=src.sn_coordinate(img+'.fits',_interactive)
        else:
            xx0=''
    if xx0 and not _interactive:
        x=xx0
        y=yy0
    else:
      repeat='y'
      while repeat=='y':
        print "_____________________________________________"
        print "  MARK SN REGION WITH - x -, EXIT  - q -"
        try:
            x,y,value=string.split(iraf.imexamine(img+'.fits',1,wcs='logical',xformat='',yformat='',Stdout=1)[0])
            src.delete("tmplabel")
            ff = open('tmplabel','w')
            ff.write(str(x)+' '+str(y)+' -1'+' \n')
            ff.close()
            iraf.tvmark(1,'tmplabel',autol='no',mark="cross",inter='no',label='no',txsize=4)
            repeat = raw_input('### repeat selection ? [y/n] ? [n] ')
            if not repeat: repeat='n'
            elif repeat=='yes': repeat='y' 
            elif repeat=='YES': repeat='y' 
        except:
            x=y=value=0
            print '### WARNING: SN REGION NOT SELECTED !!!'
            repeat = raw_input('### repeat selection ? [y/n] ? [n] ')
            if not repeat:    repeat='n'

            if repeat in ['Y','y','YES','yes','Yes']:
                repeat='y' 
            else:
                sys.exit()
    
    if _telescope=='prompt':
             prompttele = pyfits.getheader(img+'.fits')['OBSERVAT']
             if string.count(prompttele,'2'): _prompt='prompt2'
             elif string.count(prompttele,'3'): _prompt='prompt3'
             elif string.count(prompttele,'4'): _prompt='prompt4'
             elif string.count(prompttele,'5'): _prompt='prompt5'
             else: 
                 print 'WARNING: telescope not found !! '
                 _prompt=''
             imgin_cor='home$coordinate_std/'+_telescope+'/'+_prompt+'.fits'

    if not _interactive or _size:
         if _size: size=_size
	 else: size=_size0
         x1 = int(float(x)-size*fwhm0)
         if x1<=0: x1=0
         x2 = int(float(x)+size*fwhm0)
         y1 = int(float(y)-size*fwhm0)
         if y1<=0: y1=0
         y2 = int(float(y)+size*fwhm0)
         src.delete("original.fits")    
         iraf.imcopy(imglong+"["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]","original.fits")
    else:
      repeat='n'
      while repeat=='n':            
         size = raw_input('Size of the cut frame (fwhm) ['+str(_size0)+'] ? ')
         if not size : size = _size0
         else: size=int(size)
         x1 = int(float(x)-size*fwhm0)
         if x1<=0: x1=0
         x2 = int(float(x)+size*fwhm0)
         y1 = int(float(y)-size*fwhm0)
         if y1<=0: y1=0
         y2 = int(float(y)+size*fwhm0)
    
         src.delete("original.fits")
         iraf.imcopy(imglong+"["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]","original.fits")
         iraf.set(stdimage='imt512')
         _tmp1,_tmp2,goon=src.display_image('original.fits',1,'' ,'' , False, _xsize=.5, _ysize=.5)
         #sss=iraf.display('original',1,xsize=.5,ysize=.5,fill='yes', Stdout=1)
         repeat = raw_input('### ok ? [y/n] ? [y] ')
         if not repeat: repeat='y'
         elif repeat=='no': repeat='n'

    if _telescope=='prompt' and _prompt:
             _trimsection=src.trimsec(imglong,_header,_telescope)
             if _trimsection:
                 print _trimsection
                 src.delete("prompt_trimmed.fits")
                 iraf.imcopy (imgin_cor+_trimsection,"prompt_trimmed.fits")
                 
                 imgin_cor="prompt_trimmed.fits"
             src.delete("prompt_original.fits")
             iraf.imcopy(imgin_cor+"["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]","prompt_original.fits")
    _z11,_z22,good=src.display_image('original.fits', 1, '', '', False, _xsize=.5, _ysize=.5)       
    
    if _interactive:
      answ = 'y'
      answ= raw_input(">>>>> Cuts OK [y/n] [y]?")
      if not answ: answ='y'
      elif answ=='no': answ='n' 

      while answ == 'n':  
         z11 = raw_input('>>> z1 = ? ['+str(_z11)+'] ? ')
         z22 = raw_input('>>> z2 = ? ['+str(_z22)+'] ? ')
         if not z11: z11=_z11
         else: z11=float(z11)
         if not z22: z22=_z22
         else: z22=float(z22)
         #sss=iraf.display('original',1,fill='yes', zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
         _z11,_z22,goon=src.display_image('original.fits',1, z11, z22, False)
         answ= raw_input(">>>>> Cuts OK [y/n] [y]?")
         if not answ: answ='y'
         elif answ=='no': answ='n' 

    z11=float(_z11)
    z22=float(_z22)

    if not _interactive and xx0:
        _dimension=string.split((string.split(iraf.imheader('original',Stdout=1)[0],']')[0]),'[')[1]
        aa,bb=string.split(_dimension,',')
        aa,bb=float(aa)/2,float(bb)/2
        _vec=str(aa)+'  '+str(bb)+'  1'
        vector=[_vec]
        ff = open('tmplabel','w')
        for i in vector:
            ff.write(i+' \n')
        ff.close()
	iraf.tvmark(1,'tmplabel',autol='no',mark="circle", radii=10, inter='no',label='no', number='yes',\
	pointsize=20, txsize=2,color =204)
        os.system('cp tmplabel '+img+'.sn.coo')
    else:
      answ0='n'
      while answ0=='n':
        _tmp1,_tmp2,goon=src.display_image('original.fits',1,z11,z22,False)
        #iraf.display('original',1,fill='yes', zrange='no', zscale='no', z1=z11, z2=z22)
        print "   ",str(z11),str(z22)
        print "__________________________________________________"
        print "IDENTIFY SN AND CO-STARS(S) WITH - x -, EXIT - q -"
        print "__________________________________________________"
        print " 1 1 'ID. SN AND CO-STAR(S) WITH -x- EXIT -q-'"
        src.delete("tmplabel")
        vector=iraf.imexamine('original.fits',1,wcs='logical',xformat='',yformat='',use_display='no', Stdout=1)
        if string.count(vector[0],'z1')==1: vector=vector[1:]
 
        ff = open('tmplabel','w')
        for i in vector:
            ff.write(i+' \n')
        ff.close()
	iraf.tvmark(1,'tmplabel',autol='no',mark="circle", radii=10, inter='no',label='no', number='yes',\
	pointsize=20, txsize=2,color =204)
        os.system('cp tmplabel '+img+'.sn.coo')
        answ0= raw_input(">>>>> SN AND CO-STARS(S) IDENTIFICATIONS OK [y/n] [y]?")
        if not answ0: answ0='y'
        elif answ0=='no': answ0='n' 

############################ BACKGROUND FIT   ###############################

#    if not _interactive or _fb:  
#    else:
    print ' ************  background fit **********************'
    answ0 = 'n'
    while answ0 == 'n':
        src.delete("sky.fits,bg.fits,bgs.fits,sn.fits,residual.fits")
        src.display_image('original.fits',1,z11,z22,False, _xsize=.5, _ysize=.5)
        #iraf.display("original.fits",1,fill='yes',xsize=.5,ysize=.5, zrange='no', zscale='no', z1=z11, z2=z22)
        nax=int(pyfits.getheader('original.fits')['NAXIS1'])
        nay=int(pyfits.getheader('original.fits')['NAXIS2'])
        #iraf.imget("original","i_naxis1")
        #nax = int(iraf.imget.value)
        #iraf.imget("original","i_naxis2")
        #nay = int(iraf.imget.value)

        if len(vector)==1:
            xb,yb,value=string.split(vector[0])

            checkleng0='yes'
            while checkleng0=='yes':
	        if not _interactive or _fb:
		   if _fb: leng0=_fb
		   else:  leng0=_fb0
		else:
                   leng0 = raw_input('>>> length of square for background in units of FWHM [3] ? ')
                   if not leng0: leng0=3
                try: 
                    float(leng0)
                    checkleng0='no'
                except: 
                    print 'WARNING: the FWHM should be a number !!!!'
                    checkleng0=='yes'

            iraf.tvmark(1,img+".sn.coo",auto='no',mark="rectangle",length=int(fwhm0*float(leng0)),inter='no', color=204) 
            xb1 = int(float(xb)-fwhm0*float(leng0)/2)
            xb2 = int(float(xb)+fwhm0*float(leng0)/2) 
            yb1 = int(float(yb)-fwhm0*float(leng0)/2)
            yb2 = int(float(yb)+fwhm0*float(leng0)/2)
            sec="1 "+str(xb1)+" 1 "+str(nay)+'\n'
            sec=sec+str(xb2)+' '+str(nax)+" 1 "+str(nay)+'\n'
            sec=sec+str(xb1)+' '+str(xb2)+" 1 "+str(yb1)+'\n'
            sec=sec+str(xb1)+' '+str(xb2)+' '+str(yb2)+' '+str(nay)+'\n'
            ff = open('sec','w')
            ff.write(sec)
            ff.close()

            inp = "bg.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]"
            out = "sky.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]" 
            
            checkorder='yes'
            while checkorder=='yes':
	        if not _interactive or _xord and _yord:
		   if _xord and _yord: 
		       xbgord0=_xord
		       ybgord0=_yord
		   else:  
		       xbgord0=_xbgord0
		       ybgord0=_ybgord0
		else:	    
	          xbgord0 = raw_input('>>> Order of function in x for bg fit ['+str(_xbgord0)+'] ? ')
                  if not xbgord0: xbgord0=_xbgord0
                  else: _xbgord0=xbgord0
                  ybgord0 = raw_input('>>> Order of function in y for bg fit ? ['+str(_ybgord0)+'] ? ')
                  if not ybgord0: ybgord0=_ybgord0
                  else: _ybgord0=ybgord0
                try:
                    float(xbgord0)
                    float(ybgord0)
                    checkorder='no'
                except:
		    print 'WARNING: value not valid !!'  
		    checimsurfitkorder='yes'
            iraf.imsurfit("original","bg",xorder=xbgord0,yorder=ybgord0,regions="section",section="sec")  
        else:
            src.delete("tmplabel")
            src.delete("tmptbl")
            ff = open('tmplabel','w')
            ff.write('')
            ff.close()            
            print ">>>  Mark corners of bg-region with  >b<, exit  >q<"
            iraf.tvmark(1,"tmplabel",autol='no',mark="none",inter='no',label='yes',txsize=2, color=204)
            iraf.tvmark(1,"",logfile="tmptbl",autol='yes',mark="cross",inter='yes', color=204)
            ff = open('tmptbl','r')
            ss=ff.readlines()
            ff.close()            
            xb1=int(float(string.split(ss[-2])[0]))
            yb1=int(float(string.split(ss[-2])[1]))
            xb2=int(float(string.split(ss[-1])[0]))
            yb2=int(float(string.split(ss[-1])[1]))
            sec="1 "+str(xb1)+" 1 "+str(nay)+'\n'
            sec=sec+str(xb2)+' '+str(nax)+" 1 "+str(nay)+'\n'
            sec=sec+str(xb1)+' '+str(xb2)+" 1 "+str(yb1)+'\n'
            sec=sec+str(xb1)+' '+str(xb2)+' '+str(yb2)+' '+str(nay)+'\n'
            ff = open('sec','w')
            ff.write(sec)
            ff.close()
            
            inp = "bg.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]"
            out = "sky.fits["+str(xb1)+":"+str(xb2)+","+str(yb1)+":"+str(yb2)+"]" 

            checkorder='yes'
            while checkorder=='yes':
	       if not _interactive or _xord and _yord:
		   if _xord and _yord: 
		       xbgord0=_xord
		       ybgord0=_yord
		   else:  
		       xbgord0=_xbgord0
		       ybgord0=_ybgord0
	       else:  
                  xbgord0 = raw_input('>>> Order of function in x for bg fit ['+str(_xbgord0)+'] ? ')
                  if not xbgord0: xbgord0=_xbgord0
                  else: _xbgord0=xbgord0
                  ybgord0 = raw_input('>>> Order of function in y for bg fit ['+str(_ybgord0)+'] ? ')
                  if not ybgord0: ybgord0=_ybgord0
                  else: _ybgord0=ybgord0
               try:
                        float(xbgord0)
                        float(ybgord0)
                        checkorder='no'
               except:  
      		    print 'WARNING: value not valid !!'  
	            checkorder='yes'
	      
            iraf.imsurfit("original","bg",xorder=xbgord0,yorder=ybgord0,regions="sections",sections="sec")  

        midpt=float(iraf.imstat("bg",field="mean", Stdout=1)[1])

        iraf.imcopy("original.fits","sky.fits")
        iraf.imcopy(inp,"bgs.fits")
        iraf.imcopy("bgs.fits",out)

        #iraf.imarith("original","-","sky","sn",calctype="r",pixtype="r") 
        src.pyimarith("original.fits","-","sky.fits","sn.fits")
        src.pyimarith("sn.fits","+",midpt,"sn.fits")


        answ0 = 'y'
        print answ0

        _tmp1,_tmp2,goon=src.display_image('original.fits',1, z11, z22, False, _xcen=.25, _ycen=.25, _xsize=.3, _ysize=.3)
        #iraf.display("original",1,fill='yes',xcen=.25,ycen=.25,xsize=.3,ysize=.3,zscal='no',zrang='no',z1=z11,z2=z22)  
        s1 = 1
        s2 = -int(fwhm0)
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" ORIGINAL")
        ff.close()
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
        _tmp1,_tmp2,goon=src.display_image('sky.fits',1, z11, z22, False, _xcen=.25, _ycen=.75, _xsize=.3, _ysize=.3, _erase='no')
        #iraf.display("sky",1,erase='no',fill='yes',xcen=.25,ycen=.75,xsize=.3,ysize=.3,zscal='no', zrang='no' ,z2=z22, z1=z11) 
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" BACKGROUND_FIT") 
        ff.close()
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
        _tmp1,_tmp2,goon=src.display_image('sn.fits',1, z11, z22, False, _xcen=.75, _ycen=.25, _xsize=.3, _ysize=.3, _erase='no')
        #iraf.display("sn",1,erase='no',fill='yes',xcen=.75,ycen=.25,xsize=.3,ysize=.3,zscal='no', zrang='no', z1=z11, z2=z22) 
        src.delete("tmptbl")
        ff=open('tmptbl','w')
        ff.write(str(s1)+' '+str(s2)+" STARS") 
        ff.close()
        iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)

        if not _interactive: 
	   answ0='y'
	else:
          answ0= raw_input(">>> Background fit OK [y/n] [y] ?")
          if not answ0: answ0='y'
          elif answ0=='no': answ0='n' 
	
####################################    FITSN        ###################################  

    apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy=fitsn(img,img+'.sn.coo',_recenter,fwhm0,a4,'original','sn','residual',True)
    
    if _telescope=='prompt':
        _values=[]
        _targets=iraf.field(img+'.sn.coo', field='', Stdout=1)
        for tt in _targets:
            _xx,_yy=string.split(tt)[0:2]
            _xx,_yy=int(float(_xx)+0.5),int(float(_yy)+0.5)
            _values.append(string.split(iraf.imstat('prompt_original.fits['+str(_xx)+','+str(_yy)+']',Stdout=1)[1])[2])
            try:
                float(_values[-1])
            except:
                _values[-1]=0
        src.delete("prompt_trimmed.fits")
        src.delete("prompt_original.fits")
        #print '### telescope Prompt: correction magnitudes !!'
        #print _values
        #print '###'

#################       Iterate Beckground    ###################################

    if _interactive:
        if not _numiter:
            answ0= raw_input(">>> Iterate on background [y/n] [y] ?")
            if not answ0: answ0='y'
        elif _numiter>=1:
            answ0='y'
        else:
            answ0='n'
    else:
        if not _numiter:
            _numiter=_numiter0
        if _numiter>=1:
            answ0='y'
        else:
            answ0='n'
    _count=0
    while answ0 == 'y':
        _count=_count+1
	print '######'
	print '###### iteration number  '+str(_count)
	print '######'
        src.delete("sn.fits,residual.fits,snfit.fits,tmp.fits")
        checkorder='yes'
        while checkorder=='yes':
	   if not _interactive or _xord and _yord:
		   if _xord and _yord: 
		       xbgord0=_xord
		       ybgord0=_yord
		   else:  
		       xbgord0=_xbgord0
		       ybgord0=_ybgord0
	   else:  	
             xbgord0 = raw_input('>>> Order of function in x for bg fit ['+str(_xbgord0)+'] ? ')
             if not xbgord0: xbgord0=_xbgord0
             else: _xbgord0=xbgord0
             ybgord0 = raw_input('>>> Order of function in x for bg fit ['+str(_ybgord0)+'] ? ')
             if not ybgord0: ybgord0=_ybgord0
             else: _ybgord0=ybgord0
           try:
                    float(xbgord0)
                    float(ybgord0)
                    checkorder='no'
           except:  
  	        print 'WARNING: value not valid !!'
	        checkorder='yes'	   

        iraf.imsurfit("skyfit","tmp",regions="all",xorder=xbgord0,yorder=ybgord0)
        midpt=float(iraf.imstat("tmp",field="mean", Stdout=1)[1])
        #iraf.imarith("original","-","tmp","sn",calctype="r",pixtype="r")
        src.pyimarith("original.fits","-","tmp.fits","sn.fits")
        src.pyimarith("sn.fits","+",midpt,"sn.fits")
        src.delete("skyfit.fits")
        
        apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy=fitsn(img,img+'.sn.coo',_recenter,fwhm0,a4,'original','sn','residual',True)
	
        print _numiter ,_count
        if _interactive:
            if not _numiter:
                answ0= raw_input(">>> Iterate on background [y/n] [y] ?")
                if not answ0: answ0='y'
            elif _count>=_numiter:
                answ0='n'
            else:
                answ0='y'
        else:
            if not _numiter:
                _numiter=_numiter0
            if _count>=_numiter:
                answ0='n'
            else:
                answ0='y'

	    
    print "***************************************************************************"
    print "#id  x_ori    y_ori      x       y      ap_ori   ap_bgsub   fit_mag   err_art err_fit"
    for i in range(len(fitmag)):
        print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(truemag[i]),"  ",str(arterr),"  ",str(magerr[i])
    print "**************************************************************************"
    if _telescope=='prompt':
        print 'Telescope Prompt:'
        print 'All magnitudes in the ec file will be corrected by these quantities:'
        for i in range(len(_values)):
            print 'SN',i, _values[i]

##########            AGGIUSTAMENTO MANUALE                     ###############
    newmag=list(array(truemag))
    if not _interactive:
       answ0='n'
    else:
       answ0= raw_input(">>> Not yet happy ? Do you want to adjust manually stellar peak ? [y/n] [n] ")
       if not answ0: answ0='n'
       elif answ0=='yes': answ0='y'
    dmag0=0
    while answ0== 'y':
        checkdm='yes'
        while checkdm=='yes':
            if len(truemag)>1: print "!!!! WARNING: all components scaled accordingly !!!!"
            _dmag0=raw_input(">>> D(mag) adjustment (positive=fainter) ["+str(dmag0)+"]")
            if _dmag0: dmag0=_dmag0
            try: 
                float(dmag0)
                checkdm='no' 
            except:
                checkdm='yes'
        apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,newmag=manusn(dmag0,apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy)
        try:
            dmag0=newmag[0]-truemag[0]
        except:
            dmag0=newmag[0]
        answ0=raw_input(">>> again ? [y/n] [y] ")
        if not answ0: answ0='y'
        elif answ0=='yes': answ0='y'


    truemag=list(array(newmag))


#
#### ESPERIMENTO DI STELLE ARTIFICIALI AL VOLO ############################
#
    if not _interactive:
        answ0='y'
    else:
        answ0= raw_input(">>> Errors estimate (through artificial star experiment ?) [y/n] [y] ")
        if not answ0: answ0='y'
        elif answ0=='yes': answ0='y'

    if answ0=='y':
        leng0=4
        _arterr2,_arterr=errore(img,'artlist.coo',size,truemag,fwhm0,leng0,_interactive,_numiter)
        if _interactive:
            arterr=raw_input("arterr ? [%6.6s] " % (str(_arterr)))
            if not arterr: arterr=_arterr
        else:
            arterr=_arterr

#    while answ0== 'y':
#        answ0= raw_input("satisfied ? [y/n] [y] ")
#        if not answ0: answ0='n'
#        elif answ0=='y': answ0='n'
#        else: answ0='n'

#######################   CHIUDI TUTTO ###################################
#
#fine:
    print "***************************************************************************" 
    print "#id  x_ori   y_ori     x     y    ap_ori ap_bgsub  fit_mag  err_art  err_fit"
    print "# id   ap_original ap_bgsub  fit_mag  err_art  err_fit"#,  >> nome0//".ec"
    print "# SN_FIT  "#, >> nome0//".ec"
    print "# id ap_ori ap-bg  fit_mag"#, >> nome0//".ec"
    for i in range(len(fitmag)):
        print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(truemag[i]),"  ",str(arterr),"  ",str(magerr[i])
    print "**************************************************************************"
    if _telescope=='prompt':
        print 'Telescope Prompt:'
        print 'All magnitude in the ec file will be corrected adding these values:'
        for ii in range(len(_values)):
            print 'SN',ii, _values[ii]
            try:
                apori3[ii]=float(apori3[ii])-float(_values[ii])
            except:
                apori3[ii]='INDEF'
            try:
                apmag3[ii]=float(apmag3[ii])-float(_values[ii])
            except:
                apmag3[ii]='INDEF'
            try:
                truemag[ii]=float(truemag[ii])-float(_values[ii])
            except:
                truemag[ii]='INDEF'

    snmesu="# id   ap_original ap_bgsub  fit_mag  err_art  err_fit\n"
    snmesu=snmesu+"# SN_FIT  \n"
    snmesu=snmesu+"# id ap_ori ap-bg  fit_mag \n"
    for i in range(len(fitmag)):
        misu='SN%1.1s %6.6s %6.6s %6.6s %6.6s %6.6s\n' % (str(i+1),str(apori3[i]),str(apmag3[i]),str(truemag[i]),str(arterr),str(magerr[i]))
        #iraf.hedit(img,'QUBASN'+str(i+1),misu[3:-1],add='yes',update='yes',verify='no')
        src.updateheader(imglong,0,'QUBASN'+str(i+1),misu[3:-1])
        snmesu=snmesu+misu
        
    src.delete("apori")
    src.delete("sec")
    src.delete("skyfit.fits")
    src.delete("sn.fits")
    src.delete("bg.fits,bgs.fits")
    src.delete("tmp*")
    src.delete(img+".sn.*")
    ff=open(img+'.ec','r')
    ecfile=ff.readlines()
    ff.close()
    src.delete(img+".ec")
    ff=open(img+'.ec','w')
    for i in ecfile:
        ff.write(i)
    for i in snmesu:
        ff.write(i)
    ff.close()
  else:
    	print '####'
  	print '#### WARNING: empty space in the list !!'
	print '####'

src.close_program(logincl)
