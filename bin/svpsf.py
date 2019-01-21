#!/usr/bin/env python

from numpy import mean
from numpy import array
from numpy import compress
from numpy import std
from numpy import average
from numpy import median
from numpy import zeros

from snoopy2 import *
import getopt
import os,sys,shutil,string,re
import snoopy2
import glob
from math import log10
try:
        import pyfits
except:
        from astropy.io import fits as pyfits
        

def reject_star(filecoordinate,_good,filemag,img,fattore):
	ff=open(filecoordinate,'r')
	vec1=ff.readlines()
	ff.close()
	alllines=[]
	for i in vec1:
		if len(string.split(i))==3:
			alllines.append(i)
			#print i

	starname=[]
	for i in alllines:
		starname.append(string.split(i)[2])
		
        _id=iraf.noao.digiphot.ptools.txdump(filemag,fields="ID",expr='yes',Stdout=1)     
	
        os.system('rm '+img+'.iter.mag')
        os.system('cp '+filemag+' '+img+'.iter.mag')
        src.delete(img+".iter.ps*,"+img+".iter.sub.fits")
        src.delete(img+".iter.grp,"+img+".iter.nst")
        iraf.noao.digiphot.daophot.psf(image=img,photfile=img+".iter.mag",pstfile=img+".iter.mag",psfimage=img+".iter.psf", mkstars='yes', opstfile=img+".iter.pst",groupfil=img+".iter.psg",verbose='yes',verify="no",interac="no")
        iraf.noao.digiphot.daophot.group(image=img,photfile=img+".iter.mag",psfimage=img+".iter.psf",groupfil=img+".iter.grp",verify='no',verbose='no')
        iraf.noao.digiphot.daophot.nstar(image=img,groupfil=img+".iter.grp",psfimage=img+".iter.psf",nstarfil=img+".iter.nst",rejfile="",verify='no',verbose='no')
        iraf.noao.digiphot.ptools.txsort(img+".iter.nst","ID")
        cutlist=iraf.noao.digiphot.ptools.txdump(img+".iter.nst",fields="XCENTER,YCENTER,ID,SHARPNESS,CHI",expr='yes',Stdout=1)
        src.delete(img+'.iter.nst')
        src.delete(img+'.iter.grp')
        src.delete(img+'.iter.psg')
        src.delete(img+'.iter.pst')
        _xx,_yy,_ID,_sharpness,_chi=[],[],[],[],[]
        for line in range(0,len(cutlist)):
            xx,yy,aa,bb,cc=string.split(cutlist[line])
            if bb!='INDEF' and cc!='INDEF':
                _xx.append(float(xx))
                _yy.append(float(yy))
                _ID.append(int(aa))
                _sharpness.append(float(bb))
                _chi.append(float(cc))
        _ID_s=compress((average(_sharpness)-std(_sharpness)*fattore<array(_sharpness))&(array(_sharpness)<average(_sharpness)+fattore*std(_sharpness)),array(_ID))
        _sharpness_s=compress((average(_sharpness)-std(_sharpness)*fattore<array(_sharpness))&(array(_sharpness)<average(_sharpness)+fattore*std(_sharpness)),array(_sharpness))
        _chi_s=compress((average(_sharpness)-std(_sharpness)*fattore<array(_sharpness))&(array(_sharpness)<average(_sharpness)+fattore*std(_sharpness)),array(_chi))
        _xx_s=compress((average(_sharpness)-std(_sharpness)*fattore<array(_sharpness))&(array(_sharpness)<average(_sharpness)+fattore*std(_sharpness)),array(_xx))
        _yy_s=compress((average(_sharpness)-std(_sharpness)*fattore<array(_sharpness))&(array(_sharpness)<average(_sharpness)+fattore*std(_sharpness)),array(_yy))        
        _ID_sc=compress((average(_chi_s)-std(_chi_s)*fattore<array(_chi_s))&(array(_chi_s)<average(_chi_s)+fattore*std(_chi_s)),array(_ID_s))
        _sharpness_sc=compress((average(_chi_s)-std(_chi_s)*fattore<array(_chi_s))&(array(_chi_s)<average(_chi_s)+fattore*std(_chi_s)),array(_sharpness_s))
        _chi_sc=compress((average(_chi_s)-std(_chi_s)*fattore<array(_chi_s))&(array(_chi_s)<average(_chi_s)+fattore*std(_chi_s)),array(_chi_s))
        _xx_sc=compress((average(_chi_s)-std(_chi_s)*fattore<array(_chi_s))&(array(_chi_s)<average(_chi_s)+fattore*std(_chi_s)),array(_xx_s))
        _yy_sc=compress((average(_chi_s)-std(_chi_s)*fattore<array(_chi_s))&(array(_chi_s)<average(_chi_s)+fattore*std(_chi_s)),array(_yy_s))

        src.delete('tmp_psf.coo')
        ff = open('tmp_psf.coo','w')
        for i in _ID_sc:
		ff.write(alllines[starname.index(_good[i-1])])
        ff.close()

        src.delete(img+'.iter.mag')
        iraf.noao.digiphot.daophot.phot(image=img, output=img+'.iter2.mag', coords='tmp_psf.coo', verify='no', interactive='no')
        _good=iraf.fields('tmp_psf.coo','3',Stdout=1)
        for i in range(len(_good)): _good[i]=re.sub(' ','',_good[i])
        return _good,img+'.iter2.mag','tmp_psf.coo'


###########################  options  #################################
bb=[]
aa=glob.glob(snoopy2.__path__[0]+'/coordinate_std/*/*list')
for i in aa:
    bb.append(i[len(snoopy2.__path__[0]+'/coordinate_std/'):-5])

help ="################################################################ \n\
help = Usage:   svpsf.py filename coordinatelist                           \n\
                input              filelist iraf format         \n\
                -l --list filename  list coordinate sequence stars    \n\
     	        [-s,--system value]     Specific photometric system  0 opt 1 inf 2 sloan   \n\
                [-p,--outph]   write output file for qubaph (use also option -o)  \n\
      	        [-o,--outputformat]     Specific out format for qubaph file fit  phot  ime   \n\
                [-i,--interactive]     interactive mode (check the stars)         \n\
		[-z,--scale]      interactive scale of image display              \n\
                [-f,--factor]     rejection factor for automatic psf              \n\
                [-a,--aperture] num   use fix aperture                            \n\
                [-F, --function] daophot functions for PSF (gauss,moffat15,moffat25,auto..) \n\
                [-v, --varord]   order of variability of PSF (default value = 0)  \n\
                available list:                                                   \n\
                "+str(bb)+"                                                       \n\
################################################################"
if len(sys.argv)<=3:
    print help
    sys.exit()

    
logincl=src.open_program()    

from pyraf import iraf
iraf.set(clobber='yes')
iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)
iraf.tv.rimexam.backgrou = 'yes'
iraf.tv.rimexam.magzero = 0
#iraf.set(stdimage='imt2048')

interactive = False   # default interative mode 
out = ''
system=''
scale= False
_telescope=''
subdirectory=['optical/','infrared/','sloan/']
fattore=3.0
_outputformat='fit'
outph=False
cutstars=False
_aperture=''

varord=0
_function='gauss'
imglist=src.readlist(sys.argv[1])

coordinatelist=''

options, args = getopt.getopt(sys.argv[2:],'l:,i,s:,z,f:,o:,p,a:,F:,v:',['coordinatelist=','interactive','system=','scale','fattore','outputformat=','output','aperture=','function=','varord='])
for opt,arg in options:
    if opt in ('-l', '--list'): coordinatelist = arg
    if opt in ('-i', '--interactive'): interactive = True
    if opt in ('-p', '--outph'): outph = True
    if opt in ('-s', '--system'): system = int(arg) 
    if opt in ('-z', '--scale'): scale = True
    if opt in ('-f', '--fattore'): fattore = float(arg)
    if opt in ('-o', '--outputformat'): _outputformat = str(arg)
    if opt in ('-a', '--aperture'): _aperture = arg
    if opt in ('-F', '--function'): _function = arg
    if opt in ('-v', '--varord'): varord = arg

####################################### check ##########################
if _function not in ['gauss','moffat15','moffat25','lorentz','penny1','penny2','auto']:
	print "error: daophot function ["+str(_function)+"] doesn'exist !"
	sys.exit()

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
        print '### if you want to continue anyway,'
        print '### run "svother.py" to correct headers '
        print '####################### '
#        sys.exit()
##########################################################################


mag,airmass,exptime,fwhm,magphot,magime={},{},{},{},{},{}
airmass['U'],airmass['B'],airmass['V'],airmass['R'],airmass['I']=1,1,1,1,1
exptime['U'],exptime['B'],exptime['V'],exptime['R'],exptime['I']=1,1,1,1,1
airmass['u'],airmass['g'],airmass['r'],airmass['i'],airmass['z']=1,1,1,1,1
exptime['u'],exptime['g'],exptime['r'],exptime['i'],exptime['z']=1,1,1,1,1
airmass['J'],airmass['H'],airmass['K']=1,1,1
exptime['J'],exptime['H'],exptime['K']=1,1,1



######################## SET   IRAF   PARAMETERS  #######################
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)

from iraf import digiphot
from iraf import daophot
from iraf import ptools

if _aperture:
    iraf.tv.rimexam.iterati = 1
    iraf.tv.rimexam.radius = _aperture
else:
    iraf.tv.rimexam.iterati = 3
    iraf.tv.rimexam.radius = 7


zmag = 0 
iraf.digiphot.daophot.photpars.zmag = zmag
iraf.digiphot.daophot.photpars.zmag = zmag
iraf.digiphot.daophot.daopars.function = _function
iraf.digiphot.daophot.daopars.varord = varord


system_used=[]
system_used=[system]
############ START LOOP ON LIST IMAGES  ###################################
for imglong in imglist:
  if imglong:
    if imglong[-5:]=='.fits': 
	    img=imglong[:-5]
    else:
	    img=imglong
#########################  READ HEADER  ########
    _header=src.read_parameter(_telescope,system)
    _gain=src.gain(imglong,_header,_telescope)
    _ron=src.ron(imglong,_header,_telescope)
    _datamin=src.ccdmin(imglong,_header,_telescope)
    _datamax=src.ccdmax(imglong,_header,_telescope)
    #_datamax=80000
    iraf.noao.digiphot.daophot.datapars.readnoi = _ron
    iraf.noao.digiphot.daophot.datapars.epadu = _gain
    iraf.digiphot.daophot.datapars.datamin = _datamin
    iraf.digiphot.daophot.datapars.datamax = _datamax
    instrument=src.instrument(imglong,_header,_telescope) 
    _airmass=src.airmass(imglong,_header,_telescope)
    _exptime=src.exptime(imglong,_header,_telescope)
    _date=src.date(imglong,_header,_telescope)
    _xdimen=src.xdimen(imglong,_telescope)
    _ydimen=src.ydimen(imglong,_telescope)
    _object=src.objects(imglong,_header,_telescope)
    _filter=src.filter(imglong,_header,_telescope)
    print _gain,_ron,_datamax,_datamin
    iraf.digiphot.daophot.datapars.exposure = _header['hed_exptime']
    iraf.digiphot.daophot.datapars.airmass = _header['hed_airmass']
#    iraf.digiphot.daophot.datapars.filter = _header['hed_filter1']
#    iraf.digiphot.daophot.datapars.exposure = _exptime
#    iraf.digiphot.daophot.datapars.airmass = _airmass
    iraf.digiphot.daophot.datapars.filter = _filter
    
######################## filter #################################
    filter=src.filtername(_telescope,_filter,system)
    if filter!='unknown':
	    system0=system
	    if system not in system_used: 
		    system_used.append(system)
    else:
	    filter=src.filtername(_telescope,_filter,0)
	    if filter!='unknown':
		    if 0 not in system_used: system_used.append(0)
		    system0=0
	    if filter=='unknown':
		    filter=src.filtername(_telescope,_filter,1)
		    if filter!='unknown':
			    if 1 not in system_used: system_used.append(1)
			    system0=1
            if filter=='unknown':
                    filter=src.filtername(_telescope,_filter,2)
		    if filter!='unknown':
			    if 2 not in system_used: system_used.append(2)
			    system0=2
            if filter=='unknown':
                    print '### Warning: filter not recognised in this photometric system '
                    filter=_filter
                    system0=system
                    
    print '##########  '+str(filter)+'  ##############'
    airmass[filter]=_airmass
    exptime[filter]=_exptime
    
#############################################################
    ###         cordinate from system     #####
    src.delete(coordinatelist+".tv")
    src.delete("_templ*")
    _z1,_z2='',''
    try:
        dir_system=subdirectory[system0]
        iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','2,3,1', Stdout=coordinatelist+'.tv')
        iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1', Stdout='_templ.coo')
        iraf.wcsctran(coordinatelist+'.tv','_templ.coo2','home$coordinate_std/'+dir_system+coordinatelist+'_templ.fits',inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
        a=iraf.fields('_templ.coo','1', Stdout=1)
        b=iraf.fields('_templ.coo2','1,2', Stdout=1)[2:]
        ff = open(coordinatelist+'_templ.coo','w')
        for i in range(len(a)):
            ff.write(b[i]+'\t'+a[i]+' \n')
        ff.close()
        standard=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1', Stdout=1)

        if system0==0 or system0==2 : 
            standard2=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6,7,8',Stdout=1)
        else:
            standard2=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6',Stdout=1)
    except:
        print "WARNING: no coordinate file found in"+snoopy2.__path__[0]+'/coordinate_std/'+dir_system+'  !!! '
        src.close_program(logincl)
        
    stars=[]
    for i in standard:
        stars.append(re.sub(' ','',i))
######################################################################
##############       ASTROMETRY      ################

    src.delete("tmp*")
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    print '######### Select FRAME TILE on your DS9 !!!!!'
    _z10,_z20,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+dir_system+coordinatelist+'_templ.fits',2,'','',False)
    if not goon: src.close_program(logincl)
    iraf.tvmark(2,coordinatelist+'_templ.coo',mark="circle",number='no',label='yes',radii=10,nxoffse=5,nyoffse=5,color=214,txsize=4)
    _z1,_z2,goon=src.display_image(img+'.fits',1,_z1,_z2,scale)
    if not goon: src.close_program(logincl)
    iraf.tvmark(1,'tmp.'+filter+'.coo',mark="circle",number='no',label='yes',radii=10,nxoffse=5,nyoffse=5,color=205,txsize=4)

    if interactive:
        answ = raw_input('is the astrometry of the field good [y/n] ? [y] ')
    else:
        answ = 'y'
    if not answ: answ='y'
    elif answ=='no':answ='n'

    if answ=='n':
        try:
            src.delete('tmp.'+filter+'.coo')
            src.delete('tmp.ccdb')
            iraf.ccmap('_first.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)
            iraf.wcsctran('_first_image.tv','tmp.'+filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
            iraf.tvmark(1,'tmp.'+filter+'.coo',mark="circle",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=206,txsize=4)
            answ = raw_input('AND NOW, is the astrometry of the field good [y/n] ? [y] ')
            if not answ: answ='y'
            elif answ=='no':answ='n'
        except: pass
            
    while answ=='n':
            _z1,_z2,goon=src.display_image(img+'.fits',1,_z1,_z2,False)
            if not goon: src.close_program(logincl)

            print '>> Identify (minimum 2, preferably 3) reference stars (mark with "a")' 
            iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
            iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=214,txsize=4)
            xycoo = iraf.fields('tmp.coo','1,2,13',Stdout=1)
            print '>> Identify reference stars'
            idcat = []
            for i in range(len(xycoo)):
                idcat.append(raw_input('Star '+str(i+1)+'= ? ')) 
            ff = open(coordinatelist+'.tv','r')
            rr = ff.readlines()
            ff.close()
            gg = open('tmp.ccmap','w')
            fw = []
            for i in range(len(idcat)):
                idpos=stars.index(idcat[i])
                _rr = string.split(rr[idpos])
                _x,_y,_fw = string.split(xycoo[i])
                gg.write(_x+' '+_y+' '+_rr[0]+' '+_rr[1]+' \n')
                fw.append(float(_fw))
            gg.close()
            src.delete('_first_image.tv')
            src.delete('_first.ccmap')
            os.system('cp '+coordinatelist+'.tv _first_image.tv')
            os.system('cp tmp.ccmap  _first.ccmap')
            iraf.ccmap('tmp.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)
            _z1,_z2,goon=src.display_image(img+'.fits',1,_z1,_z2,False)
            if not goon: src.close_program(logincl)
            
            src.delete('tmp.'+filter+'.coo')
            iraf.wcsctran(coordinatelist+'.tv','tmp.'+filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
            iraf.tvmark(1,'tmp.'+filter+'.coo',mark="circle",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=205,txsize=4)
            src.delete("tmp.star")
            src.delete("tmp.ccmap")
            src.delete("tmp.coo")
            src.delete("tmp.ccdb")
            answ = raw_input('is the astrometry of the field good [y/n]?  [y] ')
            if not answ: answ='y'
            elif answ=='no':answ='n'

    src.delete("tmp.star")
    iraf.ccfind('home$coordinate_std/'+dir_system+coordinatelist+'.list','tmp.star',img,lngcolu=2,latcolu=3,lngunit='degrees',usewcs='yes')
    iraf.ccmap('tmp.star','tmp.ccdb',images=img,fitgeome='rscale',xcolum=9, ycolum=10,lngcolum=2,latcolumn=3,lngunit='degrees',update='yes',interact=False)
    src.delete('tmp.'+filter+'.coo')
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    src.delete("tmp.ccdb")
    src.delete("tmp.star")

#################   FWHM COMPUTATION      #########
            
    src.delete("tmp.imex_output")
    stars_pos=[]
    alllines=[]

    ff = open('tmp.'+filter+'.coo','r')
    for j in range(0,3): exl=ff.readline()
    for j in range(len(stars)): 
        alllines.append(ff.readline())
        xx=string.split(alllines[-1])[0:2]
	os.system('echo '+xx[0]+' '+xx[1]+'> tmp.two')
        pp=iraf.imexam(input=img, frame=1, logfile='', keeplog='no', imagecur='tmp.two', defkey = 'm', wcs='logical', use_disp='no', Stdout=1)
        if not 1.<= float(xx[0])<=float(_xdimen) or not 1.<= float(xx[1]) <=float(_ydimen):
            stars_pos.append(0)
        elif string.split(pp[1])[-1]==string.split(pp[1])[-2]:
            stars_pos.append(0)
        else:
            stars_pos.append(1)
    ff.close()
    ff = open('tmp.'+filter+'good.coo','w')
    for i in range(len(stars_pos)):
        if stars_pos[i]:
	     ff.write(alllines[i])
    ff.close()
    for i in range(len(alllines)):
        if stars_pos[i]:
            ff = open('tmp.one','w')
            xx=string.split(alllines[i])[0:2]
            ff.write(xx[0]+' '+xx[1]+'  a')
            ff.close()
            try:
                iraf.imexam(input=img, frame=1, logfile='tmp.imex_output', keeplog='yes', imagecur='tmp.one', wcs='logical', use_disp='no')
            except:
                if not os.path.isfile('tmp.imex_output'):
                    os.system("echo '# [1] "+str(img)+" - "+str(coordinatelist)+"' > tmp.imex_output")
                    os.system("echo '#   COL    LINE   COORDINATES      R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT' >> tmp.imex_output")
                    os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")
                else:
                    os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")

    _fwhm0 = iraf.fields('tmp.imex_output','15',Stdout=1)
    _mag = iraf.fields('tmp.imex_output','6',Stdout=1)
    _fwhm,_magime=[],[]
    j=0
    for i in range(len(stars)):
        if stars_pos[i]:
            try: _magime.append(float(_mag[j])+2.5*log10(_exptime))
            except: _magime.append(float(9999))
            try: _fwhm.append(float(_fwhm0[j]))
            except: _fwhm.append(float(9999))
            j=j+1
        else: 
            _magime.append(float(9999))
            _fwhm.append(float(9999))

    fwhm_ave=compress((array(_fwhm)<999),_fwhm)
    
    fwhm_ave1=compress((average(fwhm_ave)-2*std(fwhm_ave)<array(fwhm_ave))&(array(fwhm_ave)<average(fwhm_ave)+std(fwhm_ave)*2),array(fwhm_ave))
    _fwhm_ave=mean(compress((average(fwhm_ave1)-2*std(fwhm_ave1)<array(fwhm_ave1))&(array(fwhm_ave1)<average(fwhm_ave1)+std(fwhm_ave1)*2),array(fwhm_ave1)))

    fwhm_ave22=compress((average(fwhm_ave)-2*std(fwhm_ave)<array(fwhm_ave))&(array(fwhm_ave)<average(fwhm_ave)+std(fwhm_ave)*2),array(fwhm_ave))
    _fwhm_ave2=median(compress((average(fwhm_ave22)-2*std(fwhm_ave22)<array(fwhm_ave22))&(array(fwhm_ave22)<average(fwhm_ave22)+std(fwhm_ave22)*2),array(fwhm_ave22)))
    
    checkfwhm='yes'
    while checkfwhm=='yes':
        print '################ FWHM(median) = '+str(_fwhm_ave2)+'  '
        print '################ FWHM(mean) = '+str(_fwhm_ave)+'  '
        if interactive:
            fwhm_ave = raw_input('################ FWHM = ['+str(_fwhm_ave2)+'] ? ')
        else:
            fwhm_ave = str(_fwhm_ave2)
        try:
            if not fwhm_ave: fwhm_ave=_fwhm_ave2
            else: fwhm_ave=float(fwhm_ave)
            checkfwhm='no'
        except:
            print 'WARNING: FWHM not good !!!!'
            checkfwhm='yes'

    fwhm[filter]=fwhm_ave
    src.updateheader(imglong,0,'qubafwhm',fwhm_ave)

###############################################################################################
#################        magnitude cut               ###################
    if len(compress(array(_magime)<-5,array(_magime)))>=10:
	    print 'More than 10 stars are brighter than -5 '
	    _magime2=compress(array(_magime)<-5,array(_magime))
	    if len(_magime2)!=len(compress(array(_magime)<=99,array(_magime))):
			    
		    print 'Cut all the others'
		    src.delete("tmp."+filter+'good_cut.coo')	    
		    ff = open('tmp.'+filter+'good_cut.coo','w')	    
		    for i in range(len(_magime)):
			    if _magime[i] < -5:
				    ff.write(alllines[i])
		    ff.close()
		    cutstars=True
	    else:
		    print 'All stars brighter than -5 '
	    
    elif len(compress(array(_magime)<-4,array(_magime)))>=10:
	    print 'More than 10 stars brighter than -4 '
	    print 'cut all the others'
	    _magime2=compress(array(_magime)<-4,array(_magime))
	    if len(_magime2)!=len(compress(array(_magime)<=99,array(_magime))):
		    print 'cut all the others'
		    src.delete("tmp."+filter+'good_cut.coo')	    
		    ff = open('tmp.'+filter+'good_cut.coo','w')
		    for i in range(len(_magime)):
			    if _magime[i] < -4:
				    ff.write(alllines[i])
		    ff.close()
		    cutstars=True
	    else:
		    print 'All stars brighter than -4 '
		    
    elif len(compress(array(_magime)<-3,array(_magime)))>=10:
	    print 'More than 10 stars brighter than -3 '
	    print 'cut all the others'
	    _magime2=compress(array(_magime)<-3,array(_magime))
	    if len(_magime2)!=len(compress(array(_magime)<=99,array(_magime))):
		    src.delete("tmp."+filter+'good_cut.coo')	    
		    ff = open('tmp.'+filter+'good_cut.coo','w')
		    for i in range(len(_magime)):
			    if _magime[i] < -3:
				    ff.write(alllines[i])
		    ff.close()
		    cutstars=True
	    else:
		    print 'All stars brighter than -3 '
    else:
	    print 'No cut in magnitude. Not enough stars.'
	    cutstars=False
	    if len(compress(array(_magime)<=99,array(_magime)))>=10:
		    print 'WARNING: I will use stars fainter that -3 (instrumental magnitude) '
		
##########################################################################################



##############  DAOPHOT MAGNITUDE     ####################
################ AND APERTURE CORRECTION  ################# 

    _z1,_z2,goon=src.display_image(img+'.fits',1,_z1,_z2,False)
    if not goon: src.close_program(logincl)
    

    if _aperture:
	    annulus=float(_aperture)+float(fwhm_ave)
	    iraf.noao.digiphot.daophot.fitskypars.annulus = annulus
	    iraf.noao.digiphot.daophot.photpars.apertures = float(_aperture)
    else:
	    annulus=4*fwhm_ave
	    iraf.noao.digiphot.daophot.fitskypars.annulus = annulus
	    iraf.noao.digiphot.daophot.photpars.apertures = fwhm_ave*3.

    src.delete(img+".mag")
    iraf.noao.digiphot.daophot.phot(image=img, output=img+'.mag', coords='tmp.'+filter+'good.coo', verify='no', interactive='no')
    _magphot2=iraf.noao.digiphot.ptools.txdump(img+".mag",fields="mag",expr='yes',Stdout=1)
    j=0
    _magphot=[]
    for i in range(len(stars_pos)):
        if stars_pos[i]:
	        l=i+j
                try: 
                    _magphot.append(float(_magphot2[l]))
                except: _magphot.append(9999.)
	else:
	    _magphot.append(9999.)
	    j=j-1
	    
##############################    prompt correction
    if _telescope=='prompt':
        try:
                prompttele = pyfits.getheader(imglong)['OBSERVAT']
	except:
                prompttele=''
        if string.count(prompttele,'2'): _prompt='prompt2'
        elif string.count(prompttele,'3'): _prompt='prompt3'
        elif string.count(prompttele,'4'): _prompt='prompt4'
        elif string.count(prompttele,'5'): _prompt='prompt5'
        else: print 'WARNING: telescope not found !! '
        _trimsection=src.trimsec(imglong,_header,_telescope)
        if _trimsection:
            xx=string.split(string.split(_trimsection[1:-1],',')[0],':')[0]
            yy=string.split(string.split(_trimsection[1:-1],',')[1],':')[0]
        else:
            xx,yy=0,0
        _coordinate_star=iraf.noao.digiphot.ptools.txdump(img+".mag",fields="XCENTER,YCENTER",expr='yes',Stdout=1)
        _coord_origi_x,_coord_origi_y,_values=[],[],[]
	try:
		_z10,_z20,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+_telescope+'/'+_prompt+'.fits',2,'','',False)
		if not goon: src.close_program(logincl)
	except: pass
        j=0
        for ii in range(len(stars_pos)):
            if stars_pos[ii]:
                l=ii+j
                xxn,yyn=string.split(_coordinate_star[l])
                _coord_origi_x.append(int(float(xx)+float(xxn)))
                _coord_origi_y.append(int(float(yy)+float(yyn)))
                _xx=int(float(xx)+float(xxn))
                _yy=int(float(yy)+float(yyn))
		try:
			_values.append(float(string.split(iraf.imstat('home$coordinate_std/'+_telescope+'/'+_prompt+'.fits['+str(_xx)+','+str(_yy)+']',Stdout=1)[1])[2]))
		except:
			_values.append(0)
            else:
                _coord_origi_x.append(0)
                _coord_origi_y.append(0)
                _xx=int(0)
                _yy=int(0)
                _values.append(0)
                j=j-1

################################################
    iraf.noao.digiphot.daophot.daopars.psfrad = annulus
    iraf.noao.digiphot.daophot.daopars.fitrad = fwhm_ave
    iraf.noao.digiphot.daophot.daopars.fitsky = 'yes'
    iraf.noao.digiphot.daophot.daopars.sannulus = annulus
    iraf.noao.digiphot.daophot.daopars.recenter = 'yes'
    iraf.noao.digiphot.daophot.daopars.varorder = varord
    
    iraf.tvmark(1,'tmp.'+filter+'good.coo',mark="circle",number='yes',label='no',radii=10,nxoffse=5,nyoffse=5,color=203,txsize=4)

    print "-------------------------------------------------------------------------------"
    _id=iraf.noao.digiphot.ptools.txdump(img+".mag",fields="ID",expr='yes',Stdout=1)
    _good=iraf.fields('tmp.'+filter+'good.coo','3',Stdout=1)
    for i in range(len(_good)): _good[i]=re.sub(' ','',_good[i])

#################################################
    #  _id   is the identification number of stars
    # _good  is the string identification of stars in the list
    for i in range(len(_id)):
         print _id[i]+' = '+_good[i]
    
    print "-------------------------------------------------------------------------------"
    print " on GRAPHIC display  > a <  to accept star,  > d <  to delete  "
    print " on IMAGE   display  > f <  to fit PSF, > w <  to write PSF fit,  > q <  to quit"
    print "-------------------------------------------------------------------------------"

    src.delete(img+".ps*,"+img+".sub.fits")
    src.delete(img+".grp,"+img+".nst")

    if interactive:
        iraf.noao.digiphot.daophot.psf(image=img,photfile=img+".mag",pstfile=img+".mag",psfimage=img+".psf", mkstars='yes', opstfile=img+".pst",groupfil=img+".psg",verbose='yes',verify="no",interac="yes")
        apnum=iraf.noao.digiphot.ptools.txdump(img+".pst",fields="ID",expr='yes',Stdout=1)      
        iraf.noao.digiphot.daophot.group(image=img,photfile=img+".mag",psfimage=img+".psf",groupfil=img+".grp",verify='no',verbose='no')
        iraf.noao.digiphot.daophot.nstar(image=img,groupfil=img+".grp",psfimage=img+".psf",nstarfil=img+".nst",rejfile="",verify='no',verbose='no')
        iraf.noao.digiphot.daophot.substar(image=img,photfile=img+".nst",exfile="",psfimage=img+".psf",subimage=img+".sub",verify='no',verbose='no')

    else:

 ###########################################################################################################################################
  ####  choose star with one interaction
##
        if cutstars:
		src.delete(img+"cut.mag")
		iraf.noao.digiphot.daophot.phot(image=img, output=img+'cut.mag', coords='tmp.'+filter+'good_cut.coo', verify='no', interactive='no')
		ff=open('tmp.'+filter+'good_cut.coo','r')
		vec1=ff.readlines()
		ff.close()
		alllinescut=[]
		for i in vec1:
			if len(string.split(i))==3:
				alllinescut.append(i)
		starnamecut=[]
		for i in alllinescut:
			#print i
			starnamecut.append(string.split(i)[2])
		idcut=[]
		for ii in range(len(alllinescut)):
			idcut.append(str(_good.index(starnamecut[ii])+1))
		imputgood=starnamecut
		filemag=img+'cut.mag'
	else:
		imputgood=_good
		filemag=img+'.mag'
#########################################################################################

        filecoo='tmp.'+filter+'.coo'
	_z1,_z2,goon=src.display_image(img+'.fits',1,_z1,_z2,False)
	iraf.tvmark(1,'tmp.'+filter+'.coo',mark="circle",number='no',label='yes',radii=20,nxoffse=15,nyoffse=15,color=205,txsize=4)
        _good2,filemag,filecoo=reject_star(filecoo,imputgood,filemag,img,fattore)
	
	iraf.tvmark(1,filecoo,mark="circle",number='no',label='yes',radii=22,nxoffse=30,nyoffse=30,color=204,txsize=4)
        _good3,filemag2,filecoo=reject_star(filecoo,_good2,filemag,img,fattore)
	
	iraf.tvmark(1,filecoo,mark="circle",number='no',label='yes',radii=24,nxoffse=45,nyoffse=45,color=203,txsize=4)
	
        src.delete(img+".ps*,"+img+".sub.fits")
        src.delete(img+".grp,"+img+".nst")
########
        iraf.noao.digiphot.daophot.psf(image=img,photfile=filemag2,pstfile=filemag2,psfimage=img+".psf", mkstars='yes', opstfile=img+".pst",groupfil=img+".psg",verbose='yes',verify="no",interac="no")
        src.delete(img+'.nst')
        src.delete(img+'.grp')

        apnum=iraf.noao.digiphot.ptools.txdump(img+".pst",fields="ID",expr='yes',Stdout=1)
        for i in range(0,len(apnum)):
            apnum[i]=str(_good.index(_good3[i])+1)

#############################################################################################################################################            
        iraf.noao.digiphot.daophot.group(image=img,photfile=img+".mag",psfimage=img+".psf",groupfil=img+".grp",verify='no',verbose='no')
        iraf.noao.digiphot.daophot.nstar(image=img,groupfil=img+".grp",psfimage=img+".psf",nstarfil=img+".nst",rejfile="",verify='no',verbose='no')
        iraf.noao.digiphot.daophot.substar(image=img,photfile=img+".nst",exfile="",psfimage=img+".psf",subimage=img+".sub",verify='no',verbose='no')


    
    _z1,_z2,goon=src.display_image(img+".sub.fits",1,_z1,_z2,False)
    if not goon: src.close_program(logincl)
    
    iraf.noao.digiphot.ptools.txsort(img+".nst","ID")
    dddd=iraf.noao.digiphot.ptools.txdump(img+".nst",fields="ID,mag,merr",expr='yes',Stdout=1)
    dddd2=iraf.noao.digiphot.ptools.txdump(img+".nst",fields="ID",expr='yes',Stdout=1)

    j=0
    fitmag,fitmagerr=[],[]
    for i in range(len(stars)):
        if stars_pos[i]:
           try:
              indice=dddd2.index(str(i-j+1))     #   number of star in img.nst
              fitmag.append(float(string.split(dddd[indice])[1]))
              fitmagerr.append(float(string.split(dddd[indice])[2]))
           except: 
              fitmag.append(float(9999))
              fitmagerr.append(float(0.0))
	else:
	   j=j+1
	   fitmag.append(float(9999))
	   fitmagerr.append(float(0.0))


    apco,apco2=[],[]
    for i in range(len(stars)):
        if stars_pos[i]:
		try: 
			_id[_good.index(stars[i])] 
			if _id[_good.index(stars[i])] in apnum:
				if _magphot[i]>=99 or fitmag[i]>= 99:
					apco.append(9999)
				else:
					apco.append(_magphot[i]-fitmag[i])
				if _magime[i]>=99 or fitmag[i]>= 99:
					apco2.append(9999)
				else:
					apco2.append(_magime[i]-fitmag[i])
			else:
				apco.append(9999)
				apco2.append(9999)
		except:
			apco.append(9999)
			apco2.append(9999)

        else:
            apco.append(9999)
            apco2.append(9999)
    if _telescope=='prompt' and _values:
        print '#####  Telescope = prompt'
        print '#####  Applied correction  !! '
        for i in range(len(stars)):
            if _magime[i]<=99:
                _magime[i]=_magime[i]-_values[i]
            if _magphot[i]<=99:
                _magphot[i]=_magphot[i]-_values[i]
            if fitmag[i]<=99:
                fitmag[i]=fitmag[i]-_values[i]

    magime[filter]=_magime
    magphot[filter]=_magphot
    mag[filter]=fitmag
    
    print "*****************************************************"
    print "id      ime_mag   ap_mag   ph_mag  fit_mag   ap_cor  ph_cor"   
    for i in range(len(stars)):
        print '\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s' % (str(stars[i]),str(_magime[i]),str(_magphot[i]),str(fitmag[i]),str(apco2[i]),str(apco[i]))

    #raw_input('go on ?')
    #aper_cor=compress((array(fitmag)<999)&(array(_magphot)<999),array(_magphot)-array(fitmag))
    #aper_cor2=compress((array(fitmag)<999)&(array(_magime)<999),array(_magime)-array(fitmag))
    aper_cor=compress((array(apco)<999),array(apco))
    aper_cor2=compress((array(apco2)<999),array(apco2))
    print"*****************************************************"
    try:
	    print 'APERTURE CORRECTION --> imex   (ap_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(mean(aper_cor2)),str(std(aper_cor2)))
	    print 'APERTURE CORRECTION --> phot   (ph_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(mean(aper_cor)),str(std(aper_cor)))
	    print 'number of stars = '+str(len(aper_cor2))
	    print '*********************************************************************************'
    except:
	    print 'APERTURE CORRECTION --> imex   (ap_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(aper_cor2),0)
	    print 'APERTURE CORRECTION --> phot   (ph_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(aper_cor),0)
	    print 'number of stars = '+str(len(aper_cor2))
	    print '*********************************************************************************'
    try:
	    aper_cor22=compress((mean(aper_cor2)-std(aper_cor2)*3<array(aper_cor2))&(array(aper_cor2)<mean(aper_cor2)+std(aper_cor2)*3),array(aper_cor2))
	    _aper_cor2=compress((mean(aper_cor22)-std(aper_cor22)*3<array(aper_cor22))&(array(aper_cor22)<mean(aper_cor22)+std(aper_cor22)*3),array(aper_cor22))
	    aper_cor1=compress((mean(aper_cor)-std(aper_cor)<array(aper_cor))&(array(aper_cor)<mean(aper_cor)+std(aper_cor)),array(aper_cor))
	    _aper_cor=compress((mean(aper_cor1)-std(aper_cor1)*3<array(aper_cor1))&(array(aper_cor1)<mean(aper_cor1)+std(aper_cor1)*3),array(aper_cor1))
	    print 'APERTURE CORRECTION --> imex   (ap_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(mean(_aper_cor2)),str(std(_aper_cor2)))
	    print 'APERTURE CORRECTION --> phot   (ph_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(mean(_aper_cor)),str(std(_aper_cor)))
	    print 'number of stars= '+str(len(_aper_cor2))+'       (after 3 sigma rejection 2 times) '
    except:
	    aper_cor22=aper_cor2
	    _aper_cor2=aper_cor22
	    aper_cor1=aper_cor
	    _aper_cor=aper_cor1
	    print 'APERTURE CORRECTION --> imex   (ap_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(_aper_cor2),0)
	    print 'APERTURE CORRECTION --> phot   (ph_mag) = %6.6s\t+/-\t%6.6s \n ' % (str(_aper_cor),0)
	    print 'number of stars= '+str(len(_aper_cor2))+'       (after 3 sigma rejection 2 times) '

    checkapco='yes'
    while checkapco=='yes':
        if interactive :
            apcofin = raw_input('Aperture correction ?  [%6.6s]' % (str(mean(_aper_cor2))))
        else:
            apcofin=mean(_aper_cor2)
        if not apcofin: apcofin=mean(_aper_cor2)
        try: 
            float(apcofin)
            checkapco='no'
        except:
            print 'WARNING: aperture correction not valid !!!'
            checkapco='yes'
############  WRITE EC FILE ##############
#    mag[filter]=_magime


    fil = open(img+'.ec','w')
    fil.write('#file = '+str(img)+'   title = '+str(_object)+'    filter='+str(filter)+'  expt=%6.6s \n' % (str(_exptime)))
    fil.write('#airmass= %6.6s \n'% (str(_airmass)))
    fil.write('# PSF \n')
    fil.write('#    Reference FWHM= %6.6s    aperture correction =%6.6s \n' % (str(fwhm_ave),str(apcofin)))
    fil.write('# id fwhm  ph_mag  fit_mag  fiterr \n')
    for i in range(len(stars)):
        fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(stars[i]),str(_fwhm[i]),str(_magphot[i]),str(fitmag[i]),str(fitmagerr[i])))
    fil.close()
    src.updateheader(img+'.psf.fits',0,'QUBAFIL',filter)
    src.updateheader(img+'.fits',0,'QUBAAPCO',apcofin)
    src.delete(img+'.nst')
    src.delete(img+'.grp')
    src.delete(img+'.psg')
    src.delete(img+'.pst')
    src.delete(img+'.sub.fits')
    src.delete('tmp.'+filter+'.coo')
    
#    if interactive:
#        fil=open('aperture_check_'+str(filter)+'_'+str(fattore)+'_manuale.asci','w')
#    else:
#        fil=open('aperture_check_'+str(filter)+'_'+str(fattore)+'_automatico.asci','w')
#    fil.write('#\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\n' % ('stars','_magime','_magphot','fitmag','apco2','apco'))
#    for i in range(len(stars)):
#        fil.write('\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\n' % (str(stars[i]),str(_magime[i]),str(_magphot[i]),str(fitmag[i]),str(apco2[i]),str(apco[i])))
#    fil.close()

  else:
      	print '####'
  	print '#### WARNING: empty space in the list !!'
	print '####'


src.delete(coordinatelist+'_templ.coo')
src.delete(coordinatelist+'.tv')
src.delete('tmp*')
src.delete('tmp.ccdb')
src.delete('_first_image.tv')
src.delete('_first.ccmap')
src.delete('_templ.co*')
src.delete('tmp.imex_output')
src.delete('tmp.one')
src.delete('airmass.txt')
src.delete('*iter*.mag')
src.delete('*iter*.psf.fits')

#########################################
#  output for qubaph
#############################################

outputformat={}
outputformat['fit']=mag
outputformat['ime']=magime
outputformat['phot']=magphot

if _outputformat not in ['fit','phot','ime']:
    _outputformat='fit'


standard_value={}
standard_value[system]=standard2

if len(system_used)>=2:
    for ii in system_used:
        if ii!=system:
            if ii==0 or ii==2 : 
                othervalue=iraf.fields('home$coordinate_std/'+subdirectory[ii]+coordinatelist+'.list','1,4,5,6,7,8',Stdout=1)
            else:
                othervalue=iraf.fields('home$coordinate_std/'+subdirectory[ii]+coordinatelist+'.list','1,4,5,6',Stdout=1)
            standard_value[ii]=othervalue


for ii in system_used:
    if ii==0:    
        try:
            outputformat[_outputformat]['U']
            summaryU='band U -> ok\n'
        except:
            outputformat[_outputformat]['U']=list(zeros(len(standard_value[ii]))+9999) 
            summaryU='band U -> none\n'
        try:
            outputformat[_outputformat]['B']
            summaryB='band B -> ok\n'
        except:
            outputformat[_outputformat]['B']=list(zeros(len(standard_value[ii]))+9999) 
            summaryB='band B -> none\n'
        try:
            outputformat[_outputformat]['V']
            summaryV='band V -> ok\n'
        except:
            outputformat[_outputformat]['V']=list(zeros(len(standard_value[ii]))+9999) 
            summaryV='band V -> none\n'
        try:
            outputformat[_outputformat]['R']
            summaryR='band R -> ok\n'
        except:
            outputformat[_outputformat]['R']=list(zeros(len(standard_value[ii]))+9999) 
            summaryR='band R -> none\n'
        try:
            outputformat[_outputformat]['I']
            summaryI='band I -> ok\n'
        except:
            outputformat[_outputformat]['I']=list(zeros(len(standard_value[ii]))+9999) 
            summaryI='band I -> none\n'
    elif ii==2:
        try:
            outputformat[_outputformat]['u']
            summaryU='band u -> ok\n'
        except:
            outputformat[_outputformat]['u']=list(zeros(len(standard_value[ii]))+9999) 
            summaryU='band u -> none\n'
        try:
            outputformat[_outputformat]['g']
            summaryB='band g -> ok\n'
        except:
            outputformat[_outputformat]['g']=list(zeros(len(standard_value[ii]))+9999) 
            summaryB='band g -> none\n'
        try:
            outputformat[_outputformat]['r']
            summaryV='band r-> ok\n'
        except:
            outputformat[_outputformat]['r']=list(zeros(len(standard_value[ii]))+9999) 
            summaryV='band r -> none\n'
        try:
            outputformat[_outputformat]['i']
            summaryR='band i -> ok\n'
        except:
            outputformat[_outputformat]['i']=list(zeros(len(standard_value[ii]))+9999) 
            summaryR='band i -> none\n'
        try:
            outputformat[_outputformat]['z']
            summaryI='band z -> ok\n'
        except:
            outputformat[_outputformat]['z']=list(zeros(len(standard_value[ii]))+9999) 
            summaryI='band z -> none\n'
    elif ii==1:
        try:
            outputformat[_outputformat]['J']
            summaryU='band J -> ok\n'
        except:
            outputformat[_outputformat]['J']=list(zeros(len(standard_value[ii]))+9999) 
            summaryU='band J -> none\n'
        try:
            outputformat[_outputformat]['H']
            summaryB='band H -> ok\n'
        except:
            outputformat[_outputformat]['H']=list(zeros(len(standard_value[ii]))+9999) 
            summaryB='band H -> none\n'
        try:
            outputformat[_outputformat]['K']
            summaryV='band K-> ok\n'
        except:
            outputformat[_outputformat]['K']=list(zeros(len(standard_value[ii]))+9999) 
            summaryV='band K -> none\n'

if outph:
  for ii in system_used:
    out = _telescope+'_'+str(coordinatelist)+'_'+str(_date)+'_'+str(_outputformat)+'_'+str(ii)
    fil = open(out,'w')
    fil.write(str(instrument)+' '+str(_date)+'\n')
    fil.write('*** '+coordinatelist+' '+str(len(standard))+'\n')
    if ii==0:
        print ii
        fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(1),str(1),str(1),str(1),str(1)))
        fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['U']),str(airmass['B']),str(airmass['V']),str(airmass['R']),str(airmass['I'])))
        for i in range(len(standard_value[ii])):
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                      % (str(outputformat[_outputformat]['U'][i]),str(outputformat[_outputformat]['B'][i]),str(outputformat[_outputformat]['V'][i]),str(outputformat[_outputformat]['R'][i]),str(outputformat[_outputformat]['I'][i]),str(standard_value[ii][i])))
    elif ii==2:
        print ii
        fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(1),str(1),str(1),str(1),str(1)))
        fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['u']),str(airmass['g']),str(airmass['r']),str(airmass['i']),str(airmass['z'])))
        for i in range(len(standard_value[ii])):
            fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                      % (str(outputformat[_outputformat]['u'][i]),str(outputformat[_outputformat]['g'][i]),str(outputformat[_outputformat]['r'][i]),str(outputformat[_outputformat]['i'][i]),str(outputformat[_outputformat]['z'][i]),str(standard_value[ii][i])))
    elif ii==1:
        print ii
        fil.write('%6.6s\t%6.6s\t%6.6s\n' % (str(1),str(1),str(1)))
        fil.write('%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['J']),str(airmass['H']),str(airmass['K'])))
        for i in range(len(standard_value[ii])):
            fil.write('%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                      % (str(outputformat[_outputformat]['J'][i]),str(outputformat[_outputformat]['H'][i]),str(outputformat[_outputformat]['K'][i]),str(standard_value[ii][i])))
    fil.close()

if _aperture:
    iraf.tv.rimexam.iterati = 3
    iraf.tv.rimexam.radius = 7

src.close_program(logincl)
