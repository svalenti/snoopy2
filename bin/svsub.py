#!/usr/bin/env python

from numpy import mean
from numpy import array
from numpy import compress
from numpy import std
from numpy import median
from snoopy2 import *
import getopt
import os,sys,shutil,string,re
import snoopy2
import glob
try:
        import pyfits
except:
        from astropy.io import fits as pyfits

import math 

####################  READ OPTION ###########################################
bb=[]
aa=glob.glob(snoopy2.__path__[0]+'/coordinate_std/optical/*list')
for i in aa:
    bb.append(i[len(snoopy2.__path__[0]+'/coordinate_std/optical/'):-5])

help ="################################################################       \n\
help = Usage:   svsn.py filename                                              \n\
                target            filelist (iraf format)                       \n\
                template          filelist (iraf format)                       \n\
                -l --list filename  list coordinate sequence stars            \n\
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
subdirectory=['optical/','infrared/','sloan/']

if len(sys.argv)<=2:
    print help
    sys.exit() 


logincl=src.open_program()
from pyraf import iraf

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)
iraf.rimexam.magzero = 0

if sys.argv[1][0]=='@':   
   ff = open(sys.argv[1][1:])
   files = ff.readlines()
   imglist = []
   for ff in files: 
         if not ff=='\n' and ff[0]!='#':
             imglist.append(re.sub('\n','',ff))
elif ',' in sys.argv[1]: imglist = string.split(sys.argv[1],sep=',')
else:
    if 'Error reading' in iraf.imstat(sys.argv[1],Stdout=1)[1]:
        print '\n####### Error ######'
        print 'If "'+str(sys.argv[1])+'" is an image, it is corrupted ...'
        print 'or you just forgot the "@" before the list name ...\n'
        src.close_program(logincl)
    else:
      imglist = [sys.argv[1]]

if sys.argv[2][0]=='@':   
   ff = open(sys.argv[2][1:])
   files = ff.readlines()
   imglist2 = []
   for ff in files: 
         if not ff=='\n' and ff[0]!='#':
             imglist2.append(re.sub('\n','',ff))
elif ',' in sys.argv[2]: imglist2 = string.split(sys.argv[2],sep=',')
else:
    if 'Error reading' in iraf.imstat(sys.argv[2],Stdout=1)[1]:
        print '\n####### Error ######'
        print 'If "'+str(sys.argv[2])+'" is an image, it is corrupted ...'
        print 'or you just forgot the "@" before the list name ...\n'
        src.close_program(logincl)
    else:
      imglist2 = [sys.argv[2]]

psflista=''
coordinatelist = ''

options, args = getopt.getopt(sys.argv[3:],'p:,i,s:,c,z:,b:,x:,y:,n:,r,l:',['psflist=','interactive','system=','coordinate','size=','background_region=','xorder=','yorder=','iteraction=','recenter','coordinatelist='])
for opt,arg in options:
    if opt in ('-l', '--list'): coordinatelist = arg
    if opt in ('-p', '--psf'): psflista = arg
    if opt in ('-i', '--interactive'): _interactive = True
    if opt in ('-s', '--system'): 
        system = int(arg)
        if system not in [0,1,2]:
            print 'ERROR: system not recognised !!'
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

if coordinatelist=='':
        print help
        sys.exit() 

######################## SET   IRAF   PARAMETERS  #######################
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)

from iraf import digiphot
from iraf import daophot
from iraf import ptools

zmag = 0 
iraf.digiphot.daophot.photpars.zmag = zmag
iraf.digiphot.daophot.datapars.datamin = 'INDEF' #-2000
iraf.digiphot.daophot.datapars.datamax = 'INDEF' #_datamax1    
#################################    CHECK HEADER    ######################################

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

warn='##################################\n'

for imglong in imglist:
  target=''
  if imglong:
      if imglong[-5:]=='.fits': img=imglong[:-5]
      _telescope=src.telescope(imglong)
      #print '### TELESCOPE= '+_telescope
      if not system:
          system=src.check_system(_telescope,imglong,Stdout=False)
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
      check=src.check_tel(imglong,_telescope,system)
      if check==0:
              print '####################### '
              print "#### Error with the header !!!!"
              print '### telescope not in the list '
              print '### if you want to continueanyway,'
              print '### run "svother.py" to correct the header '
              print '####################### '
      else:
          _header=src.read_parameter(_telescope,system)
          _filter=src.filter(imglong,_header,_telescope)
          filter0=src.filtername(_telescope,_filter,system)
          if filter0=='unknown':
              print '### WARNING: filter not recognised in this photometric system '
              filter0=_filter
          #print '##########  '+str(filter)+'  ##############'
          template=''
          for imglong2 in imglist2:
              if imglong2:
                  if imglong2[-5:]=='.fits': img2=imglong2[:-5]
                  _telescope2=src.telescope(imglong2)
                  _header2=src.read_parameter(_telescope2,system)
                  _filter2=src.filter(imglong2,_header2,_telescope2)
                  filter2=src.filtername(_telescope2,_filter2,system)
                  if filter2=='unknown': 
                      print '### WARNING: filter not recognised in this photometric system '
                      filter2=_filter2
                  #print '##########  '+str(filter2)+'  ##############'
                  if filter2==filter0:
                      template=img2
                      target=img
  if target and template:
      fwhmv={}
      exptimev={}
      imgpsf0={}
      apcov={}
      magimestar={}
      magphotstar={}
      fitmagstar={}
      fitmagstarerr={}
      fwhm_star={}
      airmassv={}

      print ''
      print '### start template subtraciton with '
      print '###'
      print '### target   = ',target
      print '### template = ',template
      print '### filter   = ',filter0
      print '###'
      print '### STEP 1 =  geometric registration and trim'
      print '###'

     #####   target
      _telescope=src.telescope(target+'.fits')
      _header=src.read_parameter(_telescope,system)
#      _datamin=src.ccdmin(target,_header,_telescope)
      _datamin1=-2000  #  for template subtraction the counts can go negative quite a bit
      _datamax1=src.ccdmax(target+".fits",_header,_telescope)

      exptimev['t_'+target+'.fits']=src.exptime(target+'.fits',_header,_telescope)
      airmassv['t_'+target+'.fits']=src.airmass(target+'.fits',_header,_telescope)

      _object=src.objects(target+'.fits',_header,_telescope)
      _airmass=src.airmass(target+'.fits',_header,_telescope)
     #####   teplate
      _telescope2=src.telescope(template+'.fits')
      _header2=src.read_parameter(_telescope2,system)
      _airmass2=src.airmass(template+'.fits',_header2,_telescope2)

      exptimev['temp_'+target+'.fits']=src.exptime(template+'.fits',_header2,_telescope2)
      airmassv['temp_'+target+'.fits']=src.airmass(template+'.fits',_header2,_telescope2)
      _datamin2=-2000  #  for template subtraction the counts can go negative quite a bit
      _datamax2=src.ccdmax(template+'.fits',_header2,_telescope2)

      stars,tmpfiltercoo=src.register_module(target+'.fits',system,coordinatelist,_interactive,logincl,filter0)
      stars,tmpfiltercoo=src.register_module(template+'.fits',system,coordinatelist,_interactive,logincl,filter2)

      answ='n'
      if os.path.isfile('t_'+target+'.fits') and os.path.isfile('temp_'+target+'.fits'):
          _z1,_z2,goon=src.display_image('t_'+target+'.fits',1,'','',False)
          _z1,_z2,goon=src.display_image('temp_'+target+'.fits',2,'','',False)
          if _interactive:
                  answ = raw_input('images already trimmed, do you want to skip this step [y/n] [y] ?')
                  if not answ:
                          answ='y'
          else:
                  answ = 'y'
      if answ in ['n','N','No','','NO']:
              src.delete("diff_"+target+'.fit?')
              src.delete("diff_"+target+'.psf.fit?')
              src.delete("T"+target)
              iraf.sregister(template,target,'T'+target+'.fits',wcs='world')
              _z1,_z2,goon=src.display_image('T'+target+'.fits',1,'','',True)
              src.delete("t_"+target+".fit?,temp_"+target+".fit?,temp_"+target+".co?")
              while answ in ['n','N','No','','NO']:
                  print ">>>  mark corners of trim region with  >b<, exit  >q<" 
                  print "     (the first one > b < low on the left, the second up to right) "
                  src.delete('tmptbl')
                  iraf.tvmark(1,"",logfile="tmptbl",autol='yes',mark="cross",interactive='yes')
                  answ=raw_input(' >>>  trim ok ? [y] ')
                  if not answ:
                          answ='y'
              x1,x2=iraf.field('tmptbl', '',Stdout=1)
              xb1,yb1=str(int(float(string.split(x1)[0]))),str(int(float(string.split(x1)[1])))
              xb2,yb2=str(int(float(string.split(x2)[0]))),str(int(float(string.split(x2)[1])))
              print ">> Trim frames: ",xb1," ",xb2," ",yb1," ",yb2 #, >> target+".log")
              iraf.imcopy(target+'.fits'+"["+xb1+":"+xb2+","+yb1+":"+yb2+"]","t_"+target+'.fits')
              iraf.imcopy("T"+target+'.fits'+"["+xb1+":"+xb2+","+yb1+":"+yb2+"]","temp_"+target+'.fits')

              _z1,_z2,goon=src.display_image('t_'+target,1,'','',False)
              _z1,_z2,goon=src.display_image('temp_'+target,2,'','',False)
              src.delete('tmp.temp_'+target+'.coo')      
              src.delete('tmp.t_'+target+'.coo')      
              iraf.wcsctran(coordinatelist+'.tv','tmp.t_'+target+'.coo','t_'+target,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
              iraf.wcsctran(coordinatelist+'.tv','tmp.temp_'+target+'.coo','temp_'+target,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
              _fwhm_target,_magime_targ,stars_pos_targ,_fwhm_star_targ=src.fwhm_computation('t_'+target+'.fits','tmp.t_'+target+'.coo',stars,system,_interactive)
              _fwhm_template,_magime_templ,stars_pos_templ,_fwhm_star_templ=src.fwhm_computation('temp_'+target+'.fits','tmp.temp_'+target+'.coo',stars,system,_interactive)
      else:
          src.delete('tmp.temp_'+target+'.coo')      
          src.delete('tmp.t_'+target+'.coo')      
          iraf.wcsctran(coordinatelist+'.tv','tmp.t_'+target+'.coo','t_'+target,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
          iraf.wcsctran(coordinatelist+'.tv','tmp.temp_'+target+'.coo','temp_'+target,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
          _fwhm_target=float(pyfits.getheader('t_'+target+'.fits')['QUBAFWHM'])
          stars_pos_targ=''
          _fwhm_template=float(pyfits.getheader('temp_'+target+'.fits')['QUBAFWHM'])
          stars_pos_templ=''

      _z1,_z2,goon=src.display_image('t_'+target,1,'','',False)
      _z1,_z2,goon=src.display_image('temp_'+target,2,'','',False)
      iraf.tvmark(1,'tmp.t_'+target+'.coo',mark="circle",number='yes',label='no',radii=20,nxoffse=15,nyoffse=15,color=203,txsize=4)
      iraf.tvmark(2,'tmp.temp_'+target+'.coo',mark="circle",number='yes',label='no',radii=20,nxoffse=15,nyoffse=15,color=205,txsize=4)

      fwhmv['t_'+target+'.fits']=_fwhm_target
      fwhmv['temp_'+target+'.fits']=_fwhm_template
##################################### PSF ###################################
      print '###'
      print '### STEP 2 =  PSF '
      print '###'
      if psflista:
          for _imgpsf in psflista:
                  print(_imgpsf)
          imgpsf0[0]=raw_input('which is the psf file for '+target+' ['+psflista[0]+']?')
          if not imgpsf0[0]: imgpsf0[0]=psflista[0]

      if _fwhm_target<= _fwhm_template:
          print("target has the best PSF \n both PSF are needed ")
          m=2
          if not psflista:
              imgpsf0[0]='t_'+target+'.psf'
              imgpsf0[1]='temp_'+target+'.psf'
      else:
          print "template has the best PSF \n  only target PSF is needed " 
          m=1
          if not psflista:
              imgpsf0[0]='t_'+target+'.psf'

      for jj in range(0,m):
          if jj==0: 
              _file='t_'+target
              _fw=_fwhm_target
          elif jj==1:
              _file = 'temp_'+target
              _fw = _fwhm_template
          if os.path.isfile(imgpsf0[jj]+'.fits'):
              print imgpsf0[jj]+'.fits #########################'
              print '##################'
              if _interactive:
                      answ=raw_input('psf file already in the directory, use it [y/n] ? [y] ')
                      if not answ:
                              answ='y'
              else:
                      answ='y'
          else:
              answ='n'
              print 'psf file not in the directory'
          if answ in ['N','No','NO','n']:
              print 'make psf'
              src.delete('tmp.'+_file+'.coo')
              iraf.wcsctran(coordinatelist+'.tv','tmp.'+_file+'.coo',_file,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
              if jj==0: 
                  if not stars_pos_targ:
                      _fwhm_target,_magime_targ,stars_pos_targ,_fwhm_star_targ=src.fwhm_computation(_file+'.fits','tmp.'+_file+'.coo',stars,system,_interactive)
                  _fw=_fwhm_target
                  stars_pos=stars_pos_targ
                  _magime=_magime_targ
                  _fw_star=_fwhm_star_targ
              elif jj==1:
                  if not stars_pos_templ:
                      _fwhm_template,_magime_templ,stars_pos_templ,_fwhm_star_templ=src.fwhm_computation(_file+'.fits','tmp.'+_file+'.coo',stars,system,_interactive)
                  _fw = _fwhm_template
                  stars_pos=stars_pos_templ
                  _magime=_magime_templ
                  _fw_star=_fwhm_star_templ
              ff = open('tmp.'+_file+'.coo','r')
              alllines=[]
              for j in range(0,3): exl=ff.readline()
              for j in range(len(stars_pos)): 
                  alllines.append(ff.readline())
              ff.close()
              ff = open('tmp.'+_file+'good.coo','w')
              for i in range(len(stars_pos)):
                  if stars_pos[i]:
                      ff.write(alllines[i])
              ff.close()

              _telescopev=src.telescope(_file+".fits")
              _headerv=src.read_parameter(_telescopev,system)
              _exptime=src.exptime(_file+".fits",_headerv,_telescopev)
              _airmass=src.airmass(_file+".fits",_headerv,_telescopev)
              _filter1=src.filter(_file+".fits",_headerv,_telescopev)
              _filter=src.filtername(_telescopev,_filter1,system)
              if _filter=='unknown':
                  _filter=_filter1
              ititle='xxx'
              print jj, _filter, _airmass, _exptime, _telescopev

##############################################################################
##############  DAOPHOT MAGNITUDE     ####################
################ AND APERTURE CORRECTION  ################# 

              _z1,_z2,goon=src.display_image(_file+'.fits',1,'','',False)
              varord=0
              annulus=4*_fw

              iraf.digiphot.daophot.fitskypars.annulus = annulus
              iraf.digiphot.daophot.photpars.apertures = _fw*3.
              iraf.digiphot.daophot.datapars.exposure = _headerv['hed_exptime']
              iraf.digiphot.daophot.datapars.airmass = _headerv['hed_airmass']
              iraf.digiphot.daophot.datapars.filter = _headerv['hed_filter1']

              src.delete(_file+".mag")
              iraf.digiphot.daophot.phot(image=_file, output=_file+'.mag', coords='tmp.'+_file+'good.coo', verify='no', interactive='no')
              _magphot2=iraf.digiphot.ptools.txdump(_file+".mag",fields="mag",expr='yes',Stdout=1)
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
              if _telescopev=='prompt':
                  prompttele=pyfits.getheader(_file+'.fits')['OBSERVAT']
                  if string.count(prompttele,'2'): _prompt='prompt2'
                  elif string.count(prompttele,'3'): _prompt='prompt3'
                  elif string.count(prompttele,'1'): _prompt='prompt1'
                  elif string.count(prompttele,'4'): _prompt='prompt4'
                  elif string.count(prompttele,'5'): _prompt='prompt5'
                  else: print 'WARNING: telescope not found !! '
                  _trimsection=src.trimsec(_file+".fits",_headerv,_telescopev)
                  if _trimsection:
                      xx=string.split(string.split(_trimsection[1:-1],',')[0],':')[0]
                      yy=string.split(string.split(_trimsection[1:-1],',')[1],':')[0]
                  else:
                      xx,yy=0,0
                  _coordinate_star=iraf.digiphot.ptools.txdump(_file+".mag",fields="XCENTER,YCENTER",expr='yes',Stdout=1)
                  _coord_origi_x,_coord_origi_y,_values=[],[],[]
                  _z10,_z20,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+_telescopev+'/'+_prompt+'.fits',2,'','',False)
                  if not goon: src.close_program(logincl)        
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
                              _values.append(float(string.split(iraf.imstat('home$coordinate_std/'+_telescopev+'/'+_prompt+'.fits['+str(_xx)+','+str(_yy)+']',Stdout=1)[1])[2]))
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
              iraf.digiphot.daophot.daopars.psfrad = annulus
              iraf.digiphot.daophot.daopars.fitrad = _fw
              iraf.digiphot.daophot.daopars.fitsky = 'yes'
              iraf.digiphot.daophot.daopars.sannulus = annulus
              iraf.digiphot.daophot.daopars.recenter = 'yes'
              iraf.digiphot.daophot.daopars.varorder = 0
    
              iraf.tvmark(1,'tmp.'+_file+'good.coo',mark="circle",number='yes',label='no',radii=20,nxoffse=15,nyoffse=15,color=203,txsize=4)

              print "-------------------------------------------------------------------------------"
              _id=iraf.noao.digiphot.ptools.txdump(_file+".mag",fields="ID",expr='yes',Stdout=1)
              _good=iraf.fields('tmp.'+_file+'good.coo','3',Stdout=1)
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

              src.delete(_file+".ps*,"+_file+".sub.fits")
              src.delete(_file+".grp,"+_file+".nst")
              src.delete(imgpsf0[jj]+'.fits')
              if _interactive:
                  iraf.noao.digiphot.daophot.psf(image=_file,photfile=_file+".mag",pstfile=_file+".mag",psfimage=imgpsf0[jj], mkstars='yes', opstfile=_file+".pst",groupfil=_file+".psg",verbose='yes',verify="no",interac="yes")
              else:
                  iraf.noao.digiphot.daophot.psf(image=_file,photfile=_file+".mag",pstfile=_file+".mag",psfimage=imgpsf0[jj], mkstars='yes', opstfile=_file+".pst",groupfil=_file+".psg",verbose='yes',verify="no",interac="no")
              apnum=iraf.noao.digiphot.ptools.txdump(_file+".pst",fields="ID",expr='yes',Stdout=1)      
              iraf.noao.digiphot.daophot.group(image=_file,photfile=_file+".mag",psfimage=imgpsf0[jj],groupfil=_file+".grp",verify='no',verbose='no')
              iraf.noao.digiphot.daophot.nstar(image=_file,groupfil=_file+".grp",psfimage=imgpsf0[jj],nstarfil=_file+".nst",rejfile="",verify='no',verbose='no')
              iraf.noao.digiphot.daophot.substar(image=_file,photfile=_file+".nst",exfile="",psfimage=imgpsf0[jj],subimage=_file+".sub",verify='no',verbose='no')
#############################################################################################################################################            
    
              _z1,_z2,goon=src.display_image(_file+".sub",1,_z1,_z2,False)
              if not goon: src.close_program(logincl)
    
              iraf.noao.digiphot.ptools.txsort(_file+".nst","ID")
              dddd=iraf.noao.digiphot.ptools.txdump(_file+".nst",fields="ID,mag,merr",expr='yes',Stdout=1)
              dddd2=iraf.noao.digiphot.ptools.txdump(_file+".nst",fields="ID",expr='yes',Stdout=1)

              j=0
              _fitmag,_fitmagerr=[],[]
              for i in range(len(stars)):
                  if stars_pos[i]:
                      try:
                          indice=dddd2.index(str(i-j+1))     #   number of star in img.nst
                          _fitmag.append(float(string.split(dddd[indice])[1]))
                          _fitmagerr.append(float(string.split(dddd[indice])[2]))
                      except: 
                          _fitmag.append(float(9999))
                          _fitmagerr.append(float(0.0))
                  else:
                      j=j+1
                      _fitmag.append(float(9999))
                      _fitmagerr.append(float(0.0))

              _apco,_apco2=[],[]
              for i in range(len(stars)):
                  if stars_pos[i]:
                      try: 
                          _id[_good.index(stars[i])] 
                          if _id[_good.index(stars[i])] in apnum:
                              if _magphot[i]>=99 or _fitmag[i]>= 99:
                                  _apco.append(9999)
                              else:
                                  _apco.append(_magphot[i]-_fitmag[i])
                              if _magime[i]>=99 or _fitmag[i]>= 99:
                                  _apco2.append(9999)
                              else:
                                  _apco2.append(_magime[i]-_fitmag[i])
                          else:
                              _apco.append(9999)
                              _apco2.append(9999)
                      except:
                          _apco.append(9999)
                          _apco2.append(9999)
                  else:
                      _apco.append(9999)
                      _apco2.append(9999)

              if _telescope=='prompt' and _values:
                  print '#####  Telescope = Prompt'
                  print '#####  Applied correction  !! '
                  for i in range(len(stars)):
                      if _magime[i]<=99:
                          _magime[i]=_magime[i]-_values[i]
                      if _magphot[i]<=99:
                          _magphot[i]=_magphot[i]-_values[i]
                      if _fitmag[i]<=99:
                          _fitmag[i]=_fitmag[i]-_values[i]
                                  
    #  magnitute imexam, phot , psf fit
              magimestar[_file+'.fits']=_magime
              magphotstar[_file+'.fits']=_magphot
              fitmagstar[_file+'.fits']=_fitmag
              fitmagstarerr[_file+'.fits']=_fitmagerr
              fwhm_star[_file+'.fits']=_fw_star
    
              print "*****************************************************"
              print "id      ime_mag   ap_mag   ph_mag  fit_mag   ap_cor  ph_cor"   
              for i in range(len(stars)):
                  print '\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s\t%8.8s' % (str(stars[i]),str(_magime[i]),str(_magphot[i]),str(_fitmag[i]),str(_apco2[i]),str(_apco[i]))

              aper_cor=compress((array(_apco)<999),array(_apco))
              aper_cor2=compress((array(_apco2)<999),array(_apco2))
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
                  if _interactive :
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
              apcov[_file+'.fits']=apcofin
              src.updateheader(_file+".fits",0,'QUBAAPCO',apcofin)
              src.updateheader(imgpsf0[jj]+'.fits',0,'QUBAAPCO',apcofin)
               ############  WRITE EC FILE ##############
               #    mag[filter]=_magime
              fil = open(_file+'.ec','w')
              fil.write('#file = '+str(_file)+'   title = '+str(_object)+'    filter='+str(filter0)+'  expt=%6.6s \n' % (str(exptimev[_file+'.fits'])))
              fil.write('#airmass= %6.6s \n'% (str(airmassv[_file+'.fits'])))
              fil.write('# PSF \n')
              fil.write('#    Reference FWHM= %6.6s    aperture correction =%6.6s \n' % (str(_fw),str(apcov[_file+'.fits'])))
              fil.write('# id fwhm  ph_mag  fit_mag  fiterr \n')
              for i in range(len(stars)):
                  fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(stars[i]),str(fwhm_star[_file+'.fits'][i]),str(magphotstar[_file+'.fits'][i]),str(fitmagstar[_file+'.fits'][i]),str(fitmagstarerr[_file+'.fits'][i])))
              fil.close()
          else:
              apcov[_file+'.fits']=float(pyfits.getheader(imgpsf0[jj]+'.fits')['QUBAAPCO'])
              apcofin=apcov[_file+'.fits']
              src.updateheader(_file+".fits",0,'QUBAAPCO',apcofin)
              src.updateheader(imgpsf0[jj]+'.fits',0,'QUBAAPCO',apcofin)

#####################################################################################
############################### IMAGE DIFFERENCE #############################
#isis:
      print '###'
      print '### STEP 3 =  ISIS '
      print '###'
      answ='n'
      if os.path.isfile('diff_'+target+'.fits'):
              if _interactive :
                      answ=raw_input('difference image already in the directory, skip this step [y/n] ? [y] ')
                      if not answ:
                              answ='y'
              else:
                      answ = 'y'
      else:
          src.delete("diff_"+target+".fit?")
          answ='n'

      _z1,_z2,goon=src.display_image('t_'+target+'.fits',1,'','',False)
      _z1,_z2,goon=src.display_image('temp_'+target+'.fits',2,'','',False)

      exptarg = exptimev['t_'+target+'.fits']
      exptemp = exptimev['temp_'+target+'.fits']
      airtarg=airmassv['t_'+target+'.fits']
      airtemp=airmassv['temp_'+target+'.fits']

      _fwhm_target=fwhmv['t_'+target+'.fits']
      _fwhm_template=fwhmv['temp_'+target+'.fits']

      if _fwhm_target<= _fwhm_template:
          isa = "t_"+target+'.fits'
          isb=  "temp_"+target+'.fits'
          exptime = exptarg
          airmass = airtarg
          _imgpsf0=imgpsf0[1]

      else: 
          isa = "temp_"+target+'.fits'
          isb = "t_"+target+'.fits'
          exptime = exptemp
          airmass = airtemp
          _imgpsf0=imgpsf0[0]

      try:
          enviro=os.environ['IRAFARCH']
          if enviro=='macintel':
              directory_prog=snoopy2.__path__[0]+'/exec_prog/macintel'
          elif enviro=='redhat' or enviro=='linux':
              directory_prog=snoopy2.__path__[0]+'/exec_prog/redhat'
          else:
              print 'ERROR: architecture not found !!'
              directory_prog=''
      except:
          directory_prog=snoopy2.__path__[0]+'/exec_prog/redhat'
              
      _nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp=9,9,1,1,60000,1,1
#      nstamps_x,nstamps_y,sub_x,sub_y,saturation,degbg,degsp=src.isis_parameter(_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp)
      while answ=='n':
          sumkern=src.run_isis(isa,isb,_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp,directory_prog)
          src.delete("diff_"+target+".fit?")
          src.delete("diff_"+target+".sn.co?") 
          print "********************************************************" 

          if _fwhm_target<= _fwhm_template:
              radphot = _fwhm_template*3.
              os.system('mv conv.fits '"diff_"+target+'.fits')
              #iraf.imrename("conv","diff_"+target, verbose='yes')   
              print "!!!!!  Photometric reference "+isa+"    !!!!!!!"
          else: 
              radphot = _fwhm_target*3.  
              iraf.imarith("conv.fits","*","-1","diff_"+target+".fits")
              #iraf.imarit("conv","*","-1","diff_"+target)
              print "!!!!!!! Photometric reference "+isa+" !!!!!!!!"
              print "********************************************************" 
          
          src.updateheader("diff_"+target+".fits",0,_header['hed_exptime'],exptime)
          src.updateheader("diff_"+target+".fits",0,_header['hed_airmass'],airmass)
          src.updateheader("diff_"+target+".fits",0,_header['hed_filter1'],filter0)
   
          print "##Difference image diff_"+target+" on frame 3  #####" 
          _tm1,_tm2,goon=src.display_image("diff_"+target+'.fits',3,'','',False)

          if _interactive :
                  answ=raw_input(">>> are you happy of difference [y/n] ? [y] ")
                  if not answ:
                          answ='y'
          else:
                  answ='y'
                  
          if answ in ['N','No','NO','n']:
              answ2=raw_input(">>> do you want to calculate with new parameters [y/n] [y] ? ")
              if not answ2: answ2='y'
              if answ2 in ['Y','Yes','YES','yes','y']:
                  _nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp=src.isis_parameter(_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp)
                  answ='n'
              else:
                  src.close_program(logincl)
      if _fwhm_target>= _fwhm_template:
          src.delete('temp_sc.fits')
          try:
              _gain=src.gain(template+'.fits',_header2,_telescope2)
              _ron=src.ron(template+'.fits',_header2,_telescope2)
              _datamin2=src.ccdmin(template+'.fits',s_header2,_telescope2)
              _datamax2=src.ccdmax(template+'.fits',_header2,_telescope2)
              os.system('cp conv0.fits temp_sc.fits')
              ff=pyfits.getheader('temp_'+target+'.fits')
              src.updateheader('temp_sc.fits',0,_header2['hed_exptime'],ff[_header2['hed_exptime']])
              src.updateheader('temp_sc.fits',0,_header2['hed_airmass'],ff[_header2['hed_airmass']])
              src.updateheader('temp_sc.fits',0,'telescop',ff['telescop'])
              src.updateheader('temp_sc.fits',0,_header2['GAIN_value'],_gain)
              src.updateheader('temp_sc.fits',0,_header2['RON_value'],_ron)
              src.updateheader('temp_sc.fits',0,_header2['GAIN_value'],_gain)
              src.updateheader('temp_sc.fits',0,_header2['datamin_value'],_datamin2)
              src.updateheader('temp_sc.fits',0,_header2['datamax_value'],_datamax2)
          except:
              print 'ISIS not runned again, DM from sequence stars not possible'
      #iraf.delete("tmp$*,conv.fit?,conv0.fit?,kc_*.fits,kt_*.fits",verify='no')

##########################################################################
################# artificial star experiment  for calibration ############
      print '###'
      print '### STEP 4 =  STAR EXPERIMENT FOR CALIBRATION  Z_target - Z-tmeplate= DM'
      print '###'
      exptarg = exptimev['t_'+target+'.fits']
      exptemp = exptimev['temp_'+target+'.fits']
      print imgpsf0[0]
      try:
          DM=float(pyfits.getheader('diff_'+target+'.fits')['QUBADM'])
      except:
          DM=''
      answ='y'
      if DM:
          print '######  '+str(DM)
          print 'DM already measured '
          answ=raw_input('do you want to measure DM again [y,n] [n] ? ')
          if not answ: answ='n'
      if answ in ['y','Y','YES','Yes','yes']:
          if _fwhm_target<= _fwhm_template:
              _zeromag = 0
              src.updateheader('diff_'+target+".fits",0,'QUBADM',_zeromag)
              print '#####'
              print '#####  template scaled to target '
              print "#####   DM = 0 "
              print '#####'
          else:
    #   DM taked from ISIS parameter sum_kernel (allard communication)
              DM = 2.5*math.log10(sumkern)+2.5*math.log10(exptemp)-2.5*math.log10(exptarg)
              print "********************************************************"
              print '>>>   zero point mag difference = '+str(DM)
              print "********************************************************"
              src.updateheader('diff_'+target+".fits",0,'QUBADM',DM)

###########################################################################
##   SN measurement
###########################################################################
      print '###'
      print '### STEP 5 =  SN MEASUREMENT '
      print '###'
      answ='n'
      exptarg = exptimev['t_'+target+'.fits']
      exptemp = exptimev['temp_'+target+'.fits']
      airtarg = airmassv['t_'+target+'.fits']
      airtemp = airmassv['temp_'+target+'.fits']
      try:
          DM = float(pyfits.getheader('diff_'+target+'.fits')['QUBADM'])
      except:
          DM = ''
      try:
          QUBASN1 = pyfits.getheader('t_'+target+'.fits')['QUBASN1']
      except:
          QUBASN1=''
      if QUBASN1:
             print(QUBASN1)
             if _interactive :
                     answ=raw_input('magnitude already measured, do you want to skip this step [y/n] [y] ? ')
                     if not answ:
                             answ='y'
             else:
                     answ='y'
      else:
              answ='n'
      if answ in ['n','N','no','','NO','No']:
            print(imgpsf0[0])
            src.delete('diff_'+target+'.psf.fits')
            os.system('cp '+_imgpsf0+'.fits diff_'+target+'.psf.fits')
            img='diff_'+target
            if _fwhm_target<= _fwhm_template:
                fwhm0=float(_fwhm_template)
                nor=exptarg/exptemp
                apco0=apcov['temp_'+target+'.fits']
            else:
                nor=exptemp/exptarg
                fwhm0=float(_fwhm_target)
                apco0=apcov['t_'+target+'.fits']
         
            # GAIN,RON, DATAMAX DATAMIN FROM target .....to be checked !!!!!
            _telescope=src.telescope(target+'.fits')
            _header=src.read_parameter(_telescope,system)
            _gain=src.gain(target+'.fits',_header,_telescope)
            _ron=src.ron(target+'.fits',_header,_telescope)
            _datamin1=src.ccdmin(target+'.fits',_header,_telescope)
            _datamax1=src.ccdmax(target+'.fits',_header,_telescope)
            iraf.digiphot.daophot.datapars.readnoi = _ron
            iraf.digiphot.daophot.datapars.epadu = _gain
            iraf.digiphot.daophot.datapars.datamin = -2000
            iraf.digiphot.daophot.datapars.datamax = 1e9 #_datamax1
            print _gain,_ron,_datamax1    
         
            a1 = int(fwhm0+.5)
            a2 = int(2.*fwhm0+.5)
            a3 = int(3.*fwhm0+.5)
            a4 = int(4.*fwhm0+.5)
            a5 = int(5.*fwhm0+.5)
            ap = str(a1)+","+str(a2)+","+str(a3)
         
              ##  fitskypars.salgorithm = "constant"
              ##  fitskypars.skyvalue = 0
         
            iraf.daophot.fitskypars.annulus=a3
            iraf.daophot.photpars.apertures = ap
            iraf.digiphot.daophot.datapars.exposure = _header['hed_exptime']
            iraf.digiphot.daophot.datapars.airmass = _header['hed_airmass']
            iraf.digiphot.daophot.datapars.filter = _header['hed_filter1']
         
   #########################    POINT TO SN #######################################
         
            iraf.set(stdimage='imt2048')
            _z1,_z2='',''
            if _interactive:
                _z1,_z2,goon=src.display_image(img+'.fits',1,'','',True)
            else:
                _z1,_z2,goon=src.display_image(img+'.fits',1,'','',False)
              
            if coordinate:
                xx0,yy0 = src.sn_coordinate('t_'+target+'.fits',_interactive)
            else:
                xx0=''
         
            if xx0 and not _interactive:
              x=xx0
              y=yy0
            else:
              print "_____________________________________________"
              print "  MARK SN REGION WITH - x -, EXIT  - q -"
         
              repeat='y'
              while repeat=='y':
                try:
                    x,y,value=string.split(iraf.imexamine(img+'.fits',1,wcs='logical',xformat='',yformat='',Stdout=1)[0])
                    src.delete("tmplabel")
                    ff = open('tmplabel','w')
                    ff.write(str(x)+' '+str(y)+' '+str(value)+' \n')
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
                    if not repeat: repeat='n'
                    elif repeat=='yes': repeat='y' 
                    elif repeat=='YES': repeat='y' 
                    print "_____________________________________________"
                    print "  MARK SN REGION WITH - x -, EXIT  - q -"
          
            if _telescope=='prompt':
                     prompttele=pyfits.getheader(target+'.fits')['OBSERVAT']
                     if string.count(prompttele,'2'): _prompt='prompt2'
                     elif string.count(prompttele,'1'): _prompt='prompt1'
                     elif string.count(prompttele,'3'): _prompt='prompt3'
                     elif string.count(prompttele,'4'): _prompt='prompt4'
                     elif string.count(prompttele,'5'): _prompt='prompt5'
                     else: print 'WARNING: telescope not found !! '
                     imgin_cor='home$coordinate_std/'+_telescope+'/'+_prompt+'.fits'
                     _z1,_z2,goon=src.display_image(target+'.fits',1,'','',False)
                     print ' telescope prompt: mark position on the original image to correct the final magnitude'
                     if coordinate:
                         xx0,yy0=src.sn_coordinate(target+'.fits')
                     else:
                         xx0=''
         
                     if xx0 and not _interactive:
                         xprompt=xx0
                         yprompt=yy0
                     else:
                         print "_____________________________________________"
                         print "  MARK SN REGION WITH - x -, EXIT  - q -"
         
                         repeat='y'
                         while repeat=='y':
                             try:
                                 xprompt,yprompt,valuepompt=string.split(iraf.imexamine(target+'.fits',1,wcs='logical',xformat='',yformat='',Stdout=1)[0])
                                 src.delete("tmplabel")
                                 ff = open('tmplabel','w')
                                 ff.write(str(xprompt)+' '+str(yprompt)+' '+str(value)+' \n')
                                 ff.close()
                                 iraf.tvmark(1,'tmplabel',autol='no',mark="cross",inter='no',label='no',txsize=4)
                                 repeat = raw_input('### repeat selection ? [y/n] ? [n] ')
                                 if not repeat: repeat='n'
                                 elif repeat=='yes': repeat='y' 
                                 elif repeat=='YES': repeat='y' 
                             except:
                                 xprompt=yprompt=valueprompt=0
                                 print '### WARNING: SN REGION NOT SELECTED !!!'
                                 repeat = raw_input('### repeat selection ? [y/n] ? [n] ')
                                 if not repeat: repeat='n'
                                 elif repeat=='yes': repeat='y' 
                                 elif repeat=='YES': repeat='y' 
                                 print "_____________________________________________"
                                 print "  MARK SN REGION WITH - x -, EXIT  - q -"
                                 
                     _trimsection=src.trimsec(target+'.fits',_header,_telescope)
                     if _trimsection:
                         print _trimsection
                         src.delete("prompt_trimmed.fits")
                         iraf.imcopy (imgin_cor+_trimsection,"prompt_trimmed.fits")
                     else:
                         src.delete("prompt_trimmed.fit?")
                         os.system('cp '+imgin_cor+".fits prompt_trimmed.fits")
                     _values=string.split(iraf.imstat('prompt_trimmed.fits['+str(int(float(xprompt)))+','+str(int(float(yprompt)))+']',Stdout=1)[1])[2]
                     src.delete("prompt_trimmed.fit?",ve='no')
            else:
                _values=0
   ###############################################################################
         
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
                 iraf.imcopy (img+'.fits'+"["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]","original.fits")
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
                 iraf.imcopy (img+'.fits'+"["+str(x1)+":"+str(x2)+","+str(y1)+":"+str(y2)+"]","original.fits")
                 iraf.set(stdimage='imt512')
                 _tmp1,_tmp2,goon=src.display_image('original.fits',1,_z1,_z2,False,_xsize=.5,_ysize=.5)
                 repeat = raw_input('### ok ? [y/n] ? [y] ')
                 if not repeat: repeat='y'
                 elif repeat=='no': repeat='n'
         
            _z1,_z2,goon=src.display_image('original.fits', 1,'','', False, _xsize=.5, _ysize=.5)

            if _interactive:
              answ = 'y'
              answ= raw_input(">>>>> Cuts OK [y/n] [y]?")
              if not answ: answ='y'
              elif answ=='no': answ='n' 
         
              while answ == 'n':  
                 _z11=float(_z1)
                 _z22=float(_z2)
                 z11 = raw_input('>>> z1 = ? ['+str(_z11)+'] ? ')
                 z22 = raw_input('>>> z2 = ? ['+str(_z22)+'] ? ')
                 if not z11: z11=_z11
                 else: z11=float(z11)
                 if not z22: z22=_z22
                 else: z22=float(z22)
                 _tm1,_tm2,goon=src.display_image("original.fits",1,z11,z22,False)
                 answ= raw_input(">>>>> Cuts OK [y/n] [y]?")
                 if not answ: answ='y'
                 if answ in ['N','n','NO','no','No']:
                     answ='n'
                 else:
                     answ='y'
                     _z1,_z2=_tmp1,_tmp2

            z11=float(_z1)
            z22=float(_z2)
         
            if not _interactive and xx0:
                aa=float(pyfits.getheader('original.fits')['NAXIS1'])/2
                bb=float(pyfits.getheader('original.fits')['NAXIS2'])/2
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
                _tm1,_tm2,goon=src.display_image("original.fits",1,z11,z22,False)
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
         
   ############################### BACKGROUND FIT   ###############################
         
   ## ##    if not _interactive or _fb:  
   ## ##    else:
            print ' ************  background fit **********************'
            answ0 = 'n'
            while answ0 == 'n':
                src.delete("sky.*,bg.*,bgs.*,sn.*,residual.*")
                #iraf.display("original.fits",1,fill='yes',xsize=.5,ysize=.5, zrange='no', zscale='no', z1=z11, z2=z22)
                _tmp1,_tmp2,goon=src.display_image('original.fits',1, z11, z22, False, _xsize=.5, _ysize=.5)
                nax=int(pyfits.getheader('original.fits')['naxis1'])
                nay=int(pyfits.getheader('original.fits')['naxis2'])
         
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
                            checkorder='yes'
                    iraf.imsurfit("original.fits","bg.fits",xorder=xbgord0,yorder=ybgord0,regions="section",section="sec")  
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
                    
                    iraf.imsurfit("original.fits","bg.fits",xorder=xbgord0,yorder=ybgord0,regions="sections",sections="sec")  
         
                midpt=float(iraf.imstat("bg.fits",field="mean", Stdout=1)[1])
         
                iraf.imcopy("original.fits","sky.fits")
                iraf.imcopy(inp,"bgs.fits")
                iraf.imcopy("bgs.fits",out)
                iraf.imarith("original.fits","-","sky.fits","sn.fits") 
                iraf.imarith("sn.fits","+",midpt,"sn.fits")
         
                answ0 = 'y'
                print answ0
                _tmp1,_tmp2,goon=src.display_image('original.fits',1,z11,z22,False,_xsize=.3,_ysize=.3,_xcen=.25,_ycen=.25)
                
                s1 = 1
                s2 = -int(fwhm0)
                src.delete("tmptbl")
                ff=open('tmptbl','w')
                ff.write(str(s1)+' '+str(s2)+" ORIGINAL")
                ff.close()
                iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
                _tmp1,_tmp2,goon=src.display_image('sky.fits',1, z11, z22, False, _erase='no', _xsize=.3,_ysize=.3,_xcen=.25,_ycen=.75)
                
                src.delete("tmptbl")
                ff=open('tmptbl','w')
                ff.write(str(s1)+' '+str(s2)+" BACKGROUND_FIT") 
                ff.close()
                iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
                _tmp1,_tmp2,goon=src.display_image('sn.fits',1, z11, z22, False, _erase='no', _xsize=.3,_ysize=.3,_xcen=.75,_ycen=.25)
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
###           
   #######################################    FITSN        ###################################  
            if _fwhm_target<= _fwhm_template: 
                _imgpsf1 = imgpsf0[1]
            else: 
                _imgpsf1 = imgpsf0[0]
            print _imgpsf1
            print imgpsf0
            _imgpsf=raw_input('which psf ? '+str(_imgpsf1)+' ')
            if not _imgpsf: _imgpsf=_imgpsf1
            apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy=src.fitsn(_recenter,img,_imgpsf,fwhm0,apco0,z22,z11,midpt,size,nor,_values,DM)
          
   ####################       Iterate Beckground    ###################################
         
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
                src.delete("sn.*,residual.*,snfit.*,tmp.*")
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
                        if not xbgord0:
                                xbgord0=_xbgord0
                        else:
                                _xbgord0=xbgord0
                        ybgord0 = raw_input('>>> Order of function in x for bg fit ['+str(_ybgord0)+'] ? ')
                        if not ybgord0:
                                ybgord0=_ybgord0
                        else:
                                _ybgord0=ybgord0
                    try:
                          float(xbgord0)
                          float(ybgord0)
                          checkorder='no'
                    except:  
                        print 'WARNING: value not valid !!'
                        checkorder='yes'	   
         
                iraf.imsurfit("skyfit.fits","tmp.fits",regions="all",xorder=xbgord0,yorder=ybgord0)
                midpt=float(iraf.imstat("tmp.fits",field="mean", Stdout=1)[1])
                iraf.imarith("original.fits","-","tmp.fits","sn.fits")
                iraf.imarith("sn.fits","+",midpt,"sn.fits")
                src.delete("skyfit.fits")
                if _fwhm_target<= _fwhm_template: _imgpsf = imgpsf0[1]
                else: _imgpsf = imgpsf0[0]
                
                apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy=src.fitsn(_recenter,img,_imgpsf,fwhm0,apco0,z22,z11,midpt,size,nor,_values,DM)
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
                  
   ##########################################################################
         
            print "***************************************************************************"
###         print "#id  x_ori    y_ori      x       y      ap_ori   ap_bgsub   fit_mag   err_art err_fit"
            for i in range(len(fitmag)):
                print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(truemag[i]),"  ",str(arterr),"  ",str(magerr[i])
            print "**************************************************************************"
         
   #############            AGGIUSTAMENTO MANUALE                     ###############
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
                    dmag0=raw_input(">>> D(mag) adjustment (positive=fainter) ["+str(dmag0)+"]")
                    if not dmag0: dmag0=0
                    try: 
                        float(dmag0)
                        checkdm='no' 
                    except:
                        checkdm='yes'
                apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,newmag=src.manu(dmag0,apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,z22,z11,midpt,size,fwhm0,img,x1,y1,arterr)
                dmag0=newmag[0]-truemag[0]
                answ0=raw_input(">>> again ? [y/n] [y] ")
                if not answ0: answ0='y'
                elif answ0=='yes': answ0='y'
         
            truemag=list(array(newmag))
         
            print '***************************************************************************' 
            print "#id  x_ori   y_ori     x     y    ap_ori ap_bgsub  fit_mag  err_art  err_fit"
            print "# id   ap_original ap_bgsub  fit_mag  err_art  err_fit"#,  >> nome0//".ec"
            print "# SN_FIT  "#, >> nome0//".ec"
            print "# id ap_ori ap-bg  fit_mag"#, >> nome0//".ec"
            for i in range(len(fitmag)):
                print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(truemag[i]),"  ",str(arterr),"  ",str(magerr[i])
            print "**************************************************************************"
         
            for i in range(len(fitmag)):
                misu='SN%1.1s %6.6s %6.6s %6.6s %6.6s %6.6s\n' % (str(i+1),str(apori3[i]),str(apmag3[i]),str(truemag[i]),str(arterr),str(magerr[i]))
                src.updateheader('t_'+target+".fits",0,'QUBASN'+str(i+1),misu[3:-1])
                src.updateheader(target+".fits",0,'QUBASN'+str(i+1),misu[3:-1])

###############################################################
##########################################################################
################# artificial star experiment  for error ############
      print '###'
      print '### STEP 6 =  ARTIFICIAL STAR EXPERIMENT'
      print '###'
      exptarg = exptimev['t_'+target+'.fits']
      exptemp = exptimev['temp_'+target+'.fits']
      airtarg = airmassv['t_'+target+'.fits']
      airtemp = airmassv['temp_'+target+'.fits']
      print imgpsf0[0]
      try:
          QUBASN1=pyfits.getheader('t_'+target+'.fits')['QUBASN1']
      except:
          QUBASN1=''
      answ='y'
      if not QUBASN1:
          print 'warning: magnitude not measured (QUBASN1 not in the header)'
          #src.close_program(logincl)
      else:
          try:
              xxx=apori3
          except:
              truemag,apori3,apmag3,arterr,magerr=[],[],[],[],[]
              apori3.append(string.split(QUBASN1)[0])
              apmag3.append(string.split(QUBASN1)[1])
              truemag.append(string.split(QUBASN1)[2])
              arterr=string.split(QUBASN1)[3]
              magerr.append(string.split(QUBASN1)[4])

          if string.split(QUBASN1)[3]!='0':
              print QUBASN1,len(QUBASN1)
              print 'error already measured= '+string.split(QUBASN1)[3]
              answ=raw_input('do you want to measure again [y,n] [n] ? ')
              if not answ: answ='n'
          else:
              answ='y'
      if answ in ['y','Y','YES','Yes','yes']:
          _z1,_z2,goon=src.display_image('t_'+target,1,'','',True)
          if _fwhm_target<= _fwhm_template:
                nor=exptarg/exptemp
                exptime = exptarg
                airmass = airtarg
                isa = "_targart.fits"
                isb=  "temp_"+target+'.fits'
                fwhm0=float(_fwhm_template)
          else:
              nor=exptemp/exptarg
              isa = "temp_"+target+'.fits'
              isb = "_targart.fits"
              exptime = exptemp
              airmass = airtemp
              fwhm0=float(_fwhm_target)

          print ">>> Give position for artificial stars for error estimation"
          print "(> x < then > q <)"
          xy=iraf.imexamine('t_'+target+'.fits',1,wcs='logical',xformat='',yformat='',Stdout=1)
          iraf.digiphot.daophot.datapars.exposure = _header['hed_exptime']
          iraf.digiphot.daophot.datapars.airmass = _header['hed_airmass']
          iraf.digiphot.daophot.daopars.psfrad =  4.*fwhm0
          iraf.digiphot.daophot.daopars.fitrad = fwhm0
          iraf.digiphot.daophot.daopars.fitsky = 'yes'
          iraf.digiphot.daophot.daopars.recenter = 'yes'
          iraf.digiphot.daophot.daopars.varorder = 0
          iraf.digiphot.daophot.datapars.fwhmpsf = fwhm0
          iraf.digiphot.daophot.centerpars.cbox = fwhm0*2.  
          iraf.digiphot.daophot.centerpars.calgori = "centroid"
          iraf.digiphot.daophot.fitskypars.salgori = "centroid"  
          iraf.digiphot.daophot.fitskypars.annulus = 4.*fwhm0
          iraf.digiphot.daophot.photpars.apertures = fwhm0*3. 
          iraf.digiphot.daophot.datapars.datamin = 'INDEF' #-2000
          iraf.digiphot.daophot.datapars.datamax = 'INDEF' #_datamax1    
          answ='n'
          while answ=='n':
                  print ">>> magnitude for the artificial stars " 
                  print '>>> SN magnitude= '+string.split(QUBASN1)[2]
                  _artmag=raw_input(">>> mag ["+string.split(QUBASN1)[2]+"] ? ")
                  if not _artmag: _artmag=float(string.split(QUBASN1)[2])
                  else: _artmag=float(_artmag)
                  ff=open("_art.list",'w')
                  for st in range(1,len(xy)):
                      x,y,value=string.split(xy[st])
                      ff.write(str(x)+' '+str(y)+' '+str(_artmag)+' \n')
                  ff.close()
                  src.delete("_targart.*")
                  if _fwhm_target<= _fwhm_template: _imgpsf = imgpsf0[1]
                  else: _imgpsf = imgpsf0[0]
                  #   artificial stars always with target psf !!!! 
                  #iraf.addstar("t_"+target,"_art.list",_imgpsf,"_targart",simple='yes',veri='no')  
                  iraf.addstar("t_"+target+'.fits',"_art.list",imgpsf0[0],"_targart.fits",simple='yes',veri='no')  

                  sumkern=src.run_isis(isa,isb,_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp,directory_prog)
                  src.delete("diff_targart.fit?")
                  src.delete("diff_targart.sn.co?") 
                  print "********************************************************" 
                  if _fwhm_target<= _fwhm_template:
                      radphot = _fwhm_template*3.
                      os.system('mv conv.fits diff_targart.fits')
                      print "!!!!!  Photometric reference "+isa+"    !!!!!!!"
                  else: 
                      radphot = _fwhm_target*3.  
                      iraf.imarith("conv.fits","*","-1","diff_targart.fits")

                  print "!!!!!!! Photometric reference "+isa+" !!!!!!!!"
                  src.updateheader("diff_targart.fits",0,_header['hed_filter1'],filter0)
                  src.updateheader("diff_targart.fits",0,_header['hed_exptime'],exptime)
                  src.updateheader("diff_targart.fits",0,_header['hed_airmass'],airmass)
                  print "********************************************************" 
        
                  print "##Difference image diff_targart on frame 3  #####" 
                  _tm1,_tm2,goon=src.display_image("diff_targart.fits",3,'','',False)
                  src.delete("artlist.mag") 
                  iraf.phot("diff_targart","_art.list","artlist.mag",veri='no',inter='no')   
                  _artmag2=iraf.noao.digiphot.ptools.txdump("artlist.mag",fields="mag",expr='yes',Stdout=1)
                  for val in range(len(_artmag2)):
                      try:
                          _artmag2[val]=float(_artmag2[val])
                      except:
                          _artmag2[val]=float(9999)
                  _artmag3=compress(array(_artmag2)<=99,array(_artmag2))
                  if len(_artmag3)>=1:
                      print _artmag3-_artmag
                      answ='y'
                  else:
                      print '### stars magnitude not good (INDEF)'
                      answ2=raw_input(">>> do you want to run ISIS with new parameters (fit and ISIS) [y/n] [y] ? ")
                      if not answ2: answ2='y'
                      if answ2 in ['Y','Yes','YES','yes','y']:
                          annulus2=5
                          aperture2=4
                          _annulus2,_aperture2,_datamin2,_datamax2=src.daophot_parameters(_datamin1,_datamax1,aperture2,annulus2,fwhm0)

                          iraf.digiphot.daophot.fitskypars.annulus = float(_annulus2)*fwhm0
                          iraf.digiphot.daophot.photpars.apertures = float(_aperture2)*fwhm0
                          iraf.digiphot.daophot.datapars.datamin = float(_datamin2)
                          iraf.digiphot.daophot.datapars.datamax = float(1e11)#float(_datamax2) 

                          _nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp=src.isis_parameter(_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp)
                          answ='n'
                      else:
                          src.close_program(logincl)                          

          src.delete("tmp$*,conv.fit?,conv0.fit?,kc_*.fits,kt_*.fits")
          _zeromag =  mean(_artmag3)
          _zeroerr =  std(_artmag3)
        
          print '>>>   artificial stars mag = '+str(_zeromag)
          print  '>>>                 error ='+str(_zeroerr)
          print "********************************************************" 
          if len(_artmag3)>=5:
                  print ".... min-max rejection ...."
                  _artmag4=compress((mean(_artmag3)-std(_artmag3)<=_artmag3)&(_artmag3<=mean(_artmag3)+std(_artmag3)),_artmag3)
                  _zeromag2 =  mean(_artmag4)
                  _zeroerr2 =  std(_artmag4)
                  rejection=len(_artmag3)-len(_artmag4)
                  print '>>>   star rejected  = '+str(rejection)+' of '+str(len(_artmag3))+'  stars'
                  print '>>>   artificial magnitude  = '+str(_zeromag2)
                  print  '>>>                  error ='+str(_zeroerr2)
                  print "********************************************************"
                  arterr=raw_input('arterr = ['+str(_zeroerr2)+'] ? ') 
                  if not arterr: arterr=_zeroerr2
          else:
                  arterr=raw_input('arterr = ['+str(_zeroerr)+'] ? ') 
                  if not arterr: arterr=_zeroerr

          QUBASN2='%s %s %s %6.6s %s'%(string.split(QUBASN1)[0],string.split(QUBASN1)[1],string.split(QUBASN1)[2],str(arterr),string.split(QUBASN1)[4])
          src.updateheader('t_'+target+".fits",0,'QUBASN1',QUBASN2)
          src.updateheader(target+".fits",0,'QUBASN1',QUBASN2)
      print "**************************************************************************"        
      snmesu="# id   ap_original ap_bgsub  fit_mag  err_art  err_fit\n"
      snmesu=snmesu+"# SN_FIT  \n"
      snmesu=snmesu+"# id ap_ori ap-bg  fit_mag \n"
      for i in range(len(truemag)):
          misu='SN%1.1s %6.6s %6.6s %6.6s %6.6s %6.6s\n' % (str(i+1),str(apori3[i]),str(apmag3[i]),str(truemag[i]),str(arterr),str(magerr[i]))
          src.updateheader('t_'+target+".fits",0,'QUBASN'+str(i+1),misu[3:-1])
          src.updateheader(target+".fits",0,'QUBASN'+str(i+1),misu[3:-1])
          
          snmesu=snmesu+misu          
          ff=open('t_'+target+'.ec','a')
          for i in snmesu:
              ff.write(i)
          ff.close()
      
      src.delete("T"+target+'.fits')
      src.delete("apori")
      src.delete("sec")
      src.delete("skyfit.fits")
      src.delete("sn.fits")
      src.delete("bg.fits,bgs.fits")
      src.delete("tmp*")
      src.delete('t_'+target+".sn.*")

#############################################################################
  else:
      print 'Warning: no template found for '+img+' !'

src.close_program(logincl)
