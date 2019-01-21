#!/usr/bin/env python

from numpy import mean
from numpy import array
from numpy import compress
from astropy.io import fits as pyfits
#from pyfits import getheader
from snoopy2 import *
import getopt
import os,sys,shutil,string,re
import snoopy2
import glob
userPath = os.getcwd()

bb=[]
aa=glob.glob(snoopy2.__path__[0]+'/coordinate_std/optical/*list')
for i in aa:
    bb.append(i[len(snoopy2.__path__[0]+'/coordinate_std/optical/'):-5])
    
help ="################################################################ \n\
help = Usage:   svstandard.py filename  -l coordinatelist       \n\
                input              filelist iraf format         \n\
                -l,--list   file coordinate                     \n\
                [-o,--output filename] output file              \n\
                [-s,--system value]     Specific photometric system  0 opt 1 inf 2 sloan\n\
                [-z,--scale]        choose image cut             \n\
                [-p,--path filename]  working directory         \n\
                [-a,--aperture] num   use fix aperture          \n\
                available list:                                 \n\
                "+str(bb)+"                                     \n\
################################################################"
if len(sys.argv)<=3:
    print help
    sys.exit()


logincl=src.open_program()
from pyraf import iraf

interactive = True
coordinatelist=''
out = ''
_system=''
_telescope=''
_aperture=''
subdirectory=['optical/','infrared/','sloan/']
scale=False

imglist=src.readlist(sys.argv[1])

options, args = getopt.getopt(sys.argv[2:],'l:,o,s:,p:,z,a:',['coordinatelist=','out=','system=','userPath=','scale','aperture='])
for opt,arg in options:
    if opt in ('-l', '--list'): coordinatelist = arg
    if opt in ('-o', '--output'): out = arg
    if opt in ('-s', '--system'): _system = int(arg)
    if opt in ('-p', '--path'): userPath = arg
    if opt in ('-z', '--scale'): scale = True
    if opt in ('-a', '--aperture'): _aperture = arg

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)
iraf.set(stdimage='imt2048')
iraf.tv.rimexam.backgrou = 'yes'
iraf.tv.rimexam.magzero = 0

if _aperture:
    iraf.tv.rimexam.iterati = 1
    iraf.tv.rimexam.radius = _aperture
else:
    iraf.tv.rimexam.iterati = 3
    iraf.tv.rimexam.radius = 7

if not coordinatelist:
    print "Please define coordinate list !!!! [svstandard.py filename  -l coordinatelist]"

####################################### check ##########################
check=0
_telescope=src.telescope(imglist[0])
print '### TELESCOPE= '+_telescope
if _system not in [1,2,0]:
     system=src.check_system(_telescope,imglist[0],Stdout=True)
else:
    system=_system

check=src.check_tel(imglist[0],_telescope,system)
if check==0:
        print '####################### '
        print "#### Error with the header !!!!"
        print '### telescope not in the list '
        print '### if you want to continue anyway,'
        print '### run "svother.py" to correct headers '
        print '####################### '

##########################################################################

mag,airmass,exptime,fwhm={},{},{},{}
airmass['U'],airmass['B'],airmass['V'],airmass['R'],airmass['I']=1,1,1,1,1
exptime['U'],exptime['B'],exptime['V'],exptime['R'],exptime['I']=1,1,1,1,1
airmass['u'],airmass['g'],airmass['r'],airmass['i'],airmass['z']=1,1,1,1,1
exptime['u'],exptime['g'],exptime['r'],exptime['i'],exptime['z']=1,1,1,1,1
airmass['J'],airmass['H'],airmass['K']=1,1,1
exptime['J'],exptime['H'],exptime['K']=1,1,1


src.delete(coordinatelist+".tv") 
src.delete(coordinatelist+"_templ.coo")
src.delete("_templ.*")
try:
    dir_system=subdirectory[system]
    iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','2,3,1', Stdout=coordinatelist+'.tv')
    iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1', Stdout='_templ.coo')
    iraf.wcsctran(coordinatelist+'.tv','_templ.coo2','home$coordinate_std/'+dir_system+coordinatelist+'_templ.fits',inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    a=iraf.fields('_templ.coo','1', Stdout=1)
    b=iraf.fields('_templ.coo2','1,2', Stdout=1)[2:]
    ff = open(coordinatelist+'_templ.coo','w')
    for i in range(len(a)):
        ff.write(b[i]+'\t'+a[i]+' \n')
    ff.close()
    if system==0 or system==2 : 
       standard=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6,7,8',Stdout=1)
    else:
       standard=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6',Stdout=1)
    #for i in standard:
    #    print standard[i]
except:
    print "Warning: no coordinate file found in "+str(iraf.show('home',Stdout=1)[0])+'coordinate_std/'+dir_system+'  !!! '
#    print "Warning: no coordinate file found in"+snoopy2.__path__[0]+'/coordinate_std/'+dir_system+'  !!! '
    src.close_program(logincl)
    
stars=[]
for i in standard:
    nome=string.split(i)
    stars.append(nome[0])
    
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
warn='##################################\n'
for img in imglist:
  if img:
    _header=src.read_parameter(_telescope,system)        
    src.delete("tmp*") 
    _filtro=src.filter(img,_header,_telescope)
    _filter=src.filtername(_telescope,_filtro,system)
    if _filter=='unknown':
        _filter=_filtro
        print '******************'
        print 'Warning: you are mixing the filter systems....the program will crash at the end !!!'
        print '******************'
    print '##########  '+str(_filter)+'  ##############'
    instrument=src.instrument(img,_header,_telescope) 
    _airmass=src.airmass(img,_header,_telescope)
    _exptime=src.exptime(img,_header,_telescope)
    print '#########'+str(_exptime)
    date=src.date(img,_header,_telescope)
    _ut=src.UT(img,_header,_telescope)
    xdim=src.xdimen(img,_telescope)
    ydim=src.ydimen(img,_telescope)
    airmass[_filter]=_airmass
    exptime[_filter]=_exptime
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    _z1,_z2='',''
    
    if interactive:
        print '######### Select FRAME TILE on your DS9 !!!!!'
        _z10,_z20,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+dir_system+coordinatelist+'_templ.fits',2,'','',False)
        if not goon: src.close_program(logincl)
                
        iraf.tvmark(2,coordinatelist+'_templ.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=3)

        _z1,_z2,goon=src.display_image(img,1,_z1,_z2,scale)
        if not goon: src.close_program(logincl)

        iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=3)
        answ = raw_input('Is the astrometry of the field good [y/n] ? [y] ')
        if not answ: answ='y'

        if answ=='n':
            try:
                src.delete('tmp.'+_filter+'.coo')
                src.delete('tmp.ccdb')
                iraf.ccmap('_first.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)
                iraf.wcsctran('_first_image.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
                iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=206,txsize=3)
                answ = raw_input('AND NOW, is the astrometry of the field good [y/n] ? [y] ')
                if not answ: answ='y'
            except: pass
        
        while answ=='n':
            _z1,_z2,goon=src.display_image(img,1,_z1,_z2,False)
            if not goon: src.close_program(logincl)

            print '>> Identify (minimum 2, preferably 3) reference stars (mark with "a")' 
            iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
            iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=3)
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
            _z1,_z2,goon=src.display_image(img,1,_z1,_z2,False)
            if not goon: src.close_program(logincl)

            src.delete('tmp.'+_filter+'.coo')
            iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
            iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=2)
            src.delete("tmp.ccmap")
            src.delete("tmp.coo")
            src.delete("tmp.ccdb")
            answ = raw_input('Is the astrometry of the field good  [y/n] ? [y]')
            if not answ: answ='y'

    src.delete("tmp.star")
    iraf.ccfind('home$coordinate_std/'+dir_system+coordinatelist+'.list','tmp.star',img,lngcolu=2,latcolu=3,lngunit='degrees',usewcs='yes')
    iraf.ccmap('tmp.star','tmp.ccdb',images=img,fitgeome='rscale',xcolum=9, ycolum=10,lngcolum=2,latcolumn=3,lngunit='degrees',update='yes',interact=False)
    src.delete('tmp.'+_filter+'.coo')
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    src.delete("tmp.ccdb,tmp.star")

    ff = open('tmp.'+_filter+'.coo','r')
    for j in range(0,3):  exl=ff.readline()
    rr = ff.readlines()
    ff.close()
    ff = open('tmp.imexam','w')
    for i in range(len(rr)):
        xx=string.split(rr[i])
        #print xx
        if not 1.<= float(xx[0])<=xdim or not 1.<= float(xx[1]) <=ydim:
            warn= warn +'##### WARNING: standard '+str(i)+', filter '+str(_filter[-1:])+' out of field !!!!!\n'
        ff.write(xx[0]+'\t'+xx[1]+'\ta\n')
    ff.close()
    src.delete("tmp.imex_output")
    ff = open('tmp.imexam','r')
    alllines = ff.readlines()
    ff.close()
    for i in alllines:
            xx=string.split(i)
            ff = open('tmp.one','w')
            ff.write(i)
            ff.close()            
            if not 1.<= float(xx[0])<=xdim or not 1.<= float(xx[1]) <=ydim:
                if not os.path.isfile('tmp.imex_output'):
                    os.system("echo '# [1] "+str(img)+" - "+str(coordinatelist)+"' > tmp.imex_output")
                    os.system("echo '#   COL    LINE   COORDINATES      R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT' >> tmp.imex_output")
                    os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")
                else:
                    os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")
            else:
                try:
                    iraf.imexam(input=img, frame=1, logfile='tmp.imex_output', keeplog='yes', imagecur='tmp.one', wcs='logical', use_disp='no')
                except:
                    if not os.path.isfile('tmp.imex_output'):
                        os.system("echo '# [1] "+str(img)+" - "+str(coordinatelist)+"\n' > tmp.imex_output")
                        os.system("echo '#   COL    LINE   COORDINATES      R    MAG    FLUX     SKY    PEAK    E   PA BETA ENCLOSED   MOFFAT DIRECT' >> tmp.imex_output")
                        os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")
                    else:
                        os.system("echo '999.  999. 999. 999.  INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF INDEF' >> tmp.imex_output")
    _fwhm = iraf.fields('tmp.imex_output','14',Stdout=1)
    for i in range(len(standard)):
        try: _fwhm[i]=float(_fwhm[i])
        except: _fwhm[i]=9999
	
    print '################### '+str(len(_fwhm))+' '+str(len(standard))
    fwhm[_filter]=_fwhm
    #iraf.hedit(img,'Seeing',mean(compress((array(_fwhm)<999),_fwhm)),add='yes',update='yes',verify='no')
    src.updateheader(img,0,'Seeing',mean(compress((array(_fwhm)<999),_fwhm)))
##############################    prompt correction
    if _telescope=='prompt':
        try: prompttele = pyfits.getheader(img)['OBSERVAT']
	except: prompttele=''        
        if string.count(prompttele,'2'): _prompt='prompt2'
        elif string.count(prompttele,'3'): _prompt='prompt3'
        elif string.count(prompttele,'4'): _prompt='prompt4'
        elif string.count(prompttele,'5'): _prompt='prompt5'
        else: print 'WARNING: telescope not found !! '
        print 'Telescope prompt: Magnitude Correction !!!!!'
        _trimsection=src.trimsec(img,_header,_telescope)
        if _trimsection:
            xxx=string.split(string.split(_trimsection[1:-1],',')[0],':')[0]
            yyy=string.split(string.split(_trimsection[1:-1],',')[1],':')[0]
        else:
            xxx,yyy=0,0
        _coordinate_star = iraf.fields('tmp.imex_output','3,4',Stdout=1)
        _coord_origi_x,_coord_origi_y,_values=[],[],[]
        try:
            _z10,_z20,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+_telescope+'/'+_prompt+'.fits',2,'','',False)
        except:
            pass
        for ii in _coordinate_star:
            xxn,yyn=string.split(ii)
            _coord_origi_x.append(int(float(xxx)+float(xxn)))
            _coord_origi_y.append(int(float(yyy)+float(yyn)))
            _xx=int(float(xxx)+float(xxn))
            _yy=int(float(yyy)+float(yyn))
#            _values.append(0.0)
            try:
                _values.append(string.split(iraf.imstat('home$coordinate_std/'+_telescope+'/'+_prompt+'.fits['+str(_xx)+','+str(_yy)+']',Stdout=1)[1])[2])
            except:
                _values.append(0.0)

################################################

    _mag = iraf.fields('tmp.imex_output','6',Stdout=1)
    for i in range(len(standard)):
        if _telescope=='prompt':
            try: _mag[i]=float(_mag[i])-float(_values[i])
            except: _mag[i]=9999
        else:
            try: _mag[i]=float(_mag[i])
            except: _mag[i]=9999
        if _mag[i]>=9999: warn= warn +'##### WARNING: standard '+str(i)+', filter '+str(_filter[-1:])+' INDEF !!!!!\n'
        if 2 >=_mag[i]>= -2: warn= warn +'##### WARNING: standard '+str(i)+', filter '+str(_filter[-1:])+' very faint !!!!!\n'
    mag[_filter]=_mag
  else:
      	print '####'
  	print '#### WARNING: empty space in the list !!'
	print '####'

if system==0:    
  try:
    mag['U']
    summaryU='band U -> ok\n'
  except:
    mag['U']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryU='band U -> none\n'
  try:
    mag['B']
    summaryB='band B -> ok\n'
  except:
    mag['B']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryB='band B -> none\n'
  try:
    mag['V']
    summaryV='band V -> ok\n'
  except:
    mag['V']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryV='band V -> none\n'
  try:
    mag['R']
    summaryR='band R -> ok\n'
  except:
    mag['R']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryR='band R -> none\n'
  try:
    mag['I']
    summaryI='band I -> ok\n'
  except:
    mag['I']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryI='band I -> none\n'
elif system==2:
  try:
    mag['u']
    summaryU='band u -> ok\n'
  except:
    mag['u']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryU='band u -> none\n'
  try:
    mag['g']
    summaryB='band g -> ok\n'
  except:
    mag['g']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryB='band g -> none\n'
  try:
    mag['r']
    summaryV='band r-> ok\n'
  except:
    mag['r']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryV='band r -> none\n'
  try:
    mag['i']
    summaryR='band i -> ok\n'
  except:
    mag['i']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryR='band i -> none\n'
  try:
    mag['z']
    summaryI='band z -> ok\n'
  except:
    mag['z']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryI='band z -> none\n'
elif system==1:
  try:
    mag['J']
    summaryU='band J -> ok\n'
  except:
    mag['J']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryU='band J -> none\n'
  try:
    mag['H']
    summaryB='band H -> ok\n'
  except:
    mag['H']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryB='band H -> none\n'
  try:
    mag['K']
    summaryV='band K-> ok\n'
  except:
    mag['K']=list(array(_fwhm)-array(_fwhm)+9999)
    summaryV='band K -> none\n'


if not out: out = _telescope+'_'+str(coordinatelist)+'_'+str(date)+'_'+str('0'+string.split(_ut,':')[0])[-2:]+'_'+str(system)


src.delete(out)


fil = open(out,'w')
fil.write(str(instrument)+' '+str(date)+'_ime\n')
fil.write('*** '+coordinatelist+' '+str(len(standard))+'\n')
if system==0:
  fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(exptime['U']),str(exptime['B']),str(exptime['V']),str(exptime['R']),str(exptime['I'])))
  fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['U']),str(airmass['B']),str(airmass['V']),str(airmass['R']),str(airmass['I'])))
  for i in range(len(standard)):
    fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                 % (str(mag['U'][i]),str(mag['B'][i]),str(mag['V'][i]),str(mag['R'][i]),str(mag['I'][i]),str(standard[i])))
elif system==2:
  fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(exptime['u']),str(exptime['g']),str(exptime['r']),str(exptime['i']),str(exptime['z'])))
  fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['u']),str(airmass['g']),str(airmass['r']),str(airmass['i']),str(airmass['z'])))
  for i in range(len(standard)):
    fil.write('%6.6s\t%6.6s\t%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                 % (str(mag['u'][i]),str(mag['g'][i]),str(mag['r'][i]),str(mag['i'][i]),str(mag['z'][i]),str(standard[i])))
elif system==1:
  fil.write('%6.6s\t%6.6s\t%6.6s\n' % (str(exptime['J']),str(exptime['H']),str(exptime['K'])))
  fil.write('%6.6s\t%6.6s\t%6.6s\n' % (str(airmass['J']),str(airmass['H']),str(airmass['K'])))
  for i in range(len(standard)):
    fil.write('%6.6s\t%6.6s\t%6.6s\t%60.60s\n' \
                 % (str(mag['J'][i]),str(mag['H'][i]),str(mag['K'][i]),str(standard[i])))
fil.close()
if system==0:
   print '##############################\n'+summaryU+summaryB+summaryV+summaryR+summaryI+warn+'##############################\n'
elif system==1:
   print '##############################\n'+summaryU+summaryB+summaryV+warn+'##############################\n'
elif system==2:
   print '##############################\n'+summaryU+summaryB+summaryV+summaryR+summaryI+warn+'##############################\n'

src.delete(coordinatelist+".tv")
src.delete(coordinatelist+"_templ.coo")
src.delete("tmp*,_templ.co*,_first_image.tv,_first.ccmap")

if _aperture:
    iraf.tv.rimexam.iterati = 3
    iraf.tv.rimexam.radius = 7

src.close_program(logincl)
