#!/usr/bin/env python

from snoopy2 import *
import getopt
import os,sys,shutil,string,re
import snoopy2
import glob
userPath=os.getcwd()
from astropy.io import fits as pyfits
from astropy import wcs as pywcs
#import pywcs

bb=[]
aa=glob.glob(snoopy2.__path__[0]+'/coordinate_std/optical/*list')
for i in aa:
    bb.append(i[len(snoopy2.__path__[0]+'/coordinate_std/optical/'):-5])
    
help ="################################################################ \n\
help = Usage:   svstandard.py filename  -l coordinatelist           \n\
                input              filelist iraf format             \n\
                -l,--list   file coordinate                         \n\
                [-o,--output filename] output file                  \n\
                [-t,--transform]   consider geometric transformation\n\
                [-r,--rejection]   type of rejection (none,minmax,  \n\
                                         ccdclip,sigclip,avsigclip)  \n\
                [-c,--combine] type of combination (average,median,sum)\n\
                [-z,--zero region]    region to scale the sky \n\
                [-f, --factor]  size factor  (working with transform) \n\
                available list:                                     \n\
                "+str(bb)+"                                         \n\
################################################################"
if len(sys.argv)<=3:
    print help
    sys.exit()


def enlargeimage(img,outputimage,_system,factor):
    import snoopy2
    from snoopy2 import src
    _telescope=src.telescope(img)
    _header=src.read_parameter(_telescope,_system)
    _xdimen=src.xdimen(img,_telescope)
    _ydimen=src.ydimen(img,_telescope)
    S=int(_xdimen/factor)
    xtot=S+int(_xdimen)+S
    ytot=S+int(_ydimen)+S
    xstart=S
    xstop=S+int(_xdimen)
    ystart=S
    ystop=S+int(_ydimen)
    from iraf import ctio
    src.delete('testim.fits')
    iraf.imcreate('testim.fits',naxis=2,naxis1=xtot,naxis2=ytot,header='copy',referen=img)
    mm = iraf.imstat(img,fields='mean',nclip=3,lsigma=5,usigma=2,Stdout=1)[1]
    src.delete('testim2.fits')
    iraf.imarith(operand1='testim.fits', op='+', operand2=mm,result='testim2.fits')
    src.delete(outputimage)
    print 'testim2.fits['+str(xstart)+':'+str(xstop)+','+str(ystart)+':'+str(ystop)+']'
    iraf.imcopy(input=img, output='testim2.fits['+str(xstart)+':'+str(xstop)+','+str(ystart)+':'+str(ystop)+']')
    iraf.imrename('testim2.fits',outputimage)
    return outputimage

logincl=src.open_program()
    
from pyraf import iraf

iraf.astcat(_doprint=0)
iraf.imcoords(_doprint=0)
iraf.set(stdimage='imt1024')
iraf.tv.rimexam.backgrou = 'yes'

_rejection='none'
interactive = False
coordinatelist=''
out = ''
_system=''
_zero=False
_combine='average'
_telescope=''
subdirectory=['optical/','infrared/','sloan/']
_transform=False
_factor=''

imglist=src.readlist(sys.argv[1])

options, args = getopt.getopt(sys.argv[2:],'l:,o:,s:,i,t,r:,c:,z,f:',['coordinatelist=','output=','system=','interactive=','transform=','rejection=','combine=','zero','factor='])
for opt,arg in options:
    if opt in ('-l', '--list'): coordinatelist = arg
    if opt in ('-i', '--interactive'): interactive = True
    if opt in ('-t', '--transform'): _transform = True
    if opt in ('-o', '--output'): out = arg
    if opt in ('-r', '--rejection'): _rejection = arg
    if opt in ('-c', '--combine'): _combine = arg
    if opt in ('-s', '--system'): _system = int(arg)
    if opt in ('-z', '--zero'): _zero = True
    if opt in ('-f', '--factor'): _factor = int(arg)

if not coordinatelist:
    print "Please define coordinate list !!!! [svstandard.py filename  -l coordinatelist]"

####################################### check ##########################

check=0
_telescope=src.telescope(imglist[0])
print '### TELESCOPE= '+_telescope
if _system not in [0,1,2]:
     _system=src.check_system(_telescope,imglist[0],Stdout=True)
else:
    print '### system = '+str(_system)

check=src.check_tel(imglist[0],_telescope,_system)
if check==0:
        print '####################################### '
        print "#### Error with the header !!!!"
        print '### telescope not in the list '
        print '### if you want to continue anyway,'
        print '### run "svother.py" to correct header '
        print '######################################## '

##########################################################################

src.delete(coordinatelist+".tv") 
src.delete(coordinatelist+"_templ.coo")
src.delete("_templ.*")
try:
    dir_system=subdirectory[_system]
    iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','2,3,1', Stdout=coordinatelist+'.tv')
    iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1', Stdout='_templ.coo')
    iraf.wcsctran(coordinatelist+'.tv','_templ.coo2','home$coordinate_std/'+dir_system+coordinatelist+'_templ.fits',inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    a=iraf.fields('_templ.coo','1', Stdout=1)
    b=iraf.fields('_templ.coo2','1,2', Stdout=1)[2:]
    ff = open(coordinatelist+'_templ.coo','w')
    for i in range(len(a)):
        ff.write(b[i]+'\t'+a[i]+' \n')
    ff.close()
    if _system==0 or _system==2 : 
       standard=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6,7,8',Stdout=1)
    else:
       standard=iraf.fields('home$coordinate_std/'+dir_system+coordinatelist+'.list','1,4,5,6',Stdout=1)

except:
    print "WARNING: no coordinate file found in "+snoopy2.__path__[0]+'/coordinate_std/'+dir_system+'  !!! '
#    src.close_program(logincl)
    
stars=[]
for i in standard:
    nome=string.split(i)
    stars.append(nome[0])
    
print imglist
print len(imglist)

if _transform and _factor:
    enlargeimage(imglist[0],'pippo.fits',_system,_factor)
    imglist.append('pippo.fits')

iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
warn='##################################\n'
filtri=[]
filerejected=[]
for img in imglist:
  print img,'#######'
  if img:
#    if not system:
#        system=src.check_system(_telescope,img,Stdout=False)
    if _system=='99':
        print 'Warning: system not known '
    _header=src.read_parameter(_telescope,_system)
    _gain=src.gain(img,_header,_telescope)
    _ron=src.ron(img,_header,_telescope)
    if not _gain: _gain=1
    if not _ron: _ron=1
    src.delete('tmp2.coo')
    src.delete("tmp*") 
    _filtro=src.filter(img,_header,_telescope)
    _filter=src.filtername(_telescope,_filtro,_system)
    if _filter=='unknown':
        _filter=_filtro
    if _filter not in filtri:
        filtri.append(_filter)
    print '##########  '+str(_filter)+'  ##############'
    instrument=src.instrument(img,_header,_telescope) 
    _date=src.date(img,_header,_telescope)
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')

    if interactive:
        print '######### Select FRAME TILE on your DS9 !!!!!'
        _z1,_z2,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+dir_system+coordinatelist+'_templ.fits',2,'','',False)
        if not goon: src.close_program(logincl)
                
        iraf.tvmark(2,coordinatelist+'_templ.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=4)

        _z1,_z2,goon=src.display_image(img,1,'','',False)
        if not goon: src.close_program(logincl)
        
        iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=2)
        answ = raw_input('is the astrometry of the field good ? \n  <d> to cancel the image from the list \n [y/n/d] ? [y] ')
        if not answ: answ='y'
        
        if answ=='d':  filerejected.append('d')
        else: filerejected.append('a')
        if answ=='n':
            try:
                src.delete('tmp.'+_filter+'.coo')
                src.delete('tmp.ccdb')
                iraf.ccmap('_first.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)
                iraf.wcsctran('_first_image.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
                iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='no',label='yes',radii=5,nxoffse=5,nyoffse=5,color=206,txsize=4)
                answ = raw_input('AND NOW, is the astrometry of the field good [y/n] ? [y] ')
                if not answ: answ='y'
            except: pass
        
        while answ=='n':

            _z1,_z2,goon=src.display_image(img,1,'','',False)
            if not goon: src.close_program(logincl)

            print '>> Identify (minimum 2, preferably 3) reference stars (mark with "a")' 
            iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
            iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=4)
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
            
            _z1,_z2,goon=src.display_image(img,1,'','',False)
            if not goon: src.close_program(logincl)
            
            src.delete('tmp.'+_filter+'.coo')
            iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
            iraf.tvmark(1,'tmp.'+_filter+'.coo',mark="circle",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=205,txsize=4)
            src.delete("tmp.ccmap")
            src.delete("tmp.coo")
            src.delete("tmp.ccdb")
            answ = raw_input('is the astrometry of the field good  [y/n] ? [y]')
            if not answ: answ='y'
    
 
    src.delete("tmp.star")
    iraf.ccfind('home$coordinate_std/'+dir_system+coordinatelist+'.list','tmp.star',img,lngcolu=2,latcolu=3,lngunit='degrees',usewcs='yes')
    iraf.ccmap('tmp.star','tmp.ccdb',images=img,fitgeome='rscale',xcolum=9, ycolum=10,lngcolum=2,latcolumn=3,lngunit='degrees',update='yes',interact=False)
    src.delete('tmp.'+_filter+'.coo')
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+_filter+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    src.delete("tmp.ccdb")
    src.delete("tmp.star")



imglistnew=[]
if filerejected:
    for i in range(len(filerejected)):
        if filerejected[i]!='d':
            imglistnew.append(imglist[i])
else:
    imglistnew=imglist

if _transform and _factor:
    imglistnew.remove('pippo.fits')

for fil in filtri:
    out2 = 'merge_'+fil+'.fits'
    tmplist=[]
    for img in imglistnew:
        _filtro=src.filter(img,_header,_telescope)
        _filter=src.filtername(_telescope,_filtro,_system)
        if _filter=='unknown':
            _filter=_filtro
        if _filter==fil:
            tmplist.append(img)
            system0=src.check_system(_telescope,img,Stdout=False)
            print img

    src.delete("listtmp")
    if len(tmplist)>=2:
        dd=open('listtmp','w')
        for i in tmplist:
            dd.write(i+'\n')
        dd.close()
        src.delete(out2)
        if _transform==False:
            if _zero:
                src.delete('tmp_imcomb.fits')
                iraf.imcombine('@listtmp','tmp_imcomb.fits',combine='sum', reject='none', offsets='wcs',rdnoise=_ron, gain=_gain)
                _z1,_z2,goon=src.display_image('tmp_imcomb.fits',1,'','',False)
                ########################
                print 'use sky to scale image'
                print '>> Identify the region to measure the background with >x< and >q< for quit ' 
                src.delete('tmp.coo')
                src.delete('tmp2.coo')
                _z1,_z2,goon=src.display_image('tmp_imcomb.fits',1,'','',True)
                iraf.imexamine('tmp_imcomb.fits',1,logfile='tmp.coo',keeplog='yes')
                iraf.tvmark(1,'tmp.coo',mark="cross",number='yes',label='no',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=4)
                xycoo = iraf.fields('tmp.coo','1,2,3',Stdout=1)
                
                x0,y0,z0=string.split(xycoo[0])
                hdulist = pyfits.open('tmp_imcomb.fits')
                wcs = pywcs.WCS(hdulist[0].header)
                pix = wcs.wcs_pix2sky(float(x0),float(y0), 1)
                ff=open('statfile','w')
                m=0
                for i in range(0,len(tmplist)):
                    hdulist2 = pyfits.open(tmplist[i])
                    wcs2 = pywcs.WCS(hdulist2[0].header)
                    xx,yy=wcs2.wcs_sky2pix(pix[0],pix[1], 1)
                    xx,yy=int(xx[0]),int(yy[0])
                    _statset='['+str(xx-30)+':'+str(xx+30)+','+str(yy-30)+':'+str(yy+30)+']'
                    print _statset
                    print 'here2'
                    #val=string.split(iraf.imstat(tmplist[m]+_statset,Stdout=1)[-1])[2]
                    val = iraf.imstat(img+_statset,fields='mean',nclip=3,lsigma=5,usigma=2,Stdout=1)
                    if i==0: val0=float(val[1])
                    print val
                    ff.write(str(float(val[1])/val0)+'\n')
                    m=m+1
                ff.close()
                iraf.imcombine('@listtmp', out2 , combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain, lthresh=-9000, statsec='', scale='@statfile', \
                 zero='none', weight='none', nlow=1, nhigh=1,logfile='')

#                src.delete('tmp.tv')
#                src.delete('tmp2.coo')
#                iraf.wcsctran('tmp.coo','tmp.tv','tmp_imcomb.fits',inwcs='logical',units='degrees degrees',outwcs='world',columns='1 2',formats='%10.1f %10.1f')
#                iraf.wcsctran('tmp.tv','tmp2.coo','@listtmp',inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%6.2f %6.2f')
#                xycoo2 = iraf.fields('tmp2.coo','1,2,3',Stdout=1)
#                ff=open('statfile','w')
#                m=0
#                for i in range(0,len(xycoo2)):
#                    if xycoo2[i]:
#                        x0=int(float(string.split(xycoo2[i])[0]))
#                        y0=int(float(string.split(xycoo2[i])[1]))
#                        z0=int(float(string.split(xycoo2[i])[2]))
#                        print x0,y0
#                        _statset='['+str(x0-30)+':'+str(x0+30)+','+str(y0-30)+':'+str(y0+30)+']'
#                        print _statset
#                        print 'here2'
#                        val=string.split(iraf.imstat(tmplist[m]+_statset,Stdout=1)[-1])[2]
#                        ff.write(val+'\n')
#                        m=m+1
#                ff.close()
#                iraf.imcombine('@listtmp',out2,combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain, lthresh=-10, statsec='', scale='@statfile', zero='none')
#                ########################
#                iraf.imcombine('@listtmp',out2,combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain)
            else:
                iraf.imcombine('@listtmp',out2,combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain)
        else:
            dd=open('listtmp2','w')
            for i in tmplist:
                dd.write('T'+i+'\n')
                src.delete('T'+i)
            dd.close()
            
#            iraf.sregister('@listtmp',tmplist[0],'@listtmp2',wcs='world', boundar='constant', constant=-10000)
            if _factor:
                iraf.sregister('@listtmp','pippo.fits','@listtmp2',wcs='world', boundar='constant', constant=-10000)
            else:
                iraf.sregister('@listtmp',tmplist[0],'@listtmp2',wcs='world', boundar='constant', constant=-10000)
            if _zero:
                src.delete('tmp_imcomb.fits')
                iraf.imcombine('@listtmp2','tmp_imcomb.fits',combine='sum', reject='none', offsets='wcs',rdnoise=_ron, gain=_gain)
                _z1,_z2,goon=src.display_image('tmp_imcomb.fits',1,'','',False)
                ########################
                print 'use sky to scale image'
                print '>> Identify the region to measure the background with >x< and >q< for quit ' 
                src.delete('tmp.coo')
                _z1,_z2,goon=src.display_image('tmp_imcomb.fits',1,'','',True)
                iraf.imexamine('tmp_imcomb.fits',1,logfile='tmp.coo',keeplog='yes')
                iraf.tvmark(1,'tmp.coo',mark="cross",number='no',label='no',radii=5,nxoffse=5,nyoffse=5,color=214,txsize=4)
                xycoo = iraf.fields('tmp.coo','1,2,3',Stdout=1)
                print xycoo
                print 'here'
                x0=int(float(string.split(xycoo[0])[0]))
                y0=int(float(string.split(xycoo[0])[1]))
                _statsec='['+str(x0-30)+':'+str(x0+30)+','+str(y0-30)+':'+str(y0+30)+']'
                iraf.imcombine('@listtmp2',out2,combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain, lthresh=-9000, statsec=_statsec, zero= 'median')
                raw_input('go on ? ')
                ########################
            else:
                iraf.imcombine('@listtmp2',out2,combine=_combine, reject=_rejection, offsets='wcs',rdnoise=_ron, gain=_gain, lthresh=-9000)
            src.delete('listtmp*')
            #for i in tmplist:
            #    src.delete('T'+i)

    elif len(tmplist)==1:
        src.delete(out2)
        iraf.imcopy(tmplist[0],out2)
        src.delete('listtmp*')

    if glob.glob(out2):
        _date=src.date(out2,_header,_telescope)
        _ut=src.UT(out2,_header,_telescope)
        if not out:
            out=src.objects(out2,_header,_telescope)
            out=re.sub('-','',out)
            if not out:
                out='SN'
        src.delete(out+'_'+_date+'_'+str('0'+string.split(_ut,':')[0])[-2:]+'_'+fil+str(system0)+'.fits')
        iraf.imrename(out2,out+'_'+_date+'_'+str('0'+string.split(_ut,':')[0])[-2:]+'_'+fil+str(system0)+'.fits',verbose='yes')

src.delete("tmp*")
src.delete("_templ.co*")
src.delete(coordinatelist+"_templ.coo")
src.delete(coordinatelist+".tv")

src.close_program(logincl)




#x0,y0,z0=string.split(xycoo[0])
#hdulist = pyfits.open('tmp_imcomb.fits')
#wcs = pywcs.WCS(hdulist[0].header)
#pix = wcs.wcs_pix2sky(x0,y0, 1)

#hdulist2 = pyfits.open(tmplist[i])
#wcs2 = pywcs.WCS(hdulist2[0].header)
#xx,yy=wcs2.wcs_sky2pix(pix[0],pix[1], 1)
#xx,yy=xx[0],yy[0]
