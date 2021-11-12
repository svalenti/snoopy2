def register_module(img,_system,coordinatelist,interactive,logincl,_filter):
    from snoopy2 import src
    import snoopy2
    import string,os

#    logincl=src.open_program()
    from pyraf import iraf
    iraf.astcat(_doprint=0)
    iraf.imcoords(_doprint=0)
    iraf.set(stdimage='imt1024')
    iraf.tv.rimexam.backgrou = 'yes'

    subdirectory=['optical/','infrared/','sloan/']
    iraf.delete(coordinatelist+".tv",verify='no') 
    iraf.delete(coordinatelist+"_templ.coo",verify='no')
    iraf.delete("_templ.*",verify='no')
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
    except Exception as e:
        print e
        print "WARNING: no coordinate  " + coordinatelist + " file found in " + snoopy2.__path__[0] + '/coordinate_std/'+dir_system+'  !!! '
        src.close_program(logincl)
    stars=[]
    for i in standard:
        nome=string.split(i)
        stars.append(nome[0])

#    _filter='XXX'
    src.delete("tmp."+img+".coo") 
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+img+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')

    if interactive:
        print '######### Select FRAME TILE on your DS9 !!!!!'
        _z1,_z2,goon=src.display_image(snoopy2.__path__[0]+'/coordinate_std/'+dir_system+coordinatelist+'_templ.fits',2,'','',False)
        if not goon: src.close_program(logincl)
                
        iraf.tvmark(2,coordinatelist+'_templ.coo',mark="circle",number='no',label='yes',radii=20,nxoffse=15,nyoffse=15,color=214,txsize=4)

        _z1,_z2,goon=src.display_image(img,1,'','',False)
        if not goon: src.close_program(logincl)
        
        iraf.tvmark(1,'tmp.'+img+'.coo',mark="circle",number='no',label='yes',radii=20,nxoffse=15,nyoffse=15,color=205,txsize=2)
        answ = raw_input('is the astrometry of the field good ? \n [y/n] ? [y] ')
        if not answ: answ='y'

        if answ=='n':
            try:
                src.delete('tmp.'+img+'.coo')
                src.delete('tmp.ccdb')
                iraf.ccmap('_first.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)
                iraf.wcsctran('_first_image.tv','tmp.'+img+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
                iraf.tvmark(1,'tmp.'+img+'.coo',mark="circle",number='no',label='yes',radii=20,nxoffse=15,nyoffse=15,color=206,txsize=4)
                answ = raw_input('AND NOW, is the astrometry of the field good [y/n] ? [y] ')
                if not answ: answ='y'
            except: pass
        
        while answ=='n':

            _z1,_z2,goon=src.display_image(img,1,'','',False)
            if not goon: src.close_program(logincl)

            print '>> Identify (minimum 2, preferably 3) reference stars (mark with "a")' 
            iraf.delete('tmp.coo',verify='no')
            iraf.imexamine(img,1,logfile='tmp.coo',keeplog='yes')
            iraf.tvmark(1,'tmp.coo',mark="circle",number='yes',label='no',radii=15,nxoffse=15,nyoffse=15,color=214,txsize=4)
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
            iraf.delete('_first_image.tv',verify='no')
            iraf.delete('_first.ccmap',verify='no')
            os.system('cp '+coordinatelist+'.tv _first_image.tv')
            os.system('cp tmp.ccmap  _first.ccmap')
            iraf.ccmap('tmp.ccmap','tmp.ccdb',images=img,fitgeome='rscale',lngunit='degrees',update='yes',interact=False)            
            _z1,_z2,goon=src.display_image(img,1,'','',False)
            if not goon: src.close_program(logincl)
            
            iraf.delete('tmp.'+img+'.coo',verify='no')
            iraf.wcsctran(coordinatelist+'.tv','tmp.'+img+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
            iraf.tvmark(1,'tmp.'+img+'.coo',mark="circle",number='yes',label='no',radii=20,nxoffse=15,nyoffse=15,color=205,txsize=4)
            iraf.delete("tmp.ccmap",verify='no')
            iraf.delete("tmp.coo",verify='no')
            iraf.delete("tmp.ccdb",verify='no')
            answ = raw_input('is the astrometry of the field good  [y/n] ? [y]')
            if not answ: answ='y'
     
    iraf.delete("tmp.star",verify='no')
    iraf.ccfind('home$coordinate_std/'+dir_system+coordinatelist+'.list','tmp.star',img,lngcolu=2,latcolu=3,lngunit='degrees',usewcs='yes')
    iraf.ccmap('tmp.star','tmp.ccdb',images=img,fitgeome='rscale',xcolum=9, ycolum=10,lngcolum=2,latcolumn=3,lngunit='degrees',update='yes',interact=False)
    iraf.delete('tmp.'+img+'.coo',verify='no')
    iraf.wcsctran(coordinatelist+'.tv','tmp.'+img+'.coo',img,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
    iraf.delete("tmp.ccdb",verify='no')
    iraf.delete("tmp.star",verify='no')
    iraf.delete("tmp.coo",verify='no')
    return stars,'tmp.'+img+'.coo'

#######################################################################################
##########################        DEFINE FITSN  #########################################

def fitsn(_recenter,img,imgpsf,fwhm0,apco0,z22,z11,midpt,size,nor,_values,DM):
    from pyraf import iraf
    import string,os,sys
    from numpy import log10

    a1 = int(fwhm0+.5)
    a2 = int(2.*fwhm0+.5)
    a3 = int(3.*fwhm0+.5)
    a4 = int(4.*fwhm0+.5)
    a5 = int(5.*fwhm0+.5)
    ap = str(a1)+","+str(a2)+","+str(a3)

    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)

    from iraf import digiphot
    from iraf import daophot
    from iraf import ptools

    zmag = 0 
    iraf.noao.digiphot.daophot.photpars.zmag = zmag

    iraf.delete("apori",verify='no')
    iraf.delete(img+".sn.mag",verify='no')
    iraf.noao.digiphot.daophot.phot("original",img+".sn.coo","apori",veri='no')   
    iraf.noao.digiphot.daophot.phot("sn",img+".sn.coo",img+".sn.mag",veri='no')   
    iraf.noao.digiphot.daophot.daopars.psfrad = a4
    iraf.noao.digiphot.daophot.daopars.fitrad = fwhm0  
    iraf.noao.digiphot.daophot.daopars.fitsky = 'no'
    iraf.noao.digiphot.daophot.daopars.sannulus = int(a4)

    if _recenter:
        answ = raw_input(">>> recentering for targets [yes/no] ? [yes] ")
        if not answ:answ='yes'
    else:
        answ='yes'
    iraf.noao.digiphot.daophot.daopars.recenter = answ
    iraf.noao.digiphot.daophot.daopars.fitsky = 'yes'
    iraf.delete(img+".sn.als",verify='no')
    iraf.allstar("sn",img+".sn.mag",imgpsf,img+".sn.als","","residual",veri='no',verb='no')
    iraf.delete("snfit.fits",verify='no')
    iraf.imarith("sn","-","residual","snfit")
    iraf.delete("skyfit.fits",verify='no')
    iraf.imarith("original","-","snfit","skyfit")
    iraf.txsort(img+".sn.als","ID")
    tmptbl = iraf.txdump(img+".sn.als","mag,merr,xcenter,ycenter",expr='yes', Stdout=1)
    magerr,fitmag,centx,centy=[],[],[],[]
    for i in tmptbl:
        try:
            fitmag.append(float(string.split(i)[0]))#+2.5*log10(nor))
        except:
            fitmag.append('INDEF')
        try:
            magerr.append(float(string.split(i)[1]))
        except:
            magerr.append('INDEF')

        centx.append(float(string.split(i)[2]))
        centy.append(float(string.split(i)[3]))
    tmptbl=iraf.txdump("apori","mag",expr='yes', Stdout=1)
    apori1,apori2,apori3=[],[],[]
    for i in tmptbl:
        try:
            apori1.append(float(string.split(i)[0])-float(_values)-float(DM))#+2.5*log10(nor)
        except:
            apori1.append('INDEF')
        try:            
            apori2.append(float(string.split(i)[1])-float(_values)-float(DM))#+2.5*log10(nor)
        except:
            apori2.append('INDEF')
        try:
            apori3.append(float(string.split(i)[2])-float(_values)-float(DM))#+2.5*log10(nor)
        except:
            apori3.append('INDEF')

           
    iraf.txsort(img+".sn.mag","YCENTER")
    tmptbl=iraf.txdump(img+".sn.mag","mag,magerr",expr='yes', Stdout=1) 

    print "********************************************************************"
    print "ID <apmag on original>  <apmag on bgsubt> fitmag truemag err_fit"         
    print "     ",a1,"       ",a2,"      ",a3,"        ",a1,"     ",a2,"     ",a3 


    apmag1,apmag2,apmag3,truemag=[],[],[],[]
    for i in range(len(tmptbl)):
        try:
            apmag1.append(float(string.split(tmptbl[i])[0])-float(_values)-float(DM))#+2.5*log10(nor)
        except:
            apmag1.append('INDEF')
        try:
            apmag2.append(float(string.split(tmptbl[i])[1])-float(_values)-float(DM))#+2.5*log10(nor)
        except:
            apmag2.append('INDEF')
        try:
            apmag3.append(float(string.split(tmptbl[i])[2])-float(_values)-float(DM))#+2.5*log10(no)
        except:
            apmag3.append('INDEF')
        try:
            truemag.append(float(fitmag[i])+float(apco0)+2.5*log10(nor)-float(_values)-float(DM))
        except:
            truemag.append('INDEF')
        print i,apori1[i],apori2[i],apori3[i],apmag1[i],apmag2[i],apmag3[i],fitmag[i],truemag[i],magerr[i]

    print "********************************************************************"
    
    iraf.display("original",1,fill='yes',xcen=.25,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11) 

    z01 = z11-midpt
    z02 = z22-midpt 
    s1 = 1
    s2 = -int(fwhm0)
    iraf.delete("tmptbl",ve='no')
    ff=open('tmptbl','w')
    ff.write(str(s1)+' '+str(s2)+" ORIGINAL")
    ff.close()    
    iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
    
    iraf.display("snfit",1,erase='no',fill='yes',xcen=.25,ycen=.75,xsize=.3,ysize=.3, zscal='no',zrang='no', z2=z02, z1=z01)
    iraf.delete("tmptbl",ve='no')
    tmptbl0=iraf.txdump(img+".sn.als","xcen,ycen",expr='yes',Stdout=1)
    ff=open('tmptbl','w')
    for i in tmptbl0:
        ff.write(i+'\n')
    ff.close()    
    lra = int((2*float(size)*float(fwhm0))*2)
    iraf.tvmark(1,"tmptbl",autol='no',mark="circle", number='yes',nyoffset=lra,radi=a2,txsize=2,inter='no')
    s1 = 1
    s2 = -1*int(fwhm0)
    iraf.delete("tmptbl",ve='no')
    ff=open('tmptbl','w')
    ff.write(str(s1)+' '+str(s2)+" FITTED")
    ff.close()    
    iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
    
    iraf.display("skyfit",1,erase='no',fill='yes',xcen=.75,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11)
    s1 = 1
    s2 = -1*int(fwhm0)
    iraf.delete("tmptbl",ve='no')
    ff=open('tmptbl','w')
    ff.write(str(s1)+' '+str(s2)+" RESIDUAL")
    ff.close()    
    iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
    print '###'
    print '### magnitude scaled to target exposure time using: '+str(2.5*log10(nor))
    print '### fit magnitude corrected with aperture correction: '+str(float(apco0))
    print '### magnitude scaled to target using DM: '+str(DM)
    print '###'
    return apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy

###########################  DEFINE MANUAL CHANGING  #########################################

def manu(dmag0,apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,z22,z11,midpt,size,fwhm0,img,x1,y1,arterr):
   from pyraf import iraf
   import string,os,sys
   from numpy import array 

   a1 = int(fwhm0+.5)
   a2 = int(2.*fwhm0+.5)
   a3 = int(3.*fwhm0+.5)
   a4 = int(4.*fwhm0+.5)
   a5 = int(5.*fwhm0+.5)
   ap = str(a1)+","+str(a2)+","+str(a3)

   fdmag = 10**(-0.4*float(dmag0))
   iraf.delete(",_snfit.fit?,skyfit.fit?",ve='no')
   iraf.imarith("snfit","*",fdmag,"_snfit")   
   iraf.imarith("original","-","_snfit","skyfit")
   iraf.display("original",1,fill='yes',xcen=.25,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1= z11) 
   z01 = z11-midpt
   z02 = z22-midpt 
   s1 = 1
   s2 = -int(fwhm0)
   iraf.delete("tmptbl",ve='no')
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" ORIGINAL") 
   ff.close()
   iraf.tvmark(1,"tmptbl",autol='yes',mark="none",inter='no',label='yes',txsize=2)

   iraf.display("_snfit",1,erase='no',fill='yes',xcen=.25,ycen=.75,xsize=.3,ysize=.3, zscal='no', zrang='no', z2=z02, z1=z01 )
   iraf.delete("tmptbl",ve='no')
   tmptbl=iraf.txdump(img+".sn.als","xcen,ycen",expr='yes', Stdout=1)
   ff=open('tmptbl','w')
   for i in tmptbl:
       ff.write(i) 
   ff.close()  
   lra = int((2*size*fwhm0)*2)
   iraf.tvmark(1,"tmptbl",autol='no',mark="circle", number='yes',nyoffset=lra,radi=a2,txsize=2,inter='no')
   s1 = 1
   s2 = -int(fwhm0)
   iraf.delete("tmptbl",ve='no')
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" FITTED") 
   ff.close()  
   iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)
   iraf.display("skyfit",1,erase='no',fill='yes',xcen=.75,ycen=.25,xsize=.3,ysize=.3, zscal='no', zrang='no' ,z2=z22, z1=z11) 
   s1 = 1
   s2 = -int(fwhm0)
   iraf.delete("tmptbl",ve='no')
   ff=open('tmptbl','w')
   ff.write(str(s1)+' '+str(s2)+" RESIDUAL") 
   ff.close()
   iraf.tvmark(1,"tmptbl",autol='no',mark="none",inter='no',label='yes',txsize=2)

#   newmag=truemag
#   for i in range(0,len(truemag)):
#       try:
#           newmag[i]=float(truemag[i])+float(dmag0)
#       except:
#           newmag[i]=float(dmag0)
   newmag = list(array(truemag)+float(dmag0))   

   print "***************************************************************************" 
   print "#id  x_ori   y_ori     x     y    ap_ori ap_bgsub  fit_mag  err_art  err_fit" 
   for i in range(len(fitmag)):
       print "SN",i,str(centx[i]+x1),str(centy[i]+y1),str(centx[i]),str(centy[i]),"  ",str(apori3[i]),"  ",str(apmag3[i]),"  ",str(newmag[i]),"  ",str(arterr),"  ",str(magerr[i])
   print "**************************************************************************"

   return apori1,apori2,apori3,apmag1,apmag2,apmag3,fitmag,truemag,magerr,centx,centy,newmag


###################   DEFINE ERROR     ######################################

def errore(size,truemag,fwhm0,leng0):
    from pyraf import iraf
    import string,os,sys

    dartf=100
    while dartf>=size-1:
        artfac0 =raw_input('>>> Dispersion of artificial star positions (in units of FWHM) [1] ')
        if not artfac0: artfac0=1
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
    iraf.delete("tmpar?",ve='no')
    i=0
    tmpart=[]
    while i<=8:
        iraf.delete("reserr.fit?,artbg.fit?,artstar.fit?,artres.fit?,artfit.fit?,artskyfit.fit?",ve='no')
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
        print xbb,ybb

        iraf.delete("artlist.coo",ve='no')
        ff=open("artlist.coo",'w')
        ff.write(str(xbb)+'  '+str(ybb)+'  '+str(truemag[0])+"  1")
        ff.close()
                        
        xb1 = int(float(xbb)-fwhm0*float(leng0)/2)
        xb2 = int(float(xbb)+fwhm0*float(leng0)/2) 
        yb1 = int(float(ybb)-fwhm0*float(leng0)/2)
        yb2 = int(float(ybb)+fwhm0*float(leng0)/2)
        sec="1 "+str(xb1)+" 1 "+str(nay)+'\n'
        sec=sec+str(xb2)+' '+str(nax)+" 1 "+str(nay)+'\n'
        sec=sec+str(xb1)+' '+str(xb2)+" 1 "+str(yb1)+'\n'
        sec=sec+str(xb1)+' '+str(xb2)+' '+str(yb2)+' '+str(nay)+'\n'
        ff = open('sec','w')
        ff.write(sec)
        ff.close()
        
        iraf.delete("reserr.art",ve='no')
        iraf.delete("artlist.mag",ve='no')
        iraf.delete("artlist.als",ve='no')
        iraf.addstar("skyfit","artlist.coo",imgpsf,"reserr",nstar=1,veri='no',simple='yes')
        iraf.imsurfit("reserr","artbg",xorder=xbgord0,yorder=ybgord0,regions="section",section="sec")  
        iraf.imarith("reserr","-","artbg","artstar")
        iraf.phot("artstar","artlist.coo","artlist.mag",veri='no')
        iraf.allstar("artstar","artlist.mag",imgpsf,"artlist.als","","artres",veri='no',verb='no')
        iraf.imarith("artstar","-","artres","artfit")
        iraf.imarith("reserr","-","artfit","artskyfit")
        if i==0: era='yes'
        else:    era='no'
        artx = .5+.25*artx
        arty = .5+.25*arty
        iraf.display("artskyfit",1,fill='yes',erase=era,xcen=artx,ycen=arty,xsize=.25,ysize=.25)
	try:
	   tmpart.append(float(iraf.txdump("artlist.als","mag",expr='yes', Stdout=1)[0]))
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
    iraf.delete("reserr.fit?,artbg.fit?,artstar.fit?,artres.fit?,artfit.fit?,artskyfit.fit?",ve='no')
    iraf.delete("reserr.art",ve='no')
    iraf.delete("artlist.*",ve='no')
    return arterr2,arterr

#################################################################################################################
#################   FWHM COMPUTATION      #########
def fwhm_computation(img,tmpfiltercoo,stars,system,interactive):
    from snoopy2 import src
    import snoopy2
    import string,os
    from numpy import compress
    from numpy import array
    from numpy import average
    from numpy import mean
    from numpy import median
    from numpy import std
    from pyraf import iraf
    from numpy import log10
    _telescope=src.telescope(img)
    _header=src.read_parameter(_telescope,system)
    _xdimen=src.xdimen(img,_telescope)
    _ydimen=src.ydimen(img,_telescope)
    _filter=src.filter(img,_header,_telescope)
    try:
        filter=src.filtername(_telescope,_filter,system)
    except:
        filter=_filter
    _exptime=src.exptime(img,_header,_telescope)

    iraf.delete("tmp.imex_output",verify='no')
    stars_pos=[]
    alllines=[]

    ff = open(tmpfiltercoo,'r')
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
    ff = open('tmp.'+img+'good.coo','w')
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

    #fwhm[filter]=fwhm_ave
    iraf.hedit(img,'qubafwhm',fwhm_ave,add='yes',update='yes',verify='no')
    return fwhm_ave,_magime,stars_pos,_fwhm

###########################################################################################
def run_isis(isa,isb,nstamps_x,nstamps_y,sub_x,sub_y,saturation,degbg,degsp,_direc):
    import os,sys,string
    norm=''
    os.system('rm -rf default_config')
    default_config="nstamps_x     "+str(nstamps_x)+" /*** Number of stamps along X axis ***/     \n"
    default_config=default_config+"nstamps_y     "+str(nstamps_y)+" /*** Number of stamps along Y axis***/  \n"
    default_config=default_config+"sub_x "+ str(sub_x)+" /*** Number of sub_division of the image along X axis ***/ \n"
    default_config=default_config+"sub_y "+ str(sub_y)+" /*** Number of sub_division of the image along Y axis ***/ \n"
    default_config=default_config+"half_mesh_size  "+str(10)+" /*** Half kernel size ***/ \n"
    default_config=default_config+"half_stamp_size "+str(15)+" /*** Half stamp size ***/ \n"
    default_config=default_config+"deg_bg "+str(degbg)+" /** degree to fit differential bakground variations **/    \n"
    default_config=default_config+"saturation  "+str(saturation)+" /** degree to fit background varaitions **/ \n"
    default_config=default_config+"pix_min "+str(5.0)+" /*** Minimum vaue of the pixels to fit ****/             \n"
    default_config=default_config+"min_stamp_center "+str(130)+" /*** Minimum value for object to enter kernel fit *****/ \n"
    default_config=default_config+"ngauss "+str(3)+" /*** Number of Gaussians ****/ \n"
    default_config=default_config+"deg_gauss1 "+str(6)+"      /*** Degree associated with 1st Gaussian ****/    \n"
    default_config=default_config+"deg_gauss2 "+str(4)+"      /*** Degree associated with 2nd Gaussian ****/    \n"
    default_config=default_config+"deg_gauss3 "+str(2)+"      /*** Degree associated with 3rd Gaussian ****/    \n"
    default_config=default_config+"sigma_gauss1 "+str(0.7)+" /*** Sigma of 1st Gaussian ****/                   \n"
    default_config=default_config+"sigma_gauss2 "+str(1.5)+" /*** Sigma of 2nd Gaussian ****/                   \n"
    default_config=default_config+"sigma_gauss3 "+str(2.5)+" /*** Sigma of 3rd Gaussian ****/                   \n"
    default_config=default_config+"deg_spatial "+str(degsp)+" /*** Degree of the fit of the spatial variations of the Kernel ****/  \n"
    ff = open('default_config','w')
    ff.write(default_config)
    ff.close()
    ### RUN ISIS
    os.system(_direc+"/mrj_phot "+isa+"  "+isb+" -c ./default_config > aaaa ")# | csh
    f=open('aaaa','r')
    a=f.readlines()
    f.close()
    for i in a:
        if string.count(i,'sum_kernel')>=1:
            norm=float(string.split(i)[-1])
    if not norm:
            norm='INDEF'
    return norm

##########################################################################################

def isis_parameter(_nstamps_x,_nstamps_y,_sub_x,_sub_y,_saturation,_degbg,_degsp):
        answ='n'
        while answ=='n':
            nstamps_x=raw_input(">>> nstamp x ["+str(_nstamps_x)+"] ? ")
            try:
                if not nstamps_x: nstamps_x=float(_nstamps_x)
                else: nstamps_x=float(nstamps_x)
                answ='y'
            except:
                print 'Warning value for nstamps_x not valid'
                answ='n'
        answ='n'
        while answ=='n':
            nstamps_y=raw_input(">>> nstamp y ["+str(_nstamps_y)+"] ? ")
            try:
                if not nstamps_y: nstamps_y=float(_nstamps_y)
                else: nstamps_y=float(nstamps_y)
                answ='y'
            except:
                print 'Warning value for nstamps_y not valid'
                answ='n'    
        answ='n'
        while answ=='n':
            sub_x=raw_input(">>> sub x ["+str(_sub_x)+"] ? ")
            try:
                if not sub_x: sub_x=float(_sub_x)
                else: sub_x=float(sub_x)
                answ='y'
            except:
                print 'Warning value for sub_x not valid'
                answ='n'
        answ='n'
        while answ=='n':
            sub_y=raw_input(">>> sub y ["+str(_sub_y)+"] ? ")
            try:
                if not sub_y: sub_y=float(_sub_y)
                else: sub_y=float(sub_y)
                answ='y'
            except:
                print 'Warning value for sub_y not valid'
                answ='n'
        answ='n'
        while answ=='n':
            degbg=raw_input(">>> degbg  ["+str(_degbg)+"] ? ")
            try:
                if not degbg: degbg=float(_degbg)
                else: degbg=float(degbg)
                answ='y'
            except:
                print 'Warning value for degbg not valid'
                answ='n'
        answ='n'
        while answ=='n':
            degsp=raw_input(">>> degsp ["+str(_degsp)+"] ? ")
            try:
                if not degsp: degsp=float(_degsp)
                else: degsp=float(degsp)
                answ='y'
            except:
                print 'Warning value for degsp not valid'
                answ='n'
        return nstamps_x,nstamps_y,sub_x,sub_y,_saturation,degbg,degsp

#########################################################################

def daophot_parameters(_datamin,_datamax,aperture2,annulus2,_fwhm_target):
    print '### try to change the fit parameters'
    print 'fwhm target= '+str(_fwhm_target) 
    print 'default annulus = 4 '
    print 'default aperture= 3 '
    print 'try to take larger annulus and aperture or change datamin and datamax'
    _datamin2=raw_input('datamin ? ['+str(_datamin)+'] ')
    if not _datamin2:_datamin2=_datamin
    _datamax2=raw_input('datamax ? ['+str(_datamax)+'] ')
    if not _datamax2:_datamax2=_datamax
    _aperture2=raw_input('aperture ? ['+str(aperture2)+'] ')
    if not _aperture2:_aperture2=aperture2
    _annulus2=raw_input('annulus ? ['+str(annulus2)+'] ')
    if not _annulus2:_annulus2=annulus2
    return _annulus2,_aperture2,_datamin2,_datamax2
