

def ReadAscii3(ascifile):
    f=file(ascifile,'r')
    s=f.readlines()
    f.close()
    vecb1,vecb2,vecb3 =  [],[],[]
    for x in s:
        if x[0]!="#":
         c1,c2,c3 = x.split()
         vecb1.append(c1)
         vecb2.append(c2)
         vecb3.append(c3)
    return vecb1,vecb2,vecb3

def ReadAscii33(ascifile):
    import re,string
    f=file(ascifile,'r')
    s=f.readlines()
    f.close()
    vecb1,vecb3 =  [],[]
    for x in s:
        if x[0]!="#":
         c1,c3 = x.split('=')
         vecb1.append(re.sub(' ','',re.sub('\n','',re.sub('\t','',c1))))
         vecb3.append(re.sub('\n','',re.sub('\t','',c3)))
    return vecb1,'ss',vecb3


#######################################################################

def read_parameter(telescopio,system):
    import snoopy2
    subdirectory=['optical/','infrared/','sloan/']
    dir_system=subdirectory[system]
    directory =snoopy2.__path__[0]+'/telescope/'
    read_file = directory +dir_system+'keywords_'+str(telescopio)+'.txt'
    try:
        keyword,pp,name_keyword = ReadAscii33(read_file)
        header={}
        for i in range(len(keyword)): header[keyword[i]]=name_keyword[i]
    except:
        header=''
    return header

############################################################################

def  close_program(logincl):
   import os,sys    
   os.system('rm -rf login.cl')
   if logincl:
      f=open('login.cl','w')
      for i in logincl:
         f.write(i)
      f.close()
   if os.path.isdir('uparm'):
       os.system('rm -rf uparm/')
   sys.exit()    


def  open_program(directory='',interactive=True):
   import snoopy2
   if not directory:
       directory=snoopy2.__path__[0]
   import os,sys,shutil,string
   parent_dir = os.getcwd()+'/'
   logincl=''
   if interactive:
       if os.path.isfile('login.cl') and not os.path.isdir('uparm'):
           answ=raw_input('WARNING: There is a login.cl file in this directory. Can I cancel it [y/n] [y]?')
           if not answ: answ='y'
       elif os.path.isfile('login.cl') and os.path.isdir('uparm'):
           answ=raw_input('WARNING: There is a login.cl file and uparm dir in this directory. Can I cancel them [y/n] [y]?')
           if not answ: answ='y'
       elif not os.path.isfile('login.cl') and os.path.isdir('uparm'):
           answ=raw_input('WARNING: There is an uparm dir in this directory. Can I cancel it [y/n] [y]?')
           if not answ: answ='y'
       else:
           answ='y'
   else:
       answ='y'
   
   if answ in ['y','Y','YES','Yes']:
       if os.path.isfile('login.cl'):
           f=open('login.cl','r')
           logincl=f.readlines()
           f.close()    
           os.system('rm -rf login.cl')
       else:
           logincl=''

       if os.path.isdir('uparm'):
           os.system('rm -rf uparm/')

       filelogin=directory+'/login.cl'
       f=open(filelogin,'r')
       ss=f.readlines()
       f.close()

       ff=open('login.cl','w')
       for i in ss:
           if string.count(i,'uparm')==2:
               ff.write('set uparm = "'+parent_dir+'uparm/"\n')
           elif string.count(i,'set home = ')==1:
               ff.write('set home = "'+snoopy2.__path__[0]+'/"\n')
           else:
               ff.write(i)
       ff.close()
       os.system('mkdir uparm/')
       os.system('cp '+directory+'/uparm/* '+parent_dir+'uparm/')           
   else:
       print 'ERROR: the procedure can not run in the main iraf directory'
       close_program(logincl)
   return logincl

#########################################################################      

def open_ds9():
    import snoopy2
    from snoopy2 import src
    import os
    try:
        enviro=os.environ['IRAFARCH']
        if enviro=='macintel':
            directory_prog=snoopy2.__path__[0]+'/exec_prog/macintel'
            os.system(directory_prog+'/SAOImage\ DS9.app/Contents/MacOS/ds9&')
        elif enviro=='redhat' or enviro=='linux':
            directory_prog=snoopy2.__path__[0]+'/exec_prog/redhat'
            os.system(directory_prog+'/ds9&')
        else:
            print 'ERROR: architecture not found !!'
            os.system('ds9&')
    except:
        try:
            directory_prog=snoopy2.__path__[0]+'/exec_prog/redhat'
            os.system(directory_prog+'/ds9&')
        except:
            print 'ERROR: architecture not found !!'
            os.system('ds9&')      

def display_image(img,frame,_z1,_z2,scale,_xcen=0.5,_ycen=0.5,_xsize=1,_ysize=1,_erase='yes'):
    goon='True'
    import glob
    if glob.glob(img):
       import snoopy2
       from snoopy2 import src
       from pyraf import iraf
       import string,os
       import time
       #img=src.replace_directory(img)
       if _z2: 
          try:
              sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2,Stdout=1)
              #sss=iraf.display(img, frame, fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2,Stdout=1)
          except:
              try:
                  src.open_ds9()
                  time.sleep(5)
                  sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2,Stdout=1)
              except:
                  print ''
                  print '### ERROR: PROBLEM OPENING DS9'
                  print ''
                  goon='False'                 
       else:
        try:  
            sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', Stdout=1)
        except:
            try:
                src.open_ds9()
                time.sleep(5)
                sss=iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, fill='yes', Stdout=1)
            except:
                print ''
                print '### ERROR: PROBLEM OPENING DS9'
                print ''
                goon='False'                 
        if scale and goon:
          answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
          if not answ0: answ0='y'
          elif answ0=='no' or answ0=='NO': answ0='n' 

          while answ0=='n':
              _z11=float(string.split(string.split(sss[0])[0],'=')[1])
              _z22=float(string.split(string.split(sss[0])[1],'=')[1])
              z11 = raw_input('>>> z1 = ? ['+str(_z11)+'] ? ')
              z22 = raw_input('>>> z2 = ? ['+str(_z22)+'] ? ')
              if not z11: z11=_z11
              else: z11=float(z11)
              if not z22: z22=_z22
              else: z22=float(z22)
              print z11,z22
              sss=iraf.display(img,frame,fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase, zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
              answ0 = raw_input('>>> Cuts OK ? [y/n] ? [y] ')
              if not answ0: answ0='y'
              elif answ0=='no' or answ0=='NO': answ0='n'
       if goon:
          _z1,_z2=string.split(string.split(sss[0])[0],'=')[1],string.split(string.split(sss[0])[1],'=')[1]
    else:
        print 'Warning: image '+str(img)+' not found in the directory '
    return _z1,_z2,goon

##########################################################################
def readlist(listfile):
    import string,os,sys,re
    from snoopy2 import src
    from astropy.io import fits as pyfits
    
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:
        try:
            hdulist= pyfits.open(listfile)
        except:
            hdulist=[]
            if listfile[-5:]!='.fits': 
                listfile=listfile+'.fits'
                try:
                    hdulist = pyfits.open(listfile)
                except:
                    hdulist=[]
        if not hdulist:
#        if 'Error reading' in iraf.imstat(listfile,Stdout=1)[1]:
            print '\n####### Error ######'
            print 'If "'+str(listfile)+'" is an image, it is corrupted ...'
            print 'or you just forgot the "@" before the list name ...\n'
            sys.exit()
        else:
            imglist = [listfile]
    return imglist


def delete(listfile):
    import os,string,re,glob
    if listfile[0]=='@':   
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files: 
            ff=re.sub(' ','',ff)
            if not ff=='\n' and ff[0]!='#':
                ff=re.sub('\n','',ff)
                imglist.append(ff)
    elif ',' in listfile: imglist = string.split(listfile,sep=',')
    else:
        imglist=[listfile]
    
    lista=[]
    for _file in imglist:
        lista=lista+glob.glob(_file)
    if lista:
        for _file in lista:
            try:
                os.system('rm '+_file)
            except:
                pass
 ###############################################################

def pyimarith(img1,operation,img2,imgout):
    from astropy.io import fits as pyfits

#    from pyfits import writeto
    from snoopy2 import src
    import os,sys
    hdus1 = pyfits.open(img1)
    hdus1[0].verify('fix')
    if os.path.isfile(str(img2)):
        hdus2 = pyfits.open(img2)
        hdus2[0].verify('fix')
        if operation=='+':
            hdus1[0].data=hdus1[0].data+hdus2[0].data
        elif operation=='-':
            hdus1[0].data=hdus1[0].data-hdus2[0].data
        elif operation=='*':
            hdus1[0].data=hdus1[0].data*hdus2[0].data
        elif operation=='/':
            hdus1[0].data=hdus1[0].data/hdus2[0].data
        else:
            print 'warning: operation not recognised '
    elif float(img2):
        if operation=='+':
            hdus1[0].data=hdus1[0].data+float(img2)
        elif operation=='-':
            hdus1[0].data=hdus1[0].data-float(img2)
        elif operation=='*':
            hdus1[0].data=hdus1[0].data*float(img2)
        elif operation=='/':
            hdus1[0].data=hdus1[0].data/float(img2)
        else:
            print 'warning: operation not recognised '
    else:
        print 'warning: second imput not an image and not a float '
    if os.path.isfile(imgout):
        os.system('rm '+imgout)
    hdus1.writeto(imgout) # write current form to new file

def pycopy(img1,imgout):
    from astropy.io import fits as pyfits
    from pyfits import writeto
    
    import os,sys,string
    from numpy  import asarray
    from numpy import zeros
    # read input
    if string.count(img1,'[')== 1:   #and string.count(imgout,'[')== 0 :
        wind=string.split(string.split(img1,'[')[1],']')[0]
        name=string.split(img1,'[')[0]
        hdus1 = pyfits.open(name)
        a,c=string.split(wind,',')
        aa,bb=string.split(a,':')
        cc,dd=string.split(c,':')
        aa=int(aa)-1
        bb=int(bb)
        cc=int(cc)-1
        dd=int(dd)
        hdus1[0].data=hdus1[0].data[cc:dd,aa:bb]
    else:
        name=img1
        hdus1 = pyfits.open(name)
    if string.count(imgout,'[')== 1:
        imgoutcut=string.split(imgout,'[')[0]
    else:
        imgoutcut=imgout
    # read output
    if not os.path.isfile(imgoutcut):
        hdus1.writeto(imgoutcut) # write current form to new file
    else:
        if string.count(imgout,'[')== 1:
            wind2=string.split(string.split(imgout,'[')[1],']')[0]
            name2=string.split(imgout,'[')[0]
            hdus2 = pyfits.open(name2)
            a2,c2=string.split(wind2,',')
            aa2,bb2=string.split(a2,':')
            cc2,dd2=string.split(c2,':')
            aa2=int(aa2)-1
            bb2=int(bb2)
            cc2=int(cc2)-1
            dd2=int(dd2)
            cut='yes'
        else:
            name2=imgout
            hdus2 = pyfits.open(name2)
            cut='no'
        os.system('rm '+imgoutcut)
        if cut=='no':
            hdus1.writeto(imgout) # write current form to new file
        else:
            hdus2[0].data[cc2:dd2,aa2:bb2]=hdus1[0].data
            hdus2.writeto(imgoutcut) # write current form to new file


def correctcard(img):
    from astropy.io import fits as pyfits

    from numpy  import asarray
    import re
    hdulist=pyfits.open(img)
    a=hdulist[0]._verify('fix')    
    _header=hdulist[0].header
    for i in range(len(a)):
        if not a[i]:
            a[i]=['']
    ww=asarray([i for i in range(len(a)) if (re.sub(' ','',a[i][0])!='')])
    if len(ww)>0:
        newheader=[]
        headername=[]
        for j in _header.items():
            headername.append(j[0])
            newheader.append(j[1])
        hdulist.close()
        imm=pyfits.open(img,mode='update')
        _header=imm[0].header
        for i in ww:
            if headername[i]:
                try:
                    _header.update(headername[i],newheader[i])
                except:
                    _header.update(headername[i],'xxxx')
        imm.flush()
        imm.close()


def deg2HMS(ra='', dec='', round=False):
      import string
      RA, DEC= '', ''
      if dec:
          if string.count(str(dec),':')==2:
              dec00=string.split(dec,':')
              dec0,dec1,dec2=float(dec00[0]),float(dec00[1]),float(dec00[2])
              if '-' in str(dec0):       DEC=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
              else:                      DEC=(dec2/60.+dec1)/60.+dec0
          else:
              if str(dec)[0]=='-':      dec0=(-1)*abs(int(dec))
              else:                     dec0=abs(int(dec))
              dec1=int((abs(dec)-abs(dec0))*(60))
              dec2=((((abs(dec))-abs(dec0))*60)-abs(dec1))*60
              DEC='00'[len(str(dec0)):]+str(dec0)+':'+'00'[len(str(dec1)):]+str(dec1)+':'+'00'[len(str(int(dec2))):]+str(dec2)
      if ra:
          if string.count(str(ra),':')==2:
              ra00=string.split(ra,':')
              ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
              RA=((ra2/60.+ra1)/60.+ra0)*15.
          else:
              ra0=int(ra/15.)
              ra1=int(((ra/15.)-ra0)*(60))
              ra2=((((ra/15.)-ra0)*60)-ra1)*60
              RA='00'[len(str(ra0)):]+str(ra0)+':'+'00'[len(str(ra1)):]+str(ra1)+':'+'00'[len(str(int(ra2))):]+str(ra2)
      if ra and dec:          return RA, DEC
      else:                   return RA or DEC
