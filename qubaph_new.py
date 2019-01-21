#!/usr/bin/env python

from numpy import array
from numpy import sort
from numpy import arange
from numpy import argsort

import os,sys,string,re,shutil
import getopt
import math
from snoopy2 import *
import snoopy2
import sys

help ="################################################################ \n\
help = Usage:   ph_new.py filename             \n\
                Input format: the output file of svstandards.py (E.g. lp_ru149_20080212)    \n\
                [-r  list] list format: the input file of svstandard.py (E.g. lista_lp_20080212)   \n\
                -i  interactive       ** FOR UBUNTU USERS **  \n\
                -f  fix colour term  \n\
                -s  system num    (landolt [0],infrared [1], sloan [2])   \n\
                -w  write constant                                        \n\
################################################################"

if len(sys.argv)==1:
    print help
    sys.exit()

def updatecolorlog(line):
    from snoopy2 import src
    import snoopy2
    import string,os
    from numpy import compress
    listacolor=snoopy2.__path__[0]+'/standard/photlog.txt'
    f=open(listacolor,'r')
    liststd=f.readlines()
    f.close()
    _telescope0=string.split(line)[0]
    _date0=string.split(line)[1]
    _ut0=string.split(line)[2]
    _airmass0=string.split(line)[3]
    _band0=string.split(line)[4]
    _color0=string.split(line)[6]
    _telescope,_date,_ut,_airmass,_const=[],[],[],[],[]
    _band,_color=[],[]
    _recon=[]
    for i in range(0,len(liststd)):
        if liststd[i][0]!='#':
            _recon.append([string.split(liststd[i])[0],string.split(liststd[i])[2],string.split(liststd[i])[4],string.split(liststd[i])[6]])
            _telescope.append(string.split(liststd[i])[0])
            _date.append(string.split(liststd[i])[1])
            _ut.append(string.split(liststd[i])[2])
            _airmass.append(string.split(liststd[i])[3])
            _band.append(string.split(liststd[i])[4])
            _color.append(string.split(liststd[i])[6])
            if _recon[-1]==[_telescope0,_ut0,_band0,_color0]:
#                print _recon[-1],[_telescope0,_ut0,_band0,_color0]
                liststd[i]=line
                break

    _date1=compress((array(_color)==_color0)&(array(_telescope)==_telescope0)&(array(_band)==_band0)&(array(_date)==_date0),array(_date))
    if len(_date1)==0:
        liststd.append(line)
    else:
        print 'already there'
    aa=os.umask(0)
    try:
        f=open(listacolor,'w')
        for i in liststd:
            f.write(i)
        f.close()
    except:
        print 'Warning: it seems you do not have permission to write in '+snoopy2.__path__[0]+'/standard/photlog.txt '
    bb=os.umask(aa)
    return line



def readqubacons(img):
    import string
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _qubacons= geth['QUBACONS']
        print _qubacons
        try:
            _zeropoint=float(string.split(_qubacons)[1])
            _colore=string.split(_qubacons)[2]
            _colorevalue=float(string.split(_qubacons)[3])
        except:
            _zeropoint=''
            _colore=''
            _colorevalue=''
        try:
            _zeroerr=float(string.split(_qubacons)[5])
            _colorerr=float(string.split(_qubacons)[6])
        except:
            _zeroerr=0
            _colorerr=0
    except:
            _zeropoint=''
            _colore=''
            _colorevalue=''
            _zeroerr=0
            _colorerr=0
    return _zeropoint,_colore,_colorevalue,_zeroerr,_colorerr

#######################################################################
_system=''
p0=[0,0,0,0]
parent_dir = os.getcwd()+'/'

fixcolor=False
lista=''
fisso=0
interactive=False
writecon=False

options, args = getopt.getopt(sys.argv[2:],'r:,i,f,s:,w',['list=','interative','fixcolor','system=','write'])
for opt,arg in options:
    if opt in ('-r', '--list'): lista= arg
    if opt in ('-i', '--interactive'): interactive=True
    if opt in ('-f', '--fixcolor'): fixcolor = True
    if opt in ('-s', '--system'): _system = int(arg)
    if opt in ('-w', '--write'): writecon=True

print lista,'####'
if lista:
    imglist=src.readlist(lista)

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
          print '############################################'
          print "### Error with the header !!!!             #"
          print '### telescope not in the list              #'
          print '### if you want to continue anyway,        #'
          print '### run "svother.py" to correct headers    #'
          print '########################################### '
          src.close_program(logincl)

    _header=src.read_parameter(_telescope,system)
    _zero,_color={},{}
    _zeroerr,_colorerr={},{}
    for img in imglist:
        print img
        _filter=src.filter(img,_header,_telescope)
        filter=src.filtername(_telescope,_filter,system)
        if filter!='unknown':
            _zeropoint,_colore,_colorevalue,zeroerr,colorerr=readqubacons(img)
            _zero[filter]=_zeropoint
            _color[filter]=_colorevalue
            _zeroerr[filter]=zeroerr
            _colorerr[filter]=colorerr
        else:
                 print 'WARNING: The image form is not recognised.'

##########################################################################
else:
    imglist=''
    system='none'
    while system=='none':
       system=raw_input('Which photometric system? UBVRI [0], JHK [1], ugriz [2] ')
       if system=='0' or system=='1' or system=='2':
           system=int(system)
       else:
	   system='none' 
#k=src.atmospheric_site(system)
k=src.atmospheric_site2()
qubacont={}
again='yes'

while again=='yes':
   if system==0:            
       try:
           mu,mb,mv,mr,mi,st,V,BV,UB,VR,RI,VI,field,name_field=src.Readphfile_o(sys.argv[1],k)
       except:
           print '#### WARNING: File format not correct!'
   
       _y={}
       _y['U']=mu
       _y['B']=mb
       _y['V']=mv
       _y['R']=mr
       _y['I']=mi
 
       band=raw_input('Define band [B]')
       if not band: band='B'
       y=_y[band]
       color=raw_input('Define colour [BV]')
       if not color: color='BV'

       U,B,R,I={},{},{},{}
       x={}
       for i in st:
            B[i]=list(array(V[i])+array(BV[i]))
            R[i]=list(array(V[i])-array(VR[i]))
            U[i]=list(array(B[i])+array(UB[i]))
            I[i]=list(array(R[i])-array(RI[i]))
            if color=='BV':    x[i]=list(array(B[i])-array(V[i]))
            elif color=='VR':  x[i]=list(array(V[i])-array(R[i]))
            elif color=='RI':  x[i]=list(array(R[i])-array(I[i]))
            elif color=='UB':  x[i]=list(array(U[i])-array(B[i]))
            elif color=='VI':  x[i]=list(array(V[i])-array(I[i]))
            elif color=='BR':  x[i]=list(array(B[i])-array(R[i]))
            elif color=='UR':  x[i]=list(array(U[i])-array(R[i]))
            elif color=='BI':  x[i]=list(array(B[i])-array(I[i]))
            elif color=='UV':  x[i]=list(array(U[i])-array(V[i]))

       _yy=[]
       _xx=[]
       _position=[]
       yy=[]
       xx=[]
       position=[]
       _label=[]
       label=[]
       for i in range(len(x)):
         for j in range(len(x[i+1])):
            if abs(y[i+1][j])<=40 and abs(x[i+1][j])<=8:
              _yy.append(y[i+1][j])
              _xx.append(x[i+1][j])
              _position.append([i,j])
  	      _label.append(st[i+1][j]) 
       
       xx=sort(_xx)
       xx_elem=argsort(_xx)
       for i in range(len(xx)):
              yy.append(_yy[xx_elem[i]])
              position.append(_position[xx_elem[i]])
              label.append(_label[xx_elem[i]])
       try:
             phmin=min(xx)-(max(xx)-min(xx))/10.
             phmax=max(xx)+(max(xx)-min(xx))/10.
       except:
             phmin=-1
             phmax=2

       phmin=min(xx)-(max(xx)-min(xx))/10.
       phmax=max(xx)+(max(xx)-min(xx))/10.
       if (max(yy)-min(yy))>=1:
           mm=min(yy)-(max(yy)-min(yy))/1.
           mma=max(yy)+(max(yy)-min(yy))/1.
       else:
           mm=min(yy)-.5
           mma=max(yy)+.5

       if fixcolor:
               try:
                   print '######' 
                   print band,_zero[band],_color[band]
                   _fisso= _color[band]
               except: 
                   print '### WARNING: No color term in the header '
                   print '#### ',band,'   ',color
                   _fisso='xxx'
               fisso=raw_input('### What is the value of the color term you want to use ['+str(_fisso)+'] ? ')
               if not fisso: fisso=_fisso
               fisso=float(fisso)
               sss=band
               f=color
               a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,fisso)
               value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
               qubacont[band]=value
       else:
           sss=band
           f=color
           a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,'')
           value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
           qubacont[band]=value
   elif system==1:
      try:
         mj,mh,mk,st,J,JH,HK,JK,field,name_field=src.Readphfile_i(sys.argv[1],k)
      except:
          print '#### WARNING: file format not correct !'

      _y={}
      _y['J']=mj
      _y['H']=mh
      _y['K']=mk
      band=raw_input('Define band [J]')
      if not band: band='J'
      y=_y[band]
      color=raw_input('Define colour [JH]')
      if not color: color='JH'

      H,K={},{}
      x={}
      for i in st:
           H[i]=list(array(J[i])-array(JH[i]))
           K[i]=list(array(H[i])-array(HK[i]))
           if color=='JH':    x[i]=list(array(J[i])-array(H[i]))
           elif color=='HK':  x[i]=list(array(H[i])-array(K[i]))
           elif color=='JK':  x[i]=list(array(J[i])-array(K[i]))

      _yy=[]
      _xx=[]
      _position=[]
      yy=[]
      xx=[]
      position=[]
      _label=[]
      label=[]
      for i in range(len(x)):
         for j in range(len(x[i+1])):
            #if abs(y[i+1][j])<=40:
            if abs(y[i+1][j])<=40 and abs(x[i+1][j])<=8:
               _yy.append(y[i+1][j])
               _xx.append(x[i+1][j])
               print _xx[-1]
               _position.append([i,j])
	       _label.append(st[i+1][j])

      xx=sort(_xx)
      xx_elem=argsort(_xx)
      for i in range(len(xx)):
         yy.append(_yy[xx_elem[i]])
         position.append(_position[xx_elem[i]])
         label.append(_label[xx_elem[i]])
         try:
             phmin=min(xx)-(max(xx)-min(xx))/10.
             phmax=max(xx)+(max(xx)-min(xx))/10.
         except:
             phmin=-1
             phmax=2
      if (max(yy)-min(yy))>=1:
          mm=min(yy)-(max(yy)-min(yy))/1.
          mma=max(yy)+(max(yy)-min(yy))/1.
      else:
          mm=min(yy)-.5
          mma=max(yy)+.5

      if fixcolor:
               try:
                   print '######' 
                   print band,_zero[band],_color[band]
                   _fisso= _color[band]
               except: 
                   print '### WARNING: No color term in the header '
                   print '#### ',band,'   ',color
                   _fisso='xxx'
               fisso=raw_input('### What is the value of the color term you want to use ['+str(_fisso)+'] ? ')
               if not fisso: fisso=_fisso
               fisso=float(fisso)
               sss=band
               f=color
               a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,fisso)
               value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
               qubacont[band]=value
      else:
          sss=band
          f=color
          a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,'')
          value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
          qubacont[band]=value

   elif system==2:
     try:
        mu,mg,mr,mi,mz,st,r,gr,ug,ri,iz,gi,field,name_field=src.Readphfile_s(sys.argv[1],k)
     except:
        print '#### WARNING: File format incorrect !'

     _y={}
     _y['u']=mu
     _y['g']=mg
     _y['r']=mr
     _y['i']=mi
     _y['z']=mz

     band=raw_input('define band [g]')
     if not band: band='g'
     y=_y[band]
     color=raw_input('define color [gr]')
     if not color: color='gr'

     u,g,i,z={},{},{},{}
     x={}
     for ii in st:
           g[ii]=list(array(r[ii])+array(gr[ii]))
           u[ii]=list(array(g[ii])+array(ug[ii]))
           i[ii]=list(array(r[ii])-array(ri[ii]))
           z[ii]=list(array(i[ii])-array(iz[ii]))
           if color=='ug':    x[ii]=list(array(u[ii])-array(g[ii]))
           elif color=='gr':  x[ii]=list(array(g[ii])-array(r[ii]))
           elif color=='ri':  x[ii]=list(array(r[ii])-array(i[ii]))
           elif color=='iz':  x[ii]=list(array(i[ii])-array(z[ii]))
           elif color=='gi':  x[ii]=list(array(g[ii])-array(i[ii]))
           elif color=='gz':  x[ii]=list(array(g[ii])-array(z[ii]))
           elif color=='rz':  x[ii]=list(array(r[ii])-array(z[ii]))
           elif color=='ur':  x[ii]=list(array(u[ii])-array(r[ii]))

     _yy=[]
     _xx=[]
     _position=[]
     yy=[]
     xx=[]
     position=[]
     _label=[]
     label=[]
     for i in range(len(x)):
         for j in range(len(x[i+1])):
             if abs(y[i+1][j])<=40 and abs(x[i+1][j])<=8:
             #if abs(y[i+1][j])<=40:
               _yy.append(y[i+1][j])
               _xx.append(x[i+1][j])
               _position.append([i,j])
	       _label.append(st[i+1][j])
       
     xx=sort(_xx)
     xx_elem=argsort(_xx)
     for i in range(len(xx)):
        yy.append(_yy[xx_elem[i]])
        position.append(_position[xx_elem[i]])
        label.append(_label[xx_elem[i]])

     try:
             phmin=min(xx)-(max(xx)-min(xx))/10.
             phmax=max(xx)+(max(xx)-min(xx))/10.
     except:
             phmin=-1
             phmax=2
             
     phmin=min(xx)-(max(xx)-min(xx))/10.
     phmax=max(xx)+(max(xx)-min(xx))/10.
     if (max(yy)-min(yy))>=1:
         mm=min(yy)-(max(yy)-min(yy))/1.
         mma=max(yy)+(max(yy)-min(yy))/1.
     else:
         mm=min(yy)-.5
         mma=max(yy)+.5
         
     if fixcolor:
               try:
                   print '######' 
                   print band,_zero[band],_color[band]
                   _fisso= _color[band]
               except: 
                   print '### WARNING: No color term in the header '
                   print '#### ',band,'   ',color
                   _fisso='xxx'
               fisso=raw_input('### What is the value of the color term you want to use ['+str(_fisso)+'] ? ')
               if not fisso: fisso=_fisso
               fisso=float(fisso)
               sss=band
               f=color
               a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,fisso)
               value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
               qubacont[band]=value
     else:
         a,sa,b,sb=qubaphdef.fitcol(xx,yy,band,color,'')
         value='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(band,a,color,b,sa,sb)
         qubacont[band]=value
#################################################
   again=raw_input('again [yes/no]? [no] ')
   if not again: agina='no'
   if again in ['Yes','YES','yes','Y','y']:
       again='yes'

if imglist:
   _header=src.read_parameter(_telescope,system)
   for band in qubacont:
      for img in imglist:
         try:
            _filter=src.filter(img,_header,_telescope)
            filtern=src.filtername(_telescope,_filter,system)
            if filtern=='unknown':
                filtern=_filter
            if filtern==band:
               _airmass=src.airmass(img,_header,_telescope)
               _date=src.date(img,_header,telescope)
               _JD=src.JD(img,_header,telescope) 
               src.updateheader(img,0,'qubacons',qubacont[band])
               print '###################################'
               print '#'
               print '# ',_telescope,_date,_JD,_airmass,qubacont[band]
               line='%8s\t%8s\t%8.8s\t%5.5s\t%30s\n' % (str(_telescope),str(_date),str(_JD+0.005),str(_airmass),str(qubacont[band]))
               if writecon:
                   xxx=updatecolorlog(line)
                   print 'update zeropoints with  '+line
               print '#'
               print '##################################'
         except:
            print 'WARNING: The image form is not recognised.'
else:
      print '###################################'
      print '#'
      for band in qubacont:
             print '# ',qubacont[band]
      print '#'
      print '##################################'
      if writecon:
          _telescope=raw_input('which telescope [NOT,TNG,NTT,ekar,TNG,FTS,FTN] [FTS]? ')
          if not _telescope: _telescope='FTS'

          _date=raw_input('which date  [20091123] ? ')
          if not _date: _date='20091123'

          _airmass=raw_input('which airmass  [0] ? ')
          if not _airmass: _airmass=0

          _hour=12.
          _min=0.0
          _year=int(_date[0:4])
          _month=int(_date[4:6])
          _day=int(_date[6:8])
          _JD=src.julday(_year,_month,_day,_hour,_min)
          for band in qubacont:
              line='%8s\t%8s\t%8.8s\t%5.5s\t%30s\n' % (str(_telescope),str(_date),str(_JD+0.005-2400000),str(_airmass),str(qubacont[band]))
              xxx=updatecolorlog(line)
              print 'update zeropoints with  '+line
