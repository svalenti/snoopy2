#!/usr/bin/env python

import getopt
import os,sys,shutil,string,re
from snoopy2 import *
import snoopy2
import glob
from pyfits import getheader

help ="################################################################   \n\
help = Usage:   svother.py filelist                                       \n\
                input                filelist iraf format                 \n\
      Program identifies header keywords for relevant keywords, useful for\n\
      preparing files to be used with qubaprephot                         \n\
      ################################################################"
if len(sys.argv)<=1:
    print help
    sys.exit()


imglist=src.readlist(sys.argv[1])
img=imglist[0]

userPath=os.getcwd()
options, args = getopt.getopt(sys.argv[2:],'p:',['userPath='])
for opt,arg in options:
    if opt in ('-p', '--path'): userPath = arg

######################################################################

src.correctcard(img)
aa=getheader(img)
print aa

check=0
while check==0:
    hedairmass=''
    print ' type [+] if keyword is missing '
    hedairmass = raw_input('Enter the header keyword for airmass [AIRMASS] ')
    if not hedairmass: hedairmass='AIRMASS'
    elif hedairmass =='+':
        src.updateheader(img,0,'AIRMASS',1)
        hedairmass='AIRMASS'
    try:
        _airmass=float(getheader(img)[hedairmass])
        print 'AIRMASS= '+str(_airmass)
        check=1
    except:
        print 'header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()

print aa        
check=0
while check==0:
    hedexptime=''
    print ' type [+] if it is missing '
    hedexptime = raw_input('Enter the header keyword for the exposure time ? [EXPTIME] ')
    if not hedexptime: 
        hedexptime='EXPTIME'
    elif hedexptime =='+':
        src.updateheader(img,0,'EXPTIME',1)
        hedexptime='EXPTIME'
    try:
        _exptime=float(getheader(img)[hedexptime])
        check=1
        print 'EXPTIME= '+str(_exptime)
    except:
        print 'Header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()

print aa
check=0
while check==0:
    hedobject=''
    print ' type [+] if it is missing '
    hedobject = raw_input('Enter the header keyword for object name? [OBJECT] ')
    if not hedobject: hedobject='OBJECT'
    elif hedobject =='+':
        _object = raw_input('Enter the object name? [SN] ')
        if not _object: _object='SN'
        src.updateheader(img,0,'OBJECT',_object)
        hedobject='OBJECT'
    try:
        _object=str(getheader(img)[hedobject])
        check=1
        print 'OBJECT= '+str(_object)
    except:
        print 'header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()


print aa
check=0
while check==0:
    heddate=''
    print ' type [+] if it is missing '
    heddate = raw_input('Enter the header keyword for the date [DATE-OBS] ')
    if not heddate: heddate='DATE-OBS'
    elif heddate =='+':
        date = raw_input('What DATE? [2009-01-01] ')
        if not date: date='2009-01-01'
        src.updateheader(img,0,'DATE-OBS',date)
        heddate='DATE-OBS'
    try:
        _date=str(getheader(img)[heddate])
        check=1
        print 'DATE= '+str(_date)
    except:
        print 'Header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else:
            sys.exit()


print aa
check=0
while check==0:
    hedinst=''
    print ' type [+] if keyword is missing '
    hedinst = raw_input('Enter the header keyword for the instrument [INSTRUME] ')
    if not hedinst: hedinst='INSTRUME'
    elif hedinst =='+':
        inst = raw_input('Enter the header keyword for the instrument ? [CCD] ')
        if not inst: inst='CCD'
        src.updateheader(img,0,'INSTRUME',inst)
        hedinst='INSTRUME'
    try:
        _instrument=str(getheader(img)[hedinst])
        check=1
        print 'INSTRUMENT= '+str(_instrument)
    except:
        print 'Header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()


print aa
check=0
while check==0:
    hedUT=''
    print ' type [+] if keyword is missing '
    hedUT = raw_input('Enter the header keyword for UT [UT] ')
    if not hedUT: hedUT='UT'
    elif hedUT =='+':
        UT = raw_input('What is the UT value? [12.0] ')
        if not UT: UT='12.0'
        src.updateheader(img,0,'UT',UT)
        hedUT='UT'
    try:
        _UT=str(getheader(img)[hedUT])
        check=1
        print 'UT= '+str(_UT)
    except:
        print 'Header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()


print aa
check=0
while check==0:
    hedJD=''
    print ' type [+] if keyword is missing '
    hedJD = raw_input('Enter the header keyword for JD [JD] ')
    if not hedJD: 
        hedJD='JD'
    elif hedJD =='+':
        hedJD='JD'
        check=1
        #try:
        #    print _UT,_date
        #    uh,um=string.split(str(_UT),':')
        #    if len(_date)==8:
        #        y=int(_date[0:4])
        #        m=int(_date[4:6])
        #        g=int(_date[6:8])
        #    else:
        #        if int(_date[4:6])>50:
        #            y=int('19'+str(_date[4:6]))
        #        else:
        #            y=int('20'+str(_date[4:6]))
        #        m=int(_date[2:4])
        #        g=int(_date[0:2])
        #    _JD=src.julday(y,m,g,uh,um)
        #    if str(_JD)[0:2]=='24':_JD=_JD-2400000
        #    print '#####  JD suggested = '+str(_JD)
        #except:
        #    _JD='0'
        #    print "WARNING: JD keyword has not been found in "+str(img)+" !!!!"
        #JD2 = raw_input('Enter the JD value ['+str(_JD)+'] ')
        #if not JD2: JD2=_JD
        #src.updateheader(img,0,'JD',JD2)
    if check==0:
        try:
            _JD=float(getheader(img)[hedJD])
            print 'JD= '+str(_JD)
            check=1
        except:
            print 'Header not found'
            again = raw_input('try again [y/n] [y] ? ')
            if not again: again='y'
            if again=='y': 
                check=0
            else: 
                sys.exit()


        
print aa
check=0
while check==0:
    hedfilter=''
    print ' type [+] if keyword is missing '
    hedfilter = raw_input('Enter the header keyword for the filter id[FILTER] ')
    if not hedfilter: hedfilter='FILTER'
    elif hedfilter=='+':
        _filter = raw_input('Select filter from this list [UBVRIugrizJHK] ? [V] ')
        if not _filter: _filter='V'
        src.updateheader(img,0,'FILTER',_filter)
        hedfilter='FILTER'
    try:
        _filter=str(getheader(img)[hedfilter])
        print 'FILTER= '+str(_filter)
        check=1
    except:
        print 'Header not found'
        again = raw_input('try again [y/n] [y] ? ')
        if not again: again='y'
        if again=='y': check=0
        else: 
            sys.exit()


print aa
check=0
while check==0:
    chose = raw_input('(RON) Header or Value [H/V] ')
    if chose=='H':
        hedron=''
        hedron = raw_input('Enter the header keyword for read-noise[RON] ')
	if not hedron: 
            hedron='RON'
        try:
            _ron=float(getheader(img)[hedron])
            check=1
            print 'RON= '+str(_ron)
        except:
            print 'Header not found'
            again = raw_input('try again [y/n] [y] ? ')
            if not again: again='y'
            if again=='y': check=0
            else: 
                sys.exit()
    elif chose=='V':
        _ron = raw_input('Enter the value for read-noise[1] ')
	if not _ron: _ron='1'
        check=1
    else:
        print '####'
        print "WARNING: Error in choice !!!"
	print '####'

print aa
check=0
while check==0:
    chose = raw_input('(GAIN) Header or Value [H/V] ')
    if chose=='H':
        hedgain = raw_input('Enter the header keyword for the gain [GAIN] ')
	if not hedgain: hedgain='GAIN'
        try:
            _gain=float(getheader(img)[hedgain])
            check=1
            print 'GAIN= '+str(_gain)
        except:
            print 'Header not found'
            again = raw_input('try again [y/n] [y] ? ')
            if not again: again='y'
            if again=='y': check=0
            else: 
                sys.exit()
    elif chose=='V':
        _gain=''
        _gain = raw_input('Enter the gain value[1] ')
	if not _gain: _gain='1'
        check=1
    else:
        print '####'
        print "WARNING: Error in choice !!!"
	print '####'

print aa
check=0
while check==0:
    chose = raw_input('DATAMIN (min CCD range) header or value [H/V] ')
    if chose=='H':
        hedmin = raw_input('Enter the header keyword for the datamin[DATAMIN] ')
	if not hedmin: hedmin='DATAMIN'
        try:
            _datamin=float(getheader(img)[hedmin])
            check=1
            print 'DATAMIN= '+str(_datamin)
        except:
            print 'Header not found'
            again = raw_input('try again [y/n] [y] ? ')
            if not again: again='y'
            if again=='y': check=0
            else: 
                sys.exit()
    elif chose=='V':
        _datamin = raw_input('Enter the datamin value [-100] ')
	if not _datamin: _datamin='-100'
        check=1
    else:
        print '####'
        print "WARNING: Error in choice!!"
	print '####'

print aa
check=0
while check==0:
    chose = raw_input('DATAMAX (max CCD range) header or value [H/V] ')
    if chose=='H':
        hedmax = raw_input('Enter the header keyword for datamax [DATAMAX] ')
        if not hedmax: hedmax='DATAMAX'
        try:
            _datamax=str(getheader(img)[hedmax])
            check=1
            print 'DATAMAX= '+str(_datamax)
        except:
            print 'header not found'
            again = raw_input('try again [y/n] [y] ? ')
            if not again: again='y'
            if again=='y': check=0
            else: 
                sys.exit()
    elif chose=='V':
        _datamax = raw_input('What is the datamax value [60000] ')
        if not _datamax: _datamax='60000'
        check=1
    else:
        print '####'
        print "WARNING: Error in choice !!!"
	print '####'

filterlist=[]
for img in imglist:
   if img:
       _filter=str(getheader(img)[hedfilter])
       if string.count(str(filterlist),_filter)==0:
           filterlist.append(_filter)
       else:
           print '####'
           print '### WARNING: more images with the same filter !!!! '
	   print '####'
   else:
    	print '####'
  	print '#### WARNING: empty space in the list !!'
	print '####'   	   

filterlist2=[]
for i in filterlist:
    filter2 = raw_input('Select the filter "'+i+'" amongst [UBVRIugrizJHK] ? ')
    filterlist2.append(filter2) 
        
for img in imglist:
    _airmass=float(getheader(img)[hedairmass])
    _exptime=float(getheader(img)[hedexptime])
    _object=str(getheader(img)[hedobject])
    _date=str(getheader(img)[heddate])
    _instrument=str(getheader(img)[hedinst])
    _UT=str(getheader(img)[hedUT])
    src.updateheader(img,0,'TELESCOP','other')
    src.updateheader(img,0,'QUBATELE','other')
    src.updateheader(img,0,'QUBAAIR',_airmass)
    src.updateheader(img,0,'QUBAEXP',_exptime)
    src.updateheader(img,0,'QUBADATE',_date)
    src.updateheader(img,0,'QUBAOB',_object)
    src.updateheader(img,0,'QUBAINS',_instrument)
    src.updateheader(img,0,'QUBAUT',_UT)
    src.updateheader(img,0,'QUBAGAIN',_gain)
    src.updateheader(img,0,'QUBARON',_ron)
    src.updateheader(img,0,'QUBADMAX',_datamax)
    src.updateheader(img,0,'QUBADMIN',_datamin)
    try:
        _JD=float(getheader(img)[hedJD])
        src.updateheader(img,0,'QUBAJD',_JD)
    except:
        print _UT,_date
        if string.count(_UT,':')==1:
            uh,um=string.split(str(_UT),':')
        elif string.count(_UT,':')==2:
            uh,um,us=string.split(str(_UT),':')
            um=float(um)+float(us)/60.
        else:
            try:
                uh=int(_UT)
                um=(float(_UT)-int(_UT))*60
            except:
                try:
                    ut=float(_UT)/3600.
                    uh=int(ut)
                    um=(float(ut)-int(ut))*60
                except:
                    print 'UT ', str(_UT),' not recognised !!!'
                    uh=float(raw_input('hour ?'))
                    um=float(raw_input('minutes ?'))
            try:
                int(_date)
                if len(_date)==8:
                    y=int(_date[0:4])
                    m=int(_date[4:6])
                    g=int(_date[6:8])
                else:
                    if int(_date[4:6])>50:
                        y=int('19'+str(_date[4:6]))
                    else:
                        y=int('20'+str(_date[4:6]))
                    m=int(_date[2:4])
                    g=int(_date[0:2])
            except:
                print 'date format ', str(_date),' not recognised !!!'
                y=float(raw_input('year ?'))
                m=float(raw_input('month ?'))
                g=float(raw_input('day ?'))

            _JD=src.julday(y,m,g,uh,um)
            if str(_JD)[0:2]=='24':_JD=_JD-2400000
            print '#####  JD suggested = '+str(_JD)
            src.updateheader(img,0,'QUBAJD',_JD)
            
    _filter=str(getheader(img)[hedfilter])
    for i in range(len(filterlist)):
        if filterlist[i]==_filter:
            _filter2=filterlist2[i]
            src.updateheader(img,0,'QUBAFIL',_filter2)
            break 
