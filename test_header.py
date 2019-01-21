#!/usr/bin/env python

from snoopy2 import src
import os,glob,shutil,sys


if len(sys.argv)<=1:
     listfits=glob.glob('*fits')
     for img in listfits:
          print img
     img=raw_input('Which image do you want to test ['+str(listfits[0])+'] ? ')
     if not img: img=listfits[0]
else:
     img=sys.argv[1]

_telescope=src.telescope(img)
print '## '
print _telescope
print '## '

_system=src.check_system(_telescope,img,Stdout=True)


if _system not in [0,1,2]:
     _header=src.read_parameter(_telescope,0)
else:
     _header=src.read_parameter(_telescope,_system)
     
_imagetype=src.imagetype(img,_header,_telescope)
_object=src.objects(img,_header,_telescope)
_JD=src.JD(img,_header,_telescope)
_airmass=src.airmass(img,_header,_telescope)
_filter=src.filter(img,_header,_telescope)
_grism=src.grism(img,_header,_telescope)
_exptime=src.exptime(img,_header,_telescope)
_date=src.date(img,_header,_telescope)
_gain=src.gain(img,_header,_telescope)
_ron=src.ron(img,_header,_telescope)
_lampid=src.lampid(img,_header,_telescope)
_RA=src.RA(img,_header,_telescope)
_DEC=src.DEC(img,_header,_telescope)
_ccdmax=src.ccdmax(img,_header,_telescope)
_ccdmin=src.ccdmin(img,_header,_telescope)
_cenwav=src.cenwav(img,_header,_telescope)
_slitw=src.slitw(img,_header,_telescope)
_UT=src.UT(img,_header,_telescope)
_xdimen=src.xdimen(img,_telescope)
_ydimen=src.ydimen(img,_telescope)
_instrument=src.instrument(img,_header,_telescope)
_obsmode=src.obsmode(img,_header,_telescope)

if not _gain: _gain='########'
if not _ron: _ron='########'
if not _instrument: _instrument=='########'
if not _ydimen: _ydimen='########'
if not _xdimen: _xdimen='########'
if not _filter: _filter='########'
if not _RA: _RA='########'
if not _grism: _grism='########'
if not _slitw: _slitw='########'
if not _lampid: _lampid='#######'
if not _date: _date='#######'
if not _cenwav: _cenwav='#######'
if not _UT: _UT='#######'
if not _ccdmin:_ccdmin='#######'
if not _ccdmax:_ccdmax='#######'
if not _obsmode:_obsmode='#######'
if not _object:_object='#######'
if _system not in [0,1,2]:_system='#######'

#_slitw='ddd'
print '####################################################################'
print 'IMG                OBJECT  IMAGETYPE    EXPTIME    FILTER        GRISM      '
print str(img)+'\t'+str(_object)+'\t'+str(_imagetype)+'\t'+str(_exptime)+'\t'+str(_filter)+'\t'+str(_grism)
print '####################################################################'
print 'AIRMASS             JD             DATE          XDIM   YDIM    GAIN   RON '
print str(_airmass)+'\t'+str(_JD)+'\t'+str(_date)+'\t'+str(_xdimen)+'\t'+str(_ydimen)+'\t'+str(_gain)+'\t'+str(_ron)
print '####################################################################'
print 'LAMP_ID     slitw         RA          DEC     CCDMIN    CCDMAX   CENWAV '
print str(_lampid)+'\t'+str(_slitw)+'\t'+str(_RA)+'\t'+str(_DEC)+'\t'+str(_ccdmin)+'\t'+str(_ccdmax)+'\t'+str(_cenwav)
print '####################################################################'
print ' UT     xdimension    ydimension        instrument      SYSTEM   OBSMODE '
print str(_UT)+'\t'+str(_xdimen)+'\t'+str(_ydimen)+'\t'+str(_instrument)+'\t'+str(_system)+'\t'+str(_obsmode)
print '####################################################################'

#src.close_program(logincl)

