def replace_directory(img):
    import snoopy2
    from snoopy2 import src
    import string,re
    if string.count(img,snoopy2.__path__[0])==1:
        img2=re.sub(snoopy2.__path__[0]+'/','home$',img)
    else:
        img2=img
    return img2

####  FILTER #####
def filter(img,_header,_telescope):
    import os,re,string
    import snoopy2
    from snoopy2 import src
    from astropy.io import fits as pyfits
    #img=src.replace_directory(img)
##############################  filter with pyfits  
    geth = pyfits.getheader(img)
    try:
        _filter=geth[_header['hed_filter1']]
    except:
        _filter=''
###############################
    if not _filter:
        if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
            try:
                _filter=src.esosecondheader(img,_telescope,'hed_filter1',_header['hed_filter1'])
            except:
                _filter=''
        elif _telescope=='WHT':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
        elif _telescope=='ACAM':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter1'])
            if _filter=='CLEAR':
                _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])

    if _telescope=='lp' and _filter=='clear':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
    if _telescope=='FTN' and _filter=='air':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
            if _filter=='air':
                _filter=pyfits.open(img)[0].header.get(_header['hed_filter3'])
    if _telescope=='mer' and _filter=='air':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
            if _filter=='air':
                _filter=pyfits.open(img)[0].header.get(_header['hed_filter3'])
    if _telescope=='FTS' and _filter=='air':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
            if _filter=='air':
                _filter=pyfits.open(img)[0].header.get(_header['hed_filter3'])
    if _telescope=='lsc' and _filter=='air':
            _filter=pyfits.open(img)[0].header.get(_header['hed_filter2'])
            if _filter=='air':
                _filter=pyfits.open(img)[0].header.get(_header['hed_filter3'])
    if _telescope=='PS1':
            _filter=re.sub('.00000', '', _filter)

    if _filter:
        _filter=re.sub('-', '', _filter)
        _filter=re.sub('"', '', _filter)
        _filter=re.sub('_', '', _filter)
        _filter=re.sub(' ', '', _filter)
        _filter=re.sub('\#', '', _filter)
        _filter=re.sub('\*', '', _filter)
        _filter=re.sub('\t', '', _filter)
    return _filter


def filtername(_telescope,_filter,_system):
    filtro='unknown'
    if _system==0:
      filter_name={}
      filter_name['reference']=['U','B','V','R','I']
      filter_name['lp']=['SDSSU','BessellB','BessellV','SDSSR','SDSSI']
      filter_name['FTN']=['SDSSU','BessellB','BessellV','BessellR','BessellI']
      filter_name['mer']=['SDSSU','BessellB','BessellV','BessellR','BessellI']
      filter_name['FTS']=['SDSSU','BessellB','BessellV','BessellR','BessellI']
      filter_name['lsc']=['U','B','V','R','I']
      filter_name['ekar']=['U','B','V','R','i']
      filter_name['sampur']=['JU','JB','JV','JR','JI']
      filter_name['sch']=['U','B','V','R','I']
      filter_name['ctio']=['U','B','V','R','I']
      filter_name['montsec']=['U','B','V','R','I']
      filter_name['dutch']=['U','B','V','R','I']
      filter_name['danish']=['U','B','V','R','I']
      filter_name['trp']=['U','B','V','Rc','I']
      filter_name['hct']=['','6BesB','5BesV','4BesR','3BesI']
      filter_name['NOT']=['UBes36260','BBes440100','VBes53080','RBes650130','iint797157']
      filter_name['TNG']=['UJohn01','BJohn10','VJohn11','RJohn12','IJohn13']
      filter_name['NTT']=['U640','B639','V641','R642','i705']
      filter_name['FORS1']=['USPECIAL','BBESS','vHIGH','RBESS','IBESS']
      filter_name['FORS2']=['USPECIAL','BBESS','VBESS','RSPECIAL','IBESS']
                          #['UHIGH','bHIGH','vHIGH','RSPECIAL','IBESS']
      filter_name['WHT']=['?','harb1','harv3','harr4','?']
      filter_name['ACAM']=['?','BESB','HARV3','BESR','BESI']
      filter_name['ca']=['JohnU','JohnB','JohnV','JohnR','JohnI']
      filter_name['pennar']=['','','','','']
      filter_name['wise']=['U','B','V','R','I']
      filter_name['sub']=['?','?','SCFCFLBV01','SCFCFLBR01','SCFCFLBI01']
      filter_name['lo']=['2','3','4','5','6']
      filter_name['wend']=['UBess','Bold','Vold','Rnew(OG570+Cfx)','Inew(RG780)']
      filter_name['kait']=['u','b','v','r','i']
      filter_name['s70']=['u','b','v','r','i']
      filter_name['rem']=['u','b','v','r','i']
      filter_name['other']=['U','B','V','R','I']
      filter_name['prompt']=['U','B','V','R','I']
      filter_name['Meckering']=['U','B','V','R','I']
      filter_name['swift']=['U','B','V','','']
      if _telescope in filter_name:
        for i in range(len(filter_name[_telescope])):
 	  if filter_name[_telescope][i]==_filter:   
              filtro=filter_name['reference'][i]
    elif _system==2:
      filter_name={}
      filter_name['reference']=['u','g','r','i','z']
      filter_name['lp']=['SDSSU','SDSSG','SDSSR','SDSSI','SDSSZ']
      filter_name['mer']=['SDSSU','SDSSG','SDSSR','SDSSI','PanStarrsZ']
      filter_name['FTN']=['SDSSU','SDSSG','SDSSR','SDSSI','PanStarrsZ']
      filter_name['FTS']=['SDSSU','SDSSG','SDSSR','SDSSI','PanStarrsZ']
      filter_name['lsc']=['up','gp','rp','ip','zs']
      filter_name['NTT']=['','g782','r784','i705','z623']
      filter_name['ekar']=['','','','','']
      filter_name['dutch']=['','','','','']
      filter_name['sch']=['U','V','R','I','']
      filter_name['NOT']=['','','','','']
      filter_name['TNG']=['','','','','']
      filter_name['FORS1']=['','','','','']
      filter_name['trp']=['','','','','']
      filter_name['FORS2']=['uHIGH','VBESS','RSPECIAL','IBESS','zGUUN']
      filter_name['WHT']=['?','?','?','?','?']
      filter_name['ACAM']=['?','?','?','?','?']
      filter_name['ca']=['u','g','r','i','z']
      filter_name['pennar']=['','','','','']
      filter_name['wise']=['u','g','r','i','z']
      filter_name['sampur']=['u','g','r','i','z']
      filter_name['sub']=['?','?','?','?','?']
      filter_name['lo']=['u','g','r','i','z']
      filter_name['wend']=['u','g','r','i','z']
      filter_name['kait']=['u','g','r','i','z']
      filter_name['s70']=['u','g','r','i','z']
      filter_name['mag']=['uSloan','gSloan','rSloan','iSloan','zSloan']
      filter_name['rem']=['u','g','r','i','z']
      filter_name['PS1']=['u','g','r','i','z']
      filter_name['other']=['u','g','r','i','z']
      filter_name['prompt']=['uprime','gprime','rprime','iprime','zprime']
      filter_name['Meckering']=['uprime','gprime','rprime','iprime','zprime']
      if _telescope in filter_name:
        for i in range(len(filter_name[_telescope])):
          if filter_name[_telescope][i]==_filter:   
              filtro=filter_name['reference'][i] 
    elif _system==1:
      filter_name={}
      filter_name['reference']=['J','H','K']
      filter_name['lp']=['BarrJ','BarrH','']
      filter_name['FTN']=['BarrJ','BarrH','']
      filter_name['ekar']=['','','']
      filter_name['NOT']=['J','H','Ks']
      filter_name['TNG']=['J','H','Kprim']
      filter_name['NICS']=['J','H','Kprim']
      filter_name['NTT']=['J','H','Ks']
      filter_name['SOFI']=['J','H','Ks']
      filter_name['trp']=['','','','','']
      filter_name['CI']=['J','H','K']
      filter_name['FORS1']=['','','']
      filter_name['FORS2']=['J','H','Ks']
      filter_name['WHT']=['?','?','?']
      filter_name['ca']=['','','']
      filter_name['pennar']=['','','']
      filter_name['wise']=['','','']
      filter_name['sub']=['?','?','?']
      filter_name['lo']=['2','3','4']
      filter_name['wend']=['','','']
      filter_name['kait']=['','','']
      filter_name['s70']=['','','']
      filter_name['rem']=['J','H','K']
      filter_name['other']=['J','H','K']
      if _telescope in filter_name:
        for i in range(len(filter_name[_telescope])):
          if filter_name[_telescope][i]==_filter:   
              filtro=filter_name['reference'][i]    	         
    else:
        filtro='unknown'
    return filtro


def check_system(_telescope,img,Stdout=True):
    from snoopy2 import src
    try:
        import pyfits
    except:
        from astropy.io import fits as pyfits

    import re
    if pyfits.getheader(img).get('INSTRUME')=='SOFI' or pyfits.getheader(img).get('INSTRUME')=='NICS':
        _system=1
    else:
      aa='unknown'
      _header=src.read_parameter(_telescope,0)
      if _header:
        _filter=src.filter(img,_header,_telescope)
        aa=src.filtername(_telescope,_filter,0)
      if aa!='unknown':
        if Stdout:
            print '##################################'
            print 'System optical: UBVRI'
            print '##################################'
        _system=0
      if aa=='unknown':
        _header=src.read_parameter(_telescope,1)
        if _header:
            _filter=src.filter(img,_header,_telescope)
            aa=src.filtername(_telescope,_filter,1)
        if aa!='unknown':
            if Stdout:
                print '##################################'
                print 'System infrared: JHK'
                print '##################################'	 
            _system=1
      if aa=='unknown':
        _header=src.read_parameter(_telescope,2)
        if _header:
            _filter=src.filter(img,_header,_telescope)	   
            aa=src.filtername(_telescope,_filter,2)
        if aa!='unknown':
            if Stdout:
                print '##################################'
                print 'System sloan: ugriz'
                print '##################################'	 
            _system=2
      if aa=='unknown':
        if Stdout:
            print '##################################'
            print 'WARNING: FILTER not recognised in '+str(img)+' !!!  '
            print '##################################'      
        _system=99		        
       #      sys.exit()
    return _system

##### DATE #######
def date(img,_header,_telescope):
    import re,string
    from snoopy2 import src
    from astropy.io import fits as pyfits 
    geth = pyfits.getheader(img)
    img=src.replace_directory(img)
    _date=''
    try:
        _date=geth[_header['hed_date']]
        _date = re.sub('-', '', _date)
        _date = re.sub('/', '', _date)
        if not string.find(_date,'T')==-1:
            _date = _date[0:string.find(_date,'T')]
        if _telescope=='kait':
            print '###telescope= '+_telescope
            _date=_date[4:8]+_date[2:4]+_date[0:2]
        if len(_date)==6:
            if int(_date[4:6])>50:
                y=str('19'+str(_date[4:6]))
            else:
                y=str('20'+str(_date[4:6]))
            m=str(_date[2:4])
            g=str(_date[0:2])
            _date=str(y)+str(m)+str(g)
    except:
        try: 
            JD=src.JD(img,_header,_telescope)
            month,day,year,hour,_min=jd2ut(JD+2400000)
            if len(str(month))==1: month='0'+str(month)
            if len(str(day))==1: day='0'+str(day)

            _date=str(year)+str(month)+str(day)
            print month,day,year,hour,_min
            src.updateheader(img,0,_header['hed_date'],_date)
        except:
            print "WARNING: DATE keyword has not been found in "+str(img)+" !!!!"
#           sys.exit()
    return _date

def jd2ut(jd):
   import math
   ijd = 0
   fjd = 0.0
   a = 0
   b = 0
   c = 0
   d = 0
   e = 0
   g = 0
   day = 0.0
   month = 0
   year = 0

   jd += 0.5
   ijd = math.floor(jd)
   fjd = jd - ijd
   if ijd > 2299160:   # Oct. 15, 1582
      a = math.floor((ijd - 1867216.25) / 36524.25)
      b = ijd + 1 + a - math.floor(a / 4.0)
   else: 
      b = ijd
   
   c = b + 1524
   d = math.floor((c - 122.1) / 365.25)
   e = math.floor(d * 365.25)
   g = math.floor((c - e) / 30.6001)
   day = c - e + fjd - math.floor(g * 30.6001)
   if  g < 13.5:
         month = g - 1
   else: 
      month = g - 13
   
   if  month > 2.5: 
      year = d - 4716
   else: 
      year = d - 4715
   hour=math.floor((day-math.floor(day))*24)
   _min=(((day-math.floor(day))*24)-math.floor((day-math.floor(day))*24))*60
   day=int(day)
   year=int(year)
   month=int(month)
   return month,day,year,hour,_min

###### JD ######
def julday(year,month,day,hour,min):
   import math
   univTime = float(hour)+(float(min)/60.)
   if ((100*float(year))+float(month)-190002.5) >= 0:
      sign = 1
   else:
      sign = -1
   part1 = 367 * float(year)
   part2 = math.floor((7*(float(year)+math.floor((float(month)+9)/12)))/4)
   part3 = float(day)+math.floor((275*float(month))/9)
   part4 = 1721013.5+(univTime/24)
   part5 = 0.5*sign
   jd = part1-part2+part3+part4-part5+0.5
   return jd

def JD(img,_header,_telescope):
    import string
    from snoopy2 import src
    from astropy.io import fits as pyfits 
    #img=src.replace_directory(img)

    geth = pyfits.getheader(img)
    _system=src.check_system(_telescope,img,Stdout=False)
    if _system not in [0,1,2]: _system=0
    try:
        if _telescope=='wise':
            try:
                _JD=float(geth[_header['hed_JulianD']])
            except:
                _JD=float(geth[_header['hed_JulianD2']])
        else:    
            _JD=float(geth[_header['hed_JulianD']])
        if _telescope=='lp': _JD=_JD+0.5
        if _telescope=='FTN': _JD=_JD+0.5
        if _telescope=='FTS': _JD=_JD+0.5
        if _telescope=='FORS1' or _telescope=='FORS2': _JD=_JD+0.5        
        if _telescope=='CI': _JD=_JD+0.5
        if _telescope=='ps1': _JD=_JD+0.5
        if _telescope=='dutch': _JD=_JD+0.5
        if _telescope=='danish': _JD=_JD+0.5
	if _telescope=='prompt' or _telescope=='Meckering' or _telescope=='NTT' : _JD=_JD+0.5        
        if str(_JD)[0:2]=='24':_JD=_JD-2400000
        if _telescope=='NICS' and _system==1: _JD=_JD+50000
    except:
        try:
            _UT=UT(img,_header,_telescope)
            uh,um=string.split(str(_UT),':')
            _date=date(img,_header,_telescope)
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
            _JD=julday(y,m,g,uh,um)
            if str(_JD)[0:2]=='24':_JD=_JD-2400000
            updateheader(img,0,_header['hed_JulianD'],_JD)
            #iraf.hedit(img,_header['hed_JulianD'],_JD,add='yes',update='yes',verify='no')
        except:
            _JD=''
            print "WARNING: JD keyword has not been found in "+str(img)+" !!!!"
            #sys.exit()
    return _JD


def UT(img,_header,_telescope):
    from snoopy2 import src
    from astropy.io import fits as pyfits
    import string,re
    img = src.replace_directory(img)

    geth = pyfits.getheader(img)
    try:
        _UT=str(geth[_header['hed_UT']])
        _UT=re.sub(' ','',_UT)
        _UT=re.sub('"','',_UT)
        _UT=re.sub("'","",_UT)
        if not _UT:
            if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
                _UT=src.esosecondheader(img,_telescope,'hed_UT',_header['hed_UT'])
	if _telescope=='prompt':_UT = string.split(_UT,'T')[1]
        if _telescope=='Meckering':_UT = string.split(_UT,'T')[1]
        if _telescope=='PS1':_UT = string.split(_UT,'T')[1]
        if _telescope=='wise':_UT = string.split(_UT,'T')[1]
        if _telescope=='trp':_UT = string.split(_UT,'T')[1]
        if _telescope=='swift':_UT = string.split(_UT,'T')[1]
        if _telescope=='sch':_UT = string.split(_UT,'T')[1]
        if string.count(_UT,':')==1:
            uh=int(string.split(_UT,':')[0])
            um=(float(string.split(_UT,':')[1]))
        elif string.count(_UT,':')==2:
            uh=int(string.split(_UT,':')[0])
            um=(float(string.split(_UT,':')[1])+(float(string.split(_UT,':')[2])/60.))
            _UT=str(uh)+':'+str(um)
        elif string.count(_UT,':')==0:
	    if int(float(_UT))>=24: _UT=float(_UT)/3600. 
            uh=int(float(_UT))
            um=(float(_UT)-int(float(_UT)))*60.
            _UT=str(uh)+':'+str(um)
    except:
        _UT=''
        print "WARNING: UT keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _UT

def objects(img,_header,_telescope):
    _object=''
    import re
    from snoopy2 import src
    from astropy.io import fits as pyfits
    img=src.replace_directory(img)

    geth = pyfits.getheader(img)

    try:
        _object=geth[_header['hed_object']]
        if not _object:
            if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
                _object=src.esosecondheader(img,_telescope,'hed_object',_header['hed_object'])
        _object = re.sub("'", "", _object)
        _object = re.sub(" ", "", _object)
        _object = re.sub('"', '', _object)
        _object = re.sub('/', '_', _object)
        _object = re.sub("\+", "_", _object)
        _object = re.sub("\.", "_", _object)
        _object = re.sub("\,", "_", _object)
        _object = re.sub("\#", "", _object)
    except:
        print "WARNING: object keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _object
    
def imagetype(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    import re

    from astropy.io import fits as pyfits

    geth = pyfits.getheader(img)
    try:
        _imagetype=geth[_header['hed_imagetyp']]
    except:
        _imagetype=''
    
    if not _imagetype:
            if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
                try:
                    _imagetype=src.esosecondheader(img,_telescope,'hed_imagetyp',_header['hed_imagetyp'])
                except:
                    _imagetype=''
    if _imagetype:
        _imagetype = re.sub("'", "", _imagetype)
        _imagetype = re.sub('"', '', _imagetype)
    else:
        print "WARNING: imagetype keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _imagetype

def exptime(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth=pyfits.getheader(img)

    try:
        _exptime=float(geth[_header['hed_exptime']])
    except:
        _exptime=''
        print "WARNING: exposure time keyword has not been found in "+str(img)+" !!!!"
    return _exptime
    
def cenwav(img,_header,telescope):
    from snoopy2 import src
    img=src.replace_directory(img)

    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)

    _cenwav=''
    try:
        _cenwav=float(geth[_header['hed_cenwav']])
    except:
        print "WARNING: central wavlength keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _cenwav
    
def slitw(img,_header,_telescope):
    import snoopy2
    from snoopy2 import src
    img=src.replace_directory(img)
    import string,re

    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)

    _slitw=''
    if _telescope=='ca':
        _slitw='slit'
    else:
      try:
        _slitw=float(geth[_header['hed_slitw']])
	if _telescope=='WHT': _slitw='%.1f' % (_slitw) #_slitw=((round_(float(_slitw)*2,0))/2)
	if _telescope=='ekar': _slitw='%.2f' % (_slitw)
      except:
          try:
              _slitw=geth[_header['hed_slitw']]
          except:
              _slitw=''
      if not _slitw:
          if _telescope=='NTT' or _telescope=='FORS1' or _telescope=='FORS2':
              try:
                  _slitw=src.esosecondheader(img,_telescope,'hed_slitw',_header['hed_slitw'])
              except:
                  _slitw=''
      if _slitw:
        if _telescope=='TNG':
            try:
                for i in string.split(_slitw,'_'):
                    try: tmpslit=float(i)
                    except: pass
                    _slitw=str(tmpslit)+'_arcsec'
            except: pass
        _slitw=re.sub('"', '', _slitw)
        _slitw=re.sub('#', '', _slitw)
      else:
            print "WARNING: slit width keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _slitw
 
def xdimen(img,_telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    _xdimen=''
    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)
    try:
        _xdimen=float(geth['NAXIS1'])
    except:
        print "WARNING: XDIMEN keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _xdimen

def ydimen(img,telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    _ydimen=''
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _ydimen=float(geth['NAXIS2'])
    except:
        print "WARNING: YDIMEN keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _ydimen

def ron(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)

    try:
        _ron=float(_header['RON_value'])
    except:
        try:
            _ron=float(geth[_header['RON_value']])
        except:
            _ron=0
            print "WARNING: RON keyword has not been found in "+str(img)+" !!!!"
    return _ron

def gain(img,_header,_telescope):
    print '#######',img
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)
    try:
        _gain=float(_header['GAIN_value'])
    except:
        try:
            _gain=float(geth[_header['GAIN_value']])
        except:
            _gain=[]
            print "WARNING: GAIN keyword has not been found in "+str(img)+" !!!!"
    return _gain

def ccdmin(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)
    try:
        _ccdmin=float(_header['datamin_value'])
    except:
        try:
            _ccdmin=float(geth[_header['datamin_value']])
        except:
            _ccdmin=[]
            print "WARNING: datamin keyword has not been found in "+str(img)+" !!!!"
    return _ccdmin

def ccdmax(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)
    try:
        _ccdmax=float(_header['datamax_value'])
    except:
        try:
            _ccdmax=float(geth[_header['datamax_value']])
        except:
            _ccdmax=[]
            print "WARNING: datamin keyword has not been found in "+str(img)+" !!!!"
    return _ccdmax

def grism(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    import string,re
    from snoopy2 import src
    import snoopy2
    from astropy.io import fits as pyfits
    geth= pyfits.getheader(img)
    from time import sleep    
    rec_grism2=['R158R','R300B','R316R','R1200R','R600B','LR-B','LR-R','VHR-R','VHR-V','Grism_4','Grism_8','Grism_7','Grism_5','GR02','GR04','gr02','gr04','VPH4','gr08','GR08','Gr_4Hori', 'green-200', 'red-200', 'blue-200', 'blue-100', 'red-100']
    _grism=''
    try:
        _grism=geth[_header['hed_grism']]
    except:
        _grism=''
    if not _grism:
            if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
                try:
                    _grism=src.esosecondheader(img,_telescope,'hed_grism',_header['hed_grism'])
                except:
                    _grism=''
    if _grism:
	if _telescope=='ekar':
            _filter=filter(img,_header,_telescope)
            if string.count(rec_grism2,_filter):
                _grism=_filter
	if _grism[0]=='"' and _grism[-1]=='"':
  	        _grism=eval(_grism)
	_grism=re.sub('"', '', _grism)
	_grism=re.sub("'", "", _grism)
	_grism=re.sub("#", "", _grism)
	_grism=re.sub(" ", "", _grism)
    else:
        print "WARNING: GRISM keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _grism

def instrument(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    import string,re
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _instrum=geth[_header['hed_instrume']]
        _instrum = re.sub("'", "", _instrum)
    except:
        _instrum=''
        print "WARNING: INSTRUMENT keyword has not been found in "+str(img)+" !!!!"
    return _instrum
    
def dichroic(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _dichroic=geth[_header['hed_dichroic']]
    except:
        print "WARNING: Dichroic keyword has not been found !!!!"
    return _dichroic

def xbin(img,_header):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _xbin=geth[_header['hed_xbin']]
        _xbin=int(_xbin)
    except:
        _xbin=0
        #print "WARNING: X Binning keyword has not been found !!!!"
        #sys.exit()
    return _xbin

def ybin(img,_header):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _ybin=geth[_header['hed_ybin']]
        _ybin=int(_ybin)
    except:
        _ybin=0
        #print "WARNING: Y Binning keyword has not been found !!!!"
        #sys.exit()
    return _ybin
   
def biassec(img,header,telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _biassec=geth[_header['hed_biassec']]
    except:
        print "WARNING: biassec keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _biassec
    
def trimsec(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    if _telescope == 'NOT': 
        try:
	  _trimsec=_header['hed_trimsec']
	except:  
	  print "WARNING: trimsec keyword has not been found in "+str(img)+" !!!!"
	  _trimsec=[]
         #sys.exit()
    else:
       try:
           _trimsec=geth[_header['hed_trimsec']]
       except:  
         print "WARNING: trimsec keyword has not been found in "+str(img)+" !!!!"
	 _trimsec=[]
         #sys.exit()    
    return _trimsec
    
def obsmode(img,_header,_telescope):
    from snoopy2 import src
    img=src.replace_directory(img)
    _obsmode=''
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _obsmode=geth[_header['hed_OBSMODE']]
    except:
        print "WARNING: obsmode keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _obsmode
    
def lampid(img,_header,_telescope):
    _lampid=''
    from snoopy2 import src
    #img=src.replace_directory(img)
    import string,re
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    rec_HGCD=['HgCd','hgcd','hg','cd','Hg','HG','CD','CD']
    rec_NE=['Ne','NE','ne','Neon','NEON']
    if _telescope=='NOT':
       try:
            _lampid=geth[_header['hed_lampid']]
       except:
            _lampid=''
       if not _lampid:
         try:
           clamp={}
           clamp[src.CLAMPNM1(img,_header,_telescope)]=src.CLAMP1(img,_header,_telescope)
           clamp[src.CLAMPNM2(img,_header,_telescope)]=src.CLAMP2(img,_header,_telescope)
           clamp[src.CLAMPNM3(img,_header,_telescope)]=src.CLAMP3(img,_header,_telescope)
           clamp[src.CLAMPNM4(img,_header,_telescope)]=src.CLAMP4(img,_header,_telescope)
           for i in clamp:
              if clamp[i]!=0 and i!='Halogen':
                   _lampid=_lampid+i
           updateheader(img,0,_header['hed_lampid'],_lampid)
         except:
             updateheader(img,0,_header['hed_lampid'],'XXXXXX')
    elif _telescope=='ekar':
         try:
             _lampid=geth[_header['hed_lampid']]
         except:
             _lampid=''
         if not _lampid:
	     _object=src.objects(img,_header,_telescope)
	     for i in rec_HGCD:
	        if string.count(_object,i):
	            _lampid='HgCd'
	     for i in rec_NE:
	        if string.count(_object,i):
	            _lampid='Ne'
             if not _lampid:
	         print 'Warning: no lamp reconized for '+str(img)+' !!!'
	         _imagetype=src.imagetype(img,_header,_telescope)
	         print _object,_imagetype
                 try:
                     iraf.display(img,1)
                 except:
                     pass
                 _lampid=raw_input('Write which is the lamp [Ne],[HgCd] or press enter. ')
                 if not _lampid:_lampid='xxx'
             updateheader(img,0,_header['hed_lampid'],_lampid)
    elif _telescope=='NTT':
        _lampid=''
    else:
       try:
           _lampid=geth[_header['hed_lampid']]
       except:
          print "WARNING: lampid keyword has not been found in "+str(img)+" !!!!"    
    _lampid=re.sub('/', '', _lampid)
    _lampid=re.sub(' ', '', _lampid)
    _lampid=re.sub('"', '', _lampid)
    return _lampid

def catlog(img,_header,_telescope):
    from snoopy2 import src
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    img=src.replace_directory(img)
    try:
        _catlog=geth[_header['hed_catlog']]
    except:
        print "WARNING: catalogue keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _catlog

def RA(img,_header,_telescope):
    import string,re
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    _RA=''
    try:
       _RA=float(geth[_header['hed_RA']])
       if _RA:
           if _telescope=='NOT' or _telescope=='NTT' or _telescope=='FORS1' or _telescope=='FORS2' \
                   or _telescope=='ca' or _telescope=='rem' or _telescope=='swift' or _telescope=='dutch' \
                   or _telescope=='TNG' or _telescope=='SOFI':
               print 'change RA(degree) in RA(h)'
               _RA=_RA/15. 
    except:
       try:
          _RA=geth[_header['hed_RA']]
          if string.count(_RA,'"'): _RA=eval(_RA)
          if _telescope=='trp':     _RA=string.join(string.split(_RA),':')
	  if string.count(_RA,':')==2:
              if _telescope=='PS1':
                  _RA=((float(string.split(_RA,':')[2])/60+float(string.split(_RA,':')[1]))/60)+(float(string.split(_RA,':')[0])/15.)
              else:
                  _RA=((float(string.split(_RA,':')[2])/60+float(string.split(_RA,':')[1]))/60)+float(string.split(_RA,':')[0])
          
       except:
           print "WARNING: RA keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _RA
    
def DEC(img,_header,_telescope):
    import string
    from snoopy2 import src
    img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    _DEC=''
    try:
        _DEC=float(geth[_header['hed_DEC']])
    except:
       try: 
          _DEC=geth[_header['hed_DEC']]
          if string.count(_DEC,'"'): _DEC=eval(_DEC)
          if _telescope=='trp':     _DEC=string.join(string.split(_DEC),':')
	  if string.count(_DEC,':')==2:
            if string.count(str(string.split(_DEC,':')[0]),'-')==0:
          	    _DEC=((float(string.split(_DEC,':')[2])/60+float(string.split(_DEC,':')[1]))/60)+float(string.split(_DEC,':')[0])
	    else:
        	    _DEC=(-1)*(((abs(float(string.split(_DEC,':')[2])/60)+float(string.split(_DEC,':')[1]))/60)+abs(float(string.split(_DEC,':')[0])))
       except:
        print "WARNING: DEC keyword has not been found in "+str(img)+" !!!!"
        #sys.exit()
    return _DEC    
    
    
def CLAMP1(img,_header,_telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clamp1=int(geth[_header['hed_clamp1']])
    except:
         _clamp1=[]
         print "WARNING: CLAMP1 keyword has not been found in "+str(img)+" !!!!"
    return _clamp1    
    
def CLAMP2(img,_header,_telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clamp2=int(geth[_header['hed_clamp2']])
    except:
         _clamp2=[]
         print "WARNING: CLAMP2 keyword has not been found in "+str(img)+" !!!!"
    return _clamp2    
    
def CLAMP3(img,_header,_telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clamp3=int(geth[_header['hed_clamp3']])
    except:
        _clamp3=[]
        print "WARNING: CLAMP3 keyword has not been found in "+str(img)+" !!!!"
    return _clamp3    
    
def CLAMP4(img,_header,_telescope):
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clamp4=int(geth[_header['hed_clamp4']])
    except:
            _clamp4=[]
            print "WARNING: CLAMP4 keyword has not been found in "+str(img)+" !!!!"
    return _clamp4        
    


def CLAMPNM1(img,_header,_telescope):
    import re
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clampnm1=geth[_header['hed_clampnm1']]
	_clampnm1=re.sub(' ', '', _clampnm1)
        if _clampnm1[0]=='"' and _clampnm1[-1]=='"':
  	       _clampnm1=eval(_clampnm1)
    except:
        _clampnm1=[]
        print "WARNING: CLAMPNM1 keyword has not been found in "+str(img)+" !!!!"
    return _clampnm1    
    
def CLAMPNM2(img,_header,_telescope):
    import re
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
         _clampnm2=geth[_header['hed_clampnm2']]
	 _clampnm2=re.sub(' ', '', _clampnm2)
	 if _clampnm2[0]=='"' and _clampnm2[-1]=='"':
  	        _clampnm2=eval(_clampnm2)	 
    except:
         _clampnm2=[]
         print "WARNING: CLAMPNM2 keyword has not been found in "+str(img)+" !!!!"
    return _clampnm2    
    
def CLAMPNM3(img,_header,_telescope):
    import re
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
         _clampnm3=geth[_header['hed_clampnm3']]
	 _clampnm3=re.sub(' ', '', _clampnm3)
	 if _clampnm3[0]=='"' and _clampnm3[-1]=='"':
  	       _clampnm3=eval(_clampnm3)
    except:
            _clampnm3=[]
            print "WARNING: CLAMPNM3 keyword has not been found in "+str(img)+" !!!!"
    return _clampnm3    
    
def CLAMPNM4(img,_header,_telescope):
    import re
    from snoopy2 import src
    #img=src.replace_directory(img)
    from astropy.io import fits as pyfits
    geth = pyfits.getheader(img)
    try:
        _clampnm4=geth[_header['hed_clampnm4']]
        _clampnm4=re.sub(' ', '', _clampnm4)
        if _clampnm4[0]=='"' and _clampnm4[-1]=='"':
  	        _clampnm4=eval(_clampnm4)
    except:
            _clampnm4=[]
            print "WARNING: CLAMPNM4 keyword has not been found in "+str(img)+" !!!!"
    return _clampnm4        

#############################################################################################

def airmass(img,_header,_telescope):
  import string,os,sys
  import snoopy2
  from snoopy2 import src
  img=src.replace_directory(img)
  from astropy.io import fits as pyfits
  geth = pyfits.getheader(img)
  _checkairmass=1
  _airmass=''
  try: 
    _airmass=float(geth[_header['hed_airmass']])
  except:
      pass

  if not _airmass:
        if _telescope=='FORS1' or _telescope=='FORS2' or _telescope=='NTT':
            try:
                _airmass=src.esosecondheader(img,_telescope,'hed_airmass',_header['hed_airmass'])
            except:
                _airmass=''
  if not _airmass:  #or _airmass>=99:
       print 'WARNING: Airmass header keyword not found !!'
       _airmass=0.00001

       _RA=src.RA(img,_header,_telescope)
       _DEC=src.DEC(img,_header,_telescope)  
       _UT=src.UT(img,_header,_telescope)
       _date=src.date(img,_header,_telescope)
       if not _RA or not _DEC:
           _checkairmass=0
       try:
           _UT2=string.split(_UT,':')[0]+':'+string.split(string.split(_UT,':')[1],'.')[0]+':'+str(int((float(string.split(_UT,':')[1])-float(string.split(string.split(_UT,':')[1],'.')[0]))*60))
           _date2=_date[0:4]+'-'+_date[4:6]+'-'+_date[6:8]
       except:
           print 'WARNING: UT or DATE keywords not found in header, airmmass check not possible.'
           _checkairmass=0

       if _telescope=='WHT' or _telescope=='TNG' or _telescope=='NOT' or _telescope=='lp':	
           _observatory='lapalma'
       elif _telescope=='ekar':
           _observatory='ekar'
       elif _telescope=='NTT' or _telescope=='dutch' or _telescope=='danish':
           _observatory='lasilla'
       elif _telescope=='FORS1' or _telescope=='FORS2':
           _observatory='paranal'
       elif _telescope=='montsec':
           _observatory='montsec'
       else:
           print 'WARNING: observatory not found, airmass check not possible.'
           _observatory=''
           _checkairmass=0

       if _observatory and _checkairmass:
           import pyraf
           from pyraf import iraf
           from iraf import astutil
           f = file('airmass.txt','w')
           f.write('mst = mst ("'+str(_date2)+'",'+str(_UT2)+', obsdb ("'+str(_observatory)+'", "longitude"))\n')
           f.write('air = airmass ('+str(_RA)+','+str(_DEC)+',mst, obsdb ("'+str(_observatory)+'", "latitude"))\n')
           f.write('print(air)\n')
           f.close()

           _air=iraf.astcalc(image=img, command="airmass.txt",Stdout=1)[0]
           if ((float(_air)-float(_airmass))**2)>0.01:
               if _telescope=='ekar' or _telescope=='dutch' or _telescope=='danish' or _telescope=='montsec':
                   src.updateheader(img,0,_header['hed_airmass'],_air)
                   #iraf.hedit(img,_header['hed_airmass'],_air,add='yes',update='yes',verify='no')
           else:
               print 'Airmass value: ok '
               _airmass=float(geth[_header['hed_airmass']])
       else:
           print 'WARNING: airmass check not done '    

  output=_airmass
  return output

###############################################################################################

def sn_coordinate(imgl,_interactive):
       from math import sqrt
       from snoopy2 import src
       import snoopy2
       from pyraf import iraf
       import string
       import os
       imgl=src.replace_directory(imgl)
       if imgl[-5:]!='.fits':
           imgllong=imgl+'.fits'
       else:
           imgllong=imgl
       listastandard=snoopy2.__path__[0]+'/standard/fluxstandard/supernovaelist.txt'
       f=open(listastandard,'r')
       liststd=f.readlines()
       f.close()
       star,ra,dec=[],[],[]
       for i in liststd:
          vector=string.split(i)
          star.append(vector[0])
          #print vector
	  ra.append(float(vector[1])+((float(vector[2])+(float(vector[3])/60.))/60.))
          if float(vector[4])>0:
         	  dec.append(float(vector[4])+((float(vector[5])+(float(vector[6])/60.))/60.))	  
          else:
                  aa= -1*((abs(float(vector[4])))+((float(vector[5])+(float(vector[6])/60.))/60.))	  
               	  dec.append(aa)
	  print ra[-1],dec[-1],star[-1]
       _telescope=src.telescope(imgllong)
       _system=src.check_system(_telescope,imgllong,Stdout=True)
       _header=src.read_parameter(_telescope,_system)
       if _telescope=='other':
           print 'other'
           _xdimen=src.xdimen(imgllong,_telescope)
           _ydimen=src.ydimen(imgllong,_telescope)
           xcenter=int(float(_xdimen)/2)
           ycenter=int(float(_ydimen)/2)
           f=open('_tmp.tv','w')
           f.write(str(xcenter)+'  '+str(ycenter)+'\n')
           f.close()
           iraf.delete('tmp.coo')
           iraf.wcsctran('_tmp.tv','tmp.coo',imgl,inwcs='logical',units='degrees degrees',outwcs='world',columns='1 2',formats='%10.8f %10.8f')
           sss=iraf.fields('tmp.coo','1,2', Stdout=1)
           _RA,_DEC=string.split(sss[2])
           _RA=float(_RA)/15.
           _DEC=float(_DEC)
       else:
           _RA=src.RA(imgllong,_header,_telescope)
           _DEC=src.DEC(imgllong,_header,_telescope)
       print _RA,_DEC
       ra0=''
       if _RA and _DEC:
          for st in range(len(star)):
              distance=sqrt((float(_RA)-float(ra[st]))**2+(float(_DEC)-float(dec[st]))**2)
              print st,distance
              if distance<=0.5:
	   		print str(star[st]),str(_RA),str(_DEC),str(distance)
                        refstar=string.split(star[st])[0]
                        ra0,dec0=ra[st],dec[st]
                        if _interactive:
                            print ''
                            print '############################################'
                            question=raw_input('Is this the right object ? [y/n] [y]')
                            print '############################################'
                            print ''
                            if not question: question='y'
                            if question in ['Yes','Y','y','YES','yes']: 
                                break
                            else:
                                ra0,dec0='',''
                        else:
                            print '############################################'
                            print ' object found '+str(refstar)
                            print '############################################'
                            break
       if ra0:
           os.system('rm -rf tmp.*')
           ff=open('tmp.tv','w')
           ff.write(str(ra0*15.)+' '+str(dec0))
           ff.close()
           iraf.wcsctran('tmp.tv','tmp.pix',imgl,inwcs='world',units='degrees degrees',outwcs='logical',columns='1 2',formats='%10.1f %10.1f')
           iraf.tvmark(1,'tmp.pix',mark="circle",number='yes',radii=10,nxoffse=5,nyoffse=5,color=214,txsize=2)
           xx,yy=string.split(iraf.fields('tmp.pix','1,2',Stdout=1)[2])
           print xx,yy
           os.system('rm -rf tmp.*')
       else:
           xx,yy='',''
           print '### WARNING: no object in the list'
       return xx,yy

############################################################################################

def updateheader(image,dimension,_headername,_value):
    from astropy.io import fits as pyfits

    try:
        imm = pyfits.open(image,mode='update')
        _header=imm[dimension].header
        _header[_headername] = _value
        #_header.update(_headername,_value)
        imm.flush()
        imm.close()
    except Exception as e:
        print e
        from snoopy2 import src
        print 'warning: problem to update header, try to correct header format ....'
        src.correctcard(image)
        import sys
        try:
            imm=pyfits.open(image,mode='update')
            _header=imm[dimension].header
            _header.update(_headername,_value)
            imm.flush()
            imm.close()
        except Exception as e:
            print e
            print 'error: not possible update header'
            sys.exit()

######################################
