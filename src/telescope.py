
######################################################################
def telescope(img):
    import string
    from astropy.io import fits as pyfits

    
    telescopelist=['ekar','pennar','TNG','NOT','ACAM','WHT','lo','wise','ca','FORS1','FORS2','NTT','lp','wise','rem','CI','sub','mag','FTN','danish','SOFI','NICS','ctio','montsec','hct','trp']
    geth = pyfits.getheader(img)
    try:
        _telescope = geth['TELESCOP']
        if string.find(_telescope,'Ekar')!=-1:
            _telescope='ekar'
        elif string.find(_telescope,'Reflector')!=-1:
            _telescope='ekar'
        elif string.find(_telescope,'Schmidt')!=-1:
            _telescope='sch'
        elif string.find(_telescope,'ASIAGO')!=-1:
            _telescope='pennar'
        elif string.find(_telescope,'TNG')!=-1:
            if string.find(geth['INSTRUME'],'NICS')!=-1:
                _telescope='NICS'
            else:
                _telescope='TNG'
        elif string.find(_telescope,'NOT')!=-1:
            _telescope='NOT'
        elif string.find(_telescope,'Danish')!=-1:
            _telescope='danish'
        elif string.find(_telescope,'WHT')!=-1:
            _instrum=geth['INSTRUME']
            if _instrum=='ACAM':
                _telescope='ACAM'
            else:
                _telescope='WHT'
        elif string.find(_telescope,'Orzale')!=-1:
            _telescope='lo'
        elif string.find(_telescope,'40 INCH')!=-1:
            _telescope='wise'
	elif string.find(_telescope,'CA-2.2')!=-1:
            _telescope='ca'    
        elif _telescope=='ca':
            _telescope='ca'    
        elif string.find(_telescope,'Subaru')!=-1:
            _telescope='sub'
        elif string.find(_telescope,'K.A.I.T.')!=-1:
            _telescope='kait'
        elif string.find(_telescope,'Clay_Mag')!=-1:
            _telescope='mag'
        elif string.find(_telescope,'ESO-VLT-U1')!=-1:
            _telescope='FORS2'
        elif string.find(_telescope,'ESO-VLT-U2')!=-1:
            _telescope='FORS1'
        elif string.find(_telescope,'ESO-NTT')!=-1:
            if string.find(geth['INSTRUME'],'SOFI')!=-1:
                _telescope='SOFI'
            else:
                _telescope='NTT'
        elif string.find(_telescope,'Wendelstein')!=-1:
            _telescope='wend'
        elif string.find(_telescope,'Liverpool')!=-1:
            _telescope='lp'
        elif string.find(_telescope,'REM')!=-1:
            _telescope='rem'
        elif string.find(_telescope,'AZT-24')!=-1:
            _telescope='CI'	
	elif string.find(_telescope,'Prompt')!=-1:
            _telescope='prompt'
	elif string.find(_telescope,'Faulkes')!=-1:
            if string.find(_telescope,'North')!=-1:
                _instrume=geth['instrume']
                if _instrume=='EM01':
                    _telescope='mer'
                else:
                    _telescope='FTN'
            elif string.find(_telescope,'South')!=-1:
                _telescope='FTS'
        elif string.find(_telescope,'PS1')!=-1:
            _telescope='PS1'	        
        elif string.find(_telescope,'SWIFT')!=-1:
            _telescope='swift'	        
        elif string.find(_telescope,'Dutch')!=-1:
            _telescope='dutch'
	elif string.find(_telescope,'ct13m')!=-1:
            _telescope='ctio'	        
	elif string.find(_telescope,'Montsec')!=-1:
            _telescope='montsec'	        
	elif string.find(_telescope,'TRAPPIST')!=-1:
            _telescope='trp'	        
	elif string.find(_telescope,'Sampur')!=-1:
            _telescope='sampur'
	elif string.find(_telescope,'other')!=-1:
            _telescope='other'
	elif _telescope in ['1m0-04','1m0-05','1m0-08','1m0-09']:
            _telescope='lsc'
        elif _telescope in telescopelist:
            pass
        else:
	   try:
               _telescope=geth['OBSERVAT']
               if string.find(_telescope,'Prompt')!=-1: 
                   _telescope='prompt'
               elif string.find(_telescope,'Meckering')!=-1:
                   _telescope='Meckering'
               else:
                   _telescope=''
	   except: _telescope=''
    except: _telescope=''

    if _telescope=='':
        try:
            _instrume=geth['INSTRUME']
            if string.find(_instrume,'HFOSC')!=-1: 
                   _telescope='hct'
            else: _telescope=''
        except:  _telescope=''
    if _telescope=='':
        try:
            _telescope=geth['QUBATELE']
            if string.find(_telescope,'other')!=-1:
                _telescope='other'
            else: _telescope=''
        except: 
            if pyfits.open(img)[0].header.get('PSCAMERA')=='GPC1':
                _telescope='PS1'
            else:
                _telescope=''
    if _telescope=='':
            print 'WARNING: Telescope not defined !!!'
    return _telescope


def check_tel(img,_aa,_system):
    from snoopy2 import src
    try:
        _header = src.read_parameter(_aa,_system)
        filter = src.filter(img,_header,_aa)
        instrument = src.instrument(img,_header,_aa)
        _airmass = src.airmass(img,_header,_aa)
        _exptime = src.exptime(img,_header,_aa)
        _date = src.date(img,_header,_aa)
        xdim = src.xdimen(img,_aa)
        ydim = src.ydimen(img,_aa)
        _object = src.objects(img,_header,_aa)
        _julday = src.JD(img,_header,_aa)
        _UT = src.UT(img,_header,_aa)
        print '#### Header check ......ok '
        check=1 
    except Exception as e:
        print e
        check=0
    return check

#os.system('rm -rf login.cl')
#if logincl:
#   f=open('login.cl','w')
#   for i in logincl:
#      f.write(i)
#   f.close()

###########################################################################

