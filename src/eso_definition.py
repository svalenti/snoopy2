#HIERARCH ESO INS MODE = MOD
#HIERARCH ESO INS FILT3 NAME = FILTER3
#HIERARCH ESO INS FILT4 NAME = FILTER4
#HIERARCH ESO INS SLIT2 NAME = SLIT
#HIERARCH ESO INS GRAT NAME  = GRAT
#HIERARCH ESO INS GRAT2 NAME = GRAT2
#HIERARCH ESO INS GRAT WLEN  = WLEN
#HIERARCH ESO INS GRAT1 WLEN = WLEN1
#HIERARCH ESO INS GRAT2 WLEN = WLEN2
#HIERARCH ESO INS GRAT ORDER = ORDER
#HIERARCH ESO DET OUT1 GAIN   = GAIN
#HIERARCH ESO DET OUT1 RON    = RON

def hierarch28(img,_telescope,header_name,headerupdate):
    import snoopy2
    from snoopy2 import src
    import os,sys,time,string,re
    parent_dir=os.getcwd()+'/'
    if img[0]=='/':
        os.chdir(re.sub(string.split(img,'/')[-1],'',img))
        img=string.split(img,'/')[-1]
        print 'attenzione'
    if img[0:5]=='home$':
        os.chdir(re.sub(string.split(img,'/')[-1],'',img))
        img=string.split(img,'/')[-1]
    try:
        sss=os.environ['IRAFARCH']
    except:
        sss='linux'
    if sss=='macintel':
        directory_prog=snoopy2.__path__[0]+'/exec_prog/macintel'
    elif sss=='redhat' or sss=='linux':
        directory_prog=snoopy2.__path__[0]+'/exec_prog/redhat'
    else:
       print 'ERROR: architecture not reconized !!!'
       sys.exit()
    HIERARCH={}
    HIERARCH['hed_object']= 'HIERARCH ESO DPR TECH = '+str(headerupdate)+' '
    HIERARCH['hed_imagetyp']= 'HIERARCH ESO DPR CATG = '+str(headerupdate)+' '
    HIERARCH['hed_filter1']= 'HIERARCH ESO INS FILT1 NAME = '+str(headerupdate)+' '
    HIERARCH['hed_grism']= 'HIERARCH ESO INS GRIS1 NAME = '+str(headerupdate)+' '
    HIERARCH['hed_slitw']= 'HIERARCH ESO INS SLIT1 NAME = '+str(headerupdate)+' '
    HIERARCH['hed_airmass']= 'HIERARCH ESO TEL AIRM START = '+str(headerupdate)+' '
    HIERARCH['hed_UT']= 'UTC = '+str(headerupdate)+' '
    HIERARCH['pippo']= 'pippovecchio    = '+str(headerupdate)+' '
    if not HIERARCH[header_name]:
         print 'WARNING:  no header in hierarch28 table ' 
         sys.exit()
    else:
        print 'run hierarch28 on image:'+img+'  for header = '+str(headerupdate)
        os.system('rm -rf table.conv')
        f=open('table.conv','w')
        f.write(HIERARCH[header_name])
        f.close()
        os.system(directory_prog+'/hierarch28 '+img+' table.conv > _tmp')
        os.system('cp '+img+' _'+img)
        os.system('rm -rf '+img)
        os.system('cp _'+img+' '+img)
        os.system('rm -rf _'+img)
        os.system('rm -rf _tmp')
        os.chdir(parent_dir)
    return


def esosecondheader(img,_telescope,header_name,headerupdate):
    from pyfits import getheader
    import string,re
    if img[-5:]!='.fits':
        img=img+'.fits'
    a=getheader(img)
    HIERARCH={}
    HIERARCH['hed_object']= 'HIERARCH ESO DPR TECH'
    HIERARCH['hed_imagetyp']= 'HIERARCH ESO DPR CATG'
    HIERARCH['hed_filter1']= 'HIERARCH ESO INS FILT1 NAME'
    HIERARCH['hed_grism']= 'HIERARCH ESO INS GRIS1 NAME'
    HIERARCH['hed_slitw']= 'HIERARCH ESO INS SLIT1 NAME'
    HIERARCH['hed_airmass']= 'HIERARCH ESO TEL AIRM START'
    HIERARCH['hed_UT']= 'UTC'
    if not HIERARCH[header_name]:
         print 'WARNING:  no header in hierarch28 table ' 
         sys.exit()
         return ''
    else:
        return a[HIERARCH[header_name]]
