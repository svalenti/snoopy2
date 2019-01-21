
######################################## leggi file ottico #########################
def Readphfile_o(ascifile,k):   # read ascii file
    import string
    from numpy import math
    _field=0
    f=file(ascifile,'r')
    tit=f.readline()
    s=f.readlines()
    f.close()
    fieldin={}
    for i in range(len(s)):
       if s[i][0:3]=='***':
             _field=_field+1
             fieldin[_field]=i
    field=raw_input('How many standard fields ['+str(_field)+']? ')
    if not field: field=_field
    else: field=int(field)
    mu,mb,mv,mr,mi,st={},{},{},{},{},{}
    V,BV,UB,VR,RI,VI={},{},{},{},{},{}
    name_field=[]
    tu,tb,tv,tr,ti=[],[],[],[],[]
    au,ab,av,ar,ai=[],[],[],[],[]
    numstar=[]
    for j in range(1,field+1):
       kk=fieldin[j]
       xxx,_name_field,num=string.split(s[kk])
       numstar.append(num)
       tu,tb,tv,tr,ti=string.split(s[kk+1])
       au,ab,av,ar,ai=string.split(s[kk+2])
       _u,_b,_v,_r,_i=[],[],[],[],[]
       _V,_BV,_UB,_VR,_RI,_VI,_st=[],[],[],[],[],[],[]
       for ii in range(kk+3,kk+3+int(num)):
          varu,varb,varv,varr,vari,varst,varV,varBV,varUB,varVR,varRI=string.split(s[ii])
          if au!='999':  varu=float(varu)+2.5*math.log10(float(tu))-k['U']*float(au)
          if ab!='999':  varb=float(varb)+2.5*math.log10(float(tb))-k['B']*float(ab)
          if av!='999':  varv=float(varv)+2.5*math.log10(float(tv))-k['V']*float(av)
          if ar!='999':  varr=float(varr)+2.5*math.log10(float(tr))-k['R']*float(ar)
          if ai!='999':  vari=float(vari)+2.5*math.log10(float(ti))-k['I']*float(ai)

          _u.append(-float(varu)+(float(varV)+float(varBV)+float(varUB)))
          _b.append(-float(varb)+(float(varV)+float(varBV)))
          _v.append(-float(varv)+float(varV))
          _r.append(-float(varr)+(float(varV)-float(varVR)))
          _i.append(-float(vari)+(float(varV)-float(varVR)-float(varRI)))

          _st.append(varst)
          _V.append(float(varV))
          _BV.append(float(varBV))
          _UB.append(float(varUB))
          _VR.append(float(varVR))
          _VI.append(float(varVR)+float(varRI))
          _RI.append(float(varRI))
       name_field.append(_name_field)
       mu[j],mb[j],mv[j],mr[j],mi[j],st[j],V[j],BV[j],UB[j],VR[j],RI[j]=_u,_b,_v,_r,_i,_st,_V,_BV,_UB,_VR,_RI
       VI[j]=_VI
    return   mu,mb,mv,mr,mi,st,V,BV,UB,VR,RI,VI,field,name_field

##################################### leggi file infrarosso #############################
def Readphfile_i(ascifile,k):   # read ascii file
    import string
    from numpy import math
    _field=0
    f=file(ascifile,'r')
    tit=f.readline()
    s=f.readlines()
    f.close()
    fieldin={}
    for i in range(len(s)):
       if s[i][0:3]=='***':
             _field=_field+1
             fieldin[_field]=i
    field=raw_input('How many standard fields ['+str(_field)+']? ')
    if not field: field=_field
    else: field=int(field)
    mj,mh,mk,st={},{},{},{}
    J,JH,HK,JK={},{},{},{}
    name_field=[]
    tj,th,tk=[],[],[]
    aj,ah,ak=[],[],[]
    numstar=[]
    for jj in range(1,field+1):
       kk=fieldin[jj]
       xxx,_name_field,num=string.split(s[kk])
       numstar.append(num)
       tj,th,tk=string.split(s[kk+1])
       aj,ah,ak=string.split(s[kk+2])
       _j,_h,_k=[],[],[]
       _J,_JH,_HK,_JK,_st=[],[],[],[],[]
       for ii in range(kk+3,kk+3+int(num)):
          varj,varh,vark,varst,varJ,varJH,varHK=string.split(s[ii])
          if aj!='999':  varj=float(varj)+2.5*math.log10(float(tj))-k['J']*float(aj)
          if ah!='999':  varh=float(varh)+2.5*math.log10(float(th))-k['H']*float(ah)
          if ak!='999':  vark=float(vark)+2.5*math.log10(float(tk))-k['K']*float(ak)

          _j.append(-float(varj)+float(varJ))
          _h.append(-float(varh)+(float(varJ)-float(varJH)))
          _k.append(-float(vark)+(float(varJ)-float(varJH)-float(varHK)))

          _st.append(varst)
          _J.append(float(varJ))
          _JH.append(float(varJH))
          _JK.append(float(varJH)+float(varHK))
          _HK.append(float(varHK))
       name_field.append(_name_field)
       mj[jj],mh[jj],mk[jj],st[jj],J[jj],JH[jj],HK[jj]=_j,_h,_k,_st,_J,_JH,_HK
       JK[jj]=_JK
    return   mj,mh,mk,st,J,JH,HK,JK,field,name_field


################################  leggi file sloan  ######################################

def Readphfile_s(ascifile,k):   # read ascii file
    import string
    from numpy import math
    _field=0
    f=file(ascifile,'r')
    tit=f.readline()
    s=f.readlines()
    f.close()
    fieldin={}
    for i in range(len(s)):
       if s[i][0:3]=='***':
             _field=_field+1
             fieldin[_field]=i
    field=raw_input('How many standard fields ['+str(_field)+']? ')
    if not field: field=_field
    else: field=int(field)
    mu,mg,mr,mi,mz,st={},{},{},{},{},{}
    r,gr,ug,ri,iz,gi={},{},{},{},{},{}
    name_field=[]
    tu,tg,tr,ti,tz=[],[],[],[],[]
    au,ag,ar,ai,az=[],[],[],[],[]
    numstar=[]
    for j in range(1,field+1):
       kk=fieldin[j]
       xxx,_name_field,num=string.split(s[kk])
       numstar.append(num)
       tu,tg,tr,ti,tz=string.split(s[kk+1])
       au,ag,ar,ai,az=string.split(s[kk+2])
       _u,_g,_r,_i,_z=[],[],[],[],[]
       _rr,_ug,_gr,_ri,_iz,_gi,_st=[],[],[],[],[],[],[]
       for ii in range(kk+3,kk+3+int(num)):
          var_u,var_g,var_r,var_i,var_z,varst,varr,vargr,varug,varri,variz=string.split(s[ii])
          if au!='999':  var_u=float(var_u)+2.5*math.log10(float(tu))-k['u']*float(au)
          if ag!='999':  var_g=float(var_g)+2.5*math.log10(float(tg))-k['g']*float(ag)
          if ar!='999':  var_r=float(var_r)+2.5*math.log10(float(tr))-k['r']*float(ar)
          if ai!='999':  var_i=float(var_i)+2.5*math.log10(float(ti))-k['i']*float(ai)
          if az!='999':  var_z=float(var_z)+2.5*math.log10(float(tz))-k['z']*float(az)

          _u.append(-float(var_u)+(float(varr)+float(vargr)+float(varug)))
          _g.append(-float(var_g)+(float(varr)+float(vargr)))
          _r.append(-float(var_r)+float(varr))
          _i.append(-float(var_i)+(float(varr)-float(varri)))
          _z.append(-float(var_z)+(float(varr)-float(varri)-float(variz)))
          _st.append(varst)
          _rr.append(float(varr))
          _ug.append(float(varug))
          _gr.append(float(vargr))
          _ri.append(float(varri))
          _gi.append(float(vargr)+float(varri))
          _iz.append(float(variz))
       name_field.append(_name_field)
       mu[j],mg[j],mr[j],mi[j],mz[j],st[j],r[j],ug[j],gr[j],ri[j],iz[j]=_u,_g,_r,_i,_z,_st,_rr,_ug,_gr,_ri,_iz
       gi[j]=_gi
    return   mu,mg,mr,mi,mz,st,r,gr,ug,ri,iz,gi,field,name_field

###################################################################################


#############################  cancella punti ############################################
def cancella(mm,mma,phmin,phmax,xx,yy,position,xxcut,yycut,positioncut,st,name_field,_color,_band,label,_p0,x_elim,fixcolor,fisso,interactive):
    from numpy import math
    from numpy import array
    from numpy import arange
    from pylab import plot
    from pylab import gca
    from pylab import setp
    from pylab import title
    from pylab import xlabel
    from pylab import ylabel
    from pylab import getp
    from pylab import legend
    from pylab import annotate
    from pylab import ginput
    from pylab import xlim
    from pylab import ylim
    from pylab import clf
    import socket,os,re,string,glob
    from pylab import rcParams
    from matplotlib.font_manager import FontProperties

    if interactive:
        from pylab import ion
        ion()
    else:
        from pylab import show
        from pylab import draw

    clf()
    xlim(phmin,phmax)
    ylim(mm,mma)
    pointx_elim=[]
    pointy_elim=[]
    if x_elim:
        for i in x_elim:
            for j in range(len(xx)):
                if list(position[j])==i:
                    pointx_elim.append(xx[j])
                    pointy_elim.append(yy[j])

 
    def SetPlot():           ############################ set plot defaults
       yticklabels = getp(gca(),'yticklabels')
       xticklabels = getp(gca(),'xticklabels')
       setp(xticklabels,fontsize='20')    
       setp(yticklabels,fontsize='20')    
       font = FontProperties( size="smaller" )
       legend(numpoints=1,markerscale=.8,loc=(1.01,0.0),ncol=1,prop=font ,fancybox=True)
       #legend(numpoints=1,markerscale=1.5)
       leg = gca().get_legend()
       ltext = leg.get_texts()
       setp(ltext,fontsize=15)

    def SY(plcolor):  ###################################    define plot symbol
       sym = 'osDd1234hHpx+<>^vosDd1234hHpx+<>^v' # simboli
       col = 'bgrcmk'            # colori
       _sy = []
       if plcolor: colrange = col
       else:     colrange = 'k'
       for s in sym:
           for c in colrange:      
               _sy.append(c+s)
       return _sy

    def residuals(p,y,xx):
      err = (y-pval(xx, p))
      return err

    def pval(_xx, p):
      _y=+p[0]+p[1]*_xx
      return _y

    def residuals0(p,y,xx):
      err = (y-pval0(xx, p))
      return err

    def pval0(_xx, p):
      _y=+p[0]+0*_xx
      return _y

    def residualsfix(p,y,xx):
       err = (y-pvalfix(xx, p))
       return err

    def pvalfix(_xx, p):
       _y=+p[0]+fisso*_xx
       return _y

    def new_fit(xx,yy,fixvalue): 
        from numpy.linalg import lstsq
        from numpy import average
        from numpy import std
        from numpy import compress
        from numpy import mean
        from numpy import array
        from numpy import ones
        from numpy import sqrt
        if fixvalue or fixvalue==0.0:
            xx*fixvalue
            zero=mean(yy-xx*fixvalue)
            slope=fixvalue
            sloperr=0
            zerr=std(yy-xx*fixvalue)
        else:
#            f=open('test','w')
#            for i in range(0,len(xx)):
#               f.write(str(xx[i])+' '+str(yy[i])+' \n')
#            f.close()
            xx=array(xx)
            yy=array(yy)
            A = ones((len(yy), 2), dtype=float)
            A[:,0] = xx
            result=lstsq(A,yy)
            zero=result[0][1]
            slope=result[0][0]
            if len(xx)>2:
                dd=(1./(len(xx)))+(((mean(xx))**2)/(sum((xx-mean(xx))**2)))
                cc=1/(sum((xx-mean(xx))**2))
                zerr=sqrt((result[1][0]/(len(xx)-2)))*sqrt(dd)
                sloperr=sqrt((result[1][0]/(len(xx)-2)))*sqrt(cc)
            else:
                zerr=0
                sloperr=0
            #Print result
        return zero,slope,zerr,sloperr

##########################################################################

    titolo='%3.3s  %6.6s (%4.4s)  %3.3s  %6.6s (%4.4s)\n' % (_band,str(_p0[0]),str(_p0[2]),_color,str(_p0[1]),str(_p0[3]))
    title(titolo)
    xlabel(_color);
    ylabel(_band);

    n=len(name_field)
    for i in range(0,n):
      xxx,yyy,lab=[],[],[]
      for j in range(len(yy)):
          if position[j][0]==i:
                xxx.append(xx[j])
                yyy.append(yy[j])
                lab.append(label[j])
      xlim(phmin,phmax)
      ylim(mm,mma) 
      pp = plot(xxx,yyy,SY(plcolor=False)[position[i][0]],label=' ') #name_field[i])
      xlim(phmin,phmax)
      ylim(mm,mma) 
      setp(pp,markerfacecolor='none')
      xlim(phmin,phmax)
      ylim(mm,mma) 
      setp(pp,markersize=10)
#      for ii in range(len(xxx)):
#          m=i+1
#          annotate(lab[ii], xy=(xxx[ii], yyy[ii]),  xycoords='data',
#                   xytext=(-15, 10), textcoords='offset points',
#                   horizontalalignment='right', verticalalignment='bottom')

    if pointx_elim:
        xlim(phmin,phmax)
        ylim(mm,mma) 
        pp = plot(pointx_elim,pointy_elim,'pb')
        xlim(phmin,phmax)
        ylim(mm,mma) 
        setp(pp,markersize=10)
    
    xlim(phmin,phmax)
    ylim(mm,mma) 
    SetPlot()
    if not interactive:
        draw()
        draw()

    print 'MARK POINT '
    fff=ginput(1,timeout=100)
    xlim(phmin,phmax)
    ylim(mm,mma)


    xw=fff[0][0]
    yw=fff[0][1]
    xlim(phmin,phmax)
    ylim(mm,mma)
    pp = plot([xw],[yw])
    xlim(phmin,phmax)
    ylim(mm,mma)
    setp(pp,markersize=5)
    
    dd=[]
    for j in range(len(xxcut)):
            dd.append(math.sqrt((xw-xxcut[j])**2+(yw-yycut[j])**2))
         
    for i in range(len(dd)):
        if dd[i]==min(dd):
            x_elim.append(list(positioncut[i]))
            pointx_elim.append(xxcut[i])
            pointy_elim.append(yycut[i])

    for i in x_elim:
        for j in range(len(xxcut)):
                if positioncut[j]==i:
                    xxcut=array(list(xxcut[0:j])+list(xxcut[j+1:]))
                    yycut=list(yycut[0:j])+list(yycut[j+1:])
                    positioncut=list(positioncut[0:j])+list(positioncut[j+1:])
                    print 'leng '+str(len(xxcut))+' '+str(len(positioncut))
                    break

    print pointx_elim,pointy_elim
    if pointx_elim:
        xlim(phmin,phmax)
        ylim(mm,mma)
        pp = plot(pointx_elim,pointy_elim,'pr')
        setp(pp,markersize=5)


    if fixcolor:
       if len(xxcut)==1:
           p=[yycut[0],0,0,0]
           xfit=arange(phmin,phmax+0.099,.1)
           yfit=pval0(xfit, p)
       else:
           _p0[1]=fisso
           xfit=arange(phmin,phmax+0.099,.1)
           zero,slope,zeroerr,sloperr=new_fit(xxcut,yycut,fisso)
           p=[zero,slope,zeroerr,sloperr]
#           plsq = leastsq(residualsfix, _p0, args=(yycut,xxcut))
#           xfit=arange(phmin,phmax+0.099,.1)
#           p=plsq[0]
           yfit=pvalfix(xfit, p)
    else:
       if len(xxcut)==1:
           p=[yycut[0],0,0,0]
           xfit=arange(phmin,phmax+0.099,.1)
           yfit=pval0(xfit, p)
       else:
#            plsq = leastsq(residuals, _p0, args=(yycut,xxcut))
            zero,slope,zeroerr,sloperr=new_fit(xxcut,yycut,'')
            p=[zero,slope,zeroerr,sloperr]
            xfit=arange(phmin,phmax+0.099,.1)
#            p=plsq[0]
            yfit=pval(xfit, p)
    
    print 'zero point and color term = '+str(p)
    try:
        titolo='%3.3s  %6.6s (%4.4s)  %3.3s  %6.6s (%4.4s)\n' % (_band,str(p[0]),str(p[2]),_color,str(p[1]),str(p[3]))
        #titolo='%5.5s   %6.6s   %6.6s   %6.6s\n' % (_band,str(p[0]),_color,str(p[1]))
    except:
        titolo='%3.3s  %6.6s (%4.4s)  %3.3s  %6.6s (%4.4s)\n' % (_band,str(_p0[0]),str(_p0[2]),_color,str(_p0[1]),str(_p0[3]))
#        titolo='%5.5s   %6.6s   %6.6s   %6.6s\n' % (_band,str(_p0[0]),_color,str(_p0[1]))
    xlim(phmin,phmax)
    ylim(mm,mma)
    print titolo
    title(titolo)
    pp = plot(xfit,yfit,'-b',label='fit')
    xlim(phmin,phmax)
    ylim(mm,mma)
    setp(pp,markersize=5)
    xlim(phmin,phmax)
    ylim(mm,mma)
    SetPlot()
    if not interactive:
       draw()
       draw()

    return xx,yy,position,xxcut,yycut,positioncut,x_elim,p
 
#################################################################################################
############################## pylab plot  ############################################
def grafico_phnew(_band,_color,_p,phmax,phmin,mma,mm,xx,yy,st,position,label,name_field,xfit,yfit,fixcolor,fisso,interactive):

   from pylab import plot
   from pylab import gca
   from pylab import setp
   from pylab import title
   from pylab import xlabel
   from pylab import ylabel
   from pylab import getp
   from pylab import legend
   from pylab import annotate
   from pylab import close
   from pylab import figure
   from pylab import ginput
   from pylab import xlim
   from pylab import ylim
   from pylab import clf
   import socket,os,re,string,glob

   if interactive:
       from pylab import ion
       ion()
       print 'using interactive mode for ubuntu !!!!!'
   else:
       from pylab import show
       from pylab import draw

   from pylab import rcParams
   from pylab import subplot
   from pylab import subplots_adjust
   from matplotlib.font_manager import FontProperties

   rcParams['figure.figsize'] = 11, 5
   figure()
   subplot(111)
   subplots_adjust(right=0.8)
   def SetPlot():           ############################ set plot defaults
       yticklabels = getp(gca(),'yticklabels')
       xticklabels = getp(gca(),'xticklabels')
       setp(xticklabels,fontsize='20')
       setp(yticklabels,fontsize='20')    
       font = FontProperties( size="smaller" )
       legend(numpoints=1,markerscale=.8,loc=(1.01,0.0),ncol=1,prop=font ,fancybox=True)
       #legend(numpoints=1,markerscale=1.5)
       leg = gca().get_legend()
       ltext = leg.get_texts()
       setp(ltext,fontsize=15)

   def SY(plcolor):  ###################################    define plot symbol
       sym = 'osDd1234hHpx+<>^vosDd1234hHpx+<>^v' # simboli
       col = 'bgrcmk'            # colori
       _sy = []
       if plcolor: colrange = col
       else:     colrange = 'k'
       for s in sym:
           for c in colrange:      
               _sy.append(c+s)
       return _sy

   titolo='%3.3s  %6.6s (%4.4s)  %3.3s  %6.6s (%4.4s)\n' % (_band,str(_p[0]),str(_p[2]),_color,str(_p[1]),str(_p[3]))
   title(titolo)
   xlabel(_color)
   ylabel(_band)

   n=len(name_field)
   for i in range(0,n):
      xxx,yyy,lab=[],[],[]
      for j in range(len(yy)):
          if position[j][0]==i:
                xxx.append(xx[j])
                yyy.append(yy[j])
                lab.append(label[j])
      xlim(phmin,phmax)
      ylim(mm,mma) 
      pp = plot(xxx,yyy,SY(plcolor=False)[position[i][0]],label=name_field[i])
      xlim(phmin,phmax)
      ylim(mm,mma)
      setp(pp,markerfacecolor='none')
      xlim(phmin,phmax)
      ylim(mm,mma)
      setp(pp,markersize=10)
      xlim(phmin,phmax)
      ylim(mm,mma)
      if len(xxx)<=30:
          for ii in range(len(xxx)):
              xlim(phmin,phmax)
              ylim(mm,mma)
              annotate(lab[ii], xy=(xxx[ii], yyy[ii]),  xycoords='data',
                       xytext=(-4, 4), textcoords='offset points',fontsize='10',
                       horizontalalignment='right', verticalalignment='bottom')
   xlim(phmin,phmax)
   ylim(mm,mma) 
   ppp=plot(xfit,yfit)
   xlim(phmin,phmax)
   ylim(mm,mma)
   setp(ppp,markerfacecolor='none')
   xlim(phmin,phmax)
   ylim(mm,mma)
   SetPlot()
   xlim(phmin,phmax)
   ylim(mm,mma)

   if not interactive:     
       draw()
       draw()

   fff=ginput(1,timeout=0.1,show_clicks=False)

   print ''
   print '################################################'
   answ=raw_input('Do you want to cancel a point ? [y/n] [n]')
   if not answ: answ='n'
   xxcut,yycut,positioncut=xx,yy,position
   x_elim=[]
   while answ=='y':
      xx,yy,position,xxcut,yycut,positioncut,x_elim,p=cancella(mm,mma,phmin,phmax,xx,yy,position,xxcut,yycut,positioncut,st,name_field,_color,_band,label,_p,x_elim,fixcolor,fisso,interactive)


      ######################################################
      ###  added to avoid crash on elisa fedora 11 
      fff=ginput(1,timeout=0.1,show_clicks=False)
      ######################################################
      if len(xxcut)>=2:
         answ=raw_input('Do you want to cancel a point ? [y/n]  [y]')
      else:
         answ='n'
      if not answ: answ='y'
   try:
       qubacont='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(_band,p[0],_color,p[1],p[2],p[3])
   except:
       qubacont='%3.3s  %6.6s  %3.3s  %6.6s | %5.5s  %5.5s'%(_band,_p[0],_color,_p[1],_p[2],_p[3])

   close()
   close()

   return qubacont

######################################################################################

def atmospheric_site2():
   import sys,string
   k={}
   kk={}
   site=raw_input('site: e-so lasilla, p-aranal, a-siago, c-alar alto, s-iding spring, m-auna kea, r-oque, k-ait, t-ctio, ms-montsec, h-ct (India 2.0), st-(India 1.0) or b-ao ? [] ')
   bands='UuBgVrRiIzJHK'
   kk['telescopse']=['U','u','B','g','V','r','R','i','I','z','J','H','K']
   kk['e'] = [0.46, 0.46, 0.27, 0.20, 0.12, 0.09, 0.09, 0.02, 0.02, 0.03, 0.07, 0.023, 0.045]  
   kk['p'] = [0.43, 0.43, 0.22, 0.18, 0.11, 0.07, 0.07, 0.05, 0.05, 0.06, 0.11, 0.06, 0.07]
   kk['r'] = [0.46, 0.46, 0.22, 0.16, 0.12, 0.08, 0.08, 0.04, 0.04, 0.06, 0.12, 0.06, 0.09]
   kk['c'] = [0.47, 0.47, 0.35, 0.30, 0.23, 0.17, 0.17, 0.09, 0.09, 0.10, 0.102, 0.07, 0.09]  
   kk['s'] = [0.63, 0.70, 0.32, 0.26, 0.18, 0.150, 0.13, 0.08, 0.07, 0.06, 0.0, 0.0, 0.0]
   kk['m'] = [0.45, 0.48, 0.21, 0.16, 0.12, 0.09, 0.07,0.04, 0.03, 0.03, 0.0, 0.0, 0.0]
   kk['a'] = [0.58, 0.58, 0.29, 0.20, 0.16, 0.12, 0.12, 0.08, 0.08, 0.09, 0.16, 0.09, 0.19]
   kk['k'] = [0.56, 0.56, 0.28, 0.22, 0.17, 0.13, 0.13, 0.07, 0.07, 0.09, 0.0, 0.0, 0.0]
   kk['b'] = [0.76, 0.76, 0.41, 0.35, 0.30, 0.28, 0.28, 0.13, 0.13, 0.15, 0.0, 0.0, 0.0]
   kk['w'] = [0.70, 0.70, 0.60, 0.55, 0.52, 0.42, 0.42, 0.27, 0.27, 0.30, 0.0, 0.0, 0.0]
   kk['h'] = [0.36, 0.0 ,0.21, 0.0, 0.12, 0.0, 0.09, 0.0, 0.05, 0.0, 0.0, 0.0,  0.0]
   kk['st'] = [0.57, 0.0 ,0.28, 0.0, 0.17, 0.0, 0.11, 0.0, 0.07, 0.0, 0.0, 0.0,  0.0]
   kk['t'] = [0.50, 0.516, 0.232, 0.203, 0.144, 0.116, 0.083, 0.077, 0.056, 0.04, 0.0, 0.0, 0.0]
   kk['ms'] = [0.50, 0.516, 0.23, 0.23, 0.15, 0.11, 0.09, 0.04, 0.04, 0.04, 0.0, 0.0, 0.0]
   kk['0'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   for band in bands:
       try:
           k[band]=kk[site][string.find(bands,band)]
       except:
             print '#################################'
             print 'WARNING:  band not found !!!!'
             print '#################################'
             k[band]=0.0
             #sys.exit()
   return k   
##########################################################################

def grafico_exa(xx,yy,position,zz,_label,p,banda,xfit,yfit,phmin,phmax,mm,mma,interactive,pnew):

   from pylab import plot
   from pylab import gca
   from pylab import setp
   from pylab import title
   from pylab import xlabel
   from pylab import ylabel
   from pylab import getp
   from pylab import legend
   from pylab import annotate
   from pylab import close
   from pylab import figure
   from pylab import ginput
   from pylab import xlim
   from pylab import ylim
   from pylab import clf
   import socket,os,re,string,glob

   if interactive:
       from pylab import ion
       ion()
       print 'using interactive mode for ubuntu !!!!!'
   else:
       from pylab import show
       from pylab import draw

   from pylab import rcParams
   from pylab import subplot
   from pylab import subplots_adjust
   from matplotlib.font_manager import FontProperties

   rcParams['figure.figsize'] = 11, 5
   figure()
   subplot(111)
   subplots_adjust(right=0.8)


   def SetPlot():           ############################ set plot defaults
       yticklabels = getp(gca(),'yticklabels')
       xticklabels = getp(gca(),'xticklabels')
       setp(xticklabels,fontsize='20')
       setp(yticklabels,fontsize='20')
       font = FontProperties( size="smaller" )
       #font= FontProperties(size='x-small');
       legend(numpoints=1,markerscale=.8,loc=(1.01,0.0),ncol=1,prop=font ,fancybox=True)
       
       leg = gca().get_legend()
       ltext = leg.get_texts()
       setp(ltext,fontsize=10)

   def SY(plcolor):  ###################################    define plot symbol
       sym = 'osDd1234hHpx+<>^vosDd1234hHpx+<>^v' # simboli
       col = 'bgrcmk'            # colori
       _sy = []
       if plcolor: colrange = col
       else:     colrange = 'k'
       for s in sym:
           for c in colrange:      
               _sy.append(c+s)
       return _sy

   media,error=p[0],p[1]
   titolo='%5.5s  +-  %6.6s\t\n' % (str(media),str(error))
   title(titolo)
   xlabel('NIGHTS')
   ylabel(banda)

   xlim(phmin-1,phmax+1)
   ylim(mm,mma)

   for ii in range(len(xx)):
            pp = plot([xx[ii]],[yy[ii]],SY(plcolor=False)[zz[ii]],label=_label[ii])
            if zz[ii]!=1:
                setp(pp,markerfacecolor='none')
            setp(pp,markersize=10)
   ppp=plot(xfit,yfit,'-r',label='fit')
   setp(ppp,markerfacecolor='none')
   xlim(phmin-1,phmax+1)
   ylim(mm,mma)
   SetPlot()
   xlim(phmin-1,phmax+1)
   ylim(mm,mma)

   if not interactive: 
      draw()
      draw()

   fff=ginput(1,timeout=0.1,show_clicks=False)

   print ''
   print '################################################'
   answ=raw_input('Do you want to cancel a point ? [y/n] [n]')
   if not answ: answ='n'
   xxcut,yycut,positioncut=[],[],[]
   print xx,yy,position
   xxcut,yycut,positioncut=xx,yy,position
   x_elim=[]
   while answ=='y':
      print 'cancella'
      xx,yy,position,xxcut,yycut,positioncut,x_elim,pnew=cancella_exa(mm,mma,phmin,phmax,xx,yy,zz,position,xxcut,yycut,positioncut,banda,_label,p,pnew,x_elim,interactive)

      ######################################################
      ###  added to avoid crash on elisa fedora 11 
      fff=ginput(1,timeout=0.1,show_clicks=False)
      ######################################################
      if len(xxcut)>=2:
         answ=raw_input('Do you want to cancel a point ? [y/n]  [y]')
      else:
         answ='n'
      if not answ: answ='y'

   print pnew
   close()
   close()

   return pnew

################################################################################################
def cancella_exa(mm,mma,phmin,phmax,xx,yy,zz,position,xxcut,yycut,positioncut,banda,_label,p,pnew,x_elim,interactive):
    from numpy import math
    from numpy import array
    from numpy import arange
    from numpy import mean
    from numpy import std
    from pylab import plot
    from pylab import gca
    from pylab import setp
    from pylab import title
    from pylab import xlabel
    from pylab import ylabel
    from pylab import figure
    from pylab import getp
    from pylab import legend
    from pylab import annotate
    from pylab import ginput
    from pylab import xlim
    from pylab import ylim
    from pylab import clf
    from pylab import rcParams
    from matplotlib.font_manager import FontProperties
    import socket,os,re,string,glob

    if interactive:
        from pylab import ion
        ion()
    else:
        from pylab import show
        from pylab import draw

    media,error=pnew[0],pnew[1]
    clf()
    xlim(phmin-1,phmax+(phmax-phmin))
    ylim(mm,mma)
    pointx_elim=[]
    pointy_elim=[]
    if x_elim:
        for i in x_elim:
            for j in range(len(xx)):
                if list(position[j])==i:
                    pointx_elim.append(xx[j])
                    pointy_elim.append(yy[j])

 
    def SetPlot():           ############################ set plot defaults
       yticklabels = getp(gca(),'yticklabels')
       xticklabels = getp(gca(),'xticklabels')
       setp(xticklabels,fontsize='20')    
       setp(yticklabels,fontsize='20')    
       font = FontProperties( size="smaller" )
       legend(numpoints=1,markerscale=.8,loc=(1.01,0.0),ncol=1,prop=font ,fancybox=True)
       leg = gca().get_legend()
       ltext = leg.get_texts()
       setp(ltext,fontsize=10)


    def SY(plcolor):  ###################################    define plot symbol
       sym = 'osDd1234hHpx+<>^vosDd1234hHpx+<>^v' # simboli
       col = 'bgrcmk'            # colori
       _sy = []
       if plcolor: colrange = col
       else:     colrange = 'k'
       for s in sym:
           for c in colrange:      
               _sy.append(c+s)
       return _sy

    def pval0(_xx, p):
      _y=+p[0]+0*_xx
      return _y
##########################################################################
    titolo='%5.5s  +-  %6.6s\t\n' % (str(media),str(error))
    title(titolo)
    xlabel('NIGHTS')
    ylabel(banda)

    for ii in range(len(xx)):
            pp = plot([xx[ii]],[yy[ii]],SY(plcolor=False)[zz[ii]],label=_label[ii])
            if zz[ii]!=1:
                setp(pp,markerfacecolor='none')
            setp(pp,markersize=10)
    
    xfit=arange(phmin,phmax,.1)
    yfit=pval0(xfit, p)
    ppp=plot(xfit,yfit,'-r',label='fit')
    setp(ppp,markerfacecolor='none')
    xlim(phmin-1,phmax+1)
    ylim(mm,mma)

    if pointx_elim:
        xlim(phmin-1,phmax+1)
        ylim(mm,mma) 
        pp = plot(pointx_elim,pointy_elim,'pb')
        xlim(phmin-1,phmax+1)
        ylim(mm,mma) 
        setp(pp,markersize=10)
    
    xlim(phmin-1,phmax+1)
    ylim(mm,mma) 
    SetPlot()
    if not interactive:
        draw()
        draw()

    print 'MARK POINT '
    fff=ginput(1,timeout=100)
    xlim(phmin-1,phmax+1)
    ylim(mm,mma)

    xw=fff[0][0]
    yw=fff[0][1]
    xlim(phmin-1,phmax+1)
    ylim(mm,mma)
    pp = plot([xw],[yw])
    xlim(phmin-1,phmax+1)
    ylim(mm,mma)
    setp(pp,markersize=5)
    SetPlot()
    
    dd=[]
    for j in range(len(xxcut)):
            dd.append(math.sqrt((xw-xxcut[j])**2+(yw-yycut[j])**2))
         
    for i in range(len(dd)):
        if dd[i]==min(dd):
            x_elim.append(list(positioncut[i]))
            pointx_elim.append(xxcut[i])
            pointy_elim.append(yycut[i])

    for i in x_elim:
        for j in range(len(xxcut)):
                if positioncut[j]==i:
                    xxcut=array(list(xxcut[0:j])+list(xxcut[j+1:]))
                    yycut=list(yycut[0:j])+list(yycut[j+1:])
                    positioncut=list(positioncut[0:j])+list(positioncut[j+1:])
                    print 'leng '+str(len(xxcut))+' '+str(len(positioncut))
                    break

    print pointx_elim,pointy_elim
    if pointx_elim:
        xlim(phmin-1,phmax+1)
        ylim(mm,mma)
        pp = plot(pointx_elim,pointy_elim,'pr')
        setp(pp,markersize=5)
        
        if len(yycut)>=2:
               media=mean(array(yycut))
               error=std(array(yycut))
        else:
               media=yycut[0]
               error=0
        pnew=[media,error]
        xfit=arange(phmin,phmax,.1)
        yfit=pval0(xfit, pnew)

    print 'zero point and color term = '+str(pnew)
    try:
           titolo='%5.5s  +-  %6.6s\t\n' % (str(pnew[0]),str(pnew[1]))
    except:
           titolo='%5.5s  +-  %6.6s\t\n' % (str(p[0]),str(p[1]))
           pnew=p
    print titolo
    title(titolo)
    pp = plot(xfit,yfit,'-b',label='new fit ')
    setp(pp,markersize=5)
    xlim(phmin-1,phmax+1)
    ylim(mm,mma)
    SetPlot()

    if not interactive:
       draw()
       draw()

    return xx,yy,position,xxcut,yycut,positioncut,x_elim,pnew

#######################################################################################################################
