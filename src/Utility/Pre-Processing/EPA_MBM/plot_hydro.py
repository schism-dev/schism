#!/usr/bin/env python3
from pylib import *
close("all")

#------------------------------------------------------------------------------
#inputs
#------------------------------------------------------------------------------
#runs=['RUN05b','RUN05g']; StartTs=[datenum(2004,1,1),datenum(2004,1,1)]
runs=['RUN09c']; StartTs=[datenum(1991,1,1),]

#plot option for elev, temp, salinity,sediment
iplot=[1,1,1,1]  #1: plot; 0:skip
ishow=[0,0,0,0]  #1: show figures

#databases
bdir='/sciclone/data10/wangzg/CBP/'
elev_obs=bdir+'../Database/elev/noaa_elev_msl.npz'
elev_obs_stainfo=bdir+'../Database/elev/stainfo.npz'
wqdata=bdir+'setup_files/CBP_WQData.npz'
grd='../grid.npz'

#------------------------------------------------------------------------------
#functions: statistics
#------------------------------------------------------------------------------
STS=zdata() #used to store statistics
def stat(mti,myi,oti,oyi,dt=0):
    #get data pairs
    fp=(oti>=mti[0])*(oti<=mti[-1]); oti,oyi=oti[fp],oyi[fp]; fyi=[]
    if dt!=0:
        for otii,oyii in zip(oti,oyi):
            y=myi[abs(mti-otii)<=dt]; ay=abs(y-oyii); k=nonzero(ay==ay.min())[0][0]; fyi.append(y[k])
        fyi=array(fyi)
    else:
        fyi=interpolate.interp1d(mti,myi)(oti)

    #compute statistics
    if std(oyi)==0:
        s=zdata(); s.R=nan
        dx=fyi-oyi; s.ME=mean(dx); s.MAE=mean(abs(dx)); s.RMSD=sqrt((dx**2).mean())
    else:
        s=get_stat(fyi,oyi)
    return s

#------------------------------------------------------------------------------
#elevation: subtidal signal, tidal
#------------------------------------------------------------------------------
if fexist('elev.npz') and iplot[0]==1:
    if not fexist('elev'): os.mkdir('elev')
    stations=[8573364,8575512,8571892,8635750,8636580,8637689,8632200,8638610] #stations

    #read results
    models=[loadz(bdir+'{}/results/elev.npz'.format(i)) for i in runs]; mts=[i.time+t0 for i,t0 in zip(models,StartTs)]
    C=loadz(elev_obs); S=loadz(elev_obs_stainfo); sdict=dict(zip(S.station,S.staname))

    xts,xls=get_xtick(xts=[1990,2024])
    xm=[min([i.min() for i in mts])-1,max([i.max() for i in mts])] #range of model runs
    mt0=max([i.min()+10 for i in mts]); mt1=mt0+15   #range for plotting tide
    htm=[max([i.min() for i in mts]),min([i.max() for i in mts])] #range for HA

    cs=['g','b','k']; snames=[sdict[i] for i in stations]
    #plot tidal and subtidal
    H=zdata(); H.MA,H.MP,H.OA,H.OP=[],[],[],[];  H.sname=snames; H.station=stations
    E=zdata(); E.R,E.ME,E.MAE,E.RMSD=[],[],[],[]; E.sname=snames; E.station=stations
    for mm in arange(2):
        figure(figsize=[19.5,9.4])
        for m,[station,sname] in enumerate(zip(stations,snames)):
            subplot(4,2,m+1); lstr=[]

            MAi,MPi,OAi,OPi=[],[],[],[]
            #obs
            fp=(C.station==station)*(C.time>=xm[0])*(C.time<=xm[1]); oti,oyi=C.time[fp],C.elev[fp]

            if len(oti)!=0:
               sind=argsort(oti); oti,oyi=oti[sind],oyi[sind]; foyi=lpfilt(oyi, 1/24, 0.1); toyi=oyi-foyi
               if mm==0: plot(oti,foyi,'r-')
               if mm==1: fpt=(oti>=mt0)*(oti<=mt1);  plot(oti[fpt]-mt0,toyi[fpt],'r-');
               lstr.append('Obs: {}'.format(sname))
               hti=arange(max(htm[0],oti.min()),min(htm[1],oti.max()),1/24); hti=hti[:24*365]
               if mm==0: #HA
                   hoyi=interpolate.interp1d(oti,oyi)(hti)
                   T=harmonic_analysis(hoyi,1/24,t0=hti[0]);  H.OA.append(T.amplitude); H.OP.append(T.phase)
            else:
               if mm==0: H.OA.append(-99*ones(9)); H.OP.append(-99*ones(9))

            #model
            for n,run in enumerate(runs):
                M=models[n]; sid=nonzero(M.bp.station==str(station))[0][0]
                mti=M.time+StartTs[n]; myi=M.elev[sid]; fmyi=lpfilt(myi, 1/24, 0.1); tmyi=myi-fmyi
                if mm==0: plot(mti,fmyi,color=cs[n])
                if mm==1: fpt=(mti>=mt0)*(mti<=mt1); plot(mti[fpt]-mt0,tmyi[fpt],color=cs[n]);
                lstr.append(run)

                if mm==0: #HA analysis
                    hmyi=interpolate.interp1d(mti,myi)(hti)
                    T=harmonic_analysis(hmyi,1/24,t0=hti[0]);  MAi.append(T.amplitude); MPi.append(T.phase)
                    H.tidal_name=T.tidal_name; H.freq=T.freq
                    st=get_stat(hmyi,hoyi); E.R.append(st.R); E.ME.append(st.ME); E.MAE.append(st.MAE); E.RMSD.append(st.RMSD)
            if mm==0: H.MA.append(MAi); H.MP.append(MPi)

            #note
            hl=legend(lstr,loc=1); hl.set_alpha(1.0)
            if mm==0:
                setp(gca(),xticks=xts,xticklabels=xls,xlim=xm,ylim=[-0.5,0.7])
                if m%2==0: ylabel('subtidal elev. (m)')
                gca().xaxis.grid('on')
                rtext(0.25,0.94,'R={:6.3f}, ME={:7.4f}, MAE={:7.4f}, RMSD={:7.4f}'.format(st.R,st.ME,st.MAE,st.RMSD))
            elif mm==1:
                txts=arange(15); txls=[num2date(mt0).strftime('%m/%d,%Y') if i==0 else str(i) for i in txts]
                setp(gca(),xticks=txts,xticklabels=txls,xlim=[0,15])

        gcf().tight_layout(); mvfig()
        if mm==0: savefig('elev/elev_subtidal_1')
        if mm==1: savefig('elev/elev_tidal_1')
        if ishow[0]==0: close()

    #write statistics
    H.OA=array(H.OA); H.OP=array(H.OP); H.MA=array(H.MA); H.MP=array(H.MP)
    sdata=array([['subtide',*E.sname],  ['R',*E.R], ['RMSD',*E.RMSD], ['ME',*E.ME], ['MAE',*E.MAE]],dtype='O')
    write_excel('statistics',sdata,sht='elev')
    for m in arange(8):
        ma=H.MA[:,0,m+1]; mp=H.MP[:,0,m+1]*180/pi; oa=H.OA[:,m+1]; op=H.OP[:,m+1]*180/pi; da=ma-oa; dp=mod((mp-op),360); dp[dp>180]=dp[dp>180]-360
        sdata=array([['tide: '+H.tidal_name[m+1],'', *snames], ['Amplitude','model',*ma], ['', 'Obs',*oa], ['','diff',*da], ['Phase','model', *mp], ['','Obs',*op],  ['','diff',*dp] ],dtype='O')
        write_excel('statistics',sdata,sht='elev',indy=6+8*m)
    STS.subtide=E; STS.tide=H

    #plot Harmonic Analysis
    figure(figsize=[19,9.4])
    xi=arange(len(stations)); xls=['Yorktwon' if i.startswith('Yorktown') else i for i in snames]
    for m in arange(2):
        noff=0 if m==0 else 4;
        dw=0.25; dws=[-dw/2,dw/2] if len(runs)==1 else [-dw,0,dw]

        for n in arange(4):
            #amplitude
            subplot(4,4,n+1+2*noff)
            bar(xi+dws[0],H.OA[:,n+1+noff],  width=dw,facecolor='r',label='Obs')
            bar(xi+dws[1],H.MA[:,0,n+1+noff],width=dw,facecolor='g', label=runs[0])
            if len(runs)==2: bar(xi+dws[2],H.MA[:,1,n+1+noff],width=dw,facecolor='b', label=runs[1])

            setp(gca(),xticks=xi,xticklabels=[],xlim=[-0.5,7.5])
            setp(gca(),ylim=[0,0.5]) if H.tidal_name[n+1+noff]=='M2' else setp(gca(),ylim=[0,0.1])
            rtext(0.3,0.95,'Amplitude: {}'.format(H.tidal_name[n+1+noff]))
            if n==0: ylabel('Amplitude (m)')
            gca().xaxis.grid('on')
            if n==0 and m==0: legend(loc=1)

            #phase
            subplot(4,4,n+5+2*noff)
            bar(xi+dws[0],H.OP[:,n+1+noff]*180/pi,width=dw,facecolor='r',label='Obs')
            bar(xi+dws[1],H.MP[:,0,n+1+noff]*180/pi,width=dw,facecolor='g', label=runs[0])
            if len(runs)==2: bar(xi+dws[2],H.MP[:,1,n+1+noff]*180/pi,width=dw,facecolor='b', label=runs[1])

            plot(arange(-10,10),zeros(20),'k:')
            setp(gca(),xticks=xi,xticklabels=[],xlim=[-0.5,7.5])
            setp(gca(),yticks=arange(-180,200,90),ylim=[-180,190])
            rtext(0.3,0.95,'Phase: {}'.format(H.tidal_name[n+1+noff]))
            if n==0: ylabel('Phase (degree)')
            gca().xaxis.grid('on')
            if m==1:
                setp(gca(),xticks=xi,xticklabels=xls,xlim=[-0.5,7.5])
                xticks(fontsize=10,rotation=-45)

    gcf().tight_layout(); mvfig()
    savefig('elev/elev_Harmonics_1')
    if ishow[0]==0: close()

#------------------------------------------------------------------------------
#temperautre: at main channels
#------------------------------------------------------------------------------
if fexist('salt_temp.npz') and iplot[1]==1:
    if not fexist('temp'): os.mkdir('temp')
    regions=['channel','upper_WT','patuxent','potomac','rappahannock','york','james','ET','sample',] #main stations
    pstations=(
        ['CB1.1','CB2.1','CB2.2','CB3.1','CB3.2','CB3.3C','CB4.1C','CB4.2C','CB4.3C','CB4.4','CB5.1','CB5.2','CB5.3','CB5.4','CB5.5','CB6.1','CB6.2','CB6.3','CB6.4','CB7.1','CB7.2','CB7.3','CB7.4','CB8.1'], #channel
        ['CB1.1','WT1.1','XJG4337','WT2.1','WT3.1','WT4.1','WT5.1','XIE5748','XIE2581','WT6.1','WT7.1','WT8.1'], #upper WT
        ['TF1.5','TF1.6','TF1.7','RET1.1','LE1.1','XDE4587','LE1.2','LE1.3','LE1.4','CB5.1W','CB5.1'], # patuxent
        ['TF2.1','TF2.2','TF2.3','TF2.4','RET2.1','RET2.2','RET2.4','LE2.2','LE2.3','CB5.3'], # potomac
        ['TF3.2','TF3.2A','RPP057.00','TF3.3','RET3.1','RET3.2','LE3.1','LE3.2','LE3.3','LE3.4','LE3.6','LE3.7'], # rappahannock
        ['TF4.2','TF4.4','RET4.1','RET4.2','RET4.3','YRK028.58','LE4.1','YRK015.09','LE4.2','YRK005.40','LE4.3','WE4.2'], # york
        ['TF5.2','TF5.2A','TF5.3','TF5.4','TF5.5','TF5.5A','TF5.6','JMS050.74','RET5.2','LE5.1','LE5.2','LE5.3','LE5.4','CB8.1'], # james
        ['ET2.1','ET2.3','ET2.2','ET4.1','XIH0077','XHH4916','ET4.2','EE1.1','ET5.1','ET5.2','EE2.1','EE2.2','ET6.2','ET7.1','EE3.0','EE3.1','ET8.1','ET9.1','EE3.2','EE3.4'], #ET
        # ['CB3.3C','CB4.4','CB5.4'], #sample stations
    )

    #read data
    models=[read(bdir+'{}/results/salt_temp.npz'.format(i)) for i in runs]; C=read(wqdata)
    [exec('i.temp=i.temp_elem') for i in models if ((not hasattr(i,'temp')) and hasattr(i,'temp_elem'))]

    #labels
    depths=[0,1000]; dtags=['surf','bot']  #depth
    cs=array([['g','k'],['b','m'],['c','r']]); ms=array([['-','-'],['-','-'],['-','-']]); aps=array([[0.9,0.9],[0.6,0.6],[0.95,0.95]])
    xts,xls=get_xtick(xts=[1990,2024])
    xm=[min([i.time.min()+t0 for i,t0 in zip(models,StartTs)])-1,max([i.time.max()+t0 for i,t0 in zip(models,StartTs)])+30]
    if diff(xm)>365*10: xls=[i[-2:] for i in xls] 

    #plots
    TS=zdata(); TS.stations=pstations; TS.regions=regions; TS.stat=[]
    for m,[stations,region] in enumerate(zip(pstations,regions)):
        # if m!=8: continue
        if m==0: figure(figsize=[25,9.8]); fsize=[6,4]
        if m in [1,2,3,4,5]: figure(figsize=[18,9]); fsize=[4,3]
        if m==6: figure(figsize=[20,9]); fsize=[4,4]
        if m==7: figure(figsize=[20,9]); fsize=[5,4]
        if m==8: figure(figsize=[15,8]); fsize=[3,1]

        T=zdata(); ds=[len(stations),2]; t=nan*zeros(ds); T.R,T.ME,T.MAE,T.RMSD=t.copy(),t.copy(),t.copy(),t.copy()
        for i,station in enumerate(stations):
            subplot(*fsize,i+1); lstr=[]
            #models
            for n,run in enumerate(runs):
                for k,[depth,dtag] in enumerate(zip(depths,dtags)):
                    M=models[n]; mti=M.time+StartTs[n]; zi=abs(M.bp.z-depth)
                    sid=nonzero((M.bp.station==station)*(zi==zi.min()))[0][0]; myi=M.temp[sid]
                    plot(mti,myi,ms[n,k],color=cs[n,k],alpha=aps[n,k]); lstr.append('{}: {}'.format(run,dtag))
            #Obs
            fps=(C.station==station)*(C.var=='WTEMP')*(C.time>=xm[0])*(C.time<=xm[1])
            oti=C.time[fps]; oyi=C.data[fps]; plot(oti,oyi,'r.',ms=3)
            lstr.append('Obs. WTEMP')

            #statistics
            sp=M.bp.station; zs=M.bp.z; layers=['S','B']
            for k in arange(2):
                sid=nonzero((sp==station)*(zs==depths[k]))[0][0]; myi=M.temp[sid]
                fp=(C.station==station)*(C.var=='WTEMP')*(C.time>=xm[0])*(C.time<=xm[1])*(C.layer==layers[k]); oti,oyi=C.time[fp],C.data[fp]
                if len(oti)==0: continue
                st=stat(mti,myi,oti,oyi,dt=0.5); T.R[i,k]=st.R; T.ME[i,k]=st.ME; T.MAE[i,k]=st.MAE; T.RMSD[i,k]=st.RMSD

            #note
            if m==0: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>19 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m in [1,4,5]: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>8 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==2: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>7 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==3: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>6 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==6: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>10 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==7: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>14 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==8: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm); ylabel('temperature (oC)'); legend(['Model: Surface','Model: Bottom','Obs'],loc=9)

            gca().xaxis.grid('on')
            ym=ylim(); setp(gca(),ylim=[ym[0],max(ym[1],1)])
            rtext(0.01,0.93,station,fontsize=11,fontweight='bold',color='b')
            if i==0 and m!=8: rtext(0.5,0.93,region,fontsize=11,fontweight='bold',color='b')
            if i==0 and m!=8: hl=legend(lstr,fontsize=8,loc=1)

        gcf().tight_layout()
        savefig('temp/temp_ts_'+region)
        if ishow[1]==0: close()

        #write statistics
        sdata=array([['region: '+region,'', *stations],
                     ['Surface','R',*T.R[:,0]],['','RMSD',*T.RMSD[:,0]],['','ME',*T.ME[:,0]],['','MAE',*T.MAE[:,0]],
                     ['Bottom','R',*T.R[:,1]], ['','RMSD',*T.RMSD[:,1]], ['','ME',*T.ME[:,1]], ['','MAE',*T.MAE[:,1]]],dtype='O')
        write_excel('statistics',sdata,sht='temperature',indy=10*m)
        TS.stat.append(T)
    STS.temp=TS

#------------------------------------------------------------------------------
#salinity: time series, profile
#------------------------------------------------------------------------------
if fexist('salt_temp.npz') and iplot[2]==1:
    if not fexist('salt'): os.mkdir('salt')
    regions=['channel','upper_WT','patuxent','potomac','rappahannock','york','james','ET'] #main stations
    pstations=(
        ['CB1.1','CB2.1','CB2.2','CB3.1','CB3.2','CB3.3C','CB4.1C','CB4.2C','CB4.3C','CB4.4','CB5.1','CB5.2','CB5.3','CB5.4','CB5.5','CB6.1','CB6.2','CB6.3','CB6.4','CB7.1','CB7.2','CB7.3','CB7.4','CB8.1'], #channel
        ['CB1.1','WT1.1','XJG4337','WT2.1','WT3.1','WT4.1','WT5.1','XIE5748','XIE2581','WT6.1','WT7.1','WT8.1'], #upper WT
        ['TF1.5','TF1.6','TF1.7','RET1.1','LE1.1','XDE4587','LE1.2','LE1.3','LE1.4','CB5.1W','CB5.1'], # patuxent
        ['TF2.1','TF2.2','TF2.3','TF2.4','RET2.1','RET2.2','RET2.4','LE2.2','LE2.3','CB5.3'], # potomac
        ['TF3.2','TF3.2A','RPP057.00','TF3.3','RET3.1','RET3.2','LE3.1','LE3.2','LE3.3','LE3.4','LE3.6','LE3.7'], # rappahannock
        ['TF4.2','TF4.4','RET4.1','RET4.2','RET4.3','YRK028.58','LE4.1','YRK015.09','LE4.2','YRK005.40','LE4.3','WE4.2'], # york
        ['TF5.2','TF5.2A','TF5.3','TF5.4','TF5.5','TF5.5A','TF5.6','JMS050.74','RET5.2','LE5.1','LE5.2','LE5.3','LE5.4','CB8.1'], # james
        ['ET2.1','ET2.3','ET2.2','ET4.1','XIH0077','XHH4916','ET4.2','EE1.1','ET5.1','ET5.2','EE2.1','EE2.2','ET6.2','ET7.1','EE3.0','EE3.1','ET8.1','ET9.1','EE3.2','EE3.4'], #ET
    )

    #read data
    models=[loadz(bdir+'{}/results/salt_temp.npz'.format(i)) for i in runs]
    [exec('i.salt=i.salt_elem') for i in models if ((not hasattr(i,'salt')) and hasattr(i,'salt_elem'))]
    C=loadz(wqdata)

    #labels
    depths=[0,1000]; dtags=['surf','bot']  #depth
    cs=array([['g','k'],['b','m'],['c','r']]); ms=array([['-','-'],['-','-'],['-','-']]); aps=array([[0.9,0.9],[0.6,0.6],[0.95,0.95]])
    xts,xls=get_xtick(xts=[1990,2024])
    xm=[min([i.time.min()+t0 for i,t0 in zip(models,StartTs)])-1,max([i.time.max()+t0 for i,t0 in zip(models,StartTs)])+30]
    if diff(xm)>365*10: xls=[i[-2:] for i in xls] 

    #plots
    TS=zdata(); TS.stations=pstations; TS.regions=regions; TS.stat=[]
    for m,[stations,region] in enumerate(zip(pstations,regions)):
        if m==0: figure(figsize=[25,9.8]); fsize=[6,4]
        if m in [1,2,3,4,5]: figure(figsize=[18,9]); fsize=[4,3]
        if m==6: figure(figsize=[20,9]); fsize=[4,4]
        if m==7: figure(figsize=[20,9]); fsize=[5,4]

        T=zdata(); ds=[len(stations),2]; t=nan*zeros(ds); T.R,T.ME,T.MAE,T.RMSD=t.copy(),t.copy(),t.copy(),t.copy()
        for i,station in enumerate(stations):
            subplot(*fsize,i+1); lstr=[]
            #models
            for n,run in enumerate(runs):
                for k,[depth,dtag] in enumerate(zip(depths,dtags)):
                    M=models[n]; mti=M.time+StartTs[n]; zi=abs(M.bp.z-depth)
                    sid=nonzero((array(M.bp.station)==station)*(zi==zi.min()))[0][0]; myi=M.salt[sid]
                    plot(mti,myi,ms[n,k],color=cs[n,k],alpha=aps[n,k]); lstr.append('{}: {}'.format(run,dtag))
            #Obs
            fps=(C.station==station)*(C.var=='SALINITY')*(C.time>=xm[0])*(C.time<=xm[1])
            oti=C.time[fps]; oyi=C.data[fps]; plot(oti,oyi,'r.',ms=3)
            lstr.append('Obs. SALINITY')

            #statistics
            sp=array(M.bp.station); zs=M.bp.z; layers=['S','B']
            for k in arange(2):
                sid=nonzero((sp==station)*(zs==depths[k]))[0][0]; myi=M.salt[sid]
                fp=(C.station==station)*(C.var=='SALINITY')*(C.time>=xm[0])*(C.time<=xm[1])*(C.layer==layers[k]); oti,oyi=C.time[fp],C.data[fp]
                if len(oti)==0: continue
                st=stat(mti,myi,oti,oyi,dt=0.5); T.R[i,k]=st.R; T.ME[i,k]=st.ME; T.MAE[i,k]=st.MAE; T.RMSD[i,k]=st.RMSD

            #note
            if m==0: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>19 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m in [1,4,5]: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>8 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==2: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>7 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==3: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>6 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==6: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>10 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            if m==7: setp(gca(),xticks=xts,xticklabels=xls,xlim=xm) if i>14 else setp(gca(),xticks=xts,xticklabels=[],xlim=xm)

            gca().xaxis.grid('on')
            ym=ylim(); setp(gca(),ylim=[ym[0],max(ym[1],1)])
            rtext(0.01,0.93,station,fontsize=11,fontweight='bold',color='b')
            if i==0: rtext(0.5,0.93,region,fontsize=11,fontweight='bold',color='b')
            if i==0: hl=legend(lstr,fontsize=8,loc=1)

        gcf().tight_layout()
        savefig('salt/salinity_ts_'+region)
        if ishow[2]==0: close()

        #write statistics
        sdata=array([['region: '+region,'', *stations],
                     ['Surface','R',*T.R[:,0]],['','RMSD',*T.RMSD[:,0]],['','ME',*T.ME[:,0]],['','MAE',*T.MAE[:,0]],
                     ['Bottom','R',*T.R[:,1]], ['','RMSD',*T.RMSD[:,1]], ['','ME',*T.ME[:,1]], ['','MAE',*T.MAE[:,1]]],dtype='O')
        write_excel('statistics',sdata,sht='salinity',indy=10*m)
        TS.stat.append(T)

        #salinity: plot for statistics
        figure(figsize=[11,6]); xi=arange(T.R.shape[0]); tdata=[T.R,T.RMSD,T.ME,T.MAE]
        yms=[[0,1],[0,3],[-1,3],[0,3]]; tstr=['R','RMSD','ME','MAE'];
        tts=['Correlation Coefficients (R)','Root Mean Square Deviation (PSU)','Mean Error (PSU)','Mean Absolute Error (PSU)']
        for i in arange(4):
            subplot(2,2,i+1)
            plot(xi,tdata[i][:,0],'g.-',ms=8); plot(xi,tdata[i][:,1],'r.-',ms=8)
            setp(gca(),xticks=xi,xticklabels=stations,xlim=[-0.5,xi.max()+0.5],ylim=yms[i])
            legend([tstr[i]+': Surface',tstr[i]+': Bottom']); title(tts[i])
            gca().xaxis.grid('on')
            xticks(fontsize=8,rotation=45)
        gcf().tight_layout()
        savefig('salt/salinity_stat_'+region)
    STS.salt=TS

    #------------------------------------------------------------------------------
    #salinity: profiles @channel stations (note: only for last run)
    #------------------------------------------------------------------------------
    #main channel stations
    stations=['CB2.2','CB3.1','CB3.2','CB3.3C','CB4.1C','CB4.2C','CB4.3C','CB4.4','CB5.1','CB5.2',
              'CB5.3','CB5.4','CB5.5','CB6.1','CB6.2','CB6.3','CB6.4','CB7.1','CB7.2','CB7.3','CB7.4','CB8.1']
    xms=[[0,10],[0,15],[0,20],[0,20],[5,25],[5,25],[5,25],[5,25],[5,25],[5,25],
          [5,25],[10,27],[10,30],[10,25],[13,25],[13,25],[15,30],[15,30],[15,30],[18,33],
          [20,33],[18,35]]
    yms=[[0,13],[0,13],[0,13],[0,25],[0,35],[0,30],[0,30],[0,35],[0,35],[0,35],
          [0,30],[0,35], [0,20], [0,15], [0,12], [0,12], [0,12], [0,25], [0,22],[0,17],
          [0,17],[0,10]]

    #read data
    M=loadz(bdir+'{}/results/salt_profile.npz'.format(runs[-1])); C=loadz(wqdata)
    if ((not hasattr(M,'salt')) and hasattr(M,'salt_elem')): M.salt=M.salt_elem
    gd=loadz(grd).hgrid; vd=loadz(grd).vgrid; zcor=vd.compute_zcor(gd.dp)
    pie,pip,pacor=gd.compute_acor(c_[M.bp.x,M.bp.y]); mz=(zcor[pip]*pacor[...,None]).sum(axis=1)

    #plots
    mti=M.time+StartTs[-1]; xts=arange(0,40,5); xls=xts; yts=arange(-50,5,5); yls=abs(yts)
    for m, station in enumerate(stations):
        figure(figsize=[30,9.8]); fsize=[5,20]; xm=xms[m]; ym=yms[m]

        #find profiles
        fp=(C.var=='SALINITY')*(C.station==station)*(C.time>=mti[0])*(C.time<=mti[-1])
        oti,oyi,odi=C.time[fp],C.data[fp],C.depth[fp]; uti=unique(oti)

        #get station id
        sid=nonzero(M.bp.station==station)[0][0]; mzi=mz[sid]

        #plot each profile
        for i in arange(prod(fsize)):
            subplot(*fsize,i+1)
            if i<len(uti):
                utii=uti[i]
                #model
                fp=abs(mti-utii)<0.5; myi=M.salt[sid,fp]
                plot(myi.T,mzi,'k-',lw=0.2,alpha=0.5)
                plot(myi.mean(axis=0),mzi,'g')

                #obs
                fp=oti==utii;
                plot(oyi[fp],-odi[fp],'r.')
            else:
                plot(arange(10),arange(10)*nan,'k')

            #note
            if mod(i,20)==0:
                setp(gca(),yticks=yts,yticklabels=yls); ylabel('Depth (m)')
            else:
                setp(gca(),yticks=yts,yticklabels=[])

            if i>79:
                setp(gca(),xticks=xts,xticklabels=xls)
            else:
                setp(gca(),xticks=xts,xticklabels=[])
            setp(gca(),xlim=[xm[0]-0.1,xm[1]+0.1],ylim=[-ym[1]-0.1,ym[0]+0.1,])
            gca().xaxis.grid('on'); gca().yaxis.grid('on');

            if i<len(uti): rtext(0.25,0.95,'{}: {}'.format(station,num2date(utii).strftime('%Y-%m-%d')),fontsize=7,color='b')
            # if i==0: rtext(0.4,0.9,station,fontsize=8, fontweight='bold',color='b')

        subplots_adjust(top=0.98,bottom=0.04,left=0.021,right=0.995,hspace=0.063,wspace=0.092)
        if not fexist('salt'): os.mkdir('salt')
        savefig('salt/salinity_profile_'+station+'.png')
        if ishow[2]==0: close()

#------------------------------------------------------------------------------
#sediment
#------------------------------------------------------------------------------
if fexist('sed_profile.npz') and iplot[3]==1:
    if not fexist('sed'): os.mkdir('sed')
    #read data
    models=[loadz(bdir+'{}/results/sed_profile.npz'.format(i)) for i in runs]; mts=[i.time+t0 for i,t0 in zip(models,StartTs)]
    C=loadz(wqdata)

    xts,xls=get_xtick(xts=[1990,2024])
    xm=[min([i.min() for i in mts]),max([i.max() for i in mts])] #range of model runs
    #xm=[datenum(1991,1,1),datenum(1994,5,1)]
    mt0=max([i.min()+10 for i in mts]); mt1=mt0+15   #range for plotting tide
    htm=[max([i.min() for i in mts]),min([i.max() for i in mts])] #range for HA
    if diff(xm)>365*10: xls=[i[-2:] for i in xls] 

    #------------------------------------------------------------------------------
    #sediment: time series in CB stations
    #------------------------------------------------------------------------------
    #inputs
    pstations=(['CB1.1','CB2.1','CB2.2','CB3.1','CB3.2','CB3.3C','CB4.1C','CB4.2C','CB4.3C','CB4.4','CB5.1','CB5.2','CB5.3','CB5.4','CB5.5','CB6.1', 'CB6.2','CB6.3','CB6.4','CB7.1','CB7.2','CB7.3','CB7.4','CB8.1'],
        ['CB3.3E','CB3.3W','CB4.1E','CB4.1W','CB4.2E','CB4.2W','CB4.3E','CB4.3W','CB5.1W','CB5.4W','CB7.1N','CB7.1S','CB7.2E','CB7.3E','CB7.4N','CB8.1E'])

    #plot surface and bottom TSS
    for i,stations in enumerate(pstations):

        if i==0: figure(figsize=[30,9.4]); fsize=[5,5]
        if i==1: figure(figsize=[19,9.2]); fsize=[4,4]
        for m in arange(prod(fsize)):
            subplot(*fsize,m+1)
            if m<len(stations):
                #get bottom and surface
                cs=[['g','b'],['k','m']]; station=stations[m]; lstr=[]
                for n,run in enumerate(runs):
                    M=models[n]; sid=list(M.bp.station).index(station)
                    if hasattr(M,'sed_tconc'):
                        myi=1e3*M.sed_tconc[sid]
                    else:
                        myi=1e3*(M.sed_conc1[sid]+M.sed_conc2[sid]+M.sed_conc3[sid]+M.sed_conc4[sid])
                    mti=M.time+StartTs[n]; myis=myi[:,-1]; myib=myi[:,0]
                    plot(mti,smooth(myis,3*24),color=cs[n][0],lw=0.5,alpha=1); lstr.append('{}_surf'.format(run))
                    plot(mti,smooth(myib,3*24),color=cs[n][1],lw=0.5,alpha=1); lstr.append('{}_bott'.format(run))

                # plot obs.
                #fp=(C.time>=xm[0])*(C.time<=xm[-1])*(C.var=='TSS')*(C.station==station) #*(C.depth<2.0)
                #oti=C.time[fp]; oyi=C.data[fp]; plot(oti,oyi,'r.',ms=5); lstr.append('Obs_TSS')
                fp=(C.time>=xm[0])*(C.time<=xm[-1])*(C.var=='TSS')*(C.station==station)*(C.layer=='S')
                oti=C.time[fp]; oyi=C.data[fp]; plot(oti,oyi,'k^',ms=3); lstr.append('Obs: TSS_Surf')
                fp=(C.time>=xm[0])*(C.time<=xm[-1])*(C.var=='TSS')*(C.station==station)*(C.layer=='B')
                oti=C.time[fp]; oyi=C.data[fp]; plot(oti,oyi,'r.',ms=5); lstr.append('Obs: TSS_Bott')

            else:
                plot(arange(10),zeros(10)*nan)

            #note
            ym=[0,50]
            if m<(prod(fsize)-fsize[1]):
                setp(gca(),xticks=xts,xticklabels=[],xlim=xm)
            else:
                setp(gca(),xticks=xts,xticklabels=xls,xlim=xm)
            if m%fsize[1]==0:
                setp(gca(),yticks=arange(0,100,10),ylim=ym)
                ylabel('TSS (mg/L)')
            else:
                setp(gca(),yticks=arange(0,100,10),yticklabels=[],ylim=ym)
            if m<len(stations): rtext(0.02,0.9,station,fontweight='bold',fontsize=12,color='k')
            if m==0: legend(lstr,loc=9)
            gca().xaxis.grid('on')

        gcf().tight_layout(); mvfig()
        savefig('sed/sed_TS_{}'.format(i+1))
        if ishow[3]==0: close()

    #------------------------------------------------------------------------------
    #plot sediment variables at each station
    #------------------------------------------------------------------------------
    fname=bdir+'{}/results/sed_profile_ext.npz'.format(runs[-1])
    if fexist(fname):
       stations=[*pstations[0],*pstations[1]]
       M=read(fname); mti=M.time+StartTs[-1]
       for m,station in enumerate(stations):
           figure(figsize=[16,8]); sid=list(M.bp.station).index(station)

           lstr=['TSS','sed_1','sed_2','sed_3','sed_4',['frac1','frac2','frac3','frac4'],'shear stress','roughness','depositional flux','erosion flux']
           for i in arange(10):
               if i==0: #TSS
                   myi=1e3*(M.sed_conc1[sid]+M.sed_conc2[sid]+M.sed_conc3[sid]+M.sed_conc4[sid])[:,-1]
               elif i==1: #sed_con
                       myi=M.sed_conc1[sid,:,-1]*1e3
               elif i==2: #sed_con
                       myi=M.sed_conc2[sid,:,-1]*1e3
               elif i==3: #sed_con
                       myi=M.sed_conc3[sid,:,-1]*1e3
               elif i==4: #sed_con
                       myi=M.sed_conc4[sid,:,-1]*1e3
               elif i==5: #sed frac
                   myi=c_[M.sed_frac1[sid,:],M.sed_frac2[sid,:],M.sed_frac3[sid,:],M.sed_frac4[sid,:]]
               elif i==6: #stress
                   myi=M.sed_str[sid]
               elif i==7: #rough
                   myi=M.sed_rough[sid]
               elif i==8: #dflux
                   myi=M.sed_dflux[sid]
               elif i==9: #eflux
                   myi=M.sed_eflux[sid]
               else:
                   continue
               if i<=4: myi[abs(myi)>1e5]=nan

               #surface TSS and sed
               cs='krgb'
               subplot(3,4,i+1)
               if myi.ndim==1:
                   if i>=1 and i<=4:
                       plot(mti,myi,color=cs[i-1])
                   else:
                       plot(mti,myi,'k')
               else:
                   for k in arange(myi.shape[1]): plot(mti,myi[:,k],color=cs[k])

               #note
               setp(gca(),xticks=xts,xticklabels=xls,xlim=xm)
               if i==5:
                   legend(lstr[i],loc=1)
               else:
                   legend([lstr[i]],loc=1)
               gca().xaxis.grid('on')

           gcf().tight_layout(); mvfig()
           savefig('sed/sed_station_{}.png'.format(station))
           if ishow[3]==0: close()

   #------------------------------------------------------------------------------
   #plot sediment pofiles
   #------------------------------------------------------------------------------
   # M=models[-1]; mti=mts[-1]; years=arange(num2date(mti.min()).year,num2date(mti.max()).year+1)

   # gd=loadz(grd).hgrid


