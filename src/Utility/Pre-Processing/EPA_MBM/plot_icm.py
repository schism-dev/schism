#!/usr/bin/env python3
from pylib import *
close("all")

#------------------------------------------------------------------------------
#inputs
#------------------------------------------------------------------------------
runs=['CH3D','RUN12f']; StartTs=[datenum(1991,1,1),datenum(1991,1,1),]
xm=[datenum(1991,1,1),datenum(2001,1,1)]   #time range for plotting

#variables
mvars=['temp', 'salt',    'chla','DO','TN','TP','NO3',  'NH4', 'PO4', 'DOC','POC','PON','POP','RPOC','LPOC','RPON','LPON','RPOP','LPOP','PB1', 'PB2', 'PB3', 'DON','DOP'] #model variables
ovars=['WTEMP','SALINITY','CHLA','DO','TN','TP','NO23F','NH4F','PO4F','DOC','PC', 'PN', 'PP', 'PC',  'PC',  'PN',   'PN', 'PP',  'PP',  'CHLA','CHLA','CHLA','DON','DOP'] #observation

#stations
regions=['channel','shoal','WT','ET','patuxent','potomac','rappahannock','york','james',] #main stations
pstations=(['CB1.1','CB2.1','CB2.2','CB3.1','CB3.2','CB3.3C','CB4.1C','CB4.2C','CB4.3C','CB4.4','CB5.1','CB5.2','CB5.3','CB5.4','CB5.5','CB6.1','CB6.2','CB6.3','CB6.4','CB7.1','CB7.2','CB7.3','CB7.4','CB8.1'], #channel
    ['CB3.3E','CB3.3W','CB4.1E','CB4.1W','CB4.2E','CB4.2W','CB4.3E','CB4.3W','CB5.1W','CB5.4W','CB7.1N','CB7.1S','CB7.2E','CB7.3E','CB7.4N','CB8.1E'], #shoal
    ['WT1.1','XJG4337','WT2.1','WT3.1','WT4.1','WT5.1','XIE5748','XIE2581','WT6.1','WT7.1','WT8.1'], #upper WT
    ['ET2.1','ET2.3','ET2.2','ET4.1','XIH0077','XHH4916','ET4.2','EE1.1','ET5.1','ET5.2','EE2.1','EE2.2','ET6.2','ET7.1','EE3.0','EE3.1','ET8.1','ET9.1','EE3.2','EE3.4'], #ET
    ['TF1.5','TF1.6','TF1.7','RET1.1','LE1.1','XDE4587','LE1.2','LE1.3','LE1.4','CB5.1W','CB5.1'], # patuxent
    ['TF2.1','TF2.2','TF2.3','TF2.4','RET2.1','RET2.2','RET2.4','LE2.2','LE2.3','CB5.3'], # potomac
    ['TF3.2','TF3.2A','RPP057.00','TF3.3','RET3.1','RET3.2','LE3.1','LE3.2','LE3.3','LE3.4','LE3.6','LE3.7'], # rappahannock
    ['TF4.2','TF4.4','RET4.1','RET4.2','RET4.3','YRK028.58','LE4.1','YRK015.09','LE4.2','YRK005.40','LE4.3','WE4.2'], # york
    ['TF5.2','TF5.2A','TF5.3','TF5.4','TF5.5','TF5.5A','TF5.6','JMS050.74','RET5.2','LE5.1','LE5.2','LE5.3','LE5.4','CB8.1'],) # james

#more details
layers=['S','B']  #surface and bottom
iflags=[1, 0, 1] # TS(1 var, stations), TS(vars, 1 station), statistics
bdir='/sciclone/data10/wangzg/CBP/'
#------------------------------------------------------------------------------
#read data and obs info
#------------------------------------------------------------------------------
fnames=[bdir+ ('Data/CH3D/CH3D.npz' if i=='CH3D' else i+'/results/icm.npz') for i in runs]; C=read(bdir+'Data/cbp/CBP_WQData.npz')
models=[]; svars=[]; mtimes=[]; mzs=[]; mstations=[]
for run,StartT,fname in zip(runs,StartTs,fnames):
    M=read(fname,1); svar=M.attr(); mtime=M.time[:] if run=='CH3D' else M.time[:]+StartT
    mz,mstation=[M.bp.z,M.bp.station] if M.hasattr('bp') else [M.z[:],M.station[:]]
    #add conversion coefficients
    M.c2chl=[0.75, 0.06, 0.06]; M.n2c=[0.167, 0.167, 0.167]; M.p2c=[0.0125, 0.0125, 0.0125] #default values
    for i in ['c2chl','n2c','p2c']: M.attr(i,array(M.icm.attr(i) if M.hasattr('icm') else M.param.attr(i) if M.hasattr('param') else M.attr(i))); 
    models.append(M); svars.append(svar); mtimes.append(mtime); mzs.append(mz); mstations.append(mstation)

def get_icm_var(model,mvar):
    M=model; mdata=None; 
    if M.hasattr(mvar): return M.attr(mvar) if (mvar not in ['PB1','PB2','PB3']) else M.attr(mvar)/M.c2chl[int(mvar[-1])-1] #variable exist
    #original variable not exist (for chlorophyll, POC,PON,POP,DO TN, TP, etc)
    if mvar=='chla':  mdata=M.chla if M.hasattr('chla') else M.CHLA if M.hasattr('CHLA') else (c_[M.attr(['PB1','PB2','PB3'])]/M.c2chl[:,None,None]).sum(axis=0)
    if mvar in ['POC','PON','POP']: mdata=M.attr('R'+mvar)+M.attr('L'+mvar)
    if (mvar[:4] in ['bPOC','bPON','bPOP']) and M.hasattr(mvar[:4]): mvari=mvar[:4]; im=int(mvar[4])-1;  mdata=M.attr(mvari)[...,im]
    if mvar=='DO' and fname.endswith('icm.nc'): mdata=M.DOX
    if mvar=='TN' and mdata==None: mdata=M.attr(['RPON','LPON','DON','NH4','NO3'],fmt=1).sum(axis=0)+(M.attr(['PB1','PB2','PB3'],fmt=1)*M.n2c[:,None,None]).sum(axis=0)+M.attr('SRPON',default=0)
    if mvar=='TP' and mdata==None: mdata=M.attr(['RPOP','LPOP','DOP','PO4'],fmt=1).sum(axis=0)+(M.attr(['PB1','PB2','PB3'],fmt=1)*M.p2c[:,None,None]).sum(axis=0)+M.attr('SRPOP',default=0)+M.attr('PIP',default=0)
    if mvar=='TSS': mdata=M.attr('sed_tconc')*1e3
    return mdata

#------------------------------------------------------------------------------
#plot each variable time series
#------------------------------------------------------------------------------
#time ticks, and time style
ys=[num2date(i).year for i in xm]; tm=[max([i.min() for i in mtimes]),min([i.max() for i in mtimes])]
xts,xls=xtick(ys,0 if diff(ys)>1 else 1); xls[0]=str(ys[0]); xls=[i[2:] if diff(ys)>1 else i for i in xls]
cs=array([['g','k'],['b','m'],['c','r']]); ms=array([['-','-'],['-','-'],['-','-']]); aps=array([[0.9,0.9],[0.75,0.75],[0.95,0.95]])

#------------------------------------------------------------------------------
#plot variable at stations based on regions
#------------------------------------------------------------------------------
if iflags[0]==1:
   S=zdata(); st_vars=['run','region','layer','station','var','R','ME','MAE','RMSD']; [S.attr(i,[]) for i in st_vars] #init. statistics capsule
   for m, [mvar,ovar] in enumerate(zip(mvars,ovars)):
       #read model variable, obs data
       mys=[]; icheck=0; t0=time.time()
       for n,[run,fname,M,mtime,mstation] in enumerate(zip(runs,fnames,models,mtimes,mstations)):
           my=get_icm_var(M,mvar); mys.append(my); icheck=icheck+(my is None)
       if icheck==len(runs): continue # no available data
       fp=(C.var==ovar)*((C.layer=='S')|(C.layer=='B')); otime,ostation,odata,olayer=C.time[fp],C.station[fp],C.data[fp],C.layer[fp] #obs

       #for each region and each station
       for nn,[region, stations] in enumerate(zip(regions,pstations)):
           fsizes=[[19,9.4],[19,9.4],[19,9.4],[19,9.4],[19,9.4],[19,9.4],[19,9.4],[19,9.4],[19,9.4]]
           fsubs= [[5,5],   [4,5],   [3,4],   [5,5],   [3,4],   [3,4],   [3,5],   [3,5],   [3,5]]
           hf=figure(figsize=fsizes[nn])
           print('plotting time series: {}, {}'.format(mvar,region))
           for i,station in enumerate(stations):
               subplot(*fsubs[nn],i+1); lstr=[]
               fp=(ostation==station)*(otime>=xm[0])*(otime<=xm[1])*(olayer=='S'); ots=otime[fp]; oys=odata[fp] #get surf. obs
               fp=(ostation==station)*(otime>=xm[0])*(otime<=xm[1])*(olayer=='B'); otb=otime[fp]; oyb=odata[fp] #get bott. obs
               for n,[run,mtime,my,mz,mstation] in enumerate(zip(runs,mtimes,mys,mzs,mstations)): #plot model
                   if my is None: continue
                   #get model variable
                   for k, layer in enumerate(layers):
                       sids=pindex(mstation==station); zs=mz[sids] 
                       fps=mstation==station; sids=nonzero(fps)[0]; zs=mz[fps]; sid=sids[argmin(zs) if layer=='S' else argmax(zs)]
                       myi=my[sid]; myi[myi<-1e4]=nan; tag='surface' if layer=='S' else 'bottom'
                       plot(mtime,myi,ms[n,k],color=cs[n,k],alpha=aps[n,k]); lstr.append('{}: {}'.format(run,tag))
   
                       #do statistics
                       if iflags[2]==1:
                          if (mvar in ['PB1','PB2','PB3','LPOC','RPOC','LPON','RPON','LPOP','RPOP']) or len(otime)==0: continue
                          ot,oy=[ots,oys] if layer=='S' else [otb,oyb]
                          if len(ot)!=0: fp=(ot>=tm[0])*(ot<=tm[1]); ot,oy=ot[fp],oy[fp];  fmyi=[]
                          if len(ot)==0: continue #no obs. data in this period
                          if sum(~isnan(myi))==0:
                             [S.attr(i).append(nan) for i in ['R','ME','MAE','RMSD']]
                          else:
                             for oti,oyi in zip(ot,oy): fp=(abs(mtime-oti)<=0.5)*(~isnan(myi)); myii=myi[fp]; fmyi.append(myii[argmin(abs(myii-oyi))]) #find best match
                             s=get_stat(array(fmyi),oy); [S.attr(i).append(s.attr(i)) for i in ['R','ME','MAE','RMSD']]; 
                          S.run.append(run); S.region.append(region); S.layer.append(layer); S.station.append(station); S.var.append(mvar)
               plot(ots,oys,'r^',ms=2); lstr.append('Obs: {}_surface'.format(mvar)) #plot surf. obs
               plot(otb,oyb,linestyle='None',color='orange',marker='s',ms=2); lstr.append('Obs: {}_bottom'.format(mvar)) #plot bott. obs
       
               #note
               setp(gca(),xticks=xts,xticklabels=xls,xlim=xm); gca().xaxis.grid('on')
               text(xlim()[0]+0.03*diff(xlim()),ylim()[0]+0.8*diff(ylim()),station,fontsize=12,fontweight='bold')
               if i==0: hl=legend(lstr)
           gcf().tight_layout()
           hl.set_bbox_to_anchor([0.82,0.075,0.1,0.1],transform=gcf().transFigure)
        
           #save fig
           sdir='_'.join(runs); sname=sdir+'/icm_{}_{}'.format(region,mvar); 
           if not fexist(sdir): os.mkdir(sdir)
           #show(block=False); sys.exit()
           #savefig(sname+'.pp')
           savefig(sname,dpi=300); close()
   if iflags[2]==1: 
      S.x,S.y=proj_pts(*array([[C.lon[i],C.lat[i]] for i in S.station]).T,'epsg:4326','epsg:26918')
      S.to_array(); S.save('stat_'+'_'.join(runs))
            
#------------------------------------------------------------------------------
#Plot all variables for each station
#------------------------------------------------------------------------------
if iflags[1]==1:
   for nn,[region, stations] in enumerate(zip(regions,pstations)):
       for i,station in enumerate(stations):
           figure(figsize=(28,14))
           print('plotting '+ region+','+ station)
           for m, [mvar,ovar] in enumerate(zip(mvars,ovars)):
               #read model variable, obs data
               mys=[]; icheck=0; t0=time.time()
               for n,[run,fname,M,mtime,mstation] in enumerate(zip(runs,fnames,models,mtimes,mstations)):
                   my=get_icm_var(M,mvar); mys.append(my); icheck=icheck+(my is None)
               if icheck==len(runs): continue # no available data
               fps=(C.station==station)*(C.var==ovar)*(C.time>=xm[0])*(C.time<=xm[1])*(C.layer=='S') #get obs
               fpb=(C.station==station)*(C.var==ovar)*(C.time>=xm[0])*(C.time<=xm[1])*(C.layer=='B')
               ots=C.time[fps]; oys=C.data[fps]; otb=C.time[fpb]; oyb=C.data[fpb]
   
               subplot(5,5,m+1)
               #plot model
               lstr=[]
               for n,[run,mtime,my,mz,mstation] in enumerate(zip(runs,mtimes,mys,mzs,mstations)):
                   if my is None: continue
                   #get model variable
                   for k, layer in enumerate(layers):
                       fps=mstation==station; sids=nonzero(fps)[0]; zs=mz[fps]
                       sid=sids[argmin(zs) if layer=='S' else argmax(zs)];  tag='surface' if layer=='S' else 'bottom'
                       myi = my[sid]
                       myi[myi<-1e4]=nan
   
                       plot(mtime,myi,ms[n,k],color=cs[n,k],alpha=aps[n,k], linewidth=1.5)
                       lstr.append('{}: {}'.format(run,tag))
               #plot obs
               plot(ots,oys,'r^',ms=3); lstr.append('Obs: surface')
               plot(otb,oyb,'s',color='orange',ms=3); lstr.append('Obs: bottom')
   
               #note
               setp(gca(),xticks=xts,xticklabels=xls,xlim=xm); gca().xaxis.grid('on')
               text(xlim()[0]+0.03*diff(xlim()),ylim()[0]+0.9*diff(ylim()),mvar+': '+station,fontsize=12,fontweight='bold')
               if m==0: hl=legend(lstr)
   
           gcf().tight_layout()
           hl.set_bbox_to_anchor([0.82,0.075,0.1,0.1],transform=gcf().transFigure)
    
           #save fig
           sdir='_'.join(runs); sname=sdir+'/icm_{}_station_{}.png'.format(region,station)
           if not fexist(sdir): os.mkdir(sdir)
           savefig(sname, dpi=300);  close()

#------------------------------------------------------------------------------
#write statistics excel files
#------------------------------------------------------------------------------
if iflags[2]!=0:
   S=read('stat_{}.npz'.format('_'.join(runs))); tname='statistics_'+'_'.join(runs)+'.xlsx'
   xvars=['temp','salt','chla', 'DO', 'TN', 'TP','NH4', 'NO3', 'PO4', 'DOC', 'DON', 'DOP', 'POC', 'PON', 'POP','TSS']
   vms  =[1.5,     2,     10,    2,    1,   0.1,  0.1,  0.4,   0.02 ,  2,    0.3,   0.02,  2,     0.2,   0.05, 10]
   xvars,vms=array([[i,k] for i,k in zip(xvars,vms) if i in unique(S.var)]).T; vms=vms.astype('float')
   if fexist(tname): os.remove(tname)
   hf=figure()

   #plot scatter plots
   gd=grd('..'); idy=0; fw,fh=[12.5,12]; sht='summary'
   xm=(280714.1267497268, 444350.75124651677); ym=(4073471.7060521273, 4393172.431081811)
   for m,mvar in enumerate(xvars): 
       print('summary plot: ', mvar)
       hf.clf(); hf.set_figwidth(fw); hf.set_figheight(fh); fn=hf
       for n, run in enumerate(runs):
           for k, layer in enumerate(layers):  
               fp=(S.var==mvar)*(S.run==run)*(S.layer==layer); sx,sy,RMSD,ME=S.x[fp],S.y[fp],S.RMSD[fp],S.ME[fp] 
               vm=vms[m]; srat=200/vm

               #plots
               subplot(2,len(runs),2*n+k+1)
               gd.plot_bnd(c='g',lw=0.2); 
               if sum(~isnan(RMSD))!=0 and sum(~isnan(ME))!=0:
                  hg=scatter(sx,sy,s=srat*RMSD,c=ME,vmin=-vm,vmax=vm,cmap='bwr',edgecolor='k',linewidth=0.3); 
                  hb=colorbar(aspect=50)
               setp(gca(),xticks=[],yticks=[],xlim=xm,ylim=ym)
                       
               #note                        
               title('RMSD+Bias: {}, layer={}, {}'.format(run,layer,mvar))
               unit=r'$^o$C' if mvar=='temp' else 'PSU' if mvar=='salt' else 'ug/L' if mvar=='chla' else 'mg/L'
               hl=legend(*hg.legend_elements("sizes", num=[srat*vm]),fontsize=14) #legend
               hl.texts[0].set_text('{} ({})'.format(vm,unit))  #set legend value             
               
               rtext(0.05,0.93,'mean(RMSD)={:0.3f}'.format(mean(RMSD)),fontsize=14)
               rtext(0.05,0.90,'mean(Bias)={:0.3f}'.format(mean(ME)),fontsize=14)
       gcf().tight_layout()
       #show(block=False); sys.exit()
       write_excel(tname,mvar,sht,indy=idy+25,indx=0,fontsize=18,fontweight='bold',color='k')
       write_excel(tname,fn,sht,fmt=3,indy=idy+1,indx=3,figsize=[80*fw,80*fh]); idy=idy+50

   #plot barplot
   if iflags[2]==2:
      for mvar in xvars: #adjust the order 
          idy=0; unit='oC' if mvar=='temp' else 'mg/L'
          for region,stations in zip(regions,pstations): 
              print('barplot: ', mvar,region)
              for layer in ['S','B']:
                  #get header, title
                  tag='Surface' if layer=='S' else 'Bottom'; header1='{}: {} @{} ({})'.format(region,mvar,tag,unit); header2=['Station']
                  for i in ['R2','Bias','RMSE']:
                      for run in runs: header2.append('{} ({})'.format(i,run if run=='CH3D' else 'MBM'))

                  #get statistics data
                  CC=[]; ME=[]; RMSD=[] #orangize stat data
                  for run in runs:
                      CCi=[]; MEi=[]; RMSDi=[]
                      for station in stations:
                          fp=(S.run==run)*(S.var==mvar)*(S.station==station)*(S.layer==layer)
                          if sum(fp)==0:
                             CCi.append(nan); MEi.append(nan); RMSDi.append(nan)
                          else:
                             fp=nonzero(fp)[0][0]; CCi.append(S.R[fp]); MEi.append(S.ME[fp]); RMSDi.append(S.RMSD[fp])
                      CC.append(CCi); ME.append(MEi); RMSD.append(RMSDi)
                  CC=squeeze(array(CC))**2; ME=squeeze(array(ME)); RMSD=squeeze(array(RMSD))
                  #fpn=~isnan(ME.sum(axis=0)); stations,CC,ME,RMSD=array(stations)[fpn],CC[:,fpn],ME[:,fpn],RMSD[:,fpn] #remove nan data
                  data=r_[CC,ME,RMSD].T; mdata=data.mean(axis=0)
                   
                  #write excel data
                  write_excel(tname,header1,mvar,indy=idy,indx=1,fontsize=12,fontweight='bold',color='k') 
                  write_excel(tname,header2,mvar,indy=idy+1,fontsize=10,fontweight='bold',color='b') 
                  write_excel(tname,stations,mvar,indy=idy+2,align='column')
                  write_excel(tname,'mean',mvar,indy=idy+2+len(stations),indx=0,fontweight='bold')
                  write_excel(tname,data,mvar,indy=idy+2,indx=1,number_format="0.0000")
                  write_excel(tname,mdata,mvar,indy=idy+2+len(data),indx=1,number_format="0.0000",fontweight='bold')

                  #insert figure
                  ns=len(stations); xi=arange(ns); dw=0.35
                  fs=[2*max(24*ns,400),max(22*ns,240)]; fsize=[2*max(ns*0.5,8),max(ns*0.5,5)]
                  hf.clf(); hf.set_figwidth(fsize[0]); hf.set_figheight(fsize[1]); fn=hf

                  for mm in arange(3):
                      subplot(2,2,mm+1); lstr=[]; tstr=[]
                      for m, run in enumerate(runs):
                          if mm==0: cdata=ME[m]; stag='Bias'
                          if mm==1: cdata=CC[m]; stag='R2'
                          if mm==2: cdata=RMSD[m]; stag='RMSE'
                          lstr.append(run); tstr.append('{}({})={:0.4f}'.format(stag,run,mean(cdata)))
                          bar(xi-dw/2+m*dw,cdata, width=0.35)
                      if mm==1: setp(gca(),ylim=[0,1])
                      legend(lstr,loc=1, fontsize=12)
                      if mm==0: plot([-1,100],[0,0],'k:')
                      setp(gca(),xticks=xi,xticklabels=stations,xlim=[-0.5,max(ns,8)]); xticks(fontsize=12, rotation=45)
                      ylabel(stag,fontsize=12); title(', '.join(tstr),fontsize=14)
                  gcf().tight_layout()

                  write_excel(tname,fn,mvar,fmt=3,indy=idy+2,indx=9,figsize=fs)
                  idy=idy+5+max(len(data),10)

