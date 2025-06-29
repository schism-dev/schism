#!/usr/bin/env python3
'''
Template for running MPI job on sciclone/ches 
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *
from mpi4py import MPI

#-----------------------------------------------------------------------------
#Input
#hpc: kuro, femto, bora, potomac, james, frontera, levante, stampede2
#ppn:  64,   32,    20,    12,     20,     56,      128,      48
#-----------------------------------------------------------------------------
#P7 loading versions
#p7_version ='20241101_P7beta_wPSX'; p7_num='v1'  #zip names of P7 loading, and #version number
p7_version ='20250101_P7beta'; p7_num='v2'  #zip names of P7 loading, and #version number
iunzip     = 1             #unzip files

#walltime='10:00:00'; nnode=5;  ppn=64
walltime='10:00:00'; nnode=3;  ppn=64

#optional: (frontera,levante,stampede2)
ibatch     =1              #0: serial mode;  1: parallel mode
qnode      =None           #specify node name, or default qnode based on HOST will be used
qname      =None           #partition name
account    =None           #account name
#reservation=None           #reservation information
reservation='wangzg_91'           #reservation information
scrout     ='screen.out'   #fname for outout and error
jname      ='P7_read'      #job name


#-----------------------------------------------------------------------------
#on front node: 1). submit jobs  2) running parallel jobs
#-----------------------------------------------------------------------------
bdir=os.path.abspath(os.path.curdir)
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt,'param',qname,account,reservation)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
if ibatch==0: nproc=1; myrank=0
if ibatch==1: comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
wsm='wsm_{}'.format(p7_num) #directory of loading

#unzip all the files, and re-organize
if myrank==0 and iunzip==1: 
   #rename folders
   import tarfile
   tarfile.open('{}.tar.gz'.format(p7_version)).extractall();  os.rename(p7_version,wsm)
   for fn in glob(wsm+'/*'): os.rename(fn,'{}/{}'.format(wsm,fn.split('_')[-1]))

   #remove these files 
   snames=['EL0_4562_0003','YM4_6620_0003', #different format for these two files. other files don't have shapefile segments
   'EL0_008382377', 'EL0_008383979', 'EL0_008383983', 'EL0_008384033', 'EL0_008384037', 'EL0_008384441', 'EL0_008384489', 'EL0_008384503', 'EL0_008384507', 'EL0_008384565', 'EL0_008384841', 'EL0_008385371', 'EL0_008385389', 'EL0_008385417', 'EL0_008385507', 'EL0_008385577',
   'EL0_008385841', 'EL0_008385857', 'EL0_008385859', 'EL0_008385907', 'EL0_008385933', 'EL0_008386011', 'EL0_008386015', 'EL0_008386051', 'EL0_008386081', 'EL0_008386091', 'EL0_008386105', 'EL0_008386107', 'EL0_008386259', 'EL0_008386271', 'EL0_008386277', 'EL0_008386285',
   'EL0_008386329', 'EL0_008386345', 'EL0_008387397', 'EL0_008387399', 'EL0_008387405', 'EL0_008387411', 'EL0_008387419', 'EL0_008392684', 'EL0_008392916', 'EL0_008394030', 'EL0_008394152', 'EL0_008394156', 'EL0_008394204', 'EL0_008401537', 'EL0_008404137', 'EL0_008404225',
   'EL0_008404231', 'EL0_009411788', 'EL0_009411794', 'EL0_010232417', 'EL0_010232557', 'EL0_010232783', 'EL0_010233307', 'EL2_008397316', 'EM0_004769848', 'EM0_009408464', 'EM0_009409296', 'EM0_009409440', 'EM0_009410128', 'EM0_009411644', 'EM0_009411874', 'EM0_009411888',
   'EM0_009411890', 'EM2_009408394', 'EM2_009408436', 'EM2_009415984', 'EM3_009408326', 'EM3_009408400', 'EM4_009408514', 'EM4_009408574', 'EU0_004763958', 'EU0_004765462', 'EU0_004765480', 'EU0_004765568', 'EU0_004765618', 'EU0_004766978', 'EU0_004767218', 'EU0_004767232',
   'EU0_004767384', 'EU0_004769016', 'EU0_004769350', 'EU0_004769406', 'EU0_004770068', 'EU1_004763756', 'JA5_008608001', 'JA5_008608021', 'JA5_008608041', 'JA5_008608057', 'JA5_008608093', 'JA5_008608099', 'JA5_008608111', 'JA5_008608117', 'JA5_008608127', 'JA5_008608129',
   'JA5_008608131', 'JA5_008608143', 'JA5_008608145', 'JA5_008608149', 'JA5_008608159', 'JA5_008608167', 'JA5_008608169', 'JA5_008608173', 'JA5_008608937', 'JB0_009890708', 'JB0_009891226', 'JB0_009891320', 'JB0_009891646', 'JB0_009892970', 'JB0_010053560', 'JB0_010054690',
   'JB0_010054710', 'JB0_010054806', 'JB0_010054816', 'JB0_010054898', 'JB0_010055538', 'JB0_010056112', 'JB0_010056194', 'JB0_010057870', 'JB0_010057886', 'JB0_010064659', 'JB0_010064671', 'JB0_010064707', 'JB0_010064727', 'JB0_010064729', 'JB0_010064795', 'JB0_010064833',
   'JB0_010064983', 'JB0_010065029', 'JB0_010065041', 'JB0_010065047', 'JB0_010065061', 'JB0_010065077', 'JB0_010065095', 'JB0_010065363', 'JB0_010066405', 'JB0_010067445', 'JB0_010067545', 'JB0_010067605', 'JB0_010068667', 'JB0_010069229', 'JB0_010069261', 'PL0_004529691',
   'PL0_004532013', 'PL0_004532419', 'PL0_004532461', 'PL0_004532463', 'PL0_004532465', 'PL0_004532615', 'PL0_004532895', 'PL0_004534709', 'PL0_004535901', 'PL0_004535939', 'PL0_004536011', 'PL0_022340861', 'PL7_022340595', 'RL0_008456540', 'RL0_008456844', 'RL0_008457496',
   'RL0_008457542', 'RL0_008478716', 'RL0_008478808', 'RL0_008479370', 'RL0_008479378', 'RL0_008479402', 'RL0_008479744', 'RL0_008479748', 'RL0_008479834', 'RL0_008479838', 'RL0_008479962', 'RL0_008479988', 'RL0_008480036', 'RL0_008480110', 'RL0_008480190', 'RL0_008480216',
   'RL1_008475368', 'RL1_008475372', 'RL1_008475384', 'RL1_008475394', 'RL1_008479338', 'RL5_008475278', 'RL5_008475280', 'RL5_008475354', 'RL5_008475380', 'RL5_008475458', 'RL5_008476556', 'RL5_008477114', 'RL5_008479498', 'RL5_008479554', 'RL5_008480212', 'RL5_008480218',
   'RL5_008480328', 'RL5_008480330', 'RL5_008484332', 'WL0_008376969', 'WL0_008378195', 'WM0_011689372', 'WM0_011690166', 'WU0_011690446', 'WU0_011690552', 'XL0_008378481', 'XL0_008378487', 'XL0_011908750', 'XL0_011908756', 'XL0_011909712', 'XL3_011908556', 'XL3_011908574',
   'XL3_011908634', 'XL3_011908644', 'YL0_008457024', 'YL0_008460300', 'YL0_008460410', 'YL0_008460734', 'YL0_009888064', 'YL0_009888124', 'YL0_009888518', 'YL0_009888526', 'YL0_009888574', 'YL0_009888632', 'YL0_009896928', 'YL0_009896944', 'YL0_009896946', 'YL0_009896970',
   'YL0_009897048', 'YL0_009897052', 'YL0_009897118', 'YL0_009897328', 'YL0_009897332', 'YM0_008491696', 'YM0_008491700', 'YM0_008491728', 'YP0_008503234', 'YP0_008503310', 'YP0_008503604', 'YP0_008503632', 'YP0_008503644', 'YP0_008503810', 'YP0_008503812', 'YP0_008503850',
   'YP0_008503876', 'YP0_008503944', 'YP5_008507040', 'YP5_008507042', 'nul_004529655', 'nul_004531965', 'nul_004532023', 'nul_004532031', 'nul_004532047', 'nul_004532087', 'nul_004532103', 'nul_004532105', 'nul_004532113', 'nul_004532123', 'nul_004532251', 'nul_004532341',
   'nul_004532383', 'nul_004532501', 'nul_004532537', 'nul_004532549', 'nul_004532603', 'nul_004532831', 'nul_004533007', 'nul_004534377', 'nul_004534489', 'nul_004534525', 'nul_004535927', 'nul_004535937', 'nul_004536097', 'nul_004536293', 'nul_004763760', 'nul_004765494',
   'nul_004765514', 'nul_004765630', 'nul_004765766', 'nul_004766908', 'nul_004766960', 'nul_004766988', 'nul_004767038', 'nul_004767390', 'nul_004767466', 'nul_004769188', 'nul_004769242', 'nul_004770184', 'nul_004774658', 'nul_004774848', 'nul_008376149', 'nul_008376307',
   'nul_008376309', 'nul_008376967', 'nul_008378355', 'nul_008378471', 'nul_008378475', 'nul_008380259', 'nul_008380269', 'nul_008384457', 'nul_008384547', 'nul_008384645', 'nul_008384819', 'nul_008385981', 'nul_008385985', 'nul_008386039', 'nul_008392700', 'nul_008394140',
   'nul_008394170', 'nul_008401527', 'nul_008401955', 'nul_008401977', 'nul_008402063', 'nul_008402089', 'nul_008402143', 'nul_008402229', 'nul_008402231', 'nul_008402237', 'nul_008402243', 'nul_008402245', 'nul_008402259', 'nul_008403893', 'nul_008403961', 'nul_008404023',
   'nul_008404143', 'nul_008404169', 'nul_008404187', 'nul_008405121', 'nul_008405219', 'nul_008457676', 'nul_008458268', 'nul_008460340', 'nul_008460382', 'nul_008460648', 'nul_008460686', 'nul_008462042', 'nul_008475320', 'nul_008475428', 'nul_008475482', 'nul_008476990',
   'nul_008479426', 'nul_008479574', 'nul_008479656', 'nul_008479724', 'nul_008484212', 'nul_008484256', 'nul_008484262', 'nul_008484322', 'nul_008484344', 'nul_008491522', 'nul_008491740', 'nul_008494264', 'nul_008494280', 'nul_008503786', 'nul_008508002', 'nul_008508024',
   'nul_008607643', 'nul_009408456', 'nul_009409460', 'nul_009409916', 'nul_009410168', 'nul_009410814', 'nul_009411588', 'nul_009411858', 'nul_009415824', 'nul_009415892', 'nul_009416020', 'nul_009888566', 'nul_009888794', 'nul_009891260', 'nul_009891298', 'nul_009891580',
   'nul_009891634', 'nul_009894674', 'nul_009896994', 'nul_009897120', 'nul_009897142', 'nul_009897180', 'nul_009897194', 'nul_009897298', 'nul_009897318', 'nul_009897904', 'nul_010054824', 'nul_010055160', 'nul_010055260', 'nul_010055272', 'nul_010056018', 'nul_010057972',
   'nul_010064635', 'nul_010064773', 'nul_010064825', 'nul_010065385', 'nul_010067465', 'nul_010067583', 'nul_010068687', 'nul_010069221', 'nul_010069237', 'nul_010072859', 'nul_010072917', 'nul_010072947', 'nul_010072971', 'nul_010073053', 'nul_010231727', 'nul_010232283',
   'nul_010232443', 'nul_010232463', 'nul_010232467', 'nul_010233305', 'nul_010234563', 'nul_010234597', 'nul_010235935', 'nul_010236787', 'nul_011689864', 'nul_011689938', 'nul_011690042', 'nul_011690118', 'nul_011690548', 'nul_011690550', 'nul_011908792', 'nul_011911864',
   'nul_022338105', 'nul_022339591', 'nul_022340531', 'nul_022340559', 'nul_022340645', 'nul_022340849', 'nul_022340891', 'nul_022340893', 'nul_022341049', 'EL0_008384127', 'EL0_008384129', 'EL0_008384131', 'EL0_008384133', 'EL0_008384837', 'EL0_008384839', 'EL0_008386163',
   'EL0_008391988', 'EL0_008392528', 'EL0_008393782', 'EL0_008393784', 'EL0_008393786', 'EL0_008393788', 'EL0_008393790', 'EL0_008393794', 'EL0_008401421', 'EL0_008401431', 'EL0_008401439', 'EL0_008401529', 'EL0_008401535', 'EL0_008401545', 'EL0_008401893', 'EL0_008404183',
   'EL0_008404227', 'EL0_008404229', 'EL0_010231769', 'EL0_010231777', 'EL0_010231791', 'EL0_010231793', 'EL0_010231795', 'EL0_010232163', 'EL0_010232167', 'EL0_010233123', 'EL3_008402043', 'EM0_009409064', 'EM0_009409534', 'EM0_009410078', 'EM0_009411904', 'EU0_004766492',
   'JB0_009890722', 'JB0_009890740', 'JB0_009890752', 'JB0_009890792', 'JB0_009890812', 'JB0_009890826', 'JB0_009891578', 'JB0_010053650', 'JB0_010053800', 'JB0_010054088', 'JB0_010054202', 'JB0_010054206', 'JB0_010054210', 'JB0_010054284', 'JB0_010055102', 'JB0_010055540',
   'JB0_010063989', 'JB0_010064005', 'JB0_010065045', 'JB0_010066395', 'JB0_010067115', 'JB0_010072977', 'PL0_004534123', 'RL0_008456840', 'RL0_008478792', 'RL5_008477796', 'WL0_008377427', 'XL0_008378483', 'YL0_008459858', 'YL0_008459862', 'YL0_008460470', 'YL0_008460472',
   'YL0_008460644', 'YL0_008460646', 'YL0_008460650', 'YL0_008460652', 'YL0_008460656', 'YL0_008460658', 'YL0_008460664', 'YL0_008460736', 'YL0_009896514', 'YL0_009897154', 'YL0_009897286', 'YL0_009897288', 'YL0_009897290', 'YL0_009897292', 'YL0_009897294', 'YL0_009897300',
   'YP0_008503214', 'YP0_008503238', 'nul_004529663', 'nul_008378191', 'nul_008378485', 'nul_008401543', 'nul_008402271', 'nul_008460642', 'nul_008460654', 'nul_008460660', 'nul_008460662', 'nul_008460688', 'nul_008460690', 'nul_009888790', 'nul_009888792', 'nul_009892900',
   'nul_009897122', 'nul_009897282', 'nul_009897284', 'nul_010064905', 'nul_011693992', 'nul_011908654', 'nul_022340533', 'nul_022340543']
   [os.remove(fn) for fn in glob(wsm+'/*/*') if fn.replace('/','.').split('.')[-2] in snames]
if ibatch==1: comm.Barrier()

#read terminal and tidal data
fns=glob(wsm+'/*/*_*.*'); svar=[]; sname=[]; stype=[]; data=[] 
for i,fn in enumerate(fns):
    if i%nproc!=myrank: continue
    bn=os.path.basename(fn); svari=bn.split('.')[1]; snamei=bn.split('.')[0]
    print(i,fn,myrank); sys.stdout.flush() 
    fid=open(fn,'r'); data0=array([i.strip().split(',') for i in fid.readlines()[1:]]); datai=data0[:,3].astype('float32'); fid.close()
    if len(data0)!=13149 or data0.shape[1]!=4: sys.exit('wrong:',fn, data0.shape); sys.stdout.flush()
    svar.append(svari); sname.append(snamei); stype.append(fn.split('/')[1]); data.append(datai)
C=zdata(); C.var=array(svar); C.sname=array(sname); C.stype=array(stype); C.data=array(data); C.save('data_{}'.format(myrank))

#combine data
if ibatch==1: comm.Barrier()
if myrank==0:
   #combine data 
   svar=[]; sname=[]; stype=[]; data=[]
   for myrank in arange(nproc):
       fn='data_{}.npz'.format(myrank)
       print('reading {} loading '.format(fn)); sys.stdout.flush() 
       C=read(fn); svar.extend(C.var); sname.extend(C.sname); stype.extend(C.stype)
       data.extend(C.data if C.data.ndim==2 else C.data[None,:]); os.remove(fn)
   svar=array(svar); sname=array(sname); stype=array(stype); data=array(data)  

   C=zdata(); C.var=unique(svar); cdict=C.__dict__
   C.sname,fps=unique(sname,return_index=True); C.stype=stype[fps] 
   for svari in C.var:
       fp=svar==svari; datai=data[fp]; snamei=sname[fp]
       sind=argsort(snamei); datai=datai[sind]; cdict[svari]=datai 
       if not array_equal(snamei[sind],C.sname): sys.exit('names are different')
   C.time=datenum(1985,1,1)+arange(C.flow.shape[1])
   C.unit={'chla':'ug/l', 'clay':'kg/d', 'doxx':'mg/l', 'flow':'cms',  'nh4x':'kg/d', 'no3x':'kg/d',
           'orgn':'kg/d', 'orgp':'kg/d', 'phyt':'kg/d', 'pipx':'kg/d', 'po4x':'kg/d', 'sand':'kg/d',
           'silt':'kg/d', 'tocx':'kg/d', 'totn':'kg/d', 'totp':'kg/d', 'tssx':'kg/d', 'wtmp':'dg C'}

   #add attributes of watershed
   S=read('shapefile/wsm.shp'); atn=list(S.attname)
   sname=S.attvalue[atn.index('SegmentNam')]; fps=sindex(C.sname,sname,1); C.xy=S.xy[fps]
   C.river=S.attvalue[atn.index('river')][fps]; C.cbseg=S.attvalue[atn.index('cbseg')][fps] 
   print('finished reading {} loading, {} '.format(p7_version,time.time()-t0)); sys.stdout.flush() 

   #-------------------------------------------------------------------
   #compute allocation to MBM sides
   #-------------------------------------------------------------------
   #input
   distc=100
   slink={'PAXMH4':'PAXMH','WBRTF':'PAXTF'}
   bdir='../../database/region/'
   rnames=['SL9_2720_0001','PM7_4820_0001','JL7_7070_0001','RU5_6030_0001','JA5_7480_0001','YP4_6750_0001','YM4_6620_0001','XU3_4650_0001','EM2_3980_0001'] #RIM river

   #read schism grid side info
   gd=grd('../../grid/v8f'); gd.compute_all()
   sid0=pindex(gd.isdel[:,1],-1) #side indices: sid is the boundary side index that received loads
   sid=sid0[read(bdir+'estuary.reg').inside(gd.sxy[sid0],prj=['epsg:26918','epsg:4326'])]
   xs,ys=gd.sxy[sid].T; sxy=xs+1j*ys; ind=gd.isidenode[sid]; bxy=c_[gd.x[ind[:,0]],gd.y[ind[:,0]],gd.x[ind[:,1]],gd.y[ind[:,1]]]

   #get the nearby side index for each watershed segments
   print('computing WSM loading allocation to MBM sides'); sys.stdout.flush() 
   xc,yc=array([i[~isnan(i[:,0])].mean(axis=0) for i in C.xy]).T; cxy=xc+1j*yc
   wsid=[]; dms=[]
   for i,[sname,xy] in enumerate(zip(C.sname,C.xy)):
       if i%1000==0: print(i)
       for di in [10e3,15e3,20e3,30e3]: #search nearby sides
           pid=pindex(abs(sxy-cxy[i])<di)
           if len(pid)!=0: break
       if sname=='YM4_6620_0001': pid=read(bdir+'river_head_mattaponi.reg').inside(c_[xs,ys]) #Mattaponi 
       if sname=='YP4_6750_0001': pid=read(bdir+'river_head_pamunkey.reg').inside(c_[xs,ys])  #Pamunkey
       dist=mdist(bxy[pid],xy,4); dm=dist.min() if sname not in ['YM4_6620_0001','YP4_6750_0001'] else dist.max()
       wsid.append(sid[pid[dist==dm]]); dms.append(dm)
   wsid=array(wsid,'O'); dms=array(dms)

   #replace side index for segments
   sinds=pindex(dms>distc);  sindp=[]
   for m,sind in enumerate(sinds):
       sname=C.sname[sind]; cbseg=C.cbseg[sind]; sindn=pindex((C.cbseg==cbseg)*(dms<=distc)) #sindn refer to same cbseg watersheds
       if len(sindn)==0: sindn=pindex((C.cbseg==slink[C.cbseg[sind]])*(dms<=distc))
       ds=abs(cxy[sindn]-cxy[sind]); wid=sindn[ds.argmin()]; sindp.append(wid);
       if sname in rnames: #check RIM one by one
          cid=intersect1d(wsid[sind],wsid[wid]); wsid[sind]=cid if len(cid)!=0 else wsid[sind]
       else:
          wsid[sind]=wsid[wid]
   sindp=array(sindp)

   #orangize information for watershed loading partition
   sinds=[]; [sinds.extend(i) for i in wsid]; sinds=list(unique(sinds)); sdist=gd.distj[sinds]
   C.sxy=gd.sxy[sinds]; C.sid=array([array([sinds.index(k) for k in i]) for i in wsid],'O')
   C.srat=array([sdist[i]/sdist[i].sum() for i in C.sid],'O')
   C.note='sxy: loading points; sid: mapping for each watershed to each loading point; srat: mapping ratio'

   #-------------------------------------------------------------------
   #read shoreline erosion
   #-------------------------------------------------------------------
   print('reading shoreline erosion data'); sys.stdout.flush() 
   S=read('sho_v1/shoreef.nc',2); S.sxy=c_[S.x,S.y]; S.time=datenum(1991,1,1)+arange(len(S.time))

   #save results
   W=zdata(); W.wsm=C; W.sho=S; W.save('load_p7_{}'.format(p7_num))

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
