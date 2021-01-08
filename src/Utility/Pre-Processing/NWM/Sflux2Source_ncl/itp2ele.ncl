; This script is used to interpolate sflux prc into element
;
; Input: sflux prc ncfiles, modify path in "dir"  
;        itp_X,itp_Y,itp_A from get_itpXYA.f90
;
; Output: vsource.bin
;
; Itp option: itpm=0, bilinear ; =1, curvilinear (much slower)
;             Only use itpm=0 if you are sure lon/lat are actually 1d array.
;

system("rm -f ./vsource.bin")
itp_x=asciiread("itp_X" , -1, "float")
itp_y=asciiread("itp_Y" , -1, "float")
itp_area=asciiread("itp_A" , -1, "float")

dir = "./sflux_nc/"
itpm = 0 ; 0: bilinear, 1: curvilinear (slower)
files= systemfunc("ls "+ dir + "sflux_prc_1.*nc")
fs= addfiles(files,"r")
prc=fs[:]->prate
if (itpm.eq.0) then
   prc_x=fs[0]->lon(0,:)
   prc_y=fs[0]->lat(:,0)
else
   prc_x=fs[0]->lon
   prc_y=fs[0]->lat
end if
nt0=fs[:]->time

np=dimsizes(itp_x)
nt=dimsizes(nt0)
vsource=new((/np/),"float")
outd=new((/2/),"integer")
outd(0)=np
outd(1)=nt

ic=0
do l=1,nt
   if (itpm.eq.0) then
       vsource(:)=linint2_points(prc_x,prc_y,prc(l-1,:,:),False,itp_x,itp_y,0)*itp_area(:)/1000
   else
       vsource(:)=rcm2points(prc_y,prc_x,prc(l-1,:,:),itp_y,itp_x,2)*itp_area(:)/1000
   end if
   printMinMax(vsource(:),False)
   ic=ic+1
   fbinrecwrite("vsource.bin",-1,vsource)   
end do
asciiwrite("dim.dat",outd)

