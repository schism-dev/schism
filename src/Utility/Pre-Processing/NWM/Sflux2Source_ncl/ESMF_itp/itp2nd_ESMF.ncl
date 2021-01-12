; This script is used to interpolate sflux prc into node
;
; Input: sflux prc ncfiles, modify path in "dir"
;        itp_X,itp_Y,ele.dat from get_itpXYA.f90
;
; Output: vsource_np.bin,dim.dat
;
; If sflux prc nc files are global, set Opt@SrcRegional=False
; If hgrid.gr3 is global, set Opt@DstRegional=False
;


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
itp_x=asciiread("itp_X" , -1, "float")
itp_y=asciiread("itp_Y" , -1, "float")
nv=asciiread("ele.dat",-1,"integer") ; only accept triangle
elesize=dimsizes(nv)/3
ele=onedtond(nv,(/elesize,3/))


dir = "./sflux_nc/"
files= systemfunc("ls "+ dir + "sflux_prc_1.*nc")
fs= addfiles(files,"r")
prc=fs[:]->prate ;(:,:,:)
prc_x=fs[0]->lon
prc_y=fs[0]->lat

;---Data file containing source grid

;---Get variable to regrid
    var      = prc                          ;;
    src_lat  = prc_y                        ;;
    src_lon  = prc_x                        ;;


;---Set up regridding options
    Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod      = "bilinear"        ;;

    Opt@WgtFileName       = "curv_to_unstruct.nc"

    Opt@SrcGridLat        = prc_y           ; source grid
    Opt@SrcGridLon        = prc_x
    Opt@SrcRegional       = True             ;;

    Opt@DstGridLat        = itp_y             ; destination grid
    Opt@DstGridLon        = itp_x
    Opt@DstRegional       = True             ;;
    Opt@DstGridType       = "unstructured"
    Opt@DstTriangularMesh = ele-1              ; (ele table, remember to -1)
    Opt@DstMask2D = where(.not.ismissing(itp_y).and.\
                          .not.ismissing(itp_x),1,0) ; Necessary if lat/lon

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding and output weight file

    printVarSummary(var_regrid)
    system("rm -f ./vsource_np.bin")
    ndim=dimsizes(var_regrid)
    do l=1,ndim(0)
;      printMinMax(var_regrid(l-1,:),False)
       fbinrecwrite("vsource_np.bin",-1,var_regrid(l-1,:))
    end do
    asciiwrite("dim.dat",ndim)

end
