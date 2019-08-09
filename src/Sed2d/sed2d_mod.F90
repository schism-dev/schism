!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

MODULE sed2d_mod
!--------------------------------------------------------------------
! This module contains the variables shared by sed2d routines      
!                                                                 
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)    
! Date:   06/12/2012      
!
! History:                                          
! 01/2013 - G.Dodet:  Sorting and comments
! 02/2013 - G.Dodet:  Removed ibslope and added idrag and irough
! 03/2013 - G.Dodet:  Added Cdsed and d50/d90 arrays
! 04/2013 - G.Dodet:  Added volume control area (for extrema filter),
!                     iskip,nskip, imeth, qtot_s, bed_delta_e, qb_s,
!                     qb_e, qs_s, qs_e, dpdxy_e, dpdxy_s.
! 05/2013 - G.Dodet:  Added slope_cr
! 06/2013 - G.Dodet:  - Added bedforms-associated z0 for outputs;
!                     - Added qfilter, ufilter, nskip. 
!                     - Added qramp, qsum_e,qdt_e and qav;
!                     - Added cflsed
! 07/2013 - T.Guerin: Added parameters for CL11 formula: hs0, tp0
! 09/2013 - T.Guerin: Added iasym
! 03/2014 - T.Guerin: Added variables for multi-class multi-layer
!                     approach: nb_class, nb_layer, h_inf, h_top,
!                     d50moy, d90moy, F_class, h2
! 04/2014 - T.Guerin: - Added h_lim_max and h_lim_min
!                     - Removed hs0 and tp0
! 10/2016 - T.Guerin: Modifications related to the merge of single-
!                     class and multi-class routines
! 04/2017 - T.Guerin: - Added dryslope and wetslope
!                     - Added imnp (i.e. morphological ramp variable)
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: rkind

  IMPLICIT NONE
  SAVE

  INTEGER :: dtsed2d,   & !Morphodynamic time step = dtsed2d x dt
             iasym,     & !Flag for asymetric treatment of wave orbital velocity in CL (2011)
             idrag,     & !Flag for drag coefficient formula
             ifilt,     & !Flag for filter during simulation
             imeth,     & !Flag for numerical method
             imorpho,   & !Flag for morphodynamic simulation
             ipre_filt, & !Flag for filter type at pre-processing
             ipre_flag, & !Flag for filter at pre-processing
             irough,    & !Flag for roughness length formula
             iskip,     & !Nb of iterations before starting morpho.
             islope,    & !Flag for slope effect in SVR (1997)
             itrans,    & !Flag for transport formula
             nb_class,  & !Nb of sediment classes
             nb_layer,  & !Nb of sediment layers
             qfilter,   & !Flag for total transport filtering
             qramp,     & !Ramp period (s) for total transport
             nskip,     & !Nb of iterations between filters
             ufilter,   & !Flag for velocity filtering (Shapiro)
             ised_dump    !dredging/dumping option

  REAL(rkind) :: diffac,   & !Diffusion factor (-)
                 dryslope, & !Maximum dry slope (-) if slope filter is active
                 h_inf,    & !thickness of layers 3...N, and init. thickness of layer 2 (m)
                 h_lim_max,& !Maximum thickness for sediment layer #2 (m)
                 h_lim_min,& !Minimum thickness for sediment layer #2 (m)
                 h_top,    & !Top-layer thickness (m)
                 h0_sed,   & !Minimum depth for transport (m)
                 poro,     & !Porosity (-)
                 s,        & !Ratio rhosed/rho0 (-)
                 transfac, & !Transport correction factor (-)
                 wetslope, & !Maximum wet slope (-) if slope filter is active
                 wvisco,   & !Water kinematic viscosity (m2/s)
                 morfac      !morphological acceleration factor (>0)

  REAL(rkind), ALLOCATABLE :: bc_val(:),   & !Bnd nodes num. flag (-)
                              bed_delta(:),& !Bed delta (m) at nodes; bed_delta(npa)
                              bed_del_e(:),& !... same at elements; bed_del_e(nea)
                              Cdsed(:),    & !Drag coefficient (-); Cdsed(npa)
                              cflsed(:),   & !Courant number in sed2d; cflsed(npa)
                              d_class(:),  & !Grain size for each sediment class (m); d_class(nb_class)
                              d50(:,:),    & !d50 (m) at each node and for each sediment layer; d50(npa,nb_layer)
                              d90(:,:),    & !d90 (m); d90(npa,nb_layer)
                              dpdxy(:,:),  & !Bed slope (m/m),nodes; dpdxy(npa,2)
                              dpdxy_e(:,:),& !... same at elements; dpdxy_e(nea,2)
                              dpdxy_s(:,:),& !... same at sides; dpdxy_s(nsa,2)
                              F_class(:,:,:),& !Sediment classes fraction (-) (>=0); F_class(npa,nb_class,nb_layer) (only if nb_class>1). Layer 1 is at the top.
                              h2(:),       & !layer #2 thickness (m); h2(npa). Only this layer has a variable thickness
                              imnp(:),     & !morphological ramp value (-); imnp(npa)
                              mcoefd(:,:), & !Matrix coef for JCG; mcoefd(0:mnei_p,np)
                              qav(:,:),    & !Averaged qtot (kg/m/s); qav(npa,2)
                              qb(:,:),     & !Bed load (m3/m/s),nodes; qb(npa,2)
                              qb_e(:,:),   & !... same at elements; qb_e(nea,2)
                              qb_s(:,:),   & !... same at sides; qb_s(nsa,2)
                              qdt_e(:,:),  & !Integ. qtot (m3/m); qdt_e(nea,2)
                              qs(:,:),     & !Suspended load (m3/m/s); qs(npa,2)
                              qs_e(:,:),   & !... same at elements; qs_e(nea,2)
                              qs_s(:,:),   & !... same at sides; qs_s(nsa,2)
                              qtot(:,:),   & !Total load (m3/m/s); qtot(npa,2)
                              qtot_e(:,:), & !...same at elements; qtot_e(nea,2)
                              qtot_s(:,:), & !...same at sides; qtot_s(nsa,2)
                              slope_cr(:), & !Max. slopes for filter; slope_cr(npa)
                              u2_tmp(:),   & !u2_tmp(npa)
                              v2_tmp(:),   & !v2_tmp(npa)
                              vc_area(:),  & !Volume control area (m2); vc_area(npa)
                              z0_e(:),     & !Roughness length z0 (m); z0_e(nea)
                              z0mr_e(:),   & !Mega-ripples z0 (m); z0mr_e(nea)
                              z0cr_e(:),   & !Current ripples z0 (m); z0cr_e(nea)
                              z0sw_e(:),   & !Sand- wave z0 (m); z0sw_e(nea)
                              z0wr_e(:)     !Wave ripples z0 (m); z0wr_e(nea)


  LOGICAL, ALLOCATABLE :: bc_flag(:) !Bnd nodes log. flag

  INTEGER, PARAMETER :: idsed2d = 26 !File ID for "sed2d.out"

END MODULE sed2d_mod
