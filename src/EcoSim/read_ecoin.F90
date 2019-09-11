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

!!==============================================================================!
!! February, 2009                                                               !
!!======================================================Marta Rodrigues=========!
!===============================================================================!
!===============================================================================!
! Routine to read in ecosim.in  - adapted from read_param.F90                   ! 
!===============================================================================!
!===============================================================================!
subroutine read_ecoin

     use bio_param
     use biology
     use schism_glbl, only: rkind,itr_met
     use schism_msgp, only : myrank,parallel_abort  

     implicit none
     
     logical :: tmp3
     character(len=20) :: filename
     character(len=1) :: it_char 
     integer*4 :: tmp1,i,j
     real(rkind) :: tmp2
     
! Initialize filename

     filename(1:20)='                    '     
     
!     call read_inputs('BioIter',1,tmp1,tmp2,tmp3)
!     BioIter=tmp1
     
     filename(1:20)='                    '

     call read_inputs('RtUVR_flag',0,tmp1,tmp2,tmp3)
     RtUVR_flag=tmp3
     
     filename(1:20)='                    '
         
     call read_inputs('Regen_flag',0,tmp1,tmp2,tmp3)
     Regen_flag=tmp3
     
     filename(1:20)='                    '  

     call read_inputs('NIT_flag',1,tmp1,tmp2,tmp3)
     NIT_flag=tmp1
     
     filename(1:20)='                    '
     
     call read_inputs('DENIT_flag',1,tmp1,tmp2,tmp3)
     DENIT_flag=tmp1
     
     filename(1:20)='                    '
     
     call read_inputs('REAER_flag',1,tmp1,tmp2,tmp3)
     REAER_flag=tmp1
     
              
     
!-----------------------------------------------------------------------------
! Phytoplankton group parameters.
!-----------------------------------------------------------------------------
!
     filename(1:20)='                    '
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsNO3'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsNO3(i)=tmp2
     END DO
     
     filename(1:20)='                    '

     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsNH4'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsNH4(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsSiO'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsSiO(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsPO4'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsPO4(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:4)='HsFe'
       filename(5:5)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsFe(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='GtALG_max'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       GtALG_max(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='PhyTbase'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       PhyTbase(i)=tmp2
     END DO
     
     filename(1:20)='                    '

     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:7)='PhyTfac'
       filename(8:8)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       PhyTfac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:4)='BET_'
       filename(5:5)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       BET_(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='maxC2nALG'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       maxC2nALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='minC2nALG'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       minC2nALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:12)='C2nALGminABS'
       filename(13:13)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2nALGminABS(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:10)='maxC2SiALG'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       maxC2SiALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:10)='minC2SiALG'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       minC2SiALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:13)='C2SiALGminABS'
       filename(14:14)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2SiALGminABS(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='maxC2pALG'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       maxC2pALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='minC2pALG'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       minC2pALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:12)='C2pALGminABS'
       filename(13:13)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2pALGminABS(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:10)='maxC2FeALG'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       maxC2FeALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:10)='minC2FeALG'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       minC2FeALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:13)='C2FeALGminABS'
       filename(14:14)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2FeALGminABS(i)=tmp2
     END DO 
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='qu_yld'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       qu_yld(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:7)='E0_comp'
       filename(8:8)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       E0_comp(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='E0_inhib'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       E0_inhib(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='inhib_fac'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       inhib_fac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='C2CHL_max'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2CHL_max(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxC2Cl'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxC2Cl(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_C2Cl'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_C2Cl(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxC2Cn'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxC2Cn(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_C2Cn'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_C2Cn(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='mxPacEff'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxPacEff(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='b_PacEff'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_PacEff(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxChlB'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxChlB(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_ChlB'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_ChlB(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxChlC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxChlC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_ChlC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_ChlC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='mxPSC'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxPSC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='b_PSC'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_PSC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='mxPPC'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxPPC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='b_PPC'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_PPC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxLPUb'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxLPUb(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_LPUb'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_LPUb(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='mxHPUb'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       mxHPUb(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='b_HPUb'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       b_HPUb(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='FecDOC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       FecDOC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO j=1,Nfec
       write(it_char,'(i1)')j
       filename(8:8)=it_char
       DO i=1,Nphy
         write(it_char,'(i1)')i
         filename(1:6)='FecPEL'
         filename(7:7)=it_char
         call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
         FecPEL(i,j)=tmp2
       END DO 	 
     END DO
     
     filename(1:20)='                    '
             
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:6)='FecCYC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       FecCYC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='ExALG'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       ExALG(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:2)='WS'
       filename(3:3)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       WS(i)=tmp2
       if(WS(i)>0.d0.and.itr_met/=3) call parallel_abort('ECOSIM: Phy sinking must use itr_met=3')
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsGRZ'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsGRZ(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='basalPhy'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       basalPhy(i)=tmp2
     END DO
    
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:4)='QPhy'
       filename(5:5)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       QPhy(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:7)='gamaPhy'
       filename(8:8)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       gamaPhy(i)=tmp2
     END DO  
     
     filename(1:20)='                    ' 
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='MinRefuge'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       MinRefuge(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='RefugeDep'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RefugeDep(i)=tmp2
     END DO
     
     filename(1:20)='                    '
    
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:8)='Norm_Vol'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       Norm_Vol(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='Norm_Surf'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       Norm_Surf(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsDOP'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsDOP(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:10)='C2pALKPHOS'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2pALKPHOS(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:5)='HsDON'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsDON(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nphy
       write(it_char,'(i1)')i
       filename(1:9)='C2nNupDON'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       C2nNupDON(i)=tmp2
     END DO
      
!-----------------------------------------------------------------------------
! Bacteria group parameters.
!-----------------------------------------------------------------------------
!
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:8)='HsDOC_ba'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsDOC_ba(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:9)='GtBAC_max'
       filename(10:10)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       GtBAC_max(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:8)='BacTbase'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       BacTbase(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:7)='BacTfac'
       filename(8:8)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       BacTfac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     call read_inputs('C2nBAC',2,tmp1,tmp2,tmp3)
     C2nBAC=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('C2pBAC',2,tmp1,tmp2,tmp3)
     C2pBAC=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('C2FeBAC',2,tmp1,tmp2,tmp3)
     C2FeBAC=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('BacDOC',2,tmp1,tmp2,tmp3)
     BacDOC=tmp2
     
     call read_inputs('BacPEL',2,tmp1,tmp2,tmp3)
     BacPEL=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('BacCYC',2,tmp1,tmp2,tmp3)
     BacCYC=tmp2
     
     call read_inputs('ExBAC_c',2,tmp1,tmp2,tmp3)
     ExBAC_c=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('ExBacC2N',2,tmp1,tmp2,tmp3)
     ExBacC2N=tmp2
     
     call read_inputs('Bac_Ceff',2,tmp1,tmp2,tmp3)
     Bac_Ceff=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('RtNIT',2,tmp1,tmp2,tmp3)
     RtNIT=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('HsNIT',2,tmp1,tmp2,tmp3)
     HsNIT=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('QN',2,tmp1,tmp2,tmp3)
     QN=tmp2
     
     filename(1:20)='                    '

     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:8)='basalBac'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       basalBac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:4)='QBac'
       filename(5:5)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       QBac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:5)='GEE0C'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       GEE0C(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nbac
       write(it_char,'(i1)')i
       filename(1:6)='HsBacO'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsBacO(i)=tmp2
     END DO   
     
!-----------------------------------------------------------------------------
! DOM group parameters.
!-----------------------------------------------------------------------------
!
     filename(1:20)='                    '
     
     DO i=1,Ndom
       write(it_char,'(i1)')i
       filename(1:10)='cDOCfrac_c'
       filename(11:11)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       cDOCfrac_c(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     call read_inputs('RtUVR_DIC',2,tmp1,tmp2,tmp3)
     RtUVR_DIC=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('RtUVR_DOC',2,tmp1,tmp2,tmp3)
     RtUVR_DOC=tmp2
      
!-----------------------------------------------------------------------------
! Fecal and detritus group parameters.
!-----------------------------------------------------------------------------
!
     filename(1:20)='                    '     

     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:2)='WF'
       filename(3:3)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       WF(i)=tmp2
       if(WF(i)>0.d0.and.itr_met/=3) call parallel_abort('ECOSIM: Fec sinking must use itr_met=3')
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:8)='RegTbase'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegTbase(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:7)='RegTfac'
       filename(8:8)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegTfac(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:5)='RegCR'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegCR(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:5)='RegNR'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegNR(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:5)='RegSR'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegSR(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:5)='RegPR'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegPR(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nfec
       write(it_char,'(i1)')i
       filename(1:5)='RegFR'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       RegFR(i)=tmp2
     END DO

!-----------------------------------------------------------------------------
! Zooplankton group parameters 
!-----------------------------------------------------------------------------
!          
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:6)='zoo_sp'
       filename(7:7)=it_char
       call read_inputs(trim(filename),1,tmp1,tmp2,tmp3)
       zoo_sp(i)=tmp1
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:6)='ZooDOC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       ZooDOC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO j=1,Nfec
       write(it_char,'(i1)')j
       filename(8:8)=it_char
       DO i=1,Nzoo
         write(it_char,'(i1)')i
         filename(1:6)='ZooPEL'
         filename(7:7)=it_char
         call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
         ZooPEL(i,j)=tmp2
       END DO	 
     END DO
     
     filename(1:20)='                    ' 
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:6)='ZooCYC'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       ZooCYC(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO j=1,Nphy
       write(it_char,'(i1)')j
       filename(10:10)=it_char
       DO i=1,Nzoo
         write(it_char,'(i1)')i
         filename(1:8)='DeltaZoo'
         filename(9:9)=it_char
         call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
         DeltaZoo(i,j)=tmp2
       END DO
     END DO 
     
     filename(1:20)='                    ' 
     
     DO j=1,Nphy
       write(it_char,'(i1)')j
       filename(8:8)=it_char
       DO i=1,Nzoo
         write(it_char,'(i1)')i
         filename(1:6)='EfcCap'
         filename(7:7)=it_char
         call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
         EfcCap(i,j)=tmp2
       END DO
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:5)='HsZoo'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       HsZoo(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:6)='EfcPrd'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       EfcPrd(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:5)='ExZoo'
       filename(6:6)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       ExZoo(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:2)='GZ'
       filename(3:3)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       GZ(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:8)='basalZoo'
       filename(9:9)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       basalZoo(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:4)='QZoo'
       filename(5:5)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       QZoo(i)=tmp2
     END DO
     
     filename(1:20)='                    '
     
     DO i=1,Nzoo
       write(it_char,'(i1)')i
       filename(1:6)='etaZoo'
       filename(7:7)=it_char
       call read_inputs(trim(filename),2,tmp1,tmp2,tmp3)
       etaZoo(i)=tmp2
     END DO         

!-----------------------------------------------------------------------------
! Oxygen parameters 
!-----------------------------------------------------------------------------
!
     filename(1:20)='                    '     

     call read_inputs('omegaO2C',2,tmp1,tmp2,tmp3)
     omegaO2C=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('omegaO2N',2,tmp1,tmp2,tmp3)
     omegaO2N=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('omegaS2O',2,tmp1,tmp2,tmp3)
     omegaS2O=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('QWind',2,tmp1,tmp2,tmp3)
     QWind=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('psiWind',2,tmp1,tmp2,tmp3)
     psiWind=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('alfaWind',2,tmp1,tmp2,tmp3)
     alfaWind=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('RtDenit',2,tmp1,tmp2,tmp3)
     RtDenit=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('MDenit',2,tmp1,tmp2,tmp3)
     MDenit=tmp2
     
     filename(1:20)='                    '

     call read_inputs('omegaO2NDenit',2,tmp1,tmp2,tmp3)
     omegaO2NDenit=tmp2

     filename(1:20)='                    '
     
     call read_inputs('reox_COD',2,tmp1,tmp2,tmp3)
     reox_COD=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('QCOD',2,tmp1,tmp2,tmp3)
     QCOD=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('HsCOD',2,tmp1,tmp2,tmp3)
     HsCOD=tmp2
    
     filename(1:20)='                    '
     
     call read_inputs('pCO2a',2,tmp1,tmp2,tmp3)
     pCO2a=tmp2
     
     filename(1:20)='                    '
     
     call read_inputs('pH',2,tmp1,tmp2,tmp3)
     pH=tmp2        
end subroutine read_ecoin

!--------------------------------------------------------------------------------
subroutine read_inputs (varname,vartype,ivarvalue,varvalue1,varvalue2)
  ! Get a parameter from ecosim.in
  ! Inputs:
  !        varname: parameter name (string no longer than 90)
  !        vartype: parameter value type (0: logical; 1: integer; 2: float)
  ! Outputs:
  !        ivarvalue: integer output;
  !        varvalue1: float output;
  !        varvalue2: logical output.
  ! Format rules for ecosim.in:
  ! (1) Lines beginning with "!" are comments; blank lines are ignored;
  ! (2) one line for each parameter in the format: keywords= value;
  !     keywords are case sensitive; spaces allowed between keywords and "=" and value;
  !     comments starting with "!"  allowed after value;
  ! (3) value is an integer, double, or 2-char string; for double, any of the format is acceptable:
  !     40 40. 4.e1
  !     Use of decimal point in integers is OK but discouraged.
  use schism_glbl, only : rkind,errmsg,in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : parallel_abort,myrank
  use bio_param
  use biology
  implicit real(rkind)(a-h,o-z), integer(i-n)

  character(*),intent(in) :: varname
  integer*4,intent(in) :: vartype
  integer, intent(out) :: ivarvalue
  real(rkind), intent(out) :: varvalue1
  logical, intent(out) :: varvalue2

  character(len=90) :: line_str,str_tmp,str_tmp2

  str_tmp2=adjustl(varname)
  lstr_tmp2=len_trim(str_tmp2)
!  print*, varname ,str_tmp2(1:lstr_tmp2)

  ! Scan param.in
  open(5,file=in_dir(1:len_in_dir)//'ecosim.in',status='old')
  rewind(5)
  line=0
  do
    line=line+1
    read(5,'(a)',end=99)line_str
    line_str=adjustl(line_str) !place blanks at end
    len_str=len_trim(line_str)
    if(len_str==0.or.line_str(1:1)=='!') cycle

    loc=index(line_str,'=')
    loc2=index(line_str,'!')
    if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: exclam. before =')
    
    str_tmp=''
    str_tmp(1:loc-1)=line_str(1:loc-1) !keyword
    str_tmp=trim(str_tmp)
    lstr_tmp=len_trim(str_tmp)
    
    if(str_tmp(1:lstr_tmp)==str_tmp2(1:lstr_tmp2)) then
       if(loc2/=0) then
         str_tmp2=line_str(loc+1:loc2-1)
       else
         str_tmp2=line_str(loc+1:len_str)
       endif
       str_tmp2=adjustl(str_tmp2)
       str_tmp2=trim(str_tmp2)
       if(vartype==0) then !string
         read(str_tmp2,*) varvalue2
#ifdef DEBUG
         if(myrank==0) write(99,*)varname,' = ',varvalue2
#endif
       else if(vartype==1) then !integer
	  read(str_tmp2,*) ivarvalue
  
#ifdef DEBUG
         if(myrank==0) write(99,*)varname,' = ',ivarvalue
#endif
       else if(vartype==2) then !float
          read(str_tmp2,*)varvalue1
	 
#ifdef DEBUG
         if(myrank==0) write(99,*)varname,' = ',real(varvalue1)
#endif
       else
         write(errmsg,*)'read_param: unknown type:',vartype
         call parallel_abort(errmsg)
       endif
       exit
    endif
  enddo !scane ecosim.in
  
!   print*, 'Found it on line: ',line
   close(5)
   return

99  close(5)
   write(errmsg,*)'Failed to find parameter:',varname
   call parallel_abort(errmsg)

end subroutine read_inputs
!---------------------------------------------------------------------------

