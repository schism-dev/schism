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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       
!           ______  _____  __   __  __  ______   __________
!          /_ ___/ / ___/ / /  / / I I /_ ___/  / ___*__  /
!          ____ \ / /    / /__/ / I I  ____ \  / /  *  / /
!         ____/ // /___ / /HH/ / I I  ____/ / / /  *  / /
!         /____/ \____//_/  /_/ I_I   /____/ /_/  *  /_/
!
!
!         SCHISM (Semi-implicit Cross-scale Hydroscience Integrated System Model)                        
!         A Three-Dimensional Model on Unstructured Grids for Hydroscience
!         This is a derivative work from the original SELFE model of STC-CMOP,
!         Oregon Health & Science University, Beaverton, Oregon 97006, USA         
!         A copy of the original SELFE (v3.1dc as of Dec 13, 2014) can be found in NOTICE.
!                                                                                       
!         SCHISM developers:
!                    Lead: Joseph Zhang (OHSU & VIMS)
!                    Air-sea exchange: Mike Zulauf (OHSU)
!                    Ecology: Marta Rodrigues,Anabela Oliveira (LNEC)
!                    Sediment: Guillaume Dodet, Florian Ganthy, Ligia Pinto,Andre Fortunato (LNEC)
!                    Oil Spill: Alberto Azvedo/Anabela Oliveira (LNEC)
!                    Waves: Aron Roland (Zanke & Partner),Ulrich Zanke (Zanke & Partner) 
!                    Water quality: Harry Wang (VIMS)
!                    Hydraulic structures: Eli Ateljevich (CAL-DWR)
!                    Ecosystem: Fei Chai (U. Maine)
!                                                                                       
!       The heat exchange module makes use of the bulk aerodynamic surface flux         
!       algorithm introduced by Zeng et al (1998), and the polynomial fits to           
!       saturation vapor pressure of Flatau et al (1992):                                
!       Zeng, X., M. Zhao, and R. E. Dickinson, 1998:  Intercomparison of bulk          
!       aerodynamic algorithms for the computation of sea surface fluxes using          
!       TOGA COARE and TAO data.  J. Clim., 11, 2628-2644.                              
!       Flatau, P. J., R. L. Walko and W. R. Cotton, 1992:  Polynomial fits to          
!       saturation vapor pressure.  J. Appl. Meteor., 31, 1507-1513.                    
!                                                                                       
!       Attenuation of solar radiation (and solar heating) within the water column      
!       is based upon the expression given by Paulson and Simpson (1977), for the       
!       water types defined by Jerlov (1968):                                           
!       Jerlov, N. G., Optical Oceanography, Elsevier, 1968.                            
!       Paulson, C. A., and J. J. Simpson, Irradiance measurements in the upper       
!       ocean, J. Phys. Oceanogr., 7, 952-956, 1977.                                   
!                                                                                       
!       In addition, the module must be linked with netcdf library.
!
!       The GOTM option was taken from gotm.net.
!
!       A very special thanks to Dr. Tim Campbell, the author of MPI ELCIRC. We have
!       largely followed his style in MPI SELFE/SCHISM. We are indebted to his generous
!       help throughout the process of parallelizing SELFE.
!
!                                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OMP: search for 'new21' for TODO
!===============================================================================
!===============================================================================
! SCHISM main driver 
!===============================================================================
!===============================================================================
program schism_driver
  use schism_msgp, only: parallel_init,parallel_finalize,parallel_abort
  use schism_version
  implicit none
  character*8 cli
  call get_command_argument(1,cli)
  if (cli(1:2) == "-v")then
     print*, ""
     call print_version
     stop
  endif

  call parallel_init

  !Deal with command args
  !call get_command_args
  call schism_main
  call parallel_finalize
 
end program schism_driver

subroutine schism_main
  use schism_msgp, only: myrank !! debug only
  implicit none
  integer :: it,iths,ntime
  call schism_init(0,'./',iths,ntime)

  do it=iths+1,ntime
    call schism_step(it)
  enddo !it
  call schism_finalize
end subroutine schism_main

