module netcdf_var_names
  
 !parameter (maxlen=200,long_name_len=100,dataset_name_len=8)
 integer, parameter :: maxlen=200,long_name_len=100,dataset_name_len=8 
 character(len=48)  :: outfile_name(maxlen)
 character(len=long_name_len):: var_long_name1(maxlen) ! netcdf long name for var component 1, scalar only needs this one
 character(len=long_name_len):: var_long_name2(maxlen) ! netcdf long name for var component 2, only used by vector
 character(len=long_name_len):: var_standard_name1(maxlen) ! cf standard name for var component 1, scalar only needs this one
 character(len=long_name_len):: var_standard_name2(maxlen) ! cf standard name for var component 2, only used by vector
 character(len=dataset_name_len):: dataset_name1(maxlen) ! netcdf var name for var component 1, scalar only needs this one
 character(len=dataset_name_len):: dataset_name2(maxlen) ! netcdf var name for var component 2, only used by vector
  character(len=long_name_len):: var_cf_standard_name1(maxlen) ! netcdf var name for var component 1, scalar only needs this one
 character(len=long_name_len):: var_cf_standard_name2(maxlen) ! netcdf var name for var component 2, only used by vector
 character(len=2)  :: outfile_type(10)
 character(len=4)  :: out_location(3)
 character(len=5)  :: level_location(2)
 integer           :: tracer_start_index,tracer_len,itmp
 integer           :: num_gen_tracers ! maximum number of general tracers
 integer           :: num_alge_tracers ! maximum number of general tracers
 integer           :: num_sed_tracers  ! maximum number of general tracers
 integer           :: num_eco_tracers ! maximum number of general tracers
 character(len=12) :: tracer_char
 character(len=7)  :: invalid_name
 
contains 

 subroutine fill_labels()
 
    num_gen_tracers=10
    num_alge_tracers=10
    num_sed_tracers=10
    num_eco_tracers=10
    invalid_name="INVALID"
        
     outfile_type(1) = "61"
     outfile_type(2) = "62"
     outfile_type(3) = "63"
     outfile_type(4) = "64"
     outfile_type(5) = "65"
     outfile_type(6) = "66"
     outfile_type(7) = "67"
     outfile_type(8) = "68"
     outfile_type(9) = "69"
     outfile_type(10) = "70"
 
     out_location(1)="node"
     out_location(2)="edge"
     out_location(3)="elem"
 
     level_location(1)="full"
     level_location(2)="half"
 
     outfile_name(1:maxlen)=invalid_name
     var_long_name1(1:maxlen)=invalid_name
     var_long_name2(1:maxlen)=invalid_name
     var_standard_name1(1:maxlen)=invalid_name
     var_standard_name2(1:maxlen)=invalid_name
     dataset_name1(1:maxlen)=invalid_name
     dataset_name2(1:maxlen)=invalid_name
     
     outfile_name(1)='elev'
     outfile_name(2)='pres'
     outfile_name(3)='airt'
     outfile_name(4)='shum'
     outfile_name(5)='srad'
     outfile_name(6)='flsu'
     outfile_name(7)='fllu'
     outfile_name(8)='radu'
     outfile_name(9)='radd'
     outfile_name(10)='flux'
     outfile_name(11)='evap'
     outfile_name(12)='prcp'
     outfile_name(13)='bdrc'
     outfile_name(14)='wind'
     outfile_name(15)='wist'
     outfile_name(16)='dahv'
     outfile_name(17)='vert'
     outfile_name(18)='temp'
     outfile_name(19)='salt'
     outfile_name(20)='conc'
     outfile_name(21)='tdff'
     outfile_name(22)='vdff'
     outfile_name(23)='kine'
     outfile_name(24)='mixl'
     outfile_name(25)='zcor'
     outfile_name(26)='qnon'
     outfile_name(27)='hvel' 
     outfile_name(28)='z0st'
     outfile_name(29)='z0eq'
     outfile_name(30)='z0cr'
     outfile_name(31)='z0sw'
     outfile_name(32)='z0wr'
     outfile_name(33)='bpgr'
     outfile_name(34)='wafo'
     outfile_name(35)='bdoc'
     outfile_name(36)='bnh4'
     outfile_name(37)='bno3'
     outfile_name(38)='bpo4'
     outfile_name(39)='bcod'
     outfile_name(40)='sbdo'
     outfile_name(41)='sbsa'
     outfile_name(42)='bthk'
     outfile_name(43)='bage'
     outfile_name(44)='dtbe'
     outfile_name(45)='tfu1'
     outfile_name(46)='tfu2'
     
     
     dataset_name1(1)='elev'
     dataset_name1(2)='pres'
     dataset_name1(3)='airt'
     dataset_name1(4)='shum'
     dataset_name1(5)='srad'
     dataset_name1(6)='flsu'
     dataset_name1(7)='fllu'
     dataset_name1(8)='radu'
     dataset_name1(9)='radd'
     dataset_name1(10)='flux'
     dataset_name1(11)='evap'
     dataset_name1(12)='prcp'
     dataset_name1(13)='bdrc'
     dataset_name1(14)='wind_u'
     dataset_name1(15)='wist_x'
     dataset_name1(16)='dahv_x'
     dataset_name1(17)='vert'
     dataset_name1(18)='temp'
     dataset_name1(19)='salt'
     dataset_name1(20)='conc'
     dataset_name1(21)='tdff'
     dataset_name1(22)='vdff'
     dataset_name1(23)='kine'
     dataset_name1(24)='mixl'
     dataset_name1(25)='zcor'
     dataset_name1(26)='qnon'
     dataset_name1(27)='hvel_u' 
     dataset_name1(28)='z0st'
     dataset_name1(29)='z0eq'
     dataset_name1(30)='z0cr'
     dataset_name1(31)='z0sw'
     dataset_name1(32)='z0wr'
     dataset_name1(33)='bpgr_x'
     dataset_name1(34)='wafo'
     dataset_name1(35)='bdoc'
     dataset_name1(36)='bnh4'
     dataset_name1(37)='bno3'
     dataset_name1(38)='bpo4'
     dataset_name1(39)='bcod'
     dataset_name1(40)='sbdo'
     dataset_name1(41)='sbsa'
     dataset_name1(42)='bthk'
     dataset_name1(43)='bage'
     dataset_name1(44)='dtbe'
     dataset_name1(45)='tfu1_u'
     dataset_name1(46)='tfu2_u'
     
     dataset_name2(14)='wind_v'
     dataset_name2(15)='wist_y'
     dataset_name2(16)='dahv_y'
     dataset_name2(27)='hvel_v'
     dataset_name2(33)='bpgr_y'
     dataset_name2(45)='tfu1_v'
     dataset_name2(46)='tfu2_v'
     
     tracer_start_index = 47
   
     do i=1,num_gen_tracers 
        write(tracer_char,'(i03)')i
        tracer_char=adjustl(tracer_char)  !place blanks at end
        tracer_len=len_trim(tracer_char)  !length without trailing blanks
        itmp= tracer_start_index+i-1
        outfile_name(itmp)='GEN_'//tracer_char(1:tracer_len)
        dataset_name1(itmp)=outfile_name(itmp)
        var_long_name1(itmp)='Generic tracer #'//trim(tracer_char)
     enddo !i
     
    tracer_start_index= tracer_start_index+num_gen_tracers
    do i=1,num_alge_tracers 
        write(tracer_char,'(i03)')i
        tracer_char=adjustl(tracer_char)  !place blanks at end
        tracer_len=len_trim(tracer_char)  !length without trailing blanks
        itmp= tracer_start_index+i-1
        outfile_name(itmp)='ALGE_'//tracer_char(1:tracer_len)
        dataset_name1(itmp)=outfile_name(itmp)
        var_long_name1(itmp)='ALGE tracer #'//trim(tracer_char)
    enddo !i
     
    
    tracer_start_index= tracer_start_index+num_alge_tracers
    do i=1,num_sed_tracers 
        write(tracer_char,'(i03)')i
        tracer_char=adjustl(tracer_char)  !place blanks at end
        tracer_len=len_trim(tracer_char)  !length without trailing blanks
        itmp= tracer_start_index+i-1
        outfile_name(itmp)='SED_'//tracer_char(1:tracer_len)
        dataset_name1(itmp)=outfile_name(itmp)
        var_long_name1(itmp)='SED 3D tracer #'//trim(tracer_char)
    enddo !i
    
    
    tracer_start_index= tracer_start_index+num_sed_tracers
    do i=1,num_sed_tracers 
        write(tracer_char,'(i03)')i
        tracer_char=adjustl(tracer_char)  !place blanks at end
        tracer_len=len_trim(tracer_char)  !length without trailing blanks
        itmp= tracer_start_index+i-1
        outfile_name(itmp)='SED_qbdl'//tracer_char(1:tracer_len)
        dataset_name1(itmp)=outfile_name(itmp)
        var_long_name1(itmp)='SED3D bedfr (top layer) #'//trim(tracer_char)
    enddo !i
    
    
      
    tracer_start_index= tracer_start_index+num_sed_tracers
    do i=1,num_sed_tracers 
        write(tracer_char,'(i03)')i
        tracer_char=adjustl(tracer_char)  !place blanks at end
        tracer_len=len_trim(tracer_char)  !length without trailing blanks
        itmp= tracer_start_index+i-1
        outfile_name(itmp)='SED_bfrac'//tracer_char(1:tracer_len)
        dataset_name1(itmp)=outfile_name(itmp)
        var_long_name1(itmp)='SED bedload #'//trim(tracer_char)
    enddo !i
     
   itmp = tracer_start_index+num_sed_tracers
   outfile_name(itmp)='SED_depth'
   dataset_name1(itmp)=outfile_name(itmp)
   var_long_name1(itmp)='SED depth in m'
   itmp=itmp+1
   outfile_name(itmp)='SED_bedd50'
   dataset_name1(itmp)=outfile_name(itmp)
   var_long_name1(itmp)='sed median grain size (mm)'
   itmp=itmp+1
   outfile_name(itmp)='SED_bstress.61'
   dataset_name1(itmp)=outfile_name(itmp)
   var_long_name1(itmp)='sed bottom shear stress (N.m-2)'
   itmp=itmp+1
   outfile_name(itmp)='SED_brough'
   dataset_name1(itmp)=outfile_name(itmp)
   var_long_name1(itmp)='sed bottom roughness lenght z0 (m)'
    
   
   var_long_name1(1)="water surface height above reference datum" 
   var_long_name1(17)="vertical velocity at face center at all levels"
   var_long_name1(18)="water temperature"
   var_long_name1(19)="water practical salinity" 
   var_long_name1(14)="eastward wind velocity" 
   var_long_name2(14)="northward wind velocity" 
   var_long_name1(15)="eastward wind stress" 
   var_long_name2(15)="northward wind stress" 
   var_long_name1(16)="eastward depth averaged horizontal velocity" 
   var_long_name2(16)="northward depth averaged horizontal velocity" 
   var_long_name1(25)="elevation of every node over the 3D mesh" 
   var_long_name1(27)="eastward water velocity" 
   var_long_name2(27)="northward water velocity" 
   var_long_name1(33)="eastward barotropic pressure gradient force " 
   var_long_name2(33)="northward barotropic pressure gradient force " 
   var_long_name1(44)="minimal time step allowed at each element" 
   var_long_name1(45)="eastward salt flux at face center" 
   var_long_name2(45)="northward salt flux at face center" 
   var_long_name1(46)="eastward temperature flux at face center" 
   var_long_name2(46)="northward temperature flux at face center" 
   
   var_standard_name1(1)="sea_surface_height_above_reference_ellipsoid" 
   var_standard_name1(17)="upward_sea_water_velocity" 
   var_standard_name1(18)="sea_water_temperature"
   var_standard_name1(19)="sea_water_practical_salinity" 
   var_standard_name1(14)="eastward_wind" 
   var_standard_name2(14)="northward_wind" 
   var_standard_name1(15)="eastward_wind_shear" 
   var_standard_name2(15)="northward_wind_shear" 
   var_standard_name1(16)="eastward_depth_averaged_velocity" !not from cf standard name table
   var_standard_name2(16)="northward_depth_averaged_velocity" !not from cf standard name table
   var_standard_name1(25)="elevation" !not from cf standard name table
   var_standard_name1(27)="eastward_sea_water_velocity" 
   var_standard_name2(27)="northward_sea_water_velocity" 
   var_standard_name1(33)="eastward_barotropic_pressure_gradient_force " 
   var_standard_name2(33)="northward_barotropic_pressure_gradient_force " 
   var_standard_name1(44)="minimal_time_step" !not from cf standard name table
   var_standard_name1(45)="eastward_salt_flux_face_center" 
   var_standard_name2(45)="northward_salt_flux_face_center" 
   var_standard_name1(46)="eastward_temperature_flux_face_center" 
   var_standard_name2(46)="northward_temperature_flux_face_center" 
   
 end subroutine
 
 
 subroutine get_netcdf_var_names(out_name,out_type,a_dataset_name,a_var_long_name,&
            & a_var_standard_name,a_out_location,a_level_location,a_var_dim,stat)

  implicit none
  character(*), intent(in)    :: out_name,out_type
  character(len=long_name_len),intent(out):: a_var_long_name(2) 
  character(len=long_name_len),intent(out):: a_var_standard_name(2) 
  character(len=dataset_name_len),intent(out):: a_dataset_name(2)  
  character(len=4),intent(out):: a_out_location
  character(len=5),intent(out):: a_level_location
  logical,intent(out)::stat
  integer,intent(out)::a_var_dim
  
  integer index,i,len
  character(10) type_temp
  
  
  
  a_dataset_name(1:2) = invalid_name
  a_var_long_name(1:2) = invalid_name
  out_location = invalid_name
  level_location = invalid_name
  a_var_dim=1
  
  call  fill_labels()
  index=-1
  do i=1,maxlen
      if (outfile_name(i).eq.out_name) then
         index = i
         exit
      endif
  enddo
  stat=.false.
  if (index.eq.-1) then
      return
  end if
  stat=.true.
  a_dataset_name(1) =  dataset_name1(index)
  a_dataset_name(2) =  dataset_name2(index)
  
  a_var_long_name(1) = var_long_name1(index)
  a_var_long_name(2) = var_long_name2(index)
  
  a_var_standard_name(1) = var_standard_name1(index)
  a_var_standard_name(2) = var_standard_name2(index)
  
  
  type_temp=adjustl(out_type)   
  
  if (type_temp(1:2).eq.'61') then
      a_var_dim =1
      a_out_location=out_location(1)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'62') then
      a_var_dim =2
      a_out_location=out_location(1)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'63') then
      a_var_dim =1
      a_out_location=out_location(1)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'64') then
      a_var_dim =2
      a_out_location=out_location(1)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'65') then
      a_var_dim =2
      a_out_location=out_location(2)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'66') then
      a_var_dim =1
      a_out_location=out_location(3)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'67') then
      a_var_dim =2
      a_out_location=out_location(2)
      a_level_location=level_location(1)
   else if (type_temp(1:2).eq.'68') then
      a_var_dim =2
      a_out_location=out_location(2)
      a_level_location=level_location(2)
  else if (type_temp(1:2).eq.'69') then
      a_var_dim =1
      a_out_location=out_location(3)
      a_level_location=level_location(1)
  else if (type_temp(1:2).eq.'70') then
      a_var_dim =1
      a_out_location=out_location(3)
      a_level_location=level_location(2)
  else
      a_var_dim =1
      a_out_location=out_location(1)
      a_level_location=level_location(1)
  endif
  
 end subroutine
 
 end module  netcdf_var_names
 
 
