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

!===============================================================================
!     The following routines/functions require explicit interface to 
!     use advanced FORTRAN features (e.g., optional argument, keyword)
!
!     subroutine get_param
!     subroutine schism_output_custom
!===============================================================================

      module misc_modules
      implicit none

      contains

!===============================================================================
!===============================================================================
! Routine to read in param.in-like inputs
!===============================================================================
!===============================================================================
      subroutine get_param(fname,varname,vartype,ivarvalue,varvalue1,varvalue2,ndim1,iarr1,arr1)
! Get a parameter from param.in
! Inputs:
!        fname: the file name (e.g. 'param.in')
!        varname: parameter name (string no longer than 90)
!        vartype: parameter value type (0: 2-char string; 1: integer; 2: double; 3: integer array; 4: double array)
!        optional, ndim1: actual dimension for arrays (<=10000)
! Outputs:
!        ivarvalue: integer output;
!        varvalue1: float output;
!        varvalue2: 2-char string output.
!        optional, iarr1(:),arr1(:): for array inputs (integer or double)

! Format rules for param.in:
! (1) Lines beginning with "!" are comments; blank lines are ignored;
! (2) one line for each parameter in the format: keywords= value;
!     keywords are case sensitive; spaces allowed between keywords and "=" and value;
!     comments starting with "!"  allowed after value;
! (3) value is an integer, double, 2-char string, integer| double arrays; for double, any of the format is acceptable:
!     40 40. 4.e1
!     Use of decimal point in integers is OK but discouraged. For
!     optional array outputs, make sure you use keywords; e.g.
!     call get_param('param.in','grain_size',4,itmp,tmp,stringvalue,ndim1=10,arr1=grain_size)
!     where grain_size(:) is an array declared in driver routine.

      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort,myrank
      implicit none

      character(*),intent(in) :: fname 
      character(*),intent(in) :: varname
      integer,intent(in) :: vartype
      integer,intent(out) :: ivarvalue
      real(rkind),intent(out) :: varvalue1
      character(len=2),intent(out) :: varvalue2
      integer,optional,intent(in) :: ndim1
      integer,optional,intent(out) :: iarr1(10000)
      real(rkind),optional,intent(out) :: arr1(10000)

      character(len=300) :: line_str,str_tmp,str_tmp2
      integer :: lstr_tmp,lstr_tmp2,line,len_str,loc,loc2

      str_tmp2=adjustl(varname)
      lstr_tmp2=len_trim(str_tmp2)
!      print*, varname !,str_tmp2(1:lstr_tmp2)

      ! Scan param.in
      open(15,file=fname,status='old')
      !open(15,file='param.in',status='old')

      rewind(15)
      line=0
      do
        line=line+1
        read(15,'(a)',end=99)line_str
        line_str=adjustl(line_str) !place blanks at end
        len_str=len_trim(line_str)
        if(len_str==0.or.line_str(1:1)=='!') cycle

        loc=index(line_str,'=')
        loc2=index(line_str,'!')
        if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: comments before =')
!'

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
            varvalue2=str_tmp2(1:2)
#ifdef DEBUG
            if(myrank==0) write(99,*)varname,' = ',varvalue2
#endif
          else if(vartype==1) then !integer
            read(str_tmp2,*)ivarvalue
#ifdef DEBUG
            if(myrank==0) write(99,*)varname,' = ',ivarvalue
#endif
          else if(vartype==2) then !double
            read(str_tmp2,*)varvalue1
#ifdef DEBUG
            if(myrank==0) write(99,*)varname,' = ',real(varvalue1)
#endif
          else !arrays
            if(.not.present(ndim1)) call parallel_abort('get_param: ndim1 not found')
            if(ndim1>10000) call parallel_abort('get_param: ndim1>10000')
!'

            if(vartype==3) then !integer array
              if(.not.present(iarr1)) call parallel_abort('get_param: iarr1 not found')
!'
              read(str_tmp2,*)iarr1(1:ndim1)
            else if(vartype==4) then !double array
              if(.not.present(arr1)) call parallel_abort('get_param: arr1 not found')
!'
              read(str_tmp2,*)arr1(1:ndim1)
            else
              write(errmsg,*)'get_param: unknown type:',vartype
              call parallel_abort(errmsg)
            endif
          endif
          exit
        endif
      enddo !scan param.in
  
!     print*, 'Found it on line: ',line
      close(15)
      return

99    close(15)
#ifdef USE_FABM
        if (varname(1:3)=='fab') then
          ivarvalue=1
        else
#endif
      write(errmsg,*)'Failed to find parameter:',varname
      call parallel_abort(errmsg)
#ifdef USE_FABM
      end if
#endif

      end subroutine get_param

!===============================================================================
!===============================================================================
      subroutine schism_output_custom(lwrite,i23d,ivs,ichanout,fname,idim1,idim2,outvar1,outvar2)
!-------------------------------------------------------------------------------
!     Custom outputs for generic use. Can be called from any routine, but make sure that
!     the calling routine is called inside the main time loop 
!     exactly ONCE per step! Also beware of output channel conflict - (210 - 240
!     seem to be available now.
!     These outputs share nspool and ihfskip with standard outputs.
!
!     Inputs:
!            i23d: indicates location where outputs are defined; 4 - 2D side; 
!                  5 - 2D element; 6,7 - 3D side and whole/half levels; 
!                  8,9 - 3D element and whole/half levels; 10 - 2D node; 11,12 -
!                  3D node and whole/half levels;
!            ivs: 1 - scalar; 2 - vector;
!            ichanout: output channel number (make sure there is no conflict);
!            fname: output file base name; must be 4 character in length (e.g. 'hvel');
!                   full name will be fname//'.65' etc where suffix is determined by
!                   i23d;
!            idim1,idim2: dimensions of output array(s) in the driving routine. 
!                         For 2D variables (e.g., bottom
!                         stress), idim1 must be 1; for 3D variables, idim1 must be nvrt.
!                         idim2 must be consistent with the type of output as given by
!                         i23d (e.g., idim2=nsa or ns for i23d=4);
!            outvar1(idim1,idim2): output array;
!            outvar2(idim1,idim2): optional output array for the case of ivs=2 (vectors).
!     Outputs:
!            lwrite: 0 - didn't output (not output step); 1 - output successfully.
!-------------------------------------------------------------------------------
      use schism_glbl, only: rkind,errmsg,ihfskip,nspool,time_stamp, &
     &it_main,iths_main,ifile_char,a_4,eta2,np,ne,ns,out_rkind !,fileopenformat
      use schism_msgp, only: parallel_abort,myrank
      implicit none
      integer,intent(in) :: i23d,ivs,ichanout,idim1,idim2
      character(len=4),intent(in) :: fname
      real(rkind),intent(in) :: outvar1(idim1,idim2)
      real(rkind),optional,intent(in) :: outvar2(idim1,idim2)
      integer,intent(out) :: lwrite
 
      character(len=3) :: sfix
      character(len=72) :: fgb
      character(len=140) :: fullname
      logical :: lex1,lex2
      integer :: lim_out,ifile,ifile_len,i,k,lfgb
      
!     Check optional argument
      if(ivs==2.and..not.present(outvar2)) then
        write(errmsg,*)'schism_output_custom: vector needs more info'
        call parallel_abort(errmsg)
      endif

!     Compute suffix for file name
!     Define upper limits for outer loop
      if(i23d==4) then
        sfix='.65'; lim_out=ns
      else if(i23d==5) then
        sfix='.66'; lim_out=ne
      else if(i23d==6) then
        sfix='.67'; lim_out=ns
      else if(i23d==7) then
        sfix='.68'; lim_out=ns
      else if(i23d==8) then
        sfix='.69'; lim_out=ne
      else if(i23d==9) then
        sfix='.70'; lim_out=ne
      else if(i23d==10) then
        sfix='.71'; lim_out=np
      else if(i23d==11) then
        sfix='.72'; lim_out=np
      else if(i23d==12) then
        sfix='.73'; lim_out=np
      else
        call parallel_abort('schism_output_custom: unknown name')
      endif

!     Open first stack
      if(it_main==iths_main+1) then
        ifile=(it_main-1)/ihfskip+1
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        fullname='outputs/'//(fgb(1:lfgb)//'_'//fname//sfix)

!        inquire(file=fullname,exist=lex1)
!        inquire(ichanout,exist=lex2)
!        if(.not.lex1.or.(lex1.and..not.lex2)) then
         open(ichanout,file=fullname,status='replace',form="unformatted",access="stream")
      endif !it_main

!     Return if not output step
      if(mod(it_main,nspool)/=0) then
        lwrite=0
        return
      endif
      lwrite=1

!     Write data at each step
!Error: lat/lon not working
      write(ichanout) real(time_stamp,out_rkind)
      write(ichanout) it_main
!     Additional info: ivs
      write(ichanout) ivs
      write(ichanout) (real(eta2(i),out_rkind),i=1,np)
      
      if(ivs==2.and.present(outvar2)) then
        write(ichanout) ((real(outvar1(k,i),out_rkind),real(outvar2(k,i),out_rkind),k=1,idim1),i=1,lim_out)
      else
        write(ichanout) ((real(outvar1(k,i),out_rkind),k=1,idim1),i=1,lim_out)
      endif
 
!     Open next stack to account for end of run
      if(mod(it_main,ihfskip)==0) then
        ifile=it_main/ihfskip+1
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        fullname='outputs/'//(fgb(1:lfgb)//'_'//fname//sfix)
        close(ichanout)
        open(ichanout,file=fullname,status='replace',form="unformatted",access='stream')
      endif !mod

      end subroutine schism_output_custom
!===============================================================================
!===============================================================================

      end module misc_modules
