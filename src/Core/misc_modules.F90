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

      use schism_glbl, only : rkind,errmsg,in_dir,out_dir,len_in_dir,len_out_dir
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
      open(31,file=in_dir(1:len_in_dir)//trim(fname),status='old')

      rewind(31)
      line=0
      do
        line=line+1
        read(31,'(a)',end=99)line_str
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
      close(31)
      return

99    close(31)
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

      end module misc_modules
