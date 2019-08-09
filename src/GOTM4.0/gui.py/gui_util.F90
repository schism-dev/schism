!$Id: gui_util.F90,v 1.2 2007-01-19 09:40:25 jorn Exp $

#include"cppdefs.h"

module gui_util
    implicit none
    
    logical, private :: redirected = .false.

contains

    subroutine redirectoutput(outpath,errpath)
        implicit none
        
        character(len=255), intent(in) :: outpath
        character(len=255), intent(in) :: errpath
        
        call resetoutput()
        open(unit=stdout,file=outpath)
        open(unit=stderr,file=errpath)
        redirected = .true.
    end subroutine redirectoutput

    subroutine resetoutput()
        implicit none
    
        if (redirected) then
            close(unit=stdout)
            close(unit=stderr)
            redirected = .false.
        end if
    end subroutine resetoutput
    
    subroutine getversion(versionstring)
        implicit none
        
        character(len=255), intent(out) :: versionstring
        
        versionstring = RELEASE
    end subroutine getversion
    
end module gui_util
