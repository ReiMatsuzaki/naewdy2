#include "macros.fpp"

module mod_sys
  implicit none
contains
  subroutine mkdir(path)
    character(*), intent(in) :: path
    character(100) cmd
    write(cmd, "( 'mkdir ', A )") trim(path)
    call system(cmd)
  end subroutine mkdir
  subroutine mkdir_if_not(path)
    character(*), intent(in) :: path
    character(1000) cmd
    write(cmd, "('ls ', A, ' > /dev/null 2>&1 || mkdir ', A )") trim(path), trim(path)
    call system(cmd)
  end subroutine mkdir_if_not
  subroutine mkdirp_if_not(path)
    character(*), intent(in) :: path
    character(1000) cmd
    write(cmd, "('ls ', A, ' > /dev/null 2>&1 || mkdir -p ', A )") trim(path), trim(path)
    call system(cmd)
  end subroutine mkdirp_if_not
  function path_exists(path) result(res)
    character(*), intent(in) :: path
    logical :: res
    integer, parameter :: ifile = 999999
    open(ifile, file=path, status="old", err=999)
    res = .true.
    close(ifile)
    return
    
999 continue
    res = .false.
    return
  
  end function path_exists
end module mod_sys

module Mod_csv
  ! csv writer 
  use mod_err_handle
  implicit none
  private
  integer ifile_
  integer :: num_
  character(100), allocatable :: labels_(:)
  integer col_
  public :: csv_new, csv_delete, csv_write, csv_write_row
contains
  subroutine csv_new(ifile, labels)
    integer, intent(in) :: ifile
    character(*), intent(in) :: labels(:)
    integer i
    ifile_ = ifile
    num_ = size(labels)
    allocate(labels_(num_))
    labels_(:) = labels(:)
    col_ = 0

    do i = 1, num_
       write(ifile_, "(A)", advance="no") trim(labels_(i))
       if(i .eq. num_) then
          write(ifile_,*)
       else
          write(ifile_,"(A)", advance="no") ","
       end if
    end do
    
  end subroutine csv_new
  subroutine csv_delete
    deallocate(labels_)
  end subroutine csv_delete
  subroutine csv_write_row(vals)
    double precision, intent(in) :: vals(:)
    integer i
    do i = 1, num_
       write(ifile_, "(f20.10)", advance="no") vals(i)
       if(i .eq. num_) then
          write(ifile_,*)
       else
          write(ifile_,"(A)", advance="no") ","
       end if
    end do
    
  end subroutine csv_write_row
  subroutine csv_write(v)
    double precision, intent(in) :: v

    write(ifile_, "(f20.10)", advance="no") v
    col_ = col_ + 1       
    if(col_ .eq. size(labels_)) then
       write(ifile_,*) ""
       col_ = 0
    else
       write(ifile_,"(A)",advance="no") ","
    end if    
    
  end subroutine csv_write
end module Mod_csv

