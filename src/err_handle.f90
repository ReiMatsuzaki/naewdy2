#include "macros.fpp"
!     Error handle utilities
module mod_err_handle
  implicit none  
  logical show_message_q
  integer ierr
contains
  subroutine err_handle_begin
    show_message_q = .true.
    ierr = 0
  end subroutine err_handle_begin
  subroutine err_handle_end

  end subroutine err_handle_end
  subroutine err_handle_on_message
    show_message_q = .true.
  end subroutine err_handle_on_message
  subroutine err_handle_off_message
    show_message_q = .false.
  end subroutine err_handle_off_message
  subroutine err_with_file_line(msg, ierr0, file, line)
    character(*) msg
    integer ierr0
    character(*) file
    integer line
    if(show_message_q) then
       write(0, '(A, ":", I0, ": ", A)') file, line, msg
    end if
    ierr = ierr0
  end subroutine err_with_file_line
  subroutine err_1(a)
    character(*) a
    if(show_message_q) then
       write(0, *) a
    end if
  end subroutine err_1
  subroutine err_ss(a, b)
    character(*) a, b
    if(show_message_q) then
       write(0, *) a, b
    end if
  end subroutine err_ss
  subroutine err_si(a, b)
    character(*) a
    integer b
    if(show_message_q) then
       write(0, *) a, b
    end if
  end subroutine err_si
  subroutine err_sd(a, b)
    character(*) a
    double precision :: b
    if(show_message_q) then
       write(0, *) a, b
    end if
  end subroutine err_sd
  subroutine open_w(ifile, filename)
    integer, intent(in) ::  ifile
    character(*), intent(in) ::  filename
    
    open(ifile, file=filename, status="replace", err=999)
    return
    
999 continue
    call err_with_file_line("failed to open file", 1, __FILE__, __LINE__)
    call err_ss("filename: ", filename)
    return
  end subroutine open_w
  subroutine open_r(ifile, filename)
    integer ifile
    character(*) filename
    
    open(ifile, file=filename, status='old',err=999)
    return
    
999 continue
    call err_with_file_line("failed to open file", 1, __FILE__, __LINE__)
    call err_ss("filename: ", filename)
    return
    
  end subroutine open_r
end module mod_err_handle

