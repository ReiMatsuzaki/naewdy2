module Mod_StrUtil
contains
  subroutine is_i(str, res)
    use mod_err_handle
    implicit none
  
    character(*) str
    double precision  :: xd
    integer :: xi  
    logical res
    integer ierr_d, ierr_i, idx
 
    double precision eps
    
    eps = 1.0D-10
    read(str, *, iostat=ierr_d) xd
    read(str, *, iostat=ierr_i) xi
    
    if(ierr_i .ne. 0) then
       ! -- failed to convert to integer --
       res = .false.
       return
    end if
    
    if(ierr_d .ne. 0) then
       ! -- something wrong --
       throw_err("somthing wrong in is_i", 1)
       return
    end if
    
    if(abs(xd-xi) > eps) then
       ! -- different value for xd and xi --
       res = .false.
       return
    end if
  
    idx = index(str, ".")
    if(idx .eq. 0) then
       res = .true.
       return
    else
       res = .false.
       return
    end if

  end subroutine is_i
  subroutine is_d(str, res)
    character(*) str
    logical res
    double precision a
    read(str, *, err=998) a
    res = .true.
    return
998 continue
    res = .false.
    return
  end subroutine is_d
  subroutine convert_i(str, a)
    use mod_err_handle
    character(*) str
    integer   a
    
  
    double precision d
    double precision eps
    
    eps = 1.0D-10
    read(str, *, err=999) d
    read(str, *, err=999) a
    if(abs(d-a) > eps) then
       begin_err(1)
       call err_ss("input data may be real number\n", "");
       call err_ss("str: ", str)
       end_err();
    end if
    return
999 continue
    begin_err(1)
    call err_ss("failed to convert to integer", "")
    call err_ss("str: ", str)
    end_err()
  end subroutine convert_i
  subroutine convert_d(str, a)
    use mod_err_handle
    character(*) str
    double precision    a

    read(str, *, err=999) a
    return
999 continue
    begin_err(1)
    write(0,*) "failed to convert to real"
    write(0,*) "str", str
    end_err()
  end subroutine convert_d
end module Mod_StrUtil
