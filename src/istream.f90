#include "macros.fpp"

module mod_istream
  use mod_err_handle
  implicit none
  integer, parameter :: MAX_LINE=1000
  integer, parameter :: MAX_MARKED=100
  type istream_file
     integer ifile
     integer index
     character(MAX_LINE) previous_line
     character(MAX_LINE) current_line     
  end type istream_file
  type istream_string
     integer index
     character(MAX_LINE) line
  end type istream_string
  type istream     
     character(MAX_MARKED) marked
     integer index_marked
     integer :: derived_index = 0
     type(istream_file), pointer :: pfile     
     type(istream_string), pointer :: pstr
  end type istream
contains
  ! ==== interface ====
  ! -- Constrcutor --
  subroutine istream_new_file(this, ifile, filename)
    type(istream)    :: this
    integer, intent(in) ::  ifile
    character(*)        :: filename

    if(this % derived_index .ne. 0) then
       throw_err("already allocated", 1)
    end if
    
    this % marked = ""
    this % index_marked = -1
    this % derived_index = 1    
    allocate(this % pfile)    
    call istream_file_new(this % pfile, ifile, filename)
    
  end subroutine istream_new_file
  subroutine istream_new_string(this, str)
    type(istream) :: this
    character(*)  :: str

    if(this % derived_index .ne. 0) then
       throw_err("already allocated", 1)
    end if
    
    this % marked = ""
    this % index_marked = -1
    this % derived_index = 2    
    allocate(this % pstr)    
    call istream_string_new(this % pstr, str)

  end subroutine istream_new_string
  subroutine istream_delete(this)
    type(istream) :: this

    if(this % derived_index .eq. 0) then
       throw_err("not allocated.", 1)
    end if

    if(this % derived_index .eq. 1) then
       call istream_file_delete(this % pfile)
       deallocate(this % pfile)
    else if(this % derived_index .eq. 2) then
       call istream_string_delete(this % pstr)
       deallocate(this % pstr)
    end if
    this % derived_index = 0
    
  end subroutine istream_delete
  ! -- main --
  subroutine istream_forward(this, c, is_end)
    type(istream) :: this
    character        :: c
    logical          :: is_end
    integer :: midx

    if(this % derived_index .eq. 0) then
       throw_err("not allocated.", 1)
    else if(this % derived_index .eq. 1) then
       call istream_file_forward(this % pfile, c, is_end)
    else if(this % derived_index .eq. 2) then       
       call istream_string_forward(this % pstr, c, is_end)
    end if

    if(.not. is_end) then
       ! -- work mark --
       midx = this % index_marked
       if(midx .ne. -1) then
          this % marked(midx:midx+1) = c
          this % index_marked = midx+1
       end if
    end if

  end subroutine istream_forward
  subroutine istream_return(this)
    type(istream) :: this
    integer       :: midx

    if(this % derived_index .eq. 0) then
       throw_err("not allocated.", 1)
    else if(this % derived_index .eq. 1) then
       call istream_file_return(this % pfile)
    else if(this % derived_index .eq. 2) then       
       call istream_string_return(this % pstr)
    end if

    ! -- work mark --
    if(this % index_marked .eq. 1) then
       throw_err("previous index for mark_index is nothing", 1)
       return
    else if(this % index_marked .ne. -1) then
       midx = this % index_marked
       this % index_marked = midx-1
       this % marked(midx-1:midx+1) = "  "
    end if

  end subroutine istream_return
  subroutine istream_set_mark(this)
    type(istream) :: this
    
#if defined DEBUG_ISTREAM
    write(*, *) "start istream_set_mark"
#endif    
    
    if(this % index_marked .ne. -1) then
       throw_err("do not call when mark is already activated",1)
    end if
    this % marked = ""
    this % index_marked = 1

#if defined DEBUG_ISTREAM
    write(*, *) "  (index, mark_index)", this%index, this%index_marked
    write(*, *) "end istream_set_mark"
#endif


  end subroutine istream_set_mark
  subroutine istream_get_marked(this, out_marked)
    type(istream) :: this
    character(*) :: out_marked

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_get_marked"
#endif
    
    if(this % index_marked .eq. -1) then
       throw_err("do not call get_mark when mark is not activated", 1)
    end if
    out_marked = this % marked
    this % index_marked = -1

#if defined DEBUG_ISTREAM
    write(*, *) "  (index, mark_index, out_marked)", this%index, this%index_marked,out_marked
    write(*, *) "end istream_get_marked"
#endif
    
  end subroutine istream_get_marked
  ! ==== File ====
  ! -- Constrcutors --
  subroutine istream_file_new(this, ifile, filename)
    type(istream_file)    :: this
    integer, intent(in) ::  ifile
    character(*)        :: filename

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_file_new"
#endif
    
    call open_r(ifile, filename); check_err()
    this % ifile = ifile
    this % index = -1
    this % previous_line = ""
    this % current_line = ""

  end subroutine istream_file_new
  subroutine istream_file_delete(this)    
    type(istream_file) :: this

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_file_delete"
#endif    
    
    close(this % ifile)

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_file_delete"
#endif    
  end subroutine istream_file_delete
  ! -- main --  
  subroutine istream_file_forward(this, c, is_end)
    type(istream_file) :: this
    character        :: c
    logical          :: is_end
    integer :: idx

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_file_forward"
#endif
    
    idx = this % index
    if(idx > len_trim(this%current_line) .or. idx.eq.-1) then
       call read_new_line(this, is_end)
       if(is_end) then
          return
       end if
    end if
    
    idx = this % index
    is_end = .false.
    c = this % current_line(idx:idx+1)
    this % index = idx + 1    

#if defined DEBUG_ISTREAM
    write(*, *) "  (index, mark_index, c)", this%index, this%index_marked, c
#endif
    
  end subroutine istream_file_forward
  subroutine istream_file_return(this)
    type(istream_file) :: this

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_file_return"
#endif    
    
    if(this % index < 1) then
       begin_err(1)
       write(*,*) "index must be positive."
       call err_si("index: ", this % index)
       end_err()
    else if(this % index == 1 ) then
       if(this % previous_line .eq. "") then
          throw_err("previous line is not stored. this routine fail.",1)
       end if
       this % index = len_trim(this % previous_line)
       this % current_line = this % previous_line
       this % previous_line = ""
    else
       this % index = this % index - 1
    end if
    
#if defined DEBUG_ISTREAM
    write(*, *) "end istream_file_return"
    !write(*, *) "  (index, mark_index)", this%index, this%index_marked
#endif
            
  end subroutine istream_file_return
  ! ==== String ====
  ! -- Constructors --
  subroutine istream_string_new(this, str)
    type(istream_string) :: this
    character(*)         :: str
    if(len(str) > MAX_LINE) then
       begin_err(1)
       call err_1("input string is too long")
       call err_si("MAX_LINE: ", MAX_LINE)
       call err_si("len(str): ", len(str))
       end_err()
    end if
    this % index = 1
    this % line = str
  end subroutine istream_string_new
  subroutine istream_string_delete(this)
    type(istream_string) :: this
    this % index = 0
    this % line = ""
  end subroutine istream_string_delete
  ! -- main --
  subroutine read_new_line(this, is_end)
    type(istream_file) :: this
    logical :: is_end
    
#if defined DEBUG_ISTREAM
    write(*, *) "start read_new_line"
#endif
    
    this % previous_line = this % current_line
    read(this % ifile, '(a)', end=999) this % current_line
    is_end = .false.
    this % index = 1
    return
999 continue
    is_end = .true.
    return
  end subroutine read_new_line
  subroutine istream_string_forward(this, c, is_end)
    type(istream_string) :: this
    character        :: c
    logical          :: is_end
    integer :: idx

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_string_forward"
#endif
    is_end = .false.
    idx = this % index
    if(idx > len_trim(this % line)) then
       is_end = .true.
       return
    end if
    c = this % line(idx:idx+1)
    this % index = idx + 1    

#if defined DEBUG_ISTREAM
    write(*, *) "end istream_string_forward"
    write(*, *) "  (index,c)", this%index, c
#endif
    
  end subroutine istream_string_forward
  subroutine istream_string_return(this)
    type(istream_string) :: this

#if defined DEBUG_ISTREAM
    write(*, *) "start istream_string_return"
#endif    
    
    if(this % index < 1) then
       begin_err(1)
       call err_1("index must be positive.")
       call err_si("index: ", this % index)
       end_err()
    else if(this % index == 1 ) then
       throw_err("previous line is not stored. this routine fail.",1)
    end if
    
    this % index = this % index - 1
    
#if defined DEBUG_ISTREAM
    write(*, *) "end istream_string_return"
    !write(*, *) "  (index, mark_index)", this%index, this%index_marked
#endif
            
  end subroutine istream_string_return
end module mod_istream


