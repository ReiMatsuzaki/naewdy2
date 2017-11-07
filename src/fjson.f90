#include "macros.fpp"

!
! fourtran 90 JSON parser
!      Rei Matsuzaki @Kyoto, Fukui institute
!
! ==== todo ====
! name of conversion subroutines
! more sophisticated error display.
! .   for example, print near information
! value_delete produce double free error.
! rename elements to token
! rewrite parser as object oriented programing
! remove number limitation for elements
!    as fixed parameter ?
!    or linked list
!    or vector
! remove number limitation for object
!    as fixed parameter
!    or linked list of (string, value)

module mod_fjson_parser
  use mod_err_handle
  use mod_istream
  implicit none

  type token
     integer type
     character(100) val
  end type token
  character(10), parameter :: type_names(10) = &
       (/"null   ", "}      ", "[      ", "]      ", ":      ", &
         ",      ", "number ", "string ", "{      ", "boolean"/)
  
  integer, parameter :: TOKEN_NULL=1  
  integer, parameter :: TOKEN_END_OBJ=2
  integer, parameter :: TOKEN_BEGIN_ARY=3
  integer, parameter :: TOKEN_END_ARY=4  
  integer, parameter :: TOKEN_SEP=5
  integer, parameter :: TOKEN_COMMA=6
  integer, parameter :: TOKEN_NUM=7
  integer, parameter :: TOKEN_STR=8
  integer, parameter :: TOKEN_BEGIN_OBJ=9
  integer, parameter :: TOKEN_BOOL=10
  
  integer, parameter :: MAX_ELE=1000

  type(token), save ::  tokens(MAX_ELE)
  integer, private :: ie
  integer, private, save :: status = -1

  private check_new
contains
  ! -- Constructors --
  subroutine parser_new
    integer iie

#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parser_new"
#endif
    
    ie = 1
    status = 1    
    do iie = 1, MAX_ELE       
       tokens(iie) % val = ""
       tokens(iie) % type = TOKEN_NULL
    end do

#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parser_new"
#endif    
    
  end subroutine parser_new
  subroutine parser_delete
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parser_delete"
#endif
    
    status = -1
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parser_delete"
#endif
  end subroutine parser_delete
  ! -- Error --
  subroutine check_new
    if(status .ne. 1) then
       throw_err("call parser_new first", 1)
       return
    end if
  end subroutine check_new
  ! -- Main --
  subroutine parser_parse(istream_json)
    type(istream) istream_json
    integer idx_num, idx_bool
    character c
    logical is_end
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parser_parse"
#endif    

    call check_new
    
    ! == parse ==
    is_end = .false.
    do
       call istream_forward(istream_json, c, is_end)
       if(is_end) then
          exit
       end if

#if defined DEBUG_FJSON_PARSER
       write(*, '("   ", i0, " : ", a)') ie, c
#endif       

!       if(debug_lvl > 0) then
!          write(*, *) i, c
!       end if
       idx_num = index("+-0123456789", c)
       idx_bool = index("tTfF", c)
       if(c=='{') then
          call parse_begin_obj; check_err()          
       else if(c=='}') then
          call parse_end_obj; check_err()
       else if(c=='[') then
          call parse_begin_ary; check_err()
       else if(c==']') then
          call parse_end_ary; check_err()
       else if(c==',') then
          call parse_comma; check_err()
       else if(c==':') then
          call parse_sep; check_err()
       else if(idx_num.ne.0) then
          call parse_num(istream_json); check_err()
       else if(idx_bool.ne.0) then
          call parse_bool(istream_json); check_err()
       else if(c=='"') then
          call parse_str(istream_json); check_err()
       else if(c==' ') then
       else if(iachar(c)==9) then
          ! ignore tab
       else if(c=="\n") then
       else if(trim(c)=="") then
       else
          begin_err(1)
          call err_1("invalid character")
          call err_ss("character: ", c)
          end_err()
       end if
    
    end do

#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parser_parse"
#endif    
    
  end subroutine parser_parse
  subroutine parse_begin_obj
    integer prev

#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_begin_obj"
#endif    

    
    if(ie .ne. 1) then
       prev = tokens(ie-1) % type
       if((prev.ne.TOKEN_BEGIN_ARY) .and. (prev.ne.TOKEN_SEP) .and. &
            (prev.ne.TOKEN_COMMA)) then
          begin_err(1)
          call err_1("invalid location for {")
          call err_ss("prev_val:  ", tokens(ie-1)%val)
          call err_si("prev_type: ", tokens(ie-1)%type)
          call err_si("ie: ", ie)
          end_err()
       end if
    end if
    
    tokens(ie) % type = TOKEN_BEGIN_OBJ
    tokens(ie) % val = "{"
    ie = ie+1

#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_begin_obj"
#endif        
  end subroutine parse_begin_obj
  subroutine parse_end_obj
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_end_obj"
#endif    
    
    tokens(ie) % type = TOKEN_END_OBJ
    tokens(ie) % val = "}"
    ie = ie+1

#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_begin_obj"
#endif        
  end subroutine parse_end_obj
  subroutine parse_begin_ary
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_begin_ary"
#endif        
    tokens(ie) % type = TOKEN_BEGIN_ARY
    tokens(ie) % val = "["
    ie = ie+1
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_end_ary"
#endif            
  end subroutine parse_begin_ary
  subroutine parse_end_ary
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_end_ary"
#endif            
    tokens(ie) % type = TOKEN_END_ARY
    tokens(ie) % val = "]"
    ie = ie+1
#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_end_ary"
#endif            
  end subroutine parse_end_ary
  subroutine parse_comma
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_comma"
#endif            
    tokens(ie) % type = TOKEN_COMMA
    tokens(ie) % val = ","
    ie = ie+1
#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_comma"
#endif                
  end subroutine parse_comma
  subroutine parse_sep
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_sep"
#endif                
    tokens(ie) % type = TOKEN_SEP
    tokens(ie) % val = ":"
    ie = ie+1
#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_sep"
#endif                    
  end subroutine parse_sep
  subroutine parse_num(istream_json)
    type(istream) :: istream_json
    integer idx
    character c
    logical is_end
    logical find_end_num
    character(100) :: line

#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_num"
#endif                
    
    call istream_return(istream_json)
    call istream_set_mark(istream_json)
    find_end_num = .false.
    is_end = .false.
    do
       call istream_forward(istream_json, c, is_end)
       if(is_end) then
          find_end_num = .true.
          exit
       end if
#if defined DEBUG_FJSON_PARSER
       write(*, '("   ", i0, " : ", a)') ie, c
#endif
       idx = index("+-0123456789.", c)
       if(idx.eq.0) then
          call istream_return(istream_json)
          find_end_num = .true. 
          exit
       end if
    end do

    if(.not. find_end_num) then
       throw_err("failed find end of number", 1)
    end if
    
    call istream_get_marked(istream_json, line)
    tokens(ie) % type = TOKEN_NUM
    tokens(ie) % val = line
    ie = ie + 1
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "val: ",  trim(line)
    write(*, *) "end parse_num"
#endif                             
    
  end subroutine parse_num
  subroutine parse_bool(istream_json)
    type(istream) :: istream_json
    logical       :: is_end
    character     :: c
    integer       :: i, num
    logical       :: current_res
    character(10) :: line

#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_bool"
#endif

    call istream_return(istream_json)
    call istream_set_mark(istream_json); check_err()
    call istream_forward(istream_json, c, is_end);   check_err()
    if(c .eq. "t" .or. c.eq. "T") then
       num = 4
       current_res = .true.
    else if(c .eq. "f" .or. c .eq. "F") then
       num = 5
       current_res = .false.
    else
       begin_err(1)
       call err_1("invalid first character for bool")
       call err_ss("c: ", c)
       end_err()
    end if
    
    do i = 1, num-1
       call istream_forward(istream_json, c, is_end); check_err()       
    end do

    call istream_get_marked(istream_json, line); check_err()
    if(line.eq."true" .or. line.eq."True" .or. line.eq."TRUE") then
       current_res = .true.
    else if(line.eq."false" .or. line.eq."False" .or. line.eq."FALSE") then
       current_res = .false.
    else
       begin_err(1)
       call err_1("invalid key for boolean")
       call err_ss("line: ", line)
       return
    end if

    tokens(ie) % type = TOKEN_BOOL
    if(current_res) then
       tokens(ie) % val  = "true"
    else
       tokens(ie) % val  = "false"
    end if
    ie = ie + 1
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "end parse_bool"
#endif    
  end subroutine parse_bool
  subroutine parse_str(istream_json)
    type(istream) :: istream_json
    integer idx
    logical is_end
    character(100) line
    character c
    
#if defined DEBUG_FJSON_PARSER
    write(*, *) "start parse_str"
#endif
    
    call istream_set_mark(istream_json); check_err()
    is_end = .false.
    do
       call istream_forward(istream_json, c, is_end); check_err()
       if(is_end) then
          exit
       end if

#if defined DEBUG_FJSON_PARSER
       write(*, '("   ", i0, " : ", a)') ie, c
#endif       

       idx = index('"', c)
       if(idx.ne.0) then
          call istream_return(istream_json)
          call istream_get_marked(istream_json, line)
          call istream_forward(istream_json, c, is_end)
          tokens(ie) % type = TOKEN_STR
          tokens(ie) % val = line
          ie = ie + 1
#if defined DEBUG_FJSON_PARSER
          write(*, *) "end parse_str"
#endif                
          return
       end if
    end do

    throw_err("failed to find end of string", 1)
    

    
!    do i = i+1, len_trim(str_json)
!       c = str_json(i:i+1)
!       if(c.eq.'"') then
!          tokens(ie) % type = TOKEN_STR
!          tokens(ie) % val  = str_json(i0+1:i-1)
!          ie = ie+1
!          return
!       end if
!    end do
    
  end subroutine parse_str
  ! -- IO --
  subroutine parser_dump(ifile)
    integer ifile
    integer iie
    character(100) s
    integer t
    call check_new; check_err()
    do iie =1, ie-1
       s = tokens(iie) % val
       t = tokens(iie) % type
       if(t .eq. TOKEN_STR) then
          write(ifile, '(i0, " " , a, a, a)') iie, '"', trim(s), '"'
       else if (t .eq. TOKEN_NUM) then
          write(ifile, '(i0, " " , a)') iie, trim(s)
       else if(t .eq. TOKEN_BOOL) then
          write(ifile, '(i0, " ", a)') iie, trim(s)
       else
          write(ifile, '(i0, " ", a)') iie, type_names(t)
       end if
    end do
  end subroutine parser_dump
end module mod_fjson_parser

module mod_fjson_value
  use mod_err_handle
  implicit none
  character(10), private :: json_type_name(7) = (/"N", "I", "D", "B", "S", "A", "O"/)
  integer, parameter :: TYPE_VALUE_NULL=1
  integer, parameter :: TYPE_I=2
  integer, parameter :: TYPE_D=3
  integer, parameter :: TYPE_B=4
  integer, parameter :: TYPE_S=5
  integer, parameter :: TYPE_A=6
  integer, parameter :: TYPE_O=7  
  integer, parameter :: MAX_ELE=2000
  integer, private  :: depth = 0

  type value
     integer :: type = 1
     integer val_i
     double precision  val_d
     character(100), pointer :: val_s
     logical val_b
     type(array),  pointer :: val_a
     type(object), pointer :: val_o
  end type value
  type object
     character(100) keys(MAX_ELE)
     type(value) vals(MAX_ELE)
     integer     :: size = 0
  end type object
  type array
     type(value) vals(MAX_ELE)
     integer  :: size = 0
  end type array
contains
  ! ==== utils ====
  subroutine insert_depth(ifile)
    integer :: ifile, j
    do j=1,1+depth
       write(ifile,fmt='(A)',advance='no') "   "
    end do
  end subroutine insert_depth

  ! ==== value ====
  ! -- Constructor --
  subroutine value_new_i(v, x)
    integer, intent(in) :: x
    type(value) :: v
    v % type = TYPE_I
    v % val_i = x    
  end subroutine value_new_i
  subroutine value_new_d(v, x)
    double precision, intent(in) :: x
    type(value) :: v
    v % type = TYPE_D
    v % val_d = x
  end subroutine value_new_d
  subroutine value_new_b(v, x)
    type(value)         :: v
    logical, intent(in) :: x
    v % type = TYPE_B
    v % val_b = x
  end subroutine value_new_b
  subroutine value_new_s(v, x)
    character(*) x
    type(value) :: v
    v%type = TYPE_S
    allocate(v%val_s)
    v%val_s = trim(x)
  end subroutine value_new_s
  subroutine value_new_a(v, x)
    type(value) :: v
    type(array) :: x
    v % type = TYPE_A
    allocate(v % val_a)
    !v % val_a = x
    call array_clone(v%val_a, x); check_err()
  end subroutine value_new_a
  subroutine value_new_o(v, x)
    type(value)  :: v
    type(object) :: x    
    v % type = TYPE_O
    allocate(v % val_o)
    call object_clone(v%val_o, x); check_err()
  end subroutine value_new_o
  subroutine value_clone(v, x)
    type(value) :: v, x
    select case(x%type)
    case(TYPE_VALUE_NULL)
    case(TYPE_I)
       call value_new_i(v, x%val_i); check_err()
    case(TYPE_D)
       call value_new_d(v, x%val_d); check_err()
    case(TYPE_B)
       call value_new_b(v, x%val_b); check_err()
    case(TYPE_S)
       call value_new_s(v, x%val_s); check_err()
    case(TYPE_A)
       call value_new_a(v, x%val_a); check_err()
    case(TYPE_O)
       call value_new_o(v, x%val_o); check_err()
    case default
       begin_err(1)
       write(0,*) "invalid type for x"
       write(0,*) "inputted type id:", x%type
       end_err()
    end select
  end subroutine value_clone
  recursive subroutine value_delete(v)
    type(value) :: v
    select case(v % type)
    case(TYPE_I)
    case(TYPE_D)
    case(TYPE_B)
    case(TYPE_S)
       deallocate(v % val_s)
    case(TYPE_A)
       call array_delete(v % val_a); check_err()
       deallocate(v % val_a)
    case(TYPE_O)
       call object_delete(v % val_o); check_err()
       deallocate(v % val_o)
!    case default
!       begin_err(1)
!       write(0,*) "invalid type"
!       write(0,*) "value.type:", v%type
!       end_err()
    end select
    v % type = TYPE_VALUE_NULL
  end subroutine value_delete
  ! -- Accessor --
  subroutine value_get_i(v, i)
    type(value) :: v
    integer :: i
    if(v % type .ne. TYPE_I) then
       begin_err(1)
       call err_1("value is not integer")
       call err_ss("inputed type is: ", json_type_name(v % type))
       return
    end if
    i = v % val_i
  end subroutine value_get_i
  subroutine value_get_d(v, d)
    type(value) :: v
    double precision :: d
    if(v % type .ne. TYPE_D) then
       begin_err(1)
       call err_1("value is not double")
       call err_ss("inputed type is: ", json_type_name(v % type))
       call value_dump(v, 6)
       return
    end if
    d = v % val_d
  end subroutine value_get_d
  subroutine value_get_b(v, x)
    type(value) :: v
    logical :: x
    if(v % type .ne. TYPE_B) then
       begin_err(1)
       call err_1("value is not double")
       call err_ss("inputed type is: ", json_type_name(v % type))
       return
    end if
    x = v % val_b
  end subroutine value_get_b
  subroutine value_get_s(v, x)
    type(value) :: v
    character(*) :: x
    if(v % type .ne. TYPE_S) then
       begin_err(1)
       call err_1("value is not string")
       call err_ss("inputed type is: ", json_type_name(v % type))
       return
    end if
    x = trim(v % val_s)
  end subroutine value_get_s
  subroutine value_get_a(v, x)
    type(value)  :: v
    type(array) :: x
    if(v % type .ne. TYPE_a) then
       begin_err(1)
       call err_1("value is not array")
       call err_ss("inputed type is: ", json_type_name(v % type))
       return
    end if
    call array_clone(x, v%val_a); check_err()
  end subroutine value_get_a
  subroutine value_get_o(v, o)
    type(value) :: v
    type(object) :: o
    if(v % type .ne. TYPE_o) then
       begin_err(1)
       call err_1("value is not object")
       call err_si("inputed type id is ", v % type)
       call err_ss("inputed type is: ", json_type_name(v % type))
       end_err()
    end if
    call object_clone(o, v%val_o); check_err()
  end subroutine value_get_o
  subroutine value_get_as_c(v, x)
    type(value) :: v
    complex(kind(0d0)) :: x
!    type(object) :: o
    double precision :: re, im

    if(v % type .eq. TYPE_O) then
       !    call value_get_o(v, o); check_err()
       !call object_get_d(o, "re", re); check_err()
       !call object_get_d(o, "im", im); check_err()
       call object_get_d(v%val_o, "re", re); check_err()
       call object_get_d(v%val_o, "im", im); check_err()
       x = cmplx(re, im, kind(0d0))
    else if(v % type .eq. TYPE_D) then
       call value_get_d(v, re); check_err()
       x = re
    else
       begin_err(1)       
       if(show_message_q) then
          write(*,*) "invalid type"
          write(*, '(a)', advance='no') "v = "
          call value_dump(v, 6)
          write(*, *)
       end if
       return
    end if
    
  end subroutine value_get_as_c
  ! -- IO --
  recursive subroutine value_dump(v, in_ifile)
    type(value) :: v
    integer, optional :: in_ifile
    integer :: ifile
    
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    select case(v % type)
    case (TYPE_I)
       write(ifile, "(i0)", advance='no') v % val_i
    case (TYPE_D)
       write(ifile, "(f20.10)", advance='no') v % val_d
    case (TYPE_B)
       if(v % val_b) then
          write(ifile, '(a)', advance='no') "true"
       else
          write(ifile, '(a)', advance='no') "false"
       end if
    case (TYPE_S)
       write(ifile, "(a,a,a)", advance='no') '"' ,v%val_s(:len_trim(v%val_s)),'"'
    case(TYPE_A)
       call array_dump(v % val_a, ifile); check_err()
    case (TYPE_O)
       call object_dump(v % val_o, ifile) ; check_err()
    case default
       begin_err(1)
       write(0,*) "this type is not supported now"
       write(0,*) "type:", v%type
       write(0,*) "val_i:", v%val_i
       write(0,*) "val_d:", v%val_d
       if(associated(v%val_s)) then
          write(0,*) "val_s:", v%val_s
       end if
    end select
  end subroutine value_dump
  subroutine value_dump_debug(v, ifile)
    type(value) :: v
    integer :: ifile
    write(ifile, *) "value:"
    select case(v % type)
    case (TYPE_I)
       write(ifile, "('(',i0, ', ' ,i0,')')") v % type, v % val_i
    case (TYPE_D)
       write(ifile, "('(',i0, ', ', f20.10,')')")  v % type, v % val_d
    case (TYPE_S)
       write(ifile, "('(',i0,a,a,a,')')")  v % type, '"' ,v%val_s(:len_trim(v%val_s)),'"'
    case default
       throw_err("this type is not supported now", 1)
    end select
  end subroutine value_dump_debug

  ! ==== object ====
  ! -- Constructor --
  subroutine object_new(o)
    type(object) o
    o % size = 0
  end subroutine object_new
  subroutine object_clone(this, other)
    type(object) this, other
    integer i
    this%size = other%size
    do i = 1, this%size
       call value_clone(this%vals(i), other%vals(i))
       if(get_err().ne.0) then
          begin_err(1)
          write(0,'("error at cloning ", I2, " th element")') i
          end_err()
       end if
       this%keys(i) = other%keys(i)
    end do
  end subroutine object_clone
  recursive subroutine object_delete(this)
    type(object) :: this
    integer :: i
    
    do i = 1, this % size
       call value_delete(this % vals(i)); check_err()
       !write(*,*) i, trim(this % keys(i))
       !this % keys(i) = ""       
    end do
    this % size = 0
  end subroutine object_delete
  ! -- Accessor --
  subroutine object_set(o, k, v)
    type(object) o
    character(*) k
    type(value) v
    integer i

    if(o%size < 0) then
       begin_err(1)
       call err_1("wrong state for object o")
       call err_si("size: ", o%size)
       end_err()
    end if
    
    i = o % size + 1
    o % keys(i) = trim(k)
    call value_clone(o%vals(i), v); check_err()
    o % size = i
  end subroutine object_set
  subroutine object_get(o, k, v)
    type(object) o    
    character(*) k
    type(value) v
    integer i

    if(size(o%keys) .eq. 0) then
       throw_err("no lement", 1)
    end if
    
    do i = 1, o%size
       if(o%keys(i) .eq. k) then
          !v = o%vals(i)
          call value_clone(v, o%vals(i))
          return
       end if
    end do
    begin_err(1)
    call err_1("failed to find key")
    call err_ss("key: ", k)
    end_err()
  end subroutine object_get
  function object_exist(o, k) result(res)
    type(object) o    
    character(*) k
    logical res
    integer i

    do i = 1, o%size
       if(o%keys(i) .eq. k) then
          res = .true.
          return
       end if
    end do
    res = .false.
    return

  end function object_exist
  subroutine object_set_i(o, k, i)
    type(object) o    
    character(*) k
    integer i

    type(value) v
    call value_new_i(v, i)
    call object_set(o, k, v); check_err()
  end subroutine object_set_i
  subroutine object_set_d(o, k, x)
    type(object) o    
    character(*) k
    double precision x

    type(value) v
    call value_new_d(v, x)
    call object_set(o, k, v); check_err()
  end subroutine object_set_d
  subroutine object_set_b(o, k, x)
    type(object) o    
    character(*) k
    logical x

    type(value) v
    call value_new_b(v, x)
    call object_set(o, k, v); check_err()
  end subroutine object_set_b
  subroutine object_set_s(o, k, x)
    type(object) o    
    character(*) k
    character(*) x

    type(value) v
    call value_new_s(v, x); check_err()
    call object_set(o, k, v); check_err()
    call value_delete(v); check_err()
  end subroutine object_set_s
  subroutine object_set_a(o, k, x)
    type(object) o    
    character(*) k
    type(array) x
    type(value) v
    call value_new_a(v, x); check_err()
    call object_set(o, k, v); check_err()
    call value_delete(v); check_err()
  end subroutine object_set_a
  subroutine object_set_o(o, k, x)
    type(object) o    
    character(*) k
    type(object) x
    type(value) v
    call value_new_o(v, x); check_err()
    call object_set(o, k, v); check_err()
    call value_delete(v); check_err()
  end subroutine object_set_o
  subroutine object_get_i(o, k, x)
    type(object) o    
    character(*) k
    integer x

    type(value) v
    call object_get(o, k, v); check_err()
    call value_get_i(v, x); check_err()
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       end_err()
    end if
    call value_delete(v); check_err()
    
  end subroutine object_get_i
  subroutine object_get_d(o, k, x)
    type(object) o    
    character(*) k
    double precision x
    type(value) v
    
    call object_get(o, k, v); check_err()
    call value_get_d(v, x)
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       return
    end if
    call value_delete(v); check_err()
    
  end subroutine object_get_d
  subroutine object_get_b(o, k, x)
    type(object) o    
    character(*) k
    logical :: x
    type(value) v
    
    call object_get(o, k, v); check_err()
    call value_get_b(v, x); check_err()
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       end_err()
    end if
    call value_delete(v); check_err()

  end subroutine object_get_b
  subroutine object_get_s(o, k, x)
    type(object) o    
    character(*) k
    character(*) x

    type(value) v
    call object_get(o, k, v); check_err()
    call value_get_s(v, x); check_err()
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       end_err()
    end if
    call value_delete(v); check_err()
    
  end subroutine object_get_s
  subroutine object_get_a(o, k, x)
    type(object) o    
    character(*) k
    type(array) x
    type(value) v
    call object_get(o, k, v); check_err()
    call value_get_a(v, x); check_err()
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       end_err()
    end if
    call value_delete(v); check_err()

  end subroutine object_get_a
  subroutine object_get_o(o, k, x)
    type(object) o    
    character(*) k
    type(object) x

    type(value) v
    call object_get(o, k, v); check_err()
    call value_get_o(v, x); check_err()
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on conversion")
       call err_ss("key: ", k)
       end_err()
    end if
    call value_delete(v); check_err()

  end subroutine object_get_o
  subroutine object_get_type(o, k, type)
    type(object) o    
    character(*) k
    type(value) v
    integer type, i
    
    do i = 1, o%size
       call object_get(o, k, v); check_err()
       type = v % type
       call value_delete(v); check_err()
    end do
  end subroutine object_get_type
  subroutine object_get_s_default(o, k, s0, s)
    type(object) :: o
    character(*), intent(in) :: k, s0
    character(*), intent(out) :: s
    
    if(object_exist(o, k)) then
       call object_get_s(o, k, s) ; check_err()       
    else
       s = s0
    end if
    
  end subroutine object_get_s_default
  ! -- IO --
  recursive subroutine object_dump(o, in_ifile)

    type(object) o
    integer, optional :: in_ifile
    integer i, ifile
    
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    depth = depth + 1
    
    write(ifile,fmt='(A)') "{"
    do i=1, o%size
       call insert_depth(ifile)
       write(ifile,fmt='(A,A,A)',advance='no') '"', trim(o%keys(i)), '":'
       select case(o % vals(i) % type)
       case(TYPE_O)
          call object_dump(o%vals(i)%val_o, ifile); check_err()
       case default
          call value_dump(o % vals(i), ifile); check_err()
       end select
       if(i.ne.o%size) then
          write(ifile, *) ","
       end if
    end do
    depth = depth - 1
    write(ifile, *) ""
    call insert_depth(ifile)
    write(ifile,fmt='(A)',advance='no') "}"    
    
  end subroutine object_dump

  ! ==== Array ====
  ! -- Constructor --
  subroutine array_new(this, in_size)
    type(array) :: this
    integer     :: in_size
    this % size = in_size
  end subroutine array_new
  subroutine array_clone(this, other)
    type(array) this, other
    integer i
    this%size = other%size
    do i = 1, this%size
       call value_clone(this%vals(i), other%vals(i)); check_err()
    end do
    
  end subroutine array_clone
  recursive subroutine array_delete(this)
    type(array) :: this
    integer i
    do i = 1, this % size
       call value_delete(this % vals(i)); check_err()
    end do
    this % size = 0
  end subroutine array_delete
  ! -- Accessor --
  subroutine array_get_size(this, out_size)
    type(array) :: this
    integer      :: out_size
    out_size = this % size
  end subroutine array_get_size
  function array_size(this)
    type(array) :: this
    integer :: array_size
    array_size = this % size
  end function array_size
  subroutine array_set(this, i, v)
    type(array) :: this
    integer     :: i
    type(value) :: v

    if(i < 1 .or. this % size < i) then
       begin_err(1)
       call err_1("index out of range")
       call err_si("i: ", i)
       call err_si("size: ", this % size)
       return
    end if
    call value_clone(this%vals(i), v); check_err()
  end subroutine array_set
  subroutine array_get(this, i, v)
    type(array) :: this
    integer     :: i
    type(value) :: v

    if(i < 1 .or. this % size < i) then
       begin_err(1)
       call err_1("index out of range")
       call err_si("i: ", i)
       call err_si("size: ", this % size)
       return
    end if
    !v = this % vals(i)
    call value_clone(v, this%vals(i)); check_err()
  end subroutine array_get
  subroutine array_add(this, v)
    type(array) :: this
    type(value) :: v

    this % size = this % size + 1
    !this % vals(this % size) = v
    call array_set(this, this%size, v)
  end subroutine array_add
  subroutine array_add_i(this, x)
    type(array) :: this
    type(value) :: v
    integer x    
    call value_new_i(v, x); check_err()
    call array_add(this, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_i
  subroutine array_add_d(this, x)
    type(array) :: this
    type(value) :: v
    double precision x    
    call value_new_d(v, x); check_err()
    call array_add(this, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_d
  subroutine array_add_s(this, x)
    type(array) :: this
    type(value) :: v
    character(*) :: x    
    call value_new_s(v, x); check_err()
    call array_add(this, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_s
  subroutine array_add_a(this, x)
    type(array) :: this, x
    type(value) :: v
    call value_new_a(v, x)
    call array_add(this, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_a
  subroutine array_add_o(this, x)
    type(array) :: this
    type(object) :: x
    type(value) :: v
    call value_new_o(v, x)
    call array_add(this, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_o
  subroutine array_add_as_c(this, x)
    type(array) :: this
    complex(kind(0d0)) :: x
    type(value) :: v
    type(object) :: o
    call object_set_d(o, "re", real(x)); check_err()
    call object_set_d(o, "im", aimag(x)); check_err()
    call value_new_o(v, o); check_err()
    call array_add(this, v); check_err()
    call object_delete(o); check_err()
    call value_delete(v); check_err()
  end subroutine array_add_as_c
  subroutine array_set_i(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    integer x    
    call value_new_i(v, x)
    call array_set(this, i, v); check_err()
  end subroutine array_set_i
  subroutine array_set_d(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    double precision x

    call value_new_d(v, x)
    call array_set(this, i, v); check_err()
  end subroutine array_set_d
  subroutine array_set_s(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    character(*) x

    call value_new_s(v, x)
    call array_set(this, i, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_set_s
  subroutine array_set_a(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    type(array) :: x
    call value_new_a(v, x); check_err()
    call array_set(this, i, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_set_a
  subroutine array_set_o(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    type(object) :: x
    call value_new_o(v, x)
    call array_set(this, i, v); check_err()
    call value_delete(v); check_err()
  end subroutine array_set_o
  subroutine array_add_as_darray(this, xs)
    type(array) :: this
    double precision :: xs(:)
    integer i
    do i = 1, size(xs)
       call array_add_d(this, xs(i)); check_err()
    end do
  end subroutine array_add_as_darray
  subroutine array_get_i(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    integer x

    call array_get(this, i, v)
    call value_get_i(v, x); check_err()
  end subroutine array_get_i
  subroutine array_get_d(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    double precision x
    
    call array_get(this, i, v)
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on getting element of array")
       call err_si("i: ", i)
       end_err()
    end if
    call value_get_d(v, x)
    if(ierr .ne. 0) then
       begin_err(1)
       call err_1("error on extracting double")
       call err_si("i:", i)
       end_err()
    end if
    
  end subroutine array_get_d
  subroutine array_get_s(this, i, x)
    type(array) :: this
    integer i
    type(value) :: v
    character(*) x

    call array_get(this, i, v) ; check_err()
    call value_get_s(v, x); check_err()
    call value_delete(v); check_err()
    
  end subroutine array_get_s
  subroutine array_get_o(this, i, x)
    type(array) :: this
    integer i
    type(object) :: x
    type(value) :: v    
    call array_get(this, i, v); check_err()
    call value_get_o(v, x); check_err()
    call value_delete(v); check_err()
  end subroutine array_get_o
  subroutine array_get_a(this, i, x)
    type(array) :: this
    integer i
    type(array) :: x
    type(value) :: v    
    call array_get(this, i, v); check_err()
    call value_get_a(v, x); check_err()
    call value_delete(v); check_err()
  end subroutine array_get_a
  subroutine array_get_as_iarray(this, xs)
    type(array) :: this
    integer     :: xs(:)
    integer     :: i, x
    type(value) :: v

    if(size(xs) < this % size) then
       begin_err(1)
       call err_1("xs is too small")
       call err_si("this % size: ", this % size)
       call err_si("size(xs):    ", size(xs))
       return 
    end if

    do i = 1, this % size
       call array_get(this, i, v); check_err()
       call value_get_i(v, x); check_err()
       xs(i) = x
    end do
    
  end subroutine array_get_as_iarray
  subroutine array_get_as_darray(this, xs)
    type(array) :: this
    double precision      :: xs(:)
    integer     :: i
    double precision      :: x
    type(value) :: v

    if(size(xs) < this % size) then
       begin_err(1)
       call err_1("xs is too small")
       call err_si("this % size: ", this % size)
       call err_si("size(xs):    ", size(xs))
       end_err()
    end if

    do i = 1, this % size
       call array_get(this, i, v); check_err()
       call value_get_d(v, x); check_err()
       xs(i) = x
    end do
    
  end subroutine array_get_as_darray
  subroutine array_get_as_carray(this, xs)
    type(array) :: this
    complex(kind(0d0)) :: xs(:)
    integer     :: i
    complex(kind(0d0)) :: x
    type(value) :: v

    if(size(xs) < this % size) then
       begin_err(1)
       call err_1("xs is too small")
       call err_si("this % size: ", this % size)
       call err_si("size(xs):    ", size(xs))
       end_err()
    end if

    do i = 1, this % size
       call array_get(this, i, v); check_err()
       call value_get_as_c(v, x); check_err()
       xs(i) = x
    end do
  end subroutine array_get_as_carray
  ! -- IO --
  recursive subroutine array_dump(this, ifile)
    type(array) :: this
    integer ifile
    integer i

    depth = depth + 1
    write(ifile, '(A)') "["
    do i = 1, this % size
       call insert_depth(ifile)
       call value_dump(this % vals(i), ifile); check_err()
       if(i .eq. this % size) then
          write(ifile,*) ""
       else
          write(ifile, '(a)')  ","
       end if
    end do
    depth = depth - 1
    call insert_depth(ifile)
    write(ifile, '(A)', advance='no') "]"
       

    ! -- no new line version --
!    write(ifile, '(A)', advance='no') "["
!    do i = 1, this % size
!       call value_dump(this % vals(i), ifile); check_err()
!       if(i .ne. this%size) then
!          write(ifile, '(a)', advance='no') ", "
!       end if
!    end do
!    write(ifile, '(A)', advance='no') "]"
  end subroutine array_dump
  subroutine array_dump_debug(this, ifile)
    type(array) :: this
    integer ifile
    integer i
    write(ifile, '(A)') "array object:"
    do i = 1, this % size
       write(ifile, *) i, ":"
       call value_dump_debug(this % vals(i), ifile); check_err()
    end do
    write(ifile, '(A)', advance='no') "]"
  end subroutine array_dump_debug
end module mod_fjson_value

module mod_fjson_loader
  !
  ! analyse parsed array and build json value.
  !
  ! example for json_load_obj
  ! in_ie, in_ie+1, in_ie+2, in_ie+3, in_ie+4, in_ie+5=out_ie
  !   {  ,  "xkey",   :    ,    5   ,      } ,  ????
  use mod_err_handle
  use mod_fjson_parser
  use mod_fjson_value
  implicit none
contains
  subroutine loader_error(ie)
    integer ie, iie
    character(100) line
    call err_si("ie: ", ie)

    do iie = max(1, ie-5), ie+5
       if(tokens(ie-1)%type .ne. TOKEN_NULL) then
          write(line, '("eles[", I0, "]: ", A)') iie, trim(tokens(iie)%val)
       else
          write(line, '("eles[", I0, "]: null")') iie
       end if
       call err_1(line)
    end do
    
  end subroutine loader_error
  subroutine loads_json(stream, v)
    type(istream) :: stream
    type(value)   :: v
    integer  :: out_ie
    
#if defined DEBUG_FJSON
    write(*, *) "start loader"
#endif
#if defined DEBUG_FJSON
    write(*, *) "start parse"
#endif        
    call parser_new; check_err()
    call parser_parse(stream); check_err()

#if defined DEBUG_FJSON
    write(*, *) "parse results:"
    call parser_dump(6); check_err()
#endif            
    call load_value(1, out_ie, v); check_err()

#if defined DEBUG_FJSON
    write(*, *) "end loads_json:"
    call parser_dump(6); check_err()
#endif                
  end subroutine loads_json
  subroutine loads_json_file(filename, ifile, v)
    character(*) filename
    integer ifile
    type(value) :: v
    type(istream) :: stream

    call istream_new_file(stream, ifile, filename); check_err()
    call loads_json(stream, v); check_err()
    call istream_delete(stream); check_err()
    
  end subroutine loads_json_file
  subroutine loads_json_string(str, v)
    character(*) str
    type(value) :: v
    type(istream) :: stream

    call istream_new_string(stream, str); check_err()
    call loads_json(stream, v); check_err()
    call istream_delete(stream); check_err()

  end subroutine loads_json_string
  recursive subroutine load_value(in_ie, out_ie, v)
    integer :: in_ie, out_ie
    type(value) :: v
    integer :: ie, ie_out_tmp

#if defined DEBUG_FJSON
    write(*, *) "start load_value"
#endif    
    
    ie = in_ie
    select case(tokens(ie) % type)
    case (TOKEN_NUM)
       call load_num(ie, ie_out_tmp, v); check_err()
    case (TOKEN_STR)
       call load_str(ie, ie_out_tmp, v); check_err()
    case (TOKEN_BOOL)
       call load_bool(ie, ie_out_tmp, v); check_err()
    case (TOKEN_BEGIN_OBJ)
       call load_obj(ie, ie_out_tmp, v); check_err()
    case (TOKEN_BEGIN_ARY)
       call load_ary(ie, ie_out_tmp, v); check_err()       
    case default
       begin_err(1)
       call err_1("unsupported vlaue type")
       call loader_error(ie)
       end_err()
    end select
    if(ie_out_tmp>1000) then
       begin_err(1)
       write(0,*) "out_ie out of range"
       write(0,*) "token:", tokens(ie) % type
       write(0,*) "out_ie:", ie_out_tmp
       end_err()
    end if
    out_ie = ie_out_tmp

#if defined DEBUG_FJSON
    write(*, '("  in_ie, out_ie = ",i0, ", ", i0)') in_ie, out_ie
    write(*, *) "end load_value"
#endif            
    
  end subroutine load_value
  recursive subroutine load_num(in_ie, out_ie, v)
    use Mod_math, only : is_i, is_d, convert_i, convert_d
    integer in_ie, out_ie, ie
    type(value) v
    logical res_i, res_d    
    integer val_i
    double precision  val_d

#if defined DEBUG_FJSON
    write(*, *) "start load_num"
#endif    
    
    ie = in_ie

    if(tokens(ie) % type .ne. TOKEN_NUM) then
       begin_err(1)
       call err_1("type mismatch.")
       call loader_error(ie)
       end_err()
    end if
    
    call is_i(tokens(ie) % val, res_i)
    call is_d(tokens(ie) % val, res_d)

    if(res_i) then
       call convert_i(tokens(ie)%val, val_i)
       call value_new_i(v, val_i)
    else if(res_d) then
       call convert_d(tokens(ie)%val, val_d)
       call value_new_d(v, val_d)
    else
       throw_err("conversion to number failed", 1)
       return
    end if
    ie = ie + 1
    out_ie = ie

#if defined DEBUG_FJSON
    write(*, '("  in_ie, out_ie, v = ",i0, ", ", i0,", ")', advance='no') in_ie, out_ie
    call value_dump(v, 6)
    write(*, *)
    write(*, *) "end load_num"
#endif        
        
  end subroutine load_num
  recursive subroutine load_str(in_ie, out_ie, v)
    integer in_ie, out_ie, ie
    type(value) v

#if defined DEBUG_FJSON
    write(*, *) "start load_str"
#endif    
    
    ie = in_ie
    if(tokens(ie) % type .ne. TOKEN_STR) then
       begin_err(1)
       call err_1("type must be str")       
       call loader_error(ie)
       end_err()
    end if

    call value_new_s(v, tokens(ie) % val); check_err()
    ie = ie + 1
    out_ie = ie

#if defined DEBUG_FJSON
    write(*, '("  in_ie, out_ie, v = ",i0, ", ", i0,", ")', advance='no') in_ie, out_ie
    call value_dump(v, 6)
    write(*, *)
    write(*, *) "end load_str"
#endif        
    
  end subroutine load_str
  recursive subroutine load_bool(in_ie, out_ie, v)
    integer ie, in_ie, out_ie
    type(value) v

#if defined DEBUG_FJSON
    write(*, *) "start load_bool"
#endif    
    
    ie = in_ie
    if(tokens(ie) % type .ne. TOKEN_BOOL) then
       begin_err(1)
       call err_1("type must be bool")
       call loader_error(ie)
       end_err()
    end if

    if(tokens(ie) % val .eq. "true") then
       call value_new_b(v, .true.)
    else if(tokens(ie) % val .eq. "false") then
       call value_new_b(v, .false.)
    else
       begin_err(1)
       call err_1("invalid value")
       call err_ss("val: ", tokens(ie) % val)
       end_err()
    end if
    
    out_ie = ie+1
    
#if defined DEBUG_FJSON
    write(*, '("  in_ie, out_ie, v = ",i0, ", ", i0,", ")', advance='no') in_ie, out_ie
    call value_dump(v, 6)
    write(*, *)
    write(*, *) "end load_bool"
#endif            
  end subroutine load_bool
  recursive subroutine load_obj(in_ie, out_ie, v)
    integer     :: in_ie, out_ie, ie, out_ie_tmp
    type(value) :: v, v_val, vkey
    type(object) :: o
    character(100) key

#if defined DEBUG_FJSON
    write(*, *) "start load_obj"
#endif    

    ! -- newialize --
    call object_new(o)
    ie = in_ie

    ! -- skip { --
    call load_skip_one(ie, out_ie_tmp, TOKEN_BEGIN_OBJ); check_err()
    ie = out_ie_tmp
    
    do
       ! -- get key --
       call load_str(ie, out_ie_tmp, vkey); check_err()
       call value_get_s(vkey, key); check_err()
       call value_delete(vkey); check_err()
       !key = vkey % val_s
       ie = out_ie_tmp 

       ! -- skip separator --
       call load_skip_one(ie, out_ie_tmp, TOKEN_SEP); check_err()
       ie = out_ie_tmp
       
       ! -- get value --
       call load_value(ie, out_ie_tmp, v_val); check_err()
       call object_set(o, key, v_val); check_err()
       call value_delete(v_val); check_err()
       ie = out_ie_tmp

       ! -- check final element or not
!       write(*,*) "ie:", ie
!       write(*,*) "tok:", tokens(ie)%type
       if(tokens(ie)%type .eq. TOKEN_COMMA) then
          call load_skip_one(ie, out_ie_tmp, TOKEN_COMMA); check_err()
          ie = out_ie_tmp
          cycle
       else if(tokens(ie)%type .eq. TOKEN_END_OBJ) then
          call load_skip_one(ie, out_ie_tmp, TOKEN_END_OBJ); check_err()
          ie = out_ie_tmp
          exit
       else
          begin_err(1)
          call err_1("failed to find comma or }")
          call loader_error(ie)
          end_err()
       end if
    end do
    
    ! -- convert to value --
    call value_new_o(v, o)
    out_ie = ie

    ! -- finalize --        
    call object_delete(o); check_err()

#if defined DEBUG_FJSON
    write(*, *) "  in_ie, out_ie ", in_ie, out_ie
    write(*, *) "end load_obj"
#endif        
    
  end subroutine load_obj
  recursive subroutine load_ary(in_ie, out_ie, v)
    integer :: in_ie, out_ie, ie, out_ie_tmp
    type(value) :: v, v_val
    type(array) :: a

#if defined DEBUG_FJSON
    write(*, *) "start load_ary"
#endif    

    ! -- newialize --
    call array_new(a, 0)
    ie = in_ie

    ! -- skip [ --
    call load_skip_one(ie, out_ie_tmp, TOKEN_BEGIN_ARY); check_err()
    ie = out_ie_tmp
    
    do
       ! -- get value --
       call load_value(ie, out_ie_tmp, v_val); check_err()
       call array_add(a, v_val); check_err()
       call value_delete(v_val); check_err()
       ie = out_ie_tmp

       ! -- check final or not --
       if(tokens(ie) % type .eq. TOKEN_COMMA) then
          call load_skip_one(ie, out_ie_tmp, TOKEN_COMMA); check_err()
          ie = out_ie_tmp
          cycle
       else if(tokens(ie)%type .eq. TOKEN_END_ARY) then
          call load_skip_one(ie, out_ie_tmp, TOKEN_END_ARY); check_err()
          ie = out_ie_tmp
          exit
       else
          begin_err(1)
          call err_1("failed to find comma or ]")
          call loader_error(ie)
          end_err()
       end if
    end do

    ! -- convert to value --
    call value_new_a(v, a); check_err()
    out_ie = ie

    ! -- finalize --    
    call array_delete(a)

#if defined DEBUG_FJSON
    write(*, *) "  in_ie, out_ie ", in_ie, out_ie
    write(*, *) "end load_a"
#endif            
    
  end subroutine load_ary
  recursive subroutine load_skip_one(in_ie, out_ie, ele_type)
    integer in_ie, out_ie, ele_type, ie
#if defined DEBUG_FJSON
    write(*, *) "start load_skip_one: ", type_names(tokens(in_ie) % type)
#endif

    ie = in_ie
    if(tokens(in_ie) % type .ne. ele_type) then
       begin_err(1)
       call err_1("type mismatch")
       call err_ss("desired type: ", type_names(ele_type))
       call loader_error(ie)
       return
    end if
    out_ie = in_ie + 1
    
#if defined DEBUG_FJSON
    write(*, '("  in_ie, out_ie, val = ",i0, ", ", i0, ", ", a)', advance='no') in_ie, out_ie, trim(tokens(in_ie)%val)
    write(*, *)
    write(*, *) "end load_skip_one"
#endif        
    
  end subroutine load_skip_one
  
end module mod_fjson_loader

module mod_fjson_convert
  use mod_fjson_value
  implicit none
contains
  recursive subroutine c2v(c, v)
    complex(kind(0d0)) :: c
    type(value) :: v
    type(object) :: o
    
    if(v % type .ne. TYPE_VALUE_NULL) then
       call value_delete(v); check_err()
    end if

    call object_set_d(o, "re", real(c))
    call object_set_d(o, "im", aimag(c))
    call value_new_o(v, o); check_err()

    call object_delete(o); check_err()
    
  end subroutine c2v
  subroutine v2c(v, c)
    type(value) :: v
    complex(kind(0d0)) :: c    
    type(object) :: o
    double precision       :: re, im

    if(v % type .eq. TYPE_O) then
       call value_get_o(v, o); check_err()
       call object_get_d(o, "re", re); check_err()
       call object_get_d(o, "im", im); check_err()
       c = cmplx(re, im, kind(0d0))
    else if(v % type .eq. TYPE_D) then
       call value_get_d(v, re); check_err()
       c = re
    else
       begin_err(1)
       call err_1("invalid type")
       if(show_message_q) then
          write(*, '(a)', advance='no') "v = "
          call value_dump(v, 6)
          write(*, *)
       end if
       end_err()
    end if

    call object_delete(o); check_err()
    
  end subroutine v2c
  subroutine v2d(v, d)
    type(value) :: v
    double precision d
    integer i
    if(v % type .eq. TYPE_D) then
       call value_get_d(v, d)
    else if(v % type .eq. TYPE_I) then
       call value_get_i(v, i)
       d = i
    else
    end if
  end subroutine v2d
  subroutine ivec2a(xs, a)
    integer     :: xs(:)
    type(array) :: a
    integer :: num, i

    num = size(xs)

    if(a%size.ne.num) then
       call array_new(a, num); check_err()
    end if
    
    do i = 1, num
       call array_set_i(a, i, xs(i))
    end do
    
  end subroutine ivec2a
  subroutine dvec2a(xs, a)
    double precision      :: xs(:)
    type(array) :: a
    integer :: num, i

    num = size(xs)

    if(a%size.ne.num) then
       call array_new(a, num); check_err()
    end if
    
    do i = 1, num
       call array_set_d(a, i, xs(i))
    end do
    
  end subroutine dvec2a
  subroutine dmat2a(xss, a)
    double precision      :: xss(:,:)
    type(array) :: a
    type(array) :: ai

    integer :: numi, numj, i, j

    call array_new(a, 0); check_err()
    numi = size(xss, 1)
    numj = size(xss, 2)
    do i = 1, numi
       call array_new(ai, 0); check_err()
       do j = 1, numj
          call array_add_d(ai, xss(i, j))          
       end do
       call array_add_a(a, ai); check_err()
       call array_delete(ai); check_err()
    end do
  end subroutine dmat2a
  subroutine a2dvec(a, xs)
    type(array) :: a
    double precision :: xs(:)
    integer     ::  num, i

    call array_get_size(a, num)

    if(size(xs).ne.num) then
       begin_err(1)
       write(0,*) "size mismatch"
       write(0,*) "size(a):", num
       write(0,*) "size(xs):", size(xs)
       end_err()
    end if
    
    do i = 1, num
       call array_get_d(a, i, xs(i)); check_err()
    end do
    
  end subroutine a2dvec
  subroutine a2cvec(a, xs)
    type(array) :: a
    complex(kind(0d0)) :: xs(:)
    integer num, i
    type(value) :: v
    
    num = array_size(a); check_err()
    if(num < size(xs)) then
       throw_err("size of xs is not enough", 1)
    end if

    do i = 1, num
       call array_get(a, i, v); check_err()
       call v2c(v, xs(i)); check_err()
       call value_delete(v); check_err()
    end do
    
  end subroutine a2cvec
  subroutine a2dmat(a, xss)
    type(array) :: a
    double precision :: xss(:,:)
    integer     ::  numi, numj, i, j, numtmp
    type(value) :: vi
    type(array) :: ai
    double precision :: d
    
    call array_get_size(a, numi)
    call array_get_a(a, 1, ai)
    call array_get_size(ai, numj)
    if(numi.ne.size(xss,1).or.numj.ne.size(xss,2)) then
       begin_err(1)
       write(0,*) "size mismatch"
       write(0,*) "size(xss):", size(xss,1), size(xss,2)
       write(0,*) "size(a)  :", numi, numj
       end_err()
    end if
    
    do i = 1, numi
       call array_get(a, i, vi); check_err()
       call value_get_a(vi, ai);  check_err()
       call array_get_size(ai, numtmp); check_err()
       if(numj.ne.numtmp) then
          throw_err("size mismatch", 1)
       end if
       do j = 1, numj
          call array_get_d(ai, j, d); check_err()
          xss(i, j) = d
       end do
       call array_delete(ai); check_err()
    end do
    
  end subroutine a2dmat
  subroutine a2cmat(a, xss)
    type(array) :: a
    complex(kind(0d0)) :: xss(:,:)
    integer     ::  numi, numj, i, j
    type(array) :: ai
    
    call array_get_size(a, numi)
    call array_get_a(a, 1, ai)
    call array_get_size(ai, numj)
    call array_delete(ai)
    
    if(numi.ne.size(xss,1).or.numj.ne.size(xss,2)) then
       begin_err(1)
       write(0,*) "size mismatch"
       write(0,*) "size(xss):", size(xss,1), size(xss,2)
       write(0,*) "size(a)  :", numi, numj
       end_err()
    end if
    
    do i = 1, numi
       do j = 1, numj
          call value_get_as_c(a%vals(i)%val_a%vals(j), xss(i,j))
       end do
       !call array_get(a, i, vi); check_err()
       !call value_get_a(vi, ai);  check_err()
       !call array_get_size(ai, numtmp); check_err()
       !if(numj.ne.numtmp) then
       !   throw_err("size mismatch", 1)
       !end if
       !call array_get_as_carray(ai, xss(i,:)); check_err()
       !call array_delete(ai); check_err()
    end do
    
  end subroutine a2cmat
  subroutine cvec2a(xs, a)
    complex(kind(0d0)) :: xs(:)
    type(array) :: a
    type(value) :: v

    integer :: num, i
     
    num = size(xs)
    call array_new(a, num); check_err()
    
    do i = 1, num
       call c2v(xs(i), v); check_err()
       call array_set(a, i, v); check_err()
       call value_delete(v); check_err()
    end do
    
  end subroutine cvec2a
  subroutine cmat2a(xss, a)
    complex(kind(0d0)) :: xss(:,:)
    type(array) :: a
    type(array) :: ai
    type(value) :: v
    integer :: numi, numj, i, j

    call array_new(a, 0)
    
    numi = size(xss, 1)
    numj = size(xss, 2)
    do i = 1, numi
       call array_new(ai, 0)
       do j = 1, numj
          call c2v(xss(i, j),  v)
          call array_add(ai, v)
          call value_delete(v); check_err()
       end do
       call array_add_a(a, ai)
       call array_delete(ai)
    end do
  end subroutine cmat2a
  subroutine dten2a(xsss, a)
    double precision :: xsss(:,:,:)
    type(array) :: a, aa
    integer :: numi, i

    call array_new(a, 0)
    
    numi = size(xsss, 1)
    do i = 1, numi
       call dmat2a(xsss(i,:,:), aa); check_err()
       call array_add_a(a, aa); check_err()
       call array_delete(aa); check_err()
    end do
  end subroutine dten2a
end module mod_fjson_convert

module mod_fjson  
  use mod_fjson_parser
  use mod_fjson_value
  use mod_fjson_loader
  use mod_fjson_convert
  implicit none
end module mod_fjson



