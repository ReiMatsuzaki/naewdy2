#include "macros.fpp"

!     Unit test utilities

!
! TODO
! more pretty message !
module mod_utest_check
  implicit none
  interface utest_check_near
     subroutine check_near_d(a, b, eps, a_str, b_str, file, line, stopq)
       double precision a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq       
     end subroutine check_near_d
     subroutine check_near_c(a, b, eps, a_str, b_str, file, line, stopq)
       complex(kind(0d0)) a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq
     end subroutine check_near_c
  end interface utest_check_near
  interface utest_check_eq
     subroutine check_eq_i(a, b, a_str, b_str, file, line, stopq)
       integer a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq
     end subroutine check_eq_i
     subroutine check_eq_d(a, b, a_str, b_str, file, line, stopq)
       double precision a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq
     end subroutine check_eq_d
     subroutine check_eq_c(a, b, a_str, b_str, file, line, stopq)
       complex(kind(0d0)) a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq
     end subroutine check_eq_c
     subroutine check_eq_s(a, b, a_str, b_str, file, line, stopq)
       character(*) a, b
       character(*) a_str, b_str, file
       integer line
       logical stopq
     end subroutine check_eq_s
  end interface utest_check_eq
  interface utest_check_near_0
     subroutine check_near_d_0(a, b, eps, file, line, stopq)
       double precision :: a, b, eps
       character(*) file
       integer line
       logical stopq
     end subroutine check_near_d_0
     subroutine check_near_c_0(a, b, eps, file, line, stopq)
       complex(kind(0d0)) a, b
       double precision :: eps
       character(*) file
       integer line
       logical stopq
     end subroutine check_near_c_0
  end interface utest_check_near_0
  interface utest_check_eq_0
     subroutine check_eq_i_0(a, b, file, line, stopq)
       integer a, b
       character(*) file
       integer line
       logical stopq
     end subroutine check_eq_i_0
     subroutine check_eq_d_0(a, b, file, line, stopq)
       double precision a, b
       character(*) file
       integer line
       logical stopq
     end subroutine check_eq_d_0
     subroutine check_eq_c_0(a, b, file, line, stopq)
       complex(kind(0d0)) a, b
       character(*) file
       integer line
       logical stopq
     end subroutine check_eq_c_0
     subroutine check_eq_s_0(a, b, file, line, stopq)
       character(*) a, b
       character(*) file
       integer line
       logical stopq
     end subroutine check_eq_s_0
  end interface utest_check_eq_0
end module mod_utest_check
module mod_utest
  implicit none
!  private
  integer, save, private :: num_utest
  integer, save, private :: num_failed
!  public:: utest_check_eq, utest_check_near  
contains 
  ! ==== constructors ====
  subroutine utest_begin

    num_utest = 0
    num_failed = 0
    write(*,*) 
    write(*,*) "=========================="
    write(*,*) "          UTtest          "
    write(*,*) "=========================="
    write(*,*) 
    
  end subroutine utest_begin
  subroutine utest_end
    write(*,*) "--------------------------"
    write(*,*) "          Results         "
    write(*,*) "--------------------------"
    if(num_failed .eq. 0) then
       write(*, '(I0, " tests passed")') num_utest
    else
       write(*,*) num_failed, num_utest
    end if

  end subroutine utest_end

  ! ==== utils ====
  subroutine utest_check_begin
    num_utest = num_utest + 1
  end subroutine utest_check_begin
  subroutine utest_check_fail(msg, file, line)
    character(*) msg, file
    integer line
    num_failed = num_failed + 1
    write(*, '(a, ":", i0, ": ", a)') file, line, msg
    return
  end subroutine utest_check_fail
  
end module mod_utest

! ==== main ====
subroutine check_near_d(a,b,eps,a_str,b_str,file,line,stopq)
  use mod_utest
    double precision, intent(in) :: a, b
    double precision, intent(in) :: eps
    character(*) :: a_str, b_str, file
    integer :: line
    logical :: stopq
    call utest_check_begin    
       if(abs(a-b) > eps) then
          call utest_check_fail("not equal error", file, line)
          write(*, '(a, ": ", 2E20.10)') a_str, a
          write(*, '(a, ": ", 2E20.10)') b_str, b
          write(*, '("eps: ", E20.5)') eps
          write(*, *)
       if(stopq) then
          stop
       end if
    end if       
  end subroutine check_near_d
  subroutine check_near_c(a,b,eps,a_str,b_str,file,line,stopq)
    use mod_utest
    complex(kind(0d0)), intent(in) :: a, b
    double precision, intent(in)   :: eps
    character(*) :: a_str, b_str, file
    integer :: line
    logical :: stopq

    call utest_check_begin    
    if(abs(a-b) > eps) then
       call utest_check_fail("not equal error", file, line)
       write(*, '(a, ": ", 2E20.10)') a_str, a
       write(*, '(a, ": ", 2E20.10)') b_str, b
       write(*, '("|", a, "-", a, "| = ", E20.10)') a_str, b_str, abs(a-b)
       write(*, '("eps: ", E20.5)') eps
       write(*, *)
       if(stopq) then
          stop
       end if
    end if
  end subroutine check_near_c
  subroutine check_eq_i(a, b, a_str, b_str, file, line, stopq)
    use mod_utest
    integer a, b
    character(*) a_str, b_str, file
    integer line
    logical stopq

    call utest_check_begin
    if(a .ne. b) then
       call utest_check_fail("not equal error", file, line)
       write(*, '(a, ": ", I0)') a_str, a
       write(*, '(a, ": ", I0)') b_str, b
       write(*, *)
       if(stopq) then
          stop
       end if
    end if
    return
  end subroutine check_eq_i
  subroutine check_eq_d(a, b, a_str, b_str, file, line, stopq)
    use mod_utest    
    double precision a, b
    character(*) a_str, b_str, file
    integer line
    logical stopq
    double precision, parameter :: eps = 2.0D-7
  
    call utest_check_begin
  
    if(abs(a-b) > eps) then
       call utest_check_fail("not equal error", file, line)
       write(*, '(a, ": ", E20.10)') a_str, a
       write(*, '(a, ": ", E20.10)') b_str, b
       write(*, *)
       if(stopq) then
          stop
       end if
    end if
  end subroutine check_eq_d
  subroutine check_eq_c(a, b, a_str, b_str, file, line, stopq)
    use mod_utest
    complex(kind(0d0)), intent(in) :: a, b
    character(*) :: a_str, b_str, file
    integer :: line
    logical :: stopq
    double precision, parameter :: eps = 1.0d-14
    call check_near_c(a, b, eps, a_str, b_str, file, line, stopq)
    
  end subroutine check_eq_c
  subroutine check_eq_s(a, b, a_str, b_str, file, line, stopq)
    use mod_utest
    character(*) a
    character(*) b
    character(*) a_str, b_str, file
    integer line
    logical stopq
  
    call utest_check_begin

    if(a .ne. b) then     
       num_failed = num_failed + 1
       write(*, '(A,":", I0, ":", A)') file, line, "equal test failed"  
       write(*, '(a, ": ", a)') a_str, a
       write(*, '(a, ": ", a)') b_str, b
       write(*, *)
       if(stopq) then
          stop
       end if
    end if
  end subroutine check_eq_s
  subroutine check_eq_char(a, b, a_str, b_str, file, line, stopq)
    use mod_utest
    character a
    character b
    character(*) a_str, b_str, file
    integer line
    logical stopq
    
    call utest_check_begin

    if(a .ne. b) then     
       num_failed = num_failed + 1
       write(*, '(A,":", I0, ":", A)') file, line, "equal test failed"  
       write(*, '(a, ": ", a)') a_str, a
       write(*, '(a, ": ", a)') b_str, b
       write(*, *)
       if(stopq) then
          stop
       end if
    end if
  end subroutine check_eq_char
  subroutine utest_check_true(a, a_str, file, line, stopq)
    use mod_utest
    logical a
    character(*) a_str, file
    integer line
    logical stopq
    
    call utest_check_begin
    
    if(.not. a) then
       call utest_check_fail("not true", file, line)
       write(*, '(a, ": false")') a_str
       if(stopq) then
          stop
       end if
    end if
    return
  end subroutine utest_check_true
  subroutine utest_check_false(a, a_str, file, line, stopq)
    use mod_utest
    logical a
    character(*) a_str, file
    integer line
    logical stopq
    
    call utest_check_begin
    
    if(a) then
       call utest_check_fail("not false", file, line)
       write(*, '(a, ": true")') a_str
       if(stopq) then
          stop
       end if
    end if
    return
  end subroutine utest_check_false

  ! ==== simpler version ====
  subroutine check_near_d_0(a, b, eps, file, line, stopq)
    use mod_utest
    double precision  ::  a, b
    double precision   :: eps
    character(*) file
    integer line
    logical stopq
    call check_near_d(a, b, eps, "left ", "right", file, line, stopq)    
  end subroutine check_near_d_0
  subroutine check_near_c_0(a, b, eps, file, line, stopq)
    use mod_utest
    complex(kind(0d0)) ::  a, b
    double precision   :: eps
    character(*) file
    integer line
    logical stopq
    call check_near_c(a, b, eps, "left ", "right", file, line, stopq)    
  end subroutine check_near_c_0
  subroutine check_eq_i_0(a, b, file, line, stopq)
    use mod_utest
    integer a, b
    integer line
    character(*) file
    logical stopq
    call check_eq_i(a, b, "left ", "right", file, line, stopq)
  end subroutine check_eq_i_0
  subroutine check_eq_d_0(a, b, file, line, stopq)
    use mod_utest
    double precision a, b
    character(*) file
    integer line
    logical stopq
    call check_eq_d(a, b, "left ", "right", file, line, stopq)
    
  end subroutine check_eq_d_0
  subroutine check_eq_c_0(a, b, file, line, stopq)
    use mod_utest
    complex(kind(0d0)) ::  a, b
    character(*) file
    integer line
    logical stopq
    call check_eq_c(a, b, "left ", "right", file, line, stopq)    
  end subroutine check_eq_c_0
  subroutine check_eq_char_0(a, b, file, line, stopq)
    use mod_utest
    character a
    character b
    character(*) file
    integer line
    logical stopq
    call check_eq_char(a, b, "left ", "right", file, line, stopq)

  end subroutine check_eq_char_0
  subroutine check_eq_s_0(a, b, file, line, stopq)
    use mod_utest
    character(*) a
    character(*) b
    character(*) file
    integer line
    logical stopq

    call check_eq_s(a, b, "left ", "right", file, line, stopq)

  end subroutine check_eq_s_0
  subroutine utest_check_true_0(a, file, line, stopq)
    use mod_utest
    logical a
    character(*) file
    integer line
    logical stopq
    
    call utest_check_true(a, "a", file, line, stopq)
  end subroutine utest_check_true_0
  subroutine utest_check_false_0(a, file, line, stopq)
    use mod_utest
    logical a
    character(*) file
    integer line
    logical stopq
    
    call utest_check_false(a, "a", file, line, stopq)
  end subroutine utest_check_false_0

  ! ==== using err_handle version ====
  subroutine expect_near_cmat(a, b, eps)
    use mod_utest
    use mod_err_handle
    implicit none
    
    complex(kind(0d0)), intent(in) :: a(:,:), b(:,:)
    double precision , intent(in) :: eps
    integer n, m, i, j

    call utest_check_begin
    
    if(size(a, 1) .ne. size(b, 1) .or. size(a, 2) .ne. size(b, 2)) then
       begin_err(1)
       call err_1("size mismatch")
       write(*,*) "left = ", size(a, 1), size(a, 2)
       write(*,*) "right = ", size(b, 1), size(b, 2)
       end_err()       
    end if
    
    n = size(a, 1)
    m = size(a, 2)
    
    do i = 1, n
       do j = 1, m
          if(abs(a(i,j)-b(i,j)) > eps) then
             begin_err(1)
             write(*,*) "elemnts are different"
             write(*,*) "i, j = ", i, j
             write(*,*) "left(i,j) = ", a(i,j)
             write(*,*) "right(i,j) = ", b(i,j)
             end_err()
          end if
       end do
    end do
    
  end subroutine expect_near_cmat
