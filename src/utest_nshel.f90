#include "macros.fpp"

module Mod_TestNshel
  use Mod_err_handle
  use Mod_utest_check
  use Mod_Nshel
  implicit none
contains
  subroutine TestNshel_run()
    call TestNshel_coef()
    call TestNshel_smat()
  end subroutine TestNshel_run
  subroutine TestNshel_coef()

    write(*,*)
    write(*,*) "--------------------"
    write(*,*) "TestNshel_coef begin"
    
    ![[[ 1.      0.      0.      0.      0.      0.    ]
    ![-0.2     0.5     0.      0.      0.      0.    ]
    ![ 0.54   -0.2     0.25    0.      0.      0.    ]]
    !
    ![[-0.1     0.5     0.      0.      0.      0.    ]
    ![ 0.52   -0.15    0.25    0.      0.      0.    ]
    ![-0.254   0.79   -0.125   0.125   0.      0.    ]]
    !
    ![[ 0.51   -0.1     0.25    0.      0.      0.    ]
    ![-0.202   0.775  -0.1     0.125   0.      0.    ]
    ![ 0.8154 -0.456   0.7825 -0.075   0.0625  0.    ]]]

    call expect_eq(1.0d0,  coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,0,0))
    call expect_eq(-0.2d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,1,0))
    call expect_eq(0.54d0, coef_d1(1.0d0,1.1d0,1.2d0,1.3d0, 0,2,0))

    write(*,*) "TestNshel_coef end"
    write(*,*) "--------------------"
    write(*,*) 
    
  end subroutine TestNshel_coef
  subroutine TestNshel_smat()
    type(Obj_Nshel) nshel
    integer :: ns(1,3)
    double precision :: coef_l(0:3,20)
    double precision, allocatable :: mat(:,:)

    write(*,*) "TestNshel_run begin"

    call Nshel_new(nshel, 2, 3); check_err()
    nshel%nucs%ws(1,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%ws(2,:) = (/0.0d0,0.0d0,0.0d0/);
    nshel%nucs%zs(:) = (/1.0d0, 2.0d0/)
    ns = 0
    coef_l = 0
    coef_l(0,1) = 1.0d0
    call Nshel_set(nshel, 1, (/"s"/), 1, (/1.1d0/), coef_l, 1); check_err()
    call Nshel_set(nshel, 2, (/"s"/), 1, (/1.3d0/), coef_l, 2); check_err()
    call Nshel_set(nshel, 3, (/"s"/), 1, (/1.4d0/), coef_l, 2); check_err()
    call Nshel_setup(nshel); check_err()

    !    call Nshel_dump(nshel); check_err()

    allocate(mat(nshel%nbasis, nshel%nbasis))
    call Nshel_s(nshel, mat); check_err()
    !write(*,*) "mat:"
    !write(*,*) mat

    call Nshel_delete(nshel); check_err()

    write(*,*) "TestNshel_run end"
    
  end subroutine TestNshel_smat
end module Mod_TestNshel

program main
  use Mod_err_handle
  use Mod_utest
  use Mod_utest_check
  use Mod_TestNshel

  call utest_begin
  call err_handle_begin

  write(*,*) "Hello"
  call TestNshel_run

  call err_handle_end
  call utest_end
  
end program main
