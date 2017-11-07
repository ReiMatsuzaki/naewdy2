#include "macros.fpp"

module Mod_gauss
  use Mod_const, only : ii
  implicit none
  double precision :: R0_, P0_
  complex(kind(0d0)) :: alpha_
  complex(kind(0d0)) :: A_
contains
  subroutine gauss_at(x, y)
    double precision, intent(in) :: x
    complex(kind(0d0)), intent(out) :: y
    y = A_*exp(-alpha_*(x-R0_)**2 + ii*P0_*(x-R0_))
  end subroutine gauss_at
end module Mod_gauss

module Mod_math
  use mod_err_handle
  use mod_const
  implicit none
contains
  function tdot(a, b) result(res)
    complex(kind(0d0)), intent(in) :: a(:), b(:)
    complex(kind(0d0)) res
    integer i, n
    res = 0.0d0
    n = min(size(a), size(b))
    do i = 1, n
       res = res + a(i)*b(i)
    end do
  end function tdot
  function vmv(a, S, b) result(res)
    complex(kind(0d0)), intent(in) :: a(:), S(:,:), b(:)
    complex(kind(0d0)) res
    res = dot_product(a, matmul(S, b))
  end function vmv
  subroutine poly_fit(xs, ys, n, cs)
    double precision, intent(in) :: xs(:)
    double precision, intent(in) :: ys(:)
    integer, intent(in)          :: n
    double precision, intent(out) :: cs(0:)
    double precision :: mat(0:n, 0:n)
    double precision :: val
    integer :: k, j, i, m, info
    integer :: ipiv(n+1)
    
    m = size(xs)
    if(m .ne. size(ys)) then
       throw_err("size mismatch", 1)
    end if
    
    do k = 0, n
       do j = 0, n
          val = 0.0d0
          do i = 1, m
             val = val + xs(i)**(k+j)
          end do
          mat(k, j) = val
       end do
    end do
    do k = 0, n
       val = 0.0d0
       do i = 1, m
          val = val + ys(i)*xs(i)**k
       end do
       cs(k) = val
    end do
    
    ! -- call lapack for solving linear problem --
    call dgesv(n+1, 1, mat, n+1, ipiv, cs, n+1, info)
    if(info .ne. 0) then
       throw_err("dgesv failed", 1)
    end if
    
  end subroutine poly_fit
  subroutine schmidt_orthogonalize(S, U)
    !
    !  from overlap matrix of nonorthgozalized basis S, gives
    !  transformation matrix to orthogonalized basis
    !
    !  Inputs
    !  ------
    !   S : overlap matrix of non orthogonalized basis
    !  Returns
    !  -------
    !   U : transformation matrix
    !
    !
    !  u : primitive basis
    !  v : orthogonalize basis
    !  v_k = sum_{l=1}^k U(l,k)u_l
    !
    !  v_1 = u_1
    !  v_k = u_k -sum_{l=1}^{k-1} <v_l,u_k>/<v_l,v_l> v_l
    !      = u_k -sum_{l=1}^{k-1} <v_l,u_k>/<v_l,v_l> sum_{m=1}^l U(m,l)u_l
    
    complex(kind(0d0)), intent(in) :: S(:, :)
    complex(kind(0d0)), intent(out) :: U(:, :)
    integer k, l, m, num
    complex(kind(0d0)) :: vl_uk, vl_vl
    
    if(size(S, 1) .ne. size(S, 2)) then
       throw_err("S must be square matrix", 1)
    end if
    if(size(U, 1) .ne. size(U, 2)) then
       throw_err("U must be square matrix", 1)
    end if
    if(size(U, 1) .ne. size(S, 2)) then
       throw_err("S and U must be same shape", 1)
    end if

    num = size(S, 1)
    U(:, :) = 0.0d0
    
    U(1, 1) = 1.0d0

    do k = 2, num
       U(:, k) = 0.0d0
       U(k, k) = 1.0d0
       do l = 1, k-1
          vl_uk = dot_product(U(:, l), S(:, k))
          vl_vl = dot_product(U(:, l), matmul(S(:,:), U(:, l)))
          do m = 1, l
             U(m, k) = U(m, k) - U(m, l)*vl_uk/vl_vl
          end do
       end do
    end do

    do k = 1, num
       vl_vl = dot_product(U(:, k), matmul(S(:,:), U(:, k)))
       U(:, k) = U(:, k) / sqrt(vl_vl)
    end do
    
  end subroutine schmidt_orthogonalize
  function expect_dot(A, c) result(res)
    complex(kind(0d0)) A(:,:)
    complex(kind(0d0)) c(:)
    complex(kind(0d0)) :: res
    integer n, i, j

    res = 0.0d0
    n = size(c)
    if(size(A, 1).ne.n .or. size(A,2).ne.n) then
       begin_err(1)
       call err_1("size mismatch")
       call err_si("size(c):", n)
       call err_si("size(A,1):", size(A, 1))
       call err_si("size(A,2):", size(A, 2))
       end_err()
    end if
    res = 0.0d0
    do i = 1, n
       do j = 1, n
          res = res + conjg(c(i)) * A(i,j) * c(j)
       end do
    end do
  end function expect_dot
  subroutine dump_dmat(x, ifile, in_eps)
    double precision, intent(in) :: x(:,:)
    integer, intent(in) :: ifile
    double precision, optional, intent(in) :: in_eps
    integer i, j
    double precision :: eps

    if(present(in_eps)) then
       eps = in_eps
    else
       eps = -1.0d0
    end if
    
    write(ifile,'("i,j,val")') 
    do i = 1, size(x,1)
       do j = 1, size(x,2)
          if(abs(x(i,j)) > eps) then
             write(ifile,'(i0,",",i0,",",f20.10)') i,j,x(i,j)
          end if
       end do
    end do
  end subroutine dump_dmat
  subroutine dump_dvec(x, ifile)
    double precision, intent(in) :: x(:)
    integer, intent(in) :: ifile
    integer i
    write(ifile,'("i,val")') 
    do i = 1, size(x,1)
       write(ifile,'(i0,",",f20.10)') i,x(i)
    end do
  end subroutine dump_dvec
  subroutine t2s(t,s)
    ! Copied from mangan4 written by K.Yamamoto
    ! Matrix transformation: packet strage -> symmetric.
    implicit none
    real(kind(0d0)),intent(in) :: t(:)
    real(kind(0d0)),intent(out) :: s(:,:)
    integer :: irun,i,j
    i = 1
    j = 1
    s(:,:) = 0d0
    do irun = 1,SIZE(t)
       s(i,j) = t(irun)
       s(j,i) = t(irun)
       j = j + 1
       if (j > i) then
          j = 1
          i = i + 1
       end if
    end do
  end subroutine t2s
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
  subroutine print_mat(name, X, ifile)
    character(*), intent(in) :: name
    complex(kind(0d0)), intent(in) :: X(:, :)
    integer,  optional :: ifile
    integer n, i, ifi

    if(present(ifile)) then
       ifi = ifile
    else
       ifi = 6
    end if
    
    n = size(X, 1)
    write(ifi,*) "mat:", name
    do i = 1, n
       write(ifi,'(I4)', advance='no') i
       write(ifi,*) X(i, :)
    end do
  end subroutine print_mat
  subroutine print_dmat(name, X)
    character(*), intent(in) :: name
    double precision, intent(in) :: X(:, :)
    integer n, m, i, j
    n = size(X, 1)
    m = size(X, 2)
    write(*,*) "mat:", name
    do i = 1, n
       do j = 1, m
          write(*,'(2I5, F20.10)') i, j, X(i, j)
       end do
    end do
  end subroutine print_dmat
end module Mod_math

module Mod_Sparse
  use mod_err_handle
  implicit none
  type Sparse
     integer :: n, size, capacity
     integer, allocatable :: Is(:)
     integer, allocatable :: Js(:)
     double precision , allocatable :: vs(:)
  end type Sparse
contains
  function Sparse_allocated(this) result(res)
    type(Sparse) :: this
    logical :: res
    res = allocated(this%Is)
  end function Sparse_allocated
  subroutine Sparse_new(this, n, capa)
    type(Sparse) this
    integer, intent(in) :: n, capa
    this%n=n
    this%size = 0
    this%capacity = capa
    allocate(this%Is(capa))
    allocate(this%Js(capa))
    allocate(this%vs(capa))
  end subroutine Sparse_new
  subroutine Sparse_delete(this)
    type(Sparse) this
    deallocate(this%Is)
    deallocate(this%Js)
    deallocate(this%vs)
  end subroutine Sparse_delete
  subroutine Sparse_dump(this)
    type(Sparse) :: this
    integer r
    write(*,*) "naewady:sparse:dump begin"
    write(*,*) "n:", this%n
    write(*,*) "size:", this%size
    do r = 1, this%size
       write(*,*) this%Is(r), this%Js(r), this%vs(r)
    end do
    write(*,*) "naewady:sparse:dump end"
  end subroutine Sparse_dump
  subroutine Sparse_clear(this)
    type(Sparse) :: this
    this%size=0
  end subroutine Sparse_clear
  subroutine Sparse_set_capa(this, capa)
    type(Sparse) :: this
    integer, intent(in) :: capa
    integer :: n
    n = this%n
    if(this%capacity<capa) then
       call Sparse_delete(this)
       call Sparse_new(this, n, capa)
    end if
  end subroutine Sparse_set_capa
  subroutine Sparse_add(this, I, J, v)
    type(Sparse) :: this
    integer, intent(in) :: I, J
    double precision, intent(in) :: v
    integer idx
    
    idx = this%size+1

    if(idx>size(this%Is)) then
       write(*,*) "size exceed"
       write(*,*) "idx", idx
       write(*,*) "size(Is)", size(this%Is)
       write(*,*) "capa", this%capacity
       stop
    end if

    this%Is(idx) = I
    this%Js(idx) = J
    this%vs(idx) = v

    this%size=idx
    
  end subroutine Sparse_add
  subroutine Sparse_axpy(this, alpha, x, beta, y)
    type(Sparse) :: this
    complex(kind(0d0)) , intent(in) :: alpha, beta
    complex(kind(0d0)), intent(in) :: x(this%n)
    complex(kind(0d0)), intent(inout) :: y(this%n)
    
    integer idx, I, J
    do idx = 1, this%size
       I = this%Is(idx)
       J = this%Js(idx)
       y(I) = alpha*this%vs(idx)*x(J) + beta*y(I)
    end do
  end subroutine Sparse_axpy
  
end module Mod_Sparse

