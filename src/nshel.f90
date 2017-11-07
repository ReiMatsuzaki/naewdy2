#include "macros.fpp"

module Mod_Shel
  use Mod_err_handle
  implicit none
  type Obj_Shel
     double precision, allocatable :: zeta(:)
     double precision, allocatable :: coef(:,:)
     integer, allocatable :: ns(:,:) ! ns(i,:) = (nx,ny,nz) for i th basis
     integer :: ng, num, maxn
     double precision :: w(3)
  end type Obj_Shel
contains
  subroutine Shel_new(this, ng, num)
    type(Obj_Shel) :: this
    integer, intent(in) :: ng, num

    this%ng = ng
    this%num = num
    allocate(this%zeta(ng))
    allocate(this%coef(num,ng))
    allocate(this%ns(num,3))

    this % zeta(:) = -1
    this % coef(:,:) = -1
    this % ns(:,:) = -1
    this % maxn = -1
    this % w(:) = 0
    

    !L = -1
    !do i = 1, num
    !   if(L < sum(ns(i,:))) then
    !      L = sum(ns(i,:))
    !   end if
    !end do
    !this%maxn = L
    !
    !if(size(coef_l,1)<L+1 ) then
    !   throw_err("size mismatch",1)
    !end if
    !if(size(ns,1)<num .or. size(ns,2)<ng) then
    !   throw_err("size mismatch",1)
    !end if
    !if(size(zeta)<ng .or. size(coef_l,2)<ng) then
    !   throw_err("size mismatch",1)
    !end if    
!
!    this%zeta(:) = zeta(1:ng)
!    do i = 1, num
!       L = sum(ns(i,:))
!       this%coef(i,:) = coef_l(L,1:ng)
!    end do
!    this%ns(:,:) = ns(1:num,1:3)
!
!    this%w(:) = w(:)
    
  end subroutine Shel_new
  subroutine Shel_delete(this)
    type(Obj_Shel) :: this
    deallocate(this%zeta, this%coef, this%ns)
  end subroutine Shel_delete
end module Mod_Shel

module Mod_Nucs
  use Mod_err_handle
  implicit none
  type Obj_Nucs
     integer :: num
     double precision, allocatable :: ws(:,:)
     double precision, allocatable :: zs(:)
  end type Obj_Nucs
contains
  subroutine Nucs_new(this, num)
    type(Obj_Nucs) :: this
    integer, intent(in) :: num
    this%num = num
    allocate(this%ws(num,3))
    allocate(this%zs(num))
  end subroutine Nucs_new
  subroutine Nucs_delete(this)
    type(Obj_Nucs) :: this
    deallocate(this%ws, this%zs)
  end subroutine Nucs_delete
!  subroutine Nucs_add(this, i, w, z, anum)
!    type(Obj_Nucs) :: this
!    integer, intent(in) :: i
!    double precision, intent(in) :: w(3)
!    double precision, intent(in) :: z
!    integer, intent(in) :: anum
!
!    if(i<1 .or. this%num<i) then
!       throw_err("out of range", 1)
!    end if
!
!    this%ws(i,:) = w
!    this%zs(i) = z
!    this%anum(i) = anum
!    
!  end subroutine Nucs_add
end module Mod_Nucs

Module Mod_Nshel
  use Mod_err_handle
  use Mod_Shel
  use Mod_Nucs  
  implicit none
  type Obj_Nshel
     integer :: num
     type(Obj_Shel), allocatable :: shels(:)
     integer, allocatable :: j0s(:)
     type(Obj_Nucs) :: nucs
     integer :: nbasis
  end type Obj_Nshel
contains
  subroutine Nshel_new(this, natom, num)
    type(Obj_Nshel) :: this
    integer, intent(in) :: num, natom
    this%num = num
    allocate(this%shels(num))
    allocate(this%j0s(num))
    call Nucs_new(this%nucs, natom); check_err()

    this%j0s(:) = -1
    
  end subroutine Nshel_new
  subroutine Nshel_delete(this)
    type(Obj_Nshel) :: this
    integer :: i
    do i = 1, this%num
       call Shel_delete(this%shels(i))
    end do
    deallocate(this%shels)
    call Nucs_delete(this%nucs)
  end subroutine Nshel_delete
  subroutine Nshel_dump(this, in_ifile)
    type(Obj_Nshel) :: this
    integer, intent(in) , optional :: in_ifile
    integer ifile, js, ia, j
    
    if(present(in_ifile)) then
       ifile = in_ifile
    else
       ifile = 6
    end if

    write(ifile,*)
    write(ifile,*) "==== Nshel ===="
    write(ifile,*) "---- nucs  ----"
    write(ifile,'("ia",4A10)') "x", "y", "z", "Q"
    do ia = 1, this%nucs%num
       write(ifile,'(i0,4x,4f10.5)') ia, this%nucs%ws(ia,:), this%nucs%zs(ia)
    end do
    do js = 1, this%num
       write(ifile,*) 
       write(ifile,'(A,I0,A)') " ---- jshel = ", js, " ----"
       write(ifile, '("w    = ", 3f10.5)') this%shels(js)%w(:)
       write(ifile, '("zeta = ", f10.5)') this%shels(js)%zeta(:)
        write(ifile,*) "j     n         coef"
       do j = 1, this%shels(js)%num
          write(ifile,'(I0,2x,3i3,2x,f10.5)') j, this%shels(js)%ns(j,:), this%shels(js)%coef(j,:)
       end do
    end do
    
  end subroutine Nshel_dump
  subroutine Nshel_set(this, js, ntypes, ng, zeta, coef_l, ia)
    type(Obj_Nshel) :: this
    integer, intent(in) :: js, ng
    character(*), intent(in) :: ntypes(:)
    double precision, intent(in) :: zeta(:)
    double precision, intent(in) :: coef_l(0:,:)
    integer, intent(in) :: ia
    integer it, num, j, jj

    num = 0
    do it = 1, size(ntypes)
       if(is_in_s(ntypes, "s")) then
          num = num + 1
       else if(is_in_s(ntypes, "p")) then
          num = num + 3
       else if(is_in_s(ntypes, "d")) then
          num = num + 6
       else
          throw_err("unsupported", 1)
       end if
    end do

    call Shel_new(this%shels(js), ng, num); check_err()

    this%shels(js)%zeta(:) = zeta(1:ng)
    j = 1
    do it = 1, size(ntypes)
       if(is_in_s(ntypes, "s")) then
          this%shels(js)%ns(j,:) = (/0,0,0/)
          this%shels(js)%coef(j,:) = coef_l(0,1:ng)
          j = j + 1
       else if(is_in_s(ntypes, "p")) then
          this%shels(js)%ns(j+0,:) = (/1,0,0/)
          this%shels(js)%ns(j+1,:) = (/0,1,0/)
          this%shels(js)%ns(j+2,:) = (/0,0,1/)
          do jj = j, j+2
             this%shels(js)%coef(jj,:) = coef_l(1,1:ng)
          end do
          j = j + 3
       else if(is_in_s(ntypes, "d")) then
          this%shels(js)%ns(j+0,:) = (/2,0,0/)
          this%shels(js)%ns(j+1,:) = (/0,2,0/)
          this%shels(js)%ns(j+2,:) = (/0,0,2/)
          this%shels(js)%ns(j+3,:) = (/1,1,0/)
          this%shels(js)%ns(j+4,:) = (/0,1,1/)
          this%shels(js)%ns(j+5,:) = (/1,0,1/)
          do jj = j, j+2
             this%shels(js)%coef(jj,:) = coef_l(2,1:ng)
          end do
          j = j + 6
       else
          throw_err("unsupported", 1)
       end if
    end do
    
    this%shels(js)%maxn = maxval(this%shels(js)%ns(:,:))
    this%shels(js)%w(:) = this%nucs%ws(ia,:)
    
  end subroutine Nshel_set
  subroutine Nshel_setup(this, in_normalize)
    type(Obj_Nshel) :: this
    logical ,intent(in), optional :: in_normalize
    integer :: j0, i, j, js, jj
    double precision, allocatable :: smat(:,:)
    logical normalize

    j0 = 0
    do i = 1, this%num
       this%j0s(i) = j0
       j0 = j0 + this%shels(i)%num
       write(*,*) this%shels(i)%num
    end do

    this%nbasis = j0

    if(present(in_normalize)) then
       normalize = in_normalize
    else
       normalize = .true.
    end if
    if(normalize) then
       allocate(smat(this%nbasis, this%nbasis))
       call Nshel_s(this, smat); check_err()
       do js = 1, this%num
          do jj = 1, this%shels(js)%num
             j = this%j0s(js) + jj
             this%shels(js)%coef(jj,:) = this%shels(js)%coef(jj,:)/sqrt(smat(j,j))
          end do
       end do
    end if
    
  end subroutine Nshel_setup
  subroutine Nshel_s(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    double precision :: mat(:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), i, j, k
    double precision :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    double precision :: d(3,0:5,0:5,0:10), acc, coef

    write(*,*) "wa!"
    do js = 1, this%num
       do ks = 1, this%num    
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk,0, d); check_err()

                do jj = 1, this%shels(js)%num
                   do kk = 1, this%shels(ks)%num                      
                      nj(:) = this%shels(js)%ns(jj,:)
                      nk(:) = this%shels(ks)%ns(kk,:)
                      acc = 1
                      do i = 1, 3
                         acc = acc * d(i,nj(i),nk(i),0)
                      end do
                      coef = cp * this%shels(js)%coef(jj,jg) &
                           * this%shels(ks)%coef(kk,kg)
                      j = this%j0s(js) + jj
                      k = this%j0s(js) + kk
                      mat(j,k) = mat(j,k) + coef*acc
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine Nshel_s
  subroutine Nshel_t(this, mat)
    use Mod_const, only : pi
    type(Obj_Nshel) :: this
    double precision :: mat(:,:)
    integer js, ks, jg, kg, maxnj, maxnk, jj, kk, nj(3), nk(3), ir, jr, j, k, nkp(3)
    double precision :: wj(3), wk(3), d2, zj, zk, zp, wp(3), ep, cp
    double precision :: d(3,0:5,0:5,0:10), acc, coef, cum

    do js = 1, this%num
       do ks = 1, this%num
          wj(:) = this%shels(js)%w(:)
          wk(:) = this%shels(ks)%w(:)
          d2 = dot_product(wj-wk,wj-wk)
          maxnj = this%shels(js)%maxn
          maxnk = this%shels(ks)%maxn

          do jg = 1, this%shels(js)%ng
             do kg = 1, this%shels(ks)%ng
                zj = this%shels(js)%zeta(jg)
                zk = this%shels(ks)%zeta(kg)
                zp = zj+zk
                wp(:) = (zj*wj(:) + zk*wk(:))/zp
                ep = exp(-zj*zk/zp*d2)
                cp = ep*(pi/zp)**(1.5)
                call coef_d(zp,wp,wj,wk,maxnj,maxnk,0, d); check_err()

                do jj = 1, this%shels(js)%num
                   do kk = 1, this%shels(ks)%num
                      nj(:) = this%shels(js)%ns(jj,:)
                      nk(:) = this%shels(ks)%ns(kk,:)

                      cum = 0
                      acc = 1
                      do ir = 1, 3
                         acc = acc * d(ir,nj(ir),nk(ir),0)
                      end do
                      cum = cum + 2*zk*(2*sum(nk)+3)*acc

                      do jr = 1, 3
                         nkp(:) = nk(:); nkp(jr) = nkp(jr)+2
                         acc = 1
                         do ir = 1, 3
                            acc = acc * d(ir,nj(ir),nkp(ir),0)
                         end do
                         cum = cum + 4*zk*zk*acc
                         
                         if(nk(jr)<2) continue
                         nkp(:) = nk(:); nkp(jr) = nkp(jr)-2
                         acc = 1
                         do ir = 1, 3
                            acc = acc * d(ir,nj(ir),nkp(ir),0)
                         end do
                         cum = cum + nk(jr)*(nk(jr)-1)*acc
                         
                      end do
                      
                      coef = cp * this%shels(js)%coef(jj,jg) &
                           * this%shels(ks)%coef(kk,kg)
                      j = this%j0s(js) + jj
                      k = this%j0s(ks) + kk
                      mat(j,k) = mat(j,k) + coef*acc
                   end do
                end do
             end do
          end do
       end do
    end do    
  end subroutine Nshel_t
  ! == calculation functions == 
  recursive function coef_d1(zp,wp,wj,wk,nj,nk,n) result(res)
    double precision, intent(in) :: zp, wp, wj, wk
    integer, intent(in) :: nj, nk, n
    double precision :: res

    if(nj==0 .and. nk==0 .and. n==0) then
       res = 1.0
    else if(n<0 .or. n>nj+nk) then
       res = 0.0
    else if(nj>0) then
       res = 1/(2*zp) * coef_d1(zp,wp,wj,wk,nj-1,nk,n-1) + &
            (wp-wj)   * coef_d1(zp,wp,wj,wk,nj-1,nk,n)   + &
            (n+1)     * coef_d1(zp,wp,wj,wk,nj-1,nk,n+1)
    else
       res = 1/(2*zp) * coef_d1(zp,wp,wj,wk,nj,nk-1,n-1) + &
            (wp-wk)   * coef_d1(zp,wp,wj,wk,nj,nk-1,n)   + &
            (n+1)     * coef_d1(zp,wp,wj,wk,nj,nk-1,n+1)
    end if
    
  end function coef_d1
  subroutine coef_d(zp,wp,wj,wk,maxnj,maxnk,maxn,res)
    double precision, intent(in) :: zp, wp(3), wj(3), wk(3)
    integer, intent(in) :: maxnj, maxnk, maxn
    double precision, intent(out) :: res(:,0:,0:,0:)
    integer i, nj, nk, n

    do i = 1, 3
       do nj = 0, maxnj
          do nk = 0, maxnk
             do n = 0, maxn
                res(i,nj,nk,n) = coef_d1(zp,wp(i),wj(i),wk(i),nj,nk,n)
             end do
          end do
       end do
    end do
    
  end subroutine coef_d
  recursive function mole_gammainc1(m, z) result(res)
    integer, intent(in) :: m
    double precision, intent(in) :: z
    double precision :: res
    res = z**m
  end function mole_gammainc1
  subroutine mole_gammainc(maxm, z, res)
    integer, intent(in) :: maxm
    double precision, intent(in) :: z
    double precision :: res(0:)
    integer m
    do m = 0, maxm
       res(m) = mole_gammainc1(m, z)
    end do
  end subroutine mole_gammainc
  recursive function coef_R1(zp,wp,wc,n,j) result(res)
    double precision :: zp, wp(3), wc(3)
    integer, intent(in) :: n(3), j
    double precision :: res, wpc(3), d2, Fj(0:10)

    wpc(:) = wp(:)-wc(:)

    call mole_gammainc(j,zp*d2, Fj); check_err()
    
    if(all(n==0)) then
       d2 = dot_product(wpc, wpc)
       res = (-2*zp)**j * Fj(0)
    end if
    
  end function coef_R1
  ! == utils ==
  function is_in_s(ss, s) result(res)
    character(*), intent(in) :: ss(:)
    character(*) :: s
    logical :: res
    integer i
    res = .false.
    do i = 1, size(ss)
       if(trim(ss(i)) .eq. trim(s)) then
          res = .true.
       end if
    end do    
  end function is_in_s
end Module Mod_Nshel
