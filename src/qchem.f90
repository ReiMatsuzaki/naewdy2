#include "macros.fpp"
#define ACTIVE_GAMESS 1

module Mod_gamess_common
  integer, PARAMETER :: MXATM=2000
  integer, parameter :: MXAO=8192
  integer NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, IAN
  REAL(KIND(0D0)) :: ZAN, C
  COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,&
       &                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
  INTEGER :: NWKS,NORBS,IB,JB,IUWK,JUWK,NUWK,NLWK,IHAI,LDGUGA
  double precision :: acof
  COMMON /LOOPS1/ ACOF,NWKS,NORBS,IB,JB,IUWK,JUWK,NUWK,NLWK, IHAI,LDGUGA
  double precision ZMASS
  COMMON /MASSES/ ZMASS(MXATM)  
end module Mod_gamess_common

module Mod_aij  
  implicit none
  integer, parameter :: MAXNUM=2**12
  integer :: ihb_  !! start memory location of aij in gamess
  ! -- intermediate data --
  integer, allocatable :: iwork_(:,:), ijsort_(:), conf_(:,:)
  double precision, allocatable :: rwork_(:)
  integer :: idx_, nex_
  logical :: onaij_
contains
  subroutine extend
    integer n0,n
    integer, allocatable :: iwt(:,:)
    double precision, allocatable :: wt(:)
    n0 = MAXNUM*nex_
    allocate(iwt(4,n0),wt(n0))
    iwt(:,:) = iwork_(:,:)
    wt(:) = rwork_(:)
    nex_ = nex_+1
    n = MAXNUM*nex_
    deallocate(iwork_,rwork_)
    allocate(iwork_(4,n),rwork_(n))
    iwork_(:,:) = 0
    rwork_(:) = 0
    iwork_(:,:n0) = iwt(:,:n0)
    rwork_(:n0) = wt(:n0)    
  end subroutine extend
  subroutine aij_new
    !
    !! ijsort(ijguga) : ijmo for ijguga
    allocate(ijsort_(MAXNUM), iwork_(3, MAXNUM), rwork_(MAXNUM))
    call aij_reset
    onaij_ = .true.
    nex_ = 1
  end subroutine aij_new
  subroutine aij_reset
    idx_ = 0
    onaij_ = .true.
  end subroutine aij_reset
  subroutine aij_setihb(ihb)
    integer, intent(in) :: ihb
    ihb_ = ihb
  end subroutine aij_setihb
  subroutine aij_append(ii, jj, ijguga, acoef)
    use Mod_err_handle
    integer, intent(in) :: ii, jj, ijguga
    double precision, intent(in) :: acoef
    idx_ = idx_ + 1
    if(idx_ >= MAXNUM*nex_) then
       call extend
    end if
    if(size(iwork_,2)<idx_ .or. size(iwork_,1)<3) then
       begin_err(1)
       write(0,*) "size of iwork_ is too small"
       write(0,*) "size(iwork_):", size(iwork_,1), size(iwork_,2)
       write(0,*) "idx_:", idx_
       write(0,*) "ierr:", ierr
       end_err()
    end if
    if(size(rwork_)<idx_) then
       begin_err(1)
       write(0,*) "size of rwork_ is too small"
       write(0,*) "size(rwork_):", size(rwork_)
       write(0,*) "idx_:", idx_
       write(0,*) "ierr:", ierr
       end_err()
    end if
    iwork_(1, idx_) = ii - ihb_ + 1
    iwork_(2, idx_) = jj - ihb_ + 1
    iwork_(3, idx_) = ijguga
    rwork_(idx_) = acoef
 !   write(*,*) "append:", size(rwork_), idx_
  end subroutine aij_append
  subroutine aij_ijsort(ijmo, ijguga)
    integer, intent(in) :: ijmo, ijguga
    !write(*,*) "aij_ijsort", size(ijsort_) = ijmo
    if(size(ijsort_)<ijguga) then
       write(*,*) "aij_ijsort out of range", size(ijsort_), ijguga
       call abort()
    end if
    ijsort_(ijguga) = ijmo
!    write(*,*) "ijsort:", size(ijsort_), ijguga
  end subroutine aij_ijsort
  subroutine aij_setconf(iwks, ncore, npart, ieconf, num, nwks)
    integer, intent(in) :: iwks, ncore, npart, ieconf(num), num, nwks
    integer :: i

    if(allocated(conf_)) then
       if(size(conf_,1)<num .or. size(conf_,2)<nwks) then
          deallocate(conf_)
       end if
    end if
    
    if(.not.allocated(conf_)) then
       allocate(conf_(num, nwks))
       conf_ = 0
    end if

    do i = 1, num
       if(i <= ncore) then
          conf_(i, iwks) = 2
       else if(i <= ncore + npart) then
          conf_(i, iwks) = ieconf(i)
       else
          conf_(i, iwks) = 0
       end if
    end do    
    
  end subroutine aij_setconf
  subroutine it2s(ij, i, j)
    implicit none
    integer, intent(in) :: ij
    integer, intent(out) :: i, j
    double precision :: vij, vi
    vij = DBLE(ij)
    vi = (1 + SQRT(1.0d0 + 8.0d0 * (vij -1.0d0))) / 2.0d0
    i = INT(AINT(vi))
    j = ij - i*(i-1) / 2
  end subroutine it2s
  subroutine aij_calc(ia, va, na)
    integer, allocatable , intent(out) :: ia(:,:)
    double precision, allocatable, intent(out) :: va(:)
    integer, intent(out) :: na
    integer :: i,icsf,jcsf,imo,jmo,nmo,ncsf
    integer,allocatable :: conf1(:),conf2(:)
    onaij_ = .false.
    na = idx_
    if(.not.allocated(ia)) then
       allocate(ia(4, na), va(na))
    end if
    nmo = SIZE(conf_(:,1))
    ncsf = SIZE(conf_(1,:))
    allocate(conf1(nmo),conf2(nmo))
    do i = 1,idx_
       icsf = iwork_(1,i)
       jcsf = iwork_(2,i)
       call it2s(ijsort_(iwork_(3,i)),imo,jmo)
       ! Formal anihiration operations.
       conf1(:) = conf_(:,icsf)
       conf1(imo) = conf1(imo) - 1
       conf2(:) = conf_(:,jcsf)
       conf2(jmo) = conf2(jmo) - 1
       conf1(:) = abs(conf1(:) - conf2(:))
       if (maxval(conf1) > 0) then
          call it2s(ijsort_(iwork_(3,i)),jmo,imo)
       end if
       ia(1,i) = icsf
       ia(2,i) = jcsf
       ia(3,i) = imo
       ia(4,i) = jmo
       va(i) = rwork_(i)
    end do
    deallocate(conf1, conf2)
  end subroutine aij_calc
  subroutine aij_delete
    deallocate(iwork_)
    onaij_ = .false.
  end subroutine aij_delete
end module Mod_aij
      
module Mod_nshel
  implicit none
  INTEGER,PARAMETER :: MXSH=5000, MXGTOT=20000
  DOUBLE PRECISION :: EX,CS,CP,CD,CF,CG,CH,CI
  INTEGER :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL
  COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),&
        CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT), &
        KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),&
        KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
  integer, PARAMETER :: MXATM=2000
  integer, parameter :: MXAO=8192
  integer NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, IAN
  REAL(KIND(0D0)) :: ZAN, C
  COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,&
       &                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
contains
  subroutine object_set_ivec(o, name, x)
    use Mod_fjson
    type(object) :: o
    character(*) :: name
    integer :: x(:)
    type(array) :: a
    call ivec2a(x, a); check_err()
    call object_set_a(o, name, a); check_err()
    call array_delete(a); check_err()
  end subroutine object_set_ivec
  subroutine object_set_dvec(o, name, x)
    use Mod_fjson
    type(object) :: o
    character(*) :: name
    double precision  :: x(:)
    type(array) :: a
    call dvec2a(x, a); check_err()
    call object_set_a(o, name, a); check_err()
    call array_delete(a); check_err()
  end subroutine object_set_dvec
  subroutine Nshel_dump(ifile)
    use Mod_fjson
    integer, intent(in) :: ifile
    type(value) :: v
    type(object) :: o
    type(array) :: a
    integer ng, i1(NAT)
    double precision :: cc(3,NAT), c1(NAT)

    ng = sum(KNG(:NSHELL))
    call object_set_i(o, "ng", ng); check_err()
    call object_set_i(o, "nshell", NSHELL); check_err()
    call object_set_dvec(o, "ex", EX(:ng)); check_err()
    call object_set_dvec(o, "cs", CS(:ng)); check_err()
    call object_set_dvec(o, "cp", CP(:ng)); check_err()
    call object_set_dvec(o, "cd", CD(:ng)); check_err()
    call object_set_dvec(o, "cf", CF(:ng)); check_err()
    call object_set_dvec(o, "cg", CG(:ng)); check_err()
    call object_set_dvec(o, "ch", CH(:ng)); check_err()
    call object_set_dvec(o, "ci", CI(:ng)); check_err()
    
    call object_set_ivec(o, "kstart", KSTART(:NSHELL)); check_err()
    call object_set_ivec(o, "katom",  KATOM(:NSHELL)); check_err()
    call object_set_ivec(o, 'ktype',KTYPE(:NSHELL)); check_err()
    call object_set_ivec(o, 'kng',  KNG(:NSHELL)); check_err()
    call object_set_ivec(o, 'kloc', KLOC(:NSHELL)); check_err()
    call object_set_ivec(o, 'kmin', KMIN(:NSHELL)); check_err()
    call object_set_ivec(o, 'kmax', KMAX(:NSHELL)); check_err()

    cc(:,:) = C(1:3,1:NAT)
    call dmat2a(cc, a); check_err()
    call object_set_a(o, "c", a); check_err()
    call array_delete(a); check_err()

    c1(:) = ZAN(1:NAT)
    call dvec2a(c1, a); check_err()
    call object_set_a(o, "zan", a); check_err()
    call array_delete(a);; check_err()

    i1(:) = IAN(1:NAT)
    call ivec2a(i1, a); check_err()
    call object_set_a(o, "ian", a); check_err()
    call array_delete(a); check_err()

    call value_new_o(v, o); check_err()
    call value_dump(v, ifile); check_err()

    call value_delete(v); check_err()
    call object_delete(o); check_err()
    
  end subroutine Nshel_dump
end module Mod_nshel

module Mod_asgn
  implicit none
  private
  real(kind(0d0)),parameter :: TOL=1d-20
  !---------------------------------------------------------------
  public asgn_asgn
  !---------------------------------------------------------------  
contains
  subroutine asgn_asgn(s,c,e,c0)
    use Mod_math
    implicit none
    real(kind(0d0)),intent(in) :: c0(:,:)
    real(kind(0d0)),intent(inout) :: c(:,:),e(:)
    real(kind(0d0)),intent(out) ::s(:,:)
    integer :: n
    n = SIZE(e)
    call mm('T',c,c0,s)   ! s <- c^T.c0    
    call decouple(e,s)   
    call orthnorm(s)    
    call transf(s,c,e)
  end subroutine asgn_asgn
  subroutine orthnorm(s)
    implicit none
    real(kind(0d0)),intent(inout) :: s(:,:)
    real(kind(0d0)) :: tmp
    integer :: i,j
    call vv(s(:,1),s(:,1),tmp)
    s(:,1) = s(:,1) / SQRT(tmp)
    do i = 2,SIZE(s(:,1))
       do j = 1,i-1
          call vv(s(:,i),s(:,j),tmp)
          s(:,i) = s(:,i) - tmp * s(:,j)
       end do
       call vv(s(:,i),s(:,i),tmp)
       s(:,i) = s(:,i) / SQRT(tmp)
    end do
  end subroutine orthnorm
  subroutine transf(s,c,e)
    implicit none
    real(kind(0d0)),intent(in) :: s(:,:)
    real(kind(0d0)),intent(inout) :: e(:),c(:,:)
    integer :: nb,ns
    real(kind(0d0)),allocatable :: et(:),ct(:,:)
    nb = SIZE(c(:,1))
    ns = SIZE(e)
    allocate(et(ns),ct(nb,ns))
    call mv('T',s**2,e,et)
    call mm('T',TRANSPOSE(c),s,ct)
    e(:) = et(:)
    c(:,:) = ct(:,:)
  end subroutine transf
  subroutine decouple(e,s)
    implicit none
    real(kind(0d0)),intent(in) :: e(:)
    real(kind(0d0)),intent(inout) :: s(:,:)
    real(kind(0d0)) :: vmax
    integer :: ns,i,j,jmax
    integer,allocatable :: jls(:)
    ns = SIZE(e)
    allocate(jls(ns))
    jls(:) = 0
    jmax = 0
    do i = 1,ns
       vmax = 0.0d0
       do j = 1,ns
          if (is_used(j,jls)) then
             cycle
          end if
          if (s(i,j)**2 >= vmax) then
             vmax = s(i,j)**2
             jmax = j
          end if
       end do
       jls(i) = jmax
       do j = 1,ns
          if (j /= jmax .and. (e(jmax)-e(j))**2 > TOL) then
             s(i,j) = 0.0d0
          end if
       end do
    end do
  end subroutine decouple
  logical function is_used(j,jls)
    implicit none
    integer,intent(in) :: j,jls(:)
    integer :: i
    is_used = .false.
    do i = 1,SIZE(jls)
       if (jls(i) == 0) then
          exit
       end if
       if (jls(i) == j) then
          is_used = .true.
          exit
       end if
    end do
  end function is_used
  subroutine vv(x,y,z)
    ! Scalar product.
    implicit none
    real(kind(0d0)),intent(in) :: x(:),y(:)
    real(kind(0d0)),intent(out) :: z
    z = SUM(x(:) * y(:))
  end subroutine vv
  subroutine mm(t,a,b,c)
    ! matrix x matrix product.
    implicit none
    character,intent(in) :: t
    real(kind(0d0)),intent(in) :: a(:,:),b(:,:)
    real(kind(0d0)),intent(out) :: c(:,:)
    real(kind(0d0)),parameter :: ONE=1d0,ZERO=0d0
    integer :: m,n,k
    m = SIZE(c(:,1))
    n = SIZE(c(1,:))
    k = SIZE(b(:,1))
    call DGEMM(t,'N',m,n,k,ONE,a,k,b,k,ZERO,c,m)
  end subroutine mm
  subroutine mv(t,a,x,y)
    ! matrix x vector product.
    implicit none
    character,intent(in) :: t
    real(kind(0d0)),intent(in) :: a(:,:),x(:)
    real(kind(0d0)),intent(out) :: y(:)
    real(kind(0d0)),parameter :: ONE=1d0,ZERO=0d0
    integer :: m,n
    m = SIZE(a(:,1))
    n = SIZE(a(1,:))
    call DGEMV(t,m,n,ONE,a,m,x,1,ZERO,y,1)
  end subroutine mv
end module Mod_asgn

module Mod_QChem
  use Mod_err_handle
  use Mod_aij
  implicit none
  integer :: nao_, nmo_, nwks_, nstate_
  logical :: get_hIJ_
  double precision, allocatable :: hIJ_(:,:)
  character(10) :: type_hmat_ = ""
  integer na_ ! number of a coef
  integer, allocatable          :: ia_(:,:) ! index of a coeff
  double precision, allocatable :: va_(:)   ! value of a coeff  
  ! -- for prev data --
  logical :: use_mo_asgn_, done_mo_asgn_
  double precision, allocatable :: mo_sc0_(:,:), mo_asgn_(:,:)
  logical ::  use_ci_asgn_, done_ci_asgn_
  double precision, allocatable :: ci_c0_(:,:), ci_asgn_(:,:)
contains
  ! == Constructors ==
  subroutine QChem_new
    use Mod_gamess_common, only : NAT, NUM
    nao_ = NUM     ! number of one particle orbital(AO and MO)
    nmo_ = NUM     ! number of one particle orbital(AO and MO)
    nwks_ = 0       ! number of the walks (CSFs)
    nstate_ = 0     ! number of CI state
    type_hmat_ = "dense"
    use_mo_asgn_ = .false.; done_mo_asgn_ = .false.
    use_ci_asgn_ = .false.; done_ci_asgn_ = .false.

    call aij_new
    call QChem_nameio
    
  end subroutine QChem_new
  subroutine QChem_nameio
    integer, parameter :: NNAM=13, MXAO=8192
    double precision :: SGUESS, QNAM(NNAM), GUESS
    integer ::  KQNAM(NNAM), mix
    integer JRET,NORB,NORDER,IORDER(MXAO),JORDER(MXAO),insorb
    double precision tolz,tole,prtmo,punmo,symden,purify
    integer IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
    DATA SGUESS/8HGUESS   /
    DATA QNAM/8HGUESS   ,8HNORB    ,8HNORDER  ,8HIORDER  ,8HJORDER  ,&
              8HTOLZ    ,8HTOLE    ,8HMIX     ,8HPRTMO   ,8HPUNMO   ,&
              8HSYMDEN  ,8HPURIFY  ,8HINSORB  /
    DATA KQNAM/5,1,1,-1,-1,3,3,0,0,0,0,0,1/

    CALL NAMEIO(IR,JRET,SGUESS,NNAM,QNAM,KQNAM, &
                 GUESS,NORB,NORDER,IORDER,JORDER,TOLZ,TOLE,MIX, &
                 PRTMO,PUNMO,SYMDEN,PURIFY,INSORB, &
            0, &
        0,0,0,0,0,  0,0,0,0,0, &
        0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0, &
        0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0)

    nmo_ = norb
    
  end subroutine QChem_nameio
  subroutine QChem_set_mo0(sc0)
    double precision :: sc0(:,:)
    if(nao_*nmo_.eq.0) then
       throw_err("nao_ or nmo_ are illegal value", 1)
    end if
    if(nao_.ne.size(sc0,1) .or. nmo_.ne.size(sc0,2)) then
       throw_err("sc0 is invalid size", 1)
    end if
    use_mo_asgn_ = .true.
    allocate(mo_sc0_(nao_,nmo_),mo_asgn_(nmo_,nmo_))
    mo_sc0_(:,:) = sc0(:,:)
  end subroutine QChem_set_mo0
  subroutine QChem_set_ci0(c0)
    double precision :: c0(:,:)
    integer :: nstate, nwks
    nwks   = size(c0, 1)
    nstate = size(c0, 2)
    use_ci_asgn_ = .true.
    allocate(ci_c0_(nwks,nstate),ci_asgn_(nstate,nstate))
    ci_c0_(:,:) = c0(:,:)
  end subroutine QChem_set_ci0
  ! == Utils ==
  subroutine mm(t,a,b,c)
    ! matrix x matrix product.
    implicit none
    character,intent(in) :: t
    real(kind(0d0)),intent(in) :: a(:,:),b(:,:)
    real(kind(0d0)),intent(out) :: c(:,:)
    real(kind(0d0)),parameter :: ONE=1d0,ZERO=0d0
    integer :: m,n,k
    m = SIZE(c(:,1))
    n = SIZE(c(1,:))
    k = SIZE(b(:,1))
    call DGEMM(t,'N',m,n,k,ONE,a,k,b,k,ZERO,c,m)
  end subroutine mm
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
  ! == For inserting gamess code --
  subroutine QChem_get_csfh(ht, n, eshift)
    integer,intent(in) :: n
    real(kind(0d0)),intent(in) :: ht(:), eshift
    integer I, J, idx

    get_hIJ_ = .true.
    
    select case(type_hmat_)
    case("dense")
       
       if(allocated(hIJ_)) then
          if( size(hIJ_,1)<n) then
             deallocate(hIJ_)
             !write(*,*) "QChem_get_csfh n1 failed"
             !write(*,*) "size(hIJ_), n:", size(hIJ_,1), size(hIJ_,2), n
             !stop
          end if
       end if
       
       if(.not.allocated(hIJ_)) then
          allocate(hIJ_(n,n))
       end if
       
       idx = 1
       do J = 1, n
          do I = 1, J
             !if(size(ht)<idx) then
             !   write(*,*) "QChem_get_csfh idx failed"
             !   write(*,*) "size(ht):", size(ht)
             !   write(*,*) "idx:", idx
             !   stop
             !end if
             hIJ_(I,J) = ht(idx)
             if(I.eq.J) then
                hIJ_(I,J) = hIJ_(I,J) + eshift
             else
                hIJ_(J,I) = hIJ_(I,J)
             end if
             idx = idx+1          
          end do
       end do
    case default
       throw_err("not impled",1)
    end select

    write(*,*) "QChem_get_csfh end"
    
  end subroutine QChem_get_csfh
  subroutine QChem_assign_mo() ! <- WFN:gamess.src
    use Mod_asgn, only : asgn_asgn
    ! Save the MO information,
    ! that is called from the GAMESS regular prcess.
    double precision, allocatable :: c_nao(:,:), c(:,:), e_nao(:), e(:)

    write(*,*) "QChem_assign_mo begin"
    if(.not.use_mo_asgn_) then
       write(*,*) "no assign "
       write(*,*) "QChem_assign_mo end"
       return
    end if    

    write(*,*) "assign begin"    
    
    allocate(c_nao(nao_,nao_), c(nao_,nmo_), e(nmo_), e_nao(nao_))
    call QChem_daread2(15,nao_**2,c_nao)
    call QChem_daread(17,nao_,e_nao)

    e(:) = e_nao(1:nmo_)
    c(:,:) = c_nao(:,1:nmo_)
    call asgn_asgn(mo_asgn_,c,e,mo_sc0_)
    e_nao(1:nmo_) = e(:)
    c_nao(:,1:nmo_) = c(:,:)
    
    call QChem_dawrite2(15,nao_**2,c_nao)
    call QChem_dawrite(17,nao_,e_nao)

    done_mo_asgn_ = .true.

    write(*,*) "QChem_assign_mo end"

  end subroutine QChem_assign_mo
  subroutine QChem_assign_ci()
    use Mod_asgn, only : asgn_asgn
    double precision :: t1, t2
    integer nstate, nwks, kstat    
    double precision, allocatable :: c(:,:)
    integer, parameter :: MXRT=100
    double precision :: ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
    COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP    
    integer NFT11, NFT12, NFT13, NFT14, NFT15, NFT16, IDAF20,  NEMEMX
    COMMON /CIFILS/ NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX
    integer :: IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)

    write(*,*) "QChem_assign_ci begin"
    throw_err("Do not use this subroutine. Instead, call asgn_asgn directory in naewdy.f90", 1)

    if(use_ci_asgn_) then
       write(*,*) "assign begin"

       ! -- read CI vector --
       write(*,*) "read CI vector"
       call SEQREW(NFT12)
       READ (NFT12) nstate, nwks,t1,t2
       write(*,*) "nstate=", nstate
       write(*,*) "nwks=", nwks
       allocate(c(NWKS,NSTATE))
       do kstat = 1, NSTATE
          call SQREAD(NFT12,c(:,kstat),NWKS)
       end do
       call DIRCLO(NFT12,'DELETE')

       ! -- read CI energy --
       write(*,*) "read CI energy"
       call DAREAD(IDAF, IODA, ENUCR, MXRT+15, 2, 0)

       ! -- assign --
       write(*,*) "assign"
       call asgn_asgn(ci_asgn_, c, estate(1:nstate), ci_c0_); check_err()

       ! -- write CI vector --
       write(*,*) "write CI vector"
       CALL SEQOPN(NFT12,'CIVECTR','UNKNOWN',.FALSE.,'UNFORMATTED')
       WRITE(NFT12) nstate,nwks,t1,t2
       do kstat = 1, nstate
          call sqwrit(nft12, c(:,kstat), nwks)
       end do

       ! -- write CI energy --
       write(*,*) "write CI energy"
       call DAWRIT(IDAF, IODA, ENUCR, MXRT+15, 2, 0)

       done_ci_asgn_ = .true.
       
    else
       write(*,*) "no assign"
    end if
    
    write(*,*) "QChem_assign_ci end"
    
  end subroutine QChem_assign_ci
  ! == common interface ==
  subroutine QChem_calc
    use Mod_gamess_common, only : NWKS
    use Mod_math, only : dump_dmat
    integer, parameter :: MXRT=100
    double precision ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
    COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
    double precision GUESS, CIDRT
    DATA CIDRT/8HCIDRT   /
    character(10) :: type_calc
    type_calc = "energx"
    !type_calc = "direct"
    
    call aij_reset; check_err()
    get_hIJ_ = .false.

    select case(type_calc)
    case("direct")
       call ONEEI
       call JANDK
       call GUESMO(GUESS)
       call calc_orbital
       call QChem_assign_mo

       if(done_mo_asgn_) then
          call open_w(261, "asgn_mo.csv"); check_err()
          call dump_dmat(mo_asgn_, 261, 0.1d0); check_err()
          close(261)
       end if
       
       call DRTGEN(-5,CIDRT)
       call TRFMCX(-5,0,0,0,.FALSE.,.TRUE.,&
            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
            .FALSE.,.FALSE.,.true.)
       CALL GUGSRT(0,.false.)
       CALL GUGAEM(0)
       
       write(*,*) "QChem.gugaem begin"
       CALL GUGADG(0)
       
       write(*,*) "QChem.gugadg begin"
       !       call QChem_assign_ci
       CALL GUGADM(0)
       !       call GUG2DM(0)
    case("energx")
       call ENERGX
    case default
       throw_err("not supported",1)
    end select

    write(*,*) "end of QChem.calc"
    if(.not.get_hIJ_) then
       throw_err("failed to get H_CSF.",1)
    end if
    
    call aij_calc(ia_, va_, na_); check_err()

    NWKS_ = NWKS
    NSTATE_ =  INT(STATN)

    !call SEQREW(NFT12)
    !READ (NFT12) NSTATE_,NWKS_,tit1,tit2
    !call DIRCLO(NFT12,'KEEP')    
    
  end subroutine QChem_calc
  subroutine calc_orbital
    implicit none
    double precision RRSHFT,EXTTOL,DMPTOL,VSHTOL
    integer IEXTIN
    COMMON /ACONV / RRSHFT,EXTTOL,DMPTOL,VSHTOL,IEXTIN      
    double precision ACURCY,EN,ETOT,EHF,EHF0,DIFF
    integer ITER,ICALP,ICBET
    COMMON /CONV  / ACURCY,EN,ETOT,EHF,EHF0,DIFF,ITER,ICALP,ICBET
    integer MXATM, MXRT
    PARAMETER (MXATM=2000, MXRT=100)
    double precision ENUCR,EELCT,ETOT2,SZ2,SZZ2,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
    COMMON /ENRGYS/ ENUCR,EELCT,ETOT2,SZ2,SZZ2,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
    double precision E, EG
    COMMON /FUNCT / E,EG(3,MXATM)
    integer NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
    double precision ZAN,C
    COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
         ZAN(MXATM),C(3,MXATM),IAN(MXATM)
    integer IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
    double precision zero, ENUC
    PARAMETER (ZERO=0.0D+00)

!     convergence condition
    IEXTIN = 4
    EXTTOL = 1.0D-03
    DMPTOL = 1.0D-04
    VSHTOL = 0.4D+00
    call TRMAT
      
!     init energy
    E = ZERO
    ETOT = ZERO
    EN=ENUC(NAT,ZAN,C)
    CALL RHFCL
    
!     save energies
    E = ETOT
    ENUCR = EN
    EELCT = EHF
    ETOT2 = ETOT
    ECORE = ZERO
    ESCF  = ETOT
    CALL DAWRIT(IDAF,IODA,ENUCR,MXRT+15,2,0)
  end subroutine calc_orbital
  subroutine QChem_daread(id, n, data)
    integer, intent(in) :: id, n
    double precision, intent(out) :: data(:)
    integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)

    call DAREAD(IDAF, IODA, data, n, id, 0)
    
  end subroutine QChem_daread
  subroutine QChem_daread2(id, n, data)
    integer, intent(in) :: id, n
    double precision, intent(out) :: data(:,:)
    integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)

    call DAREAD(IDAF, IODA, data, n, id, 0)
    
  end subroutine QChem_daread2
  subroutine QChem_dawrite(id, n, data)
    integer, intent(in) :: id, n
    double precision, intent(out) :: data(:)
    integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
    call DAWRIT(IDAF,IODA,data,n,id,0)
  end subroutine QChem_dawrite
  subroutine QChem_dawrite2(id, n, data)
    integer, intent(in) :: id, n
    double precision, intent(out) :: data(:,:)
    integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
    call DAWRIT(IDAF,IODA,data,n,id,0)
  end subroutine QChem_dawrite2
  subroutine QChem_dump_aij(ifile)
    integer, intent(in) :: ifile
    integer :: i, icsf, jcsf, imo, jmo
    double precision :: aij
    write(ifile, '("I,J,i,j,val")') 
    do i = 1, na_
       icsf = ia_(1,i)
       jcsf = ia_(2,i)
       imo  = ia_(3,i)
       jmo  = ia_(4,i)
       aij  = va_(i)
       write(ifile,'(I10,",",I10,",",I0,",",I10,",",f20.10)') icsf, jcsf, imo, jmo, aij
    end do
  end subroutine QChem_dump_aij
  subroutine QChem_dump_Hcsf(ifile, eps)
    use Mod_math, only : dump_dmat
    integer, intent(in) :: ifile
    double precision, intent(in) :: eps

    if(nwks_.eq.0) then
       throw_err("nwks is illegal value", 1)
    end if

    call dump_dmat(hIJ_, ifile, eps); check_err()

  end subroutine QChem_dump_Hcsf
end module Mod_QChem
