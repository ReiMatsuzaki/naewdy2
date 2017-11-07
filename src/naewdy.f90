#include "macros.fpp"
#include "err_handle.f90"
#include "istream.f90"
#include "fjson.f90"
#include "sys.f90"
#include "const.f90"
#include "math.f90"
#include "qchem.f90"

module Mod_GugaPrint
  use Mod_err_handle
  use Mod_QChem
  implicit none
  private
  public run
contains
  ! ==== Driver ====
  subroutine run
    use Mod_fjson
    use Mod_sys, only : mkdirp_if_not
    use Mod_nshel, only : Nshel_dump
    use Mod_math
    use Mod_asgn, only : asgn_asgn
    character(20) :: fn
    integer :: ifile = 1234
    double precision, allocatable :: T(:), M(:,:), e(:), c(:,:), eci(:), cmo(:,:)
    double precision :: tit1, tit2
    integer i0, i1, kstat
    integer NFT11, NFT12, NFT13, NFT14, NFT15, NFT16, IDAF20,  NEMEMX
    COMMON /CIFILS/ NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX
    type(object) :: o
    
    write(*,*)
    write(*,*) "    --------------"
    write(*,*) "NAEWDY v2.00"
    write(*,*) "    --------------"
    write(*,*)

    ! -- initialize --
    call err_handle_begin
    call QChem_new; check_err()
    call new_json
    call mkdirp_if_not("out"); check_err()
    
    ! -- run gamess --    
    call QChem_calc; check_err()    
    allocate(T(nao_*(nao_+1)/2), M(nao_,nao_))
    
    call object_set_i(o, "nao", nao_); check_err()
    call object_set_i(o, "nmo", nmo_); check_err()
    call object_set_i(o, "nwks", NWKS_); check_err()
    call object_set_i(o, "nstate", nstate_); check_err()

    ! -- AO --
    fn = "out/"//"nshel.json"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call Nshel_dump(ifile); check_err()
    close(ifile)

    ! -- AO matrix --        
    call QChem_daread(12, nao_*(nao_+1)/2, T(:)); call t2s(T,M)
    fn = "out/"//"s.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dmat(M, ifile); check_err()
    close(ifile)

    call QChem_daread(13, nao_*(nao_+1)/2, T(:)); call t2s(T,M)
    fn = "out/"//"t.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dmat(M, ifile); check_err()
    close(ifile)

    call QChem_daread(11, nao_*(nao_+1)/2, T(:)); call t2s(T,M)
    fn = "out/"//"h.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dmat(M, ifile); check_err()
    close(ifile)

    ! -- MO --
    allocate(e(nmo_), cmo(nao_,nmo_))
    call QChem_daread2(15, nao_*nmo_, cmo(:,:)); check_err()    
    fn = "out/"//"cmo.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dmat(cmo, ifile); check_err()
    close(ifile)

    call QChem_daread(17, nmo_, e(:)); check_err()
    fn = "out/"//"emo.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dvec(e, ifile); check_err()
    close(ifile)

    ! -- CSF --
    fn = "out/"//"aij.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call QChem_dump_aij(ifile); check_err()
    close(ifile)

    fn = "out/"//"hcsf.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call QChem_dump_Hcsf(ifile, 1.0d-10); check_err()

    ! -- CI --    
    call SEQREW(NFT12)
    READ (NFT12) i0, i1,tit1,tit2    
    allocate(c(NWKS_,NSTATE_), eci(NSTATE_))
    do kstat = 1, NSTATE_
       call SQREAD(NFT12,c(:,kstat),NWKS_)
    end do
    call DIRCLO(NFT12,'DELETE')        

    call read_ci_energy(eci(:)); check_err()

    if(use_ci_asgn_) then       
       call asgn_asgn(ci_asgn_, c, eci, ci_c0_); check_err()
       done_ci_asgn_ = .true.
    end if

    fn = "out/"//"cci.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dmat(c, ifile)
    close(ifile)
    
    fn = "out/"//"eci.csv"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call dump_dvec(eci, ifile); check_err()
    close(ifile)

    ! -- asgn --
    if(done_mo_asgn_) then
       fn = "out/asgn_mo.csv"; ifile=ifile+1
       call open_w(ifile, fn); check_err()
       call dump_dmat(mo_asgn_, ifile, 0.1d0); check_err()
       close(ifile)
    end if

    if(done_ci_asgn_) then
       fn = "out/asgn_ci.csv"; ifile=ifile+1
       call open_w(ifile, fn); check_err()
       call dump_dmat(ci_asgn_, ifile, 0.1d0); check_err()
       close(ifile)
    end if    

    ! -- other --
    fn = "out/"//"common.json"; ifile=ifile+1
    call open_w(ifile, fn); check_err()
    call object_dump(o, ifile); check_err()
    close(ifile)
    
  end subroutine run
  ! ==== Init ====
  subroutine new_json
    use Mod_fjson
    character(100) :: prev_out, comment
    integer, parameter :: ifile = 2345121
    type(value) :: v
    type(object) :: o, oo
    logical :: use_mo, use_ci
   
    
    call loads_json_file("naewdy.in.json", ifile, v); check_err()
    call value_get_o(v, o); check_err()

    call object_get_s(o, "comment", comment); check_err()
    write(*,*) "comment:", comment
    
    if(object_exist(o, "assign")) then
       call object_get_o(o, "assign", oo); check_err()
       call object_get_s(oo, "dir", prev_out); check_err()
       write(*,*) "dir:", prev_out       
       call object_get_b(oo, "use_mo", use_mo); check_err()
       write(*,*) "use_mo:", use_mo
       call object_get_b(oo, "use_ci", use_ci); check_err()
       write(*,*) "use_ci:", use_ci
       call object_delete(oo); check_err()
       
       if(use_mo) then
          call new_set_mo0(prev_out); check_err()
       end if

       if(use_ci) then
          call new_set_ci0(prev_out); check_err()
       end if
    end if
    
    call object_delete(o); check_err()
    call value_delete(v); check_err()
    
  end subroutine new_json
  subroutine new_set_mo0(prev_out)
    use Mod_fjson
    character(*), intent(in) :: prev_out
    character(100) :: fn
    integer, parameter :: ifile = 12325
    type(object) :: o
    type(value)  :: v
    integer :: idx, i, j, nao, nmo
    double precision :: val
    double precision, allocatable :: s(:,:), c0(:,:), sc0(:,:)

    write(*,*) "new_set_mo0 begin"

    fn = trim(prev_out) // "/common.json"
    write(*,*) "reading file:", fn
    call loads_json_file(fn, ifile, v); check_err()
    call value_get_o(v, o); check_err()
    call object_get_i(o, "nao", nao); check_err()
    call object_get_i(o, "nmo", nmo); check_err()
    write(*,*) "nao:", nao
    write(*,*) "nmo:", nmo
    call value_delete(v);
    call object_delete(o);
    close(ifile)

    allocate(c0(nao, nmo), s(nao,nao), sc0(nao,nmo))
    
    fn = trim(prev_out) // "/cmo.csv"
    write(*,*) "reading file:", fn
    call open_r(ifile, trim(fn)); check_err()    
    read(ifile,*)
    do idx = 1, nao_*nmo_
       read(ifile,*) i, j, val
       c0(i,j) = val
    end do
    close(ifile)

    fn = trim(prev_out) // "/s.csv"
    call open_r(ifile, fn); check_err()
    read(ifile,*)
    do idx = 1, nao_*nao_
       read(ifile,*) i, j, val
       s(i,j) = val
    end do
    close(ifile)

    sc0(:,:) = matmul(s(:,:), c0(:,:))

    call QChem_set_mo0(sc0)

    deallocate(s,c0,sc0)

  end subroutine new_set_mo0
  subroutine new_set_ci0(prev_out)
    use Mod_fjson
    character(*), intent(in) :: prev_out
    character(100) :: fn
    integer, parameter :: ifile = 12325
    type(object) :: o
    type(value)  :: v
    integer :: idx, i, j, nstate, nwks
    double precision :: val
    double precision, allocatable :: c0(:,:)

    write(*,*) "new_set_ci0 begin"

    fn = trim(prev_out) // "/common.json"
    write(*,*) "reading file:", fn
    call loads_json_file(fn, ifile, v); check_err()
    call value_get_o(v, o); check_err()
    call object_get_i(o, "nstate", nstate); check_err()
    write(*,*) "nstate:", nstate
    call object_get_i(o, "nwks", nwks); check_err()
    write(*,*) "nwks:", nwks
    call value_delete(v);
    call object_delete(o);
    close(ifile)

    allocate(c0(nwks, nstate))
    
    fn = trim(prev_out) // "/cci.csv"
    write(*,*) "reading file:", fn
    call open_r(ifile, trim(fn)); check_err()    
    read(ifile,*)
    do idx = 1, nstate*nwks
       read(ifile,*) i, j, val
       c0(i,j) = val
    end do
    close(ifile)

    call QChem_set_ci0(c0); check_err()

    deallocate(c0)

    write(*,*) "new_set_ci0 end"
    
  end subroutine new_set_ci0
  ! ==== Utils ====
  subroutine read_ci_energy(e)
    integer, parameter :: MXRT=100
    double precision :: e(:)
    double precision :: ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
    COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2, &
         VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
    integer :: IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
    COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)

    write(*,*) "read_ci_energy begin"
    call DAREAD(IDAF, IODA, ENUCR, MXRT+15, 2, 0)

    e(:) = ESTATE(1:size(e))

  end subroutine read_ci_energy
end module Mod_GugaPrint

subroutine naewdyx
  use Mod_GugaPrint, only : run
  call run
end subroutine naewdyx

