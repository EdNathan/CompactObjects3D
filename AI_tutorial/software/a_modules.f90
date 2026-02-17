module shared_arrays
  integer ngti1,ngti2,ngti3,ngti,mm
  double precision, allocatable :: gtistart1(:),gtiend1(:)
  double precision, allocatable :: gtistart2(:),gtiend2(:)
  double precision, allocatable :: gtistart3(:),gtiend3(:)
  double precision, allocatable :: gtistart(:),gtiend(:)
  integer, allocatable :: segend(:),segment_no(:)
  real, allocatable :: lc(:),lcsub(:),lcref(:),lcphi(:,:),lcphi3(:,:),phihist(:)
  real, allocatable :: lc_lo_12(:),lc_hi_3(:),lc_lo_3(:),lc_hi_12(:),dphihist(:)
  real, allocatable :: Ilc(:),varIlc(:),qlc(:),varQlc(:),ulc(:),varUlc(:)
  integer nevts1,nevts2,nevts3
  integer nf,nevts_sub,nevts_ref
  integer, allocatable :: np(:),iar(:)
  real, allocatable :: far(:)
  double precision pi
  integer nchn
  real, allocatable :: ECHN(:),AEmuEi(:),Aeffi(:),muEi(:)
  parameter (pi=3.141592653589793)
end module shared_arrays
