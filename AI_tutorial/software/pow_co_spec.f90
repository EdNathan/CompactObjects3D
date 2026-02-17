include 'a_modules.f90'
include 'a_subroutines.f90'

program pow_co_spec
! Extract a light curve for IXPE event mode data and sort into segments
!
! Mac Studio
! gfortran -L/opt/homebrew/Cellar/cfitsio/4.6.2/lib -lcfitsio pow_co_spec.f90 
!
! iMac
! gfortran -L/usr/local/Cellar/cfitsio/4.4.1/lib -lcfitsio pow_co_spec.f90 
!
! Laptop
! gfortran -L/opt/homebrew/Cellar/cfitsio/4.2.0/lib -lcfitsio pow_co_spec.f90 
!
! General
! heainit
! gfortran -L$HEADAS/lib -lcfitsio pow_co_spec.f90
  
  use shared_arrays
  implicit none
  integer pilo,pihi,nseg
  integer evtunit1,evtunit2,evtunit3
  integer i,j
  double precision dt
  double precision tstart,tstop
  character (len=500) evtfile1,evtfile2,evtfile3
  integer nn
  real Elo,Ehi
  double precision tlive,ton,Telapse,Tseg,tbeg,tend,total_exp
  real, allocatable :: P(:),dP(:),Pco(:),dPco(:)
  real rate,wn,df,f,c,fbin,dfbin
  real, allocatable :: Pbin(:),dPbin(:),Pcobin(:),dPcobin(:)
  integer config
  real wn_est,tau,rate_in,wiggles

! Input parameters -------------------------------------
  Elo  = 2.0            !lower end of energy range
  Ehi  = 8.0            !upper end of energy range
  dt   = 1.0/64.0       !duration of each time bin
  nseg = 2**11          !number of time bins in each segment
  c    = 1.1            !geometric rebinning constant
  tbeg = 0.0            !user-defined start time (secs since MJDref)
  tend = 1e30           !user-defined stop  time (secs since MJDref)
  config = 3            !the DU that is used as the reference band
  tau    = 1.24e-3      !deadtime (s)
  !Input event files:
  evtfile1 = '../CygX1data/03010101/event_l2/ixpe03010101_det1_evt2_v01.fits'
  evtfile2 = '../CygX1data/03010101/event_l2/ixpe03010101_det2_evt2_v01.fits'
  evtfile3 = '../CygX1data/03010101/event_l2/ixpe03010101_det3_evt2_v01.fits'
! ------------------------------------------------------

! Open output files
  open(79,file='Praw.dat')
  open(81,file='Pbin.dat')
  open(82,file='fPbin.dat')
  open(87,file='xPbin.dat')
  
! Derived quantities
  Tseg = dt * nseg
  write(*,*)"nseg=",nseg
  write(*,*)"Tseg (s) =",Tseg
  df = 1.0 / Tseg
  pilo = max( 0 , ceiling(Elo/0.04)+1 )
  pihi = min( 374 , floor(Ehi/0.04) )

! Open event files with read only access
  call open_evts(evtunit1,evtfile1,evtunit2,evtfile2,evtunit3,evtfile3)
  
! Get number of events in the event files
  call getnevts(evtunit1,nevts1)
  call getnevts(evtunit2,nevts2)
  call getnevts(evtunit3,nevts3)
  write(*,*)"Events in each DU:",nevts1,nevts2,nevts3
  write(*,*)"Total number of events=",nevts1+nevts2+nevts3

! Get tstart and tstop
  call get_tstart(evtunit1,tstart,tstop,tlive,ton)
  tstart = max( tstart , tbeg )
  tstop  = min( tstop  , tend )

! Get a master GTI list
  call GTIwrangler(evtunit1,evtunit2,evtunit3,tstart,tstop,Telapse)

! Calculate number of bins in simple light curve
  nn = ceiling( Telapse / dt )
  write(*,*)"Number of time bins in the observation=",nn

! Sort GTIs into segments
  call seg_wrangler(tstart,dt,nseg,total_exp)

! Extract full band light curves from different combinations of DUs
  call fullbandlc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,tstart,dt,config)

! Make power spectrum of summed DU light curve
  allocate(  P(nseg/2) )
  allocate( dP(nseg/2) )
  call mypow(nseg,lc,nn,real(dt),mm,segend,1,P,dP,rate,wn)
  write(*,*)"Total count rate (c/s) = ",rate
  rate_in = rate / ( 1.0 - tau * rate )
  write(*,*)"Total intrinsic count rate (c/s) = ",rate_in
  rate_in = rate_in / 3.0
  write(*,*)"Intrinsic count rate per DU (c/s) = ",rate_in
  
! Plot power spectrum and subtract Poisson noise
  write(79,*)"read serr 1 2"
  write(79,*)"skip on"
  do i = 1,nseg/2
     f = i * df
     wn_est = wn * wiggles(f,tau,rate_in)
     write(79,*)f,0.5*df,P(i),dP(i),wn_est
     P(i) = P(i) - wn_est
  end do
  
! Calculate co-spectrum
  allocate(  Pco(nseg/2) )
  allocate( dPco(nseg/2) )
  call myco(nseg,lcsub,lcref,nn,real(dt),mm,segend,1,Pco,dPco)
  
! Set up binning scheme
  call geobin_scheme(nseg,df,c)

! Bin power spectrum
  allocate( Pbin(nf) )
  allocate( dPbin(nf) )
  call geobin_bin(nseg,P,dP,Pbin,dPbin)
  
! Bin co-spectrum
  allocate( Pcobin(nf) )
  allocate( dPcobin(nf) )
  call geobin_bin(nseg,Pco,dPco,Pcobin,dPcobin)

! Plot binned co-spectrum
  write(81,*)"skip on"
  write(81,*)"read serr 1 2 3"
  do i = 1,nf
     fbin  = 0.5 * ( far(i) + far(i-1) )
     dfbin = 0.5 * ( far(i) - far(i-1) )
     write(81,*)fbin,dfbin,Pcobin(i),dPcobin(i),Pbin(i),dPbin(i)
     write(87,*)fbin-dfbin,fbin+dfbin,2.0*Pcobin(i)*dfbin,2.0*dPcobin(i)*dfbin
  end do
  
! Plot binned co-spectrum and power spectrum in f * P form
  write(82,*)"skip on"
  write(82,*)"read serr 1 2 3"
  do i = 1,nf
     fbin  = 0.5 * ( far(i) + far(i-1) )
     dfbin = 0.5 * ( far(i) - far(i-1) )
     write(82,*)fbin,dfbin,fbin*Pcobin(i),fbin*dPcobin(i),fbin*Pbin(i),fbin*dPbin(i)
  end do

! Write summary of outputs
  write(*,*)"-------------------------------------"
  write(*,*)"Outputs:"
  write(*,*)"Praw.dat: raw unbinned power with Poisson noise"
  write(*,*)"Pbin.dat: binned noise subtracted copower (black) and power (red)"
  write(*,*)"fPbin.dat: binned noise subtracted f*copower (black) and f*power (red)"
  write(*,*)"xPbin.dat: to convert this into an xspec file xPbin.pha use..."
  write(*,*)"...flx2xsp infile=xPbin.dat phafile=xPbin.pha rspfile=xPbin.rsp"
  write(*,*)"-------------------------------------"
  
! Write plot commands
  write(79,*)"log"
  write(79,*)"la y Power (rms\u2\d/Hz)"
  write(79,*)"la x Frequency (Hz)"
  write(81,*)"log"
  write(81,*)"la y Power (rms\u2\d/Hz)"
  write(81,*)"la x Frequency (Hz)"
  write(82,*)"log"
  write(82,*)"la y Frequency x Power (rms\u2\d)"
  write(82,*)"la x Frequency (Hz)"

  
end program pow_co_spec



