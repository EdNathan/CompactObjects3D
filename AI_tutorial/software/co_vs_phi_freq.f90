include 'a_modules.f90'
include 'a_subroutines.f90'

program co_vs_phi
! Caculate rms and phase as a function of frequency and
! modulation angle. This code will stack over multiple observations.
!
! Fits the sine function to output G_Q and G_U as a function of frequency 
!
! gfortran -L$HEADAS/lib -lcfitsio co_vs_phi_freq.f90
!  
  use shared_arrays
  implicit none
  integer pilo,pihi,nseg
  integer evtunit1,evtunit2,evtunit3
  integer i,j,k
  double precision dt
  double precision tstart,tstop
  character (len=500) evtfile1,evtfile2,evtfile3,evtlist
  integer nn,phibins,config,status
  real Elo,Ehi,f,Psn,phideg,dphideg,lagcyc,dlagcyc,c,delf
  double precision tlive,ton,Telapse,Tseg,tbeg,tend,dphi,Texp
  real, allocatable :: P(:),dP(:),Pco(:),dPco(:)
  real df,rate,wn,nd
  real, allocatable :: Pbin(:),dPbin(:),Pcobin(:),dPcobin(:)
  real, allocatable :: ReG(:),dReG(:),ImG(:),dImG(:)
  real, allocatable :: ReGphi(:,:),dReGphi(:,:),ImGphi(:,:),dImGphi(:,:),dG(:,:)
  real, allocatable :: rms(:,:),drms(:,:),lag(:,:),dlag(:,:)
  real, allocatable :: Pr(:),dPr(:),Prco(:),dPrco(:),Ps(:,:),dPs(:,:)
  real Psub,dPsub
  real, allocatable :: ReGav(:),ImGav(:),dGav(:),phi(:)
  real, allocatable :: ReGtot(:,:),ImGtot(:,:),r_sub(:),r_subtot(:),dr_subtot(:)
  real, allocatable :: ReG_mod(:),ImG_mod(:),hmod(:),r_tot(:),dr_tot(:)
  real, allocatable :: Pr_tot(:),Pco_tot(:),Ps_tot(:,:)
  real, allocatable :: ReG_Q(:),ImG_Q(:),ReG_U(:),ImG_U(:)
  real, allocatable :: dReG_Q(:),dImG_Q(:),dReG_U(:),dImG_U(:)
  integer Nobs,obs,listunit,mmtot
  real r_ref,r_reftot,Ttot
  integer mmax,ma
  parameter (mmax=20)
  integer ia(mmax)
  real chisq,da(mmax),dyda(mmax),a_re(mmax),a_im(mmax),chisq_full,lag_mod
  real a0(mmax),qav,uav,a_null(mmax),chisq_null,rms_mod,r_all,r_alltot
  real da_re(mmax),da_im(mmax),null,dnull
  double precision N1,Q1,U1,N1tot,Q1tot,U1tot
  
! Input parameters -------------------------------------  
  Elo      = 2.0        !lower end of energy range
  Ehi      = 8.0        !upper end of energy range
  dt       = 1.0/64.0   !duration of each time bin
  nseg     = 2**10      !number of time bins in each segment
  tbeg     = 0.0        !user-defined start time (secs since MJDref)
  tend     = 1e30       !user-defined stop  time (secs since MJDref)
  phibins  = 20         !number of modulation angle bins
  config   = 3          !DU configuration (cross the other two DUs with DU config)
  c        = 1.3        !geometric rebinning constant
  Nobs     = 1          !number of of observations to stack over
  evtlist  = '../CygX1data/cygx1_hard_obsids.txt' !File with list of event files (3 per obsid)
! ------------------------------------------------------
  
! Derived quantities
  Tseg = dt * nseg
  write(*,*)"nseg=",nseg
  write(*,*)"Tseg=",Tseg
  df = 1.0 / Tseg
  pilo = ceiling(Elo/0.04)+1
  pihi = floor(Ehi/0.04)
  write(*,*)"pilo=",pilo
  write(*,*)"pihi=",pihi
  allocate( phi(phibins) )
  dphi = pi/dble(phibins)
  dphideg = dphi * 180.0 / pi
  do j = 1,phibins
     phi(j) = -0.5*pi + (dble(j)-0.5) * dphi
  end do

! Open files
  open(77,file='ReImGU_vs_freq.dat')
  open(79,file='ReImGQ_vs_freq.dat')
  open(78,file='P_vs_freq.dat')
  
! Set up binning scheme
  call geobin_scheme(nseg,df,c)
  write(*,*)"nf=",nf

! Allocate arrays with constant size
  allocate(  P(nseg/2) )
  allocate( dP(nseg/2) )  
  allocate(  Pr(nf) )
  allocate( dPr(nf) )
  allocate(  Pr_tot(nf) )
  allocate(  Prco(nf) )
  allocate( dPrco(nf) )
  allocate(  Pco_tot(nf) )
  
  write(78,*)"read serr 1 2"
  write(78,*)"skip on"
  
  allocate(  ReG(nseg/2) )
  allocate( dReG(nseg/2) )
  allocate(  ImG(nseg/2) )
  allocate( dImG(nseg/2) )  
  allocate(  ReGphi(nf,phibins) )
  allocate( dReGphi(nf,phibins) )
  allocate(  ImGphi(nf,phibins) )
  allocate( dImGphi(nf,phibins) )
  allocate( dG(nf,phibins) )
  allocate(  lag(nf,phibins) )
  allocate( dlag(nf,phibins) )
  allocate(  rms(nf,phibins) )
  allocate( drms(nf,phibins) )
  allocate(  Pco(nseg/2) )
  allocate( dPco(nseg/2) )
  
  allocate(  Ps(nf,phibins) )
  allocate( dPs(nf,phibins) )
  allocate( Ps_tot(nf,phibins) )

  allocate( ReGtot(nf,phibins) )
  allocate( ImGtot(nf,phibins) )
  
  allocate( r_sub(0:phibins) )
  allocate( r_subtot(0:phibins) )
  allocate( dr_subtot(0:phibins) )
  allocate( r_tot(phibins) )
  allocate( dr_tot(phibins) )

  allocate( ReG_Q(nf) )
  allocate( ImG_Q(nf) )
  allocate( ReG_U(nf) )
  allocate( ImG_U(nf) )
  allocate( dReG_Q(nf) )
  allocate( dImG_Q(nf) )
  allocate( dReG_U(nf) )
  allocate( dImG_U(nf) )
  
! Initialize running sums
  mmtot    = 0
  Pr_tot   = 0.0
  Pco_tot  = 0.0
  Ps_tot   = 0.0
  ReGtot   = 0.0
  ImGtot   = 0.0
  r_reftot = 0.0
  r_subtot = 0.0
  N1tot    = 0.0
  Q1tot    = 0.0
  U1tot    = 0.0
  r_alltot = 0.0
  
! Open event list file
  status = 0
  call ftgiou(listunit,status)
  open(listunit,file=evtlist)
  
! Loop through all observations
  do obs = 1,nobs

     write(*,*)"--------------------------------------------"
     write(*,*)"Observation number = ",obs
     
     !Read the names of the event files
     read(listunit,'(a)')evtfile1
     read(listunit,'(a)')evtfile2
     read(listunit,'(a)')evtfile3
     
     !Open event files with read only access
     call open_evts(evtunit1,evtfile1,evtunit2,evtfile2,evtunit3,evtfile3)
  
     !Get number of events in the event files
     call getnevts(evtunit1,nevts1)
     call getnevts(evtunit2,nevts2)
     call getnevts(evtunit3,nevts3)
     write(*,*)"Total number of events=",nevts1+nevts2+nevts3

     !Calculate Q and U to put into the model
     call stokes(evtunit1,nevts1,pilo,pihi,N1,Q1,U1)
     
     !Get tstart and tstop
     call get_tstart(evtunit1,tstart,tstop,tlive,ton)
     tstart = max( tstart , tbeg )
     tstop  = min( tstop  , tend )
  
     !Get a master GTI list
     call GTIwrangler(evtunit1,evtunit2,evtunit3,tstart,tstop,Telapse)
  
     !Define more accurate values for tstart and tstop
     tstart  = gtistart(1)
     tstop   = gtiend(ngti)
     Telapse = tstop - tstart
  
     !Calculate number of bins in simple light curve
     nn = ceiling( Telapse / dt )

     !Sort GTIs into segments
     call seg_wrangler(tstart,dt,nseg,Texp)
     !write(*,*)"Number of realisations = ",mm*(ihi-ilo+1)
  
     !Extract full band light curves from different combinations of DUs
     write(*,*)"Calculating fullband light curve..."
     call fullbandlc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,tstart,dt,config)
     write(*,*)"...finished calculating fullband light curve"
     
     !Extract phibins light curves from two DUs
     write(*,*)"Calculating phi-binned light curves..."
     call multiphilc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,tstart,dt,config)
     write(*,*)"...finished calculating phi-binned light curves"
     
     !Calculate and correct for spurious polarisation
     call spur_corr(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,phi,nseg,config)

     !Calculate mean subject and reference band count rates
     call mean_rate(nseg,lcsub,nn,mm,segend,r_sub(0))
     write(*,*)"count rate (sub)=",r_sub(0)
     call mean_rate(nseg,lcref,nn,mm,segend,r_ref)
     write(*,*)"count rate (ref)=",r_ref
     r_all = r_sub(0) + r_ref

     !Adjust count rates to scale to all DUs
     lcphi = lcphi * r_all / r_sub(0)
     lcref = lcref * r_all / r_ref
     lcsub = lcsub * r_all / r_sub(0)
     
     !Make power spectrum of reference band light curve
     call mypow(nseg,lcref,nn,real(dt),mm,segend,2,P,dP,rate,wn)

     !Rebin the power spectrum
     call geobin_bin(nseg,P,dP,Pr,dPr)

     !Calculate co-spectrum
     call myco(nseg,lcsub,lcref,nn,real(dt),mm,segend,2,Pco,dPco)

     !Rebin the co-spectrum
     call geobin_bin(nseg,Pco,dPco,Prco,dPrco)
     
     do i = 1,nf
        f    = 0.5 * ( far(i) + far(i-1) )
        delf = 0.5 * ( far(i) - far(i-1) )
        write(78,*)f,delf,Prco(i),dPrco(i)
     end do
     write(78,*)"no no"

     do i = 1,nf
        f    = 0.5 * ( far(i) + far(i-1) )
        delf = 0.5 * ( far(i) - far(i-1) )
        write(78,*)f,delf,Pr(i),dPr(i)
     end do
     write(78,*)"no no"
     
     !Calculate a cross spectrum for each phi bin
     write(*,*)"Calculating cross spectra..."
     do j = 1,phibins
        !Calculate cross spectra
        call mycross(nseg,lcphi(:,j),lcref,nn,real(dt),mm,segend,2,ReG,dReG,ImG,dImG)
        !write(*,*)"phi range (rad) = ",phi(j)-0.5*dphi,phi(j)+0.5*dphi
        !Calculate mean count rate for each phi bin
        call mean_rate(nseg,lcphi(:,j),nn,mm,segend,r_sub(j))
        !Bin cross spectra
        call geobin_bin(nseg,ReG,dReG,ReGphi(:,j),dReGphi(:,j))
        call geobin_bin(nseg,ImG,dImG,ImGphi(:,j),dImGphi(:,j))
        !Calculate power spectrum
        call mypow(nseg,lcphi(:,j),nn,real(dt),mm,segend,2,P,dP,rate,wn)
        !Bin power spectrum
        call geobin_bin(nseg,P,dP,Ps(:,j),dPs(:,j))        
     end do
     write(*,*)"...finished calculating cross spectra"

     !Add to running sums
     mmtot    = mmtot    + mm
     Pr_tot   = Pr_tot   + mm * Pr
     Pco_tot  = Pco_tot  + mm * Prco
     Ps_tot   = Ps_tot   + mm * Ps
     ReGtot   = ReGtot   + mm * ReGphi
     ImGtot   = ImGtot   + mm * ImGphi
     r_reftot = r_reftot + mm * r_ref
     r_subtot = r_subtot + mm * r_sub
     N1tot    = N1tot + mm * N1
     Q1tot    = Q1tot + mm * Q1
     U1tot    = U1tot + mm * U1
     r_alltot = r_alltot + mm * r_all
     
     write(*,*)"--------------------------------------------"

     !Close event files
     close(evtunit1)
     call ftfiou(evtunit1, status)
     close(evtunit2)
     call ftfiou(evtunit2, status)
     close(evtunit3)
     call ftfiou(evtunit3, status)

     !Deallocate arrays
     deallocate( gtistart1 )
     deallocate( gtiend1 )
     deallocate( gtistart2 )
     deallocate( gtiend2 )
     deallocate( gtistart3 )
     deallocate( gtiend3 )
     deallocate( gtistart )
     deallocate( gtiend )
     deallocate( segend )
     deallocate( lc )
     deallocate( lcsub )
     deallocate( lcref )
     deallocate( lcphi )
     
  end do

! Finish averages
  Pr_tot    = Pr_tot   / real(mmtot)
  Pco_tot   = Pco_tot  / real(mmtot)
  Ps_tot    = Ps_tot   / real(mmtot)
  ReGtot    = ReGtot   / real(mmtot)
  ImGtot    = ImGtot   / real(mmtot)
  r_reftot  = r_reftot / real(mmtot)
  r_subtot  = r_subtot / real(mmtot)
  Ttot      = real( mmtot * nseg * dt )
  dr_subtot = sqrt( r_subtot / Ttot )
  N1tot     = N1tot / real(mmtot) / Ttot
  Q1tot     = Q1tot / real(mmtot) / Ttot
  U1tot     = U1tot / real(mmtot) / Ttot
  r_alltot  = r_alltot / real(mmtot)

! Transfer count rates to a suitable array
  do j = 1,phibins
     phideg    = phi(j) * 180.0 / pi
     r_tot(j)  = r_subtot(j)
     dr_tot(j) = sqrt( r_subtot(j) / Ttot * r_alltot / r_subtot(0) )
  end do

  
! Calculate errors and convert to rms and lag
  do j = 1,phibins
     do i = 1,nf
        nd = np(i) * mmtot
        !Calculate error from Ingram (2019) formula
        call ingram2019_err_b0(Pr_tot(i),Ps_tot(i,j),Pco_tot(i),ReGtot(i,j),ImGtot(i,j),nd,dG(i,j))
        !Convert to rms and phase lag
        call rms_and_lag(ReGtot(i,j),ImGtot(i,j),dG(i,j),Pco_tot(i),far(i-1),far(i),rms(i,j),drms(i,j),lag(i,j),dlag(i,j))
        !Convert rms to fractional
        rms(i,j)  =  rms(i,j) / r_tot(j)
        drms(i,j) = drms(i,j) / r_tot(j)
     end do
  end do
  
! ======================================================================
! Fit sine models
! ======================================================================
! Set number of model parameters
  ma = 3

! Allocate arrays
  allocate(     hmod(phibins) )
  allocate( ReG_mod(phibins) )
  allocate( ImG_mod(phibins) )
  
! Fit I, Q U model to r_tot
  
  !Parameters
  ia   = 1
  a0(1) = r_alltot / real(phibins)         !I dphi/pi
  a0(2) = 3.*Q1tot / real(phibins)         !Q dphi/pi
  a0(3) = 3.*U1tot / real(phibins)         !U dphi/pi  
  !Fit the model
  call dofit(phi,r_tot,dr_tot,phibins,a0,ia,ma,IQUsinfunc,chisq,da)
  write(*,*)"I, Q, U fit chisq/dof = ",chisq,"/",phibins-3
  !Write out the model
  do j = 1,phibins
     phideg = phi(j) * 180.0 / pi
     call IQUsinfunc(phi(j),a0,hmod(j),dyda,ma)
  end do
  qav = a0(2)/a0(1)
  uav = a0(3)/a0(1)
  write(*,*)"qav=",qav
  write(*,*)"uav=",uav
  write(*,*)"PA (deg) =",0.5*atan2(uav,qav)*180./pi

  write(79,*)"read serr 1 2 3"
  write(77,*)"read serr 1 2 3"
  
! Fit full model to cross spectra for each frequency bin
  do i = 1,nf
  
     !Fit for ReGtot
     ia      = 1
     a_re(1) = 100.0   !P_I dphi/pi
     a_re(2) = 4.0     !ReG_Q dphi/pi
     a_re(3) = 1.0     !ReG_U dphi/pi
     call dofit(phi,ReGtot(i,:),dG(i,:),phibins,a_re,ia,ma,IQUsinfunc,chisq,da_re)
     chisq_full = chisq
  
     !Fit for ImGtot
     ia      = 1
     ia(1)   = 0
     a_im(1) = 0.0     !ImP_I dphi/pi (zero)
     a_im(2) = 4.0     !ImG_Q dphi/pi
     a_im(3) = 1.0     !ImG_U dphi/pi
     call dofit(phi,ImGtot(i,:),dG(i,:),phibins,a_im,ia,ma,IQUsinfunc,chisq,da_im)
     chisq_full = chisq_full + chisq

     !Chi squared for full model
     write(*,*)"Full model:"
     write(*,*)"chisq/dof=",chisq_full,"/",2*phibins-5

     !Set G_Q and G_U
     ReG_Q(i)  =  a_re(2) * phibins
     dReG_Q(i) = da_re(2) * phibins
     ImG_Q(i)  =  a_im(2) * phibins
     dImG_Q(i) = da_im(2) * phibins
     ReG_U(i)  =  a_re(3) * phibins
     dReG_U(i) = da_re(3) * phibins
     ImG_U(i)  =  a_im(3) * phibins
     dImG_U(i) = da_im(3) * phibins

     !Write out
     f    = 0.5 * ( far(i) + far(i-1) )
     delf = 0.5 * ( far(i) - far(i-1) )
     write(79,*)f,delf,ReG_Q(i),dReG_Q(i),ImG_Q(i),dImG_Q(i)
     write(77,*)f,delf,ReG_U(i),dReG_U(i),ImG_U(i),dImG_U(i)
     
  end do


! Null hypothesis (PA and PD) constant:
! G(phi) = Pco/J + q*Pco/J * cos(2phi) + u*Pco/J * sin(2phi)
! G(phi) = 1/J [ Pco + G_Q cos(2phi) + G_U sin(2phi) ]
! =>  
! G_Q = q*Pco; G_U = u*Pco
! =>
! ReG_Q = q*Pco; ImG_Q = 0; ReG_Q = u*Pco; ImG_U = 0.

  write(79,*)"no no"
  write(77,*)"no no"
  
  do i = 1,nf

     f    = 0.5 * ( far(i) + far(i-1) )
     delf = 0.5 * ( far(i) - far(i-1) )

     null  = qav* Pco_tot(i)
     dnull = 0.0
     write(79,*)f,delf,null,dnull,0.0,0.0

     null  = uav* Pco_tot(i)
     dnull = 0.0
     write(77,*)f,delf,null,dnull,0.0,0.0
     
  end do  
  
! Write out commands
  write(77,*)"skip on"
  write(77,*)"log x"
  write(77,*)"la x Frequency (Hz)"
  write(77,*)"la y G\dU\u"
  write(77,*)"li s on 3,4"
  write(77,*)"co 1 on 3"
  write(77,*)"co 2 on 4"
  
  write(78,*)"log"
  write(78,*)"la x Frequency (Hz)"
  write(78,*)"la y Power (RMS\u2\d/Hz)"

  write(79,*)"skip on"
  write(79,*)"log x"
  write(79,*)"la x Frequency (Hz)"
  write(79,*)"la y G\dQ\u"
  write(79,*)"li s on 3,4"
  write(79,*)"co 1 on 3"
  write(79,*)"co 2 on 4"
  
  write(*,*)"----------------------------------------------------------"
  write(*,*)"Outputs:"
  write(*,*)"ReImGQ_vs_freq.dat: Black=ReG_Q, Red=ImG_Q"
  write(*,*)"ReImGU_vs_freq.dat: Black=ReG_U, Red=ImG_U"  
  write(*,*)"P_vs_freq.dat: co-power (black) and power (red) vs freq"
  write(*,*)"----------------------------------------------------------"

! Close files
  close(77)
  close(78)
  close(79)
  
end program co_vs_phi


