include 'a_modules.f90'
include 'a_subroutines.f90'

program co_vs_phi
! Caculate rms and phase integrated over a given frequency range as a function
! of modulation angle. This code will stack over multiple observations.
!
!
! gfortran -L$HEADAS/lib -lcfitsio co_vs_phi.f90
!  
  use shared_arrays
  implicit none
  integer pilo,pihi,nseg
  integer evtunit1,evtunit2,evtunit3
  integer i,j,k
  double precision dt
  double precision tstart,tstop
  character (len=500) evtfile1,evtfile2,evtfile3,evtlist
  integer nn,phibins,ilo,ihi,config,status
  real Elo,Ehi,flo,fhi,f,Psn,phideg,dphideg,lagcyc,dlagcyc
  double precision tlive,ton,Telapse,Tseg,tbeg,tend,dphi,Texp
  real, allocatable :: P(:),dP(:),Pco(:),dPco(:)
  real df,Pr,dPr,rate,wn,Prco,dPrco,dPs,nd
  real, allocatable :: Pbin(:),dPbin(:),Pcobin(:),dPcobin(:)
  real, allocatable :: ReG(:),dReG(:),ImG(:),dImG(:)
  real, allocatable :: ReGphi(:),dReGphi(:),ImGphi(:),dImGphi(:),dG(:)
  real, allocatable :: rms(:),drms(:),lag(:),dlag(:),dGopp(:)
  real Psub,dPsub
  real, allocatable :: ReGopp(:),dReGopp(:),ImGopp(:),dImGopp(:)
  real, allocatable :: ReGav(:),ImGav(:),dGav(:),phi(:),Ps(:),Ps_tot(:)
  real, allocatable :: ReGtot(:),ImGtot(:),r_sub(:),r_subtot(:),dr_subtot(:)
  real, allocatable :: ReG_mod(:),ImG_mod(:),hmod(:),r_tot(:),dr_tot(:)
  integer Nobs,obs,listunit,mmtot
  real Pr_tot,Pco_tot,r_ref,r_reftot,Ttot
  integer mmax,ma
  parameter (mmax=20)
  integer ia(mmax)
  real chisq,da(mmax),dyda(mmax),a_re(mmax),a_im(mmax),chisq_full,lag_mod
  real a0(mmax),qav,uav,a_null(mmax),chisq_null,rms_mod,r_all,r_alltot
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
  flo      = 1.0        !low end of frequency range averaged over
  fhi      = 4.0        !high end of frequency range averaged over
  Nobs     = 1                 !Number of of observations to stack over
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
  ilo  = ceiling( flo / df )
  ihi  = floor(   fhi / df )
  allocate( phi(phibins) )
  dphi = pi/dble(phibins)
  dphideg = dphi * 180.0 / pi
  do j = 1,phibins
     phi(j) = -0.5*pi + (dble(j)-0.5) * dphi
  end do

! Open output files
  open(86,file='ReG_ImG_vs_phi.dat')
  open(85,file='rms_lag_vs_phi.dat')
  open(180,file='xReG_ImG_vs_phi.dat')
  open(181,file='xrms_lag_vs_phi.dat')
  open(80,file='rate_vs_phi.dat')
  
! Allocate arrays with constant size
  allocate(  P(nseg/2) )
  allocate( dP(nseg/2) )
  allocate(  ReG(nseg/2) )
  allocate( dReG(nseg/2) )
  allocate(  ImG(nseg/2) )
  allocate( dImG(nseg/2) )
  allocate(  ReGphi(phibins) )
  allocate( dReGphi(phibins) )
  allocate(  ImGphi(phibins) )
  allocate( dImGphi(phibins) )
  allocate( dG(phibins) )
  allocate(  lag(phibins) )
  allocate( dlag(phibins) )
  allocate(  rms(phibins) )
  allocate( drms(phibins) )
  allocate(  Pco(nseg/2) )
  allocate( dPco(nseg/2) )
  allocate( Ps(phibins) )
  allocate( Ps_tot(phibins) )
  allocate( ReGtot(phibins) )
  allocate( ImGtot(phibins) )
  allocate( r_sub(0:phibins) )
  allocate( r_subtot(0:phibins) )
  allocate( dr_subtot(0:phibins) )
  allocate( r_tot(phibins) )
  allocate( dr_tot(phibins) )

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
     
     !Make power spectrum of reference band light curve and average over the frequency range
     call mypow(nseg,lcref,nn,real(dt),mm,segend,2,P,dP,rate,wn)
     call freqav(nseg,P,dP,ilo,ihi,Pr,dPr)
     write(*,*)"Pr=",Pr,"±",dPr

     !Calculate co-spectrum and average over the frequency range
     call myco(nseg,lcsub,lcref,nn,real(dt),mm,segend,2,Pco,dPco)
     call freqav(nseg,Pco,dPco,ilo,ihi,Prco,dPrco)
     write(*,*)"Prco=",Prco,"±",dPrco
     write(*,*)"rms = ",sqrt( Prco * (fhi-flo) / r_all**2 )
     
     !Calculate a cross spectrum for each phi bin
     write(*,*)"Calculating cross spectra..."
     do j = 1,phibins
        !Calculate cross spectra
        call mycross(nseg,lcphi(:,j),lcref,nn,real(dt),mm,segend,2,ReG,dReG,ImG,dImG)
        !write(*,*)"phi range (rad) = ",phi(j)-0.5*dphi,phi(j)+0.5*dphi
        !Calculate mean count rate for each phi bin
        call mean_rate(nseg,lcphi(:,j),nn,mm,segend,r_sub(j))
        !Average over frequency range
        call freqav(nseg,ReG,dReG,ilo,ihi,ReGphi(j),dReGphi(j))
        call freqav(nseg,ImG,dImG,ilo,ihi,ImGphi(j),dImGphi(j))
        !Calculate power spectrum and average
        call mypow(nseg,lcphi(:,j),nn,real(dt),mm,segend,2,P,dP,rate,wn)
        call freqav(nseg,P,dP,ilo,ihi,Ps(j),dPs)
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

! Scaling count rates from 2 DUs to 3 DUs
  write(80,*)"read serr 1 2"
  write(80,*)"skip on"
  write(*,*)"r_alltot / r_subtot(0) = ",r_alltot / r_subtot(0)
  do j = 1,phibins
     phideg    = phi(j) * 180.0 / pi
     r_tot(j)  = r_subtot(j)
     dr_tot(j) = sqrt( r_subtot(j) / Ttot * r_alltot / r_subtot(0) )
     write(80,*)phideg,0.5*dphideg,r_tot(j),dr_tot(j)
  end do

  
! Calculate errors and convert to rms and lag
  nd   = real( mmtot*(ihi-ilo+1) )
  !fac  = ( r_subtot(0) + r_reftot )**2 / ( r_subtot(0) * r_reftot )
  do j = 1,phibins
     !Calculate error from Ingram (2019) formula
     call ingram2019_err_b0(Pr_tot,Ps_tot(j),Pco_tot,ReGtot(j),ImGtot(j),nd,dG(j))
     !Convert to rms and phase lag
     call rms_and_lag(ReGtot(j),ImGtot(j),dG(j),Pco_tot,flo,fhi,rms(j),drms(j),lag(j),dlag(j))
     !Convert rms to fractional
     rms(j)  =  rms(j) / r_tot(j)
     drms(j) = drms(j) / r_tot(j)
  end do

! Plot the flo-fhi cross spectrum vs phi
  write(86,*)"read serr 1 2 3"
  write(86,*)"skip on"
  write(85,*)"read serr 1 2 3"
  write(85,*)"skip on"
  do j = 1,phibins
     phideg =  phi(j) * 180.0 / pi
     write(86,*)phideg,0.5*dphideg,ReGtot(j),dG(j),ImGtot(j),dG(j)
     write(85,*)phideg,0.5*dphideg,rms(j),drms(j),lag(j)/(2.0*pi),dlag(j)/(2.0*pi)
  end do
  
! Write out in xspec format
! (need to plot from 180 to 540 degrees to get around xspec problem with negative x-axis)
  do j = 1,phibins
     phideg = phi(j) * 180.0 / pi
     if( phideg .gt. 0.0 )then
        lagcyc  =  lag(j) / (2.0*pi)
        dlagcyc = dlag(j) / (2.0*pi)
        write(180,*)phideg-0.5*dphideg,phideg+0.5*dphideg,rms(j)*dphideg,drms(j)*dphideg
        write(181,*)phideg-0.5*dphideg,phideg+0.5*dphideg,lagcyc*dphideg,dlagcyc*dphideg
     end if
  end do


  write(*,*)"-------------------------------------"
  write(*,*)"Total number of segments = ",mmtot
  write(*,*)"Number of frequencies = ",ihi-ilo+1
  write(*,*)"Lowest frequency = ",ilo,ilo*df
  write(*,*)"Highest frequency = ",ihi,ihi*df
  write(*,*)"Total realisations = ",nd
  write(*,*)"-------------------------------------"
  
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
  write(80,*)"no no"
  do j = 1,phibins
     phideg = phi(j) * 180.0 / pi
     call IQUsinfunc(phi(j),a0,hmod(j),dyda,ma)
     write(80,*)phideg,0.5*dphideg,hmod(j),0.0
  end do
  qav = a0(2)/a0(1)
  uav = a0(3)/a0(1)

! Fit full model to cross spectra

  !Fit for ReGtot
  ia      = 1
  a_re(1) = 100.0   !P_I dphi/pi
  a_re(2) = 4.0     !ReG_Q dphi/pi
  a_re(3) = 1.0     !ReG_U dphi/pi
  call dofit(phi,ReGtot,dG,phibins,a_re,ia,ma,IQUsinfunc,chisq,da)
  chisq_full = chisq
  
  !Fit for ImGtot
  ia      = 1
  ia(1)   = 0
  a_im(1) = 0.0     !ImP_I dphi/pi (zero)
  a_im(2) = 4.0     !ImG_Q dphi/pi
  a_im(3) = 1.0     !ImG_U dphi/pi
  call dofit(phi,ImGtot,dG,phibins,a_im,ia,ma,IQUsinfunc,chisq,da)
  chisq_full = chisq_full + chisq

  !Chi squared for full model
  write(*,*)"Full model:"
  write(*,*)"chisq/dof=",chisq_full,"/",2*phibins-5
  
  !Write out the model
  write(86,*)"no no"
  write(85,*)"no no"
  do j = 1,phibins
     phideg = phi(j) * 180.0 / pi
     call IQUsinfunc(phi(j),a_re,ReG_mod(j),dyda,ma)
     call IQUsinfunc(phi(j),a_im,ImG_mod(j),dyda,ma)
     write(86,*)phideg,0.5*dphideg,ReG_mod(j),0.0,ImG_mod(j),0.0
     !Convert to rms and phase lag
     rms_mod = sqrt( ReG_mod(j)**2 + ImG_mod(j)**2 )
     rms_mod = rms_mod / hmod(j)
     rms_mod = rms_mod * sqrt( (fhi-flo) / Pco_tot )
     lag_mod = atan2( ImG_mod(j) , ReG_mod(j) ) / (2.0*pi)
     write(85,*)phideg,0.5*dphideg,rms_mod,0.0,lag_mod,0.0
  end do
  write(85,*)"no no"

  write(*,*)"rms^2 = ",a_re(1)*phibins*(fhi-flo)/r_alltot**2
  write(*,*)"rms = ",sqrt(a_re(1)*phibins*(fhi-flo)/r_alltot**2)

  
! Fit null hypothesis to cross spectra
  !Fit for ReGtot
  ia        = 0
  ia(1)     = 1
  a_null(1) = 100.0   !P_I dphi/pi
  a_null(2) = qav     !q
  a_null(3) = uav     !u
  call dofit(phi,ReGtot,dG,phibins,a_null,ia,ma,NULLsinfunc,chisq,da)
  chisq_null = chisq

  !Write out the model
  write(86,*)"no no"
  chisq = 0.0
  do j = 1,phibins
     phideg = phi(j) * 180.0 / pi
     call NULLsinfunc(phi(j),a_null,ReG_mod(j),dyda,ma)
     ImG_mod(j) = 0.0
     write(86,*)phideg,0.5*dphideg,ReG_mod(j),0.0,ImG_mod(j),0.0
     chisq = chisq + ( ImGtot(j) / dG(j) )**2
     !Convert to rms and phase lag
     rms_mod = sqrt( ReG_mod(j)**2 + ImG_mod(j)**2 )
     rms_mod = rms_mod / hmod(j)
     rms_mod = rms_mod * sqrt( (fhi-flo) / Pco_tot )
     lag_mod = atan2( ImG_mod(j) , ReG_mod(j) ) / (2.0*pi)
     write(85,*)phideg,0.5*dphideg,rms_mod,0.0,lag_mod,0.0
  end do
  chisq_null = chisq_null + chisq
  write(*,*)"Null hypothesis model:"
  write(*,*)"chisq/dof=",chisq_null,"/",2*phibins-1
  
  
! Write out commands

  write(80,*)"skip on"
  write(80,*)"li s on 2"
  
  write(86,*)"wi 1"
  write(86,*)"v .15 .55 .9 .95"
  write(86,*)"yplot 1,3,5"
  write(86,*)"wi 2"
  write(86,*)"v .15 .15 .9 .55"
  write(86,*)"yplot 2,4,6"
  write(86,*)"la x Modulation Angle (deg)"
  write(86,*)"la y ImG"
  write(86,*)"r y"
  write(86,*)"wi 1"
  write(86,*)"la y ReG"
  write(86,*)"la nx off"
  write(86,*)"r y"
  write(86,*)"ma 17 on 1,2"
  write(86,*)"ma size 2 on 1,2"
  write(86,*)"li s on 3,4,5,6"
  write(86,*)"lw 5"
  write(86,*)"co 1 on 1,2"
  write(86,*)"co 2 on 3,4"
  write(86,*)"co 4 on 5,6"
  
  write(85,*)"wi 1"
  write(85,*)"v .15 .55 .9 .95"
  write(85,*)"yplot 1,3"
  write(85,*)"wi 2"
  write(85,*)"v .15 .15 .9 .55"
  write(85,*)"yplot 2,4,6"
  write(85,*)"la x Modulation Angle (deg)"
  write(85,*)"la y Phase Lag (cycles)"
  write(85,*)"r y"
  write(85,*)"wi 1"
  write(85,*)"la y rms"
  write(85,*)"la nx off"
  write(85,*)"r y"
  write(85,*)"li s on 3,4"
  write(85,*)"co 1 on 1,2"
  write(85,*)"co 2 on 3,4"
  write(85,*)"co 15 on 5,6"

  

  write(*,*)"----------------------------------------------------------"
  write(*,*)"Outputs:"
  write(*,*)"ReG_ImG_vs_phi.dat: ReG and ImG vs modulation angle"
  write(*,*)"rms_lag_vs_phi.dat: rms and lag vs modulation angle"
  write(*,*)"xReG_ImG_vs_phi.dat: XSPEC format rms vs modulation angle"
  write(*,*)"xrms_lag_vs_phi.dat: XSPEC format lag vs modulation angle"
  write(*,*)"rate_vs_phi.dat: Count rate vs modulation angle"
  write(*,*)"----------------------------------------------------------"
  
! Close files
  close(86)
  close(85)
  close(180)
  close(181)
  close(80)
  
end program co_vs_phi


