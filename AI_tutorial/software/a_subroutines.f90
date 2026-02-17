
!-----------------------------------------------------------------------
subroutine mean_rate(nseg,lc,nn,mm,segend,mu)
!
! Calculates mean count rate of a light curve, only accounting for
! segments that are used for the timing analysis
!
! INPUTS:
! nseg        Segment length
! lc(nn)      Light curve (expected units count rate)
! nn          Total length of light curve
! mm          Number of segments  
! segend(mm)  List of time bins at end of segments
! OUTPUTS:
! mu          Mean count rate
!
  implicit none
  integer nseg,nn,mm,segend(mm)
  real lc(nn),lcseg(nseg),mu
  integer i,j,m

  mu = 0.0
  do m = 1,mm
     do i = 1,nseg
        j        = i + segend(m) - nseg
        lcseg(i) = lc(j)
        mu       = mu + lcseg(i)
     end do
  end do
  mu = mu / real(nseg*mm)

  return
end subroutine mean_rate     
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine calc_rms_lag(nphase,phibins,Rpar,PDpar,PApar,mu,rms,lag)
! Given parameters for I, PD and PA modulations, this routine calculates
! rms and lag vs modulation angle.
! rms is fractional
! lag is in radians
!
! Inputs:
! nphase      Number of QPO phases
! phibins     Number of modulation angle bins
! Rpar(mmax)  Parameters of count rate modulation model
! PDpar(mmax) Parameters of PD modulation model
! PApar(mmax) Parameters of PA modulation model
! mu          Modulation factor
!
! Outputs:
! rms(nphase/2,phibins)  Fractional rms of each harmonic for each modulation angle bin
! lag(nphase/2,phibins)  Phase lag of each harmonic for each modulation angle bin
!
  implicit none
  !Inputs
  integer nphase,phibins,mmax
  parameter (mmax=20)
  double precision Rpar(mmax),PDpar(mmax),PApar(mmax),mu
  !Outputs
  real rms(nphase/2,phibins),lag(nphase/2,phibins)
  !Internal
  real R_t(nphase),PD_t(nphase),PA_t(nphase),ReR(0:nphase/2),ImR(0:nphase/2)
  real mt(nphase,phibins)

! Calculate rate, PD and PA vs phase from parameters
  call mkmods(nphase,Rpar,PDpar,PApar,R_t,PD_t,PA_t)

! Fourier transform rate (which will be the reference band)
  call doFFT(1.0/real(nphase),nphase,R_t,ReR,ImR)

! Generate count rate vs phase for each phi bin
  call mklcs(nphase,phibins,real(mu),R_t,PD_t,PA_t,mt)

! Calculate rms and lag wrt rate for each of the phibins light curves
  call mklagamp(nphase,phibins,ReR,ImR,mt,rms,lag)

  return
end subroutine calc_rms_lag
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine lagamp_frommods(nphase,phibins,R_t,PD_t,PA_t,mu,rms,lag)
! Given I, PD and PA modulations (in nphase phases), this routine
! calculates
! rms and lag vs modulation angle.
! rms is fractional
! lag is in radians
!
! Inputs:
! nphase       Number of QPO phases
! phibins      Number of modulation angle bins
! Rt(nphase)   Count rate as a function of phase
! PD_t(nphase) PD as a function of phase
! PA_t(nphase) PA as a function of phase
! mu          Modulation factor
!
! Outputs:
! rms(nphase/2,phibins)  Fractional rms of each harmonic for each modulation angle bin
! lag(nphase/2,phibins)  Phase lag of each harmonic for each modulation angle bin
!
  implicit none
  !Inputs
  integer nphase,phibins
  real R_t(nphase),PD_t(nphase),PA_t(nphase),mu
  !Outputs
  real rms(nphase/2,phibins),lag(nphase/2,phibins)
  !Internal
  real ReR(0:nphase/2),ImR(0:nphase/2),mt(nphase,phibins)

! Fourier transform rate (which will be the reference band)
  call doFFT(1.0/real(nphase),nphase,R_t,ReR,ImR)

! Generate count rate vs phase for each phi bin
  call mklcs(nphase,phibins,mu,R_t,PD_t,PA_t,mt)

! Calculate rms and lag wrt rate for each of the phibins light curves
  call mklagamp(nphase,phibins,ReR,ImR,mt,rms,lag)

  return
end subroutine lagamp_frommods  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine mklagamp(nphase,phibins,ReR,ImR,mt,rms,lag)
!
! Inputs:
! nphase                 Number of phase bins
! phibins                Number of phi bins, where phi is modulation angle.
! ReR(0:nphase/2)        Real part of FT of count rate
! ImR(0:nphase/2)        Imaginary part of FT of count rate
! mt(nphase,phibins)     mt(k,j) is the count rate at the kth phase in the jth phi bin
!
! Outputs:
! rms(nphase/2,phibins) rms of each harmonic for each phi (FRACTION)
! lag(nphase/2,phibins) phase lag of each harmonic for each phi (RADIANS)
! 
  implicit none
  !Inputs
  integer nphase,phibins
  real ReR(0:nphase/2),ImR(0:nphase/2),mt(nphase,phibins)
  !Outputs
  integer j,k
  real ReM(0:nphase/2,phibins),ImM(0:nphase/2,phibins),S(phibins)
  real ReG,ImG,rms(nphase/2,phibins),lag(nphase/2,phibins)

  do j = 1,phibins
     !Calculate FFT
     call doFFT(1.0/real(nphase),nphase,mt(:,j),ReM(:,j),ImM(:,j))
     S(j) = ReM(0,j)  !mean count rate in each phi bin
     do k = 1,nphase/2
        rms(k,j) = sqrt( ReM(k,j)**2 + ImM(k,j)**2 ) / S(j)
        ReG      = ReM(k,j) * ReR(k) + ImM(k,j)*ImR(k)
        ImG      = ImM(k,j) * ReR(k) - ReM(k,j)*ImR(k)
        lag(k,j) = atan2( ImG , ReG )
     end do
  end do
  
  return
end subroutine mklagamp     
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mklagampct(nphase,phibins,ReR,ImR,mt,rms,lag,S)
!
! Inputs:
! nphase                 Number of phase bins
! phibins                Number of phi bins, where phi is modulation angle.
! ReR(0:nphase/2)        Real part of FT of count rate
! ImR(0:nphase/2)        Imaginary part of FT of count rate
! mt(nphase,phibins)     mt(k,j) is the count rate at the kth phase in the jth phi bin
!
! Outputs:
! rms(nphase/2,phibins)  rms of each harmonic for each phi (FRACTION)
! lag(nphase/2,phibins)  phase lag of each harmonic for each phi (RADIANS)
! S(phibins)             count rate for each phi
! 
  implicit none
  !Inputs
  integer nphase,phibins
  real ReR(0:nphase/2),ImR(0:nphase/2),mt(nphase,phibins)
  !Outputs
  integer j,k
  real ReM(0:nphase/2,phibins),ImM(0:nphase/2,phibins),S(phibins)
  real ReG,ImG,rms(nphase/2,phibins),lag(nphase/2,phibins)

  do j = 1,phibins
     !Calculate FFT
     call doFFT(1.0/real(nphase),nphase,mt(:,j),ReM(:,j),ImM(:,j))
     S(j) = ReM(0,j)  !mean count rate in each phi bin
     do k = 1,nphase/2
        rms(k,j) = sqrt( ReM(k,j)**2 + ImM(k,j)**2 ) / S(j)
        ReG      = ReM(k,j) * ReR(k) + ImM(k,j)*ImR(k)
        ImG      = ImM(k,j) * ReR(k) - ReM(k,j)*ImR(k)
        lag(k,j) = atan2( ImG , ReG )
     end do
  end do
  
  return
end subroutine mklagampct
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mklcs(nphase,phibins,mu,R_t,PD_t,PA_t,mt)
! Routine to make phibins light curves with nphase steps.
! Each light curve is for a different phi bin, where phi is modulation angle.
  implicit none
  integer nphase,phibins,j,k
  real mu,R_t(nphase),PD_t(nphase),PA_t(nphase)
  real mt(nphase,phibins),phi,dphi,pi,f,fphi
  pi    = acos(-1.0)
  dphi  = pi / real(phibins)
  do j = 1,phibins
     phi = -0.5*pi + ( real(j) - 0.5 ) * dphi
     do k = 1,nphase
        !Calculate modulation function
        f = fphi(phi,mu,PD_t(k),PA_t(k))
        !Calculate counts in the phase and phi bin
        mt(k,j) = R_t(k) * f * dphi
     end do
  end do
  return
end subroutine mklcs
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
function fphi(phi,mu,PD,PA)
! phi & PA in radians
  implicit none
  real fphi,phi,mu,PD,PA
  real pi
  pi   = acos(-1.0)
  fphi = 1. + mu*PD * cos( 2.*(PA-phi) )
  fphi = fphi / pi
  return
end function fphi
!-----------------------------------------------------------------------
      

!-----------------------------------------------------------------------
subroutine mkmods(nphase,Rpar,PDpar,PApar,R_t,PD_t,PA_t)
! Calculates count rate, PD and PA vs QPO (or pulse) phase for nphase phases.
! Input: nphase, Rpar(mmax), PDpar(mmax), PApar(mmax)
! Output: R_t(nphase), PD_t(nphase), PA_t(nphase)
  implicit none
  !Input
  integer nphase,mmax
  parameter (mmax=20)
  double precision Rpar(mmax),PDpar(mmax),PApar(mmax)
  !Output
  real R_t(nphase),PD_t(nphase),PA_t(nphase)
  !Internal
  integer k
  double precision phase,R,PD,PA,dyda(mmax)
  
  do k = 1,nphase
     phase = real(k) / real(nphase)
     call dwffunc(phase,Rpar,R,dyda,5)
     call dwffunc(phase,PDpar,PD,dyda,5)
     PD = max( PD , 0.d0 )
     call dwffunc(phase,PApar,PA,dyda,5)
     R_t(k)  = R
     PD_t(k) = PD
     PA_t(k) = PA
  end do

  return
end subroutine mkmods
!-----------------------------------------------------------------------




!------------------------------------------------------------------------
subroutine doFFT(dt,n,at,ReA,ImA)
  implicit none
  integer n,j
  real at(n),ReA(0:n/2),ImA(0:n/2)
  real data(2*n),dt,mu
  mu = 0.0
  do j = 1,n
     data(2*j-1) = at(j)
     data(2*j)   = 0.0
     mu = mu + at(j)
  end do
  mu = mu / real(n)
  call ourfour1(data,n,1)
  do j = 1, n/2
     ReA(j) = data(2*j+1) * sqrt( 2. * dt / real(n) )
     ImA(j) = data(2*j+2) * sqrt( 2. * dt / real(n) )
  end do
  ReA(0) = mu
  ImA(0) = 0.0
  return
end subroutine doFFT
!------------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine mysteppar(ma,j,kmax,par,dpar,ipar,chisqmin,pmin,pmax,ndata,x,y,dy,funcs,unit)
! Makes a steppar plot in unit number "unit". 
! ma        = number f model parameters
! j         = parameter number to do the steppar over
! kmax      = number of steps
! par(ma)   = best fitting parameters
! dpar(ma)  = uncertainty estimate on best fitting parameters
! ipar(ma)  = free or fixed array
! chisqmin  = chi squared of best fitting model
! pmin      = minimum value of parameter in steppar
! pmax      = maximum value of parameter in steppar
! ndata     = number of data points
! x(ndata)  = x-values of data points
! y(ndata)  = y-values of data points
! dy(ndata) = y-errors of data points
! funcs     = name of function to fit
! unit      = unit number to plot results in
  implicit none
  !Inputs
  integer ma,j,kmax,unit,ndata
  integer ipar(ma)
  real par(ma),dpar(ma),chisqmin,pmin,pmax,x(ndata),y(ndata),dy(ndata)
  external funcs
  !Internal
  integer ipar_s(ma),k
  real par_s(ma),dpar_s(ma),chisq
  
! Transfer parameter array to internal array
  par_s     = par
  dpar_s    = dpar
  ipar_s    = ipar

! Freeze parameter to be stepped over
  ipar_s(j) = 0

! Step over the parameter
  do k = 1,kmax
     par_s(j) = pmin + (pmax-pmin) * real(k-1)/real(kmax-1)
     call dofit(x,y,dy,ndata,par_s,ipar_s,ma,funcs,chisq,dpar_s)
     write(unit,*)par_s(j),chisq,chisqmin+1.0
  end do

  return
end subroutine mysteppar  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine multilorefunc(x,a,ymod,dyda,ma)
!
! Sum of N lorentzians. N = ma / 3
!
! Inputs:
! x          = frequency
! a(1:ma)    = array of parameters
! ma         = number of parameters
! Outputs:
! ymod       = model evaluation at x
! dyda(1:ma) = derivative wrt a
!
! Parameters:
! a(1)    = Lore 1 nu0 (centroid)
! a(2)    = Lore 1 Delta (HWHM)
! a(3)    = Lore 1 norm (squared fractional rms)
! a(4)    = Lore 2 nu0 (centroid)
! .
! a(3N)   = Lore N norm (squared fractional rms)
!  
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax),pi
  real renorm,fac,fac2,nu0,HWHM,norm
  integer N,i
  pi = acos(-1.0)
  N  = ma / 3
  ymod = 0.0
  do i = 1,N
     !Set Lorentzian parameters
     nu0  = a( (i-1)*3 + 1 )
     HWHM = a( (i-1)*3 + 2 )
     norm = a( (i-1)*3 + 3 )
     !Derived quantities     
     renorm = 1.0 / ( pi/2. + atan( nu0/HWHM ) )
     fac    = 1.0 / ( ( x-nu0 )**2. + HWHM**2. )
     fac2   = 1.0 / ( (nu0/HWHM)**2 + 1.0 )
     !Model evaluation
     ymod = ymod + 2.0 * norm * HWHM * renorm * fac
     !Derivatives
     !dL/dnu0
     dyda((i-1)*3 + 1) = 2.0*HWHM * ( x-nu0 ) * fac**2 * renorm
     dyda((i-1)*3 + 1) = dyda((i-1)*3 + 1) - renorm**2 * fac * fac2
     dyda((i-1)*3 + 1) = dyda((i-1)*3 + 1) * 2.0 * norm
     !dL/dDelta
     dyda((i-1)*3 + 2) = nu0 * fac2 * fac * renorm**2 / HWHM
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) - 2.0*HWHM**2 * fac**2 * renorm
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) + fac * renorm
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) * 2.0 * norm
     !dL/dN
     dyda((i-1)*3 + 3) = 2.0 * HWHM * renorm * fac
  end do
  return
end subroutine multilorefunc
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine multilorerfunc(x,a,ymod,dyda,ma)
!
! Sum of N lorentzians but norm is fractional rms,
! not fractional squared rms. N = ma / 3
!
! Inputs:
! x          = frequency
! a(1:ma)    = array of parameters
! ma         = number of parameters
! Outputs:
! ymod       = model evaluation at x
! dyda(1:ma) = derivative wrt a
!
! Parameters:
! a(1)    = Lore 1 nu0 (centroid)
! a(2)    = Lore 1 Delta (HWHM)
! a(3)    = Lore 1 fractional rms
! a(4)    = Lore 2 nu0 (centroid)
! .
! a(3N)   = Lore N fractional rms
!  
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax),pi
  real renorm,fac,fac2,nu0,HWHM,norm,rms
  integer N,i
  pi = acos(-1.0)
  N  = ma / 3
  ymod = 0.0
  do i = 1,N
     !Set Lorentzian parameters
     nu0  = a( (i-1)*3 + 1 )
     HWHM = a( (i-1)*3 + 2 )
     rms  = a( (i-1)*3 + 3 )
     norm = rms**2
     !Derived quantities     
     renorm = 1.0 / ( pi/2. + atan( nu0/HWHM ) )
     fac    = 1.0 / ( ( x-nu0 )**2. + HWHM**2. )
     fac2   = 1.0 / ( (nu0/HWHM)**2 + 1.0 )
     !Model evaluation
     ymod = ymod + 2.0 * norm * HWHM * renorm * fac
     !Derivatives
     !dL/dnu0
     dyda((i-1)*3 + 1) = 2.0*HWHM * ( x-nu0 ) * fac**2 * renorm
     dyda((i-1)*3 + 1) = dyda((i-1)*3 + 1) - renorm**2 * fac * fac2
     dyda((i-1)*3 + 1) = dyda((i-1)*3 + 1) * 2.0 * norm
     !dL/dDelta
     dyda((i-1)*3 + 2) = nu0 * fac2 * fac * renorm**2 / HWHM
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) - 2.0*HWHM**2 * fac**2 * renorm
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) + fac * renorm
     dyda((i-1)*3 + 2) = dyda((i-1)*3 + 2) * 2.0 * norm
     !dL/drms
     dyda((i-1)*3 + 3) = 4.0 * rms * HWHM * renorm * fac
  end do
  return
end subroutine multilorerfunc
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine ingram2019_err_b0(Pr,Ps,Prco,ReG,ImG,nd,dG)
! Calculates Ingram (2019) error on the energy-dependent cross spectrum.
! This is equation 18 of Ingram (2019) with b=0.
! Use this error formula if:
! 1) You are calculating multiple cross spectra, all with the same
!    reference band (otherwise use Bendat & Piersol errors);
! 2) The number of realisations averaged over is nd > 500 (otherwise
!    you need to subtract off a bias term from the cross spectrum).
!
! Inputs:
! Pr, Ps     Reference and subject band power spectra (include Poisson noise)
! Prco       Reference band Poisson noise-subtracted power spectrum
! ReG, ImG   Complex cross-spectrum
! nd         Number of realisations averaged over
! Outputs:
! dG         Error on cross-spectrum: dG = dReG = dImG = d|G|.
  implicit none
! Inputs
  real Pr,Ps,Prco,ReG,ImG,nd
! Outputs
  real dG
! Calculation
  dG = Pr*Ps - (Pr/Prco) * ( ReG**2 + ImG**2 )     
  dG = dG / ( 2.0 * nd )
  dG = sqrt( dG )
  return
end subroutine ingram2019_err_b0
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine rms_and_lag(ReG,ImG,dG,Prco,flo,fhi,rms,drms,lag,dlag)
! Calculates rms and phase in radians and their associated errors.
! For use when real and imaginary parts have the same error.
! ReG and ImG are in units of rms^2/Hz.
! Inputs:
! ReG, ImG:    Real and imaginary parts of the cross spectrum.
! dG           dG = dReG = dImG
! Prco         Poisson noise-subtracted reference band power spectrum
! flo, fhi     Lower and upper bounds of frequency range
! Outputs:
! rms, drms    In the same units as G
! lag, dlag    In radians  
  implicit none
! Inputs
  real ReG,ImG,dG,Prco,flo,fhi
! Outputs
  real rms,drms,lag,dlag
! Calculation
  rms  = sqrt( (ReG**2+ImG**2) ) * sqrt( (fhi-flo) / Prco )
  drms = dG * sqrt( (fhi-flo) / Prco )
  lag  = atan2( ImG , ReG )
  dlag = dG / sqrt( ReG**2 + ImG**2 )
  return
end subroutine rms_and_lag
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine rms_and_lag2(ReG,dReG,ImG,dImG,QPOrms,rms,drms,lag,dlag)
! Calculates rms and phase in radians and their associated errors.
! For use when real and imaginary parts of the cross spectrum have been
! integrated (i.e. via a Lorentzian fit). Therefore ReG and ImG are in
! units of rms^2, NOT rms^2/Hz.
! Inputs:
! ReG, ImG:    Real and imaginary parts of the integrated cross spectrum.
! dReG, dImG   Error on ReG and ImG
! QPOrms       rms of the QPO in the reference band
! Outputs:
! rms, drms    Fractional rms
! lag, dlag    In radians  
  implicit none
! Inputs
  real ReG,dReG,ImG,dImG,QPOrms
! Outputs
  real rms,drms,lag,dlag
! Calculation
  rms  = sqrt( (ReG**2+ImG**2) ) / QPOrms
  drms = ( (ReG*dReG)**2 + (ImG*dImG)**2 ) / ( ReG**2 + ImG**2 )
  drms = sqrt(drms) / QPOrms
  lag  = atan2( ImG , ReG )
  dlag = sqrt( (ImG*dReG)**2 + (ReG*dImG)**2 ) / ( ReG**2 + ImG**2 )
  return
end subroutine rms_and_lag2
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine freqav(nseg,G,dG,ilo,ihi,Gphi,dGphi)
  implicit none
  integer nseg,ilo,ihi,i
  real G(nseg/2),dG(nseg/2),Gphi,dGphi
  Gphi  = 0.0
  dGphi = 0.0
  do i = ilo,ihi
     Gphi  = Gphi  +  G(i)
     dGphi = dGphi + dG(i)**2
  end do
  Gphi  = Gphi / real( ihi-ilo+1 )
  dGphi = sqrt( dGphi ) / real( ihi-ilo+1 )
  return
end subroutine freqav
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine geobin_scheme(nseg,df,c)
  use shared_arrays
  implicit none
  !Inputs
  integer nseg
  real df,c
  !Outputs
  !in shared_arrays :: far(0:nf), nf, np(nf)
  !Internal
  integer i,j,remain,bins
  
! First calculate the number of new frequency bins
  i = 0
  j = 0
  do while( i .lt. nseg/2 )
     j      = j + 1
     remain = nseg/2 - i
     bins   = min( floor( c**j ) , remain )
     i      = i + bins
  end do
  nf      = j

! Now allocate arrays
  if( allocated( np ) ) deallocate( np )
  allocate( np(nf) )
  if( allocated( far ) ) deallocate( far )
  allocate( far(0:nf) )
  if( allocated( iar ) ) deallocate( iar )
  allocate( iar(0:nf) )

! Now fill arrays
  i = 0
  j = 0
  far(0) = 0.5 * df
  iar(0) = 0
  do while( i .lt. nseg/2 )
     j      = j + 1
     remain = nseg/2 - i
     np(j)  = min( floor( c**j ) , remain )
     i      = i + np(j)
     far(j) = ( i + 0.5 ) * df
     iar(j) = i
   end do
   far(nf) = ( nseg/2 + 0.5 ) * df
   iar(nf) = nseg/2

   return
 end subroutine geobin_scheme  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine geobin_bin(nseg,P,dP,Pbin,dPbin)
  use shared_arrays
  implicit none
  !Inputs
  integer nseg
  real P(nseg),dP(nseg)
  !Outputs
  real Pbin(nf),dPbin(nf)
  !Internal
  integer j,i
  do j = 1,nf
     Pbin(j)  = 0.0
     dPbin(j) = 0.0
     do i = iar(j-1)+1,iar(j)
        Pbin(j)  = Pbin(j) + P(i)
        dPbin(j) = dPbin(j) + dP(i)**2
     end do
     Pbin(j)  = Pbin(j) / real( np(j) )
     dPbin(j) = sqrt( dPbin(j) ) / real( np(j) )
  end do
  return
end subroutine geobin_bin
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine calc_spur_corr(evtunit1,evtunit2,pilo,pihi,nn,phibins,phi,nseg,corr)
! Calculates correction for spurious polarisation needed for each
! modulation angle bin. Also plots modulation function before and after correction.
  use shared_arrays
  implicit none
  integer evtunit1,evtunit2,pilo,pihi,nn,phibins,nseg
  real phi(phibins),corr(phibins)
  double precision Qsp,Usp
  double precision Nsp1,Qsp1,Usp1,Nsp2,qsp2,usp2,Nsp
  integer m,i,j,k
  real rate
  
! Calculate spurious polarisation Stokes parameters for DU1 and DU2
  call sp_stokes(evtunit1,nevts1,pilo,pihi,Nsp1,Qsp1,Usp1)
  call sp_stokes(evtunit2,nevts2,pilo,pihi,Nsp2,Qsp2,Usp2)

! Sum up
  Nsp = Nsp1 + Nsp2
  Qsp = Qsp1 + Qsp2
  Usp = Usp1 + Usp2
  
! Normalise
  qsp = Qsp / Nsp
  usp = Usp / Nsp
  write(*,*)"++qsp=",qsp
  write(*,*)"++usp=",usp
  
! Calculate count rate in the DU1+DU2 light curve
  rate = 0.0
  do m = 1,mm
     do i = 1,nseg
        j    = i + segend(m) - nseg
        rate = rate + lcsub(j)
     end do
  end do
  rate = rate / real(nseg*mm)
  write(*,*)"++DU1+DU2 count rate=",rate

! Calculate spurious modulation
  do j = 1,phibins
     corr(j)  = qsp * cos( 2.0*phi(j) ) + usp * sin( 2.0*phi(j) )
     corr(j)  = corr(j) * rate
     corr(j)  = corr(j) / real(phibins)
  end do
     
  return
end subroutine calc_spur_corr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine calc_spur_corr_opp(evtunit3,pilo,pihi,nn,phibins,phi,nseg,corr_opp)
! Calculates correction for spurious polarisation needed for each
! modulation angle bin. Also plots modulation function before and after correction.
  use shared_arrays
  implicit none
  integer evtunit3,pilo,pihi,nn,phibins,nseg
  real phi(phibins),corr_opp(phibins)
  double precision Qsp,Usp,Nsp
  integer m,i,j,k
  real rate
  
! Calculate spurious polarisation Stokes parameters for DU3
  call sp_stokes(evtunit3,nevts3,pilo,pihi,Nsp,Qsp,Usp)
  
! Normalise
  qsp = Qsp / Nsp
  usp = Usp / Nsp
  write(*,*)"++qsp=",qsp
  write(*,*)"++usp=",usp
  
! Calculate count rate in the DU3 light curve
  rate = 0.0
  do m = 1,mm
     do i = 1,nseg
        j    = i + segend(m) - nseg
        rate = rate + lcref(j)
     end do
  end do
  rate = rate / real(nseg*mm)
  write(*,*)"++DU3 count rate=",rate

! Calculate spurious modulation
  do j = 1,phibins
     corr_opp(j)  = qsp * cos( 2.0*phi(j) ) + usp * sin( 2.0*phi(j) )
     corr_opp(j)  = corr_opp(j) * rate
     corr_opp(j)  = corr_opp(j) / real(phibins)
  end do
     
  return
end subroutine calc_spur_corr_opp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine spur_corr(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,phi,nseg,config)
! Corrects the light curves as a function of modulation angle phi
! for a constant spurious polarisation.
  use shared_arrays
  implicit none
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,nseg,config
  real phi(phibins)
  double precision Qsp,Usp
  double precision Nsp1,Qsp1,Usp1,Nsp2,qsp2,usp2,Nsp
  double precision Q,U
  double precision N1,Q1,U1,N2,q2,u2,N
  integer m,i,j,k
  real rate,corr,rate_phi,model
  
! Extract spurious polarisation Stokes parameters for the two DUs
! used for the phi bin light curves
  if( config .eq. 1 )then
     call sp_stokes(evtunit3,nevts3,pilo,pihi,Nsp1,Qsp1,Usp1)
     call sp_stokes(evtunit2,nevts2,pilo,pihi,Nsp2,Qsp2,Usp2)
     call stokes(evtunit3,nevts3,pilo,pihi,N1,Q1,U1)
     call stokes(evtunit2,nevts2,pilo,pihi,N2,Q2,U2)
  else if( config .eq. 2 )then
    call sp_stokes(evtunit1,nevts1,pilo,pihi,Nsp1,Qsp1,Usp1)
    call sp_stokes(evtunit3,nevts3,pilo,pihi,Nsp2,Qsp2,Usp2)
    call stokes(evtunit1,nevts1,pilo,pihi,N1,Q1,U1)
    call stokes(evtunit3,nevts3,pilo,pihi,N2,Q2,U2)
  else if( config .eq. 3 )then
     call sp_stokes(evtunit1,nevts1,pilo,pihi,Nsp1,Qsp1,Usp1)
     call sp_stokes(evtunit2,nevts2,pilo,pihi,Nsp2,Qsp2,Usp2)
     call stokes(evtunit1,nevts1,pilo,pihi,N1,Q1,U1)
     call stokes(evtunit2,nevts2,pilo,pihi,N2,Q2,U2)
  end if
  
! Sum up
  Nsp = Nsp1 + Nsp2
  Qsp = Qsp1 + Qsp2
  Usp = Usp1 + Usp2
  N   = N1   + N2
  Q   = Q1   + Q2
  U   = U1   + U2
  
! Normalise
  qsp = Qsp / Nsp
  usp = Usp / Nsp
  write(*,*)"qsp=",qsp
  write(*,*)"usp=",usp

  q = Q / N
  u = U / N
  write(*,*)"q (level 2) = ",q
  write(*,*)"u (level 2) = ",u
  
! Calculate count rate in the summed DU light curve
  rate = 0.0
  do m = 1,mm
     do i = 1,nseg
        j    = i + segend(m) - nseg
        rate = rate + lcsub(j)
     end do
  end do
  rate = rate / real(nseg*mm)
  if( config .eq. 1 )then
     write(*,*)"DU2+DU3 count rate=",rate
  else if( config .eq. 2 )then
     write(*,*)"DU1+DU3 count rate=",rate
  else if( config .eq. 3 )then
     write(*,*)"DU1+DU2 count rate=",rate
  end if
     
! Subtract spurious modulation
  do j = 1,phibins

     !Calculate SP correction
     corr  = qsp * cos( 2.0*phi(j) ) + usp * sin( 2.0*phi(j) )
     corr  = corr * rate
     corr  = corr / real(phibins)

     ! !Diagnostic output
     ! rate_phi = 0.0
     ! do m = 1,mm
     !    do i = 1,nseg
     !       k    = i + segend(m) - nseg
     !       rate_phi = rate_phi + lcphi(k,j)
     !    end do
     ! end do
     ! rate_phi = rate_phi / real(nseg*mm)
     ! model    = 1.0 + q*cos(2.0*phi(j)) + u*sin(2.0*phi(j))
     ! model    = model * rate / real(phibins)
     ! write(87,*)phi(j)*180./pi,rate_phi,corr,rate_phi-corr,model

     !Subtract SP correction
     lcphi(:,j) = lcphi(:,j) - corr
     
  end do
     
  return
end subroutine spur_corr
!-----------------------------------------------------------------------








!-----------------------------------------------------------------------
subroutine modhist(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,tstart,dt,phibins,phi,Texp)
! Output: phihist(phibins), dphihist(phibins) via shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins
  double precision tstart,dt,Texp
  real phi(phibins)
  real rate,corr
  double precision Nsp,Nsp1,Nsp2,Nsp3,Qsp,Qsp1,Qsp2,Qsp3,Usp,Usp1,Usp2,Usp3
  integer j

! Allocate
  allocate(  phihist(phibins) )
  allocate( dphihist(phibins) )

! Initialise
  phihist = 0.0

! Extract
  call simple_hist(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,phibins,phihist)
  call simple_hist(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,phibins,phihist)
  call simple_hist(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,phibins,phihist)
  
! Set errors
  dphihist = sqrt( phihist )

! Convert to count rate
  phihist  =  phihist / Texp
  dphihist = dphihist / Texp

! Sum the entire histogram to get the total count rate
  rate = 0.0
  do j = 1,phibins
     rate = rate + phihist(j)
  end do
  write(*,*)"Total count rate (c/s) = ",rate
  
! Extract spurious polarisation (SP) Stokes parameters
  call sp_stokes(evtunit1,nevts1,pilo,pihi,Nsp1,Qsp1,Usp1)
  call sp_stokes(evtunit2,nevts2,pilo,pihi,Nsp2,Qsp2,Usp2)
  call sp_stokes(evtunit3,nevts3,pilo,pihi,Nsp3,Qsp3,Usp3)

! Sum SP Stokes parameters
  Nsp = Nsp1 + Nsp2 + Nsp3
  Qsp = Qsp1 + Qsp2 + Qsp3
  Usp = Usp1 + Usp2 + Usp3

! Normalise SP Stokes parameters
  qsp = Qsp / Nsp
  usp = Usp / Nsp

! Subtract SP
  do j = 1,phibins

     !Calculate SP correction
     corr  = qsp * cos( 2.0*phi(j) ) + usp * sin( 2.0*phi(j) )
     corr  = corr * rate
     corr  = corr / real(phibins)

     !Subtract SP correction
     phihist(j) = phihist(j) - corr
     
  end do
  
  return
end subroutine modhist
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine fullband_sr(evtunit_sub,evtunit_ref,pilo,pihi,nn,tstart,dt)
! Outputs: lc(nn), lcsub(nn), lcref(nn) in shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit_sub,evtunit_ref,pilo,pihi,nn
  double precision tstart,dt
  
! Allocate arrays
  allocate( lc(nn)   )
  allocate( lcsub(nn)   )
  allocate( lcref(nn)   )

! Initialise arrays
  lc    = 0.0
  lcsub = 0.0
  lcref = 0.0

! Extract light curves
  call simplest_lc(evtunit_sub,nevts_sub,tstart,dt,pilo,pihi,nn,lcsub)
  call simplest_lc(evtunit_ref,nevts_ref,tstart,dt,pilo,pihi,nn,lcref)
  lc    = lcsub + lcref

! Convert to count rate
  lc    = lc / dt
  lcsub = lcsub / dt
  lcref = lcref / dt  

  return
end subroutine fullband_sr
!-----------------------------------------------------------------------




  
!-----------------------------------------------------------------------
subroutine fullbandlc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,tstart,dt,config)
! config determines what DUs are correlated.
! config = 1: DU3+DU2 vs DU1
! config = 2: DU3+DU1 vs DU2
! config = 3: DU2+DU1 vs DU3
!
! Outputs: lc(nn), lcsub(nn), lcref(nn) in shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,config
  double precision tstart,dt

! Allocate arrays
  allocate( lc(nn)   )
  allocate( lcsub(nn)   )
  allocate( lcref(nn)   )

! Initialise arrays
  lc    = 0.0
  lcsub = 0.0
  lcref = 0.0

! Extract light curves
  if( config .eq. 1 )then
     call simplest_lc(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,lc)
     call simplest_lc(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,lc)
     lcsub = lc
     call simplest_lc(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,lcref)
     lc    = lc + lcref
  else if( config .eq. 2 )then
     call simplest_lc(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,lc)
     call simplest_lc(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,lc)
     lcsub = lc
     call simplest_lc(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,lcref)
     lc    = lc + lcref
  else if( config .eq. 3 )then
     call simplest_lc(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,lc)
     call simplest_lc(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,lc)
     lcsub = lc
     call simplest_lc(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,lcref)
     lc    = lc + lcref
  end if
  
! Convert to count rate
  lc    = lc / dt
  lcsub = lcsub / dt
  lcref = lcref / dt
  
  return
end subroutine fullbandlc
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine fullbandlc_qu(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,tstart,dt,won)
! Outputs: lc(nn), qlc(nn), dqlc(nn), ulc(nn), dulc(nn) in shared_arrays
! nchn,ECHN(nchn),AEmuEi(nchn) and Aeffi(nchn) passed via shared_arrays
  
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,won
  double precision tstart,dt

! Allocate arrays
  allocate(  lc(nn)  )
  allocate( Ilc(nn)  )
  allocate( varIlc(nn) )
  allocate(  Qlc(nn) )
  allocate( varQlc(nn) )
  allocate(  Ulc(nn) )
  allocate( varUlc(nn) )
  
! Initialise arrays
  lc    = 0.0
  Ilc   = 0.0
  Qlc   = 0.0
  Ulc   = 0.0

! Extract light curves
  call simplest_lc_qu(evtunit1,nevts1,tstart,dt,pilo,pihi,nchn,AEmuEi,Aeffi,nn,won,&
       lc,ilc,qlc,ulc,varIlc,varQlc,varUlc)
  call simplest_lc_qu(evtunit2,nevts2,tstart,dt,pilo,pihi,nchn,AEmuEi,Aeffi,nn,won,&
       lc,ilc,qlc,ulc,varIlc,varQlc,varUlc)
  call simplest_lc_qu(evtunit3,nevts3,tstart,dt,pilo,pihi,nchn,AEmuEi,Aeffi,nn,won,&
       lc,ilc,qlc,ulc,varIlc,varQlc,varUlc)
  
! Convert to count rate
  lc   = lc   / dt
  
  return
end subroutine fullbandlc_qu
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine multiphisr(evtunit_sub,pilo,pihi,nn,phibins,tstart,dt)
! Extracts phibins light curves from the sub event file
! Outputs: lcphi(nn,phibins) in shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit_sub,pilo,pihi,nn,phibins
  double precision tstart,dt
  
! Allocate arrays
  allocate( lcphi(nn,phibins)   )

! Initialise arrays
  lcphi = 0.0

! Extract light curves
  call simplest_lc_phi(evtunit_sub,nevts_sub,tstart,dt,pilo,pihi,nn,phibins,lcphi)

! Convert to count rate
  lcphi = lcphi / dt
  
  return
end subroutine multiphisr
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine multiphilc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,tstart,dt,config)
! Extracts phibins light curves from the two DUs that are not DU(config)
! Outputs: lcphi(nn,phibins) in shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,config
  double precision tstart,dt

! Allocate arrays
  allocate( lcphi(nn,phibins)   )

! Initialise arrays
  lcphi = 0.0

! Extract light curves
  if( config .eq. 1 )then
     call simplest_lc_phi(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,phibins,lcphi)
     call simplest_lc_phi(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,phibins,lcphi)
  else if( config .eq. 2 )then
     call simplest_lc_phi(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,phibins,lcphi)
     call simplest_lc_phi(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,phibins,lcphi)
  else if( config .eq. 3 )then
     call simplest_lc_phi(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,phibins,lcphi)
     call simplest_lc_phi(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,phibins,lcphi)
  end if

! Convert to count rate
  lcphi = lcphi / dt
  
  return
end subroutine multiphilc
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine multiphilc_eff(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,tstart,dt)
! Extracts phibins DU1+DU2 light curves
! Outputs: lcphi(nn,phibins) in shared_arrays
! Uses effective modulation angle: phi = 0.5 * atan2(u,q), where u and q are level 2
! Stokes parameters.
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins
  double precision tstart,dt

! Allocate arrays
  allocate( lcphi(nn,phibins)   )

! Initialise arrays
  lcphi = 0.0

! Extract light curves
  call simplest_lc_phieff(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,phibins,lcphi)
  call simplest_lc_phieff(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,phibins,lcphi)
  
! Convert to count rate
  lcphi = lcphi / dt
  
  return
end subroutine multiphilc_eff
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine multiphilc_opp(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins,tstart,dt)
! Extracts phibins DU3 light curves
! Outputs: lcphi3(nn,phibins) in shared_arrays
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn,phibins
  double precision tstart,dt

! Allocate arrays
  allocate( lcphi3(nn,phibins)   )

! Initialise arrays
  lcphi3 = 0.0

! Extract light curves
  call simplest_lc_phi(evtunit3,nevts1,tstart,dt,pilo,pihi,nn,phibins,lcphi3)
  
! Convert to count rate
  lcphi3 = lcphi3 / dt
  
  return
end subroutine multiphilc_opp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine twophi_lc(evtunit1,evtunit2,evtunit3,pilo,pihi,nn,psi0,tstart,dt)
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3,pilo,pihi,nn
  double precision psi0,tstart,dt

! Allocate arrays
  allocate( lc_lo_12(nn)   )
  allocate( lc_hi_3(nn)   )
  allocate( lc_lo_3(nn)   )
  allocate( lc_hi_12(nn)   )

! Initialise arrays
  lc_lo_12 = 0.0
  lc_hi_3  = 0.0
  lc_lo_3  = 0.0
  lc_hi_12 = 0.0

! Extract light curves
  !DU1+2 light curves
  call simplest_lc_twophi(evtunit1,nevts1,tstart,dt,pilo,pihi,nn,psi0,lc_hi_12,lc_lo_12)
  call simplest_lc_twophi(evtunit2,nevts2,tstart,dt,pilo,pihi,nn,psi0,lc_hi_12,lc_lo_12)
  !DU3 light curves
  call simplest_lc_twophi(evtunit3,nevts3,tstart,dt,pilo,pihi,nn,psi0,lc_hi_3,lc_lo_3)

! Convert to count rate
  lc_lo_12 = lc_lo_12 / dt
  lc_hi_3  = lc_hi_3  / dt
  lc_lo_3  = lc_lo_3  / dt
  lc_hi_12 = lc_hi_12 / dt

  return
end subroutine twophi_lc
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine seg_wrangler(tstart,dt,nseg,total_exp)
  use shared_arrays
  implicit none
  !Inputs
  double precision tstart,dt,total_exp
  integer nseg
  !Internal
  integer gti,kstart,kend,m,segs,j
  
! Calculate the number of segments in the light curve
  mm = 0
  total_exp = 0.0
  do gti = 1,ngti
     !The last time bin not in this GTI
     kstart = ceiling( ( gtistart(gti) - TSTART ) / dt )
     !The last time bin in this GTI
     kend   = floor(   ( gtiend(gti)   - TSTART ) / dt )
     kend   = max( kstart , kend )
     !Add the number of full segments in this GTI
     mm   = mm + floor( real(kend-kstart) / real(nseg) )
     !Add to sum for total exposure time
     total_exp = total_exp + gtiend(gti) - gtistart(gti)
  end do
  write(*,*)"Total exposure time (s) = ",total_exp
  write(*,*)"Number of good segments = ",mm
  write(*,*)"Exposure per segment (s) = ",nseg*dt
  write(*,*)"Useful exposure time (s) = ",mm*nseg*dt
  write(*,*)"Fraction of time used = ",mm*nseg*dt / total_exp

! Now find the time bin at the end of each segment
  allocate( segend(mm) )
  m = 0
  do gti = 1,ngti
     !The last time bin not in this GTI
     kstart = ceiling( ( gtistart(gti) - TSTART ) / dt )
     !The last time bin in this GTI
     kend   = floor(   ( gtiend(gti)   - TSTART ) / dt )
     kend   = max( kstart , kend )
     !Number of full segments in this GTI
     segs   = floor( real(kend-kstart) / real(nseg) )
     !Loop through these segments
     do j = 1,segs
        m = m + 1
        segend(m) = kstart + j*nseg
     end do
  end do
  
  return
end subroutine seg_wrangler
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine GTIwrangler(evtunit1,evtunit2,evtunit3,tstart,tstop,Telapse)
  use shared_arrays
  implicit none
  !Inputs
  integer evtunit1,evtunit2,evtunit3
  double precision tstart,tstop,Telapse
  !Internal
  integer master_ngti,gti,status,gtiunit
  
! Read in GTIs
  call getngti(evtunit1,ngti1)
  allocate( gtistart1(ngti1) )
  allocate( gtiend1(ngti1) )
  call getgti(evtunit1,ngti1,gtistart1,gtiend1)

  call getngti(evtunit2,ngti2)
  allocate( gtistart2(ngti2) )
  allocate( gtiend2(ngti2) )
  call getgti(evtunit2,ngti2,gtistart2,gtiend2)

  call getngti(evtunit3,ngti3)
  allocate( gtistart3(ngti3) )
  allocate( gtiend3(ngti3) )
  call getgti(evtunit3,ngti3,gtistart3,gtiend3)

! Trim GTIs for floating point accuracy
  gtistart1 = gtistart1 - tstart
  gtistart2 = gtistart2 - tstart
  gtistart3 = gtistart3 - tstart
  gtiend1   = gtiend1   - tstart
  gtiend2   = gtiend2   - tstart
  gtiend3   = gtiend3   - tstart
  Telapse   = tstop     - tstart

  !Calculate number of GTIs in the master list
  ngti = master_ngti(ngti1,ngti2,ngti3,gtistart1,gtiend1,gtistart2,&
       gtiend2,gtistart3,gtiend3,Telapse,1d-10)
  !write(*,*)"ngti (master)=",ngti
  !Allocate master GTI arrays
  allocate( gtistart(ngti) )
  allocate(   gtiend(ngti) )
  !Calculate master GTI list
  call get_master_gtis(ngti,ngti1,ngti2,ngti3,gtistart1,gtiend1,&
       gtistart2,gtiend2,gtistart3,gtiend3,Telapse,1d-10,gtistart,gtiend)

! Write out GTIs to study in graphical form
  !write(*,*)"tstart=",tstart
  !call GTIplot

! Add tstart back on
  gtistart = gtistart + tstart
  gtiend   = gtiend   + tstart
  
  return
end subroutine GTIwrangler
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine get_master_gtis(ngti,ngti1,ngti2,ngti3,gtistart1,gtiend1,gtistart2,gtiend2,gtistart3,gtiend3,&
     Telapse,tacc,gtistart,gtiend)
  implicit none
  integer ngti,ngti1,ngti2,ngti3
  double precision gtistart1(ngti1),gtiend1(ngti1),gtistart2(ngti2),gtiend2(ngti2)
  double precision gtistart3(ngti3),gtiend3(ngti3),Telapse,tacc
  double precision gtistart(ngti),gtiend(ngti),tlatest
  integer i1,i2,i3,counter
! Shift to GTIs that all start after t = 0
  i1 = 1
  do while( gtistart1(i1) .lt. 0.d0 )
     i1 = i1 + 1
  end do  
  i2 = 1
  do while( gtistart2(i2) .lt. 0.d0 )
     i2 = i2 + 1
  end do
  i3 = 1
  do while( gtistart3(i3) .lt. 0.d0 )
     i3 = i3 + 1
  end do
! Loop through master GTI list
  do counter = 1,ngti
     !Set the start of the current master GTI
     gtistart(counter) = max( gtistart1(i1), gtistart2(i2), gtistart3(i3) )
     !Set the end of the current master GTI
     gtiend(counter) = Telapse
     if( gtiend1(i1) .gt. gtistart(counter) ) gtiend(counter) = min( gtiend1(i1) , gtiend(counter) )
     if( gtiend2(i2) .gt. gtistart(counter) ) gtiend(counter) = min( gtiend2(i2) , gtiend(counter) )
     if( gtiend3(i3) .gt. gtistart(counter) ) gtiend(counter) = min( gtiend3(i3) , gtiend(counter) )
     !Move to a GTI that ends after the current master GTI ends
     do while( gtiend1(i1) - gtiend(counter) .lt. 1e-6 .and. i1 .lt. ngti1 )
        i1 = i1 + 1
     end do
     do while( gtiend2(i2) - gtiend(counter) .lt. 1e-6 .and. i2 .lt. ngti2 )
        i2 = i2 + 1
     end do
     do while( gtiend3(i3) - gtiend(counter) .lt. 1e-6 .and. i3 .lt. ngti3 )
        i3 = i3 + 1
     end do
     !Find the latest starting of the currently selected GTIs
     tlatest = max( gtistart1(i1), gtistart2(i2), gtistart3(i3) )
     !Move to a GTI that ends after tlatest
     do while( gtiend1(i1) - tlatest .lt. tacc .and. i1 .lt. ngti1 )
        i1 = i1 + 1
     end do
     do while( gtiend2(i2) - tlatest .lt. tacc .and. i2 .lt. ngti2 )
        i2 = i2 + 1
     end do
     do while( gtiend3(i3) - tlatest .lt. tacc .and. i3 .lt. ngti3 )
        i3 = i3 + 1
     end do
  end do
  return
end subroutine get_master_gtis
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function master_ngti(ngti1,ngti2,ngti3,gtistart1,gtiend1,gtistart2,gtiend2,&
     gtistart3,gtiend3,Telapse,tacc)
  implicit none
  integer master_ngti,ngti1,ngti2,ngti3
  double precision gtistart1(ngti1),gtiend1(ngti1),gtistart2(ngti2),gtiend2(ngti2)
  double precision gtistart3(ngti3),gtiend3(ngti3),Telapse,tacc
  integer i1,i2,i3,counter
  logical loop
  double precision gtistart,gtiend,tlatest
! Shift to GTIs that all start after t = 0
  i1 = 1
  do while( gtistart1(i1) .lt. 0.d0 )
     i1 = i1 + 1
  end do  
  i2 = 1
  do while( gtistart2(i2) .lt. 0.d0 )
     i2 = i2 + 1
  end do
  i3 = 1
  do while( gtistart3(i3) .lt. 0.d0 )
     i3 = i3 + 1
  end do
! Count number of master GTIs
  counter = 0
  loop    = .true.
  do while(loop)
     counter = counter + 1
     !Set the start of the current master GTI
     gtistart = max( gtistart1(i1), gtistart2(i2), gtistart3(i3) , 0.0d0 )
     !Set the end of the current master GTI
     gtiend = Telapse
     if( gtiend1(i1) .gt. gtistart ) gtiend = min( gtiend1(i1) , gtiend )
     if( gtiend2(i2) .gt. gtistart ) gtiend = min( gtiend2(i2) , gtiend )
     if( gtiend3(i3) .gt. gtistart ) gtiend = min( gtiend3(i3) , gtiend )
     !Decide if we've reached the end of the loop
     if( Telapse - gtiend .lt. tacc ) loop = .false.
     if( i1 .eq. ngti1 .and. i2 .eq. ngti2 .and. i3 .eq. ngti3 ) loop = .false.
     !Move to a GTI that ends after the current master GTI ends
     do while( gtiend1(i1) - gtiend .lt. tacc .and. i1 .lt. ngti1 )
        i1 = i1 + 1
     end do
     do while( gtiend2(i2) - gtiend .lt. tacc .and. i2 .lt. ngti2 )
        i2 = i2 + 1
     end do
     do while( gtiend3(i3) - gtiend .lt. tacc .and. i3 .lt. ngti3 )
        i3 = i3 + 1
     end do
     !Find the latest starting of the currently selected GTIs
     tlatest = max( gtistart1(i1), gtistart2(i2), gtistart3(i3) )
     !Move to a GTI that ends after tlatest
     do while( gtiend1(i1) - tlatest .lt. tacc .and. i1 .lt. ngti1 )
        i1 = i1 + 1
     end do
     do while( gtiend2(i2) - tlatest .lt. tacc .and. i2 .lt. ngti2 )
        i2 = i2 + 1
     end do
     do while( gtiend3(i3) - tlatest .lt. tacc .and. i3 .lt. ngti3 )
        i3 = i3 + 1
     end do
  end do
  master_ngti = counter
! Catch instances where the user-defined end time is before the last GTI
  if( gtistart .gt. gtiend ) master_ngti = counter - 1
  return
end function master_ngti  
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine GTIplot
  use shared_arrays
  integer status,gtiunit,gti
! Open GTI unit
  status = 0
  call ftgiou(gtiunit,status)
  open(gtiunit,file='gtis.dat')
  write(gtiunit,*)"skip on"
!Write GTIs
  do gti = 1,ngti1
     write(gtiunit,*)gtistart1(gti),1.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti1
     write(gtiunit,*)gtiend1(gti),1.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti2
     write(gtiunit,*)gtistart2(gti),2.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti2
     write(gtiunit,*)gtiend2(gti),2.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti3
     write(gtiunit,*)gtistart3(gti),3.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti3
     write(gtiunit,*)gtiend3(gti),3.0
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti
     write(gtiunit,*)gtistart(gti),4.0     
  end do
  write(gtiunit,*)"no no"
  do gti = 1,ngti
     write(gtiunit,*)gtiend(gti),4.0
  end do
  write(gtiunit,*)"no no"
  write(*,*)"GTIs plotted in gtis.dat. Use @GTIplot02.pco"
! Close GTI unit
  call ftclos(gtiunit, status)
  call ftfiou(gtiunit, status)
  return
end subroutine GTIplot  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine open_subref(evtunit_sub,evtunit_ref,evtfile_sub,evtfile_ref)
  implicit none
  integer evtunit_sub,evtunit_ref
  character (len=500) evtfile_sub,evtfile_ref
  integer status,readwrite,blocksize  
! Open event files with read-only access
  status = 0
  call ftgiou(evtunit_sub,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening event file sub..."
  readwrite = 0
  call ftopen(evtunit_sub,evtfile_sub,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened event file:",trim(evtfile_sub)
  status = 0
  call ftgiou(evtunit_ref,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening event file ref..."
  readwrite = 0
  call ftopen(evtunit_ref,evtfile_ref,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened event file:",trim(evtfile_ref)
  return
end subroutine open_subref
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine open_evts(evtunit1,evtfile1,evtunit2,evtfile2,evtunit3,evtfile3)
  implicit none
  integer evtunit1,evtunit2,evtunit3
  character (len=500) evtfile1,evtfile2,evtfile3
  integer status,readwrite,blocksize
! Open event files with read-only access
  status = 0
  call ftgiou(evtunit1,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(evtunit1,evtfile1,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"Opened event file:",trim(evtfile1)
  status = 0
  call ftgiou(evtunit2,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(evtunit2,evtfile2,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"Opened event file:",trim(evtfile2)
  status = 0
  call ftgiou(evtunit3,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(evtunit3,evtfile3,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"Opened event file:",trim(evtfile3)
  return
end subroutine open_evts
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine sp_stokes(unit,nevts,pilo,pihi,Nsp0,Qsp0,Usp0)
! Inputs:
! unit, nevts, pilo, pihi
! Outputs:
! Qsp0, Usp0   = spurious polarisation Stokes parameters
! Stokes parameters are in units of counts.
! Stokes parameters ARE initialised within this subroutine.  
  implicit none
  !Inputs
  integer unit,nevts,pilo,pihi
  !Outputs
  double precision Nsp0,Qsp0,Usp0
  !Internal
  integer status,qcol,ucol,qualcol,picol,k,quality,PIchn
  logical exact,anynull
  double precision qsp,usp
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'QSP',qcol,status)
  call ftgcno(unit,exact,'USP',ucol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'PI',picol,status)

! Sum up Stokes parameters
  Nsp0 = 0.0
  Qsp0 = 0.0
  Usp0 = 0.0
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Read in Stokes parameters
           call ftgcvd(unit,qcol,k,1,1,-1.0,qsp,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,usp,anynull,status)
           !Add to running sum
           Nsp0 = Nsp0 + 1.0
           Qsp0 = qsp0 + qsp
           Usp0 = usp0 + usp
        end if
     end if
  end do
  !write(*,*)"Qsp=",Qsp0
  !write(*,*)"Usp=",Usp0
  !write(*,*)"Nsp=",Nsp0
  !write(*,*)"nevts=",nevts

  return
end subroutine sp_stokes
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine stokes(unit,nevts,pilo,pihi,N0,Q0,U0)
! Inputs:
! unit, nevts, pilo, pihi
! Outputs:
! Q0, U0   = level 2 Stokes parameters
! Stokes parameters are in units of counts.
! Stokes parameters ARE initialised within this subroutine.  
  implicit none
  !Inputs
  integer unit,nevts,pilo,pihi
  !Outputs
  double precision N0,Q0,U0
  !Internal
  integer status,qcol,ucol,qualcol,picol,k,quality,PIchn
  logical exact,anynull
  double precision q,u
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'Q',qcol,status)
  call ftgcno(unit,exact,'U',ucol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'PI',picol,status)

! Sum up Stokes parameters
  N0 = 0.0
  Q0 = 0.0
  U0 = 0.0
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Read in Stokes parameters
           call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
           !Add to running sum
           N0 = N0 + 1.0
           Q0 = q0 + q
           U0 = u0 + u
        end if
     end if
  end do
  !write(*,*)"Q (level 2)=",Q0
  !write(*,*)"U (level 2)=",U0
  !write(*,*)"N =",N0
  !write(*,*)"nevts=",nevts

  return
end subroutine stokes
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine simplest_lc_phi(unit,nevts,t0,dt,pilo,pihi,nn,phibins,lcphi)
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, phibins
! Outputs:
! lcphi(nn,phibins)   = counts in each time bin in each phi bin
! ** lcphi is NOT initialised so we can easily add DUs together **
  implicit none
  integer unit,nevts,pilo,pihi,nn,phibins
  real lcphi(nn,phibins)
  double precision t0,dt
  integer status,tcol,picol,k,PIchn,i,phicol,qcol,quality,j
  logical exact,anynull
  double precision time,tevt,phi,pi
  character (len=50) EXNAME,comment
  pi = acos(-1.d0)
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
  
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'PHI',phicol,status)
  call ftgcno(unit,exact,'QUAL',qcol,status)
  
! Create light curve
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           i    = ceiling( tevt / dble(dt) )
           if( i .gt. 0 .and. i .le. nn )then
              !Determine phi bin
              call ftgcvd(unit,phicol,k,1,1,-1.0,phi,anynull,status)
              !j = ceiling( (phi+pi) * dble(phibins)/(2.0*pi) )
              j = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
              if( j .lt. 1 ) write(*,*)"Shit! j = ",j
              if( j .gt. phibins ) write(*,*)"Shit! j = ",j
              !Add to correct light curve
              lcphi(i,j)   = lcphi(i,j)   + 1.0
           end if
        end if
     end if
  end do
  
  return
end subroutine simplest_lc_phi
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine simplest_lc(unit,nevts,t0,dt,pilo,pihi,nn,lc)
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn  
! Outputs:
! lc(1:nn)   = counts in each time bin
! ** lc is NOT initialised so we can easily add DUs together **
  implicit none
  integer unit,nevts,pilo,pihi,nn
  real lc(nn)
  double precision t0,dt
  integer status,tcol,picol,k,PIchn,i
  logical exact,anynull
  double precision time,tevt
  character (len=50) EXNAME,comment
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  
! Create light curve
  do k = 1,nevts
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     if( pichn .ge. pilo .and. pichn .le. pihi )then
        call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
        tevt = time - t0
        i    = ceiling( tevt / dble(dt) )
        if( i .gt. 0 .and. i .le. nn )then
           lc(i)   = lc(i)   + 1.0
        end if
     end if
  end do
  
  return
end subroutine simplest_lc
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine simplest_lc_qu(unit,nevts,t0,dt,pilo,pihi,nchn,AEmuEi,Aeffi,nn,won,&
     lc,ilc,qlc,ulc,varIlc,varQlc,varUlc)
! Calculates I, Q and U light curves
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nchn, AEmuEi(nchn), Aeffi(nchn), nn , won
! Outputs:
! lc(1:nn)     = counts in each time bin
! ilc(1:nn)    = Stokes I for each time bin (units of counts / cm^2)
! Qlc(1:nn)    = Q for each time bin (units of counts / cm^2)
! Ulc(1:nn)    = U for each time bin (units of counts / cm^2)
! varQlc(1:nn) = Variance of I
! varQlc(1:nn) = Variance of Q
! varUlc(1:nn) = Variance of U
! ** lc, ilc, qlc and ulc (and errors) are NOT initialised so we can easily add DUs together **
  implicit none
! Inputs
  integer unit,nevts,pilo,pihi,nchn,nn,won
  real AEmuEi(nchn),Aeffi(nchn)
  double precision t0,dt
! Outputs
  real lc(nn),ilc(nn),qlc(nn),ulc(nn),varIlc(nn),varQlc(nn),varUlc(nn)
! Internal
  integer status,tcol,picol,qcol,ucol,k,PIchn,i,wcol
  real weight
  logical exact,anynull
  double precision time,tevt,q,u
  character (len=50) EXNAME,comment
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'Q',qcol,status)
  call ftgcno(unit,exact,'U',ucol,status)
  call ftgcno(unit,exact,'W_MOM',wcol,status)

! Catch mistakes on the weight option
  won = min( won , 1 )
  won = max( won , 0 )
  
! Create light curve
  write(*,*)"Calclating simplest_lc_qu light curve..."
  do k = 1,nevts
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     if( pichn .ge. pilo .and. pichn .le. pihi )then
        call ftgcve(unit,wcol,k,1,1,-1.0,weight,anynull,status)
        weight = weight * won + 1 - won
        call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
        tevt = time - t0
        i    = ceiling( tevt / dble(dt) )
        if( i .gt. 0 .and. i .le. nn )then
           !Add to I sum
           lc(i)     = lc(i)   + 1.0
           ! ilc(i)    = ilc(i)  + weight / Aeffi(PIchn+1)
           ! varIlc(i) = varIlc(i) + ( weight / Aeffi(PIchn+1) )**2

           ilc(i)    = ilc(i)  + weight
           varIlc(i) = varIlc(i) + weight**2
           
           !Add to Q and U sums
           call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
           ! qlc(i)    = qlc(i)    + weight   * q / AEmuEi(PIchn+1)
           ! varQlc(i) = varQlc(i) + ( weight * q / AEmuEi(PIchn+1) )**2
           ! ulc(i)    = ulc(i)    + weight   * u / AEmuEi(PIchn+1)
           ! varUlc(i) = varUlc(i) + ( weight * u / AEmuEi(PIchn+1) )**2

           qlc(i)    = qlc(i)    + weight*Aeffi(PIchn+1) * q / AEmuEi(PIchn+1)
           varQlc(i) = varQlc(i) + ( weight*Aeffi(PIchn+1) * q / AEmuEi(PIchn+1) )**2
           ulc(i)    = ulc(i)    + weight*Aeffi(PIchn+1)   * u / AEmuEi(PIchn+1)
           varUlc(i) = varUlc(i) + ( weight*Aeffi(PIchn+1) * u / AEmuEi(PIchn+1) )**2
           
        end if
     end if
  end do
  write(*,*)"...finished calculating I, Q, U light curve"
  
  return
end subroutine simplest_lc_qu
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
subroutine filter_count(unit,nevts,pilo,pihi,ngood)
!
! Counts how many events there will be after filtering
!
! Inputs:
! unit, nevts, pilo, pihi
!
! Outputs:
! ngood
!
  implicit none
  integer unit,nevts,pilo,pihi,nn,ngood
  integer status,picol,k,PIchn,qualcol
  integer quality
  logical exact,anynull
  integer ilo,ihi
  character (len=50) EXNAME,comment
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  
! Count up number of events after filtering
  ngood = 0
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           ngood = ngood + 1
        end if
     end if
  end do
  
  return
end subroutine filter_count
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine filter(unit,nevts,t0,dt,pilo,pihi,nn,ngood,tbins,pichns,qsps,usps)
!
! Filters the events and outputs a great big list of the time bins
! each event belongs to.
!
! MUST have already called filter_count
!  
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, ngood
!
! Outputs:
! tbins(ngood)
!
  implicit none
  integer unit,nevts,pilo,pihi,nn,ngood
  integer tbins(ngood),pichns(ngood)
  double precision qsps(ngood),usps(ngood)
  double precision t0,dt,qsp,usp
  integer status,tcol,picol,k,PIchn,i,qualcol
  integer quality,qcol,ucol
  logical exact,anynull
  double precision time,tevt
  integer ilo,ihi,j
  character (len=50) EXNAME,comment
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'QSP',qcol,status)
  call ftgcno(unit,exact,'USP',ucol,status)
  
! Filter events
  i = 0
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt      = time - t0
           j         = ceiling( tevt / dble(dt) )
           i         = i + 1
           tbins(i)  = j
           PIchns(i) = pichn
           !Read in spurious polarisation Stokes parameters
           call ftgcvd(unit,qcol,k,1,1,-1.0,qsp,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,usp,anynull,status)
           qsps(i) = qsp
           usps(i) = usp
        end if
     end if
  end do
  
  return
end subroutine filter
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine filterp(unit,nevts,t0,dt,pilo,pihi,nn,ngood,tbins,ta,pichns,qsps,usps)
!
! Filters the events and outputs a great big list of the time bins
! each event belongs to.
!
! MUST have already called filter_count
!  
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, ngood
!
! Outputs:
! tbins(ngood),ta(ngood),pichns(ngood),qsps(ngood),usps(ngood)
!
  implicit none
  integer unit,nevts,pilo,pihi,nn,ngood
  integer tbins(ngood),pichns(ngood)
  double precision ta(ngood),qsps(ngood),usps(ngood)
  double precision t0,dt,qsp,usp
  integer status,tcol,picol,k,PIchn,i,qualcol
  integer quality,qcol,ucol
  logical exact,anynull
  double precision time,tevt
  integer ilo,ihi,j
  character (len=50) EXNAME,comment
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'QSP',qcol,status)
  call ftgcno(unit,exact,'USP',ucol,status)
  
! Filter events
  i = 0
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt      = time - t0
           j         = ceiling( tevt / dble(dt) )
           i         = i + 1
           tbins(i)  = j
           ta(i)     = tevt
           PIchns(i) = pichn
           !Read in spurious polarisation Stokes parameters
           call ftgcvd(unit,qcol,k,1,1,-1.0,qsp,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,usp,anynull,status)
           qsps(i) = qsp
           usps(i) = usp
        end if
     end if
  end do
  
  return
end subroutine filterp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine fake_lc_phi(unit,nevts,t0,dt,pilo,pihi,nn,phibins,nseg,phase,&
     mu,PDpar,PApar)
! nseg, phase(nseg,mm), q, u
! input a model of q and u vs QPO phase, but first set them to constant  
  use shared_arrays
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, phibins
! Outputs:
! lcphi(nn,phibins)   = counts in each time bin in each phi bin
! ** lcphi is NOT initialised so we can easily add DUs together **
  implicit none
  integer mmax
  parameter (mmax=20)
  integer unit,nevts,pilo,pihi,nn,phibins,nseg
  real phase(nseg,mm)
  double precision t0,dt,q,u,dgenphi
  double precision mu,PDpar(mmax),PApar(mmax),dyda(mmax)
  double precision qsp,usp,CDF,xacc,PD,PA
  integer status,tcol,picol,k,PIchn,i,phicol,qcol,ucol,qualcol
  integer quality,j,m,maxit,idum
  logical exact,anynull
  double precision time,tevt,phi,phase_t
  integer ilo,ihi
  double precision tlo,Xlo,Xhi,Ylo,Yhi,DelY,DelX,Delg,omega,mup,psi,qtot,utot
  real ran2,getphase_t
  character (len=50) EXNAME,comment
  data idum/-28572/
  save idum
  
  maxit = 500
  xacc  = 1d-5
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'PHI',phicol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'QSP',qcol,status)
  call ftgcno(unit,exact,'USP',ucol,status)
  
! Create light curve
  write(*,*)"Calclating fake light curves..."
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           j    = ceiling( tevt / dble(dt) )
           !Determine segment
           m = segment_no(j)
           if( m .gt. 0 )then !If we're in a segment
              !Determine the time bin within the segment
              i = j - segend(m) + nseg
              !Determine instantaneous PD and PA
              phase_t = getphase_t(i,j,m,dt,nseg,mm,phase,tevt)
              call dwffunc(phase_t,PDpar,PD,dyda,5)
              PD = max( PD , 0.d0 )
              call dwffunc(phase_t,PApar,PA,dyda,5)
              q  = mu * PD * cos( 2.0*PA )
              u  = mu * PD * sin( 2.0*PA )              
              !Select a modulation angle
              !Add spurious polarisation to Stokes parameters
              call ftgcvd(unit,qcol,k,1,1,-1.0,qsp,anynull,status)
              call ftgcvd(unit,ucol,k,1,1,-1.0,usp,anynull,status)
              qtot = q + qsp
              utot = u + usp
              !Select a modulation angle
              CDF  = dble( ran2(idum) )
              mup  = min( 0.95d0 , sqrt( qtot**2 + utot**2 ) )
              psi  = 0.5 * atan2(utot,qtot)
              phi  = dgenphi(CDF,mup,psi,xacc,maxit)
              !Determine phi bin
              i = ceiling( (phi+0.5*pi) * dble(phibins)/pi )              
              !Add to correct light curve
              lcphi(j,i)   = lcphi(j,i)   + 1.0
           end if
        end if
     end if
  end do
  write(*,*)"...finished calculating fake light curves"
  
  return
end subroutine fake_lc_phi
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine pol_per_evt(dt,nseg,ngood,tbins,ta,nn,PDpar,PApar,phase,PD_t,PA_t)
! Calculates instantaneous PD and PA for each event.
! Inputs: unit,t0,dt,nseg,ngood,tbins(1:ngood),nn,PDpar(1:mmax),PApar(1:mmax),phase(nseg,mm)
! Outputs: PD_t(1:ngood), PA_t(ngood)
  use shared_arrays  !segment_no(nn),segend(mm),mm
  implicit none
  !Inputs
  integer nseg,ngood
  integer tbins(ngood),nn,mmax
  parameter (mmax=20)
  double precision dt,PDpar(mmax),PApar(mmax),ta(ngood)
  real phase(nseg,mm)
  !Outputs
  double precision PD_t(ngood),PA_t(ngood)
  !Internal
  integer k,j,m,i
  double precision phase_t,time,tevt,PD,PA,dyda(mmax)
  real getphase_t
  
! Create light curve
  write(*,*)"Calclating instantaneous PD and PD..."
  do k = 1,ngood
     j = tbins(k)
     !Determine segment
     m = segment_no(j)
     if( m .gt. 0 )then !If we're in a segment

        !Get event time
        tevt = ta(k)
        
        !Determine the time bin within the segment
        i = j - segend(m) + nseg
        
        !Determine instantaneous PD and PA
        phase_t = getphase_t(i,j,m,dt,nseg,mm,phase,tevt)
        
        call dwffunc(phase_t,PDpar,PD,dyda,5)
        PD = max( PD , 0.d0 )
        call dwffunc(phase_t,PApar,PA,dyda,5)
        PD_t(k) = PD
        PA_t(k) = PA

        ! if( m .eq. 1 ) write(701,*)tevt,phase_t,PD,PA

     else
        !If not in a segment, give average PD and PA
        PD_t(k) = PDpar(1)
        PA_t(k) = PApar(1)
     end if
  end do
  write(*,*)"...finished calculating instantaneous PD and PD"

  return

end subroutine pol_per_evt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function getphase_t(i,j,m,dt,nseg,mm,phase,tevt)
!Input: i,j,dt,nseg,mm,phase(1:nseg,1:mm),tevt
!Output: getphase_t
  implicit none
  !Inputs
  integer i,j,m,nseg,mm
  double precision dt,tevt
  real phase(nseg,mm)
  !Outputs
  real getphase_t
  !Internal
  integer ihi,ilo
  double precision tlo,Xlo,Xhi,Ylo,Yhi,DelY,DelX,Delg,omega,pi
  parameter (pi=3.141592653589793)
  real phase_t
  !Determine QPO phase
  !Find time bins that end after and before tevt
  if( i .gt. 1 )then
     ihi = i
     ilo = i-1
     tlo = (j-1) * dble(dt)
  else
     ihi = i+1
     ilo = i
     tlo = j * dble(dt)
  end if
  !Determine QPO phase difference between ilo and ihi
  !in a way that is immune to the phase ambiguity
  Xlo = cos( 2.0*pi * phase(ilo,m) )
  Xhi = cos( 2.0*pi * phase(ihi,m) )
  Ylo = sin( 2.0*pi * phase(ilo,m) )
  Yhi = sin( 2.0*pi * phase(ihi,m) )
  DelY = Yhi*Xlo - Xhi*Ylo
  DelX = Ylo*Yhi + Xhi*Xlo
  Delg = atan2( DelY , DelX ) / (2.0*pi)
  !Interpolate phase and keep on 0 to 1 interval
  omega = Delg / dt
  phase_t = phase(ilo,m) + omega * ( tevt - tlo )
  if( phase_t .le. 0.0 ) phase_t = phase_t + 1.0
  if( phase_t .gt. 1.0 ) phase_t = phase_t - 1.0
  getphase_t = phase_t
  return
end function getphase_t
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine fake_phi_filt(ngood,tbins,pichns,qsps,usps,nn,phibins,PD_t,PA_t,&
     spinc,nchn,muEi,rar,lcfilt)
!
! Calculates phibins light curves lcfilt(1:nn,1:phibins) from the
! filtered event lists and input PD and PA modulations
! by generating fake modulation angles.
! lcfilt is NOT initialised so we can easily add DUs together **
!
! rar(1:ngood) is an array of random numbers
!  
  implicit none
  integer mmax
  parameter (mmax=20)
  !Input
  integer ngood,nn,phibins,nchn
  real muEi(nchn),rar(ngood)
  integer tbins(ngood),pichns(ngood)
  double precision qsps(ngood),usps(ngood)
  double precision PD_t(ngood),PA_t(ngood)
  logical spinc
  !Output
  real lcfilt(nn,phibins)
  !Internal
  integer k,j,i,maxit
  double precision q,u,mu,mup,psi,dgenphi,qsp,usp,CDF,xacc
  double precision qtot,utot,phi,pi
  pi = acos(-1.d0)

  maxit = 1000
  xacc  = 1d-5
  
  do k = 1,ngood

     !Determine time bin
     j = tbins(k)

     !Determine mu
     mu = muEi( pichns(k) )
     
     !Set Stokes parameters
     q  = mu * PD_t(k) * cos( 2.0*PA_t(k) )
     u  = mu * PD_t(k) * sin( 2.0*PA_t(k) )
     
     !Add spurious polarisation to Stokes parameters
     if( spinc )then
        qtot = q + qsps(k)
        utot = u + usps(k)
     else
        qtot = q
        utot = u
     end if
     
     !Select a modulation angle
     CDF  = dble( rar(k) )
     mup  = min( 0.95d0 , sqrt( qtot**2 + utot**2 ) )
     psi  = 0.5 * atan2(utot,qtot)
     phi  = dgenphi(CDF,mup,psi,xacc,maxit)
     
     !Determine phi bin
     i = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
     
     !Add to correct light curve
     lcfilt(j,i)   = lcfilt(j,i)   + 1.0
     
  end do

  return
end subroutine fake_phi_filt     
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine fake_phi_filt_mul(ngood,tbins,pichns,qsps,usps,nn,phibins,PD_t,PA_t,&
     spinc,nchn,muEi,rar,lcfilt)
!
! Calculates phibins light curves lcfilt(1:nn,1:phibins) from the
! filtered event lists and input PD and PA modulations
! by generating fake modulation angles.
! Uses new, non-union multiplicative modulation function.
! lcfilt is NOT initialised so we can easily add DUs together **
!
! rar(1:ngood) is an array of random numbers
!  
  implicit none
  integer mmax
  parameter (mmax=20)
  !Input
  integer ngood,nn,phibins,nchn
  real muEi(nchn),rar(ngood)
  integer tbins(ngood),pichns(ngood)
  double precision qsps(ngood),usps(ngood)
  double precision PD_t(ngood),PA_t(ngood)
  logical spinc
  !Output
  real lcfilt(nn,phibins)
  !Internal
  integer k,j,i,maxit
  double precision q,u,mu,dgenphimul,qsp,usp,CDF,xacc
  double precision phi,pi
  pi = acos(-1.d0)

  maxit = 1000
  xacc  = 1d-5
  
  do k = 1,ngood

     !Determine time bin
     j = tbins(k)

     !Determine mu
     mu = muEi( pichns(k) )
     
     !Set Stokes parameters
     q  = mu * PD_t(k) * cos( 2.0*PA_t(k) )
     u  = mu * PD_t(k) * sin( 2.0*PA_t(k) )
     
     !Include spurious polarisation
     if( spinc )then
        qsp = qsps(k)
        usp = usps(k)
     else
        qsp = 0.0
        usp = 0.0
     end if
     
     !Select a modulation angle
     CDF  = dble( rar(k) )     
     phi  = dgenphimul(CDF,q,u,qsp,usp,xacc,maxit)
     
     !Determine phi bin
     i = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
     
     !Add to correct light curve
     lcfilt(j,i)   = lcfilt(j,i)   + 1.0
     
  end do

  return
end subroutine fake_phi_filt_mul
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dgenphimul(CDF,q,u,qsp,usp,xacc,maxit)
!
! Generates random modulation angle phi for given muPD=mu*PD and PA.
! INPUTS:
! CDF     = number between 0 and 1 (cummulative distribution function)
! q       = Stokes q, q = mu PD cos(2PA) (fraction)
! u       = Stokes u, u = mu PD sin(2PA) (fraction)
! qsp     = spurious Stokes q (fraction)
! usp     = spurious Stokes u (fraction)
! xacc    = target accuracy
! maxit   = maximum number of iterations before giving up
! OUTPUTS:
! dgenphi = modulation angle (radians)
!
! Solves within xacc=1d-4, maxit=500, for muPD up to 0.97.
! Would need more iterations for the same accuracy for higher muPD.
!
  implicit none
  integer maxit
  double precision dgenphimul,CDF,q,u,qsp,usp,xacc
  integer j
  double precision pi,phi,phi_prev,reinterval,qtot,utot,norm
  double precision temp
  pi   = acos(-1.d0)
  norm = 1 / ( 2.0 + q*qsp + u*usp )
  qtot = q + qsp
  utot = u + usp
  
!First guess
  phi      = pi * (CDF-0.5)
  phi_prev = pi

! Run iteration loop to get true phi
  j = 1
  do while( abs(phi-phi_prev) .gt. xacc .and. j .lt. maxit )
     phi_prev = phi
     temp = pi*CDF - 0.5*pi - norm*qtot*sin(2.d0*phi)
     temp = temp + norm*utot*( cos(2.d0*phi) + 1.d0 )
     temp = temp - 0.25*norm*(q*qsp-u*usp) * sin(4.d0*phi)
     temp = temp + 0.25*norm*(q*usp+u*qsp) * ( cos(4.d0*phi) - 1.d0 )
     phi  = temp
     j = j + 1
  end do
  
! Catch failures
  if(  abs(phi-phi_prev) .gt. xacc )then
     write(*,*)"Warning! dgenphi failed to converge!"
     write(*,*)"Error is (rad) = ",abs(phi-phi_prev)
     write(*,*)"phi(rad)=",phi
  end if
  
! Convert back to -pi/2 to pi/2 interval
  phi = reinterval(phi,-0.5d0*pi,0.5d0*pi)
  
! Pass to the function name
  dgenphimul = phi

  return
end function dgenphimul
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine null_phi_filt(ngood,tbins,pichns,qsps,usps,nn,phibins,PD,PA,&
     spinc,nchn,muEi,rar,lcfilt)
!
! Calculates phibins light curves lcfilt(1:nn,1:phibins) from the
! filtered event lists by generating fake modulation angles.
! lcfilt is NOT initialised so we can easily add DUs together **
!
! rar(1:ngood) is an array of random numbers
!  
  implicit none
  !Input
  integer ngood,nn,phibins,nchn
  real muEi(nchn),rar(ngood)
  integer tbins(ngood),pichns(ngood)
  double precision qsps(ngood),usps(ngood),PD,PA
  logical spinc
  !Output
  real lcfilt(nn,phibins)
  !Internal
  integer k,j,i,maxit
  double precision q,u,mu,mup,psi,dgenphi,qsp,usp,CDF,xacc
  double precision qtot,utot,phi,pi
  pi = acos(-1.d0)

  maxit = 1000
  xacc  = 1d-5
  
  do k = 1,ngood
     j = tbins(k)
     !Determine mu
     mu = muEi( pichns(k) )
     !Set Stokes parameters
     q  = mu * PD * cos( 2.0*PA )
     u  = mu * PD * sin( 2.0*PA )              
     !Select a modulation angle
     !Add spurious polarisation to Stokes parameters
     if( spinc )then
        qtot = q + qsps(k)
        utot = u + usps(k)
     else
        qtot = q
        utot = u
     end if
     !Select a modulation angle
     CDF  = dble( rar(k) )  !dble( ran1(idum) )
     mup  = min( 0.95d0 , sqrt( qtot**2 + utot**2 ) )
     psi  = 0.5 * atan2(utot,qtot)
     phi  = dgenphi(CDF,mup,psi,xacc,maxit)
     !Determine phi bin
     i = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
     !Add to correct light curve
     lcfilt(j,i)   = lcfilt(j,i)   + 1.0
  end do

  return
end subroutine null_phi_filt     
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine null_lc_phi(unit,nevts,t0,dt,pilo,pihi,nn,phibins,PD,PA,spinc)
!
! Extracts light curves in different phi-bins by generating phi assuming
! a constant PD and PA.
!
!
  use shared_arrays
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, phibins, PD, PA, spinc
! via shared_arrays: nchn,muEi(1:nchn)  
! Outputs:
!
! ** lcphi is NOT initialised so we can easily add DUs together **
! via shared_arrays: lcphi(nn,phibins) = counts in each time bin in each phi bin
!
  implicit none
  integer unit,nevts,pilo,pihi,nn,phibins
  logical spinc
  double precision t0,dt,q,u
  double precision mu,mup,psi,dgenphi
  double precision qsp,usp,CDF,xacc,PD,PA
  integer status,tcol,picol,k,PIchn,i,phicol,qcol,ucol,qualcol
  integer quality,j,m,maxit,idum
  logical exact,anynull
  double precision time,tevt,phi
  integer ilo,ihi
  double precision qtot,utot
  real phase_t,ran2
  character (len=50) EXNAME,comment
  data idum/-28572/
  save idum
  
  maxit = 1000
  xacc  = 1d-5
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'PHI',phicol,status)
  call ftgcno(unit,exact,'QUAL',qualcol,status)
  call ftgcno(unit,exact,'QSP',qcol,status)
  call ftgcno(unit,exact,'USP',ucol,status)
  
! Create light curve
  write(*,*)"Calclating fake light curves..."
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qualcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           j    = ceiling( tevt / dble(dt) )
           !Determine mu
           mu = muEi(pichn)
           !Set Stokes parameters
           q  = mu * PD * cos( 2.0*PA )
           u  = mu * PD * sin( 2.0*PA )              
           !Select a modulation angle
           !Add spurious polarisation to Stokes parameters
           if( spinc )then
              call ftgcvd(unit,qcol,k,1,1,-1.0,qsp,anynull,status)
              call ftgcvd(unit,ucol,k,1,1,-1.0,usp,anynull,status)
              qtot = q + qsp
              utot = u + usp
           else
              qtot = q
              utot = u
           end if
           !Select a modulation angle
           CDF  = dble( ran2(idum) )
           mup  = sqrt( qtot**2 + utot**2 )
           psi  = 0.5 * atan2(utot,qtot)
           phi  = dgenphi(CDF,mup,psi,xacc,maxit)
           !Determine phi bin
           i = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
           !Add to correct light curve
           lcphi(j,i)   = lcphi(j,i)   + 1.0           
        end if
     end if
  end do
  write(*,*)"...finished calculating fake light curves"
  
  return
end subroutine null_lc_phi
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
function dmod_angle(CDF,mu,p0,psi0,xacc,maxit)
!
! DO NOT USE THIS FUNCTION. I THINK IT IS WRONG.
!
! Given CDF, a number between 0 and 1, the function returns the
! corresponding value of modulation angle, psi. An iteration scheme is
! used with a target accuracy xacc and maximum number of iterations
! maxit.
  implicit none
  integer maxit,j
  double precision dmod_angle,CDF,mu,p0,psi0,xacc
  double precision psi,psiprev,pi
  pi = acos(-1.0)
  !Initial guess for psi-psi0
  psi = 2.0 * pi * CDF
  psiprev = 2.0 * psi
  !Now iterate for psi-psi0
  j = 1
  do while( abs(psi-psiprev) .gt. xacc .and. j .lt. maxit )
     psiprev = psi
     psi = 2.0*pi*CDF - 0.5*mu*p0 * sin( 2.0 * psi )
     j = j + 1
  end do
  dmod_angle = psi + psi0
  if( abs(psi-psiprev) .gt. xacc )then
     write(*,*)"Warning! mod_angle did not converge"
     write(*,*)"Error is (rad) = ",abs(psi-psiprev)
  end if
  if( dmod_angle .gt. pi  ) dmod_angle = dmod_angle - 2.0*pi
  if( dmod_angle .le. -pi ) dmod_angle = dmod_angle + 2.0*pi
  return
end function dmod_angle
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dgenphi(CDF,muPD,PA,xacc,maxit)
!
! Generates random modulation angle phi for given muPD=mu*PD and PA.
! INPUTS:
! CDF     = number between 0 and 1 (cummulative distribution function)
! mu      = modulation factor (fraction)
! PD      = polarization degree (fraction)
! PA      = polarization angle (radians)
! xacc    = target accuracy
! maxit   = maximum number of iterations before giving up
! OUTPUTS:
! dgenphi = modulation angle (radians)
!
! Solves within xacc=1d-4, maxit=500, for muPD up to 0.97.
! Would need more iterations for the same accuracy for higher muPD.
!
  implicit none
  integer maxit
  double precision dgenphi,CDF,muPD,PA,xacc
  integer j
  double precision pi,phi,phi_prev,reinterval
  pi = acos(-1.d0)
  
!First guess
  phi      = pi * (CDF-0.5)
  phi_prev = pi

! Run iteration loop to get true phi
  j = 1
  do while( abs(phi-phi_prev) .gt. xacc .and. j .lt. maxit )
     phi_prev = phi
     phi      = pi*CDF - 0.5*pi - 0.5*muPD*sin(2.d0*phi)
     j = j + 1
  end do

! Catch failures
  if(  abs(phi-phi_prev) .gt. xacc )then
     write(*,*)"Warning! dgenphi failed to converge!"
     write(*,*)"Error is (rad) = ",abs(phi-phi_prev)
     write(*,*)"phi(rad)=",phi," muPD=",muPD," PA(rad)=",PA
  end if
  
! Add on PA and convert back to -pi/2 to pi/2 interval
  phi = phi + PA
  phi = reinterval(phi,-0.5d0*pi,0.5d0*pi)
  
! Pass to the function name
  dgenphi = phi

  return
end function dgenphi
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
function reinterval(phi,phimin,phimax)
!
! Places angle phi into the interval phimin -> phimax
!
  implicit none
  double precision reinterval,phi,phimin,phimax
  double precision range
  integer intervals
  reinterval = phi
  range      = phimax - phimin
  intervals  = max( 0 , ceiling( (phi-phimax) / range ) )
  reinterval = reinterval - intervals*range
  intervals  = max( 0 , ceiling( (phimin-phi) / range ) )
  reinterval = reinterval + intervals*range
  return
end function reinterval
!-----------------------------------------------------------------------




















!-----------------------------------------------------------------------
function PDFphi(muPD,PA,phi)
! Probability distribution function.
! muPD = modulation factor * polarization degree (fraction)
! PA   = polarization angle (radians, on interval -pi/2 to +pi/2)
! phi  = modulation angle (radians)  
  implicit none
  double precision PDFphi,muPD,PA,phi
  double precision pi
  parameter (pi=3.141592653589793)
  PDFphi = 1.0 + muPD * cos( 2.0*PA - 2.0*phi )
  PDFphi = PDFphi / pi
  return
end function PDFphi
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine simple_hist(unit,nevts,t0,dt,pilo,pihi,nn,phibins,hist)
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, phibins
! Outputs:
! hist(phibins)   = counts in each phi bin
! ** hist is NOT initialised so we can easily add DUs together **
! t0, dt and nn are input so I can make time cuts
  implicit none
  integer unit,nevts,pilo,pihi,nn,phibins
  real hist(phibins)
  double precision t0,dt
  integer status,tcol,picol,k,PIchn,i,phicol,qcol,quality,j
  logical exact,anynull
  double precision time,tevt,phi,pi
  character (len=50) EXNAME,comment
  pi = acos(-1.d0)
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'PHI',phicol,status)
  call ftgcno(unit,exact,'QUAL',qcol,status)
  
! Create light curve
  write(*,*)"Calclating histogram..."
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Filter for time interval
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           i    = ceiling( tevt / dt )
           if( i .gt. 0 .and. i .le. nn )then
              !Determine phi bin
              call ftgcvd(unit,phicol,k,1,1,-1.0,phi,anynull,status)
              j = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
              if( j .lt. 1 ) write(*,*)"Shit! j = ",j
              if( j .gt. phibins ) write(*,*)"Shit! j = ",j
              !Add to correct histogram bin
              hist(j)   = hist(j)   + 1.0
           end if
        end if
     end if
  end do
  write(*,*)"...finished calculating histogram"
  
  return
end subroutine simple_hist
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
subroutine simplest_lc_phieff(unit,nevts,t0,dt,pilo,pihi,nn,phibins,lcphi)
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, phibins
! Outputs:
! lcphi(nn,phibins)   = counts in each time bin in each phi bin
! ** lcphi is NOT initialised so we can easily add DUs together **
! Uses *effective* modulation angle, phi = 0.5 * atan2(u,q), where q and u
! are level 2 Stokes parameters.
  implicit none
  integer unit,nevts,pilo,pihi,nn,phibins
  real lcphi(nn,phibins)
  double precision t0,dt
  integer status,tcol,picol,k,PIchn,i,j,qcol,ucol
  logical exact,anynull
  double precision time,tevt,phi,pi,q,u
  character (len=50) EXNAME,comment
  pi = acos(-1.d0)
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'Q',qcol,status)
  call ftgcno(unit,exact,'U',ucol,status)
  
! Create light curve
  write(*,*)"Calclating light curves..."
  do k = 1,nevts
     !Filter for energy
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     if( pichn .ge. pilo .and. pichn .le. pihi )then
        !Determine time bin
        call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
        tevt = time - t0
        i    = ceiling( tevt / dble(dt) )
        if( i .gt. 0 .and. i .le. nn )then
           !Determine phi bin
           call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
           phi = 0.5 * atan2(u,q)
           j = ceiling( (phi+0.5*pi) * dble(phibins)/pi )
           !Add to correct light curve
           lcphi(i,j)   = lcphi(i,j)   + 1.0
        else
           write(*,*)"Warning! Not in a time bin!!!"
        end if
     end if
  end do
  write(*,*)"...finished calculating light curves"
  
  return
end subroutine simplest_lc_phieff
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine simplest_lc_twophi(unit,nevts,t0,dt,pilo,pihi,nn,psi0,lc_hi,lc_lo)
! Inputs:
! unit, nevts, t0, dt, pilo, pihi, nn, psi0
! Outputs:
! lc_hi(nn)   = counts in each time bin of the high pol light curve
! lc_lo(nn)   = counts in each time bin of the low pol light curve  
! ** lc_hi and lc_lo are NOT initialised so we can easily add DUs together **
  implicit none
  integer unit,nevts,pilo,pihi,nn
  real lc_hi(nn),lc_lo(nn)
  double precision t0,dt,psi0
  integer status,tcol,picol,k,PIchn,i,phicol,qcol,quality,j
  logical exact,anynull
  double precision time,tevt,phi,pi
  character (len=50) EXNAME,comment
  pi = acos(-1.d0)
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'PHI',phicol,status)
  call ftgcno(unit,exact,'QUAL',qcol,status)
  
! Create light curve
  write(*,*)"Calclating light curves..."
  do k = 1,nevts
     !Check quality flag
     call ftgcvj(unit,qcol,k,1,1,-1.0,quality,anynull,status)
     if( quality .eq. 1 )then
        !Filter for energy
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           !Determine time bin
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           i    = ceiling( tevt / dble(dt) )
           if( i .gt. 0 .and. i .le. nn )then
              !Read phi
              call ftgcvd(unit,phicol,k,1,1,-1.0,phi,anynull,status)
              !Convert to interval psi0 - pi/2 -> psi0 + pi/2
              if( phi .le. psi0-pi/2. ) phi = phi + pi
              if( phi .gt. psi0+pi/2. ) phi = phi - pi
              !Populate light curves
              if( phi .le. psi0+pi/4. .and. phi .gt. psi0-pi/4. )then
                 lc_hi(i) = lc_hi(i) + 1.0
              else
                 lc_lo(i) = lc_lo(i) + 1.0
              end if
           else
              write(*,*)"Warning! Not in a time bin!!!"
           end if
        end if
     end if
  end do
  write(*,*)"...finished calculating light curves"
  
  return
end subroutine simplest_lc_twophi
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
function wiggles(f,tau,rin)
! Function returns the wiggly part of the deadtime power spectrum, w(f).    
! For Leahy normalisation, the Poisson noise is then P(f) = 2.0 * w(f).
! Inputs:
! f    = Fourier frequency
! tau  = Deadtime
! rin  = Intrinsic count rate (i.e. count rate in the absence of deadtime)
  implicit none
  real wiggles,f,tau,rin
  real omega,pi,cosfac,sinfac,w
  pi = acos(-1.0)
  !Useful definitions
  omega = 2.0*pi*f
  cosfac = 1.0 - cos(omega*tau)
  sinfac = sin(omega*tau)
  !Calculate the function
  w = cosfac + omega/rin * sinfac
  w = w / ( cosfac**2 + ( sinfac+omega/rin)**2 )
  w = 2.0 * w
  w = 1.0 - w
  !Result
  wiggles = w
  return
end function wiggles
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
      SUBROUTINE ourfour1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!-----------------------------------------------------------------------

      
      

!-----------------------------------------------------------------------
subroutine mypow(nseg,C,nn,dt,mm,segend,norm,P,dP,mu,wn)
! INPUT:
! nseg        Number of time bins in a segment
! C(nn)       Time series to transform
! nn          Length of time series
! dt          Time interval
! mm          Number of segments
! segend(mm)  Last time bin in each segment
! norm        |norm|=1: fractional rms normalisation
!             |norm|=2: absolute rms normalisation
!             norm<0: subtract white noise; norm>0: don't subtract
! OUTPUT:
! P(nseg/2)   Output power spectrum
! dP(nseg/2)  Error on power spectrum
! mu          Mean count rate
! wn          White noise level
  implicit none
  integer nseg,nn,mm,norm,segend(mm)
  real C(nn),dt,P(nseg/2),dP(nseg/2)
  real df,mu,Cseg(nseg),wn,per(nseg/2),f
  integer i,j,m
  df = 1. / ( real(nseg) * dt )
  !Calculate average (non white noise subtracted) power spectrum
  !in absolute rms normalisation
  P  = 0.0
  mu = 0.0
  do m = 1,mm
     do i = 1,nseg
        j       = i + segend(m) - nseg
        Cseg(i) = C(j)
        mu      = mu + Cseg(i)
     end do
     call nperiodogram(dt,nseg,Cseg,per)
     P = P + per
  end do
  mu = mu / float(nseg*mm)
  wn = 2. * mu
  P  = P / float(mm)
  dP = P / sqrt( float(mm) )
  !Subtract white noise
  if( norm .lt. 0 ) P = P - wn
  !Convert to fractional rms normalisation
  if( abs(norm) .eq. 1 )then
    P  = P  / mu**2.0
    dP = dP / mu**2.0
    wn = wn / mu**2.0
  end if
  RETURN
END subroutine mypow
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine myco(nseg,Ca,Cb,nn,dt,mm,segend,norm,Pco,dPco)
! INPUT:
! nseg          Number of time bins in a segment
! Ca(nn)        Time series from detector combination A
! Cb(nn)        Time series from detector combination B
! nn            Length of time series
! dt            Time interval
! mm            Number of segments
! segend(mm)    Last time bin in each segment
! norm          |norm|=1: fractional rms normalisation
!               |norm|=2: absolute rms normalisation
! OUTPUT:
! Pco(nseg/2)   Output co-spectrum
! dPco(nseg/2)  Error on co-spectrum
! mua           Mean count rate from detector combination A
! mub           Mean count rate from detector combination B
  implicit none
  integer nseg,nn,mm,norm,segend(mm)
  real Ca(nn),Cb(nn),dt,Pco(nseg/2),dPco(nseg/2)
  real mua,mub,Caseg(nseg),Cbseg(nseg)
  real Pa(nseg/2),Pb(nseg/2),rc(nseg/2),ic(nseg/2)
  integer i,j,m
  !Calculate average (non white noise subtracted) power spectrum
  !in absolute rms normalisation
  Pco = 0.0
  Pa  = 0.0
  Pb  = 0.0
  mua = 0.0
  mub = 0.0
  do m = 1,mm
     do i = 1,nseg
        j        = i + segend(m) - nseg
        Caseg(i) = Ca(j)
        mua      = mua + Caseg(i)
        Cbseg(i) = Cb(j)
        mub      = mub + Cbseg(i)        
     end do
     call ncperiodogram(dt,nseg,Caseg,Cbseg,rc,ic)
     Pco = Pco + rc
     call ncperiodogram(dt,nseg,Caseg,Caseg,rc,ic)
     Pa  = Pa + rc
     call ncperiodogram(dt,nseg,Cbseg,Cbseg,rc,ic)
     Pb  = Pb + rc
  end do
  mua  = mua / real(nseg*mm)
  mub  = mub / real(nseg*mm)
  Pco  = Pco / real(mm)
  Pa   = Pa  / real(mm)
  Pb   = Pb  / real(mm)
  dPco = sqrt( Pa * Pb / real(2*mm) )
  !Convert to fractional rms normalisation
  if( abs(norm) .eq. 1 )then
    Pco  = Pco  / (mua*mub)
    dPco = dPco / (mua*mub)
  end if
  RETURN
end subroutine myco
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mybiphase(nseg,Cs,Cr,nn,dt,mm,segend,nuqpo,hwhm,nit,biphase,dbiphase)
  implicit none
! Inputs
  integer nseg,nn,mm,segend(mm),nit
  real Cs(nn),Cr(nn),dt,nuqpo,hwhm
! Outputs
  real biphase,dbiphase
! Internal
  integer m,i,j,ilo,ihi
  real df,Cs_seg(nseg),Cr_seg(nseg),pi
  complex B(nseg/4),Bqpo,Bm(nseg/4,mm),Bit(nseg/4)
  real ran2,biphase_it
  integer idum,it
  pi = acos(-1.0)

! Set frequency bins to average over
  df  = 1.0 / ( real(nseg)*dt )
  ilo = nint( ( nuqpo-hwhm) / df )
  ihi = nint( ( nuqpo+hwhm) / df )
  
! Calculate estimate
  !Calculate bispectrum
  B = cmplx( 0.0 , 0.0 )
  do m = 1,mm
     do i = 1,nseg
        j        = i + segend(m) - nseg
        Cs_seg(i) = Cs(j)
        Cr_seg(i) = Cr(j)
     end do
     call cmplx_abispec(dt,nseg,Cs_seg,Cr_seg,Bm(:,m))
     B(:) = B(:) + Bm(:,m)
  end do
  !Average over QPO frequency range
  Bqpo = cmplx( 0.0 , 0.0 )
  do i = ilo,ihi
     Bqpo = Bqpo + B(i)
  end do
  !Take the argument
  biphase = atan2( aimag(Bqpo) , real(Bqpo) )
  
! Bootstrap errors
  !Initialise random number chain
  idum = -387573
  !Run through iterations
  dbiphase = 0.0
  do it = 1,nit
     !Calculate auto-bispectrum
     Bit = 0.0
     do j = 1,mm
        m = ceiling( real(mm) * ran2(idum) )
        m = min( m , mm )
        Bit(:) = Bit(:) + Bm(:,m)
     end do
     !Average over QPO frequency range
     Bqpo = cmplx( 0.0 , 0.0 )
     do i = ilo,ihi
        Bqpo = Bqpo + Bit(i)
     end do
     !Take the argument
     biphase_it = atan2( aimag(Bqpo) , real(Bqpo) )
     dbiphase = dbiphase + ( biphase - biphase_it )**2
  end do
  dbiphase = sqrt( dbiphase / real(nit) )
  
  return
end subroutine mybiphase
!-----------------------------------------------------------------------



!------------------------------------------------------------------------
subroutine cmplx_abispec(dt,n,st,rt,B)
  implicit none
! Input
  integer n
  real dt,st(n),rt(n)
! Output
  complex B(n/4)
! Internal
  integer j
  real datas(2*n),datar(2*n)
  complex Rf(n/2),Sf(n/2)

! Pack light curves into data array
  do j = 1,n         
     datas(2*j-1) = st(j)
     datas(2*j)   = 0.0
     datar(2*j-1) = rt(j)
     datar(2*j)   = 0.0
  end do

! Fourier transform
  call ourfour1(datas,n,1)
  call ourfour1(datar,n,1)

! Pack data array into complex FT arrays
  do j = 1,n/2
     Rf(j) = cmplx( datar(2*j+1) , datar(2*j+2) )
     Sf(j) = cmplx( datas(2*j+1) , datas(2*j+2) )
  end do

! Calculate bispectrum
  do j = 1,n/4
     B(j) = Rf(j) * Rf(j) * conjg( Sf(2*j) )
     B(j) = B(j) * 2. * dt / real(n)
  end do

  return
end subroutine cmplx_abispec
!------------------------------------------------------------------------







!-----------------------------------------------------------------------
subroutine mycross(nseg,Ca,Cb,nn,dt,mm,segend,norm,ReG,dReG,ImG,dImG)
! INPUT:
! nseg          Number of time bins in a segment
! Ca(nn)        Time series from detector combination A
! Cb(nn)        Time series from detector combination B
! nn            Length of time series
! dt            Time interval
! mm            Number of segments
! segend(mm)    Last time bin in each segment
! norm          |norm|=1: fractional rms normalisation
!               |norm|=2: absolute rms normalisation
! OUTPUT:
! ReG(nseg/2)   Output real part of the cross-spectrum
! dReG(nseg/2)  Error on real part
! ImG(nseg/2)   Output imaginary part of the cross-spectrum
! dImG(nseg/2)  Error on imaginary part
! mua           Mean count rate from detector combination A
! mub           Mean count rate from detector combination B
  implicit none
  integer nseg,nn,mm,norm,segend(mm)
  real Ca(nn),Cb(nn),dt,ReG(nseg/2),dReG(nseg/2)
  real ImG(nseg/2),dImG(nseg/2)
  real mua,mub,Caseg(nseg),Cbseg(nseg)
  real rc(nseg/2),ic(nseg/2)
  integer i,j,m
  !Calculate average (non white noise subtracted) power spectrum
  !in absolute rms normalisation
  ReG = 0.0
  ImG = 0.0
  mua = 0.0
  mub = 0.0
  do m = 1,mm
     do i = 1,nseg
        j        = i + segend(m) - nseg
        Caseg(i) = Ca(j)
        mua      = mua + Caseg(i)
        Cbseg(i) = Cb(j)
        mub      = mub + Cbseg(i)        
     end do
     call ncperiodogram(dt,nseg,Caseg,Cbseg,rc,ic)
     ReG = ReG + rc
     ImG = ImG + ic
  end do
  mua  = mua / real(nseg*mm)
  mub  = mub / real(nseg*mm)
  ReG  = ReG / real(mm)
  ImG  = ImG / real(mm)
  !Now error
  dReG = 0.0
  dImG = 0.0
  do m = 1,mm
     do i = 1,nseg
        j        = i + segend(m) - nseg
        Caseg(i) = Ca(j)
        Cbseg(i) = Cb(j)
     end do
     call ncperiodogram(dt,nseg,Caseg,Cbseg,rc,ic)
     dReG = dReG + ( ReG - rc )**2
     dImG = dImG + ( ImG - ic )**2
  end do
  dReG = sqrt(dReG) / real(mm)
  dImG = sqrt(dImG) / real(mm)
  !Convert to fractional rms normalisation
  if( abs(norm) .eq. 1 )then
     ReG  = ReG  / (mua*mub)
     ImG  = ImG  / (mua*mub)
     dReG = dReG / (mua*mub)
     dImG = dImG / (mua*mub)
  end if
  RETURN
end subroutine mycross
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine nperiodogram(dt,n,at,per)
! Calculates the periodogram of the time series at(n) in absolute rms
! normalisation (uses four1 from Press et al 1992)
  implicit none
  integer n,j
  real at(n),per(n/2)
  real data(2*n),dt
  do j = 1,n
    data(2*j-1) = at(j)
    data(2*j)   = 0.0
  end do
  call ourfour1(data,n,1)
  do j = 1, n/2
    per(j) = data(2*j+1)**2.0 + data(2*j+2)**2.0
    per(j) = per(j) * 2. * dt / float(n)
  end do
  return
end subroutine nperiodogram
!------------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mycrossf(nseg,Cj,Cn,nnmax,dt,mm,segend,norm,Cross,sig)
! c INPUT:
! c nseg          Number of time bins in a segment
! c Cj(nnmax)     One time series to transform
! c Cn(nnmax)     The other time series to transform
! c nnmax         Maximum length of time series
! c dt            Time interval
! c mm            Number of segments
! c segend(mm)    Last time bin in each segment
! c norm          norm=1: fractional rms; norm=2: absolute rms
! c OUTPUT:
! c Cross(nseg/2) Complex cross-spectrum coefficient (+ve phase means j lags n)
! c sig(nseg/2)   Error on the cross-spectrum
  implicit none
  integer nseg,nnmax,mm,segend(nnmax/nseg),norm
  real Cj(nnmax),Cn(nnmax),dt,sig(nseg/2)
  complex Cross(nseg/2)
  integer m,i,j
  real Pj(nseg/2),Pn(nseg/2),Cjseg(nseg),Cnseg(nseg)
  real rc(nseg/2),ic(nseg/2),muj,mun
  Cross = 0.0
  Pj    = 0.0
  Pn    = 0.0
  sig   = 0.0
  muj   = 0.0
  mun   = 0.0
  !Calculate cross-spectrum in absolute rms normalisation
  do m = 1,mm
     do i = 1,nseg
        j = i + segend(m) - nseg
        Cjseg(i) = Cj(j)
        Cnseg(i) = Cn(j)
        muj      = muj + Cj(j)
        mun      = mun + Cn(j)
     end do
     call ncperiodogram(dt,nseg,Cjseg,Cnseg,rc,ic)
     do i = 1,nseg/2
        Cross(i) = Cross(i) + complex( rc(i) , ic(i) )
     end do
     call ncperiodogram(dt,nseg,Cjseg,Cjseg,rc,ic)
     Pj = Pj + rc
     call ncperiodogram(dt,nseg,Cnseg,Cnseg,rc,ic)
     Pn = Pn + rc
  end do
  muj   = muj   / float( nseg*mm )
  mun   = mun   / float( nseg*mm )
  write(*,*)"FPMA count rate=",muj
  write(*,*)"FPMB count rate=",mun
  Cross = Cross / float( mm )
  Pj    = Pj    / float( mm )
  Pn    = Pn    / float( mm )
  sig   = sqrt( Pn * Pj / ( 2.0 * float(mm) ) )
  !Convert to fractional rms normalisation
  if( abs(norm) .eq. 1 )then
     Cross = Cross / ( muj * mun )
     sig   = sig   / ( muj * mun )
  end if
  RETURN
END subroutine mycrossf
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dec_cross(nseg,mm,FFT_sub,FFT_ref,dt,jlo,jhi,ReG,ImG)
! "Deconstructed" cross spectrum.
! FFTs are entered and this code just multiplies together and normalises.
! No error is calculated, and it is automatically in fractional rms norm.
! +ve lag is subject lagging reference.
  implicit none
  integer nseg,mm,jlo,jhi
  real FFT_sub(2*nseg,mm),FFT_ref(2*nseg,mm),dt,ReG,ImG
  integer m,j,bins
  real rc,ic,meanref,meansub
! Average ReG and ImG over segments
  bins = jhi - jlo + 1
  ReG     = 0.0
  ImG     = 0.0
  meanref = 0.0
  meansub = 0.0
  do m = 1,mm
     meanref = meanref + FFT_ref(1,m) / real(nseg)
     meansub = meansub + FFT_sub(1,m) / real(nseg)
     do j = jlo,jhi
        rc  = FFT_sub(2*j+1,m)*FFT_ref(2*j+1,m)
        rc  = rc + FFT_sub(2*j+2,m)*FFT_ref(2*j+2,m)
        ReG = ReG + rc        
        ic  = FFT_sub(2*j+2,m)*FFT_ref(2*j+1,m)
        ic  = ic - FFT_sub(2*j+1,m)*FFT_ref(2*j+2,m)
        ImG = ImG + ic
     end do
  end do
  meanref = meanref / real(mm)
  meansub = meansub / real(mm)
  ReG  = ReG / real( mm * bins )
  ReG  = ReG * 2. * dt / real(nseg)
  ImG  = ImG / real( mm * bins )
  ImG  = ImG * 2. * dt / real(nseg)
! Convert to fractional rms^2
  ReG = ReG / ( meansub*meanref )
  ImG = ImG / ( meansub*meanref )
  return
end subroutine dec_cross  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dec_crossf(nseg,mm,FFT_sub,FFT_ref,dt,ReG,dReG,ImG,dImG)
! "Deconstructed" cross spectrum vs frequency.
! FFTs are entered and this code multiplies together and normalises.
! Error is calculated, and it is automatically in fractional rms norm.
! +vs lag is subject lagging reference.
  implicit none
  integer nseg,mm
  real FFT_sub(2*nseg,mm),FFT_ref(2*nseg,mm),dt,ReG(nseg/2),dReG(nseg/2)
  real ImG(nseg/2),dImG(nseg/2)
  integer m,j
  real rc,ic,meanref,meansub
  
! Average ReG and ImG over segments
  ReG     = 0.0
  ImG     = 0.0
  meanref = 0.0
  meansub = 0.0
  do m = 1,mm
     meanref = meanref + FFT_ref(1,m) / real(nseg)
     meansub = meansub + FFT_sub(1,m) / real(nseg)
     do j = 1,nseg/2
        rc  = FFT_sub(2*j+1,m)*FFT_ref(2*j+1,m)
        rc  = rc + FFT_sub(2*j+2,m)*FFT_ref(2*j+2,m)
        ReG(j) = ReG(j) + rc
        ic  = FFT_sub(2*j+2,m)*FFT_ref(2*j+1,m)
        ic  = ic - FFT_sub(2*j+1,m)*FFT_ref(2*j+2,m)        
        ImG(j) = ImG(j) + ic
     end do
  end do
  meanref = meanref / real(mm)
  meansub = meansub / real(mm)
  ReG  = ReG / real(mm)
  ImG  = ImG / real(mm)
  
! Calculate error
  dReG = 0.0
  dImG = 0.0
  do m = 1,mm
     do j = 1,nseg/2
        rc  = FFT_sub(2*j+1,m)*FFT_ref(2*j+1,m)
        rc  = rc + FFT_sub(2*j+2,m)*FFT_ref(2*j+2,m)
        dReG(j) = dReG(j) + ( ReG(j) - rc )**2
        ic  = FFT_sub(2*j+2,m)*FFT_ref(2*j+1,m)
        ic  = ic - FFT_sub(2*j+1,m)*FFT_ref(2*j+2,m)
        dImG(j) = dImG(j) + ( ImG(j) - ic )**2
     end do
  end do
  dReG = sqrt(dReG) / real(mm)
  dImG = sqrt(dImG) / real(mm)

! Convert to fractional rms^2
  ReG  = ReG * 2. * dt / real(nseg)
  ReG  = ReG / ( meansub*meanref )
  dReG = dReG * 2. * dt / real(nseg)
  dReG = dReG / ( meansub*meanref )
  ImG  = ImG * 2. * dt / real(nseg)
  ImG  = ImG / ( meansub*meanref )
  dImG = dImG * 2. * dt / real(nseg)
  dImG = dImG / ( meansub*meanref )
  
  return
end subroutine dec_crossf
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dec_pow(nseg,mm,FFT,dt,jlo,jhi,P)
! "Deconstructed" power spectrum.
! FFT is entered and this code just takes the mod squared and normalises.
! No error is calculated, and it is automatically in fractional rms norm.
! +vs lag is subject lagging reference.
  implicit none
  integer nseg,mm,jlo,jhi
  real FFT(2*nseg,mm),dt,P
  integer m,j,bins
  real per,mean
! Average ReG and ImG over segments
  bins = jhi - jlo + 1
  P    = 0.0
  mean = 0.0
  do m = 1,mm
     mean = mean + FFT(1,m) / real(nseg)
     do j = jlo,jhi
        per = FFT(2*j+1,m)**2 + FFT(2*j+2,m)**2
        P   = P + per
     end do
  end do
  mean = mean / real(mm)
  P    = P / real( mm * bins )
  P    = P * 2. * dt / real(nseg)
! Convert to fractional rms^2
  P    = P / mean**2
  return
end subroutine dec_pow
!-----------------------------------------------------------------------


!------------------------------------------------------------------------
subroutine FFTsegs(nseg,nn,lc,mm,segend,FFTa)
!
! Calculates FFT of each segment of light curve lc(nn)
! Output FFTa(nseg,mm) is in four1 format, so odd entries are real
! and even entries are imaginary.
!
  implicit none
  integer nseg,nn,mm
  integer segend(mm)
  real lc(nn),FFTa(2*nseg,mm)
  integer m,i,j
  real lc_seg(nseg)
  do m = 1,mm
     !Chop into current segment
     do i = 1,nseg
        j         = i + segend(m) - nseg
        lc_seg(i) = lc(j)
     end do
     !Calculate FFT of this segment
     call justFFT(nseg,lc_seg,FFTa(:,m))
  end do
  return
end subroutine FFTsegs
!------------------------------------------------------------------------


!------------------------------------------------------------------------
subroutine justFFT(n,at,FFTa)
! Fourier transforms at(1:n) and outputs in four1 format as FFTa(2*n).
! Even entries and imaginary and odd entries are real.
  implicit none
  integer n
  real at(n),FFTa(2*n)
  integer j
  !Fill FFTa() array
  do j = 1,n
     FFTa(2*j-1) = at(j)
     FFTa(2*j)   = 0.0
  end do
  !Do FFT
  call ourfour1(FFTa,n,1)
  return
end subroutine justFFT  
!------------------------------------------------------------------------



!------------------------------------------------------------------------
subroutine ncperiodogram(dt,n,ht,st,rc,ic)
! Calculates a cross spectrum between the time series ht(n) and st(n)
! In absolute rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! From Press et al (1992) DFT definition, this is H(\nu)S^*(\nu)
! p.g. 491 "time shifting" property: x(t-t0) = X(f) exp(i 2pi f t0)
! s(t) = delta(t) => s is a sharp flare at t=0
! h(t) = s(t-t0)  => h is a sharp flare at t=t0
! => H(f)S^*(f) = S(f) exp(i 2pi f t0) S^*(f)
! => H(f)S^*(f) = |S(f)|^2 exp(i 2pi f t0)
! which is a positive lag
  implicit none
  integer n,j
  real ht(n),st(n),dt
  real datah(2*n),datas(2*n),rc(n/2),ic(n/2)
  do j = 1,n         
     datah(2*j-1) = ht(j)
     datah(2*j)   = 0.0
     datas(2*j-1) = st(j)
     datas(2*j)   = 0.0
  end do
  call ourfour1(datah,n,1)
  call ourfour1(datas,n,1)
  do j = 1, n/2
     rc(j) = datah(2*j+1)*datas(2*j+1)
     rc(j) = rc(j) + datah(2*j+2)*datas(2*j+2)
     rc(j) = rc(j) * 2. * dt / real(n)     
     ic(j) = datah(2*j+2)*datas(2*j+1)
     ic(j) = ic(j) - datah(2*j+1)*datas(2*j+2)
     ic(j) = ic(j) * 2. * dt / real(n)
  end do  
  return
end subroutine ncperiodogram
!------------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine rebinE(earx,px,nex,ear,p,ne)
!General rebinning scheme, should be nice and robust - BUT IT FUCKING ISN'T
!i,nex,earx,px = input
!j,ne,ear,p    = output
  implicit none
  integer i,nex,j,ne,ilo,ihi
  real earx(0:nex),ear(0:ne),px(nex),p(ne),Ehigh,Elow,upper,lower
  real FRAC,Ej,Ei,pi,Ei2,pi2,grad,cons,Ehi,Elo,phi,plo
  logical interp
  ilo = 1
  do j = 1,ne
     do while( earx(ilo) .le. ear(j-1) .and. ilo .lt. nex )
        ilo = ilo + 1
     end do
     ihi = ilo
     do while( earx(ihi) .le. ear(j) .and. ihi .lt. nex )
        ihi = ihi + 1
     end do
     if( ihi .gt. ilo )then
        p(j) = 0.0
        do i = ilo,ihi
           lower = MAX( earx(i-1) , ear(j-1)  )
           upper = MIN( earx(i)   , ear(j)    )
           p(j) = p(j) + px(i) * ( upper - lower )
        end do
        p(j) = p(j) / ( ear(j) - ear(j-1) )
     else
        !Interpolate (or extrapolate)
        i = ilo
        Ei = 0.5 * ( earx(i) + earx(i-1) )
        if( Ei .gt. Ej ) i = ilo - 1
        i = max( i , 2     )
        i = min( i , nex-1 )
        Ej  = 0.5 * ( ear(j) + ear(j-1) )
        Ehi = 0.5 * ( earx(i+1) + earx(i)   )
        Elo = 0.5 * ( earx(i)   + earx(i-1) )
        phi = px(i+1)
        plo = px(i)
        p(j) = plo + (phi-plo)*(Ej-Elo)/(Ehi-Elo)
     end if
     if( ilo .gt. 1 ) ilo = ilo - 1
  end do
  RETURN   
END subroutine rebinE
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine channel(Elo,Ehi,nchn,echn,PIlo,PIhi)
! Convert energy min and max to PI channel min and max
  implicit none
  integer nchn,PIlo,PIhi
  real Elo,Ehi,echn(0:nchn)
  integer k
  k = 1
  do while( Elo .gt. ECHN(k-1) )
     k = k + 1
  end do
  PIlo = k - 1
  k = 1
  do while( ECHN(k) .le. Ehi )
     k = k + 1
  end do
  PIhi = k - 2
  return
end subroutine channel
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getNCHN(rmffile,nchn)
  implicit none
  character (len=500) rmffile
  integer nchn
  integer status,rmfunit,readwrite,blocksize
  character (len=500) comment
  
! Open rmf file with read-only access
  status = 0
  call ftgiou(rmfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(rmfunit,rmffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open rmf file'

! Shift to  extension "EBOUNDS"
  status = 0
  call ftmnhd(rmfunit,2,'EBOUNDS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EBOUNDS'

! Get number of rows in the table
  call ftgkyj(rmfunit,'NAXIS2',NCHN,comment,status)
  if(status .ne. 0) stop 'Cannot determine No of rows'

! Close mrf file and free up unit number
  call ftclos(rmfunit,status)
  call ftfiou(rmfunit,status)

  return
end subroutine getNCHN
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getECHN(rmffile,NCHN,ECHN)
  implicit none
  integer nchn
  character (len=500) rmffile
  real ECHN(0:NCHN)
  integer status,rmfunit,readwrite,blocksize
  character (len=500) comment
  integer Elocol,Ehicol,k
  logical exact,anynull

! Open rmf file with read-only access
  status = 0
  call ftgiou(rmfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(rmfunit,rmffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open rmf file'

! Shift to  extension "EBOUNDS"
  status = 0
  call ftmnhd(rmfunit,2,'EBOUNDS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EBOUNDS'

! Work out the correct columns
  call ftgcno(rmfunit,exact,'E_MIN',Elocol,status)
  call ftgcno(rmfunit,exact,'E_MAX',Ehicol,status)

! Read in ECHN()
  do k = 1,NCHN
     call ftgcve(rmfunit,Elocol,k,1,1,-1.0,ECHN(k-1),anynull,status)
     call ftgcve(rmfunit,Ehicol,k,1,1,-1.0,ECHN(k),anynull,status)
  end do

! Close mrf file and free up unit number
  call ftclos(rmfunit,status)
  call ftfiou(rmfunit,status)
  
  return
end subroutine getECHN
!-----------------------------------------------------------------------




  
!-----------------------------------------------------------------------
subroutine getmuE(mrffile,arffile,nenerg,E,muE)
! Read modulation factor vs E from mrf and arf files
  implicit none
  integer nenerg
  character (len=500) mrffile,arffile
  real muE(nenerg),E(0:nenerg)
  integer status,mrfunit,readwrite,blocksize,arfunit
  character (len=500) comment
  integer Elocol,Ehicol,respcol,k
  real Aeff(nenerg)
  logical exact,anynull
  
! Open mrf file with read-only access
  status = 0
  call ftgiou(mrfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(mrfunit,mrffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open mrf file'

! Shift to  extension "SPECRESP"
  status = 0
  call ftmnhd(mrfunit,2,'SPECRESP',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECRESP'

! Work out the correct columns
  call ftgcno(mrfunit,exact,'ENERG_LO',Elocol,status)
  call ftgcno(mrfunit,exact,'ENERG_HI',Ehicol,status)
  call ftgcno(mrfunit,exact,'SPECRESP',respcol,status)

! Read in mu(E)*Aeff(E)
  do k = 1,NENERG
     call ftgcve(mrfunit,Elocol,k,1,1,-1.0,E(k-1),anynull,status)
     call ftgcve(mrfunit,Ehicol,k,1,1,-1.0,E(k),anynull,status)
     call ftgcve(mrfunit,respcol,k,1,1,-1.0,muE(k),anynull,status)
  end do

! Close mrf file and free up unit number
  call ftclos(mrfunit,status)
  call ftfiou(mrfunit,status)
  
! Open arf file with read-only access
  status = 0
  call ftgiou(arfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(arfunit,arffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open mrf file'

! Shift to  extension "SPECRESP"
  status = 0
  call ftmnhd(arfunit,2,'SPECRESP',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECRESP'

! Work out the correct column
  call ftgcno(arfunit,exact,'SPECRESP',respcol,status)

! Read in Aeff(E)
  do k = 1,NENERG
     call ftgcve(arfunit,respcol,k,1,1,-1.0,Aeff(k),anynull,status)
     muE(k) = muE(k) / Aeff(k)
  end do

! Close arf file and free up unit number
  call ftclos(arfunit,status)
  call ftfiou(arfunit,status)

  return
end subroutine getmuE  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getNENERG(mrffile,nenerg)
  implicit none
  character (len=500) mrffile
  integer nenerg
  integer status,mrfunit,readwrite,blocksize
  character (len=500) comment
  
! Open mrf file with read-only access
  status = 0
  call ftgiou(mrfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(mrfunit,mrffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open mrf file'

! Shift to  extension "SPECRESP"
  status = 0
  call ftmnhd(mrfunit,2,'SPECRESP',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECRESP'

! Get number of rows in the table
  call ftgkyj(mrfunit,'NAXIS2',NENERG,comment,status)
  if(status .ne. 0) stop 'Cannot determine No of rows'

! Close mrf file and free up unit number
  call ftclos(mrfunit,status)
  call ftfiou(mrfunit,status)

  return
end subroutine getNENERG
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getnevts(unit,nevts)
  implicit none
  integer unit,status,nevts
  character (len=500) comment
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'

! Get number of rows in the table: rows
  call ftgkyj(unit,'NAXIS2',nevts,comment,status)
  if(status .ne. 0) stop 'Cannot determine No of rows'
      
  return
end subroutine getnevts
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine getscal(unit,bscale,ascale)
! Read keywords AREASCAL and BACKSCAL
  implicit none
  integer unit,status
  real ascale,bscale
  character (len=500) comment
      
! Shift to  extension "SPECTRUM"
  status = 0
  call ftmnhd(unit,2,'SPECTRUM',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECTRUM'

! Get ascale
  call ftgkye(unit,'AREASCAL',ascale,comment,status)
  if(status .ne. 0) stop 'Cannot read AREASCAL'

! Get bscale
  call ftgkye(unit,'BACKSCAL',bscale,comment,status)
  if(status .ne. 0) stop 'Cannot read BACKSCAL'
  
  return
end subroutine getscal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getMJDref(unit,MJDref)
  implicit none
  integer unit,status
  double precision MJDref
  integer MJDrefI
  real MJDrefF
  character (len=500) comment
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'

! Get integer part of MJDref: MJDrefI
  call ftgkyj(unit,'MJDREFI',MJDrefI,comment,status)
  if(status .ne. 0) stop 'Cannot read MJDrefI'
  write(*,*)"MJDrefI=",MJDrefI

! Get fractional part of MJDref: MJDrefF
  call ftgkye(unit,'MJDREFF',MJDrefF,comment,status)
  if(status .ne. 0) stop 'Cannot read MJDrefF'
  write(*,*)"MJDrefI=",MJDrefF

! Add together
  MJDref = dble( MJDrefI ) + dble( MJDrefF )
      
  return
end subroutine getMJDref
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine get_tstart(unit,tstart,tstop,tlive,ton)
  implicit none
  integer unit,status
  double precision tstart,tstop,tlive,ton
  character (len=500) comment
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'

! Get tstart
  call ftgkyd(unit,'TSTART',tstart,comment,status)
  if(status .ne. 0) stop 'Cannot determine TSTART'
  !write(*,*)"tstart=",tstart

! Get tstop
  call ftgkyd(unit,'TSTOP',tstop,comment,status)
  if(status .ne. 0) stop 'Cannot determine TSTOP'
  !write(*,*)"tstop=",tstop

! Get Livetime
  call ftgkyd(unit,'LIVETIME',tlive,comment,status)
  if(status .ne. 0) stop 'Cannot determine LIVETIME'
  !write(*,*)"tlive=",tlive

! Get Ontime
  call ftgkyd(unit,'ONTIME',ton,comment,status)
  if(status .ne. 0) stop 'Cannot determine ONTIME'
  !write(*,*)"ton=",ton
  
  return
end subroutine get_tstart
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine getspec(unit,nevts,nchn,x0,y0,r0,Ispec)
! Ispec is NOT initialised to allow me to add DUs
  implicit none
  integer unit,nevts,nchn
  real x0,y0,r0,Ispec(nchn)
  integer status,xcol,ycol,picol,k,PIchn,i
  logical exact,anynull
  character (len=50) EXNAME,comment
  real x,y
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'X',xcol,status)
  call ftgcno(unit,exact,'Y',ycol,status)
  call ftgcno(unit,exact,'PI',picol,status)

! Extract spectrum
  do k = 1,nevts
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     i = PIchn + 1
     call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
     call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
     if( (x-x0)**2 + (y-y0)**2 .le. r0**2 )then
        Ispec(i) = Ispec(i) + 1.0
     end if
  end do

  return
end subroutine getspec
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine getcounts(unit,nevts,reg,x0,y0,rin,rout,counts)
! counts is NOT initialised to allow me to add DUs
! Extraction region:
! reg = 1: circle of radius rin (rout redundant)
! reg = 2: annulus inner and outer radii rin and rout
  implicit none
  integer unit,nevts,reg,counts
  real x0,y0,rin,rout
  integer status,xcol,ycol,k
  logical exact,anynull
  character (len=50) EXNAME,comment
  real x,y
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'X',xcol,status)
  call ftgcno(unit,exact,'Y',ycol,status)

! Extract events
  if( reg .eq. 1 )then
     !Circular region
     do k = 1,nevts
        call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
        call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
        if( (x-x0)**2 + (y-y0)**2 .le. rin**2 )then
           counts = counts + 1
        end if
     end do
  else if( reg .eq. 2 )then
     !Annulus
     do k = 1,nevts
        call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
        call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
        if( (x-x0)**2 + (y-y0)**2 .ge. rin**2 )then
           if( (x-x0)**2 + (y-y0)**2 .le. rout**2 )then
              counts = counts + 1
           end if
        end if
     end do
  else
     write(*,*)"Bad reg option in getcounts subroutine"
  end if
     
  return
end subroutine getcounts
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine readpha(unit,nchn,Ispec,dIspec)
! Reads spectrum from the pha file
  implicit none
  integer unit,nchn
  real Ispec(nchn),dIspec(nchn)
  integer status,rcol,drcol,k
  logical exact,anynull
  character (len=50) EXNAME,comment
  real x,y
  
! Shift to  extension "SPECTRUM"
  status = 0
  call ftmnhd(unit,2,'SPECTRUM',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECTRUM'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'RATE',rcol,status)
  call ftgcno(unit,exact,'STAT_ERR',drcol,status)

! Extract spectrum
  do k = 1,nchn
     call ftgcve(unit,rcol ,k,1,1,-1.0,Ispec(k) ,anynull,status)
     call ftgcve(unit,drcol,k,1,1,-1.0,dIspec(k),anynull,status)
  end do

  return
end subroutine readpha
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getngti(unit,ngti)
  implicit none
  integer unit,status,ngti
  character (len=500) comment
      
! Shift to  extension "GTI"
  status = 0
  call ftmnhd(unit,2,'GTI',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension GTI'

! Get number of rows in the table: rows
  call ftgkyj(unit,'NAXIS2',ngti,comment,status)
  if(status .ne. 0) stop 'Cannot determine No of rows'
  !write(*,*)"No of GTIs=",ngti
      
  return
end subroutine getngti
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine getgti(unit,ngti,gtistart,gtiend)
!
  implicit none
  integer unit,ngti
  double precision gtistart(ngti),gtiend(ngti)
  integer status,bcol,ecol,k
  logical exact,anynull
      
! Shift to  extension "GTI"
  status = 0
  call ftmnhd(unit,2,'GTI',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension GTI'
      
! Get the column number for required inputs
  exact=.false.
  status = 0
  call ftgcno(unit,exact,'start',bcol,status)
  call ftgcno(unit,exact,'stop' ,ecol,status)
  if( status .ne. 0 ) stop 'cannot find column numbers for GTIs'
  
! Read in GTIs
  do k = 1,ngti
     call ftgcvd(unit,bcol,k,1,1,-1.0,gtistart(k),anynull,status)
     call ftgcvd(unit,ecol,k,1,1,-1.0,gtiend(k)  ,anynull,status)
  end do
  
  return
end subroutine getgti
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getlivetime(ngti,gtistart,gtiend,tstart,tstop,dt,nn,livetime)
  ! Calculates the amount of in-GTI time for each time bin
  implicit none
  integer ngti,nn
  double precision gtistart(ngti),gtiend(ngti),tstart,tstop
  double precision livetime(nn)
  real dt
  integer k,gti
  double precision tmin,tmax,tlo,thi

  do k = 1,nn
     tmin = tstart + (k-1)*dt
     tmax = tstart +  k   *dt
     livetime(k) = 0.0
     do gti = 1,ngti
        tlo = max( tmin , gtistart(gti) )
        thi = min( tmax , gtiend(gti)   )
        livetime(k) = livetime(k) + max( 0.d0 , thi-tlo )
     end do
  end do
  
  return
end subroutine getlivetime
!-----------------------------------------------------------------------












!-----------------------------------------------------------------------
subroutine simple_lightcurve(unit,nevts,t0,dt,won,pilo,pihi,nchn,muEi,reg,x0,y0,rin,rout,nn,lc,ilc,qlc,ulc,W2lc)
! Extraction region:
! reg = 1: circle radius rin (rout redundant)
! reg = 2: annulus between rin and rout
! Weights:
! won = 0: don't use weights
! won = 1: use weights
! Outputs:
! lc(1:nn)   = counts in each time bin
! ilc(1:nn)  = weighted counts in each time bin (same as lc for won=0)
! qlc(1:nn)  = Stokes Q of time bin in terms of (weighted) counts
! ulc(1:nn)  = Stokes U of time bin in terms of (weighted) counts
! W2lc(1:nn) = Sum of square of weightings for each time bin
! lc, ilc, qlc, ulc and W2 are NOT initialised so that I can add detectors together.
! Calculates light curve from event file already opened under unit "unit"
! Curve is extracted for the PI channels pilo->pihi (inclusive)
! Returns I, Q and U light curves in *counts*, not count rate.
  implicit none
  integer unit,nevts,pilo,pihi,nn,nchn,reg,won
  real x0,y0,rin,rout,lc(nn),ilc(nn),dt,qlc(nn),ulc(nn),muEi(nchn),W2lc(nn)
  double precision t0
  integer status,tcol,xcol,ycol,picol,k,PIchn,i
  integer qcol,ucol,wcol
  logical exact,anynull
  double precision time,tevt,q,u
  character (len=50) EXNAME,comment
  real x,y,weight
      
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'X',xcol,status)
  call ftgcno(unit,exact,'Y',ycol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'Q',qcol,status)
  call ftgcno(unit,exact,'U',ucol,status)
  call ftgcno(unit,exact,'W_MOM',wcol,status)

! Catch mistakes on the weight option
  won = min( won , 1 )
  won = max( won , 0 )
  
! Create light curve
  write(*,*)"Calculating light curve..."
  if( reg .eq. 1 )then
     !Circular region
     do k = 1,nevts
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
           call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)           
           if( (x-x0)**2 + (y-y0)**2 .le. rin**2 )then
              call ftgcve(unit,wcol,k,1,1,-1.0,weight,anynull,status)
              weight = weight * won + 1 - won
              call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
              tevt = time - t0
              i    = ceiling( tevt / dble(dt) )
              if( i .gt. 0 .and. i .le. nn )then
                 lc(i)   = lc(i)   + 1.0
                 ilc(i)  = ilc(i)  + weight
                 W2lc(i) = W2lc(i) + weight**2
                 call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
                 call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
                 qlc(i) = qlc(i) + weight * q / muEi(PIchn+1)
                 ulc(i) = ulc(i) + weight * u / muEi(PIchn+1)
              end if
           end if
        end if
     end do
  else if( reg .eq. 2 )then
     !Annular region
     do k = 1,nevts
        call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
        if( pichn .ge. pilo .and. pichn .le. pihi )then
           call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
           call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
           if( (x-x0)**2 + (y-y0)**2 .ge. rin**2 )then
              if( (x-x0)**2 + (y-y0)**2 .le. rout**2 )then
                 call ftgcve(unit,wcol,k,1,1,-1.0,weight,anynull,status)
                 weight = weight * won + 1 - won
                 call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
                 tevt = time - t0
                 i    = ceiling( tevt / dble(dt) )
                 if( i .gt. 0 .and. i .le. nn )then
                    lc(i)   = lc(i)   + 1.0
                    ilc(i)  = ilc(i)  + weight
                    W2lc(i) = W2lc(i) + weight**2
                    call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
                    call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
                    qlc(i) = qlc(i) + weight * q / muEi(PIchn+1)
                    ulc(i) = ulc(i) + weight * u / muEi(PIchn+1)
                 end if
              end if
           end if
        end if
     end do
  else
     write(*,*)"Bad reg option in simple_lightcurve"
  end if
  write(*,*)"...finished calculating light curve"
  
  return
end subroutine simple_lightcurve
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine getmueff(unit,nevts,won,pilo,pihi,nchn,muEi,x0,y0,r0,mueff)
  !Calculates effective modulation factor
  implicit none
  integer unit,nevts,won,pilo,pihi,nchn
  real x0,y0,r0,muEi(nchn),mueff
  integer status,xcol,ycol,picol,k,PIchn
  integer qcol,ucol,wcol
  logical exact,anynull
  double precision q,u
  character (len=50) EXNAME,comment
  real x,y,qsum,usum,qwsum,uwsum,weight
  
! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'X',xcol,status)
  call ftgcno(unit,exact,'Y',ycol,status)
  call ftgcno(unit,exact,'PI',picol,status)
  call ftgcno(unit,exact,'Q',qcol,status)
  call ftgcno(unit,exact,'U',ucol,status)
  call ftgcno(unit,exact,'W_MOM',wcol,status)
      
! Read in all q and u values
  qsum  = 0.0
  usum  = 0.0
  qwsum = 0.0
  uwsum = 0.0
  do k = 1,nevts
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     if( pichn .ge. pilo .and. pichn .le. pihi )then
        call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
        call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
        if( (x-x0)**2 + (y-y0)**2 .le. r0**2 )then
           call ftgcve(unit,wcol,k,1,1,-1.0,weight,anynull,status)
           weight = weight * won + 1 - won
           call ftgcvd(unit,qcol,k,1,1,-1.0,q,anynull,status)
           call ftgcvd(unit,ucol,k,1,1,-1.0,u,anynull,status)
           !Add to sums
           qsum = qsum + weight * q
           usum = usum + weight * u
           qwsum = qwsum + weight * q / muEi(PIchn+1)
           uwsum = uwsum + weight * u / muEi(PIchn+1)
        end if
     end if
  end do

! Finally calculate mueff
  mueff = sqrt( ( qsum**2 + usum**2 ) / ( qwsum**2 + uwsum**2 ) )

  write(*,*)"mueff=",mueff

  
  return
end subroutine getmueff  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine sim_simplc(unit,nevts,t0,dt,pilo,pihi,bfrac,p0,psi0,nchn,muEi,x0,y0,r0,nn,qlc,ulc)
! qlc and ulc are NOT initialised so that I can add detectors together.
! Reads in the actual event list but generates random modulation angle
! for each event. Eventually need to subtract spurious modulation too.
! Returns I, Q and U light curves in *counts*, not count rate.
  implicit none
  integer unit,nevts,pilo,pihi,nn,nchn
  real dt,p0,psi0,qlc(nn),ulc(nn),bfrac
  real x0,y0,r0,muEi(nchn)
  double precision t0
  integer status,tcol,xcol,ycol,picol,k,PIchn,i
  logical exact,anynull
  double precision time,tevt,q,u
  character (len=50) EXNAME,comment
  integer maxit,j,idum
  real xacc,En,mu,CDF,ran2,psi,mod_angle,x,y,eta,pi
  data idum/-389539/
  save idum
  pi = acos(-1.0)

  maxit = 100
  xacc  = 1e-5

! Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension EVENTS'
      
! Get the column number for required inputs
  exact=.false.
  call ftgcno(unit,exact,'TIME',tcol,status)
  call ftgcno(unit,exact,'X',xcol,status)
  call ftgcno(unit,exact,'Y',ycol,status)
  call ftgcno(unit,exact,'PI',picol,status)
      
! Create light curve
  !write(*,*)"Calclating light curve..."
  do k = 1,nevts
     call ftgcvj(unit,picol,k,1,1,-1.0,PIchn,anynull,status)
     if( pichn .ge. pilo .and. pichn .le. pihi )then
        call ftgcve(unit,xcol,k,1,1,-1.0,X,anynull,status)
        call ftgcve(unit,ycol,k,1,1,-1.0,Y,anynull,status)
        if( (x-x0)**2 + (y-y0)**2 .le. r0**2 )then
           call ftgcvd(unit,tcol,k,1,1,-1.0,time,anynull,status)
           tevt = time - t0
           i    = ceiling( tevt / dble(dt) )
           if( i .gt. 0 .and. i .le. nn )then
              !Select a random modulation angle
              mu  = muEi(PIchn+1)
              CDF = ran2(idum)
              eta = ran2(idum)
              if( eta .lt. bfrac )then
                 !Background event
                 psi = 2.0*pi * ran2(idum) - pi
              else
                 !Source event
                 psi = mod_angle(CDF,mu,p0,psi0,xacc,maxit)
              end if
              qlc(i) = qlc(i) + 2.0 * cos( 2.0*psi ) / mu
              ulc(i) = ulc(i) + 2.0 * sin( 2.0*psi ) / mu
           end if
        end if
     end if
  end do
  return
end subroutine sim_simplc
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
function mod_angle(CDF,mu,p0,psi0,xacc,maxit)
! Given CDF, a number between 0 and 1, the function returns the
! corresponding value of modulation angle, psi. An iteration scheme is
! used with a target accuracy xacc and maximum number of iterations
! maxit.
  implicit none
  integer maxit,j
  real mod_angle,CDF,mu,p0,psi0,xacc
  real psi,psiprev,pi
  pi = acos(-1.0)
  !Initial guess for psi-psi0
  psi = 2.0 * pi * CDF
  psiprev = 2.0 * psi
  !Now iterate for psi-psi0
  j = 1
  do while( abs(psi-psiprev) .gt. xacc .and. j .lt. maxit )
     psiprev = psi
     psi = 2.0*pi*CDF - 0.5*mu*p0 * sin( 2.0 * psi )
     j = j + 1
  end do
  mod_angle = psi + psi0
  if( abs(psi-psiprev) .gt. xacc )then
     write(*,*)"Warning! mod_angle did not converge"
  end if
  return
end function mod_angle
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&     
       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
       IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
! Long period (> 2 × 10^{18}) random number generator of L’Ecuyer with Bays-Durham shuffle
! and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
! of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
! alter idum between successive deviates in a sequence. RNMX should approximate the largest
! floating value that is less than 1.
  INTEGER idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then !Initialize.
     idum=max(-idum,1) !Be sure to prevent idum = 0.
     idum2=idum
     do j=NTAB+8,1,-1 !Load the shuffle table (after 8 warm-ups).
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     end do
     iy=iv(1)
  end if
  k=idum/IQ1 !Start here when not initializing.
  idum=IA1*(idum-k*IQ1)-k*IR1     !Compute idum=mod(IA1*idum,IM1) without over-
  if (idum.lt.0) idum=idum+IM1    !flows by Schrage’s method.
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2 !Compute idum2=mod(IA2*idum2,IM2) likewise.
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV !Will be in the range 1:NTAB.
  iy=iv(j)-idum2 !Here idum is shuffled, idum and idum2 are combined to generate output.
  iv(j)=idum 
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX) !Because users don’t expect endpoint values.
  return  
END FUNCTION ran2
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
FUNCTION gasdev(idum)
  implicit none
  INTEGER idum
  REAL gasdev
!    USES ran2
  INTEGER iset
  REAL fac,gset,rsq,v1,v2,ran2
  SAVE iset,gset
  DATA iset/0/
  if (iset.eq.0) then
1    v1=2.*ran2(idum)-1.
     v2=2.*ran2(idum)-1.
     rsq=v1**2+v2**2
     if(rsq.ge.1..or.rsq.eq.0.)goto 1
     fac=sqrt(-2.*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
  else
     gasdev=gset
     iset=0
  end if
  return
END FUNCTION gasdev
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      subroutine deletefile(filename,status)

!C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!C  Simply return if status is greater than zero
      if (status .gt. 0)return

!C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
!C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!C  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mkfourchar(i,A)
! c
! c input the integer i and this outputs a character A which has
! c a length of 4 with leading zeros if required
! c
  integer i
  character (LEN=4) A
  if( i .lt. 10 )then
     write(A,'(A3,I1)')'000',i
  else if( i .lt. 100 )then
     write(A,'(A2,I2)')'00',i
  else if( i .lt. 1000 )then
     write(A,'(A1,I3)')'0',i
  else if( i .lt. 10000 )then
     write(A,'(I4)')i
  else
     write(*,*)"Integer too long in mkfourchar!"
   end if
   return
 end subroutine mkfourchar
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine mkthreechar(i,A)
! c
! c input the integer i and this outputs a character A which has
! c a length of 3 with leading zeros if required
! c
  integer i
  character (LEN=3) A
  if( i .lt. 10 )then
     write(A,'(A2,I1)')'00',i
  else if( i .lt. 100 )then
     write(A,'(A1,I2)')'0',i
  else if( i .lt. 1000 )then
     write(A,'(I3)')i
  else
     write(*,*)"Integer too long in mkthreechar!"
  end if
  return
end subroutine mkthreechar
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
function lore_f(f,param)
! Function version of my XSPEC model mylore
  implicit none
  real lore_f,f,param(3),f0,HWHM,norm,pi
  pi  = acos(-1.0)
  !Parameters
  f0   = param(1)
  HWHM = param(2)
  norm = param(3)
  !Function
  lore_f = HWHM / ( (f-f0)**2. + HWHM**2. )
  lore_f = lore_f / ( pi/2. + atan(f0/HWHM) )
  lore_f = 2.0 * norm * lore_f
  return
end function lore_f
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
function pndt_f(f,param)
! Function version of my XSPEC model for dead time affected Poisson noise.
! In fractional rms normalisation.
  implicit none
  real pndt_f,f,param(3)
  real tau,rate,norm,wiggles
  
! Parameters
  tau  = param(1)  !Dead time (sec)
  rate = param(2)  !Total intrinsic count rate per independent detector
  norm = param(3)  !Poisson noise in the absence of dead time
  
!Function
  pndt_f = norm * wiggles(f,tau,rate)

  return
end function pndt_f
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine getAEmuE(mrffile,arffile,nenerg,E,AEmuE,Aeff)
! Read modulation factor vs E from mrf and arf files
  implicit none
  integer nenerg
  character (len=500) mrffile,arffile
  real AEmuE(nenerg),E(0:nenerg)
  integer status,mrfunit,readwrite,blocksize,arfunit
  character (len=500) comment
  integer Elocol,Ehicol,respcol,k
  real Aeff(nenerg)
  logical exact,anynull
  
! Open mrf file with read-only access
  status = 0
  call ftgiou(mrfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(mrfunit,mrffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open mrf file'

! Shift to  extension "SPECRESP"
  status = 0
  call ftmnhd(mrfunit,2,'SPECRESP',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECRESP'

! Work out the correct columns
  call ftgcno(mrfunit,exact,'ENERG_LO',Elocol,status)
  call ftgcno(mrfunit,exact,'ENERG_HI',Ehicol,status)
  call ftgcno(mrfunit,exact,'SPECRESP',respcol,status)

! Read in mu(E)*Aeff(E)
  do k = 1,NENERG
     call ftgcve(mrfunit,Elocol,k,1,1,-1.0,E(k-1),anynull,status)
     call ftgcve(mrfunit,Ehicol,k,1,1,-1.0,E(k),anynull,status)
     call ftgcve(mrfunit,respcol,k,1,1,-1.0,AEmuE(k),anynull,status)
  end do

! Close mrf file and free up unit number
  call ftclos(mrfunit,status)
  call ftfiou(mrfunit,status)
  
! Open arf file with read-only access
  status = 0
  call ftgiou(arfunit,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  readwrite = 0
  call ftopen(arfunit,arffile,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open mrf file'

! Shift to  extension "SPECRESP"
  status = 0
  call ftmnhd(arfunit,2,'SPECRESP',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension SPECRESP'

! Work out the correct column
  call ftgcno(arfunit,exact,'SPECRESP',respcol,status)

! Read in Aeff(E)
  do k = 1,NENERG
     call ftgcve(arfunit,respcol,k,1,1,-1.0,Aeff(k),anynull,status)
  end do

! Close arf file and free up unit number
  call ftclos(arfunit,status)
  call ftfiou(arfunit,status)

  return
end subroutine getAEmuE  
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine optimal(nseg,df,Nlore,lore_par,Nwn,wn_par,Nqpo,fmin,fmax,coeff)
! Calculates the optimal filter given a model
  implicit none
  integer nseg,Nlore,Nwn,Nqpo
  real df,lore_par(3,Nlore),wn_par(3,Nwn),fmin,fmax
  real coeff(nseg/2)
  integer i,j
  real f,Pmodel,Pqpo,lore_f,pndt_f

  do i = 1,nseg/2
     f = i * df
     Pmodel = 0.0
     do j = 1,Nlore
        Pmodel   = Pmodel + lore_f(f,lore_par(:,j))
     end do
     do j = 1,Nwn
        Pmodel   = Pmodel + pndt_f(f,wn_par(:,j))
     end do
     Pqpo = lore_f(f,lore_par(:,Nqpo))
     coeff(i) = Pqpo / Pmodel
     if( f .gt. fmax ) coeff(i) = 0.0
     if( f .lt. fmin ) coeff(i) = 0.0
  end do

  return
end subroutine optimal  
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine FFTfilter(nseg,qseg,coeff,filt)
  implicit none
  integer nseg,j,j1,j2
  real qseg(nseg),filt(nseg),data(2*nseg),coeff(nseg/2)
! Load data array to be transformed
  do j = 1,nseg
     data(2*j-1) = qseg(j)
     data(2*j)   = 0.0
  end do
! Transform
  call ourfour1(data,nseg,1)
! Set DC components to zero
  data(1) = 0.0
  data(2) = 0.0
! Now all AC components except for Nyquist
  do j = 1,nseg/2-1
     j1 = 2*j + 1
     data(j1)   = coeff(j) * data(j1)    !REAL
     data(j1+1) = coeff(j) * data(j1+1)  !IMAGINARY
     j2 = 2*nseg + 2 - j1
     data(j2)   = data(j1)
     data(j2+1) = -data(j1+1)
  end do
! Finally Nyquist
  data(nseg+1) = coeff(nseg/2) * data(nseg+1) !Nyquist (real)
  data(nseg+2) = 0.0                          !Nyquist is always real
! Now inverse transform to get QPO light curve
  call ourfour1(data,nseg,-1)
  do j = 1,nseg
     filt(j) = data(2*j-1) / float(nseg)
  end do
  return
end subroutine FFTfilter
!-----------------------------------------------------------------------
      


!-----------------------------------------------------------------------
subroutine doHilb(nn,x,y)
  implicit none
! Apply the Hilbert transform to x(1:nn) and output it as y(1:nn)
  integer nn,j
  real x(nn),y(nn)
  complex Xf(nn/2),Yf(nn/2)      
! Now take the Fourier transform of x(k)
  call myrealFFT(nn,x,Xf)
! Do the Hilbert transform
  do j = 1,nn/2-1
     Yf(j) = Xf(j) * complex( 0. , 1. )
  end do
  Yf(nn/2) = Xf(nn/2)
! Now inverse Fourier transform, enforcing the reality condition
  call myrealinvFFT(nn,Xf,x)
  call myrealinvFFT(nn,Yf,y)
  return
end subroutine doHilb
!-----------------------------------------------------------------------

      

!------------------------------------------------------------------------
subroutine myrealFFT(n,at,Af)
  implicit none
  integer n,j
  real at(n),data(2*n)
  complex Af(n/2)
  do j = 1,n
     data(2*j-1) = at(j)
     data(2*j)   = 0.0
  end do
  call ourfour1(data,n,1)
  do j = 1, n/2
     Af(j)  = complex( data(2*j+1) , data(2*j+2) )
  end do
  return
end subroutine myrealFFT
!------------------------------------------------------------------------


!------------------------------------------------------------------------
subroutine myrealinvFFT(n,Af,at)
  implicit none
  integer n,j,i
  real at(n),data(2*n)
  complex Af(n/2)      
! +ve frequencies (including Nyquist)
  do j = 1, n/2
     data(2*j+1) =  real( Af(j) )
     data(2*j+2) = aimag( Af(j) )
  end do      
! DC component
  data(1) = 0.0
  data(2) = 0.0
! -ve frequencies
  do j = 1, n/2
     data(2*n-2*j+1) =  data(2*j+1) !Real
     data(2*n-2*j+2) = -data(2*j+2) !Imaginary
  end do      
! Inverst FFT
  call ourfour1(data,n,-1)
! Transfer to array at(n)
  do i = 1,n
     at(i) = data(2*i-1) / float(n)
  end do
  return
end subroutine myrealinvFFT
!------------------------------------------------------------------------




!-----------------------------------------------------------------------
function X2const(n,x,dx,xm)
! Calculates chisquared of data set x(1:n), dx(1:n) wrt constant
! model m.
  implicit none
  integer n
  real X2const,x(n),dx(n),xm
  integer j
  X2const = 0.0
  do j = 1,n
     X2const = X2const + ( x(j) - xm )**2 / dx(j)**2
  end do
  return
end function X2const  
!-----------------------------------------------------------------------
     
     
!-----------------------------------------------------------------------
function wav(n,x,dx)
! Returns error-weighted averaged of x(1:n)
  implicit none
  integer n
  real wav,x(n),dx(n)
  integer j
  real num,den
  num = 0.0
  den = 0.0
  do j = 1,n
     num = num + x(j) / dx(j)**2
     den = den + 1.0 / dx(j)**2
  end do
  wav = num / den
  return
end function wav
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine sterr_calc(nit,phibins,rms1,rms_av,drms_av)
  !Calculates standard error on the mean
  implicit none
  integer nit,phibins
  real rms1(phibins,nit),rms_av(phibins),drms_av(phibins)
  integer it,j
  drms_av = 0.0
  do it = 1,nit
     do j = 1,phibins
        drms_av(j) = drms_av(j) + ( rms1(j,it) - rms_av(j) )**2
     end do
  end do
  drms_av = sqrt(drms_av) / real(nit)
  return
end subroutine sterr_calc
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine meancalc(nit,phibins,rms1,rms_av)
  implicit none
  integer nit,phibins
  real rms1(phibins,nit),rms_av(phibins)
  integer it,j
  rms_av = 0.0
  do it = 1,nit
     do j = 1,phibins
        rms_av(j) = rms_av(j) + rms1(j,it)
     end do
  end do
  rms_av = rms_av / real(nit)
  return
end subroutine meancalc
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
function X2null(n,y,dy,yav)
! Calculates chisquared
  implicit none
  integer n,i
  real X2null,y(n),dy(n),yav(n)
  X2null = 0.0
  do i = 1,n
     X2null = X2null + ( y(i) - yav(i) )**2 / dy(i)**2
  end do
  return
end function X2null
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine fitplotsin(n,phi,y,dy,y0,unit)
!
! Fit data array y(n), dy(n) with a 180 degree period sine wave model
! and plot result in open unit number `unit'. phi(n) is in radians.
!
  implicit none
  integer n,unit
  real phi(n),y(n),dy(n),y0
  integer ma,mmax
  parameter (mmax=20)
  real a(mmax),chisq,da(mmax),pi,phideg,dphideg,ymod,dyda(mmax)
  integer ia(mmax),i
  external twosinfunc
  pi = acos(-1.0)
  
! Fit the model
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.0      !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  ia(4) = 0
  ia(5) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq,da)

! Write the data and model
  dphideg = ( phi(2)-phi(1) ) * 180.0 / pi
  do i = 1,n
     !Evaluate model
     call twosinfunc(phi(i),a,ymod,dyda,ma)
     !Write out model
     phideg = phi(i) * 180.0/pi
     write(unit,*)phideg,0.5*dphideg,y(i),dy(i),ymod,0.0
  end do
  write(unit,*)"no no"

  write(*,*)"chisq/dof = ",chisq,"/",real(n-3)
  write(*,*)"Best fitting parameters:"
  do i = 1,ma
    write(*,*)i,": ",a(i),"±",da(i)
  end do

  return
end subroutine fitplotsin  
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
subroutine sinebattle(n,phi,y,dy,y0,chisqww,dofww)
!
! Fit data with three models and return the chisquared and dof of the best one
! The models are:
! Model 1 = cons + 180 degree sine wave
! Model 2 = cons + 90 degree sine wave
! Model 3 = cons + 180 degree sine wave + 90 degree sine wave
!
! Inputs:
! n          number of data points
! phi(n)     array of x data points
! y(n)       array of y data points
! dy(n)      uncertainty on y
! y0         initial guess of constant y
! Outputs:
! chisqww    chisquared of winning model
! dofww      dof of winning model
!
  implicit none
  integer n
  real phi(n),y(n),dy(n),y0,chisqww,dofww
  integer ma,mmax
  parameter (mmax=20)
  real a(mmax),chisq1,da(mmax),chisq2,chisq3,chisqw,dof3,dofw
  integer ia(mmax)
  real Fstat,myftest,prob
  external twosinfunc

! Model 1: 180 degree sine wave
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.0      !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  ia(4) = 0
  ia(5) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq1,da)
! Model 2: 90 degree sine wave
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0    !Average
  a(2) = 0.0      !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.01     !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  ia(2) = 0
  ia(3) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq2,da)
! Model 3: two sine waves
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.01     !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq3,da)
! Model 1 vs Model 2
  chisqw = min( chisq1 , chisq2 )
  !Winner of that vs Model 3
  dof3  = real( n - 5 )
  dofw  = real( n - 3 )
  Fstat = myftest(chisq3,dof3,chisqw,dofw,prob)
  if( prob .lt. 0.1 )then !Accept double sine model if it is better with p<10% (2 sigma = 0.0455)
     chisqww = chisq3
     dofww   = dof3
  else
     chisqww = chisqw
     dofww   = dofw
  end if
  return
end subroutine sinebattle  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine sinebattle180(n,phi,y,dy,y0,chisqww,dofww,chisq1,dof1)
!
! Fit data with three models and return the chisquared and dof of the best one
! and the chisquared and dof of the 180 degrees model
! The models are:
! Model 1 = cons + 180 degree sine wave
! Model 2 = cons + 90 degree sine wave
! Model 3 = cons + 180 degree sine wave + 90 degree sine wave
!
! Inputs:
! n          number of data points
! phi(n)     array of x data points
! y(n)       array of y data points
! dy(n)      uncertainty on y
! y0         initial guess of constant y
! Outputs:
! chisqww    chisquared of winning model
! dofww      dof of winning model
! chisq1     chisquared of 180 degrees model
! dof1       dof of 180 degrees model
!
  implicit none
  integer n
  real phi(n),y(n),dy(n),y0,chisqww,dofww,chisq1,dof1
  integer ma,mmax
  parameter (mmax=20)
  real a(mmax),da(mmax),chisq2,chisq3,chisqw,dof3,dofw
  integer ia(mmax)
  real Fstat,myftest,prob
  external twosinfunc

! Model 1: 180 degree sine wave
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.0      !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  ia(4) = 0
  ia(5) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq1,da)
  dof1  = real( n - 3 )
! Model 2: 90 degree sine wave
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0    !Average
  a(2) = 0.0      !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.01     !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  ia(2) = 0
  ia(3) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq2,da)
! Model 3: two sine waves
  !Parameters
  ma   = 5        !Number of parameters
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 180 deg sine wave
  a(3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4) = 0.01     !Amplitude of 90 deg sine wave
  a(5) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia = 1
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,ma,twosinfunc,chisq3,da)
! Model 1 vs Model 2
  chisqw = min( chisq1 , chisq2 )
  !Winner of that vs Model 3
  dof3  = real( n - 5 )
  dofw  = real( n - 3 )
  Fstat = myftest(chisq3,dof3,chisqw,dofw,prob)
  if( prob .lt. 0.1 )then !Accept double sine model if it is better with p<10% (2 sigma = 0.0455)
     chisqww = chisq3
     dofww   = dof3
  else
     chisqww = chisqw
     dofww   = dofw
  end if
  return
end subroutine sinebattle180
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine fitfolded(n,phi,y,dy,y0,a,da,chisq,dof)
!
! Fits a model: constant + sin + sin.
! Designed to fit to QPO waveform.
!  
! Inputs:
! n          number of QPO phases
! phi(n)     array of x data points (i.e. QPO phase)
! y(n)       array of y data points (i.e. count rate / PD / PA)
! dy(n)      uncertainty on y
! y0         initial guess of constant y
! Outputs:
! a(1:5)     array of model parameters
! da(1:5)    array of uncertainties on model parameters
! chisq      chisquared of best-fitting model
! dof        dof of best-fitting model
  implicit none
! Input/Output
  integer n
  real phi(n),y(n),dy(n),y0,a(5),da(5),chisq,dof
! Internal
  integer ia(5)
  external wffunc  
  !Parameters
  a    = 0.0
  a(1) = y0       !Average
  a(2) = 0.01     !Amplitude of 1st harmonic
  a(3) = 0.0      !Phase of 1st harmonic (cycles)
  a(4) = 0.01     !Amplitude of 2nd harmonic
  a(5) = 0.0      !Phase of 2nd harmonic (cycles)
  !Free/fixed
  ia = 1
  !Fit the model
  call dofit(phi,y,dy,n,a,ia,5,wffunc,chisq,da)
  dof = real( n - 5 )
  return
end subroutine fitfolded
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine sinebattleplot(n,phi,y,dy,y0,unit)
!
! Fit data with three models and return the chisquared and dof of the best one
! The models are:
! Model 1 = cons + 180 degree sine wave
! Model 2 = cons + 90 degree sine wave
! Model 3 = cons + 180 degree sine wave + 90 degree sine wave
!
! Inputs:
! n          number of data points
! phi(n)     array of x data points
! y(n)       array of y data points
! dy(n)      uncertainty on y
! y0         initial guess of constant y
! Outputs:
! chisqww    chisquared of winning model
! dofww      dof of winning model
!
  implicit none
  integer n,unit,i
  real phi(n),y(n),dy(n),y0,chisqww,dofww
  integer ma,mmax
  parameter (mmax=20)
  real a(mmax,3),da(mmax,3)
  real chisq1,chisq2,chisq3,chisqw,dof3,dofw
  integer ia(mmax,3),winner
  real Fstat,myftest,prob,pi,phideg,dphideg
  real dyda(mmax),ymod
  external twosinfunc
  pi = acos(-1.0)

! Model 1: 180 degree sine wave
  !Parameters
  ma   = 5         !Number of parameters
  a(1,1) = y0       !Average
  a(2,1) = 0.01     !Amplitude of 180 deg sine wave
  a(3,1) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4,1) = 0.0      !Amplitude of 90 deg sine wave
  a(5,1) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia(:,1) = 1
  ia(4,1) = 0
  ia(5,1) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a(:,1),ia(:,1),ma,twosinfunc,chisq1,da(:,1))
! Model 2: 90 degree sine wave
  !Parameters
  ma   = 5         !Number of parameters
  a(1,2) = y0       !Average
  a(2,2) = 0.0      !Amplitude of 180 deg sine wave
  a(3,2) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4,2) = 0.01     !Amplitude of 90 deg sine wave
  a(5,2) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia(:,2) = 1
  ia(2,2) = 0
  ia(3,2) = 0
  !Fit the model
  call dofit(phi,y,dy,n,a(:,2),ia(:,2),ma,twosinfunc,chisq2,da(:,2))
! Model 3: two sine waves
  !Parameters
  ma   = 5        !Number of parameters
  a(1,3) = y0       !Average
  a(2,3) = 0.01     !Amplitude of 180 deg sine wave
  a(3,3) = 0.0      !Phase of 180 deg sine wave (radians)
  a(4,3) = 0.01     !Amplitude of 90 deg sine wave
  a(5,3) = 0.0      !Phase of 90 deg sine wave (radians)
  !Free/fixed
  ia(:,3) = 1
  !Fit the model
  call dofit(phi,y,dy,n,a(:,3),ia(:,3),ma,twosinfunc,chisq3,da(:,3))
! Model 1 vs Model 2
  winner = 1
  chisqw = chisq1
  if( chisq2 .lt. chisq1 )then
     winner = 2
     chisqw = chisq2
  end if
  !Winner of that vs Model 3
  dof3  = real( n - 5 )
  dofw  = real( n - 3 )
  Fstat = myftest(chisq3,dof3,chisqw,dofw,prob)
  if( prob .lt. 0.1 )then !Accept double sine model if it is better with p<10% (2 sigma = 0.0455)
     winner  = 3
     chisqww = chisq3
     dofww   = dof3
  else
     chisqww = chisqw
     dofww   = dofw
  end if

  write(*,*)"Winner is model ",winner
  write(*,*)"chisq/dof =",chisqww,"/",dofww

  write(*,*)"Best fitting parameters:"
  do i = 1,ma
    write(*,*)i,": ",a(i,winner),"±",da(i,winner)
  end do
  
! Write the data and model
  dphideg = ( phi(2)-phi(1) ) * 180.0 / pi
  do i = 1,n
     !Evaluate model
     call twosinfunc(phi(i),a(:,winner),ymod,dyda,ma)
     !Write out model
     phideg = phi(i) * 180.0/pi
     write(unit,*)phideg,0.5*dphideg,y(i),dy(i),ymod,0.0
  end do
  write(unit,*)"no no"
  
  
  return
end subroutine sinebattleplot
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine dofit_mlore(x,y,sig,ndata,a,ia,ma,chisq,da)
! Routine to fit model multilorefunc to data using mrqmin
! Inputs:
! ndata          = number of data points
! x(1:ndata)     = x data points
! y(1:ndata)     = y data points
! sig(1:ndata)   = error on y
! ma             = number of model parameters  
! a(1:ma)        = array of model parameters
! ia(1:ma)       = free (ai=1) or fixed (ai=0)
! Outputs:
! a(1:ma)        = array of best-fitting model parameters
! chisq          = minimum chisquared
! da(1:ma)       = error on parameters
!
  implicit none
  integer mmax,ndata,ma
  integer ia(ma),k,maxit
  parameter (mmax=20,maxit=100)
  real x(ndata),y(ndata),sig(ndata),a(ma),chisq,da(ma)
  external multilorefunc
  real alamda,covar(mmax,mmax),alpha(mmax,mmax)
  !Initialise
  alamda = -1.0
  call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,multilorefunc,alamda)
  !Run the fit loop
  k = 0
  do while( alamda .lt. 1e10 .and. k .lt. maxit )
     k = k + 1
     call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,multilorefunc,alamda)
  end do
  if( k .eq. maxit ) write(*,*)"Warning! Ran out of steps in fit!"
  !Calculate errors on best-fitting parameters
  alamda = 0.0
  call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,multilorefunc,alamda)
  do k = 1,ma
     da(k) = sqrt( abs(covar(k,k)) )
  end do
  return      
end subroutine dofit_mlore
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine dofit(x,y,sig,ndata,a,ia,ma,funcs,chisq,da)
! Routine to fit model defined by function funcs to data using mrqmin
! Inputs:
! ndata          = number of data points
! x(1:ndata)     = x data points
! y(1:ndata)     = y data points
! sig(1:ndata)   = error on y
! ma             = number of model parameters  
! a(1:ma)        = array of model parameters
! ia(1:ma)       = free (ai=1) or fixed (ai=0)
! funcs          = name of function that defines the model
! Outputs:
! a(1:ma)        = array of best-fitting model parameters
! chisq          = minimum chisquared
! da(1:ma)       = error on parameters
!
  implicit none
  integer mmax,ndata,ma
  integer ia(ma),k,maxit
  parameter (mmax=20,maxit=100)
  real x(ndata),y(ndata),sig(ndata),a(ma),chisq,da(ma)
  external funcs
  real alamda,covar(mmax,mmax),alpha(mmax,mmax)
  !Initialise
  alamda = -1.0
  call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,funcs,alamda)
  !Run the fit loop
  k = 0
  do while( alamda .lt. 1e10 .and. k .lt. maxit )
     k = k + 1
     call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,funcs,alamda)
  end do
  if( k .eq. maxit ) write(*,*)"Warning! Ran out of steps in fit!"
  !Calculate errors on best-fitting parameters
  alamda = 0.0
  call mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,mmax,chisq,funcs,alamda)
  do k = 1,ma
     da(k) = sqrt( abs(covar(k,k)) )
  end do
  return      
end subroutine dofit
!-----------------------------------------------------------------------


      
!-----------------------------------------------------------------------
SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,funcs,alamda)
  INTEGER ma,nca,ndata,ia(ma),MMAX
  REAL alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),sig(ndata),x(ndata),y(ndata)
  external funcs
  PARAMETER (MMAX=20)
!CU    USES covsrt,gaussj,mrqcof
  INTEGER j,k,l,m,mfit
  REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
  SAVE ochisq,atry,beta,da,mfit
  if(alamda.lt.0.)then
     mfit=0
     do j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
     end do
     alamda=0.001
     call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
     ochisq=chisq
     do j=1,ma
        atry(j)=a(j)
     end do
  end if
  j=0
  do l=1,ma
     if(ia(l).ne.0) then
        j=j+1
        k=0
        do m=1,ma
           if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
           endif
        end do
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
     endif
  end do
  call gaussj(covar,mfit,nca,da,1,1)
  if(alamda.eq.0.)then
     call covsrt(covar,nca,ma,ia,mfit)
     return
  end if
  j=0
  do l=1,ma
      if(ia(l).ne.0) then
         j=j+1
         atry(l)=a(l)+da(j)
      endif
  end do
  call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
  if(chisq.lt.ochisq)then
     alamda=0.1*alamda
     ochisq=chisq
     j=0
     do l=1,ma
        if(ia(l).ne.0) then
           j=j+1
           k=0
           do m=1,ma
              if(ia(m).ne.0) then
                 k=k+1
                 alpha(j,k)=covar(j,k)
              end if
           end do
           beta(j)=da(j)
           a(l)=atry(l)
        end if
     end do
  else
     alamda=10.*alamda
     chisq=ochisq
  end if
  return
end SUBROUTINE mrqmin
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,funcs)
  INTEGER ma,nalp,ndata,ia(ma),MMAX
  REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),y(ndata)
  EXTERNAL funcs
  PARAMETER (MMAX=20)
  INTEGER mfit,i,j,k,l,m
  REAL dy,sig2i,wt,ymod,dyda(MMAX)
  mfit=0
  do j=1,ma
     if (ia(j).ne.0) mfit=mfit+1
  end do
  do j=1,mfit
     do k=1,j
        alpha(j,k)=0.
     end do
     beta(j)=0.
  end do
  chisq=0.
  do i=1,ndata
     call funcs(x(i),a,ymod,dyda,ma)
     sig2i=1./(sig(i)*sig(i))
     dy=y(i)-ymod
     j=0
     do l=1,ma
        if(ia(l).ne.0) then
           j=j+1
           wt=dyda(l)*sig2i
           k=0
           do m=1,l
              if(ia(m).ne.0) then
                 k=k+1
                 alpha(j,k)=alpha(j,k)+wt*dyda(m)
              end if
           end do
           beta(j)=beta(j)+dy*wt
        end if
     end do
     chisq=chisq+dy*dy*sig2i
  end do
  do j=2,mfit
     do k=1,j-1
        alpha(k,j)=alpha(j,k)
     end do
  end do
  return
END SUBROUTINE mrqcof
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
  INTEGER ma,mfit,npc,ia(ma)
  REAL covar(npc,npc)
  INTEGER i,j,k
  REAL swap
  do i=mfit+1,ma
     do j=1,i
        covar(i,j)=0.
        covar(j,i)=0.
     end do
  end do
  k=mfit
  do j=ma,1,-1
     if(ia(j).ne.0)then
        do i=1,ma
           swap=covar(i,k)
           covar(i,k)=covar(i,j)
           covar(i,j)=swap
        end do
        do i=1,ma
           swap=covar(k,i)
           covar(k,i)=covar(j,i)
           covar(j,i)=swap
        end do
        k=k-1
     end if
  end do
  return
END SUBROUTINE covsrt
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
SUBROUTINE gaussj(a,n,np,b,m,mp)
  INTEGER m,mp,n,np,NMAX
  REAL a(np,np),b(np,mp)
  PARAMETER (NMAX=50)
  INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  REAL big,dum,pivinv
  do j=1,n
     ipiv(j)=0
  end do
  do i=1,n
     big=0.
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if(ipiv(k).eq.0)then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 end if
              else if (ipiv(k).gt.1) then
                 write(*,*)'singular matrix in gaussj'
              end if
           end do
        end if
     end do
     ipiv(icol)=ipiv(icol)+1
     if(irow.ne.icol)then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        end do
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        end do
     end if
     indxr(i)=irow
     indxc(i)=icol
     if(a(icol,icol).eq.0.) write(*,*)'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     end do
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     end do
     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           end do
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           end do
        end if
     end do
  end do
  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        end do
     end if
  end do
  return
END SUBROUTINE gaussj
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine wffunc(x,a,ymod,dyda,ma)
!
! Inputs:
! x         = phase (cycles)
! a(1:ma)   = array of parameters
! ma        = number of parameters (5)
! Outputs:
! ymod      = model evaluation at x
! dyda(1:ma) = derivative wrt a
!  
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax),pi
  pi = 3.141592653589793
  ymod = a(1) + a(2) * sin( 2.0*pi*(x-a(3)) ) + a(4) * sin( 4.0*pi*(x-a(5)) )
  dyda(1) = 1.0
  dyda(2) = sin( 2.0*pi*(x-a(3)) )
  dyda(3) = -2.0*pi * a(2) * cos( 2.0*pi*(x-a(3)) )
  dyda(4) = sin( 4.0*pi*(x-a(5)) )
  dyda(5) = -4.0*pi * a(4) * cos( 4.0*pi*(x-a(5)) )
  return
end subroutine wffunc  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine dwffunc(x,a,ymod,dyda,ma)
!
! Inputs:
! x         = phase (cycles)
! a(1:ma)   = array of parameters
! ma        = number of parameters (5)
! Outputs:
! ymod      = model evaluation at x
! dyda(1:ma) = derivative wrt a
!  
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  double precision x,a(mmax),ymod,dyda(mmax),pi
  pi = 3.141592653589793
  ymod = a(1) + a(2) * sin( 2.0*pi*(x-a(3)) ) + a(4) * sin( 4.0*pi*(x-a(5)) )
  dyda(1) = 1.0
  dyda(2) = sin( 2.0*pi*(x-a(3)) )
  dyda(3) = -2.0*pi * a(2) * cos( 2.0*pi*(x-a(3)) )
  dyda(4) = sin( 4.0*pi*(x-a(5)) )
  dyda(5) = -4.0*pi * a(4) * cos( 4.0*pi*(x-a(5)) )
  return
end subroutine dwffunc  
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine IQUsinfunc(x,a,ymod,dyda,ma)
!
! Inputs:
! x          = phase angle (radians)
! a(1:ma)    = array of parameters
! ma         = number of parameters
! Outputs:
! ymod       = model evaluation at x
! dyda(1:ma) = derivative wrt a
!
! Parameters:
! a(1) = dphi/pi * I
! a(2) = dphi/pi * Q
! a(3) = dphi/pi * U
!
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax)
  ymod = a(1) + a(2) * cos(2.0*x) + a(3) * sin(2.0*x)
  dyda(1) = 1.0
  dyda(2) = cos(2.0*x)
  dyda(3) = sin(2.0*x)
  return
end subroutine IQUsinfunc 
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine NULLsinfunc(x,a,ymod,dyda,ma)
!
! Inputs:
! x          = phase angle (radians)
! a(1:ma)    = array of parameters
! ma         = number of parameters
! Outputs:
! ymod       = model evaluation at x
! dyda(1:ma) = derivative wrt a
!
! Parameters:
! a(1) = dphi/pi * P_I
! a(2) = q
! a(3) = u
!
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax)
  
  ymod = a(1) * ( 1.0 + a(2) * cos(2.0*x) + a(3) * sin(2.0*x) )
  dyda(1) = 1.0 + a(2) * cos(2.0*x) + a(3) * sin(2.0*x)
  dyda(2) = a(1) * cos(2.0*x)
  dyda(3) = a(1) * sin(2.0*x)
  return
end subroutine NULLsinfunc 
!-----------------------------------------------------------------------








!-----------------------------------------------------------------------
subroutine twosinfunc(x,a,ymod,dyda,ma)
!
! Inputs:
! x          = phase angle (radians)
! a(1:ma)    = array of parameters
! ma         = number of parameters
! Outputs:
! ymod       = model evaluation at x
! dyda(1:ma) = derivative wrt a
!
  implicit none
  integer ma,mmax
  parameter (mmax=20)
  real x,a(mmax),ymod,dyda(mmax)
  ymod = a(1) + a(2) * sin( 2.0*(x-a(3)) ) + a(4) * sin( 4.0*(x-a(5)) )
  dyda(1) = 1.0
  dyda(2) = sin( 2.0*(x-a(3)) )
  dyda(3) = -2.0 * a(2) * cos( 2.0*(x-a(3)) )
  dyda(4) = sin( 4.0*(x-a(5)) )
  dyda(5) = -4.0 * a(4) * cos( 4.0*(x-a(5)) )
  return
end subroutine twosinfunc
!-----------------------------------------------------------------------

  
!-----------------------------------------------------------------------
function myftest(X2,nu2,X1,nu1,p)
! c X2  = chisquared of unrestricted (more parameters) model
! c nu2 = d.o.f. of unrestricted (more parameters) model
! c X1  = chisquared of restricted (less parameters) model
! c nu1 = d.o.f. of restricted (less parameters) model
! c ftest = value of f-statistic
! c p     = p-value (probability that X2 is only better than X1 by chance)
  implicit none
  real X2,nu2,X1,nu1,myftest,p,F,a,b,x,betai
!c Calculate f-statistic
  F   = ( X1 - X2 ) / ( nu1 - nu2 )
  F   = F / ( X2 / nu2 )
  myftest = F
!c Calculate p-value
  a = 0.5 * nu2
  b = 0.5 * ( nu1 - nu2 )
  x = nu2 / ( nu2 + ( nu1 - nu2 )*F  )
  p = betai(a,b,x)
  return
end function myftest
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
FUNCTION betai(a,b,x)
  real betai,a,b,x
!CU    USES betacf,gammln
  real bt,betacf,gammln
  if(x.lt.0..or.x.gt.1.) write(*,*)'bad argument x in betai'
  if(x.eq.0..or.x.eq.1.)then
     bt=0.
  else
     bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
  end if
  if(x.lt.(a+1.)/(a+b+2.))then
     betai=bt*betacf(a,b,x)/a
     return
  else
     betai=1.-bt*betacf(b,a,1.-x)/b
     return
  end if
END FUNCTION betai
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
FUNCTION betacf(a,b,x)
  INTEGER MAXIT
  real betacf,a,b,x,EPS,FPMIN
  PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
  INTEGER m,m2
  real aa,c,d,del,h,qab,qam,qap
  qab=a+b
  qap=a+1.
  qam=a-1.
  c=1.
  d=1.-qab*x/qap
  if(abs(d).lt.FPMIN)d=FPMIN
  d=1./d
  h=d
  do m=1,MAXIT
     m2=2*m
     aa=m*(b-m)*x/((qam+m2)*(a+m2))
     d=1.+aa*d
     if(abs(d).lt.FPMIN)d=FPMIN
     c=1.+aa/c
     if(abs(c).lt.FPMIN)c=FPMIN
     d=1./d
     h=h*d*c
     aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
     d=1.+aa*d
     if(abs(d).lt.FPMIN)d=FPMIN
     c=1.+aa/c
     if(abs(c).lt.FPMIN)c=FPMIN
     d=1./d
     del=d*c
     h=h*del
     if(abs(del-1.).lt.EPS)goto 1
  end do
  write(*,*)'a or b too big, or MAXIT too small in betacf'
1 betacf=h
  return
END FUNCTION betacf
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
FUNCTION gammln(xx)
  REAL gammln,xx
  INTEGER j
  DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
  -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  end do
  gammln=tmp+log(stp*ser/x)
  return
END FUNCTION gammln
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
FUNCTION gammq(a,x)
! USES gcf,gser
  REAL a,gammq,x
  REAL gammcf,gamser,gln
  if(x.lt.0..or.a.le.0.)write(*,*) 'bad arguments in gammq'
  if(x.lt.a+1.)then
     call gser(gamser,a,x,gln)
     gammq=1.-gamser
  else
     call gcf(gammcf,a,x,gln)
     gammq=gammcf
  end if
  return
END FUNCTION gammq
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE gcf(gammcf,a,x,gln)
! USES gammln
  INTEGER ITMAX
  REAL a,gammcf,gln,x,EPS,FPMIN
  PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
  INTEGER i
  REAL an,b,c,d,del,h,gammln
  gln=gammln(a)
  b=x+1.-a
  c=1./FPMIN
  d=1./b
  h=d
  do i=1,ITMAX
     an=-i*(i-a)
     b=b+2.
     d=an*d+b
     if(abs(d).lt.FPMIN)d=FPMIN
     c=b+an/c
     if(abs(c).lt.FPMIN)c=FPMIN
     d=1./d
     del=d*c
     h=h*del
     if(abs(del-1.).lt.EPS)goto 1
  end do
  write(*,*) 'a too large, ITMAX too small in gcf'
1 gammcf=exp(-x+a*log(x)-gln)*h
      return
    END SUBROUTINE gcf
!-----------------------------------------------------------------------    




!-----------------------------------------------------------------------
SUBROUTINE gser(gamser,a,x,gln)
! USES gammln
  INTEGER ITMAX
  REAL a,gamser,gln,x,EPS
  PARAMETER (ITMAX=100,EPS=3.e-7)
  INTEGER n
  REAL ap,del,sum,gammln
  gln=gammln(a)
  if(x.le.0.)then
     if(x.lt.0.)write(*,*) 'x < 0 in gser'
     gamser=0.
     return
  end if
  ap=a
  sum=1./a
  del=sum
  do n=1,ITMAX
     ap=ap+1.
     del=del*x/ap
     sum=sum+del
     if(abs(del).lt.abs(sum)*EPS)goto 1
  end do
  write(*,*) 'a too large, ITMAX too small in gser'
1 gamser=sum*exp(-x+a*log(x)-gln)
  return
END SUBROUTINE gser
!-----------------------------------------------------------------------

