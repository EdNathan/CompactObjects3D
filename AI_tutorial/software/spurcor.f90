include 'a_modules.f90'
include 'a_subroutines.f90'

program spurcor
! Adds back spurious polarisation to get modulation angle wrt to north.
!
! 1) Read in level 1 and level 2 event files
! 2) Match them
! 3) Calculate spurious Stokes parameters for each one
! 4) Check that the result is sensible
! 5) Create new "level 1.5" event file
!
! Compile command
! gfortran -L$HEADAS/lib -lcfitsio spurcor.f90
!  
  implicit none
  character (len=500) lvl1_DU1_evt,lvl1_DU2_evt,lvl1_DU3_evt
  character (len=500) lvl2_DU1_evt,lvl2_DU2_evt,lvl2_DU3_evt
  character (len=500) sp_dir,sp_file1,sp_file2,sp_file3,dir
  integer status,readwrite,blocksize
  integer evt1unit1,evt1unit2,evt1unit3,evt2unit1,evt2unit2,evt2unit3
  integer nevts1_l1,nevts2_l1,nevts3_l1,nevts1_l2,nevts2_l2,nevts3_l2
  integer spunit1,spunit2,spunit3
  integer deltaunit
  character (len=500) outfile1,outfile2,outfile3
  double precision tol
  
! Input parameters -------------------------------------

  !Accuracy level for spurious polarisation calc
  tol = 1e-4

  !Directory including the event files
  dir = '../CygX1data/03010101/'
  
  !Level 1 event files
  lvl1_DU1_evt = '/event_l1/ixpe03010101_det1_evt1_v01.fits'
  lvl1_DU2_evt = '/event_l1/ixpe03010101_det2_evt1_v01.fits'
  lvl1_DU3_evt = '/event_l1/ixpe03010101_det3_evt1_v01.fits'
  
  !Level 2 event files
  lvl2_DU1_evt = '/event_l2/ixpe03010101_det1_evt2_v01.fits'
  lvl2_DU2_evt = '/event_l2/ixpe03010101_det2_evt2_v01.fits'
  lvl2_DU3_evt = '/event_l2/ixpe03010101_det3_evt2_v01.fits'
  
  !Spurious polarization maps
  sp_dir   = '/Users/administrator/caldb/data/ixpe/gpd/bcf/spmod/'
  sp_file1 = 'ixpe_d1_20170101_spmod_02.fits'
  sp_file2 = 'ixpe_d2_20170101_spmod_02.fits'
  sp_file3 = 'ixpe_d3_20170101_spmod_02.fits'
  
  !New level 1.5 event files to be created
  outfile1 = '/event_1pt5/ixpe03010101_det1_evt1pt5_v01.fits'
  outfile2 = '/event_1pt5/ixpe03010101_det2_evt1pt5_v01.fits'
  outfile3 = '/event_1pt5/ixpe03010101_det3_evt1pt5_v01.fits'
  
! ------------------------------------------------------

! Add full path to file names
  lvl1_DU1_evt = trim(dir) // '/' // trim(lvl1_DU1_evt)
  lvl1_DU2_evt = trim(dir) // '/' // trim(lvl1_DU2_evt)
  lvl1_DU3_evt = trim(dir) // '/' // trim(lvl1_DU3_evt)
  lvl2_DU1_evt = trim(dir) // '/' // trim(lvl2_DU1_evt)
  lvl2_DU2_evt = trim(dir) // '/' // trim(lvl2_DU2_evt)
  lvl2_DU3_evt = trim(dir) // '/' // trim(lvl2_DU3_evt)  
  sp_file1     = trim(sp_dir) // '/' // trim(sp_file1)
  sp_file2     = trim(sp_dir) // '/' // trim(sp_file2)
  sp_file3     = trim(sp_dir) // '/' // trim(sp_file3)
  outfile1     = trim(dir) // '/' // trim(outfile1)
  outfile2     = trim(dir) // '/' // trim(outfile2)
  outfile3     = trim(dir) // '/' // trim(outfile3)
  
! Open spurious modulation files with read-only access
  write(*,*)"------------------------------------------------"
  write(*,*)"Opening suprious polarization files:"
  status = 0
  call ftgiou(spunit1,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU1 spurious modulation file..."
  readwrite = 0
  call ftopen(spunit1,sp_file1,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU1 spurious modulation file",trim(sp_file1)
  status = 0
  call ftgiou(spunit2,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU2 spurious modulation file..."
  readwrite = 0
  call ftopen(spunit2,sp_file2,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU2 spurious modulation file",trim(sp_file2)
  status = 0
  call ftgiou(spunit3,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU3 spurious modulation file..."
  readwrite = 0
  call ftopen(spunit3,sp_file3,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU3 spurious modulation file",trim(sp_file3)
  write(*,*)"------------------------------------------------"

  
! Open event files with read-only access
  write(*,*)"------------------------------------------------"
  write(*,*)"Opening level 1 event files:"
  status = 0
  call ftgiou(evt1unit1,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU1 level 1 event file..."
  readwrite = 0
  call ftopen(evt1unit1,lvl1_DU1_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU1 level 1 event file:",trim(lvl1_DU1_evt)
  status = 0
  call ftgiou(evt1unit2,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU2 level 1 event file..."
  readwrite = 0
  call ftopen(evt1unit2,lvl1_DU2_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU2 level 1 event file:",trim(lvl1_DU2_evt)
  status = 0
  call ftgiou(evt1unit3,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU3 level 1 event file..."
  readwrite = 0
  call ftopen(evt1unit3,lvl1_DU3_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU3 level 1 event file:",trim(lvl1_DU3_evt)
  write(*,*)"------------------------------------------------"
  write(*,*)"------------------------------------------------"
  write(*,*)"Opening level 2 event files"
  status = 0
  call ftgiou(evt2unit1,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU2 level 2 event file..."
  readwrite = 0
  call ftopen(evt2unit1,lvl2_DU1_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU1 level 2 event file:",trim(lvl2_DU1_evt)
  status = 0
  call ftgiou(evt2unit2,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU2 level 2 event file..."
  readwrite = 0
  call ftopen(evt2unit2,lvl2_DU2_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU2 level 2 event file:",trim(lvl2_DU2_evt)
  status = 0
  call ftgiou(evt2unit3,status)
  if( status .ne. 0 ) stop 'cannot open unused unit number'
  write(*,*)"Opening DU3 level 2 event file..."
  readwrite = 0
  call ftopen(evt2unit3,lvl2_DU3_evt,readwrite,blocksize,status)
  if( status .ne. 0 ) stop 'cannot open event file'
  write(*,*)"...opened DU3 level 2 event file:",trim(lvl2_DU3_evt)
  write(*,*)"------------------------------------------------"
 
! Get number of events in the event files
  !Level 1
  call getnevts(evt1unit1,nevts1_l1)
  call getnevts(evt1unit2,nevts2_l1)
  call getnevts(evt1unit3,nevts3_l1)
  !Level 2
  call getnevts(evt2unit1,nevts1_l2)
  call getnevts(evt2unit2,nevts2_l2)
  call getnevts(evt2unit3,nevts3_l2)
  
! Calculate the raw sky frame modulation angle for each level 2 event
! and write to a new level 1.5 event file
  call ftgiou(deltaunit,status)
  open(deltaunit,file='delta.dat')
  write(deltaunit,*)"skip on"
  write(*,*)"------------------------------------------------"
  write(*,*)"Calculating phi for each level 2 event from each DU"
  write(*,*)"------------------------------------------------"
  write(*,*)"DU1:"
  call getphi(evt1unit1,evt2unit1,nevts1_l1,nevts1_l2,spunit1,outfile1,deltaunit,tol)
  write(*,*)"output in: ",trim(outfile1)
  write(deltaunit,*)"no no"
  write(*,*)"------------------------------------------------"
  write(*,*)"DU2:"
  call getphi(evt1unit2,evt2unit2,nevts2_l1,nevts2_l2,spunit2,outfile2,deltaunit,tol)
  write(*,*)"output in: ",trim(outfile2)
  write(deltaunit,*)"no no"
  write(*,*)"------------------------------------------------"
  write(*,*)"DU3:"
  call getphi(evt1unit3,evt2unit3,nevts3_l1,nevts3_l2,spunit3,outfile3,deltaunit,tol)
  write(*,*)"output in: ",trim(outfile2)
  write(*,*)"------------------------------------------------"
  
! Close units and free up unit numbers
  call ftclos(spunit1, status)
  call ftfiou(spunit1, status)
  call ftclos(spunit2, status)
  call ftfiou(spunit2, status)
  call ftclos(spunit3, status)
  call ftfiou(spunit3, status)
  call ftclos(evt1unit1, status)
  call ftfiou(evt1unit1, status)
  call ftclos(evt1unit2, status)
  call ftfiou(evt1unit2, status)
  call ftclos(evt1unit3, status)
  call ftfiou(evt1unit3, status)
  call ftclos(evt2unit1, status)
  call ftfiou(evt2unit1, status)
  call ftclos(evt2unit2, status)
  call ftfiou(evt2unit2, status)
  call ftclos(evt2unit3, status)
  call ftfiou(evt2unit3, status)
  
end program spurcor


!-----------------------------------------------------------------------
subroutine copyhdu(inunit,outfile,outunit)
! Inputs
! inunit:  unit number of (alreay opened) input file.
! outfile: name of output file to be created
! Outputs
! outunit: unit number of output file
  implicit none
  integer status,inunit,outunit,blocksize,morekeys,hdutype
  character (len=500) outfile
  
! The STATUS parameter must always be initialized.
  status=0

! Delete the file if it already exists, so we can then recreate it
! The deletefile subroutine is listed at the end of this file.
  call deletefile(outfile,status)

! Get  unused Logical Unit Numbers to use to open the FITS file.
  call ftgiou(outunit,status)

! Create the new empty FITS file (value of blocksize is ignored)
  blocksize=1
  call ftinit(outunit,outfile,blocksize,status)

! Skip to the 2nd extension in the input file
  call ftmahd(inunit,2,hdutype,status)
  
! FTCOPY copies the current HDU from the input FITS file to the output
! file.  The MOREKEY parameter allows one to reserve space for additional
! header keywords when the HDU is created.   FITSIO will automatically
! insert more header space if required, so programmers do not have to
! reserve space ahead of time, although it is more efficient to do so if
! it is known that more keywords will be appended to the header.
  morekeys=0
  call ftcopy(inunit,outunit,morekeys,status)

! Append/create a new empty extension on the end of the output file
  call ftcrhd(outunit,status)
      
! Skip to the 3rd extension in the input file.
  call ftmahd(inunit,3,hdutype,status)

!  FTCOPY now copies the binary table from the input FITS file
!  to the output file.
   call ftcopy(inunit,outunit,morekeys,status)  

 end subroutine copyhdu
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine getphi(unit_l1,unit_l2,nevts_l1,nevts_l2,spunit,outfile,deltaunit,tol)
  implicit none
  integer unit_l1,unit_l2,nevts_l1,nevts_l2,spunit
  character (len=500) outfile
  integer status,tcol,xcol,ycol,picol,qcol,ucol,sm_chns,deltaunit
  parameter (sm_chns=6)
  logical exact,anynull
  integer k,PIchn,k1,tcol1,phicol
  double precision phi1,q1,u1,q2,u2,phi,p2,qc,uc,pc,Delta
  double precision DETQ_SM(sm_chns),DETU_SM(sm_chns),qsp,usp,tol
  real X,Y
  double precision time_l2,t0,time,time_prev,pi
  integer DETX_PX,DETY_PX
  integer row,qsmcol,usmcol,chnscol
  integer xrow,yrow,values,bad,good,outunit
  integer DET_CHN_SM(sm_chns),quality,offgrid,j,ncol
  character (len=500) comment
  integer colnum,ncols
  character (len=16) extname,ttype(5),tform(5),tunit(5)
  double precision sum_q1_sky,sum_u1_sky,sum_q2_sky,sum_u2_sky
  double precision sum_qsp_sky,sum_usp_sky,qsp_sky,usp_sky
  double precision Delq,Delu
  pi = acos( -1.d0 )
  
! Copy the level 2 event file to a new level 1pt5 event file
  call copyhdu(unit_l2,outfile,outunit)
  
! Event files: Shift to  extension "EVENTS"
  status = 0
  call ftmnhd(unit_l2,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'Level 2: cannot shift to extension EVENTS'
  status = 0
  call ftmnhd(unit_l1,2,'EVENTS',0,status)
  if( status .ne. 0 ) stop 'Level 1: cannot shift to extension EVENTS'
  
! Get number of columns in the table
  call ftgkyj(unit_l2,'TFIELDS',ncol,comment,status)

! Add 5 new columns to the end of the EVENTS binary table
  !First move to the EVENTS extension
  status = 0
  call ftmnhd(outunit,2,'EVENTS',0,status)
  if( status .ne. 0 ) write(*,*)"Couldn't move to EVENTS extension of outfile"
  
  !Insert the names and types of the new columns.
  !I can't work out how to input the units of each column.
  status = 0
  colnum = ncol + 1
  ncols  = 5
  
  ttype(1) = 'PHI'
  ttype(2) = 'QUAL'
  ttype(3) = 'DELTA'
  ttype(4) = 'QSP'
  ttype(5) = 'USP'  
  tform(1) = '1D'
  tform(2) = '1J'
  tform(3) = '1D'
  tform(4) = '1D'
  tform(5) = '1D'  
  tunit(1) = 'rad'
  tunit(2) = ' '
  tunit(3) = 'rad'
  tunit(4) = ' '
  tunit(5) = ' '

  call FTICLS(outunit,colnum,ncols,ttype,tform,status)
  if( status .ne. 0 ) write(*,*)"FTICLS fucked up"
  
! Get t0 from the level 2 event file
  status = 0
  call ftgkyd(unit_l2,'TSTART',t0,comment,status)
  if(status .ne. 0) stop 'Cannot determine TSTART'
  
! Spurious modulation files: Shift to  extension "MODULATION"
  status = 0
  call ftmnhd(spunit,2,'MODULATION',0,status)
  if( status .ne. 0 ) stop 'cannot shift to extension MODULATION'
  
! Get the column number for required inputs
  exact=.false.
  !Level 2 events
  call ftgcno(unit_l2,exact,'TIME',tcol,status)
  call ftgcno(unit_l2,exact,'PI',picol,status)
  call ftgcno(unit_l2,exact,'Q',qcol,status)
  call ftgcno(unit_l2,exact,'U',ucol,status)
  !Level 1 events
  call ftgcno(unit_l1,exact,'TIME',tcol1,status)
  call ftgcno(unit_l1,exact,'DETPHI2',phicol,status)
  call ftgcno(unit_l1,exact,'ABSX',xcol,status)
  call ftgcno(unit_l1,exact,'ABSY',ycol,status)
  !Spurious polarisation maps
  call ftgcno(spunit,exact,'DETQ_SM',qsmcol,status)
  call ftgcno(spunit,exact,'DETU_SM',usmcol,status)
  call ftgcno(spunit,exact,'PI',chnscol,status)
  
! Go through each level 2 event
  k1        = 0
  time      = -1d30
  time_prev = -1d30
  bad       = 0
  good      = 0
  offgrid   = 0

  sum_q1_sky = 0.0
  sum_u1_sky = 0.0
  sum_q2_sky = 0.0
  sum_u2_sky = 0.0
  sum_qsp_sky = 0.0
  sum_usp_sky = 0.0
  
  do k = 1,nevts_l2,1
     quality = 1
     status  = 0
     !Read in level 2 event
     call ftgcvj(unit_l2,picol,k,1,1,-1.0,PIchn,anynull,status)
     call ftgcvd(unit_l2,tcol,k,1,1,-1.0,time_l2,anynull,status)
     call ftgcvd(unit_l2,qcol,k,1,1,-1.0,q2,anynull,status)
     call ftgcvd(unit_l2,ucol,k,1,1,-1.0,u2,anynull,status)
     time_l2 = time_l2 - t0
     !Find corresponding level 1 event
     call match(unit_l1,k1,time,time_prev,time_l2,t0,nevts_l1,tcol1)
     !Read in level 1 event
     call ftgcvd(unit_l1,phicol,k1,1,1,-1.0,phi1,anynull,status)
     call ftgcve(unit_l1,xcol,k1,1,1,-1.0,X,anynull,status)
     call ftgcve(unit_l1,ycol,k1,1,1,-1.0,Y,anynull,status)
     !Calculate the correct row to read from the spurious pol map
     DETX_PX = ceiling( ( X + 7.5 ) * 300.0 / 15.0 )
     DETY_PX = ceiling( ( Y + 7.5 ) * 300.0 / 15.0 )
     row     = ( DETX_PX - 1 ) * 300 + DETY_PX
     !Set quality flag for off grid events
     if( DETX_PX .gt. 300 ) quality = 0
     if( DETX_PX .lt. 1   ) quality = 0
     if( DETY_PX .gt. 300 ) quality = 0
     if( DETY_PX .lt. 1   ) quality = 0
     if( quality .eq. 0 ) offgrid = offgrid + 1
     !Read in spurious polarisation for this pixel
     status = 0
     do j = 1,sm_chns
        call ftgcvj(spunit,chnscol,row,j,1,-1.0,DET_CHN_SM(j),anynull,status)
        call ftgcvd(spunit,qsmcol ,row,j,1,-1.0,DETQ_SM(j)   ,anynull,status)
        call ftgcvd(spunit,usmcol ,row,j,1,-1.0,DETU_SM(j)   ,anynull,status)
     end do
     !Interpolate spurious polarisation for this event
     call sm_interp(sm_chns,DET_CHN_SM,DETQ_SM,DETU_SM,PIchn,qsp,usp)
     !Subtract spurious polarisation from level 1 Stokes parameters
     q1 = 2.0 * cos( 2.0 * phi1 )
     u1 = 2.0 * sin( 2.0 * phi1 )   
     qc = q1 - qsp
     uc = u1 - usp
     !Set quality flag for events that fail the PD test
     pc = sqrt( qc**2 + uc**2 )
     p2 = sqrt( q2**2 + u2**2 )
     if( abs(pc-p2) .gt. tol ) quality = 0
     !Calculate roll angle
     Delta = atan2(uc,qc) + atan2(u2,q2)
     if( Delta .gt. pi  ) Delta = Delta - 2.0*pi
     if( Delta .lt. -pi ) Delta = Delta + 2.0*pi
     Delta = 0.5 * Delta
     !Write out roll angle
     if( time_l2 .lt. 10000.0 )then
        if( quality .eq. 1 ) write(deltaunit,*)time_l2,2.0*Delta*180.0/pi
     end if
     !Calculate raw modulation angle in the sky frame

     !Defined on the interval -pi/2 to +pi/2
     phi = 2.0*Delta - 2.0*phi1
     if( phi .gt. pi  ) phi = phi - 2.0*pi
     if( phi .lt. -pi ) phi = phi + 2.0*pi
     phi = 0.5 * phi
     
     ! !Defined on the interval -pi to +pi
     ! phi = Delta - phi1     
     ! if( phi .gt. pi  ) phi = phi - 2.0*pi
     ! if( phi .lt. -pi ) phi = phi + 2.0*pi
     
     !Rotate spurious Stokes parameters to the sky frame
     Delq = cos( 2.0*Delta )
     Delu = sin( 2.0*Delta )
     qsp_sky = Delq * qsp + Delu * usp
     usp_sky = Delu * qsp - Delq * usp
     !Count up quality flags
     if( quality .eq. 1 )then
        good = good + 1
        !Sum up Stokes parameters to check result
        sum_q1_sky = sum_q1_sky + 2.0*cos(2.0*phi)
        sum_u1_sky = sum_u1_sky + 2.0*sin(2.0*phi)
        sum_q2_sky = sum_q2_sky + q2
        sum_u2_sky = sum_u2_sky + u2
        sum_qsp_sky = sum_qsp_sky + qsp_sky
        sum_usp_sky = sum_usp_sky + usp_sky
     else
        bad  = bad + 1
     end if
     !Append to level 1pt 5 event list
     call ftpcld(outunit,ncol+1,k,1,1,phi,status)
     call ftpclj(outunit,ncol+2,k,1,1,quality,status)
     call ftpcld(outunit,ncol+3,k,1,1,Delta,status)
     call ftpcld(outunit,ncol+4,k,1,1,qsp_sky,status)
     call ftpcld(outunit,ncol+5,k,1,1,usp_sky,status)
     !Catch badly behaved rows
     if( status .ne. 0 )then
        status = 0
     end if
     if( anynull )then
        DETQ_SM = 0.0
        DETU_SM = 0.0
        anynull = .false.
     end if     
  end do
  
! Write out summary statistics
  write(*,*)"total events = ",nevts_l2
  write(*,*)"good events = ",good
  write(*,*)"on grid bad events = ",bad-offgrid
  write(*,*)"off grid events",offgrid

! Normalise summed Stokes parameters
  sum_q1_sky = sum_q1_sky   / dble(good)
  sum_u1_sky = sum_u1_sky   / dble(good)
  sum_q2_sky = sum_q2_sky   / dble(good)
  sum_u2_sky = sum_u2_sky   / dble(good)
  sum_qsp_sky = sum_qsp_sky / dble(good)
  sum_usp_sky = sum_usp_sky / dble(good)

! Write out summed Stokes parameters as a check
  write(*,*)"q2_sky=",sum_q2_sky
  write(*,*)"q1_sky-qsp_sky=",sum_q1_sky-sum_qsp_sky
  write(*,*)"u2_sky=",sum_u2_sky
  write(*,*)"u1_sky-usp_sky=",sum_u1_sky-sum_usp_sky
  
! Close out unit and free up unit number
  call ftclos(outunit, status)
  call ftfiou(outunit, status)
  
  return
end subroutine getphi
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine sm_interp(sm_chns,DET_CHN_SM,DETQ_SM,DETU_SM,PIchn,qsp,usp)
  implicit none
! Input
  integer sm_chns,DET_CHN_SM(sm_chns),PIchn
  double precision DETQ_SM(sm_chns),DETU_SM(sm_chns)
! Output
  double precision qsp,usp
! Internal
  integer j
  !Find SM channel above and below PIchn
  j = 2
  do while( DET_CHN_SM(j) .lt. PIchn .and. j .lt. sm_chns )
     j = j + 1
  end do
  !Interpolate (automatically extrapolates if needed)
  qsp = ( DETQ_SM(j) - DETQ_SM(j-1) ) * ( PIchn - DET_CHN_SM(j-1) )
  qsp = qsp / ( DET_CHN_SM(j) - DET_CHN_SM(j-1) )
  qsp = qsp + DETQ_SM(j-1)
  usp = ( DETU_SM(j) - DETU_SM(j-1) ) * ( PIchn - DET_CHN_SM(j-1) )
  usp = usp / ( DET_CHN_SM(j) - DET_CHN_SM(j-1) )
  usp = usp + DETU_SM(j-1)
  !For PIchn above the highest, just use the value for the highest
  if( PIchn .gt. DET_CHN_SM(sm_chns) )then
     qsp = DETQ_SM(j)
     usp = DETU_SM(j)
  end if
  return
end subroutine sm_interp
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine match(unit_l1,k1,time,time_prev,time_l2,t0,nevts_l1,tcol1)
! Find the row (k1) in the level 1 event list that matches the time of
! the current level 2 event.
! time,time_prev,and k1 must already be initialised.
  implicit none
  integer unit_l1,k1,nevts_l1,tcol1
  double precision time,time_prev,time_l2,t0
  integer status
  logical anynull
  double precision dt_right,dt_left,time_l1
  status = 0
  do while( time .lt. time_l2 .and. k1 .lt. nevts_l1 )
     k1        = k1 + 1
     time_prev = time
     call ftgcvd(unit_l1,tcol1,k1,1,1,-1.0,time,anynull,status)
     time      = time - t0
  end do
  dt_right = time - time_l2
  dt_left  = time_l2 - time_prev
  if( dt_right .lt. dt_left )then
     time_l1 = time
  else
     time_l1 = time_prev
     k1      = k1 - 1
  end if
  return
end subroutine match
!-----------------------------------------------------------------------



