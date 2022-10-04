!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_gritas_grids --- GRITAS dimension and grid parameters
!
! !INTERFACE:

      MODULE m_gritas_grids

! !USES:

      Implicit none

! !DESCRIPTION: 
!
!  Defines dimensions and other relevant parameters.
!
! !DEFINED PARAMETERS: 

!     Grid dimension
!     --------------
      integer, save ::     im                        ! no. zonal      gridpoints
      integer, save ::     jm                        ! no. meridional gridpoints
      integer, save ::     km                        ! max.   no. vertical levels
      integer, save ::     nlevs                     ! actual no. vertical levels
      integer, save ::     nchan                     ! no. of active channels (when applicable)

!     Longitudinal grid
!     -----------------
      real, parameter         :: LonMin = -180.      ! = lambda(1)
      real, save              :: dLon                ! mesh size (deg)
      real, allocatable, save :: glon(:)             ! zonal gridpoints

!     real, parameter :: dLon = 360. / im    ! mesh size (deg)
   

!     Latitudinal grid
!     ----------------
      real, parameter         :: LatMin = -90.       ! = phi(1)
      real, save              :: dLat                ! mesh size (deg)
      real, allocatable, save :: glat(:)             ! meridional gridpoints

!     Logitudinal slots mimicking time zones
!     --------------------------------------
      integer, parameter :: ntslots = 1 ! 360/15
      integer, allocatable, save :: lonindx(:,:)     ! 
      integer, allocatable, save :: regindx(:,:)     ! known (latitudinal) region indeces
      real, allocatable, save :: regbounds(:,:)      ! known (latitudinal) region boundaries
      integer, parameter :: mregs = 4
      character(len=7) :: regnames(mregs) = (/'global ','n.hem  ','tropics','s.hem  '/)

!     Collected statistics for scorecard
!     ----------------------------------
      integer, parameter :: nevals = 3 ! number of evaluated statistical parameters, e.g., mean, stdev, rms, efficiency, etc
      character(len= 8), parameter :: type_stats(nevals) = (/'mean','rms','stdv'/)

!     Vertical grid
!     -------------
      real, allocatable, save :: plevs(:)            ! vertical levels in hPa
      real, allocatable, save :: hlevs(:)            ! heights from hyrostatic balance (apprx; m)

!     Active Channels
!     ---------------
      integer, save :: ichan_version                 ! version=1: uses chn index
                                                     ! version=2: uses chn number
      integer, allocatable, save :: achan(:,:)       ! active channels: index, and chn number

!     Missing value
!     -------------
      real, parameter :: AMISS = 1.0E15


! !REVISION HISTORY: 
!
!  05sep97   da Silva   Original code.
!  09Jun06   Todling    Converted to module.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname = 'm_Gritas_Grids'

      interface grids_init  ; module procedure init_        ; end interface
      interface grids_clean ; module procedure clean_       ; end interface
      interface GrITAS_Pint ; module procedure GrITAS_Pint_ ; end interface
      interface GrITAS_Norm ; module procedure GrITAS_Norm_ ; end interface
      interface GrITAS_Accum; module procedure GrITAS_Accum_; end interface
      interface GrITAS_MinMax; module procedure GrITAS_MinMax_; end interface
      interface GrITAS_XCov_Accum; module procedure GrITAS_XCov_Accum_; end interface
      interface GrITAS_XCov_Norm ; module procedure GrITAS_XCov_Norm_ ; end interface
      interface GrITAS_3CH_Norm ; module procedure GrITAS_3CH_Norm_ ; end interface
      interface GrITAS_XCov_Prs_Accum ; module procedure GrITAS_XCov_Prs_Accum_ ; end interface
      interface GrITAS_CDF ; module procedure GrITAS_cdf_ ; end interface
      interface GrITAS_Strip ; module procedure GrITAS_Strip_ ; end interface
      interface GrITAS_HydroHeights ; module procedure hydro_heights_ ; end interface
      interface GrITAS_Scores
         module procedure Scores1_ 
         module procedure Scores2_ 
      end interface

      CONTAINS

      subroutine init_ (im_def,jm_def,km_def,stat)

      implicit none
      integer, intent(in)  :: im_def
      integer, intent(in)  :: jm_def
      integer, intent(in)  :: km_def
      integer, intent(out) :: stat

      integer  i, j, ier
      real :: lona, lonb

      stat = 0

      im = im_def
      jm = jm_def
      km = km_def

      allocate ( glon(im), glat(jm), plevs(km), stat=ier )
         if(ier/=0) then
            stat = 1
            return
         endif

      dLon = 360. / im
      dLat = 180. / ( jm - 1 )

!     Create horizontal grids for GFIO's sake
!     ---------------------------------------
      do i = 1, im
         glon(i) = LonMin + (i-1) * dLon
      end do
      do j = 1, jm
         glat(j) = LatMin + (j-1) * dLat
      end do

!     Create regional info
!     --------------------
      allocate(regindx(2,size(regnames)))
      allocate(regbounds(2,size(regnames)))
      ! this is to be placed in the RC file as table
      do i=1,size(regnames)
         if(trim(regnames(i))== 'global') then
            regbounds(1,i) = -90.0; regbounds(2,i) =  90.0
         endif
         if(trim(regnames(i))== 's.hem') then
            regbounds(1,i) = -80.0; regbounds(2,i) = -20.0
         endif
         if(trim(regnames(i))== 'tropics') then
            regbounds(1,i) = -20.0; regbounds(2,i) =  20.0
         endif
         if(trim(regnames(i))== 'n.hem') then
            regbounds(1,i) =  20.0; regbounds(2,i) =  80.0
         endif
      end do
      print *
      write(6,'(a)') 'Wired regions for scorecard:'
      write(6,'(a)') '----------------------------'
      do i=1,size(regnames)
         do j = 1, jm
            if(glat(j)<=regbounds(1,i)) then
               regindx(1,i) = j
            endif
            if(glat(j)>0.0 .and. glat(j)<=regbounds(2,i)) then
               regindx(2,i) = j
            endif
            if(glat(j)<0.0 .and. glat(j)<=regbounds(2,i)) then
               regindx(2,i) = j
            endif
         end do
         write(6,'(a10,2(2x,f6.2),2(2x,i4))') trim(regnames(i)), glat(regindx(1,i)), glat(regindx(2,i)), regindx(1,i), regindx(2,i)
      end do

!     Create hour slot (based on approx time-zones
      allocate(lonindx(2,ntslots))
      do j = 1, ntslots
         lona = -180 + (j-1) * 360.0
         lonb = -180 +  j    * 360.0
         lonindx(1,1)=1 
         lonindx(2,1)=im
      enddo
      do j = 1,ntslots
         lona =        (j-1) * 360.0
         lonb =         j    * 360.0
         write(6,'(i3,2(2x,f7.2),2(2x,i5))') j, lona, lonb, lonindx(1,j), lonindx(2,j)
      enddo
      end subroutine init_

      subroutine clean_(stat)
      implicit none
      integer, intent(out) :: stat
      integer ier
      stat = 0
      deallocate(lonindx)
      deallocate(regindx)
      deallocate ( glon, glat, plevs, stat=ier )
         if(ier/=0) then
            stat = 1
            return
         endif
!     if (allocated(hlevs)) deallocate(hlevs)
      end subroutine clean_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Pint_ --- Set pressure Intervals for gridding
! 
! !INTERFACE:
!
      subroutine GrITAS_Pint_ ( p1, p2 )

! !INPUT PARAMETERS: 
 
! !OUTPUT PARAMETERS:
!
      real       p1(nlevs), p2(nlevs)   ! pressure intervals for gridding
!
! !DESCRIPTION: Determine the pressure level ranges for binning
!               off-mandatory level data. Observations with pressure
!               p1(k) <= p < p2(k) will be gridded at level plevs(k). 
!               The binnings are log-pressure linear. 
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  17Oct97  G. P. Lou  Modification.
!  25Mar98  da Silva   Removed tighning of vertical level interval
!                      introduced by G. P. Lou.
!  24Jun04  Todling    Changed handling of 1-level case.
!
!EOP
!-------------------------------------------------------------------------


      character(len=*), parameter :: myname_ = myname//'GrITAS_Pint_'

      integer i, j, k

!     Check number of levels
!     ----------------------
      if ( nlevs .le. 0 ) call gritas_perror
     &   ( 'Pint: must have at least 1 vertical levels' )

!     Special single level case
!     -------------------------
      if ( nlevs .eq. 1 ) then
           p1(1) = plevs(1)
           p2(1) = plevs(1) - 0.1
           print *, '       p1         p        p2'
           print 100, p1(1), plevs(1), p2(1)
           print *
           return
      end if


!     Check monotonicity
!     ------------------
      do k = 1, nlevs-1
         if ( plevs(k) .le. plevs(k+1) ) call gritas_perror
     &      ('Pint: pressure levels do not decrease monotonically')
      end do

      if ( plevs(nlevs) .eq. 0. ) call gritas_perror
     &   ( 'Pint: cannot handle 0 hPa' )

      k = 1
      p1(1) = 2000. 
      p2(1) = exp ( 0.5 * ( log(plevs(k+1)) + log(plevs(k)) ) )

      do k = 2, nlevs-1
         p1(k) = exp ( 0.5 * ( log(plevs(k-1)) + log(plevs(k)) ) )
         p2(k) = exp ( 0.5 * ( log(plevs(k+1)) + log(plevs(k)) ) )
      end do

      k = nlevs
      p1(k) = exp ( 0.5 * ( log(plevs(k-1)) + log(plevs(k)) ) )
      p2(k) = 0.
       
!     Echo parameters
!     ---------------
      print *, '       p1         p        p2'
      do k = 1, nlevs
         print 100, p1(k), plevs(k), p2(k)
      end do
 100  format(1x,3F12.3)

      print *

      end subroutine GrITAS_Pint_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GrITAS_Accum_ --- Accumulate O-F on del grids
! 
! !INTERFACE:
!
      subroutine GrITAS_Accum_ ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                           kxkt_table, p1, p2,
     &                           ob_flag,
     &                           l2d, bias_2d, rtms_2d, xcov_2d, nobs_2d,
     &                           l3d, bias_3d, rtms_3d, xcov_3d, nobs_3d )

! !USES:
!

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lat(nobs)     ! latitude (degrees)
      real,    intent(in) ::       lon(nobs)     ! longitude (degrees)
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index

      
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      real,    intent(in) ::           p1(nlevs)
      real,    intent(in) ::           p2(nlevs)     ! pressure intervals for gridding

      integer,    intent(in) ::         l2d          ! actual number of 2d grids
      integer,    intent(in) ::         l3d          ! actual number of 3d grids

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer,    intent(inout) :: ob_flag(nobs)     ! On input, ob is rejected if
                                                     !    ob_flag is other than 0
                                                     ! On output, ob_flag = 0 if ob
                                                     !    was used to accumulate the
                                                     !    bias/RMS
                                                     ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(im,jm,l2d,nr)    !   time means
      real, intent(inout) ::  rtms_2d(im,jm,l2d,nr)    !   root time mean square
      real, intent(inout) ::  xcov_2d(im,jm,l2d)       !   x (co)variance
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                     ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(im,jm,km,l3d,nr) !  time means
      real, intent(inout) ::  rtms_3d(im,jm,km,l3d,nr) !   root time mean square
      real, intent(inout) ::  xcov_3d(im,jm,km,l3d)    !   root time mean square
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box

!
! !DESCRIPTION: Given the O-F for a given synoptic time, this routine
!               accumulates the bias and the RMS at the nearest gridpoint.
!
! !REMARKS:
!
!   1.  xrtms is redundant when nr=1, i.e., xrtms=rtms
!   2.  when nr=2, assumes residuals are from same experiment, that is,
!       no need to distinguish between nobs
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  25Mar98  da Silva   Removed "goto 10" statement; now all must be
!                      assigned a level or else a fatal error will occur
!                      (this was the original design).
!  26Mar98  da Silva   Included periodic boundary conditions in longitude.
!  14Jul99  da SIlva   Added land-surface kt's
!  28oct06  Todling    Generalized to handle cross-terms; note xrtms is 
!                      redundant when nr=1, i.e., xrtms=rtms
!  28May09  Redder     Changed code to define negative value in
!                      kxkt_table as indices for surface variables.
!  07Jul09  Redder     Added parameters ob_flag
!  07Jan22  Todling    Accommodate accum for 3CH (nr=3); in this case xcov left
!                      untouched
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_Accum_'

      integer n, i, j, k, kk, nx, l, ir
      logical surface_var, reject_ob

      if ( nr > 3 ) call gritas_perror ('Accum: could not handle '
     &                               // 'more then 2 residuals')
      if ( nr==3 ) then
         nx = 0
      else
         nx = max(1,nr)
      endif

!     Loop over observations
!     ----------------------
      do 10 n = 1, nobs
         reject_ob = ob_flag ( n ) .ne. 0
         if ( reject_ob  ) go to 10                    ! reject unwanted obs

         l = kxkt_table(kx(n),kt(n))
         surface_var = kxkt_table (kx(n),kt(n)) .lt. 0
         l = abs (l)

         if ( l .le. 0 ) then
            ob_flag ( n ) = 1
            go to 10                                   ! forget it
         end if
         
!        Horizontal gridbox: gridpoint is assumed at center of box
!        ---------------------------------------------------------
         i = 1 + nint ( (lon(n) - LonMin) / dLon )
         if ( i .eq. (im+1) ) i = 1                      ! periodic bc
         if ( i .eq.   0    ) i = im                     ! periodic bc
         j = 1 + nint ( (lat(n) - LatMin) / dLat )

         i = min(im,max(1,i))
         j = min(jm,max(1,j))
        

 
!        Surface (2D) grid
!        -----------------
         if ( surface_var )  then  ! surface variables 

!             Accumulate
!             ----------              
              do ir = 1, nr
                 bias_2d(i,j,l,ir) = bias_2d(i,j,l,ir) + del(n,ir)
                 rtms_2d(i,j,l,ir) = rtms_2d(i,j,l,ir) + del(n,ir)*del(n,ir)
              end do
              if(nx>0) xcov_2d(i,j,l)       = xcov_2d(i,j,l)    + del(n,1) *del(n,nx)
                       nobs_2d(i,j,l)       = nobs_2d(i,j,l)    + 1.


!        Upper-air (3D) grid
!        -------------------
         else

!             Find vertical level
!             -------------------
              do kk = 1, nlevs
                 if ( lev(n) .le. p1(kk) .AND. 
     &                lev(n) .gt. p2(kk) ) then
                      k = kk
                      go to 21
                 end if
              end do

              print *, 'Accum: obs level is ', lev(n), ' hPa'
              call gritas_perror ('Accum: could not find level')

 21           continue

!             Accumulate
!             ----------              
              do ir = 1, nr
                 bias_3d(i,j,k,l,ir) = bias_3d(i,j,k,l,ir) + del(n,ir)
                 rtms_3d(i,j,k,l,ir) = rtms_3d(i,j,k,l,ir) + del(n,ir)*del(n,ir)
              end do
              if(nx>0) xcov_3d(i,j,k,l)       = xcov_3d(i,j,k,l)    + del(n,1) *del(n,nx)
                       nobs_3d(i,j,k,l)       = nobs_3d(i,j,k,l)    + 1.

           end if

 10     continue



        end subroutine GrITAS_Accum_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Norm_ --- Normalize time mean statistics
! 
! !INTERFACE:
!
       subroutine GrITAS_Norm_ ( nr, nrmlz,
     &                           l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                           l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,
     &                           ntresh )
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::      nr           ! number of distinct residual series
      logical, intent(in) ::      nrmlz        ! normalize or not
      integer, intent(in) ::      l2d          ! actual number of 2d grids
      integer, intent(in) ::      l3d          ! actual number of 3d grids
      integer, OPTIONAL, intent(in) :: ntresh  ! tresh hold to compute means

! !INPUT/OUTPUT PARAMETERS:
!
                                                      ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(im,jm,l2d,nr)    !   time means
      real, intent(inout) ::  rtms_2d(im,jm,l2d,nr)    !   root time mean square
      real, intent(inout) ::  xcov_2d(im,jm,l2d,nr)    !   x-(co)variance
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                      ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(im,jm,km,l3d,nr) !   time means
      real, intent(inout) ::  rtms_3d(im,jm,km,l3d,nr) !   root time mean square
      real, intent(inout) ::  xcov_3d(im,jm,km,l3d,nr) !   x-(co)variance
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box
!
! !OUTPUT PARAMETERS:
!
      real, intent(inout) ::  stdv_2d(im,jm,l2d,nr)    ! 2D root time mean square
      real, intent(inout) ::  stdv_3d(im,jm,km,l3d,nr) ! 3D root time mean square

!
! !DESCRIPTION:  Normalizes the time mean statistics (such as
!                dividing the accumulated O-F by the number of
!  observations), set to 'missing value' the gridboxes with no
!  data.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  26Mar98  da Silva   Change normalization: n=1 is now OK for means.
!                      I needed this for being able to use the same
!                      code for synoptic time averages (GrISAS).
!  24Jun04  Todling    Added protection for nlev<2.
!  12oct05  da Silva   Added ntresh to impose a min number of obs
!  28Oct06  Todling    Generalized to handle cross-terms.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_Norm_'

      integer i, j, k, l, m, nx, ir, ntresh_
      real    tmp(nr), xtmp, n

      integer icount
      real, allocatable :: glo2dmean(:,:)
      real, allocatable :: glo3dmean(:,:,:), glo3drtms(:,:,:)

      if ( nr > 2 ) call gritas_perror ('Accum: could not handle more then 2 residuals')
      nx = max(1,nr)

      ntresh_ = 0                   ! by default no effect
      if ( present(ntresh) ) then
           ntresh_ = ntresh
      end if


!     Surface (2D) grids
!     ------------------
      m=1
      do l = 1, l2d
         do j = 1, jm
            do i = 1, im

               n = nobs_2d(i,j,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = bias_2d(i,j,l,ir) / m                                ! mean(x)
                       tmp(ir)           = rtms_2d(i,j,l,ir) / m                                ! mean(x**2)
                       rtms_2d(i,j,l,ir) = sqrt(tmp(ir))
                       tmp(ir)           = abs(tmp(ir)  - bias_2d(i,j,l,ir) *bias_2d(i,j,l,ir)) ! mean(x'**2)
                    end do
                    xtmp             = xcov_2d(i,j,l,1) / m                                    ! <x,y>
                    xcov_2d(i,j,l,1) = xtmp                                                    ! save xRMS in slot 1
                    xtmp             = xtmp  - bias_2d(i,j,l,1) *bias_2d(i,j,l,nx)             ! <x',y'>
               else
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = AMISS
                       rtms_2d(i,j,l,ir) = AMISS
                    end do
                       xcov_2d(i,j,l,1)  = AMISS
               end if

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr 
                       stdv_2d(i,j,l,ir)  = sqrt( m * tmp(ir) / (m-1) )
                    end do
                       xcov_2d(i,j,l,nx)  =       m * xtmp    / (m-1) 
               else
                    do ir = 1, nr 
                       stdv_2d(i,j,l,ir)  = AMISS
                    end do
                       xcov_2d(i,j,l,nx)  = AMISS
               end if

            end do
         end do
      end do

      if( nlevs<2 ) return

!     Upper-air (3D) grids
!     --------------------

      m=1
      do l = 1, l3d
         do k = 1, nlevs
            do j = 1, jm
               do i = 1, im
                  
               n = nobs_3d(i,j,k,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = bias_3d(i,j,k,l,ir) / m                        ! mean(x) and mean(y)
                       tmp(ir)             = rtms_3d(i,j,k,l,ir) / m                        ! mean(x**2) and mean(y**2)
                       rtms_3d(i,j,k,l,ir) = sqrt(tmp(ir))
                       tmp(ir)  = abs(tmp(ir)  - bias_3d(i,j,k,l,ir) *bias_3d(i,j,k,l,ir))  ! mean(x'**2) and mean(y'**2)
                    end do
                    xtmp               = xcov_3d(i,j,k,l,1) / m                             ! <x,y>
                    xcov_3d(i,j,k,l,1) = xtmp                                               ! save xRMS in slot 1 
                    xtmp               = xtmp  - bias_3d(i,j,k,l,1) *bias_3d(i,j,k,l,nx)    ! <x',y'>
               else
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = AMISS
                       rtms_3d(i,j,k,l,ir) = AMISS
                    end do
                       xcov_3d(i,j,k,l,1)  = AMISS
               endif

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    stdv_3d(i,j,k,l,1:nr)  = sqrt( m * tmp(1:nr) / (m-1) )
                    xcov_3d(i,j,k,l,nx)    =       m * xtmp      / (m-1)
               else
                    if (n==1) then 
                       stdv_3d(i,j,k,l,1:nr)  = sqrt( tmp(1:nr) )
                       xcov_3d(i,j,k,l,nx)    =       xtmp
                    else
                       stdv_3d(i,j,k,l,1:nr) = AMISS
                       xcov_3d(i,j,k,l,nx)   = AMISS
                    endif
               end if

               end do
            end do
         end do
      end do

      end subroutine GrITAS_Norm_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GrITAS_XCov_Accum_ --- X-covariance accumulator
! 
! !INTERFACE:
!
      subroutine GrITAS_XCov_Accum_ ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                                kxkt_table,
     &                                ob_flag,
     &                                l3d, bias_lv, xcov_lv, nobs_lv, nprof, tch )

! !USES:
!

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lat(nobs)     ! latitude (degrees)
      real,    intent(in) ::       lon(nobs)     ! longitude (degrees)
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index

      
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      integer,    intent(in) ::         l3d          ! actual number of 3d grids
      logical,    intent(in) ::         tch

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer,    intent(inout) :: ob_flag(nobs)     ! On input, ob is rejected if
                                                     !    ob_flag is other than 0
                                                     ! On output, ob_flag = 0 if ob
                                                     !    was used to accumulate the
                                                     !    bias/RMS

                                                     ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_lv(:,:,:)         !   time means
      real, intent(inout) ::  xcov_lv(nchan,nchan,l3d,nr)  !   root time mean square
      real, intent(inout) ::  nobs_lv(nchan,l3d)     !   number of obs for each case

      integer, intent(out) :: nprof
!
! !DESCRIPTION: Given two residual sequences (nr) at given synoptic time, 
!               this routine accumulates bias and RMS at to construct 
!               an interchannel (cross) covariance.
!
! !REVISION HISTORY: 
!
!  11Apr14  Todling    Adapt from original to estimate interchannel (cross)
!                      covariances "a la Desroziers method".
!  02Aug14  Todling    Revisit algorithm; now works with real channel numbers
!  08May22  Todling    Accommodate 3CH
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_XCov_Accum_'

      integer n1, n2, m1, k1, k2, ii, jj, k, kk, nx, l, ir, nrr
      logical surface_var, docalc

      logical,allocatable,dimension(:):: reject_ob
      integer,allocatable,dimension(:):: indx,ll
      real,allocatable,dimension(:)   :: rlon,rlat,rlev
      real,allocatable,dimension(:,:) :: rdel

      if ( nr > 3 ) then
          call gritas_perror (myname_//': improper dim settings')
      endif
      if (tch) then
        if (nr/=3) then
          call gritas_perror (myname_ //': TCH setting inconsistent w/ dims')
        endif
      else
        if (nr/=2) then
          call gritas_perror (myname_ //': Desroziers setting inconsistent w/ dims')
        endif
      endif

!     Check to see that a set of complete profiles is present
!     -------------------------------------------------------
      if( mod(nobs,km) == 0 ) then
        print *
        print *, 'XCov_Accum: found ', nobs/km, ' total complete profiles: '
      else
        print *, 'XCov_Accum:  nobs, profiles ', nobs, km
        call gritas_perror ('XCov_Accum: not all profiles complete, aborting ...')
      endif
      allocate(indx(nchan))
      allocate(ll(nchan))
      allocate(reject_ob(nchan))
      allocate(rlon(nchan),rlat(nchan),rlev(nchan),rdel(nchan,nr))
      indx = 0
      rlon = amiss
      rlat = amiss
      rlev = amiss
      rdel = amiss
      
!     Loop over observations
!     ----------------------
      n1 = 0; n2 = 0
      nprof = 0
      do n1 = 1, nobs, km
         n2 = n1 + km - 1
         m1 =0
         reject_ob=.false.
         do ii = n1, n2
            k=0
            do kk = 1,nchan
               if(abs(achan(kk,ichan_version)-nint(lev(ii)))<1.e-5) then
                  k=kk
                  go to 9
               endif
            enddo
    9       continue
            if(k>0) then ! ob channel among active channels 
               m1 = m1 + 1
               if(m1>nchan) 
     &            call gritas_perror ('XCov-Accum: error in channel matching')
               indx(m1) = ii  ! index for active channel
               reject_ob(m1) = ob_flag ( ii ) .ne. 0
               ll(m1) = kxkt_table(kx(ii),kt(ii))
               ll(m1) = abs (ll(m1))
               if ( ll(m1) .le. 0 ) then
                  ob_flag ( ii ) = 1
                  reject_ob(m1) = .true. ! quality mark of active channel
               end if
            endif
         enddo
         if(ANY(reject_ob)) cycle ! if profile of active channels has rejected ob, skip profile
         if(m1/=nchan) cycle      ! this shouldn't happen, nonetheless skip

!        Beyond this point, only considering active channels of profiles whose ob are fully accepted
!        -------------------------------------------------------------------------------------------
         rlon(:)   = lon(indx(:)) 
         rlat(:)   = lat(indx(:)) 
         rlev(:)   = lev(indx(:)) 
         rdel(:,:) = del(indx(:),:) 

         nprof = nprof + 1

!        Loop over observations
!        ----------------------
         do kk = 1, nchan
            l = ll(kk)
            bias_lv(kk,l,:) = bias_lv(kk,l,:)  + rdel(kk,:)

!           Accumulate
!           ----------              
            if (tch) then    ! for 3CH method
               do ii = 1, nchan
                  xcov_lv(ii,kk,l,:) = xcov_lv(ii,kk,l,:) + rdel(ii,:)*rdel(kk,:)
               enddo
            else             ! for Desroziers et al
               do ii = 1, nchan
                  xcov_lv(ii,kk,l,1) = xcov_lv(ii,kk,l,1) + rdel(ii,1)*rdel(kk,2)
! Low-level symmetric version of Desroziers
!                 xcov_lv(ii,kk,l,1) = xcov_lv(ii,kk,l,1) + 0.5*(rdel(ii,1)*rdel(kk,2)+rdel(ii,2)*rdel(kk,1))
               enddo
            endif
            nobs_lv(kk,l) = nobs_lv(kk,l) + 1.
         enddo

         indx = 0
         rlon = amiss
         rlat = amiss
         rlev = amiss
         rdel = amiss

      enddo
      print*, 'GrITAS_XCov_Accum, total profiles processed = ', nprof

      deallocate(rlon,rlat,rlev,rdel)
      deallocate(reject_ob)
      deallocate(ll)
      deallocate(indx)

      end subroutine GrITAS_XCov_Accum_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_XCov_Norm_ --- Normalize time mean statistics
! 
! !INTERFACE:
!
       subroutine GrITAS_XCov_Norm_ ( nr, nrmlz, tch, self,
     &                                l3d, bias_3d, xcov_3d, nobs_3d,
     &                                ntresh )
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::      nr           ! number of distinct residual series
      logical, intent(in) ::      nrmlz        ! normalize or not
      logical, intent(in) ::      tch          ! 3CH method
      logical, intent(in) ::      self         ! allows getting R=(o-<o>)^2, same for B and A
      integer, intent(in) ::      l3d          ! actual number of 3d grids
      integer, OPTIONAL, intent(in) :: ntresh  ! tresh hold to compute means

! !INPUT/OUTPUT PARAMETERS:
!
                                                      ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(:,:,:)          !   time means
      real, intent(inout) ::  xcov_3d(:,:,:,:)        !   x-(co)variance
      real, intent(inout) ::  nobs_3d(:,:)            !   number of obs for each case

!
! !DESCRIPTION:  Normalizes the time mean statistics (such as
!                dividing the accumulated O-F by the number of
!                observations).
!
! !REVISION HISTORY: 
!
!  11Apr14  Todling    Adapt from original to handle interchannel covariance
!                      case
!  08May22  Todling    Accommodate 3CH
!  26May22  Todling    Added self opt (has anyone looked at is in the 
!                      literature), in the context of Desroziers?
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_XCov_Norm_'

      integer k1, k2, k3, k, l, m, nx, ir, ndim, ntresh_, i1, i2, nrr
      real    tmp(nr), xtmp, n
      real,   allocatable :: var_3d(:,:,:,:)
   
      if ( nr > 3 ) then
          call gritas_perror (myname_ //': improper dim settings')
      endif
      if (tch) then
        if (nr/=3) then
          call gritas_perror (myname_ //': TCH setting inconsistent w/ dims')
        endif
      else
        if (nr/=2) then
          call gritas_perror (myname_ //': Desroziers setting inconsistent w/ dims')
        endif
      endif

      ndim = size(bias_3d,1)
      if(size(bias_3d,2)/=l3d)
     &   call gritas_perror (myname_//': error dim(bias_3d,2)')
      if(size(xcov_3d,1)/=ndim)
     &   call gritas_perror (myname_//': error dim(xcov_3d,1)')
      if(size(xcov_3d,2)/=ndim)
     &   call gritas_perror (myname_//': error dim(xcov_3d,2)')
      if(size(xcov_3d,3)/=l3d)
     &   call gritas_perror (myname_//': error dim(xcov_3d,2)')

      ntresh_ = 0                   ! by default no effect
      if ( present(ntresh) ) then
           ntresh_ = ntresh
      end if

!     Upper-air (3D) grids
!     --------------------
      nrr=1
      if(tch) nrr=3
      do k3=1,nrr
         i1=1; i2=2
         if(nrr==3)then
            i1=k3; i2=k3
         endif
         m=1
         do l = 1, l3d
            do k2 = 1, ndim
               n = nobs_3d(k2,l) 
               if ( nrmlz ) m=n
               if ( m .gt. 1 .AND. m .ge. ntresh_ ) then
                  if(tch) then
                     bias_3d(k2,l,k3) = bias_3d(k2,l,k3) / m        ! mean(x) and mean(y)
                  else
                     bias_3d(k2,l,:)  = bias_3d(k2,l,:)  / m        ! mean(x) and mean(y)
                  endif
               endif
            end do
         end do
         m=1
      enddo
      do k3=1,nrr
         i1=1; i2=2
         if(nrr==3)then
            i1=k3; i2=k3
         endif
         do l = 1, l3d
            do k2 = 1, ndim
               do k1 = 1, ndim
                  n = nobs_3d(k2,l) 
   
                  if ( nrmlz ) m=n
                  if ( m .gt. 1 .AND. m .ge. ntresh_ ) then
                       xtmp = xcov_3d(k1,k2,l,k3) / m                       ! <x,y>
                       if ( tch ) then
                         xtmp = xtmp  - bias_3d(k1,l,i1) *bias_3d(k2,l,i2)  ! <x',y'>
                       else
! Low-level symmetric version of Desroziers
                         xtmp = xtmp  - 0.5*( bias_3d(k1,l,i1) *bias_3d(k2,l,i2) 
     .                                     +  bias_3d(k2,l,i1) *bias_3d(k1,l,i2) ) ! <x',y'>
                       endif
                  endif
      
                  if ( m .gt. 1 .AND. m .ge. ntresh_ ) then
                       xcov_3d(k1,k2,l,k3)   = m * xtmp / (m-1)
                  end if

               end do
            end do
         end do
      end do

      if (.not.tch) return
      if ( self )  return  ! this returns R, B and A, for (o,b,a) corners

      allocate(var_3d(size(xcov_3d,1),size(xcov_3d,2),size(xcov_3d,3),nrr))
      var_3d = xcov_3d

!     Now calculate 3CH estimate (o,b,a)
!     --------------------------
!                                     (x-y)*2           (x-z)*2            (y-z)*2
      xcov_3d(:,:,:,1) = 0.5*(var_3d(:,:,:,1) + var_3d(:,:,:,2) - var_3d(:,:,:,3))
      xcov_3d(:,:,:,2) = 0.5*(var_3d(:,:,:,3) + var_3d(:,:,:,1) - var_3d(:,:,:,2))
      xcov_3d(:,:,:,3) = 0.5*(var_3d(:,:,:,2) + var_3d(:,:,:,3) - var_3d(:,:,:,1))

!     print *, 'debug R : ', sqrt(sum(xcov_3d(:,:,:,1)*xcov_3d(:,:,:,1)))
!     print *, 'debug B : ', sqrt(sum(xcov_3d(:,:,:,2)*xcov_3d(:,:,:,2)))
!     print *, 'debug A : ', sqrt(sum(xcov_3d(:,:,:,3)*xcov_3d(:,:,:,3)))
      deallocate(var_3d)

      end subroutine GrITAS_XCov_Norm_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GrITAS_MinMax_ --- Gather min and max vals within box
! 
! !INTERFACE:
!
      subroutine GrITAS_MinMax_ ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                            kxkt_table, p1, p2,
     &                            ob_flag,
     &                            l2d, minv_2d, maxv_2d, nobs_2d,
     &                            l3d, minv_3d, maxv_3d, nobs_3d )

! !USES:
!

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lat(nobs)     ! latitude (degrees)
      real,    intent(in) ::       lon(nobs)     ! longitude (degrees)
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index

      
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      real,    intent(in) ::           p1(nlevs)
      real,    intent(in) ::           p2(nlevs)     ! pressure intervals for gridding

      integer,    intent(in) ::         l2d          ! actual number of 2d grids
      integer,    intent(in) ::         l3d          ! actual number of 3d grids

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer,    intent(inout) :: ob_flag(nobs)     ! On input, ob is rejected if
                                                     !    ob_flag is other than 0
                                                     ! On output, ob_flag = 0 if ob
                                                     !    was used to accumulate the
                                                     !    bias/RMS
                                                     ! Surface (2D) grids:
      real, intent(inout) ::  minv_2d(im,jm,l2d,nr)    !   min value within box
      real, intent(inout) ::  maxv_2d(im,jm,l2d,nr)    !   max value within box
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                     ! Upper-air (3D) grids:
      real, intent(inout) ::  minv_3d(im,jm,km,l3d,nr) !   min value within box
      real, intent(inout) ::  maxv_3d(im,jm,km,l3d,nr) !   max value within box
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box

!
! !DESCRIPTION: Given an ungridded quantity for a given synoptic time, this routine
!               collects max and min values within the nearest gridpoint.
!
! !REMARKS:
!
! !REVISION HISTORY: 
!
!  08Jul2020  Todling    Initial code, adapted from Accum.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_MinMax_'

      integer n, i, j, k, kk, nx, l, ir
      logical surface_var, reject_ob, lfound

      if ( nr > 1 ) call gritas_perror ('MinMax: could not handle '
     &                                // 'more then 1 field')
      nx = max(1,nr)

!     Loop over observations
!     ----------------------
      do n = 1, nobs
         reject_ob = ob_flag ( n ) .ne. 0
         if ( reject_ob  ) cycle                    ! reject unwanted obs

         l = kxkt_table(kx(n),kt(n))
         surface_var = kxkt_table (kx(n),kt(n)) .lt. 0
         l = abs (l)

         if ( l .le. 0 ) then
            ob_flag ( n ) = 1
            cycle
         end if
         
!        Horizontal gridbox: gridpoint is assumed at center of box
!        ---------------------------------------------------------
         i = 1 + nint ( (lon(n) - LonMin) / dLon )
         if ( i .eq. (im+1) ) i = 1                      ! periodic bc
         if ( i .eq.   0    ) i = im                     ! periodic bc
         j = 1 + nint ( (lat(n) - LatMin) / dLat )

         i = min(im,max(1,i))
         j = min(jm,max(1,j))
        

 
!        Surface (2D) grid
!        -----------------
         if ( surface_var )  then  ! surface variables 

!             Accumulate
!             ----------              
              do ir = 1, nr
                 if( del(n,ir) < minv_2d(i,j,l,ir) ) minv_2d(i,j,l,ir) = del(n,ir)
                 if( del(n,ir) > maxv_2d(i,j,l,ir) ) maxv_2d(i,j,l,ir) = del(n,ir)
              end do
              nobs_2d(i,j,l) = nobs_2d(i,j,l) + 1.


!        Upper-air (3D) grid
!        -------------------
         else

!             Find vertical level
!             -------------------
              lfound = .false.
              do kk = 1, nlevs
                 if ( lev(n) .le. p1(kk) .AND. 
     &                lev(n) .gt. p2(kk) ) then
                      k = kk
                      lfound = .true.
                      exit
                 end if
              end do

              if (.not.lfound) then
                 print *, 'Accum: obs level is ', lev(n), ' hPa'
                 call gritas_perror ('Accum: could not find level')
              endif

!             Accumulate
!             ----------              
              do ir = 1, nr
                 if( del(n,ir) < minv_3d(i,j,k,l,ir) ) minv_3d(i,j,k,l,ir) = del(n,ir)
                 if( del(n,ir) > maxv_3d(i,j,k,l,ir) ) maxv_3d(i,j,k,l,ir) = del(n,ir)
              end do
              nobs_3d(i,j,k,l) = nobs_3d(i,j,k,l) + 1.

           end if

       end do


       end subroutine GrITAS_MinMax_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_3CH_Norm_ --- Normalize time mean statistics for 3CH
! 
! !INTERFACE:
!
       subroutine GrITAS_3CH_Norm_ ( nr, nrmlz,
     &                               l2d, bias_2d, rtms_2d, var_2d, xcov_2d, nobs_2d,
     &                               l3d, bias_3d, rtms_3d, var_3d, xcov_3d, nobs_3d,
     &                               ntresh )
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::      nr           ! number of distinct residual series
      logical, intent(in) ::      nrmlz        ! normalize or not
      integer, intent(in) ::      l2d          ! actual number of 2d grids
      integer, intent(in) ::      l3d          ! actual number of 3d grids
      integer, OPTIONAL, intent(in) :: ntresh  ! tresh hold to compute means

! !INPUT/OUTPUT PARAMETERS:
!
                                                      ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(im,jm,l2d,nr)    !   time means
      real, intent(inout) ::  rtms_2d(im,jm,l2d,nr)    !   root time mean square
      real, intent(inout) ::  xcov_2d(im,jm,l2d,nr)    !   holds 3ch results
      real, intent(inout) ::  nobs_2d(im,jm,l2d)       !   number of obs per grid box

                                                      ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(im,jm,km,l3d,nr) !   time means
      real, intent(inout) ::  rtms_3d(im,jm,km,l3d,nr) !   root time mean square
      real, intent(inout) ::  xcov_3d(im,jm,km,l3d,nr) !   holds 3ch result
      real, intent(inout) ::  nobs_3d(im,jm,km,l3d)    !   number of obs per grid box
!
! !OUTPUT PARAMETERS:
!
      real, intent(inout) ::  var_2d(im,jm,l2d,nr)    ! 2D root time mean square
      real, intent(inout) ::  var_3d(im,jm,km,l3d,nr) ! 3D root time mean square

!
! !DESCRIPTION:  Normalizes the time mean statistics (such as
!                dividing the accumulated O-F by the number of
!  observations), set to 'missing value' the gridboxes with no
!  data.
!
! !REVISION HISTORY: 
!
!  28Oct06  Todling    Implementation for 3CH (based on others)
! 
! !REMARKS:
!   1. rtms is returned but not needed for 3CH
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_3CH_Norm_'

      integer i, j, k, l, m, ir, ntresh_
      real    tmp(nr), n

      integer icount
      real, allocatable :: glo2dmean(:,:)
      real, allocatable :: glo3dmean(:,:,:), glo3drtms(:,:,:)

      if ( nr > 3 ) call gritas_perror (myname//' could not handle more then 2 residuals')

      ntresh_ = 0                   ! by default no effect
      if ( present(ntresh) ) then
           ntresh_ = ntresh
      end if


!     Surface (2D) grids
!     ------------------
      m=1
      do l = 1, l2d
         do j = 1, jm
            do i = 1, im

               n = nobs_2d(i,j,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = bias_2d(i,j,l,ir) / m                                ! mean(x)
                       tmp(ir)           = rtms_2d(i,j,l,ir) / m                                ! mean(x**2)
                       rtms_2d(i,j,l,ir) = sqrt(tmp(ir))
                       tmp(ir)           = abs(tmp(ir)  - bias_2d(i,j,l,ir) *bias_2d(i,j,l,ir)) ! mean(x'**2)
                    end do
               else
                    do ir = 1, nr
                       bias_2d(i,j,l,ir) = AMISS
                       rtms_2d(i,j,l,ir) = AMISS
                    end do
               end if

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr 
                       var_2d(i,j,l,ir)  = m * tmp(ir) / (m-1)
                    end do
               else
                    do ir = 1, nr 
                       var_2d(i,j,l,ir)  = AMISS
                    end do
               end if

            end do
         end do
      end do

!     Now calculate 3CH estimate (o,b,a)
!     --------------------------
      xcov_2d(:,:,:,1) = 0.5*(var_2d(:,:,:,1) + var_2d(:,:,:,2) - var_2d(:,:,:,3))
      xcov_2d(:,:,:,2) = 0.5*(var_2d(:,:,:,2) + var_2d(:,:,:,3) - var_2d(:,:,:,1))
      xcov_2d(:,:,:,3) = 0.5*(var_2d(:,:,:,3) + var_2d(:,:,:,1) - var_2d(:,:,:,2))


      if( nlevs<2 ) return

!     Upper-air (3D) grids
!     --------------------

      m=1
      do l = 1, l3d
         do k = 1, nlevs
            do j = 1, jm
               do i = 1, im
                  
               n = nobs_3d(i,j,k,l) 
               if ( nrmlz ) m=n
               if ( n .gt. 0. .AND. n .ge. ntresh_ ) then
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = bias_3d(i,j,k,l,ir) / m                        ! mean(x) and mean(y)
                       tmp(ir)             = rtms_3d(i,j,k,l,ir) / m                        ! mean(x**2) and mean(y**2)
                       rtms_3d(i,j,k,l,ir) = sqrt(tmp(ir))
                       tmp(ir)  = abs(tmp(ir)  - bias_3d(i,j,k,l,ir) *bias_3d(i,j,k,l,ir))  ! mean(x'**2) and mean(y'**2)
                    end do
               else
                    do ir = 1, nr
                       bias_3d(i,j,k,l,ir) = AMISS
                       rtms_3d(i,j,k,l,ir) = AMISS
                    end do
               endif

               if ( n .gt. 1. .AND. n .ge. ntresh_ ) then
                    var_3d(i,j,k,l,1:nr)  = m * tmp(1:nr) / (m-1)
               else
                    if (n==1) then 
                       var_3d(i,j,k,l,1:nr)  = tmp(1:nr)
                    else
                       var_3d(i,j,k,l,1:nr) = AMISS
                    endif
               end if

               end do
            end do
         end do
      end do

!     Now calculate 3CH estimate (o,b,a)
!     --------------------------
!                                         (x-y)*2             (x-z)*2              (y-z)*2
      xcov_3d(:,:,:,:,1) = 0.5*(var_3d(:,:,:,:,1) + var_3d(:,:,:,:,2) - var_3d(:,:,:,:,3))
      xcov_3d(:,:,:,:,2) = 0.5*(var_3d(:,:,:,:,3) + var_3d(:,:,:,:,1) - var_3d(:,:,:,:,2))
      xcov_3d(:,:,:,:,3) = 0.5*(var_3d(:,:,:,:,2) + var_3d(:,:,:,:,3) - var_3d(:,:,:,:,1))

      end subroutine GrITAS_3CH_Norm_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_cdf_ --- Calculate Chi-Square given P and DF(nobs)
! 
! !INTERFACE:
!
       subroutine GrITAS_cdf_ ( CDF, Prob,
     &                          nobs_2d, cdfval_2d,
     &                          nobs_3d, cdfval_3d )
!
! !USES:
!
      use m_die, only: die

      implicit NONE
!
! !INPUT PARAMETERS: 
!

      character(len=*), intent(in) :: CDF             ! cummulative dist function: t-dist or chisq-dist
      real,    intent(in) ::  prob                    ! probability
      real,    intent(in) ::  nobs_2d(:,:,:)          !   number of 2d obs for each case
      real,    intent(in) ::  nobs_3d(:,:,:,:)        !   number of 3d obs for each case

! !INPUT/OUTPUT PARAMETERS:
!
                                                      ! Surface (2D) grids:
      real, intent(inout) ::  cdfval_2d(:,:,:)        !   chi-square/t-student result
                                                      ! Upper-air (3D) grids:
      real, intent(inout) ::  cdfval_3d(:,:,:,:)      !   chi-square/t-student result

!
! !DESCRIPTION:  Calcualate T-Student or Chi-Square to allow assinging confidence
!                interval to mean or standard-deviation, respectively.
!
! !REVISION HISTORY: 
!
!  08Nov14  Todling    Initial code
!  31Jul20  Todling    Initial code
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_cdf_'

      integer ii,jj,kk,ll
      integer nx,ny,nz,nf,which,ier
      double precision df,p,q,bound,x2
   
      if( trim(cdf)/='t-dist' .and. trim(cdf)/='chisq-dist' ) then
         call die(myname_,': unsupported CDF: '//trim(cdf) //', aborting ...',99)
      endif

      which = 2
      q=real(prob,8)
      p= 1-q
      write(6,'(2a,2(2x,f9.3))') myname_, ': Confidence levels for q/p: ', q, p

!    2D case
      nx=size(cdfval_2d,1)
      ny=size(cdfval_2d,2)
      nf=size(cdfval_2d,3)
      do ll = 1, nf
         cdfval_2d(:,:,ll)=0.0
         df = 0.0
         do jj = 1, ny
            do ii = 1, nx
               if(nobs_2d(ii,jj,ll)==0) cycle
               df=df+real(nobs_2d(ii,jj,ll),8)
!              df=real(nobs_2d(ii,jj,ll)-1,8)
!              if(trim(cdf)=='t-dist') then
!                 call cdft   (which,p,q,x2,df,ier,bound)
!              else
!                 call cdfchi (which,p,q,x2,df,ier,bound)
!              endif
!              if (ier==0) then
!                  cdfval_2d(ii,jj,ll) = x2    
!              endif
            end do
            if (df < 1.0) cycle
            df = df - 1.0
            if(trim(cdf)=='t-dist') then
               call cdft   (which,p,q,x2,df,ier,bound)
            else
               call cdfchi (which,p,q,x2,df,ier,bound)
            endif
            if (ier==0) then
                cdfval_2d(:,:,ll) = x2    
            endif
         end do
      end do

!    3D case
      nx=size(cdfval_3d,1)
      ny=size(cdfval_3d,2)
      nz=size(cdfval_3d,3)
      nf=size(cdfval_3d,4)
      do ll = 1, nf
         cdfval_3d(:,:,:,ll)=0.0
         do kk = 1, nz
            df = 0.0
            do jj = 1, ny
               do ii = 1, nx
                  if(nobs_3d(ii,jj,kk,ll)==0) cycle
                  df=df+real(nobs_3d(ii,jj,kk,ll),8)
!                 df=real(nobs_3d(ii,jj,kk,ll)-1,8)
!                 if(trim(cdf)=='t-dist') then
!                    call cdft   (which,p,q,x2,df,ier,bound)
!                 else
!                    call cdfchi (which,p,q,x2,df,ier,bound)
!                 endif
!                 if (ier==0) then
!                     cdfval_3d(ii,jj,kk,ll) = x2    
!                 endif
               end do
            end do
            if (df < 1.0) cycle
            df = df - 1
            if(trim(cdf)=='t-dist') then
               call cdft   (which,p,q,x2,df,ier,bound)
            else
               call cdfchi (which,p,q,x2,df,ier,bound)
            endif
            if (ier==0) then
                cdfval_3d(:,:,kk,ll) = x2    
            endif
         end do
      end do

      end subroutine GrITAS_cdf_

      subroutine GrITAS_Strip_(A,n,iextra)
      implicit none
      real,intent(inout) :: A(:,:)
      integer, intent(out) :: n
      integer, intent(out), optional :: iextra(:) 
      
      real,allocatable,dimension(:,:):: B
      integer, allocatable, dimension(:) :: indx
      integer m,ii,jj,i1,j1

      n = size(A,1)
      m = size(A,2)
      if(n/=m) call gritas_perror ('Strip: matrix must be square')
      allocate(indx(n))
      indx=-1
      m=0
      ! find zero rows
      do jj=1,n
         if(sum(abs(A(:,jj))) > 1.e-10) then
           m=m+1
           indx(jj)=jj
         endif
      enddo
      allocate(B(m,m))
      ! strip off rows
      j1=0
      do jj=1,n
         if(indx(jj)<0) cycle
         j1=j1+1
         i1=0
         do ii=1,n
            if(indx(ii)<0) cycle
            i1=i1+1
            B(i1,j1) = A(indx(ii),indx(jj))
         enddo
      enddo
      A=0.0
      A(1:m,1:m)=B
      if (present(iextra)) then
         if(size(iextra)>=n) then
            iextra(1:n)=indx
         endif
      endif
      n=m

      deallocate(B)
      deallocate(indx)
      end subroutine GrITAS_Strip_

      subroutine scores1_ ( expid, nymd, nhms, hour,
     &                      lat, lon, lev, kx, kt, del, nobs, nr,
     &                      kt_list, var_list, kxkt_table, p1, p2,
     &                      ob_flag, l2d, l3d )

      use m_obs_stats, only: obs_stats_init
      use m_obs_stats, only: obs_stats_collect
      use m_obs_stats, only: obs_stats_write
      use m_obs_stats, only: obs_stats_final
      use m_obs_stats, only: obs_stats

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX

      use m_die, only: die
      implicit none
!
! !INPUT PARAMETERS: 
!

      character(len=*), intent(in) :: expid      ! experiment ID
      integer, intent(in) ::       nymd          ! date: YYYYMMDD
      integer, intent(in) ::       nhms          ! date: HHMMSS
      integer, intent(in) ::       hour          ! forecast hour 0=bkg, 24, etc
      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lat(nobs)     ! latitude (degrees)
      real,    intent(in) ::       lon(nobs)     ! longitude (degrees)
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index
      integer, intent(in) ::  ob_flag(nobs)      ! data qc flag

      
      integer, intent(in) ::  kt_list(:)         ! list of kt from RC
      character(len=*), intent(in) :: var_list(:)  ! list of variables
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      real,    intent(in) ::           p1(nlevs)
      real,    intent(in) ::           p2(nlevs)     ! pressure intervals for gridding

      integer,    intent(in) ::         l2d          ! actual number of 2d grids
      integer,    intent(in) ::         l3d          ! actual number of 3d grids


!     local variables

      character(len=*), parameter :: myname_ = myname//"*scores1_"

                                               ! Surface (2D) grids:
      real, allocatable ::  bias_2d(:,:,:,:)   !   time means
      real, allocatable ::  rtms_2d(:,:,:,:)   !   root time mean square
      real, allocatable ::  stdv_2d(:,:,:,:)   !   standard deviations
      real, allocatable ::  xcov_2d(:,:,:,:)   !   cross-covariance
      real, allocatable ::  nobs_2d(:,:,:)     !   number of obs per grid box

                                               ! Upper-air (3D) grids:
      real, allocatable ::  bias_3d(:,:,:,:,:) !   time means
      real, allocatable ::  rtms_3d(:,:,:,:,:) !   root time mean square
      real, allocatable ::  stdv_3d(:,:,:,:,:) !   standard deviations
      real, allocatable ::  xcov_3d(:,:,:,:,:) !   cross-covariance
      real, allocatable ::  nobs_3d(:,:,:,:)   !   number of obs per grid box

      character(len=80) :: fname
      real(4), allocatable :: collect_stats(:,:,:,:)
      integer, allocatable :: this_ob_flag(:)
      integer ns, ier
      logical, allocatable :: sfcvars(:)
      type(obs_stats) stats

      character(len=1), parameter :: fakevars(1)=(/'h'/)

      allocate(sfcvars(size(kt_list)))
      sfcvars = .false.
      where(kt_list<0) sfcvars=.true.

      allocate(this_ob_flag(nobs))
     
      allocate ( bias_2d(im,jm,l2d,nr),  rtms_2d(im,jm,l2d,nr),
     &           stdv_2d(im,jm,l2d,nr),
     &           xcov_2d(im,jm,l2d,nr),  nobs_2d(im,jm,l2d),
     &           stat=ier )
        if(ier/=0) call die (myname_,'Error in Alloc(2d)')

      allocate ( bias_3d(im,jm,km,l3d,nr), rtms_3d(im,jm,km,l3d,nr),
     &           stdv_3d(im,jm,km,l3d,nr),
     &           xcov_3d(im,jm,km,l3d,nr), nobs_3d(im,jm,km,l3d),
     &           stat=ier )
        if(ier/=0) call die (myname_,'Error in Alloc(3d)')

      nobs_2d = 0.0
      nobs_3d = 0.0

      bias_2d = 0.0; bias_3d = 0.0
      rtms_2d = 0.0; rtms_3d = 0.0
      stdv_2d = 0.0; stdv_3d = 0.0
      xcov_2d = 0.0; xcov_3d = 0.0

      this_ob_flag = ob_flag

!     Accumulate statistics for this date/time only
!     ---------------------------------------------
      call GrITAS_Accum_ ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                     kxkt_table, p1, p2,
     &                     this_ob_flag,
     &                     l2d, bias_2d, rtms_2d, xcov_2d, nobs_2d,
     &                     l3d, bias_3d, rtms_3d, xcov_3d, nobs_3d )

!     Normalize statistics for this date/time only
!     --------------------------------------------
      call GrITAS_Norm_ ( nr, .true.,
     &                    l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                    l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d )

!     Gather scorecard-like summary of stats for this date/time
!     ---------------------------------------------------------
      allocate(collect_stats(mregs,km,l3d+l2d,nevals),
     &           stat=ier )
      if(ier/=0) call die (myname_,': Error in Alloc(3d)')

      call scores2_ ( collect_stats,
     &                bias_2d, rtms_2d, stdv_2d, nobs_2d,
     &                bias_3d, rtms_3d, stdv_3d, nobs_3d   )

!     Write out scorecard for this date/time
!     --------------------------------------
      ns = mregs*(km*(l3d+size(fakevars))+l2d)*nevals
      call obs_stats_init ( ns, stats, expid, 'residual', 'gmao', 'gmao' )
      call obs_stats_collect (nymd,nhms,hour,var_list,sfcvars,type_stats,(/'pl'/),real(plevs,4),
     &                        regnames,real(regbounds,4),collect_stats,stats,fakevars=fakevars)
      write(fname,'(a,i3.3,a,i8.8,a,i2.2,a)') 'scorecard.', hour, 'h.', nymd, '_', nhms/10000, 'z.csv'
      call obs_stats_write ( fname, stats )
      call obs_stats_final ( stats )

      deallocate(collect_stats)
      deallocate ( bias_3d,  rtms_3d, stdv_3d,
     &             xcov_3d,  nobs_3d )
      deallocate ( bias_2d,  rtms_2d, stdv_2d,
     &             xcov_2d,  nobs_2d )
      deallocate(this_ob_flag)
      deallocate(sfcvars)
     
      end subroutine scores1_

      subroutine scores2_ ( collect_stats,
     &                      bias_2d, rtms_2d, stdv_2d, nobs_2d,
     &                      bias_3d, rtms_3d, stdv_3d, nobs_3d   )
      use m_die, only: die
      implicit none

      real(4), intent(inout) :: collect_stats(:,:,:,:)
                                               ! Surface (2D) grids:
      real, intent(inout) ::  bias_2d(:,:,:,:)   !   time means
      real, intent(inout) ::  rtms_2d(:,:,:,:)   !   root time mean square
      real, intent(inout) ::  stdv_2d(:,:,:,:)   !   standard deviations
      real, intent(inout) ::  nobs_2d(:,:,:)     !   number of obs per grid box

                                                ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_3d(:,:,:,:,:) !   time means
      real, intent(inout) ::  rtms_3d(:,:,:,:,:) !   root time mean square
      real, intent(inout) ::  stdv_3d(:,:,:,:,:) !   standard deviations
      real, intent(inout) ::  nobs_3d(:,:,:,:)   !   number of obs per grid box

      integer im,jm,km,nv2d,nv3d,nr,nnr
      integer nregs,l3d,nevals
      integer nt,ne,mr,nv,kk,lm
      real sig_omf,sig_oma

      nt   = 1
      nv2d = size(bias_2d,3)
      nv3d = size(bias_3d,4)
      nr   = size(bias_3d,5)
      km   = size(bias_3d,3)
      jm   = size(bias_3d,2)
      im   = size(bias_3d,1)
      if (nr/=1) then
        call die('GrITAS_scores_',': incorrec number of residuals, aborting ...',98)
      endif

      if(size(collect_stats,2) /=km) then
        call die('GrITAS_scores_',': inconsitent level count, aborting ...',99)
      endif
      if(size(collect_stats,3) /=nv3d+nv2d) then
        call die('GrITAS_scores_',': inconsitent variable count, aborting ...',99)
      endif
      nregs  = size(collect_stats,1)
      nevals = size(collect_stats,4)

! for test: global stdv omf only
         do ne=1,nevals
            if(trim(type_stats(ne))=='mean' .or. 
     &         trim(type_stats(ne))=='rms'  .or.
     &         trim(type_stats(ne))=='stdv') nnr = 1
            do mr=1,nregs
               ! mean
               if(trim(type_stats(ne))=='mean') then
                 nv = 0
                 do lm=1,nv2d
                    nv=nv+1
                    collect_stats(mr,1,nv,ne) = ave_nomiss(bias_2d(:,:,lm,nnr),nobs_2d(:,:,lm),lonindx(:,nt),regindx(:,mr))
                 enddo
                 do lm=1,nv3d
                    nv=nv+1
                    do kk=1,km
                       collect_stats(mr,kk,nv,ne) = ave_nomiss(bias_3d(:,:,kk,lm,nnr),nobs_3d(:,:,kk,lm),
     &                                                         lonindx(:,nt),regindx(:,mr))
                    enddo
                 enddo
               endif
               ! rms
               if(trim(type_stats(ne))=='rms') then
                 nv = 0
                 do lm=1,nv2d
                    nv=nv+1
                    collect_stats(mr,1,nv,ne) = ave_nomiss(rtms_2d(:,:,lm,nnr),nobs_2d(:,:,lm),
     &                                                     lonindx(:,nt),regindx(:,mr))
                 enddo
                 do lm=1,nv3d
                    nv=nv+1
                    do kk=1,km
                       collect_stats(mr,kk,nv,ne) = ave_nomiss(rtms_3d(:,:,kk,lm,nnr),nobs_3d(:,:,kk,lm),
     &                                                         lonindx(:,nt),regindx(:,mr))
                    enddo
                 enddo
               endif
               ! stdv
               if(trim(type_stats(ne))=='stdv') then
                 nv = 0
                 do lm=1,nv2d
                    nv=nv+1
                    collect_stats(mr,1,nv,ne) = ave_nomiss(stdv_2d(:,:,lm,nnr),nobs_2d(:,:,lm),lonindx(:,nt),regindx(:,mr))
                 enddo
                 do lm=1,nv3d
                    nv=nv+1
                    do kk=1,km
                       collect_stats(mr,kk,nv,ne) = ave_nomiss(stdv_3d(:,:,kk,lm,nnr),nobs_3d(:,:,kk,lm),
     &                                                         lonindx(:,nt),regindx(:,mr))
                    enddo
                 enddo
               endif
               ! percent contribution
!              if(trim(type_stats(ne))=='contrib') then
!                do lm=1,nv2d
!                   nv=nv+1
!                   sig_omf = ave_nomiss(stdv_2d(:,:,lm,1),nobs_2d(:,:,lm),lonindx(:,nt),regindx(:,mr))
!                   sig_oma = ave_nomiss(stdv_2d(:,:,lm,2),nobs_2d(:,:,lm),lonindx(:,nt),regindx(:,mr))
!                   if (abs(sig_omf-AMISS)>0.01.and.abs(sig_oma-AMISS)>0.01) then
!                      collect_stats(mr,1,nv,ne) = 100 * ( 1.0 - sig_oma/sig_omf )
!                   else
!                      collect_stats(mr,1,nv,ne) = AMISS
!                   endif
!                enddo
!                do lm=1,nv3d
!                   nv=nv+1
!                   do kk=1,km
!                      sig_omf = ave_nomiss(stdv_3d(:,:,kk,lm,1),nobs_3d(:,:,kk,lm),lonindx(:,nt),regindx(:,mr))
!                      sig_oma = ave_nomiss(stdv_3d(:,:,kk,lm,2),nobs_3d(:,:,kk,lm),lonindx(:,nt),regindx(:,mr))
!                      if (abs(sig_omf-AMISS)>0.01.and.abs(sig_oma-AMISS)>0.01) then
!                         collect_stats(mr,kk,nv,ne) = 100 * ( 1.0 - sig_oma/sig_omf )
!                      else
!                         collect_stats(mr,kk,nv,ne) = AMISS
!                      endif
!                   enddo
!                enddo
!              endif
            enddo
         enddo

      end subroutine scores2_

      real function ave_nomiss(x,nobs,xregion,yregion)
      real, intent(in) :: x(:,:)
      real, intent(in) :: nobs(:,:)
      integer, intent(in) :: xregion(2),yregion(2)
      integer im,jm,ii,jj
      real counter,threshold
      real add
      threshold = 0.0
      im=size(x,1)  
      jm=size(x,2)  
      add=-AMISS
      counter=0.0
      do jj=yregion(1),yregion(2)
         do ii=xregion(1),xregion(2)
            if (nobs(ii,jj)>=threshold.and.abs(x(ii,jj)-AMISS)>0.01) then
               if(add==-AMISS) then
                  add=x(ii,jj) 
               else
                  add=add+x(ii,jj) 
               endif
               counter=counter+1.0
            endif
         enddo
      enddo
      if(add==-AMISS) then
         ave_nomiss = AMISS
      else
         ave_nomiss = add/counter
      endif
      end function ave_nomiss

      subroutine hydro_heights_
 
      use m_die, only: die
      implicit none

      real, parameter ::  grav   = 9.80665 ! m/s**2
      real, parameter ::  Mair   = 0.02896 ! kg/mol
      real, parameter ::  tstd   = 288.15  ! K
      real, parameter ::  pstd   = 1013.25 ! mb (mean sea level)
      real, parameter ::  rgas   = 8.3143  ! (N*m)/(mol*K)
   
      integer k
      real const

      allocate(hlevs(size(plevs)))
   
      const = rgas * tstd / (mair*grav)
      write(6,'(a)') ' Lev Index   Pressure(mb)  Height(km)'
      do k = 1, size(plevs)
         hlevs(k) = const * log(pstd/plevs(k))
         write(6,'(i5,1x,2(f11.3))') k, plevs(k), hlevs(k)/1000.
      enddo

      end subroutine hydro_heights_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GrITAS_XCov_Prs_Accum_ --- Accumulate O-F on del grids
! 
! !INTERFACE:
!
      subroutine GrITAS_XCov_Prs_Accum_ ( lev, kx, kt, ks, del, nobs, nr,
     &                                    kxkt_table, p1, p2,
     &                                    ob_flag,
     &                                    l3d, bias_lv, xcov_lv, nobs_lv, nprof, tch )

! !USES:
!

      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : KXMAX
      use m_die, only: die
      implicit NONE
!
! !INPUT PARAMETERS: 
!

      integer, intent(in) ::       nobs          ! actual number of observations
      integer, intent(in) ::       nr            ! number of distinct residual series
      real,    intent(in) ::       lev(nobs)     ! level (hPa)
      real,    intent(in) ::       del(nobs,nr)  ! residual
      integer, intent(in) ::       kx(nobs)      ! data source index
      integer, intent(in) ::       kt(nobs)      ! data type   index
      integer, intent(in) ::       ks(nobs)      ! sound index

      
      integer, intent(in) ::        kxkt_table(KXMAX,KTMAX)  ! kx-kt table

      real,    intent(in) ::           p1(nlevs)
      real,    intent(in) ::           p2(nlevs)     ! pressure intervals for gridding

      integer,    intent(in) ::         l3d          ! actual number of 3d grids
      logical,    intent(in) ::         tch          ! TCH or Desroziers et al

!
! !INPUT/OUTPUT PARAMETERS:
!
      integer,    intent(inout) :: ob_flag(nobs)     ! On input, ob is rejected if
                                                     !    ob_flag is other than 0
                                                     ! On output, ob_flag = 0 if ob
                                                     !    was used to accumulate the
                                                     !    bias/RMS

                                                     ! Upper-air (3D) grids:
      real, intent(inout) ::  bias_lv(:,:,:)         !  time means
      real, intent(inout) ::  xcov_lv(:,:,:,:)       !   root time mean square
      real, intent(inout) ::  nobs_lv(:,:)           !   number of obs per grid box

      integer,intent(out) :: nprof

!
! !DESCRIPTION: Given the O-F for a given synoptic time, this routine
!               accumulates the bias and the RMS at the nearest gridpoint.
!
! !REMARKS:
!
!   1.  xrtms is redundant when nr=1, i.e., xrtms=rtms
!   2.  when nr=2, assumes residuals are from same experiment, that is,
!       no need to distinguish between nobs
!
! !REVISION HISTORY: 
!
!  04Jun22  Todling    Calculate pressure-level covariances
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'GrITAS_XCov_Prs_Accum_'

      integer n, m, nn, mm, k1, k2, ks1, kk, nx, l, ir, ncnt, np, ip
      logical surface_var, reject_ob
      integer, allocatable, dimension(:) :: nindx,kindx
      integer, allocatable, dimension(:) :: uks
      integer, allocatable, dimension(:) :: nfound
      real, allocatable, dimension(:) :: rlev
      real, allocatable, dimension(:,:) :: rdel

      if ( nr > 3 ) then
          call gritas_perror (myname_ //': improper dim settings')
      endif
      if (tch) then
        if (nr/=3) then
          call gritas_perror (myname_ //': TCH setting inconsistent w/ dims')
        endif
      else
        if (nr/=2) then
          call gritas_perror (myname_ //': Desroziers setting inconsistent w/ dims')
        endif
      endif

      nprof=1
      ks1=ks(1)
      do n = 1, nobs
         if(ks(n)/=ks1) then
           ks1=ks(n)
           nprof=nprof+1
         endif
      enddo
      print *, 'found this many profiles: ', nprof
      allocate(uks(nprof),nfound(nprof))
      nprof=1
      ks1 = ks(1)
      uks(1)=ks1
      nfound(1) = 0
      do n = 1, nobs
         nfound(nprof) = nfound(nprof)+1
         if(ks(n)/=ks1) then
           ks1=ks(n)
           nprof=nprof+1
           if(nprof>size(uks)) call die (myname,'crap',99)
           uks(nprof) = ks1 
           nfound(nprof) = 0
         endif
      enddo
!     do n=1,nprof
!       print *, uks(n), nfound(n)
!     enddo
      if(sum(nfound) /= nobs) then
         print *, sum(nfound), nobs
         call die(myname,'incorrect sume',99)
      endif

      do np = 1, nprof
         allocate(rlev(nfound(np)),rdel(nfound(np),nr))
         ip=0
         do n=1,nobs 
            if(ks(n)/=uks(np)) cycle
            ip=ip+1
            rlev(ip)   = lev(n)
            rdel(ip,:) = del(n,:)
         enddo
         ncnt=0
         do n=1,ip
            do kk = 1, nlevs
               if ( rlev(n) .le. p1(kk) .AND. 
     &              rlev(n) .gt. p2(kk) ) then
                    ncnt = ncnt + 1
                    go to 21
               end if
            end do
            call gritas_perror ('Accum: could not find level(1)')
 21         continue
         enddo

         allocate(nindx(ncnt),kindx(ncnt))
         nindx=-1
         kindx=-1
         ncnt = 0

!        Find vertical level
!        -------------------
         do n=1,ip
            do kk = 1, nlevs
               if ( rlev(n) .le. p1(kk) .AND. 
     &              rlev(n) .gt. p2(kk) ) then
                    ncnt = ncnt + 1
                    kindx(ncnt) = kk
                    nindx(ncnt) = n
                    go to 31
               end if
            end do
            call gritas_perror ('Accum: could not find level(2)')
 31         continue
         enddo
         if(ncnt/=ip) 
     &      call gritas_perror ('Accum: counts not matching')

!        Loop over observations
!        ----------------------
         l = 1
         do n = 1, ip
            if(nindx(n)<0) cycle
            k1=kindx(n)
            nn=nindx(n)

!           Accumulate
!           ----------              
            nobs_lv(k1,l)   = nobs_lv(k1,l)   + 1.
            bias_lv(k1,l,:) = bias_lv(k1,l,:) + rdel(nn,:)
            do m = 1, ip
               if(nindx(m)<0) cycle
               k2=kindx(m)
               mm=nindx(m)
               if (tch) then ! 3CH
                   xcov_lv(k1,k2,l,:) = xcov_lv(k1,k2,l,:) + rdel(nn,:)*rdel(mm,:)
               else          ! for Desroziers et al (low-level symmetrized)
!                  xcov_lv(k1,k2,l,1) = xcov_lv(k1,k2,l,1) + rdel(nn,1)*rdel(mm,2)
                   xcov_lv(k1,k2,l,1) = xcov_lv(k1,k2,l,1) + 0.5*(rdel(nn,1)*rdel(mm,2)+rdel(nn,2)*rdel(mm,1))
               endif
            enddo
         enddo

         deallocate(nindx,kindx)
         deallocate(rlev,rdel)
      enddo
      deallocate(uks,nfound)


      end subroutine GrITAS_XCov_Prs_Accum_


      end module m_gritas_grids
