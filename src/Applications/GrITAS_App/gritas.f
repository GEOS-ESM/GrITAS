!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOI
!
! !TITLE: An Application for Gridding Innovation\\ Time Averaged Statistics (GrITAS v2.00)
!
! !AUTHORS: Arlindo da Silva and Guang Ping Lou
!
! !AFFILIATION: Data Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: September 1, 1997 (Revised June 9, 2004)
!
! !INTRODUCTION: System Overview
!
!  GrITAS is a FORTRAN 77/90 program which reads ODS or GEOS-1 style
!  del files and produces global grids of O-F (observation minus forecast),
!  O-A (observation minus analysis) or  O (observed value) with time mean,
!  standard deviation and number of observations in each grid box.
!  Grids are produced for user specified combination of data types
!  (e.g., u, v-winds, heights) and data sources (e.g., radiosondes,
!  TOVS retrievals, ships, etc.) 
!
!
!EOI
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS()
! 
! !DESCRIPTION: Driver for the {\em Gridding of Innovation Time Average
!               Statistics} system (GrITAS).
!
! !INTERFACE:
!
       Program GrITAS

! !USES:
       Use m_MergeSorts

       Use m_odsmeta, only : KTMAX
       Use m_odsmeta, only : KXMAX
       Use m_odsmeta, only : X_PRE_BAD  ! GSI  passive data
       Use m_odsmeta, only : X_PASSIVE  ! PSAS passive data

       Use m_gritas_grids, only : im, jm, km
       Use m_gritas_grids, only : GrITAS_Pint
       Use m_gritas_grids, only : GrITAS_Accum
       Use m_gritas_grids, only : GrITAS_MinMax
       Use m_gritas_grids, only : GrITAS_Norm
       Use m_gritas_grids, only : GrITAS_XCov_Accum
       Use m_gritas_grids, only : GrITAS_XCov_Prs_Accum
       Use m_gritas_grids, only : GrITAS_XCov_Norm
       Use m_gritas_grids, only : GrITAS_3CH_Norm
       Use m_gritas_grids, only : GrITAS_Strip
       Use m_gritas_grids, only : GrITAS_CDF
       Use m_gritas_grids, only : GrITAS_Scores
       Use m_gritas_grids, only : plevs
       Use m_gritas_grids, only : nlevs
       Use m_gritas_grids, only : nchan
       Use m_gritas_grids, only : achan
       Use m_gritas_grids, only : mregs
       Use m_gritas_grids, only : AMISS
 
       Use m_gritas_masks, only: gritas_maskout

       Use m_gfio_output, only: GFIO_Output

       Use m_ods

!      Use m_obs, only : Reducer

       Use m_StrTemplate

       Use m_stdio, only : stdout
       Use m_die

! !REMARKS:
!
!   1) In case of calculating statistics of (specified) sigo and obias the standard deviation 
!      output is replaced by root-mean-square; that is to say that position t=1 in the grads
!      file will have the mean as it normally has, and t=2 will have the RMS instead of the STDV.
!
!   2) In case of Desroziers et al. diagnostics, the mean and standard deviations output 
!      are replaced by cross-variance and the cross standard deviations, respectively.
!
! !REVISION HISTORY: 
!
!  01Sep97  da Silva   Initial design & code.
!  17Oct97  G.P. Lou   Modifications. (RT: What modifications?)
!  15Mar98  da Silva   Minor touch ups (comments only)
!  11Dec98  da Silva   Included ODS input along with O-A/obs option.
!  15Dec98  da Silva   Added GFIO output.
!  23dec00  da Silva   Made part of fvDAS.
!  09Jun04  Todling    - Using m_ods
!                      - Resolution independent executable
!                      - Changed action when no-obs found in syn time
!  24Jun04  Todling    Fixed time frequency of output (RUC compatible).
!  30Jul04  Todling    Added reducer capability
!  12oct05  da Silva   Clean ODS each time
!  07Nov06  Todling    - Implemented Desroziers et al. diagnostics
!                      - Implemented Jo option.
!  23Jan08  Todling    Added -obsimp
!  12Sep08  Todling    Added -nonorm to facilitate obs impact calculation
!  07Jul09  Redder     Added ob_flag as an argument in the call to the
!                      routine GrITAS_Accum
!  20Jul09  Redder     Added the code to generate the descriptions for
!                      for each entry in the kt list.
!  21Jul09  Redder     Fixed bug in call to the routine, GrITAS_Init.
!  04Aug11  Ravi       Fixed a bug in Call to the routine GFIO_Output.
!                      (added var_descr)
!  10May13  Todling    - Add -fobloc to calc fcst at ob-location
!                      - Add opt to read extra ODS file and calc ens spread 
!  10Apr14  Todling    Add opt -xcov; combined w/ say esigo it allows for
!                      estimation of full R error cov matrix
!  12Apr14  Todling    Add opt -mask to allow specification of mask 
!                      to exclude data over certain regions
!  18Oct14  Todling    Multiple fixes related to ens spread and 
!                      "prescribed" innov covariance calculation
!  31Oct14  Todling    Add opt rtms to replace stddev w/ rtms fields
!
!EOP
!-------------------------------------------------------------------------
!BOC
      Implicit NONE

      character(len=*), parameter :: myname = 'GrITAS'

      include 'gritas.h'


!     User defined vertical levels
!     ----------------------------
      real, allocatable :: p1(:), p2(:)  ! pressure intervals for gridding


!     Workspace for Holding observations for each synoptic time
!     Note: lat/lon/etc cannot be pointers because ods is r8 vs lat_r4
!     ----------------------------------------------------------------
      integer                  nobs      ! actual number of obs for this synoptic time
      integer                  nobs_good ! no. of observations after reducer operation
      type     ( ods_vect )    ods
      type     ( ods_vect )    mods
      integer, allocatable ::  kx (:)
      integer, allocatable ::  kt (:)
      integer, allocatable ::  ks (:)
      real,    allocatable ::  lat(:)
      real,    allocatable ::  lon(:)
      real,    allocatable ::  lev(:)
      real,    allocatable ::  qc (:)          ! not used, but needed
      real,    allocatable ::  del(:,:)
      integer, allocatable ::  is(:)           ! ordering index
      integer, allocatable ::  ob_flag (:)     ! = 0 if ob was used to accumulate
                                               !   bias/RMS, = nonzero otherwise

!     Workspace for holding the surface gridded fields
!     ------------------------------------------------
      real, allocatable ::  bias_2d(:,:,:,:)   ! first  time mean
      real, allocatable ::  rtms_2d(:,:,:,:)   ! root time mean square
      real, allocatable ::  stdv_2d(:,:,:,:)   ! standard deviation
      real, allocatable ::  xcov_2d(:,:,:,:)   ! x-rms and x-variance
      real, allocatable ::  nobs_2d(:,:,:)     ! number of obs per grid box
      real, allocatable ::  chsq_2d(:,:,:,:)   ! chi-square for stddev confidence
      integer list_2d(l2d_max)                 ! corresponding list item

!      equivalence ( rtms_2d(1,1,1), stdv_2d(1,1,1) )


!     Workspace for holding the upper-air gridded fields
!     --------------------------------------------------
      real, allocatable ::  bias_3d(:,:,:,:,:) ! first  time mean
      real, allocatable ::  rtms_3d(:,:,:,:,:) ! root time mean square
      real, allocatable ::  stdv_3d(:,:,:,:,:) ! standard deviation
      real, allocatable ::  xcov_3d(:,:,:,:,:) ! x-standard deviation
      real, allocatable ::  nobs_3d(:,:,:,:)   ! number of obs per grid box
      real, allocatable ::  chsq_3d(:,:,:,:,:) ! chi-square for stddev confidence
      integer list_3d(l3d_max)                 ! corresponding list item

!     Workspace for holding (channel) x-covariances
!     ---------------------------------------------
      real, allocatable ::  stdv_lv(:)       ! x-standard deviation
      real, allocatable ::  bias_lv(:,:,:)   ! first  time mean
      real, allocatable ::  xcov_lv(:,:,:,:) ! x-covariance
      real, allocatable ::  nobs_lv(:,:)     ! number of obs for each case

!     Collected statistics for scorecard
!     ----------------------------------
!     real(4), allocatable :: levels(:)

!     Data structure for user selection of data type/sources
!     ------------------------------------------------------
      integer       listsz               ! size of the list (< lmax)
      integer       kt_list(lmax)        ! each grid has one data type (kt)  
      integer       kx_list(kxmax,lmax)  ! but several data sources (kx)
      character(len=11)  var_list(lmax)       ! variable name for output 
                                         ! (e.g., uraob)
      character(len=256) var_descr(lmax)      ! ... and its description.
!     Note: Negative kt values will be treated as surface variables

!     Data type/data source table. Gives the grid number for each (kx,kt)
!     This table is derived from the user lists above. Notice that 
!     surface grids are associated with kt<4.
!     -------------------------------------------------------------------
      integer       kxkt_table(kxmax,ktmax)

      integer         l2d                      ! actual number of 2d grids
      integer         l3d                      ! actual number of 3d grids

!     Input File names
!     ----------------
      integer          ndel                      ! actual number of del-files
      character*(255)  del_ifname(NDEL_MAX)      ! O-F (del) data files
      character*(255)  odsmean_fntmpl            ! template of mean of ods files (for EnADAS only)

!     Output File names
!     -----------------
      character*(255)  grid_ofname               ! gridded O-F (del) file


!     Local variables
!     ---------------
      integer idel, isyn, ksyn, ioptn, iopt, nstat, nr, i, j, ir, nrout
      integer ii,jj,kk,ll,k1,k2,iii,jjj
      integer ns,ntimes
      integer t1, t2
      integer mdim
      integer info, nprof, nprof_total
      integer nymd, nhms, timinc, nsyn, ier, rc
      integer,pointer :: fid(:)
      integer nymd_ens, nhms_ens, nredundant
      integer nymd_min, nymd_max, nymd_out
      integer scale_res ! <=0 do not scale
                        ! 1 = scale by sigo
                        ! 2 = scale by o
                        ! 3 = scale by f
      character(len=30) expid
      character(len=8) ftype, otype
      character(len=10) dtype
      character(len=255) RC_red, fname, lb_tallyfile, covname, lwimask
      logical DO_red, first
      logical ncf, split_tgroups, meanres
      logical iampassive, iwantpassive, nrmlz
      logical xcov, self, pxcov
      logical timeselect
      logical lonselect
      logical lminmax
      logical rtms
      logical tch
      integer scorecard
      integer islot
      integer symmetrize
      real lona, lonb
      real v1, v2
      real confidence
      logical raobspd
      real,allocatable,dimension(:,:) :: aux

      nymd_min =  99999999
      nymd_max = -99999999
      nymd_out = -1         !  date of output
      nymd_ens = -1         !  redundant date (as from ensemble members)
      nhms_ens = -1         !  redundant date (as from ensemble members)

      nprof_total = 0
      lwimask='NONE'
      first = .true.
      nrmlz = .true.
      t1  =  999            ! DEFAULT: all times
      t2  = -999
      lona  = -999          ! DEFAULT: all longitudes
      lonb  =  999
      scorecard = -99       ! hours
      expid = 'NULL'
      islot = 1
      symmetrize = .true.

!...........................................................................


!     User interface
!     --------------
      call GrITAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                   ftype, dtype, otype, ioptn, ncf, nr,
     &                   split_tgroups, lb_tallyfile,
     &                   nymd_out, DO_red, RC_red,
     &                   ndel, del_ifname, grid_ofname, odsmean_fntmpl,
     &                   kt_list, kx_list, var_list, var_descr,
     &                   listsz, nrmlz, xcov, pxcov, self, rtms, t1, t2, lwimask, confidence, 
     &                   scale_res, lona, lonb, lminmax, expid, scorecard, tch, symmetrize,
     &                   raobspd )

      iopt = abs(ioptn)

      if (confidence<0.0) then
          nstat = 3     ! writes out mean, std, and nobs count
      else
          nstat = 6     ! writes out mean, std, nobs count, stdv_chisqr
      endif
      allocate(fid(1))
      fid=-1

!     Create kx-kt table from user defined list
!     -----------------------------------------
      call GrITAS_KxKt ( KTMAX, KXMAX, L2D_MAX, L3D_MAX,
     &                   listsz, kt_list, kx_list, 
     &                   kxkt_table, list_2d, list_3d, l2d, l3d )

!     Allocate space for working arrays
!     ---------------------------------
      if( xcov .or. pxcov ) then

         if(nchan>0) then
           mdim=nchan
         else
           mdim=km
         endif
         allocate ( bias_lv(mdim,l3d,nr), stdv_lv(mdim),
     &              xcov_lv(mdim,mdim,l3d,nr), nobs_lv(mdim,l3d),
     &              stat=ier )
           if(ier/=0) call die (myname,'Error in Alloc(lv)')

!        Initialize grids
!        ----------------
         bias_lv = 0.0
         xcov_lv = 0.0
         nobs_lv = 0.0

      else

         allocate ( bias_2d(im,jm,l2d,nr),  rtms_2d(im,jm,l2d,nr),
     &              stdv_2d(im,jm,l2d,nr),
     &              xcov_2d(im,jm,l2d,nr),  nobs_2d(im,jm,l2d),
     &              chsq_2d(im,jm,l2d,3),
     &              stat=ier )
           if(ier/=0) call die (myname,'Error in Alloc(2d)')

         allocate ( bias_3d(im,jm,km,l3d,nr), rtms_3d(im,jm,km,l3d,nr),
     &              stdv_3d(im,jm,km,l3d,nr),
     &              xcov_3d(im,jm,km,l3d,nr), nobs_3d(im,jm,km,l3d),
     &              chsq_3d(im,jm,km,l3d,3),
     &              stat=ier )
           if(ier/=0) call die (myname,'Error in Alloc(3d)')

!        Initialize grids
!        ----------------
         bias_2d(1:im,1:jm,     1:l2d,1:nr) = 0.0
         rtms_2d(1:im,1:jm,     1:l2d,1:nr) = 0.0
         stdv_2d(1:im,1:jm,     1:l2d,1:nr) = 0.0
         xcov_2d(1:im,1:jm,     1:l2d,1:nr) = 0.0
         nobs_2d(1:im,1:jm,     1:l2d)      = 0.0
         chsq_2d(1:im,1:jm,     1:l2d,1:3)  = 0.0
         bias_3d(1:im,1:jm,1:km,1:l3d,1:nr) = 0.0
         rtms_3d(1:im,1:jm,1:km,1:l3d,1:nr) = 0.0
         stdv_3d(1:im,1:jm,1:km,1:l3d,1:nr) = 0.0
         xcov_3d(1:im,1:jm,1:km,1:l3d,1:nr) = 0.0
         nobs_3d(1:im,1:jm,1:km,1:l3d)      = 0.0
         chsq_3d(1:im,1:jm,1:km,1:l3d,1:3)  = 0.0
         if (lminmax) then
            bias_2d= AMISS; stdv_2d=-AMISS
            bias_3d= AMISS; stdv_3d=-AMISS
         endif

      endif

      allocate ( p1(km), p2(km),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(pressure)')


!     Set pressure Intervals for gridding
!     -----------------------------------
      call GrITAS_Pint ( p1, p2 )


!     Figure out range of synoptic time loop
!     --------------------------------------
      nredundant = 1
      ksyn = 32767
      if(ncf) ksyn = 1

!     Loop over del files
!     -------------------
      do idel = 1, ndel

!         Loop over synoptic time on this del file...
!         -------------------------------------------
          do 10 isyn = 1, ksyn     
             meanres = .false.


!            Read del file for this synoptic time
!            ------------------------------------
             if ( trim(ftype) .eq. 'del' ) then

                allocate ( lat(NOBSMAX), lon(NOBSMAX), lev(NOBSMAX),
     &                     qc (NOBSMAX), del(NOBSMAX,nr),
     &                     kx(NOBSMAX),  kt(NOBSMAX), ks(NOBSMAX),
     &                     ob_flag(NOBSMAX),
     &                     stat=ier )
                   if(ier/=0) call die (myname,
     &                                ' Obs Arrays Error, alloc(del)')

                call Read_Del ( trim(del_ifname(idel)), NOBSMAX, ioptn,
     &                          lat, lon, lev, del(1,1), kx, kt, qc, nobs,
     &                          nymd, nhms, ier )
                ob_flag ( : nobs ) = 0
             else

                nymd = -1           ! get data for next synoptic time on file
                nhms =  0
                call ODSNxTime ( trim(del_ifname(idel)), nymd, nhms )
                if ( nymd .eq. -1 ) then
!                    print *, 'End-Of-File'
                     exit     ! skip to next file
                end if
                if ( first ) then
                     first = .false.
                else
                   call ods_clean ( ods, ier )
                end if
                call ODS_Get ( trim(del_ifname(idel)), nymd, nhms, ftype, ods, ier, ncf=ncf )

                if ( ier .gt. 0 ) then
                     print *, 'ODS_Get error: ier = ', ier
                     exit     ! skip to next file
                else
                     write(stdout,'(3a,i8,a,i6.6)') 
     &                    'Read ', trim(del_ifname(idel)), ' on ', nymd, ' at ', nhms
                end if

!               If so, read ODS file with mean residuals
!               ----------------------------------------
                if (trim(ODSmean_fntmpl) /= 'NONE' ) then
                    call strTemplate ( fname, ODSmean_fntmpl, 'GRADS', 'NULL', nymd, nhms, ier )
                    call ODS_Get ( trim(fname), nymd, nhms, ftype, mods, ier, ncf=ncf )
                    if ( ier .gt. 0 ) then
                         print *, 'Failed to read ODS with mean residuals:', trim(fname)
                         print *, 'Aborting, ODS_Get error: ier = ', ier
                         call exit(1)
                    else
                         write(stdout,'(3a,i8,a,i6.6)') 
     &                        'Read ', trim(fname), ' on ', nymd, ' at ', nhms
                    endif
                    meanres = .true.
                endif

!               If creating speed out of u/v components ...
!               -------------------------------------------
                if (raobspd) then
                   call generate_raob_speed(ods) ! from this point on, u holds speed 
                                                 ! and all v components have qcx=31
                endif

!               Set number of observations found in the file
!               --------------------------------------------
                nobs   =  ods%data%nobs

                if ( nobs .eq. 0 ) then
                     print *, 'No data for this synoptic time'
                     call ods_clean ( ods, ier )
                     if (trim(ODSmean_fntmpl) /= 'NONE' ) then
                         call ods_clean ( mods, ier )
                     endif
                     cycle    ! skip to next synoptic time
                end if

!               In case there is redundancy in files (ensemble case, when there
!               are multiple files valid at the same time), determine how many 
!               redundant files ...
!               ---------------------------------------------------------------
                if(nymd_ens==nymd .and. nhms_ens==nhms) then
                   nredundant = nredundant + 1
                else
                   nymd_ens=nymd; nhms_ens=nhms
                   nredundant = 1
                endif

!               Reduce input data as requested
!               ------------------------------
!               if ( DO_red ) then
!                    call Reducer ( nymd, nhms, nobs, ods, nobs_good, rcfile=RC_red )
!                    ods%data%nobs = nobs_good
!                    nobs          = nobs_good
!                    print *, 'Remaining obs after Reducer: ', nobs
!               end if

                allocate ( lat(nobs), lon(nobs), lev(nobs), qc (nobs), 
     &                     del(nobs,nr),   kx(nobs), kt (nobs), ks(nobs),
     &                     ob_flag (nobs), stat=ier )
                   if(ier/=0) call die (myname,
     &                                ' Obs Arrays Error, alloc(ods)')

                allocate ( is(nobs) )
                call IndexSet  ( nobs, is )
                call IndexSort ( nobs, is, ods%data%qcexcl(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%   lat(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%   lon(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%   lev(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%  time(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%    kt(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%    ks(1:nobs), descend=.false. )
                call IndexSort ( nobs, is, ods%data%    kx(1:nobs), descend=.false. )

                ods%data%kt    (1:nobs) = ods%data%kt    ( (/ (is(i), i=1,nobs) /) )
                ods%data%kx    (1:nobs) = ods%data%kx    ( (/ (is(i), i=1,nobs) /) )
                ods%data%ks    (1:nobs) = ods%data%ks    ( (/ (is(i), i=1,nobs) /) )
                ods%data%lon   (1:nobs) = ods%data%lon   ( (/ (is(i), i=1,nobs) /) )
                ods%data%lat   (1:nobs) = ods%data%lat   ( (/ (is(i), i=1,nobs) /) )
                ods%data%lev   (1:nobs) = ods%data%lev   ( (/ (is(i), i=1,nobs) /) )
                ods%data%time  (1:nobs) = ods%data%time  ( (/ (is(i), i=1,nobs) /) )
                ods%data%obs   (1:nobs) = ods%data%obs   ( (/ (is(i), i=1,nobs) /) )
                ods%data%OmF   (1:nobs) = ods%data%OmF   ( (/ (is(i), i=1,nobs) /) )
                ods%data%OmA   (1:nobs) = ods%data%OmA   ( (/ (is(i), i=1,nobs) /) )
                ods%data%xm    (1:nobs) = ods%data%xm    ( (/ (is(i), i=1,nobs) /) )
                ods%data%qcexcl(1:nobs) = ods%data%qcexcl( (/ (is(i), i=1,nobs) /) )
                ods%data%qchist(1:nobs) = ods%data%qchist( (/ (is(i), i=1,nobs) /) )
                ods%data%Xvec  (1:nobs) = ods%data%Xvec  ( (/ (is(i), i=1,nobs) /) )

                if ( meanres ) then ! the following assumes the members of the
                                    ! ensemble have been ordered as the mean
                                    ! (i.e., each member has gone thru odsmatch w/ mean)
                   mods%data%kt    (1:nobs) = mods%data%kt    ( (/ (is(i), i=1,nobs) /) )
                   mods%data%kx    (1:nobs) = mods%data%kx    ( (/ (is(i), i=1,nobs) /) )
                   mods%data%ks    (1:nobs) = mods%data%ks    ( (/ (is(i), i=1,nobs) /) )
                   mods%data%lon   (1:nobs) = mods%data%lon   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%lat   (1:nobs) = mods%data%lat   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%lev   (1:nobs) = mods%data%lev   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%time  (1:nobs) = mods%data%time  ( (/ (is(i), i=1,nobs) /) )
                   mods%data%obs   (1:nobs) = mods%data%obs   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%OmF   (1:nobs) = mods%data%OmF   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%OmA   (1:nobs) = mods%data%OmA   ( (/ (is(i), i=1,nobs) /) )
                   mods%data%xm    (1:nobs) = mods%data%xm    ( (/ (is(i), i=1,nobs) /) )
                   mods%data%qcexcl(1:nobs) = mods%data%qcexcl( (/ (is(i), i=1,nobs) /) )
                   mods%data%qchist(1:nobs) = mods%data%qchist( (/ (is(i), i=1,nobs) /) )
                   mods%data%Xvec  (1:nobs) = mods%data%Xvec  ( (/ (is(i), i=1,nobs) /) )
                endif

                lat     = ods%data%lat   (1:nobs)
                lon     = ods%data%lon   (1:nobs)
                lev     = ods%data%lev   (1:nobs)
                kt      = ods%data%kt    (1:nobs)
                kx      = ods%data%kx    (1:nobs)
                ks      = ods%data%ks    (1:nobs)
                qc      = real(ods%data%qcexcl(1:nobs))
                ob_flag ( : nobs ) = 0

                deallocate(is)

!               Force y2K compliance
!               --------------------
                if ( nymd .lt. 19000000 ) nymd = nymd + 19000000
                nsyn = ods%meta%nsyn

!               Collect requested residual (field) (O-F, O-A, obs, etc)
!               -------------------------------------------------------
                if ( iopt .eq. 0 .or. iopt .eq. 1 ) then
                     del(:,1) = ods%data%omf(1:nobs)
                     if(scale_res==1) then
                       del(:,1) = del(:,1)/sqrt(ods%data%xvec(1:nobs))
                     endif
                     if(scale_res==2) then
                       del(:,1) = del(:,1)/ods%data%obs(1:nobs)
                     endif
                     if(scale_res==3) then
                       del(:,1) = del(:,1)/(ods%data%obs(1:nobs)-ods%data%omf(1:nobs))
                     endif
                     if(xcov.or.pxcov)
     .               del(:,2) = del(:,1) ! ods%data%omf(1:nobs)
                else if ( iopt .eq. 2 .or. iopt .eq. 3 ) then
                     del(:,1) = ods%data%oma(1:nobs)
                     if(scale_res==1) then
                       del(:,1) = del(:,1)/sqrt(ods%data%xvec(1:nobs))
                     endif
                     if(scale_res==2) then
                       del(:,1) = del(:,1)/ods%data%obs(1:nobs)
                     endif
                     if(scale_res==3) then
                       del(:,1) = del(:,1)/(ods%data%obs(1:nobs)-ods%data%omf(1:nobs))
                     endif
                     if(xcov.or.pxcov)  
     .               del(:,2) = del(:,1) ! ods%data%oma(1:nobs)
                else if ( iopt .eq. 4 .or. iopt .eq. 5 ) then
                     del(:,1) = ods%data%obs(1:nobs)
                else if ( iopt .eq. 6 .or. iopt .eq. 7 ) then    ! prescribed sigo
                     del(:,1) = ods%data%xvec(1:nobs)
                else if ( iopt .eq. 8 .or. iopt .eq. 9 ) then    ! observation bias
                     del(:,1) = ods%data%xm(1:nobs)
                else if ( iopt .eq.10 .or. iopt .eq.11 ) then    ! X-cov <AmF,OmF> ~ HBH'
                     if (tch) then                        ! (x,y,z) -> (o,b,a)
                                                          ! though the signs are
                                                          ! irrelant this follows
                                                          ! the actual rotation
                                                          ! of indexes in
                                                          ! Todling et al. (2022)
                       del(:,1) =  ods%data%oma(1:nobs) - ods%data%omf(1:nobs) ! (y-z)
                       del(:,2) = -ods%data%omf(1:nobs)   ! (y-x)
                       del(:,3) = -ods%data%oma(1:nobs)   ! (z-x)
                      !del(:,1) = ods%data%omf(1:nobs)
                      !del(:,2) = ods%data%oma(1:nobs)
                      !del(:,3) = ods%data%omf(1:nobs) - ods%data%oma(1:nobs)
                      !islot = 2
                     else
                       del(:,1) = ods%data%omf(1:nobs) - ods%data%oma(1:nobs)
                       del(:,2) = ods%data%omf(1:nobs)
                     endif
                else if ( iopt .eq.12 .or. iopt .eq.13 ) then    ! X-cov <AmF,OmA> ~ HPaH'
                     if (tch) then                        ! (x,y,z) -> (o,b,a)
                       del(:,1) = -ods%data%oma(1:nobs)   ! (z-x)
                       del(:,2) =  ods%data%omf(1:nobs) - ods%data%oma(1:nobs) ! (z-y)
                       del(:,3) =  ods%data%omf(1:nobs)   ! (x-y)
                      !del(:,1) = ods%data%oma(1:nobs)
                      !del(:,2) = ods%data%omf(1:nobs)
                      !del(:,3) = ods%data%omf(1:nobs) - ods%data%oma(1:nobs)
                      !islot = 3
                     else
                       del(:,1) = ods%data%omf(1:nobs) - ods%data%oma(1:nobs)
                       del(:,2) = ods%data%oma(1:nobs)
                     endif
                else if ( iopt .eq.14 .or. iopt .eq.15 ) then    ! X-cov <OmA,OmF> ~ R
                     if(tch) then                         ! (x,y,z) -> (o,b,a)
                       if ( self ) then
                          del(:,1) = ods%data%obs(1:nobs)    ! (x)
                          del(:,2) = ods%data%obs(1:nobs) - ods%data%omf(1:nobs) ! (y)
                          del(:,3) = ods%data%obs(1:nobs) - ods%data%oma(1:nobs) ! (z)
                       else
                          del(:,1) = ods%data%omf(1:nobs)    ! (x-y)
                          del(:,2) = ods%data%oma(1:nobs)    ! (x-z)
                          del(:,3) = ods%data%oma(1:nobs) - ods%data%omf(1:nobs) ! (y-z)
                       endif
                     else
                       del(:,1) = ods%data%oma(1:nobs)
                       del(:,2) = ods%data%omf(1:nobs)
                     endif
                else if ( iopt .eq.16 .or. iopt .eq.17 ) then    ! Jo(OmF)/p
                     del(:,1) = ods%data%omf(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                     nstat = 1
                else if ( iopt .eq.18 .or. iopt .eq.19 ) then    ! Jo(OmA)/p
                     del(:,1) = ods%data%oma(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                     nstat = 1
                else if ( iopt.eq.20 .or. iopt.eq.21 ) then      ! original observation
                     del(:,1) = ods%data%obs(1:nobs) + ods%data%xm(1:nobs)
                else if ( iopt.eq.22 .or. iopt.eq.23 ) then      ! observation impacts
                     del(:,1) = ods%data%xvec(1:nobs)
                else if ( iopt.eq.24 .or. iopt.eq.25 ) then      ! background at obs location
                     if (meanres) then ! del = f_i - <f>
                        ! make sure the match is proper: don't vectorize this loop
                        do iii=1,nobs
                           if(ods%data%obs(iii)==mods%data%obs(iii)) then
                              del(iii,1) = mods%data%omf(iii) -  ods%data%omf(iii)
                           else
                              ods%data%qcexcl(iii) = X_PASSIVE
                           endif
                        enddo
                     else              ! del = f_i
                         del(:,1) = ods%data%obs(1:nobs) - ods%data%omf(1:nobs)
                     endif
                else if ( iopt.eq.26 .or. iopt.eq.27 ) then      ! background at obs location
                     if (meanres) then ! del(1) = f_i - <f>; del(2) = sqrt(R)
                        ! make sure the match is proper: don't vectorize this loop
                        do iii=1,nobs
                           if(ods%data%obs(iii)==mods%data%obs(iii)) then
                              del(iii,1) = mods%data%omf(iii) -  ods%data%omf(iii) ! 
                              del(iii,2) =  ods%data%xvec(iii)   ! sigO slot
                           else
                              ods%data%qcexcl(iii) = X_PASSIVE
                           endif
                        enddo
                     endif
                else if ( iopt .eq. 28 .or. iopt .eq. 29 ) then ! h(xf)
                     del(:,1) = ods%data%obs(1:nobs) - ods%data%omf(1:nobs) + ods%data%xm(1:nobs)
                else if ( iopt .eq. 30 .or. iopt .eq .31 ) then ! analysis at obs location
                     if (meanres) then ! del = a_i - <a>
                        ! make sure the match is proper: don't vectorize this loop
                        do iii=1,nobs
                           if(ods%data%obs(iii)==mods%data%obs(iii)) then
                              del(iii,1) = mods%data%oma(iii) -  ods%data%oma(iii)
                           else
                              ods%data%qcexcl(iii) = X_PASSIVE
                           endif
                        enddo
                     else              ! del = a_i
                         del(:,1) = ods%data%obs(1:nobs) - ods%data%oma(1:nobs)
                     endif
                else if ( iopt .eq. 32 .or. iopt .eq .33 ) then ! analysis at obs location
                     del(:,1) = ods%data%omf(1:nobs)-ods%data%oma(1:nobs)
                     del(:,2) = ods%data%oma(1:nobs)
                     where (ods%data%xvec(1:nobs) > 0.0)
                        del(:,1) = del(:,1)/ods%data%xvec(1:nobs)
                        del(:,2) = del(:,2)/ods%data%xvec(1:nobs)
                     elsewhere
                        del(:,1) = amiss
                        del(:,2) = amiss
                     endwhere
                else
                     call die (myname, 'iopt not valid')
                endif

!               Now that all atrributes have been read, eliminate bad obs
!                and those that have passed QC flag, if user desires. 
!                Notice that passive data types are kept if ioptn>0.
!               --------------------------------------------------------
                if ( mod(iopt,2) .ne. 0 ) then
                   j = 0
                   do i = 1, nobs
                      iampassive   = ods%data%qcexcl(i) .eq. X_PASSIVE 
     &                          .or. ods%data%qcexcl(i) .eq. X_PRE_BAD
                      iwantpassive = ioptn .gt. 0 .and. iampassive
                      if (t2<t1) then
                          timeselect = .TRUE.
                      elseif (t1<=ods%data%time(i) .and. ods%data%time(i)<=t2) then
                          timeselect = .TRUE.
                      else
                          timeselect = .FALSE.
                      end if
                      if (lona<-180.0.and.lonb>180.0) then
                          lonselect = .TRUE.
                      elseif (lona<=ods%data%lon(i) .and. ods%data%lon(i)<=lonb) then
                          lonselect = .TRUE.
                      else
                          lonselect = .FALSE.
                      end if
                      if ( timeselect .and. lonselect .and.
     &                    (ods%data%qcexcl(i) .eq. 0 .or.
     &                     iwantpassive) ) then
                         ob_flag ( i ) = 0
                      else 
                         ob_flag ( i ) = 1
                      end if
                   end do
                   if(iopt==17 .or. iopt==19) del(1:nobs,1) = del(1:nobs,1)/del(1:nobs,2)  ! omf/sigo for Joa/Job calculation
                end if

             end if

!            Apply mask
!            ----------
             call gritas_maskout (lat,lon,ob_flag,lwimask)

!            Enforce Y2K compliance, etc.
!            ----------------------------
             if ( nymd .lt. 19000000 ) nymd = nymd + 19000000
             nymd_max = max ( nymd_max, nymd )
             nymd_min = min ( nymd_min, nymd )

!            Accumulate grids
!            ----------------
             if (xcov.or.pxcov) then
                if (xcov) then
                  call GrITAS_XCov_Accum ( lat, lon, lev, kx, kt, del(1:nobs,1:nr), nobs, nr,
     &                                     kxkt_table,
     &                                     ob_flag,
     &                                     l3d, bias_lv, xcov_lv, nobs_lv, nprof, tch )
                else
                  call GrITAS_XCov_Prs_Accum ( lev, kx, kt, ks, del(1:nobs,1:nr), nobs, nr,
     &                                         kxkt_table, p1, p2,
     &                                         ob_flag,
     &                                         l3d, bias_lv, xcov_lv, nobs_lv, nprof, tch )
                endif
                nprof_total = nprof_total + nprof
                print *, ' nprof_total = ', nprof_total
             else
               if (lminmax) then
                  call GrITAS_MinMax( lat, lon, lev, kx, kt, del(1:nobs,1:nr), nobs, nr,
     &                                kxkt_table, p1, p2,
     &                                ob_flag,
     &                                l2d, bias_2d, stdv_2d, nobs_2d,
     &                                l3d, bias_3d, stdv_3d, nobs_3d )
               else
                  call GrITAS_Accum ( lat, lon, lev, kx, kt, del(1:nobs,1:nr), nobs, nr,
     &                                kxkt_table, p1, p2,
     &                                ob_flag,
     &                                l2d, bias_2d, rtms_2d, xcov_2d(:,:,:,1),   nobs_2d,
     &                                l3d, bias_3d, rtms_3d, xcov_3d(:,:,:,:,1), nobs_3d )

                  if (scorecard>=0) then
                      call GrITAS_Scores ( expid, nymd, nhms, scorecard,
     &                                     lat, lon, lev, kx, kt, del, nobs, nr,
     &                                     kt_list(1:listsz), var_list(1:listsz), kxkt_table, p1, p2,
     &                                     ods%data%qcexcl, l2d, l3d )
                  endif
               endif
             endif
             if ( ftype        .ne. 'del' .and.
     &            lb_tallyfile .ne.  ' '  .and.
     &            nobs         .gt.  0 ) then
                call GrITAS_LBTally ( del_ifname(idel), nymd, nhms,
     &                                ob_flag, ods,  rc )
             end if

             deallocate ( lat, lon, lev, qc, del, kx, kt, ks, ob_flag )

!            Clean up mean-carrying ODS structure
!            ------------------------------------
             if (meanres) then
                call ods_clean ( mods, ier )
             endif
 10       continue
 11       continue

      end do

!     Three-corner hat (but not x-cov)
!     --------------------------------
      if ((.not.xcov).and.(.not.pxcov)) then
         if (tch) then
            call GrITAS_3CH_Norm ( nr, nrmlz,
     &                             l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                             l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d)
         endif
      endif

!     Normalize and post process rms and bias into covariance
!     -------------------------------------------------------
      if ( (xcov .or. pxcov) .and.  nprof_total>0 ) then

!           Normalize and remove bias term
!           ------------------------------
            call GrITAS_XCov_Norm ( nr, .true., tch, self, l3d, bias_lv, xcov_lv, nobs_lv )

!           Symmetrize since cross-cov not really cov
!           -----------------------------------------
            nrout = 1
            if (tch) nrout=nr
            if (symmetrize) then
               print *, 'Xcov is being symmetrized'
               do kk=1,nrout
                  do ll=1,l3d
                     xcov_lv(:,:,ll,kk) = 0.5 *( xcov_lv(:,:,ll,kk) + transpose(xcov_lv(:,:,ll,kk)) )
                  enddo
               enddo
            endif

!           Write out covariance matrix
!           ---------------------------
            do kk=1,nrout
               write(covname,'(2a,i1,a)') trim(grid_ofname), 'cov', kk, '.bin'
               open(32,file=trim(covname),form='unformatted',convert='little_endian')
               write(32) real(xcov_lv(1:,1:,1,kk),4)
               close(32)
            enddo

!           Write obs count 
!           ---------------
            write(covname,'(2a)') trim(grid_ofname), 'nobs.bin'
            open(32,file=trim(covname),form='unformatted',convert='little_endian')
            write(32) real(nobs_lv(:,1),4)
            close(32)

!           For GSI: always symmetrize
!           --------------------------
            if (.not.symmetrize) then
               do kk=1,nrout
                  do ll=1,l3d
                     xcov_lv(:,:,ll,kk) = 0.5 *( xcov_lv(:,:,ll,kk) + transpose(xcov_lv(:,:,ll,kk)) )
                  enddo
               enddo
            endif

!           GSI never care about A or B 
!           (The following clearly assumes nr=1 carries R)
!           ----------------------------------------------
            covname = trim(grid_ofname) // 'cov4gsi.bin' 
            open(32,file=trim(covname),form='unformatted',convert='little_endian')
            if (nchan>0) then
               write(32) nchan, 8
               write(32) achan(:,1)
            else
               write(32) nlevs, 8
               write(32) plevs
            endif
            write(32) real(xcov_lv(1:,1:,1,1),8)
            close(32)

      else

            if (.not. lminmax .and. .not.tch) then
               call GrITAS_Norm ( nr, nrmlz,
     &                            l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                            l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d)
            endif

      endif

      if (.not. xcov .and. .not.pxcov) then ! <xcov>

!         If desired, save RTMS of omf or oma instead of stdv
!         ---------------------------------------------------
          if ( rtms .and. (
     &         iopt .eq. 0  .or. iopt .eq.  1 .or.
     &         iopt .eq. 2  .or. iopt .eq.  3 ) ) then
               stdv_2d = rtms_2d
               stdv_3d = rtms_3d
               print *
               print *, 'Replace Std-Dev with root time mean square (RTMS)'
               print *
          endif

!         STDV not always makes sense ...
!         -------------------------------
          if ( iopt==6 .or.iopt==7 .or. iopt==8 .or.iopt==9 ) then
               stdv_2d = rtms_2d
               stdv_3d = rtms_3d
               print *
               print *, 'Warning: Std-Dev slot holds root time mean square (RTMS) instead'
               print *
          endif

!         In case of x-covariances, write x-rms in the bias slot, and x-variances in stdv slot
!         ------------------------------------------------------------------------------------
          if ( iopt==10 .or.iopt==11 .or. iopt==12 .or.iopt==13 .or. iopt==14 .or.iopt==15 ) then
               if (tch) then
                  stdv_2d = xcov_2d 
                  stdv_3d = xcov_3d
               else
                  do ll=1,l2d
                  do jj=1,jm
                  do ii=1,im
                     bias_2d(ii,jj,ll,1) = sqrt(max(0.,xcov_2d(ii,jj,ll,1)))
                     enddo
                  enddo
                  enddo
                  stdv_2d(:,:,:,  1) = xcov_2d(:,:,:,  2)
                  do ll=1,l3d
                  do kk=1,km
                  do jj=1,jm
                  do ii=1,im
                     bias_3d(ii,jj,kk,ll,1) = sqrt(max(0.,xcov_3d(ii,jj,kk,ll,1)))
                  enddo
                  enddo
                  enddo
                  enddo
                  stdv_3d(:,:,:,:,1) = xcov_3d(:,:,:,:,2)
               endif
          endif

!         In case of Jo(OmA) and Jo(OmF) store Jo in bias array
!         -----------------------------------------------------
          if ( iopt==16 .or.iopt==17 .or. iopt==18 .or.iopt==19 ) then
               bias_2d(:,:,:,  1) = rtms_2d(:,:,:,  1)
               bias_3d(:,:,:,:,1) = rtms_3d(:,:,:,:,1)
          endif
       
          if ( nredundant>1 ) then

!            In this case:
!            the first  entry (nr=1) holds the forecast spread
!            the second entry (nr=2) holds the sigo
!            ---------------------------------------------------------
             if ( iopt==26 .or. iopt==27 ) then
                v1=0.
                v2=0.
                do ll=1,l2d
                do jj=1,jm
                do ii=1,im
                   if (abs(stdv_2d(ii,jj,ll,1)-AMISS)>0.01 .and. abs(bias_2d(ii,jj,ll,2)-AMISS)>0.01)  then 
                      v1 = stdv_2d(ii,jj,ll,1)*stdv_2d(ii,jj,ll,1)
                      v2 = bias_2d(ii,jj,ll,2)*bias_2d(ii,jj,ll,2)  ! stdv slot of sigo carries nill since all members use same sigo
                      stdv_2d(ii,jj,ll,1) = sqrt(v1+v2)
                   else
                      stdv_2d(ii,jj,ll,1) = AMISS
                   endif
                enddo
                enddo
                enddo 
                v1=0.
                v2=0.
                do ll=1,l3d
                do kk=1,km
                do jj=1,jm
                   do ii=1,im
                      if (abs(stdv_3d(ii,jj,kk,ll,1)-AMISS)>0.01 .and. abs(bias_3d(ii,jj,kk,ll,2)-AMISS)>0.01)  then 
                      v1 = stdv_3d(ii,jj,kk,ll,1)*stdv_3d(ii,jj,kk,ll,1)
                      v2 = bias_3d(ii,jj,kk,ll,2)*bias_3d(ii,jj,kk,ll,2) ! stdv slot of sigo carries (nearly) nill
                                                                         ! since all members use (nearly) same sigo
                      stdv_3d(ii,jj,kk,ll,1) = sqrt(v1+v2)
                   else
                      stdv_3d(ii,jj,kk,ll,1) = AMISS   
                   endif
                enddo           
                enddo           
                enddo           
                enddo
             endif    
             print *
             print *, ' [] Number of redundant samples (ensemble): ', nredundant

          endif
       
          if (lminmax) then
             where(abs(stdv_2d)==AMISS) stdv_2d=AMISS
             where(abs(stdv_3d)==AMISS) stdv_3d=AMISS
          endif

          nobs_2d = nobs_2d / nredundant
          nobs_3d = nobs_3d / nredundant

!         Calculate Chi-square for specified confidence interval
!         ------------------------------------------------------      
          if (confidence>0.0) then
             call GrITAS_CDF ( 'chisq-dist', 1.0-0.5*confidence, nobs_2d, chsq_2d(:,:,:,1), nobs_3d, chsq_3d(:,:,:,:,1) )
             call GrITAS_CDF ( 'chisq-dist',     0.5*confidence, nobs_2d, chsq_2d(:,:,:,2), nobs_3d, chsq_3d(:,:,:,:,2) )
             call GrITAS_CDF (     't-dist',     0.5*confidence, nobs_2d, chsq_2d(:,:,:,3), nobs_3d, chsq_3d(:,:,:,:,3) )
          endif

!         Write del grids to file
!         -----------------------
          if ( trim(grid_ofname) /= "NONE" ) then
             if ( nymd_out < 0 ) then
                nymd = ( nymd_min + nymd_max ) / 2
             else
                nymd = nymd_out ! user specified
             endif
             nhms = 120000
             timinc = 240000 / nsyn ! frequency of ouput
             if ( trim(otype) .eq. 'hdf'  .or. trim(otype) .eq. 'nc4') then
                if ( split_tgroups ) nstat = -1 * nstat
                call GFIO_Output   ( grid_ofname,otype,
     &                               var_list, var_descr, listsz, timinc,
     &                               l2d, list_2d, bias_2d(:,:,:,1)  , stdv_2d(:,:,:,islot)  , nobs_2d,
     &                               l3d, list_3d, bias_3d(:,:,:,:,1), stdv_3d(:,:,:,:,islot), nobs_3d,
     &                               nstat, nymd, nhms, fid, chsq_2d=chsq_2d, chsq_3d=chsq_3d ) 

             else

                grid_ofname = trim(grid_ofname)//'.dat' 
                call GrADS_Output (  trim(grid_ofname),
     &                               var_list, listsz,
     &                               l2d, list_2d, bias_2d(1,1,1,1)  , stdv_2d(1,1,1,islot)  , nobs_2d,
     &                               l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,islot), nobs_3d, nstat)

             end if
          end if

      endif  ! <xcov>

!     Summary
!     -------
      if (nprof_total>0) then
         print *, 'GrITAS, total overall profiles processed = ', nprof_total
      endif

!     Clean up
!     --------
      if (xcov.or.pxcov) then
          deallocate ( bias_lv, xcov_lv, nobs_lv, stat=ier )
            if(ier/=0) call die (myname,'Error in Alloc(lv)')
      else
          deallocate ( bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d, chsq_3d, stat=ier )
            if(ier/=0) call die (myname,'Error in Alloc(3d)')

          deallocate ( bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d, chsq_2d, stat=ier )
            if(ier/=0) call die (myname,'Error in Dealloc(2d)')
      endif

      call GrITAS_Clean ()

!     All done.
!     --------
      if ( trim(otype) .eq. 'hdf' .and. trim(otype) .eq. 'nc4' ) then
          call GFIO_Close ( fid(1), rc )   ! close GFIO file
          deallocate(fid) 
      end if

      stop
      end

!EOC

!...........................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_LBTally --- Write summary of obs not included in statistics
! 
! !INTERFACE:
!

      subroutine GrITAS_LBTally ( ods_file, date, time,
     &                            ob_flag,  ods,  rc )
!
! !USES:
      use m_ods,   only : ods_vect, 
     &                    ODS_Merge, ODS_Select, ODS_Tally, ODS_Clean
      use m_stdio, only : stdout
!
! !INPUT PARAMETERS:
      character (len=*), intent (in)  ::
     &   ods_file
      integer,           intent (in)  ::
     &   date,
     &   time,
     &   ob_flag (*)
      type (ods_vect),   intent (in)  ::
     &   ods
!
! !OUTPUT PARAMETERS:
      integer,           intent (out) ::
     &   rc

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!     10Jul09  Redder   Initial code.
!EOP
!-------------------------------------------------------------------------

      type ( ods_vect ) :: ods_out
      integer :: iob, nobs, nsel, lu
      integer, dimension (:), pointer :: qchist

      lu = stdout
      nobs   =  ods     % data % nobs

!     Transfer obs to local space
!     ---------------------------
      call ODS_Merge ( (/ ods /), 1, ods_out, rc )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Merge.' )
         return
      end if

!     Put flag values to local ods structure
!     --------------------------------------
      nobs   =  ods     % data % nobs
      qchist => ods_out % data % qchist ( : nobs )
      ods_out % data % qchist ( : nobs ) = ob_flag ( : nobs )

!     Select only the obs not included in the statistical computations
!     ----------------------------------------------------------------
      call ODS_Select ( ods_out, nobs, nsel, rc,
     &                  qchist = 1 )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Select.' )
         return
      end if

!     Write header and ...
!     --------------------
      write (lu, '(4(a/)a,i8/a,i6)',iostat = rc)
     &    'ods_tally:===========================================',
     &    'ods_tally:',
     &    'ods_tally: Talley of obs not used in the computation '
     &                               //  ' of statistics',
     &    'ods_tally:   ods file = ' // trim ( ods_file ),
     &    'ods_tally:   date     = ', date,
     &    'ods_tally:   time     = ', time
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from write statement.' )
         return
      end if

!     ... summary of selected obs
!     ---------------------------
      call ODS_Tally ( lu, ods_out, nsel, rc )
      if ( rc .ne. 0 ) then
         call gritas_perror ( 'GrITAS_LBTally: error status '
     .                     // 'returned from ODS_Tally.' )
         return
      end if

!     Clean up
!     --------
      call ODS_Clean ( ods_out, rc )
 
      return
      end 

!...........................................................................


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrITAS_Init --- Command line and resource user interface
! 
! !INTERFACE:
!
      subroutine GrITAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                         ftype, dtype, otype, iopt, ncf, nr,
     &                         split_tgroups, lb_tallyfile,
     &                         nymd_out, DO_red, RC_red,
     &                         ndel, del_ifname, grid_ofname, odsmean_fntmpl,
     &                         kt_list, kx_list, var_list, var_descr,
     &                         listsz, nrmlz, xcov, pxcov, self, rtms, t1, t2, lwimask,
     &                         confidence, scale_res,
     &                         lona, lonb, lminmax, expid, scorecard, tch,
     &                         symmetrize, raobspd )


! !USES:
 
      Use m_gritas_grids, only : grids_init
      Use m_gritas_grids, only : im, jm, km
      Use m_gritas_grids, only : nlevs
      Use m_gritas_grids, only : plevs
      Use m_gritas_grids, only : hlevs
      Use m_gritas_grids, only : nchan
      Use m_gritas_grids, only : achan
      Use m_gritas_grids, only : GrITAS_hydroheights
      Use m_gritas_grids, only : ichan_version
      Use m_gritas_masks, only : gritas_ini_masks
      Use m_inpak90
      Use m_stdio, only : stderr
      Use m_die

!
! !INPUT PARAMETERS:
!
      implicit        NONE
      integer         NDEL_MAX             ! max number of del files    
      integer         LMAX                 ! max list size
      integer         KXMAX                ! max number of kx's

!
! !OUTPUT PARAMETERS:
!

      character(len=*) ftype                 ! File type: 'ods' or 'del'
      character(len=*) dtype                 ! Data type: 'omf', 'oma' or 'obs'
      character(len=*) otype                 ! Output type: 'hdf' or 'eee'
      integer          nr                    ! number of distinct residual series
      integer          iopt                  ! =1 for O-F
                                             ! =3 for OMA 
                                             ! =5 for obs
      integer          nymd_out              ! date for output file
      logical          ncf                   ! Non-compliant format (to allow
                                             !   handling of GSI-output files)
      logical         split_tgroups          ! = true to split the tgroups,
                                             ! bias, stdv and nobs into 
                                             ! separate files.
      character(len=*) lb_tallyfile          ! name of file containing the
                                             !   summary of ODS obs not
                                             !   included in the summary.
                                             !   If blank name then no
                                             !   summary is printed.  Other-
                                             !   wise the filenam is set to
                                             !   stdout for standard output
                                             !   for this implementation.
      integer         ndel                   ! actual number of del files
      character*255   del_ifname(NDEL_MAX)   ! del file names
      character*255   grid_ofname            ! output grid file name

      integer         listsz                 ! actual list size (<LMAX)
      integer         kt_list(LMAX)          ! data type for the grid
      integer         kx_list(KXMAX,LMAX)    ! data sources for the grid
      character*11    var_list(LMAX)         ! variable name for the grid
      character*256   var_descr(LMAX)        ! ... and its description
      character(len=*) lwimask               ! land/ocean mask
      logical          DO_red                ! Defines whether to apply reducer or not
      logical          xcov                  ! knob controlling estimation of R and B fill covs (chn)
      logical          pxcov                 ! knob controlling estimation of R and B fill covs (lev)
      logical          self                  ! get R, B, A full covariance within xcov opt
      logical          rtms                  ! controls write out of RTMS instead of STDV
      character(len=*) RC_red                ! Recucer resource filename
      character(len=*) ODSmean_fntmpl        ! Template of filename with mean ODS files (for EnADAS only)
      real             confidence            ! Confidence interval (used in Chi-sqr check)
      integer          scale_res             ! Scale residual by sqrt(R), or o, or f
      logical          lminmax               ! Get min/max values within box

      logical          tch                   ! Enable three-corner hat calculation
      logical          symmetrize            ! Symmetrize x-cov estimates of cov
      logical          raobspd               ! Calculate RAOB speed

      integer,intent(inout) :: t1, t2        ! time interval of obs to be considered in stats
      real   ,intent(inout) :: lona,lonb     ! longitude interval of obs to be considered in stats
      character(len=*),intent(inout) :: expid! experiment identifier
      integer,intent(inout) :: scorecard     ! generate info for scorecard

! !INPUT/OUTPUT PARAMETERS:

       logical         nrmlz                 ! Controls normalization of accumalated

! !DESCRIPTION: This routine parses the command line for the name of
!               the input/output file names and read the resource file
!  file 'gritas.rc' for defining the combination of data type/data
!  sources for each of the output grids.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  ???????  G. P. Lou  Modifications.  
!  15Dec98  da Silva   Added ods/gfio, o-a, obs support.
!  25mar99  da Silva   Added -rc option
!  16nov99  da Silva   Added -nopassive
!  09Jun04  Todling    Added -res; updated to use m_gritas_grids
!  14Jun04  Todling    Changed reading of RC file kt/kx table to allow ":"
!  30Jul04  Todling    Added reducer option; bug fix: added i90_release
!  28Oct06  Todling    Added -sigo, and -obias options
!  07Nov06  Todling    Added various opts: joa/job/esigo/hbh/hah
!  23Jan08  Todling    Added -obsimp
!  12Sep08  Todling    Added -nonorm
!  28May09  Redder     Added capability to define kt as a surface variable
!                      by added the prefix, sfckt_mkr, to its name in the
!                      input rc file
!  08Jun09  Redder     Added the command line option, -ospl, and
!                      the subroutine output parameter, split_tgroups
!  10Jul09  Redder     Added the command line option, -lb and
!                      the subroutine output parameter, lb_tallyfile
!  20Jul09  Redder     Added the subroutine output parameter, var_descr.
!  21Jul09  Redder     Fixed bug in generating the list, var_descr
!  23Apr10  Redder     Added feature count the levels listed in the rc
!                      file before initializing the grid and to list the
!                      range of levels using the character ":"
!  10May13  Todling    Add -fobloc; bug fix: -o was in conflict w/ other opts
!                      starting with "o"; add rc file opt for ens-spread
!  06Feb14 El Akkraoui Fixed bug in vertical loop, km instead of jm. 
!  10Feb14 El Akkraoui Added -promf to get prescribed ensemble innovations using 
!                      ensemble spread and sigo (sqrt(S+R)). 
!                      Fixed bug by Modifying dimensions of allocated variables 
!                      to use l2d/l3d instead of l2d_max/l3d_max. 
!  11Apr14  Todling    Mild changes to estimate inter-channel correlations using
!                      "Desrosizers" method.
!  01Aug14  Todling    Allow for use of real channel numbers (as well as indexes). 
!  09Nov14  Todling    Add chi-square diagnostic for standard deviations
!  16Sep15  Todling    Scale residual by sqrt(R)
!  08Jul20  Todling    Add min/max option
!  06Jan22  Todling    Option to produce three-corner hat estimates
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=*), parameter :: myname = 'GrITAS_Init'

!     Prefix for variable name denoting a surface variable
!     ----------------------------------------------------
      character (len=*), parameter :: sfckt_mkr = '_'

!     Temp storage for vertical levels
!     --------------------------------
      real, allocatable, save :: plevs2(:) 

      integer        iarg, argc, iargc
      character(len=255) argv, res, SS, rc_ifname, mask_fname
      character(len=2000) token
      integer        i, j, nkx, lt, ii, jj, ilist, isfckt_mkr, ic
      integer        iret, ios, iscan, idummy, ichan_indx, ichan_numb
      integer        kx1, kx2, kxnext, ilev_add, nlevs_add, nchan_add
      real           p, p1, p2, incr
      real           swap
      logical passive, long_name_found, long_name



      passive = .true.  ! include passive data in calc of stats
      ncf   = .false.
      ftype = 'ods'
      dtype = 'omf'
      otype = 'hdf'
      otype = 'nc4'
      iopt  = 1
      nr    = 1
      im=72; jm=46; km=26  ! default resolution
      DO_red = .false.
      RC_red = 'reducer.rc'
      ODSmean_fntmpl = 'NONE'
      nymd_out = -1
      xcov=.false.
      pxcov=.false.
      self=.false.
      rtms=.false.
      nchan=0
      confidence=-999.
      scale_res=-1
      lminmax=.false.
      tch = .false.
      symmetrize = .true.
      raobspd = .false.

!     Parse command line
!     ------------------
      rc_ifname = 'gritas.rc'
      grid_ofname = 'gritas'      
      split_tgroups = .false.
      mask_fname = 'NONE'
      lb_tallyfile  = ' '
      argc =  iargc()
!      print *, "begin to read command-line input"
      if ( argc .lt. 1 ) call gritas_usage()
      ndel = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
!      print *, "iarg = ", iarg
!      print *, "argc = ", argc
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if      ( trim(argv) .eq. '-del' )  then
            ftype = 'del'
         else if ( trim(argv) .eq. '-nopassive' )  then
            passive = .false.
         else if ( trim(argv) .eq. '-omf' )  then
            dtype = 'omf'
            iopt  = 1
         else if ( trim(argv) .eq. '-oma' )  then
            dtype = 'oma'
            iopt  = 3
         else if ( trim(argv) .eq. '-obs' )  then
            dtype = 'obs'
            iopt  = 5
         else if ( trim(argv) .eq. '-sigo' )  then
            dtype = 'sigo'
            iopt  = 7
         else if ( trim(argv) .eq. '-obias' )  then
            dtype = 'obias'
            iopt  = 9
         else if ( trim(argv) .eq. '-hbh' )  then
            dtype = 'hbh'
            iopt  = 11
            nr    = 2
         else if ( trim(argv) .eq. '-hah' )  then
            dtype = 'hah'
            iopt  = 13
            nr    = 2
         else if ( trim(argv) .eq. '-esigo' )  then
            dtype = 'esigo'
            iopt  = 15
            nr    = 2
         else if ( trim(argv) .eq. '-job' )  then
            dtype = 'job'
            iopt  = 17
            nr    = 2
         else if ( trim(argv) .eq. '-joa' )  then
            dtype = 'joa'
            iopt  = 19
            nr    = 2
         else if ( trim(argv) .eq. '-oobs' )  then
            dtype = 'oobs'
            iopt  = 21
         else if ( trim(argv) .eq. '-obsimp' )  then
            dtype = 'obsimp'
            iopt  = 23
         else if (trim(argv) .eq. '-fobloc' )  then
            dtype = 'fobloc'
            iopt  = 25
         else if (trim(argv) .eq. '-aobloc' )  then
            dtype = 'aobloc'
            iopt  = 31
         else if (trim(argv) .eq. '-dfs' )  then
            dtype = 'dfs'
            iopt  = 33
            nr    = 2
         else if ( trim(argv) .eq. '-raobspd' )  then
            raobspd = .true.
         else if ( trim(argv) .eq. '-xcov' )  then
            xcov = .true.
         else if ( trim(argv) .eq. '-pxcov' )  then
            pxcov = .true.
         else if ( trim(argv) .eq. '-rtms' )  then
            rtms = .true.
         else if ( trim(argv) .eq. '-resbysigo' )  then
            scale_res = 1
         else if ( trim(argv) .eq. '-resbyo' )  then
            scale_res = 2
         else if ( trim(argv) .eq. '-resbyf' )  then
            scale_res = 3
         else if ( trim(argv) .eq. '-nosym' )  then
            symmetrize = .false.
         else if (trim(argv) .eq. '-promf' )  then
            dtype = 'promf'
            iopt  = 27
            nr    = 2
         else if ( trim(argv) .eq. '-hxb' )  then
            dtype = 'hxb'
            iopt  = 29
         else if ( trim(argv) .eq. '-3ch' )  then
            tch = .true.
         else if ( trim(argv) .eq. '-self' )  then
            self = .true.
         else if (trim(argv) .eq. '-scorecard' )  then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) expid
            if ( ios /= 0 ) then
                 print *, 'cannot parse expid ...'
                 call exit(1)
            end if
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) scorecard
            if ( ios /= 0 ) then
                 print *, 'cannot parse scorecard ...'
                 call exit(1)
            end if
            dtype = 'score'
         else if ( trim(argv) .eq. '-reduce' )  then
            DO_red = .true.
         else if ( trim(argv) .eq. '-minmax' )  then
            lminmax = .true.
         else if ( trim(argv) .eq. '-nonorm' )  then
            nrmlz = .false.
         else if ( trim(argv) .eq. '-ncf' )  then
            ncf = .true.
         else if ( trim(argv) .eq. '-hdf' )  then
            otype = 'hdf'
         else if ( trim(argv) .eq. '-ieee' )  then
            otype = 'eee'
         else if ( trim(argv) .eq. '-o' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, grid_ofname )
         else if ( trim(argv) .eq. '-ospl' ) then
            split_tgroups = .true.
         else if ( trim(argv) .eq. '-lb' ) then
            lb_tallyfile = 'stdout'
         else if ( trim(argv) .eq. '-rc' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, rc_ifname )
         else if ( trim(argv) .eq. '-mask' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, mask_fname )
            iarg = iarg + 1
            call GetArg ( iArg, lwimask )
         else if ( trim(argv) .eq.'-rcred' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, RC_red )
         else if ( trim(argv) .eq. '-conf' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) confidence
            if ( ios /= 0 ) then
                 print *, 'cannot parse confidence interval ...'
                 call exit(1)
            end if
         else if ( trim(argv) .eq. '-nlevs' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) km
            if ( ios /= 0 ) then
                 print *, 'cannot parse nlevs ...'
                 call exit(1)
            end if
         else if ( trim(argv) .eq. '-date' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) nymd_out
            if ( ios /= 0 ) then
                 print *, 'cannot parse date ...'
                 call exit(1)
            end if
         else if ( trim(argv) .eq. '-res' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()
            iarg = iarg + 1
            call GetArg ( iArg, res )
            if ( res(1:1) == "a" ) then
                 im=72; jm=46
             else if ( res(1:1) == "b" ) then
                 im=144; jm=91
             else if ( res(1:1) == "c" ) then
                 im=288; jm=181
             else if ( res(1:1) == "d" ) then
                 im=576; jm=361
             else if ( res(1:1) == "D" ) then
                 im=540; jm=361
             else if ( res(1:1) == "e" ) then
                 im=1080; jm=721
             else
                 print *, 'unknown resolution ... taking default'
             endif
         elseif (trim(argv) .eq. '-lon') then
            if ( iarg+1 .gt. argc ) call gritas_usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is lona or lona:lonb
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore lona
                read(SS,*) lona
                lonb = lona
            else                   ! colon, therefore lona:lonb
                read(SS(1:ic-1) ,*) lona
                read(SS(ic+1:lt),*) lonb
            end if
!           if (lona<lonb) then    ! let's be nice to the user..
!               swap = lonb
!               lonb = lona
!               lona = swap
!           end if
         elseif (trim(argv) .eq. '-time') then
            if ( iarg+1 .gt. argc ) call gritas_usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) t1
                t2 = t1
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) t1
                read(SS(ic+1:lt),*) t2
            end if
            if (t2<t1) then    ! let's be nice to the user..
                swap = t2
                t2   = t1
                t1   = swap
            end if
         else
           ndel = ndel + 1
           del_ifname(ndel) = argv
         end if
      end do
 111  continue

      if ( ndel .lt. 1 ) call gritas_usage()               

      if ( .not. passive ) iopt = - abs(iopt)   ! exclude passive data

!     Consistency check: del files have only O-F
!     ------------------------------------------
      if ( ftype .eq. 'del' .and. iopt .ne. 1 ) then
         print *
         print *, 'Note: DEL files only have O-F, ignoring -oma/-obs'
         dtype = 'omf'
         iopt = 1
      end if

!     If xcov chosen, force nr=2 (this way can do omf^2, oma^2, etc)
!     --------------------------
      if(xcov.or.pxcov)  nr=2
      if(tch)  then
         nr=3
         print *
         print *, ' Using 3CH method'
         print *, ' ================'
      endif
      if ( self ) then
        print *
        if (tch) then
           print *, ' will return (<v-<v>)*(v-<v>)'
           print *, ' ============================'
         else
           print *, ' self opt only work together with xcov'
           print *, ' ====================================='
           call die(myname)
        endif
      endif

!     Echo the parameters
!     -------------------
      print *
      print *, '   Residual type: ', trim(dtype)
      print *, 'Input  File type: ', trim(ftype)
      print *, 'Output File type: ', trim(otype)
      print *, '       Del files: ', ndel
      do i = 1, ndel
         print *, '     o ', trim(del_ifname(i))
      end do
      if(trim(lwimask)/='NONE') then
      print * 
      print *, 'Mask file: ',trim(mask_fname)
      print *, 'Mask type: ',trim(lwimask)
      endif
      if ( scale_res==1 ) then
      print *
      print *, 'Residuals scaled by sqrt(R)'
      print *
      endif
      if ( scale_res==2 ) then
      print *
      print *, 'Residuals scaled by observations(yo)'
      print *
      endif
      if ( scale_res==3 ) then
      print *
      print *, 'Residuals scaled by background(xb)'
      print *
      endif
      if ( lminmax ) then
      print *
      print *, 'Will generated Min/Max values within box'
      print *
      endif
      
      print *
      print *, 'Grid output file: ', trim(grid_ofname)
      print *

      if ( split_tgroups ) then
         print *
         print *, ' Bias, stdv and nobs will be placed in the files '
         print *, ' with the basenames, '
         print *, '   ' // trim (grid_ofname) // '.bias,'
         print *, '   ' // trim (grid_ofname) // '.stdv and '
         print *, '   ' // trim (grid_ofname) // '.nobs'
         print *
      end if

!     Read resource file
!     ------------------
      call i90_LoadF ( trim(rc_ifname), iret )
      if ( iret .ne. 0 ) then
         print *, 'GrITAS: cannot open resource file ',
     &            trim(rc_ifname)
         call exit(1)
      end if

!     vertical levels
!     ---------------
      call i90_label ( 'GrITAS*ODSmean_template:', iret )
      if ( iret .eq. 0 ) then
         call I90_Gtoken ( token, iret )
         if ( iret==0 ) then
              ODSmean_fntmpl = trim(token)
              print *, ' Template of ODS mean filename: ', trim(ODSmean_fntmpl)
         else
         call gritas_perror ('GrITAS: cannot read ODSmean filename template in rc file' )
         endif
      endif

      ichan_version = 1
      call i90_label ( 'GrITAS*ODS_Channel_Version:', iret )
      if ( iret .eq. 0 ) then
           ichan_version = i90_gint(iret)
           if(iret/=0) call gritas_perror ('GrITAS: missing entry Channel Version')
      endif
      print *
      if(ichan_version==1) then
         print *, 'Expects ODS level to carry indexes of satellite data '
      else if (ichan_version==2) then
         print *, 'Expects ODS level to carry channels of satellite data '
      else
         call gritas_perror ('GrITAS: invalid index option' )
      endif
      

!     Vertical levels or satellite channels.  Scan the list
!     twice.  The first scan to determine the list size, and
!     the second to extract the list levels or channels.
!     ------------------------------------------------------
      do iscan = 1, 2
         if ( iscan .eq. 2 ) then
            print *, 'Found ', nlevs, ' levels to bin to'
            allocate ( plevs2(nlevs), stat=iret )
            if(iret .ne. 0) call gritas_perror
     &            ('GrITAS: Allocation error for plevs2 ')
         end if

!        vertical levels
!        ---------------
         call i90_label ( 'GrITAS*Vertical_Levels:', iret )
         if ( iret .ne. 0 ) call gritas_perror
     &       ('GrITAS: cannot find vertical levels in the rc file ' )
 
         nlevs = 0
         do i = 1, 32767

!           get token
!           ---------
            call i90_gtoken ( token,  iret )
            if ( iret .ne. 0 ) go to 222

!           ... and process it
!           ------------------
	    ic = index(token,':')
            if (ic .eq. 0) then ! ... as a single vertical level
                                ! ------------------------------
               p1 = i90_AtoF ( token, iret )
               if ( iret .ne. 0 ) call gritas_perror
     &            ('GrITAS: Bad token for vertical level ' //
     &             'which must be a float or integer' )
               p2 = p1
            else                ! ... or as a range
                                ! -----------------
               p1 = i90_AtoF (token(1:ic-1), iret)
               if ( iret .ne. 0 ) call gritas_perror
     &            ('GrITAS: Bad value for vertical level' //
     &             'which must be a float or integer')
               p2 = i90_AtoF (token(ic+1:),  iret)
               if ( iret .ne. 0 ) call gritas_perror
     &            ('GrITAS: Bad value for vertical level' //
     &             'which must be a float or integer')
            end if
            incr = 1
            nlevs_add = incr + floor (( p2 - p1 + 0.0001 * incr) / incr )
            do ilev_add = 1, nlevs_add
               nlevs = nlevs + 1
               if ( iscan .eq. 2 )
     &            plevs2 ( nlevs ) = p1 + ( ilev_add - 1 ) * incr
            end do
         end do
 222     continue
         if ( nlevs .eq. 0 )
     &         call gritas_perror('GrITAS: no vertical levels or ' //
     &                             'satellite channels')
      end do

      if ( plevs2(1) < plevs2(2) ) plevs2(1:nlevs) = plevs2(nlevs:1:-1)
      print *
      print *, 'Vertical Levels: ', ( plevs2(i), i=1,nlevs )

!     Vertical levels or satellite channels.  Scan the list
!     twice.  The first scan to determine the list size, and
!     the second to extract the list levels or channels.
!     ------------------------------------------------------
      do iscan = 1, 2
         if ( iscan .eq. 2 ) then
            allocate ( achan(nchan,2), stat=iret )   !  columns: (1) chn index; (2) chn number 
            if(iret .ne. 0) call gritas_perror
     &            ('GrITAS: Allocation error for achan ')
         end if

!        vertical levels
!        ---------------
         nchan = 0
         call i90_label ( 'GrITAS*Satellite*Info::', iret )
         if ( iret .ne. 0 ) exit ! don't do anything when list of active channels in rc file
 
         do i = 1, 32767

            call i90_gline ( iret )
            if ( iret /=0 ) exit ! end of table

!           determine whether channel is used or not
!           ----------------------------------------
            if ( iscan .eq. 1 ) then
               idummy = i90_gint(iret)
                  if(iret/=0) call gritas_perror ('GrITAS: missing 1st entry in sat table')
               idummy = i90_gint(iret)
                  if(iret/=0) call gritas_perror ('GrITAS: missing 2nd entry in sat table')
               idummy = i90_gint(iret)
                  if(iret/=0) call gritas_perror ('GrITAS: missing 3rd entry in sat table')
               if(idummy==1) nchan=nchan+1
            endif

!           read channel index and channel number
!           -------------------------------------
            if ( iscan .eq. 2 ) then
               ichan_indx = i90_gint(iret)
               ichan_numb = i90_gint(iret)
               idummy     = i90_gint(iret)
               if(idummy==1) then
                  nchan=nchan+1
                  achan ( nchan,1 ) = ichan_indx
                  achan ( nchan,2 ) = ichan_numb
               endif
            endif
         end do
 223     continue
         if ( nchan .eq. 0 )
     &        call gritas_perror('GrITAS: satellite channels not found ' )
      end do

      if (nchan>0) then
          print *
          print *, 'Active Channels Index: '
          write(6,'(10(2x,i4))') ( achan(i,1), i=1,nchan )
          print *, 'Active Channels: '
          write(6,'(10(2x,i4))') ( achan(i,2), i=1,nchan )
      endif

!     Initialize grid
!     ---------------
      km = nlevs
      print *
      print *, '      Resolution: ', im, jm, km
      print *
      call grids_init ( im, jm, km, iret)
      if(iret/=0) call die(myname,'error from grids_init')
      plevs = plevs2
      deallocate ( plevs2 ) ! clean up

      if(nchan==0) call GrITAS_HydroHeights()

!     Initialize list
!     ---------------
      listsz = 0
      do i = 1, LMAX
         do j = 1, KXMAX
            kx_list(j,i) = -1
         end do
      end do

!     kx-kt list
!     ----------
      long_name_found = .false.
      kx_list = -1
      call i90_label ( 'GrITAS*List::', iret )
      if ( iret .ne. 0 ) call
     &   gritas_perror('GrITAS: cannot find kx-kt table on rc file')
      do 10 i = 1, 32767
         call i90_gline ( iret )
         if ( iret .eq. -1 ) 
     &        call gritas_perror ('GrITAS: premature end of table')
         if ( iret .eq. 1 ) go to 11
         listsz = listsz + 1
         if(listsz .gt. LMAX )
     &      call gritas_perror('GrITAS: kx-kt list is too large')
         kt_list(listsz) = i90_gint(iret)
         if (iret.ne.0) call gritas_perror('GrITAS: table error: kt')
         call i90_gtoken ( argv, iret )
         if (iret.ne.0) call gritas_perror('GrITAS: table error: var')
         isfckt_mkr = scan ( argv, sfckt_mkr )
         if ( isfckt_mkr .eq. 1 ) then
            kt_list(listsz) = -1 * kt_list(listsz)
            var_list(listsz) = argv(isfckt_mkr+1:)
         else
            var_list(listsz) = argv
         end if
         var_descr(listsz) = ' '

         j = 0
         do jj = 1, KXMAX
            call I90_GToken(token, iret )
            if (iret/=0) then      ! read error
                exit
            end if
               ii = index(token,':') ! token is single kx or range of kx's
               lt = len_trim(token)
               if (ii==0) then    ! no colon, therefore single kx
                  read(token,*) kx1
                  kx2 = kx1
               else                  ! colon, therefore kx1:kx2
                  read(token(1:ii-1),*) kx1
                  read(token(ii+1:lt),*) kx2
               end if
               if (kx1>kx2) then      ! range error
                  write(stderr,'(2a,i5)')
     &               myname, ': Invalid range: ', token
                  call die(myname)
               end if
               do kxnext = kx1, kx2
                  if (j==KXMAX) then    ! check space
                     write(stderr,'(2a,i5)') myname,': increase KXMAX'
                     call die(myname)
                  else if (iret==0) then
                     j = j + 1
                     kx_list(j,listsz) = kxnext
                  end if
               end do
         end do ! while

 21      continue
 10   continue
 11   continue

      if ( .not. long_name_found ) then
         do ilist = 1, listsz
            var_descr ( ilist ) = var_list ( ilist )
         end do
      end if

!     Echo the parameters
!     -------------------
      print *
      print *, 'Kx-Kt table...'
      do i = 1, listsz

         nkx = 0
         do j = 1, KXMAX
            if ( kx_list(j,i) .lt. 0 ) go to 51
            nkx = j
         end do
 51      continue


         print *
         print *, 'Variable ', i, ': ', trim(var_list(i))
         print *, '   o kt = ', abs(kt_list(i))
         print *, '   o kx = ', (kx_list(j,i),j=1,nkx)

      end do

      if ( confidence > 0.0 ) then
         print *
         print *, 'Std-dev confidence level: ', confidence
         confidence = 1.0-confidence ! this is what the cdf routine needs
      endif

!     Release resources
!     -----------------
      call I90_Release()

!     Create land-water-ice mask if required
!     --------------------------------------
      if(trim(mask_fname)/='NONE') then
         call gritas_ini_masks(trim(mask_fname),im,jm,lwimask)
      endif

!     All done
!     --------
      return
      end

!...........................................................................
      subroutine GrITAS_Clean ()
      use m_gritas_grids, only : grids_clean
      use m_gritas_masks, only: gritas_fnl_masks
      use m_die
      implicit none
      integer ier
      call gritas_fnl_masks()
      call grids_clean(ier)
        if(ier/=0) call die('GrITAS_Clean','error from grids_clean')
      end subroutine GrITAS_Clean
!...........................................................................

      subroutine gritas_usage()
      print *
      print *,'Usage:  gritas.x  [-rc fn] [-del]  [-oma|-obs]  [-ieee] [-reduce]'
      print *,'                  [-rcred fn] [-o fn] [-ncf] [-res RES] obs_file(s)'
      print *
      print *,  'where'
      print *
      print *, '-rc fname   resource file name (default: gritas.rc)'
      print *, '-conf NUM   confidence interval for stdev (e.g.: 0.95)'
      print *, '-del        input files are GEOS-1 DEL files'
      print *, '             (by default input files are ODS)'
      print *, '-nopassive  passive data types will not be included'
      print *, '-oma        produces gridded O-A residuals instead'
      print *, '              of O-F residuals'
      print *, '-obs        produces gridded observed values instead'
      print *, '              of O-F residuals'
      print *, '-obias      produces gridded observed bias values instead'
      print *, '              of O-F residuals'
      print *, '-obsimp     produces observation impacts'
      print *, '-sigo       produces observation error values'
      print *, '-esigo      produces estimate of sigo from residuals'
      print *, '-hbh        produces gridded O-A x O-B residual statistics'
      print *, '-hah        produces gridded A-F x O-A residual statistics'
      print *, '-job        produces gridded Jo(O-F)/p'
      print *, '-joa        produces gridded Jo(O-A)/p'
      print *, '-oobs       produces original obs (not bias corrected)'
      print *, '-fobloc     produces forecast (background) at ob location'
      print *, '            (see note 3)'
      print *, '-promf      produces prescribed ensemble innovation covariance'
      print *, '            at ob location (see note 4)'
      print *, '-ieee       output file is in IEEE binary format'
      print *, '              instead of HDF/GFIO format; both formats'
      print *, '              are GrADS compatible'
      print *, '-reduce     to apply reducer (default RC file: reducer.rc)'
      print *, '-resbysigo  scale residuals by sqrt(R)                    '
      print *, '-resbyf     scale residuals by forecast (good for GPS eval)'
      print *, '-resbyo     scale residuals by observations (good for GPS eval)'
      print *, '-rcred fn   specify alternative name for reducer resource file'
      print *, '-ncf        non-compliant format (apply for GSI-diag-conv files)'
      print *, '-nonorm     skips normalization of accumulated summs (good for Obs Impact)'
      print *, '-res RES    resolution definition (X=a, b, c, d, e or m,for merra) '
      print *, '            (default: a)'
      print *, '-nlevs km   Number of vertical levels'
      print *, '            (default: 26)'
      print *, '-mask mask.nc4 WHICH  allow for using lwi mask when calculating xcov'
      print *, '            mask.nc4 is a file containg fractions; WHICH is ocean/land/ice'
      print *, '-xcov       Estimate cross-covariance/correlations <oma,omb>'
      print *, '-date nymd  Date for output file; default is mean of input'
      print *, '-ospl       Split bias, stdv and nobs into separate files'
      print *, '            (only for sdf files)'
      print *, '-lb         Write tally of obs not included in'
      print *, '            computation of statistics'
      print *, '-o bname    specifies output file base name, optional '
      print *, '            (default: gritas)'
      print *, '-nosym      do not symmetrize x-cov estimate of cov (for test only)'
      print *, '-raobspd    construct RAOB speed'
      print *
      print *, 'obs_files   ODS, DEL, or DIAG_CONV file names (no default)'
      print *
      print *
      print *, ' NOTES: '
      print *
      write(6,*) ' 1) In case of calculating statistics of (specified) sigo and obias the standard deviation '
      write(6,*) '    output is replaced by root-mean-square; that is to say that position t=1 in the grads '
      write(6,*) '    file will have the mean as it normally has, and t=2 will have the RMS instead of the STDV. '
      print *
      write(6,*)  ' 2) In case of Desroziers et al. diagnostics, the mean and standard deviations output '
      write(6,*) '     are replaced with cross-variance and the cross-standard-deviations, respectively. '
      print *
      write(6,*)  ' 3) The knob -fobloc can be combined with an rc-file specific parameter to calculate'
      write(6,*)  '    the ensmeble spread = f_i - <f>, where f=bkg/fcst in observation space.'
      print *
      write(6,*)  ' 4) The knob -promf must use an rc file with specific parameter to calculate'
      write(6,*)  '    ensemble spread; Bkg error cov estimated from ensemble spread.'
      print *
      write(6,*)  ' 5) The knob -xcov combined with -esigo allows for estimation of full <oma,omb> '
      write(6,*)  '    per channel (level); this can be used to infer channels cross-correlations.'
      print *
      print *, 'Revised: 25 Apr 2014 (R. Todling)'
      print *, 'Created:  1 September 1997 (A. da Silva)'
      print *
      call gritas_perror('GrITAS: not enough input parameters')
      end

!EOC




