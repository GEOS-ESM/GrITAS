!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOI
!
! !TITLE: An Application for Gridding Innovation\\ Synoptic Averaged Statistics (GrISAS v1.00)
!
! !AUTHORS: Arlindo da Silva 
!
! !AFFILIATION: Data Assimilation Office, NASA/GSFC, Greenbelt, MD 20771
!
! !DATE: September 1, 1997 (Revised June 9, 2004)
!
! !INTRODUCTION: System Overview
!
!  GrISAS is a FORTRAN 77 program which reads innovation (observation minus
!  forecast residuals) files and produces global grids with time mean,
!  standard deviation and number of observations in each grid box,
!  for a given synoptic time.
!
!  Grids are produced for user specified combination of data types
!  (e.g., u, v-winds, heights) and data sources (e.g., radiosondes,
!  TOVS retrievals, ships, etc.) 
!
!  iopt is an option parameter that let the user decide what kind of
!  o-f are included in the statistics.
!  Currently, iopt = 1 is hardwired in this code as we focus on those
!  observations which have passed the on-line quality control.
!
!EOI
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GrISAS()
! 
! !DESCRIPTION: Driver for the {\em Gridding of Innovation Synoptic Average
!               Statistics} system (GrISAS).
!
! !INTERFACE:
!
       Program GrISAS

! !USES:

       Use m_MergeSorts

       Use m_odsmeta, only : KTMAX
       Use m_odsmeta, only : KXMAX
       Use m_odsmeta, only : X_PASSIVE  ! PSAS passive data

       Use m_gritas_grids, only : im, jm, km
       Use m_gritas_grids, only : GrITAS_Pint
       Use m_gritas_grids, only : GrITAS_Accum
       Use m_gritas_grids, only : GrITAS_Norm

       Use m_gfio_output, only : GFIO_Output
       Use m_ods

!      Use m_obs, only : Reducer
       Use m_gritas_grids, only : GrITAS_CDF

       Use m_StrTemplate

       Use m_stdio, only : stdout
       Use m_die
!
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
!  26Mar98  da Silva   Derived from GrITAS.
!  15Dec98  da Silva   Added ODS/GFIO output.
!  02feb99  da Silva   Fixed missing obs bug.
!  23dec00  da Silva   Made part of fvDAS.
!  08Jun04  Todling    - Converted to use of m_ods.
!                      - Resolution independent executable 
!  24Jun04  Todling    Fixed time frequency of output (RUC compatible).
!  30Jul04  Todling    Added reducer capability
!  12oct05  da Silva   Added ntresh to impose a minimum number of obs
!  07Nov06  Todling    Implemented Desroziers et al. diagnostics
!  23Jan08  Todling    Added -obsimp
!  07Jul09  Redder     Added ob_flag as an argument in the call of the
!                      routine GrITAS_Accum
!  10May13  Todling    - Add -fobloc to calc fcst at ob-location
!                      - Add opt to read extra ODS file and calc ens spread
!
!EOP
!-------------------------------------------------------------------------
!BOC
      Implicit NONE

      character(len=*), parameter :: myname = 'GrISAS'

      include 'gritas.h'


!     User defined vertical levels
!     ----------------------------
      real, allocatable ::  p1(:), p2(:) ! pressure intervals for gridding


!     Workspace for Holding observations for each synoptic time
!     Note: lat/lon/etc cannot be pointers because ods is r8 vs lat_r4
!     ----------------------------i-----------------------------------
      integer                  nobs      ! actual number of obs for this synoptic time
      integer                  nobs_good ! no. of observations after reducer operation
      type    ( ods_vect )     ods
      type    ( ods_vect )     mods
      integer, allocatable ::  kx(:)
      integer, allocatable ::  kt(:)
      real,    allocatable ::  lat(:)
      real,    allocatable ::  lon(:)
      real,    allocatable ::  lev(:)
      real,    allocatable ::  qc(:)   ! not used, but needed
      real,    allocatable ::  del(:,:)
      integer, allocatable ::  is(:)           ! ordering index
      integer, allocatable ::  ob_flag (:)     ! = 0 if ob was used to assumulate
                                               !   bias/RMS, = nonzero otherwise

!     Workspace for holding the surface gridded fields
!     ------------------------------------------------
      real, allocatable ::  bias_2d(:,:,:,:)   ! time means
      real, allocatable ::  rtms_2d(:,:,:,:)   ! root time mean square
      real, allocatable ::  stdv_2d(:,:,:,:)   ! standard deviation
      real, allocatable ::  xcov_2d(:,:,:,:)   ! dummy in grisas
      real, allocatable ::  nobs_2d(:,:,:)     ! number of obs per grid box
      real, allocatable ::  chsq_2d(:,:,:,:)   ! chi-square for stddev confidence
      integer list_2d(l2d_max)                 ! corresponding list item

!      equivalence ( rtms_2d(1,1,1), stdv_2d(1,1,1) )


!     Workspace for holding the upper-air gridded fields
!     --------------------------------------------------
      real, allocatable ::  bias_3d(:,:,:,:,:) ! time means
      real, allocatable ::  rtms_3d(:,:,:,:,:) ! root time mean square
      real, allocatable ::  stdv_3d(:,:,:,:,:) ! standard deviation
      real, allocatable ::  xcov_3d(:,:,:,:,:) ! dummy in grisas
      real, allocatable ::  chsq_3d(:,:,:,:,:) ! chi-square for stddev confidence
      real, allocatable ::  nobs_3d(:,:,:,:)   ! number of obs per grid box
      integer list_3d(l3d_max)                 ! corresponding list item



!     Data structure for user selection of data type/sources
!     ------------------------------------------------------
      integer       listsz               ! size of the list (< lmax)
      integer       kt_list(lmax)        ! each grid has one data type (kt)  
      integer       kx_list(kxmax,lmax)  ! but several data sources (kx)
      character*11  var_list(lmax)       ! variable name for output 
                                         ! (e.g., uraob)
      character*256 var_descr(lmax)      ! ... and its description.
!     Note: Negative kt values will be treated as surface variables

!     Data type/data source table. Gives the grid number for each (kx,kt)
!     This table is derived from the user lists above. Notice that 
!     surface grids are associated with kt<4.
!     -------------------------------------------------------------------
      integer       kxkt_table(kxmax,ktmax)

      real confidence
      integer         l2d                      ! actual number of 2d grids
      integer         l3d                      ! actual number of 3d grids

!     Input File names
!     ----------------
      integer          ndel                      ! actual number of del-files
      character*(255)  del_ifname(NDEL_MAX)      ! O-F (del) data files
      character*(255)  odsmean_fntmpl            ! template of mean of ods files (for EnADAS only)

!     Output File names
!     -----------------
      character*(255)  grid_obasen               ! gridded O-F (del) base name
      character*(255)  grid_ofname               ! gridded O-F (del) file


!     Local variables
!     ---------------
      integer,pointer :: fid(:)
      integer idel, isyn, ksyn, ioptn, iopt, ios, i, j, iii
      integer nymd, nhms, timinc, nsyn, ier, nstat, nr, rc
      character*8 ftype, dtype, otype
      character(len=255) RC_red, fname
      logical DO_red
      logical ncf, split_tgroups, first, nrmlz, meanres
      integer ntresh
      integer nymd_save, nhms_save, imiss, nmiss

      integer scale_res ! <=0 do not scale
                        ! 1 = scale by sigo
                        ! 2 = scale by o
                        ! 3 = scale by f

!...........................................................................
      
      ioptn = 1                          ! keep only obs which passed QC
      first = .true.
      nrmlz = .true.

!     User interface
!     --------------
      call GrISAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                   ftype, dtype, otype, ioptn, ncf, nr,
     &                   split_tgroups,
     &                   DO_red, RC_red,
     &                   ndel, del_ifname, grid_obasen, odsmean_fntmpl,
     &                   kt_list, kx_list, var_list, var_descr, 
     &                   listsz, ntresh, nrmlz, confidence, scale_res ) 
      iopt = abs(ioptn)

      if (confidence<0.0) then
          nstat = 3     ! writes out mean, std, and nobs count
      else
          nstat = 6     ! writes out mean, std, nobs count, stdv_chisqr
      endif
      if (split_tgroups) then
          allocate(fid(nstat))
          nstat = -1 * nstat
      else
          allocate(fid(1))
      endif
      fid=-1

!     Create kx-kt table from user defined list
!     -----------------------------------------
      call GrITAS_KxKt ( KTMAX, KXMAX, L2D_MAX, L3D_MAX,
     &                   listsz, kt_list, kx_list, 
     &                   kxkt_table, list_2d, list_3d, l2d, l3d )


!     Allocate space for working arrays
!     ---------------------------------
      allocate ( bias_2d(im,jm,l2d,nr), rtms_2d(im,jm,l2d,nr),
     &           stdv_2d(im,jm,l2d,nr), nobs_2d(im,jm,l2d),
     &           xcov_2d(im,jm,l2d,nr), chsq_2d(im,jm,l2d,3),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(2d)')

      allocate ( bias_3d(im,jm,km,l3d,nr), rtms_3d(im,jm,km,l3d,nr),
     &           stdv_3d(im,jm,km,l3d,nr), nobs_3d(im,jm,km,l3d), 
     &           xcov_3d(im,jm,km,l3d,nr), chsq_3d(im,jm,km,l3d,3),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(3d)')

      allocate ( p1(km), p2(km),
     &           stat=ier )
        if(ier/=0) call die (myname,'Error in Alloc(pressure)')


!     Set pressure Intervals for gridding
!     -----------------------------------
      call GrITAS_Pint ( p1, p2 )

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      if(ncf) ksyn = 1

!     Loop over del files
!     -------------------
      do idel = 1, ndel

!         Loop over synoptic time on this del file
!         ----------------------------------------
          do 10 isyn = 1, ksyn     
             meanres = .false.


!            Initialize grids with zeros
!            ---------------------------
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


!            Read del file for this synoptic time
!            ------------------------------------
             if ( trim(ftype) .eq. 'del' ) then

                allocate ( lat(NOBSMAX), lon(NOBSMAX), lev(NOBSMAX), 
     &                     qc (NOBSMAX), del(NOBSMAX,nr),
     &                     kx(NOBSMAX), kt(NOBSMAX), ob_flag(NOBSMAX),
     &                     stat=ier )
                   if(ier/=ier) call die (myname,
     &                                  ' Obs Arrays Error, alloc()')

                call Read_Del ( trim(del_ifname(idel)), NOBSMAX, ioptn,
     &               lat, lon, lev, del(1,1), kx, kt, qc, nobs,
     &               nymd, nhms, ier )
                ob_flag ( : nobs ) = 0

             else

                nymd = -1           ! get data for next synoptic time on file
                nhms =  0
                call ODSNxTime ( trim(del_ifname(idel)), nymd, nhms )
                if ( nymd .eq. -1 ) then
!                  print *, 'End-Of-File'
                   exit      ! skip to next file
                end if

!               PRC
!               Hard check for now...we want to ensure we write to the file
!               even if times are missing in the ods
!               Conditions: have read past end of file
!               without finding last synoptic time
!               expected...

!_RT            if ( nymd .eq. -1 ) then
!_RT              nmiss = (240000-nhms_save) / (240000/nsyn)
!_RT              print *, 'PETE: ',nymd, nhms_save, 240000/nsyn, nmiss
!_RT              do imiss = 1,nmiss
!_RT                 nymd = nymd_save
!_RT                 print *, isyn, nymd_save, nhms_save
!_RT                 nhms = nhms_save + 240000/nsyn
!_RT                 nhms_save = nhms
!_RT                 if(nhms .ge. 240000) exit  ! really are past end!
!_RT                 nobs = 0
!_RT                 print *, nhms_save, nhms, nsyn, isyn, 240000/nsyn
!_RT                 print *, 'End-Of-File'

!<-- PRC

!           Normalize grids
!           ---------------

!_RT        call GrITAS_Norm ( nr,nrmlz,
!_RT &                         l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
!_RT &                         l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,
!_RT &                         ntresh = ntresh )

!           if (nobs > 0) then

!            Write del grids to file
!            -----------------------
!_RT         nstat = 1               ! writes only means for now
!_RT         timinc = 240000 / nsyn  ! frequency of output
!_RT         if ( trim(otype) .eq. 'hdf' .or.  trim(otype) .eq. 'nc4' ) then
!_RT           print *,' PETE grisas: nymd,nhms ',nymd,nhms
!_RT
!_RT           call GFIO_Output  ( grid_obasen,otype,
!_RT &                             var_list, var_descr, listsz, timinc,
!_RT &                             l2d, list_2d, bias_2d(1,1,1,1),   stdv_2d(1,1,1,1),   nobs_2d,
!_RT &                             l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d,
!_RT &                             nstat, nymd, nhms, fid )
!_RT          else
!_RT            call Make_Fname ( grid_obasen, nymd, nhms, grid_ofname )
!_RT            call GrADS_Output ( grid_ofname,
!_RT &                             var_list, listsz,
!_RT &                             l2d, list_2d, bias_2d(1,1,1,1),   stdv_2d(1,1,1,1),   nobs_2d,
!_RT &                             l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d, nstat )
!_RT           end if

! --> PRC
!         end do
!               exit     ! skip to next file
!_RT      else
!_RT         nymd_save = nymd
!_RT         nhms_save = nhms
!_RT      end if

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
                         print *, 'Failed to read ODS with mean residuals, ODS_Get error: ier = ', ier
                         call exit(1)
                    else
                         write(stdout,'(3a,i8,a,i6.6)') 
     &                        'Read ', trim(fname), ' on ', nymd, ' at ', nhms
                    endif
                    meanres = .true.
                endif

!               Set number of observations found in the file
!               --------------------------------------------
                nobs   =  ods%data%nobs
                nsyn   = ods%meta%nsyn

!
                if ( nobs .eq. 0 ) then
                     print *, 'No data for this synoptic time'

!_RT                 nstat = 1

!_RT                 if(nsyn > 8 .or. nsyn < 1) then
!_RT                   nsyn = 4
!_RT                 endif

!_RT                 timinc = 240000 / nsyn  ! frequency of output
!_RT
!_RT                 call GrITAS_Norm ( nr,nrmlz,
!_RT &                             l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
!_RT &                             l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,
!_RT &                             ntresh = ntresh )

!_RT                 if ( trim(otype) .eq. 'hdf' .or. trim(otype) .eq. 'nc4' ) then
!_RT                    print *,' grisas: nymd,nhms,timinc ',nymd,nhms,timinc
!_RT                    call GFIO_Output  ( grid_obasen,otype,
!_RT &                                      var_list, var_descr, listsz, timinc,
!_RT &                                      l2d, list_2d, bias_2d(1,1,1,1),   stdv_2d(1,1,1,1),   nobs_2d,
!_RT &                                      l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d,
!_RT &                                      nstat, nymd, nhms, fid )
    
!_RT                 else
!_RT                   call Make_Fname ( grid_obasen, nymd, nhms, grid_ofname )
!_RT                   call GrADS_Output ( grid_ofname,
!_RT &                                     var_list, listsz,
!_RT &                                     l2d, list_2d, bias_2d(1,1,1,1),   stdv_2d(1,1,1,1),   nobs_2d,
!_RT &                                     l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d, nstat )
!_RT                 end if


                     cycle    ! skip to next synoptic time
                end if
                                                                                                                    
!               Reduce input data as requested
!               ------------------------------
!               if ( DO_red ) then
!                    call Reducer ( nymd, nhms, nobs, ods, nobs_good, rcfile=RC_red )
!                    ods%data%nobs = nobs_good
!                    nobs          = nobs_good
!               end if

                allocate ( lat(nobs), lon(nobs), lev(nobs),
     &                     qc (nobs), del(nobs,nr),  kx(nobs), kt (nobs),
     &                     ob_flag (nobs), stat=ier )
                   if(ier/=ier) call die (myname,
     &                                  ' Obs Arrays Error, alloc()')

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
                qc      = real(ods%data%qcexcl(1:nobs))
                ob_flag ( : nobs ) = 0

                deallocate(is)

!               Force y2K compliance
!               --------------------
                if ( nymd .lt. 19000000 ) nymd = nymd + 19000000 
                nsyn = ods%meta%nsyn

!               Observed value (O-F, O-A, or obs)
! 	        ---------------------------------
                if ( iopt .eq. 0 .or. iopt .eq. 1 ) then                      ! obs minus background residuals
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
                else if ( iopt .eq. 2 .or. iopt .eq. 3 ) then                 ! obs minus analysis residuals
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
                else if ( iopt .eq. 4 .or. iopt .eq. 5 ) then                 ! observations
                     del(:,1) = ods%data%obs(1:nobs)
                else if ( iopt .eq. 6 .or. iopt .eq. 7 ) then                 ! prescribed sigo
                     del(:,1) = ods%data%xvec(1:nobs)
                else if ( iopt .eq. 8 .or. iopt .eq. 9 ) then                 ! obs bias
                     del(:,1) = ods%data%xm(1:nobs)
                else if ( iopt .eq.10 .or. iopt .eq.11 ) then                 ! X-cov <AmF,OmF> ~ HBH'
                     del(:,1) = ods%data%oma(1:nobs) - ods%data%omf(1:nobs)
                     del(:,2) = ods%data%omf(1:nobs)
                else if ( iopt .eq.12 .or. iopt .eq.13 ) then                 ! X-cov <AmF,OmA> ~ HPaH'
                     del(:,1) = ods%data%oma(1:nobs) - ods%data%omf(1:nobs)
                     del(:,2) = ods%data%oma(1:nobs)
                else if ( iopt .eq.14 .or. iopt .eq.15 ) then                 ! X-cov <OmA,OmF> ~ R
                     del(:,1) = ods%data%oma(1:nobs)
                     del(:,2) = ods%data%omf(1:nobs)
                else if ( iopt .eq.16 .or. iopt .eq.17 ) then                 ! Jo(OmF)/p
                     del(:,1) = ods%data%omf(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                else if ( iopt .eq.18 .or. iopt .eq.19 ) then                 ! Jo(OmA)/p
                     del(:,1) = ods%data%oma(1:nobs)
                     del(:,2) = ods%data%xvec(1:nobs)
                else if ( iopt .eq. 20 .or. iopt .eq. 21 ) then               ! original observations
                     del(:,1) = ods%data%obs(1:nobs) + ods%data%xm(1:nobs)
                else if ( iopt .eq. 22 .or. iopt .eq. 23 ) then               !  observation impacts
                     del(:,1) = ods%data%xvec(1:nobs)
                else if ( iopt.eq.24 .or. iopt.eq.25 ) then      !  background on obs location
                     if (meanres) then ! del = f_i - <f>
                         do iii=1,nobs
                           if(ods%data%obs(iii)==mods%data%obs(iii)) then
!                             del(iii,1) = mods%data%obs(iii) - mods%data%omf(iii)
!                             del(iii,1) =  ods%data%obs(iii) - ods%data%omf(iii) - del(iii,1)
                              del(iii,1) = mods%data%omf(iii) - ods%data%omf(iii)
                           else
                              ods%data%qcexcl(iii) = X_PASSIVE
                           endif
                        enddo
                     else  ! del = f_i
                         del(:,1) = ods%data%obs(1:nobs) - ods%data%omf(1:nobs)
                     endif
                else
                     call die (myname, 'iopt not valid')
                endif

!               Now that all atrributes have been read, eliminate obs
!                which have bot passwd QC, if user desires. Notice that
!                passive data types (qcexcl=7) are kept if ioptn>0.
!               ------------------------------------------------------
                if ( mod(iopt,2) .ne. 0 ) then
                    do i = 1, nobs
                       if ( ods%data%qcexcl(i) .eq. 0 .or.
     &                     ( ioptn.gt.0 .and. ods%data%qcexcl(i) .eq. 7 ) ) then
                          ob_flag ( i ) = 0
                       else 
                          ob_flag ( i ) = 1
                       end if
                    end do
                    if(iopt==17 .or. iopt==19) del(1:nobs,1) = del(1:nobs,1)/del(1:nobs,2)  ! omf/sigo for Joa/Job calculation
               end if
             endif

!            Accumulate grids for this synoptic time
!            ---------------------------------------
             if ( nobs .gt. 0 ) then
                call GrITAS_Accum ( lat, lon, lev, kx, kt, del, nobs, nr,
     &                              kxkt_table, p1, p2,
     &                              ob_flag,
     &                              l2d, bias_2d, rtms_2d, xcov_2d(:,:,:,1),   nobs_2d,
     &                              l3d, bias_3d, rtms_3d, xcov_3d(:,:,:,:,1), nobs_3d )
             end if

            deallocate ( lat, lon, lev, qc, del, kx, kt, ob_flag )

!           Calculate Chi-square for specified confidence interval
!           ------------------------------------------------------      
            if (confidence>0.0) then
               call GrITAS_CDF ( 'chisq-dist', 1.0-0.5*confidence, nobs_2d, chsq_2d(:,:,:,1), nobs_3d, chsq_3d(:,:,:,:,1) )
               call GrITAS_CDF ( 'chisq-dist',     0.5*confidence, nobs_2d, chsq_2d(:,:,:,2), nobs_3d, chsq_3d(:,:,:,:,2) )
               call GrITAS_CDF (     't-dist',     0.5*confidence, nobs_2d, chsq_2d(:,:,:,3), nobs_3d, chsq_3d(:,:,:,:,3) )
            endif

!           Normalize grids
!           ---------------
            if ( nrmlz ) then
                 call GrITAS_Norm ( nr, nrmlz,
     &                             l2d, bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,
     &                             l3d, bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,
     &                             ntresh = ntresh )
            endif


!           In case of x-covariances, write x-rms in the bias slot, and x-variances in stdv slot
!           NOTE: in grisas, only means matter.
!           ------------------------------------------------------------------------------------
            if ( iopt==10 .or.iopt==11 .or. iopt==12 .or.iopt==13 .or. iopt==14 .or.iopt==15 ) then
                 bias_2d(:,:,:,1)   = xcov_2d(:,:,:,1)
                 bias_3d(:,:,:,:,1) = xcov_3d(:,:,:,:,1)
            endif

!           In case of Jo(OmA) and Jo(OmF) store Jo in bias array
!           -----------------------------------------------------
            if ( iopt==16 .or.iopt==17 .or. iopt==18 .or.iopt==19 ) then
                 bias_2d(:,:,:,  1) = rtms_2d(:,:,:,  1)
                 bias_3d(:,:,:,:,1) = rtms_3d(:,:,:,:,1)
            endif


!           Write del grids to file 
!           -----------------------
            nstat = 1               ! writes only means for now
            timinc = 240000 / nsyn  ! frequency of output
            if ( trim(otype) .eq. 'hdf' .or. trim(otype) .eq. 'nc4' ) then
!              call Make_Fname ( grid_obasen, nymd, nhms, grid_ofname )
               call GFIO_Output  ( grid_obasen,otype, 
     &                             var_list, var_descr, listsz, timinc,
     &                             l2d, list_2d, bias_2d(:,:,:,1),   stdv_2d(:,:,:,1),   nobs_2d,
     &                             l3d, list_3d, bias_3d(:,:,:,:,1), stdv_3d(:,:,:,:,1), nobs_3d,
     &                             nstat, nymd, nhms, fid, chsq_2d=chsq_2d, chsq_3d=chsq_3d )

            else
               call Make_Fname ( grid_obasen, nymd, nhms, grid_ofname )
               call GrADS_Output ( grid_ofname, 
     &                             var_list, listsz,
     &                             l2d, list_2d, bias_2d(1,1,1,1),   stdv_2d(1,1,1,1),   nobs_2d,
     &                             l3d, list_3d, bias_3d(1,1,1,1,1), stdv_3d(1,1,1,1,1), nobs_3d, nstat )
            end if


!           Clean up mean-carrying ODS structure
!           ------------------------------------
            if (meanres) then
                call ods_clean ( mods, ier )
            endif

 10       continue
 11       continue


       end do

!      Clean up
!      --------
       deallocate ( bias_3d, rtms_3d, stdv_3d, xcov_3d, nobs_3d,chsq_3d, stat=ier )
         if(ier/=0) call die (myname,'Error in Alloc(3d)')

       deallocate ( bias_2d, rtms_2d, stdv_2d, xcov_2d, nobs_2d,chsq_2d,stat=ier )
         if(ier/=0) call die (myname,'Error in Dealloc(2d)')

       call GrITAS_Clean ()

!      All done.
!      --------
       if ( trim(otype) .eq. 'hdf' .or. trim(otype) .eq. 'nc4' ) then
          do iii=1,size(fid)
             call GFIO_Close ( fid(iii), rc )   ! close GFIO file
          enddo
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
! !ROUTINE:  GrISAS_Init --- Command line and resource user interface
! 
! !INTERFACE:
!
      subroutine GrISAS_Init ( NDEL_MAX, LMAX, KXMAX,
     &                         ftype, dtype, otype, iopt, ncf, nr,
     &                         split_tgroups,
     &                         DO_red, RC_red,
     &                         ndel, del_ifname, grid_ofname, odsmean_fntmpl,
     &                         kt_list, kx_list, var_list, var_descr, 
     &                         listsz, ntresh, nrmlz, confidence, scale_res )

! !USES:

      Use m_gritas_grids, only : grids_init
      Use m_gritas_grids, only : im, jm, km
      Use m_gritas_grids, only : nlevs
      Use m_gritas_grids, only : plevs
      Use m_gritas_grids, only : ichan_version
      Use m_gritas_grids, only : nchan
      Use m_gritas_grids, only : achan
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
      character(len=*) dtype                 ! Data type: 'omf', 'oma', 'obs', 'amf', 'sigo', 'obias'
      character(len=*) otype                 ! Output type: 'hdf','nc4', or 'eee'
      integer          nr                    ! number of distinct residual series
      integer          iopt                  ! =1 for O-F
                                             ! =3 for OMA 
                                             ! =5 for obs
      logical          ncf                   ! Non-compliant format (to allow
                                             !   handling of GSI-output files)
      logical         split_tgroups          ! = true to split the tgroups,
                                             !   bias, stdv and nobs into 
                                             !   separate files.

      integer         ndel                   ! actual number of del files
      character*255   del_ifname(NDEL_MAX)   ! del file names
      character*255   grid_ofname            ! output grid file name
      character*255   odsmean_fntmpl         ! template of filename with mean ODS (for EnADAS)

      integer         ntresh                 ! minimum number of obs to
                                             !  compote the mean
      integer         listsz                 ! actual list size (<LMAX)
      integer         kt_list(LMAX)          ! data type for the grid
      integer         kx_list(KXMAX,LMAX)    ! data sources for the grid
      character*11    var_list(LMAX)         ! variable name for the grid
      character*256   var_descr(LMAX)        ! ... and its description
   
      logical          DO_red                ! Defines whether to apply reducer or not
      character(len=*) RC_red                ! Recucer resource filename

! !INPUT/OUTPUT PARAMETERS:

       logical         nrmlz                 ! Controls normalization of accumalated
       real            confidence            ! Confidence level
       integer         scale_res             ! opt to scale residuals

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
!  18Jan08  Todling    Add -nonorm; bug fix: -reduce opt was unset
!  28May09  Redder     Added capability to define kt as a surface variable
!                      by added the prefix, sfckt_mkr, to its name in the
!                      input rc file
!  20Jul09  Redder     Added the subroutine output parameter, var_descr.
!  10May13  Todling    Add -fobloc; bug fix: -o was in conflict w/ other opts
!                      starting with "o"; add rc file opt for ens-spread
!  15Aug20  Redder/RT  Added the command line option, -ospl, and
!                      the subroutine output parameter, split_tgroups
!
!EOP
!-------------------------------------------------------------------------
!BOC

      character(len=*), parameter :: myname = 'GrISAS_Init'

!     Prefix for variable name denoting a surface variable
!     ----------------------------------------------------
      character (len=*), parameter :: sfckt_mkr = '_'

!     Temp storage for vertical levels
!     --------------------------------
      real, allocatable, save :: plevs2(:)

      integer        iarg, argc, iargc
      character*255  argv, rc_ifname, res, token
      integer i, j, nkx, lt, ii, jj, ilist, isfckt_mkr, ic
      integer iret, ios, iscan, idummy, ichan_indx, ichan_numb
      integer kx1, kx2, kxnext, ilev_add, nlevs_add, nchan_add
      real p, p1, p2, incr
      logical passive, long_name_found, long_name

      ntresh = 0          ! Bo treshold imposition
      passive = .true.
      split_tgroups = .false.
      ncf   = .false.
      ftype = 'ods'
      dtype = 'omf'
      otype = 'hdf'
      otype = 'nc4'
      iopt = 1
      nr   = 1
      im=72; jm=46; km=26  ! default resolution
      DO_red = .false.
      RC_red = 'reducer.rc'
      ODSmean_fntmpl = 'NONE'
      confidence = -999.

!     Parse command line
!     ------------------
      rc_ifname = 'gritas.rc'
      grid_ofname = 'grisas'      
      argc =  iargc()
      if ( argc .lt. 1 ) call gritas_usage()
      ndel = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if(trim(argv) .eq. '-del')  then
            ftype = 'del'
         else if(trim(argv) .eq. '-nopassive')  then
            passive = .false.
         else if(trim(argv) .eq. '-omf')  then
            dtype = 'omf'
            iopt  = 1
         else if(trim(argv) .eq. '-oma')  then
            dtype = 'oma'
            iopt  = 3
         else if(trim(argv) .eq. '-obs')  then
            dtype = 'obs'
            iopt  = 5
         else if(trim(argv) .eq. '-sigo')  then
            dtype = 'sigo'
            iopt  = 7
         else if(trim(argv) .eq. '-obias')  then
            dtype = 'obias'
            iopt  = 9
         else if(trim(argv) .eq. '-hbh')  then
            dtype = 'hbh'
            iopt  = 11
            nr    = 2
         else if(trim(argv) .eq. '-hah')  then
            dtype = 'hah'
            iopt  = 13
            nr    = 2
         else if(trim(argv) .eq. '-esigo')  then
            dtype = 'esigo'
            iopt  = 15
            nr    = 2
         else if(trim(argv) .eq. '-job') then
            dtype = 'job'
            iopt  = 17
            nr    = 2
         else if(trim(argv) .eq. '-joa')  then
            dtype = 'joa'
            iopt  = 19
            nr    = 2
         else if(trim(argv) .eq. '-oobs')  then
            dtype = 'oobs'
            iopt  = 21
         else if(trim(argv) .eq. '-obsimp')  then
            dtype = 'obsimp'
            iopt  = 23
         else if ( trim(argv) .eq. '-resbysigo' )  then
            scale_res = 1
         else if ( trim(argv) .eq. '-resbyo' )  then
            scale_res = 2
         else if ( trim(argv) .eq. '-resbyf' )  then
            scale_res = 3
         else if(trim(argv) .eq. '-fobloc')  then
            dtype = 'fobloc'
            iopt  = 25
         else if(trim(argv) .eq. '-reduce')  then
            DO_red = .true.
         else if(trim(argv) .eq. '-nonorm')  then
            nrmlz = .false.
         else if(trim(argv) .eq. '-reduce')  then
            DO_red = .true.
         else if(trim(argv) .eq. '-ncf')  then
            ncf = .true.
         else if(trim(argv) .eq. '-ieee')  then
            otype = 'eee'
         else if ( trim(argv) .eq. '-conf' ) then
            if ( iarg+1 .gt. argc ) call gritas_usage()
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) confidence
            if ( ios /= 0 ) then
                 print *, 'cannot parse confidence interval ...'
                 call exit(1)
            end if
         else if(trim(argv) .eq. '-o ') then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, grid_ofname )
         else if ( trim(argv) .eq. '-ospl' ) then
            split_tgroups = .true.
         else if(trim(argv) .eq. '-rc') then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, rc_ifname )
         else if(trim(argv) .eq. '-rcred') then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, RC_red )
         else if(trim(argv) .eq. '-ntresh') then
            if ( iarg+1 .gt. argc ) call gritas_usage()               
            iarg = iarg + 1
            call GetArg ( iArg, token )
            read(token,*,iostat=ios) ntresh
            if ( ios /= 0 ) then
                 print *, 'cannot parse ntresh ...'
                 call exit(1)
            end if
         else if(trim(argv) .eq. '-res') then
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
             else if ( res(1:1) == "D" ) then ! this is the MERRA resolution
                 im=540; jm=361
             else if ( res(1:1) == "e" ) then
                 im=1080; jm=721
             else
                 print *, 'unknon resolution ... taking default'
             endif
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
      if ( trim(ftype) .eq. 'del' .and. iopt .ne. 1 ) then
         print *
         print *, 'Note: DEL files only have O-F, ignoring -oma/-obs'
         dtype = 'omf'
         iopt = 1
      end if


      print *
      print *, 'Grid output file: ', trim(grid_ofname)
      print *

!     Read resource file
!     ------------------
      call i90_LoadF ( trim(rc_ifname), iret )
      if ( iret .ne. 0 ) then
         print *, 'GrISAS: cannot open resource file ',
     &            trim(rc_ifname)
         call exit(1)
      end if

!     specify template for files containing ensemble mean
!     ---------------------------------------------------
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
            nlevs_add = 1 + floor (( p2 - p1 + 0.0001 * incr) / incr )
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
         call i90_label ( 'GrITAS*Satellite*Info::', iret )
         if ( iret .ne. 0 ) exit ! don't do anything when list of active channels in rc file
 
         nchan = 0
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

!     Echo the parameters
!     -------------------
      km = nlevs
      print *
      print *, '      Resolution: ', im, jm, km
      print *, '   Residual type: ', trim(dtype)
      print *, 'Input  File type: ', trim(ftype)
      print *, 'Output File type: ', trim(otype)
      print *, '       Del files: ', ndel
      do i = 1, ndel
         print *, '     o ', trim(del_ifname(i))
      end do
      
!     Initialize grid
!     ---------------
      print *
      call grids_init ( im, jm, km, iret)
      if(iret/=0) call die(myname,'error from grids_init')
      plevs = plevs2
      deallocate ( plevs2 ) ! clean up


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
!            print *, trim(token)
            if (iret/=0) then      ! read error
                exit
!               call die(myname,'I90_GToken error, iret=',iret)
            end if
!           long_name = verify ( token, '-0123456780: ') .gt. 0 .and.
!    &                  jj .eq. 1
!           if ( long_name ) then
!              var_descr ( listsz ) =  token
!              long_name_found      = .true.
!           else
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
!           end if
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
         print *, '   o kt = ', kt_list(i)
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

      if ( split_tgroups ) then
         print *
         print *, ' Bias, stdv and nobs will be placed in the files '
         print *, ' with the basenames, '
         print *, '   ' // trim (grid_ofname) // '.bias,'
         print *, '   ' // trim (grid_ofname) // '.stdv and '
         print *, '   ' // trim (grid_ofname) // '.nobs'
         print *
      end if

!     All done
!     --------
      return
      end
        
!...........................................................................
      subroutine GrITAS_Clean ()
      use m_gritas_grids, only : grids_clean
      implicit none
      integer rc
      call grids_clean(rc)
      end subroutine GrITAS_Clean
!...........................................................................

      subroutine gritas_usage()
      print *
      print *,'Usage:  grisas.x  [-rc fn] [-del]  [-oma|-obs]  [-ieee] [-reduce]'
      print *,'                  [-rcred fn] [-o fn] [-ncf] [-res RES] obs_file(s)'
      print *
      print *,  'where'
      print *
      print *, '-rc fname   resource file name (default: gritas.rc)'
      print *, '-del        input files are GEOS-1 DEL files'
      print *, '             (by default input files are ODS)'
      print *, '-nopassive  passive data types will not be included'
      print *, '-oma        produces gridded O-A residuals instead'
      print *, '              of O-F residuals'
      print *, '-obs        produces gridded observed values instead'
      print *, '              of O-F residuals'
      print *, '-sigo       produces gridded observed error values instead'
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
      print *, '-oobs       produces original observations (not bias corrected)'
      print *, '-fobloc     produces forecast (background) at ob location'
      print *, '            (see note 3)'
      print *, '-ieee       output file is in IEEE binary format'
      print *, '              instead of HDF/GFIO format; both formats'
      print *, '              are GrADS compatible'
      print *, '-reduce     to apply reducer (default RC file: reducer.rc)'
      print *, '-rcred fn   specify alternative name for reducer resource file'
      print *, '-ncf        non-compliant format (apply for GSI-diag-conv files)'
      print *, '-nonorm     skips normalization of accumulated summs (good for Obs Impact)'
      print *, '-ntresh     minimum number of obs to compute the mean; if'
      print *, '             obs<ntresh then mean will set to missing'
      print *, '             (default: 0, no action)'
      print *, '-res RES    resolution definition (X=a, b, c, d/D, or e) '
      print *, '                                   a  72x46              '
      print *, '                                   b 144x91              '
      print *, '                                   c 288x181             '
      print *, '                                   D 540x361 (MERRA)     '
      print *, '                                   d 576x361             '
      print *, '                                   e 1080x721            '
      print *, '            (default: a)'
      print *, '-o bname    specifies output file base name, optional '
      print *, '            (default: gritas)'
      print *, '-ospl       Split bias, stdv and nobs into separate files'
      print *, '            (only for sdf files)'
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
      write(6,*) '     are replaced with cross-variance and the cross standard deviations, respectively. '
      print *
      write(6,*)  ' 3) The knob -fobloc can be combined with an rc-file specific parameter to calculate'
      write(6,*)  '    the ensmeble spread = f_i - <f>, where f=bkg/fcst in observation space.'
      print *
      print *, 'Revised: 15 Aug 2020 (Todling)'
      print *, 'Created:  1 September 1997 (da Silva)'
      print *
      call gritas_perror('GrISAS: not enough input parameters')
      end

!EOC




