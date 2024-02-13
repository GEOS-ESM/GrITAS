      module m_GFIO_OutPut
! Todling: convert to module (but kept filename the same)
      private
      public  GFIO_Output
      interface GFIO_Output
         module procedure GFIO_Output_
      end interface
      contains
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_Output --- Write del grids to IEEE GrADS file
! 
! !INTERFACE:
!
       subroutine GFIO_Output_( grid_obasen,otype,
     &                          var_list,var_descr,listsz, timeinc,
     &                          l2d, list_2d, bias_2d, stdv_2d, nobs_2d,
     &                          l3d, list_3d, bias_3d, stdv_3d, nobs_3d,
     &                          nstat, punits, yyyymmdd, hhmmss, fid,
     &                          chsq_2d, chsq_3d )
!
! !USES:
!
      use m_gritas_grids, only : im, jm, km

      implicit NONE
!
! !INPUT PARAMETERS: 
!

      character(len=*)  grid_obasen      ! output file basen name
      character(len=*)  otype            ! tag
      character(len=*)  var_descr(:)     ! tag

      integer       listsz               ! size of the list (< lmax)
      character*11  var_list(listsz)     ! variable name for output 

      integer         l2d          ! actual number of 2d grids
      integer         l3d          ! actual number of 3d grids

      integer         list_2d(l2d) ! list item for 2d grids
      integer         list_3d(l3d) ! list item for 3d grids

                                   ! Surface (2D) grids:
      real  bias_2d(im,jm,l2d)     !   time mean
      real  stdv_2d(im,jm,l2d)     !   stdv
      real  nobs_2d(im,jm,l2d)     !   number of obs per grid box
      real, optional ::  chsq_2d(im,jm,l2d,3)   !   chi-square for stddev confidence

                                   ! Upper-air (3D) grids:
      real  bias_3d(im,jm,km,l3d)  !   time mean
      real  stdv_3d(im,jm,km,l3d)  !   stdv
      real  nobs_3d(im,jm,km,l3d)  !   number of obs per grid box
      real,optional :: chsq_3d(im,jm,km,l3d,3)!   chi-square for stddev confidence
      integer nstat                ! No of stats to write to file:
                                   ! nstat = 1   ...  bias
                                   ! nstat = 2   ...  bias, stdv
                                   ! nstat = 3   ...  bias, stdv, nobs
                                   ! nstat >= 4  ...  bias, stdv, nobs, chisqr/tstud
                                   ! nstat < 0   ...  write each stat in
                                   !   separate files each with the base
                                   !   name grid_obasen // '.' // stat
                                   !   where stat = 'bias', 'stdv', or 'nobs' 
      integer        yyyymmdd      ! Year-month-day, e.g., 19971003
      integer        hhmmss        ! Hour-minute-second, e.g., 120000
  
      integer        timeinc        ! output frequency (HHMMSS)
      character(len=3) ::   punits ! units of "levels", hPa or 1
!
! !OUTPUT PARAMETERS: 
!
      integer,pointer ::    fid(:)  ! GFIO file (for closing the file later)
!
! !DESCRIPTION: Write del-grids to HDF file using the GFIO interface.
!
! !REVISION HISTORY: 
!
!  03Sep97  da Silva   Initial code.
!  17OCt97  G. P. Lou  Modifications. (RT: What modifications?)
!  26Mar98  da Silva   Added nstat parameter for GrISAS's sake.
!  10Jun04  Todling    Bug fix: fid missing from argument list.
!  05Jun09  Redder     Added code to enable routine to write each stat 
!                      in separate files (by setting n < 0).  Fixed bug
!                      by generalizing the algorithm for increment time.
!  04Aug11  Ravi       Added var_descr to the call GFIO_TGroup and 
!                      to the call GFIO_PutVar.
!  09Nov14  Todling    Add chisqr output; bug fix var_desc 
!
!EOP
!-------------------------------------------------------------------------
!BOC

      character*255 BaseN
      integer indx
      integer nymd, nhms, nstat2

      nstat2 = abs ( nstat )
      if ( nstat2 .lt. 1 .or. nstat2 .gt. 6 ) then
         print *,' nstat: ',nstat,' nstat2: ',nstat2
         call gritas_perror ( 'GFIO_Output: invalid nstat' )
      endif

      BaseN = grid_obasen
      nymd = yyyymmdd
      nhms = hhmmss

      print *, 'debug : ', nstat, fid, nymd, nhms
      indx = 1
      if ( nstat .lt. 0 ) then
         BaseN = trim ( grid_obasen ) // '.bias'
         indx=1
      end if
      call    GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, bias_2d,
     &                      l3d,  list_3d, bias_3d,
     &                      punits, nymd, nhms, fid(indx) )

      if ( nstat2 .ge. 2 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.stdv'
            indx=2
         else
            call GFIO_TIncr_ ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, stdv_2d,
     &                      l3d,  list_3d, stdv_3d,
     &                      punits, nymd, nhms, fid(indx) )
      end if

      if ( nstat2 .ge. 3 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.nobs'
            indx=3
         else
            call GFIO_TIncr_ ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, nobs_2d,
     &                      l3d,  list_3d, nobs_3d,
     &                      punits, nymd, nhms, fid(indx) )
      end if

      if ( nstat2 .ge. 4 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.chsqA'
            indx=4
         else
            call GFIO_TIncr_ ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, chsq_2d(:,:,:,1),
     &                      l3d,  list_3d, chsq_3d(:,:,:,:,1),
     &                      punits, nymd, nhms, fid(indx) )
      end if
      if ( nstat2 .ge. 5 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.chsqB'
            indx=5
         else
            call GFIO_TIncr_ ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, chsq_2d(:,:,:,2),
     &                      l3d,  list_3d, chsq_3d(:,:,:,:,2),
     &                      punits, nymd, nhms, fid(indx) )
      end if
      if ( nstat2 .ge. 6 ) then
         if ( nstat .lt. 0 ) then
            BaseN = trim ( grid_obasen ) // '.tstud'
            indx=6
         else
            call GFIO_TIncr_ ( nymd, nhms, timeinc, nymd, nhms )
         end if
         call GFIO_TGroup_( BaseN,otype,
     &                      var_list,var_descr, listsz, timeinc,
     &                      l2d,  list_2d, chsq_2d(:,:,:,3),
     &                      l3d,  list_3d, chsq_3d(:,:,:,:,3),
     &                      punits, nymd, nhms, fid(indx) )
      end if



!     All done
!     --------
      return
      end subroutine GFIO_Output_

!...............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_TGroup --- Write del grids to IEEE GrADS file for one time group
! 
! !INTERFACE:
!
       subroutine GFIO_TGroup_( grid_obasen, otype,
     &                          var_list,var_descr, listsz, timeinc,
     &                          l2d, list_2d, data_2d,
     &                          l3d, list_3d, data_3d,
     &                          punits, yyyymmdd, hhmmss, fid )
!
! !USES:
!
      use m_gritas_grids, only : im, jm, km
      use m_gritas_grids, only : glon
      use m_gritas_grids, only : glat
      use m_gritas_grids, only : AMISS
      use m_gritas_grids, only : nlevs
      use m_gritas_grids, only : plevs

      implicit NONE
      include 'gritas.h'       ! lmax
!
! !INPUT PARAMETERS: 
!

      character(len=*)  grid_obasen      ! output file basen name
      character(len=*)  otype            ! Tag
      character(len=*)  var_descr(:)     ! Description

      integer       listsz               ! size of the list (< lmax)
      character*11  var_list(listsz)     ! variable name for output 

      integer         l2d          ! actual number of 2d grids
      integer         l3d          ! actual number of 3d grids

      integer         list_2d(l2d) ! list item for 2d grids
      integer         list_3d(l3d) ! list item for 3d grids

                                   ! Surface (2D) grids:
      real  data_2d(im,jm,l2d)     !   time mean
                                   ! Upper-air (3D) grids:
      real  data_3d(im,jm,km,l3d)  !   time mean
      character(len=*) :: punits   ! level units
      integer        yyyymmdd      ! Year-month-day, e.g., 19971003
      integer        hhmmss        ! Hour-minute-second, e.g., 120000
  
      integer        timeinc       ! output frequency (HHMMSS)
!
! !OUTPUT PARAMETERS: 
!
      integer        fid           ! GFIO file (for closing the file later)
!
! !DESCRIPTION: Write del-grids to HDF file using the GFIO interface.
!
! !REVISION HISTORY: 
!
!  04Jun2009  Redder  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      integer  i, j, k, l, nymd, nhms

!     GFIO work space
!     ---------------
      character*255   fname, title, source, contact
      character*12    vunits(lmax)
      integer         kmvar(lmax)
      real            valid_range(2,lmax), packing_range(2,lmax)
      integer         rc, prec

      save fname

!     character*255 PrevBasen
!     save          PrevBasen

!     data          PrevBasen / '!@#$%^&**()_+|' /   ! just garbage

!     Create the GFIO file
!     --------------------
!     fname = trim(grid_obasen)//'.hdf'
      fname = trim(grid_obasen)//'.'//trim(otype)
         
!     New outputfile, created it
!     --------------------------
!     if ( trim(PrevBasen) .ne. trim(grid_obasen) ) then            
      if ( fid < 0 ) then

!        PrevBasen = grid_obasen

         print *
         print *, ' [] Writing to new HDF/GFIO file ', trim(fname),
     &            ' on ', yyyymmdd, hhmmss/10000, 'Z'

         title = 'Gridded O-F, O-A or Obs Values'
         source = 'From GEOS/DAS DEL or ODS files'
         contact = 'data@gmao.gsfc.nasa.gov'
         do i = 1, listsz
            vunits(i) = 'none'  ! for now
            valid_range(1,i)   = amiss      
            valid_range(2,i)   = amiss      
            packing_range(1,i) = amiss      
            packing_range(2,i) = amiss      
         end do
         kmvar = -1
         do i = 1, l2d
            l = list_2d(i)
            kmvar(l) = 0
         end do
         do i = 1, l3d
            l = list_3d(i)
            kmvar(l) = nlevs
         end do
         prec = 0               ! 32 bits

         call GFIO_Create ( trim(fname), trim(title), trim(source), trim(contact), amiss,
     &                      im, jm, nlevs, glon, glat, plevs, trim(punits), 
     &                      yyyymmdd, hhmmss, timeinc,
     &                      listsz, var_list, var_descr, vunits, kmvar,
     &                      valid_range, packing_range, prec,
     &                      fid, rc )

        print *, 'debug fid = ', fid
         if ( rc .ne. 0 ) call
     & gritas_perror ( 'GFIO_TGroup: could not create HDF/GFIO file.' )

      else

         print *
         print *, ' [] Writing to existing HDF/GFIO file ', 
     &       trim(fname), ' on ', yyyymmdd, hhmmss/10000, 'Z'

      end if


!     Write the data to file
!     ----------------------
      print *

!     Means (uses GrADS t=1 slot)
!     ---------------------------
      nymd = yyyymmdd
      nhms = hhmmss
      do i = 1, l2d
         l = list_2d(i)
         call GFIO_PutVar ( fid, var_list(l), nymd, nhms,
     &                      im, jm, 0, 1, data_2d(1,1,i),  
     &                      rc )
         if ( rc .ne. 0 ) then
            print *, 'rc = ', rc, 'var = ', var_list(l)
            call gritas_perror ( 
     &           'GFIO_TGroup: cannot write 2d bias' )
         end if
      end do
      do i = 1, l3d
         l = list_3d(i)
         do k = 1, nlevs
         call GFIO_PutVar ( fid, var_list(l), nymd, nhms,
     &                      im, jm, k, 1, data_3d(1,1,k,i),  
     &                      rc )
         if ( rc .ne. 0 ) then 
            print *, 'rc = ', rc, 'var = ', var_list(l)
            call gritas_perror ( 
     &          'GFIO_TGroup: cannot write 3d bias' )
         end if
      end do
      end do
      
!     All done
!     --------
      return
      end subroutine GFIO_TGroup_
!...............................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  GFIO_TIncr --- Increment time
! 
! !INTERFACE:
!
      subroutine GFIO_TIncr_( Date, Time, TInc, Date2, Time2 ) 
!
! !USES:
!
      implicit NONE
!
! !INPUT PARAMETERS: 
      integer  Date          ! Date (in YYYYMMDD format )
      integer  Time          ! Time (in   HHMMSS format )
      integer  TInc          ! TInc (in time increment )
!
! !OUTPUT PARAMETERS: 
      integer  Date2         ! Output date
      integer  Time2         ! Output time
!
!...............................................................

      integer JInc, JTime, ODS_Julian, ODS_CalDat
      integer JTime_New, JDay_New, Date_New, Time_New

      JInc  =       ( abs ( TInc ) / 10000 ) *  60   * 60
     .      + ( mod ( abs ( TInc ),  10000 ) / 100 ) * 60
     .      +   mod ( abs ( TInc ),    100 )
      if ( TInc .lt. 0 ) JInc = - JInc  
      JTime =             ( Time   / 10000 ) * 60    * 60 
     .      + ( mod (       Time,    10000 ) / 100 ) * 60
     .      +   mod (       Time,      100 )
      JTime_New = modulo ( JTime + JInc,  24 * 60 * 60 )
      JDay_New  = ODS_Julian ( Date ) + ( JTime + JInc - JTime_New )
     .                                / ( 24 * 60 * 60 )
      
      Date_New  = ODS_CalDat    ( JDay_New )
      Time_New  = 10000 *       ( JTime_New / ( 60 * 60 ))
     .          +   100 * ( mod ( JTime_New,    60 * 60 ) / 60 )
     .          +           mod ( JTime_New,    60 )
      Date2    = Date_New
      Time2    = Time_New 
      return
      end subroutine GFIO_TIncr_
!...............................................................

!EOC

      end module m_GFIO_OutPut

