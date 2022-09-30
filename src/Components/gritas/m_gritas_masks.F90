module m_gritas_masks
! By R. Todling, 12Apr2014 - needs prologue and comments 
use m_gritas_binning, only: gritas_binning
implicit none
private
public gritas_ini_masks
public gritas_maskout
public gritas_fnl_masks

interface gritas_ini_masks; module procedure ini_masks_; end interface
interface gritas_maskout;  module procedure maskout_; end interface
interface gritas_fnl_masks; module procedure fnl_masks_; end interface

real, allocatable, dimension(:,:) :: maskfld

integer :: im_target=-1
integer :: jm_target=-1
real    :: undefa
logical :: g5grid  ! when .t. orientation as in GEOS-5 (-180,180)
logical :: initialized_=.false.

! internal only to these module
real :: lonMin_ = -180.0  ! caution: not to be confused w/ that in gritas_grids
real :: latMin_ = -90.0   ! caution: not to be confused w/ that in gritas_grids
contains

subroutine ini_masks_ (fname,im_user,jm_user,usr_mask)
  character(len=*), intent(in) :: fname ! file name where land-mask found
  integer,intent(in) :: im_user, jm_user
  character(len=*),intent(in) :: usr_mask

  integer id,rc
  integer im,jm,lm,ntime,nvars,ngatts
  integer nymd, nhms, inc_sec
  logical binning

  real, allocatable, dimension(:,:) :: fld2d
  character(len=80) maskname

  im_target = im_user
  jm_target = jm_user

  call gfio_open       ( trim(fname),1,id,rc )
  call gfio_diminquire ( id,im,jm,lm,ntime,nvars,ngatts,rc)
  call GetBegDateTime  ( id, nymd, nhms, inc_sec, rc )
  call grid_orientation_(id,im,jm,lm,ntime,nvars)
  if(.not.g5grid) then 
     call gritas_perror ('ini_landmask_: cannot handle grid orientation')
  endif

  allocate(fld2d(im,jm))

! get mask fraction
  maskname='fr'//trim(usr_mask)
  call gfio_getvar ( id,trim(maskname),nymd,nhms,im,jm,0,1,fld2d,rc )
  if(rc/=0) then
     call gritas_perror ('ini_landmask_: trouble reading mask file')
  endif
  call gfio_close ( id,rc )

! when target grid different from input grid, bin input to target
  binning = im/=im_target .or. jm/=jm_target
  allocate(maskfld(im_target,jm_target))
  if (binning) then
     call gritas_binning ( fld2d,im,jm,maskfld,im_target,jm_target,undefa,0 )
  else
     maskfld = fld2d 
  endif

  deallocate(fld2d)

  initialized_=.true.

end subroutine ini_masks_

! this is done for safety only since binning routines 
! assumes grid to be oriented from -180 to 180
subroutine grid_orientation_(id,im,jm,lm,ntime,nvars)
  ! all this mambo-jambo to find out orientation of grid
  integer,intent(in)  :: id,im,jm,lm,ntime,nvars

  real,    allocatable ::    lat(:)
  real,    allocatable ::    lon(:)
  real,    allocatable ::    lev(:)
  real,    allocatable :: vrange(:,:)
  real,    allocatable :: prange(:,:)
  integer, allocatable :: yymmdd(:)
  integer, allocatable :: hhmmss(:)
  integer, allocatable ::  kmvar(:)
  integer, allocatable ::   nloc(:)
  integer timinc,rc

  character*256  title
  character*256  source
  character*256  contact
  character*256  levunits
  character*256, allocatable ::  vname(:)
  character*256, allocatable :: vtitle(:)
  character*256, allocatable :: vunits(:)

  allocate ( lon(im) )
  allocate ( lat(jm) )
  allocate ( lev(lm) )
  allocate ( yymmdd(ntime) )
  allocate ( hhmmss(ntime) )
  allocate (  vname(nvars) )
  allocate ( vtitle(nvars) )
  allocate ( vunits(nvars) )
  allocate (  kmvar(nvars) )
  allocate ( vrange(2,nvars) )
  allocate ( prange(2,nvars) )

  call gfio_inquire ( id,im,jm,lm,ntime,nvars, &
                      title,source,contact,undefa, &
                      lon,lat,lev,levunits, &
                      yymmdd,hhmmss,timinc, &
                      vname,vtitle,vunits,kmvar, &
                      vrange,prange,rc )
  if( lon(1) < -1.0 ) then
    g5grid=.true.
  endif

  deallocate ( lon )
  deallocate ( lat )
  deallocate ( lev )
  deallocate ( yymmdd )
  deallocate ( hhmmss )
  deallocate (  vname )
  deallocate ( vtitle )
  deallocate ( vunits )
  deallocate (  kmvar )
  deallocate ( vrange )
  deallocate ( prange )

end subroutine grid_orientation_

subroutine maskout_ (lat,lon,qcx,which)
  real, intent(in)    :: lon(:)
  real, intent(in)    :: lat(:)
  integer, intent(inout) :: qcx(:)
  character(len=*),intent(in) :: which ! this could eventually be land/ocean/ice

  integer i,j,n,im,jm,nobs,nexcl
  real dlon, dlat

  im=im_target
  jm=jm_target
  nexcl=0

! these are not to be confused w/ the grid defined in gritas_grids ...
! these are internal to this module only
  dLon = 360. / im
  dLat = 180. / ( jm - 1 )

  if(.not.initialized_) return ! nothing to do

  if(trim(which)/='ocean'.and.trim(which)/='land'.and.trim(which)/='seaice'.and.&
     trim(which)/='landice'.and.trim(which)/='lake') then
     return ! unknown mask, do nothing
  endif

  nobs=size(lat)

  do n=1,nobs

!    Horizontal gridbox: gridpoint is assumed at center of box
!    ---------------------------------------------------------
     i = 1 + nint ( (lon(n) - LonMin_) / dLon )
     if ( i .eq. (im+1) ) i = 1                      ! periodic bc
     if ( i .eq.   0    ) i = im                     ! periodic bc
     j = 1 + nint ( (lat(n) - LatMin_) / dLat )

     i = min(im,max(1,i))
     j = min(jm,max(1,j))

!    Location with fraction of ocean larger than 0.9 are taken for ocean
!    -------------------------------------------------------------------
     if(maskfld(i,j)<0.99.and.qcx(n)==0) then
        nexcl=nexcl+1
        qcx(n) = 1
     endif

  enddo
  print*, 'Number of Obs Exluded due to maskout: ', nexcl 
  
end subroutine maskout_

subroutine fnl_masks_
  if(.not. initialized_) return
  deallocate(maskfld)
  initialized_=.false.
end subroutine fnl_masks_

end module m_gritas_masks

