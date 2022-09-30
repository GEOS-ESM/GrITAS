subroutine generate_raob_speed (ods)
use m_MergeSorts
use m_odsmeta, only : X_RED 
use m_ods
implicit none
type(ods_vect)  ods
integer ii,iu,iv,ndim,ired,nx,ny
integer i,j,k,jt,jj
real d
real,allocatable :: ulon(:),ulat(:),ulev(:)
integer,allocatable::ukt(:),ukx(:),utime(:)
real,allocatable :: vlon(:),vlat(:),vlev(:)
integer,allocatable::vkt(:),vkx(:),vtime(:)
real,allocatable :: obs(:),omf(:),oma(:)
integer, allocatable :: is(:), js(:)
real :: ana,ubkg,vbkg
real :: bkg,uana,vana
! count number of u/v components
nx=0;ny=0
do ii=1,ods%data%nobs
   if (ods%data%qcexcl(ii)/=0)  cycle
   if (ods%data%kx(ii)/=220) cycle
   if (ods%data%kt(ii)==4) nx=nx+1
   if (ods%data%kt(ii)==5) ny=ny+1
enddo
if(nx/=ny) then
  print*, 'Unequal numnber of u/v data, likely troubled ', nx,ny
endif
allocate(ulon(nx),ulat(nx),ulev(nx),utime(nx),ukx(nx),ukt(nx))
allocate(vlon(ny),vlat(ny),vlev(ny),vtime(ny),vkx(ny),vkt(ny))
allocate(omf(nx),oma(nx),obs(nx))
iu=0;iv=0
do ii=1,ods%data%nobs
   if (ods%data%qcexcl(ii)/=0)  cycle
   if (ods%data%kx(ii)/=220) cycle
   if (ods%data%kt(ii)==4) then 
      iu=iu+1
      ukt (iu) =ods%data%kt (ii)
      ukx (iu) =ods%data%kx (ii)
      ulon(iu) =ods%data%lon(ii)
      ulat(iu) =ods%data%lat(ii)
      ulev(iu) =ods%data%lev(ii)
      utime(iu)=ods%data%time(ii)
   endif
   if (ods%data%kt(ii)==5) then 
      iv=iv+1
      vkt (iv) =ods%data%kt (ii)
      vkx (iv) =ods%data%kx (ii)
      vlon(iv) =ods%data%lon(ii)
      vlat(iv) =ods%data%lat(ii)
      vlev(iv) =ods%data%lev(ii)
      vtime(iv)=ods%data%time(ii)
   endif
enddo
allocate ( is(nx) )
call IndexSet  ( nx, is )
call IndexSort ( nx, is, ukx  (1:nx), descend=.false. )
call IndexSort ( nx, is, ukt  (1:nx), descend=.false. )
call IndexSort ( nx, is, ulat (1:nx), descend=.false. )
call IndexSort ( nx, is, ulon (1:nx), descend=.false. )
call IndexSort ( nx, is, ulev (1:nx), descend=.false. )
call IndexSort ( nx, is, utime(1:nx), descend=.false. )
allocate ( js(ny) )
call IndexSet  ( ny, js )
call IndexSort ( ny, js, vkx  (1:ny), descend=.false. )
call IndexSort ( ny, js, vkt  (1:ny), descend=.false. )
call IndexSort ( ny, js, vlat (1:ny), descend=.false. )
call IndexSort ( ny, js, vlon (1:ny), descend=.false. )
call IndexSort ( ny, js, vlev (1:ny), descend=.false. )
call IndexSort ( ny, js, vtime(1:ny), descend=.false. )

! For matched cases ...
  i = 1
  j = 1
  k = 0
  ired=0
  do while (i .le. nx) ! for each u observations:
     ii = is(i)
     do jt = j, ny     ! look for a v match
        jj = js(jt)
        do
           d = ukt  (ii)-4        ; if (d /= 0) exit
           d = vkt  (jj)-5        ; if (d /= 0) exit
           d = ukx  (ii)-220      ; if (d /= 0) exit
           d = vkx  (jj)-220      ; if (d /= 0) exit
           d = ukx  (ii)-vkx  (jj); if (d /= 0) exit
           d = utime(ii)-vtime(jj); if (d /= 0) exit
           d = ulev (ii)-vlev (jj); if (d /= 0) exit
           d = ulon (ii)-vlon (jj); if (d /= 0) exit
           d = ulat (ii)-vlat (jj);             exit
        end do
        if (d==0) then   ! all match
            ! obs
            obs(ii) = sqrt(ods%data%obs(ii)*ods%data%obs(ii) + &
                           ods%data%obs(jj)*ods%data%obs(jj))
            ! omf
            ubkg = ods%data%obs(ii)-ods%data%omf(ii) 
            vbkg = ods%data%obs(jj)-ods%data%omf(jj) 
            bkg  = sqrt(ubkg*ubkg + vbkg*vbkg)
            omf(ii) = obs(ii)-bkg
            ! oma
            uana = ods%data%obs(ii)-ods%data%oma(ii) 
            vana = ods%data%obs(jj)-ods%data%oma(jj) 
            ana  = sqrt(uana*uana + vana*vana)
            oma(ii) = obs(ii)-ana

            j = min(jt + 1,ny)
            i = i + 1
            k = k + 1
            exit         ! stop looking for v match
        elseif (d<0 .OR. jt==ny) then  ! no match exists
            j = jt
            i = i + 1
            ods%data%qcexcl(ii)=X_RED
            exit         ! stop looking
        end if
     end do              ! loop on v wind
  end do                 ! loop on u wind
if(ired/=0) print *, 'Red-flagged this many obs: ',ired
!
  do while (i .le. nx) ! for each u observations:
     ii = is(i)
     if(ods%data%qcexcl(ii)/=0) cycle
     if(ods%data%kt(ii)/=4)     cycle
     if(ods%data%kx(ii)/=220)   cycle
     ods%data%obs(ii) = obs(ii)
     ods%data%omf(ii) = omf(ii)
     ods%data%oma(ii) = oma(ii)
  end do
!
deallocate(is,js)
deallocate(omf,oma,obs)
deallocate(vlon,vlat,vlev,vtime,vkx,vkt)
deallocate(ulon,ulat,ulev,utime,ukx,ukt)
end subroutine generate_raob_speed
