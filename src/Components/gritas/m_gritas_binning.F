      module m_gritas_binning 
! RTodling: "stolen" from GEOS_Util - Larry Takacs coding and wrapped
! RTodling: change q10x10 to be an allocatable array
!           in module here

      implicit none
      private
      public gritas_binning

      interface gritas_binning; module procedure bin_; end interface

      contains

      subroutine bin_ ( qin,im_in,jm_in,qout,im_out,jm_out,undef,msgn )
      implicit  none
      integer   im_in ,jm_in ,msgn
      integer   im_out,jm_out
      real      undef
      real    qin(im_in ,jm_in )
      real   qout(im_out,jm_out)
      real, allocatable, dimension(:,:) :: q10x10

      allocate(q10x10(360*6,180*6))
 
c Parse Arbitray Field (im,jm) to 10'x10' Variable
c ------------------------------------------------
!     call timebeg ('___bin')
      call bin_10x10_ ( qin,im_in,jm_in,q10x10 )
 
c Bin 10'x10' Variable to Output Field (im_out,jm_out)
c ----------------------------------------------------
      call averaged_10x10_ ( q10x10,qout,im_out,jm_out,undef,msgn )
!     call timeend ('___bin')

      deallocate(q10x10)

      return
      end subroutine bin_

      subroutine averaged_10x10_ ( z10x10,z,im,jm,undef,msgn )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Average a (10m X 10m) input array to an output array (im,jm)
C
C  INPUT:
C  ======
C    z10x10 ..... Input array(360*6,180*6)
C    msgn ....... Integer Flag for scalar (0) or vector (1)
C
C  OUTPUT:
C  =======
C    z .......... Output array(im,jm)
C    im ......... Longitudinal dimension of z
C    jm ......... Latitudinal  dimension of z
C
C  NOTES:
C  ======
C    Input array z10x10  represents values within a 10min X 10min grid-box.
C             Each box is referenced by the latitude and longitude of
C             its southwest corner, not its center point.  Thus,
C             the height associated with a coordinate actually
C             represents the heights centered to the northeast of that point.
C
C    Output array z(im,jm) is assumed to be on an A-grid.
C                 z(i,j)   represents the value at the center of the grid-box.
C                 z(1,j)   is located at lon=-180.
C                 z(i,1)   is located at lat=-90.
C                 z(i,jm)  is located at lat=+90.
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      implicit none
      integer im,jm,msgn
      real  z(im,jm)
      real  dlam(im), dphi(jm)
      real  z10x10(360*6,180*6)

      integer i,j,ibeg,iend,jbeg,jend
      integer ii,jj,itmp
      real    sum1,sum2
      real    zlat,zlon
      real    lon1,lon2,wx
      real    lat1,lat2,wy
      real    lonbeg,lonend,lat,coslat
      real    latbeg,latend
      real    undef 
      real    pi,dz 
      real    lon_cmp(im)
      real    lat_cmp(jm)
!     logical defined

!     call timebeg ('____ave_10x10')
      pi   = 4.*atan(1.)
      dz   = pi/(6.*180)
      dlam = 2*pi/ im
      dphi =   pi/(jm-1)

c Compute Computational Lambda's and Phi's
c ----------------------------------------
      lon_cmp(1) = -pi
      do i=2,im
      lon_cmp(i) = lon_cmp(i-1) + dlam(i-1)
      enddo
      lat_cmp(1) = -pi*0.5
      do j=2,jm-1
      lat_cmp(j) = lat_cmp(j-1) + dphi(j-1)
      enddo
      lat_cmp(jm) =  pi*0.5


c Compute average away from poles
c -------------------------------
      do j=2,jm-1
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat-dphi(j-1)/2
      latend = zlat+dphi(j)  /2
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1.+(latbeg+pi/2)/dz
      jend = 1.+(latend+pi/2)/dz

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo
      enddo

c Compute average at South Pole
c -----------------------------
         j=1
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat
      latend = zlat+dphi(j)  /2
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1
      jend = 1.+(latend+pi/2)/dz

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo

c Compute average at North Pole
c -----------------------------
         j=jm
      do i=1,im

      zlat = lat_cmp(j)
      zlon = lon_cmp(i)

      latbeg = zlat-dphi(j-1)/2
      latend = zlat
      if( i.eq.1 ) then
      lonbeg = zlon-dlam(im) /2
      else
      lonbeg = zlon-dlam(i-1)/2
      endif
      lonend = zlon+dlam(i)  /2
      
      ibeg = 1.+(lonbeg+pi)  /dz
      iend = 1.+(lonend+pi)  /dz
      jbeg = 1.+(latbeg+pi/2)/dz
      jend = 1080

      sum1 = 0
      sum2 = 0
      do jj=jbeg,jend
      lat = -pi/2+(jj-0.5)*dz
      coslat = cos(lat)
      lat1 = -pi/2  + (jj-1)*dz
      lat2 = -pi/2  +  jj   *dz
                           wy = 1.0
      if( lat1.lt.latbeg ) wy = (lat2-latbeg)/dz
      if( lat2.gt.latend ) wy = (latend-lat1)/dz

         if(ibeg.ge.1) then
           do ii=ibeg,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         else
                 itmp = 1.+(lonbeg+0.1*dz+3*pi)/dz
           do ii=itmp,360*6
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                     wx = 1.0
           if( lon1.lt.lonbeg+2*pi ) wx = (lon2-lonbeg-2*pi)/dz
           if( lon2.gt.lonend+2*pi ) wx = (2*pi+lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
           do ii=1,iend
           if( defined(z10x10(ii,jj),undef) ) then
           lon1 = -pi  + (ii-1)*dz
           lon2 = -pi  +  ii   *dz
                                wx = 1.0
           if( lon1.lt.lonbeg ) wx = (lon2-lonbeg)/dz
           if( lon2.gt.lonend ) wx = (lonend-lon1)/dz
           sum1 = sum1 + z10x10(ii,jj)*coslat*wx*wy
           sum2 = sum2 +               coslat*wx*wy
           endif
           enddo
         endif

      enddo
      if( sum2.ne.0.0 ) then
          z(i,j) = sum1/sum2
      else
          z(i,j) = undef
      endif
      enddo

c Average Pole Values
c -------------------
      if( msgn.eq.0 ) then
      sum1 = 0
         j = 0
      do i=1,im
         if( defined(z(i,1),undef) ) then
             sum1 = sum1 + z(i,1)
                j = j + 1
         endif
      enddo
      if( j.ne.0 ) then
      z(:,1) = sum1/j
      else
      z(:,1) = undef
      endif

      sum2 = 0
         j = 0
      do i=1,im
         if( defined(z(i,jm),undef) ) then
             sum2 = sum2 + z(i,jm)
                j = j + 1
         endif
      enddo
      if( j.ne.0 ) then
      z(:,jm) = sum2/j
      else
      z(:,jm) = undef
      endif

      endif

!     call timeend ('____ave_10x10')
      return
      end subroutine averaged_10x10_

      subroutine bin_10x10_ ( z,im,jm,z10x10 )
C***********************************************************************
C
C  PURPOSE:
C  ========
C    Compute a (10m X 10m) array binned from an input array (im,jm)
C
C  INPUT:
C  ======
C    z .......... Input array(im,jm)
C    im ......... Longitudinal dimension of z
C    jm ......... Latitudinal  dimension of z
C
C  OUTPUT:
C  =======
C    z10x10 ..... Output array(360*6,180*6)
C
C  NOTES:
C  ======
C    Input array z(im,jm) is assumed to be on an A-grid.
C                z(i,j)   represents the value at the center of the grid-box.
C                z(1,j)   is located at lon=-180.
C                z(i,1)   is located at lat=-90.
C                z(i,jm)  is located at lat=+90.
C
C    Output array z10x10  represents values within a 10min X 10min grid-box.
C             Each box is referenced by the latitude and longitude of
C             its southwest corner, not its center point.  Thus,
C             the height associated with a coordinate actually
C             represents the heights centered to the northeast of that point.
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

      implicit none
      integer im,jm
      real  z(im,jm)
      real z10x10(360*6,180*6)

      integer i,j,ii,jj,ibeg,iend,jbeg,jend
      real    zlatc,zlonc
      real    lonbeg,lonend,lat
      real    latbeg,latend
      real    pi,dl,dp,dz 

!     call timebeg ('____bin_10x10')
      pi = 4.*atan(1.)
      dl = 2*pi/im
      dp = pi/(jm-1)
      dz = pi/(6.*180)
      
      do j=1,180*6
      do i=1,360*6

      zlatc = -pi/2+(j-0.5)*dz  ! Latitude  at center of 10x10 box
      zlonc = -pi  +(i-0.5)*dz  ! Longitude at center of 10x10 box

c Find bounding lat and lon on IMxJM grid
c ---------------------------------------
      iend = nint( 1.+(zlonc+pi)/dl )
      lonend = -pi + (iend-1)*dl
      if( lonend.ge.zlonc ) then
      lonbeg = -pi + (iend-2)*dl
      else
      iend = iend+1
      lonbeg = lonend
      lonend = -pi + (iend-1)*dl
      endif
      ibeg = iend-1

      jend = nint( 1.+(zlatc+pi/2)/dp )
      latend = -pi/2 + (jend-1)*dp
      if( latend.ge.zlatc ) then
      latbeg = -pi/2 + (jend-2)*dp
      else
      jend = jend+1
      latbeg = latend
      latend = -pi/2 + (jend-1)*dp
      endif
      jbeg = jend-1


      if(iend.gt.im) iend=iend-im

      if( zlonc.le.lonbeg+0.5*dl ) then
      ii = ibeg
      else
      ii = iend
      endif
      if( zlatc.le.latbeg+0.5*dp ) then
      jj = jbeg
      else
      jj = jend
      endif

      z10x10(i,j) = z(ii,jj)
      
      enddo
      enddo

!     call timeend ('____bin_10x10')
      return
      end subroutine bin_10x10_

      logical function defined ( q,undef )
      implicit none
      real     q,undef
      defined = abs(q-undef).gt.0.1*abs(undef)
      end function defined

      end module m_gritas_binning
