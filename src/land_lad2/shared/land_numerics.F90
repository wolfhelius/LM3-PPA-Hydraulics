! ============================================================================
! module numerics: a collection of useful general-purpose routines
! ============================================================================
#define __ERROR__(message) \
call my_error(mod_name,message,FATAL,__FILE__,__LINE__)
#define __NOTE__(message) \
call my_error(mod_name,message,NOTE,__FILE__,__LINE__)
#define __ASSERT__(x, message) \
if(.NOT.(x))call my_error(mod_name,message,FATAL,__FILE__,__LINE__)

module land_numerics_mod

use utilities_mod, only: logunit, error_mesg, FATAL, NOTE

implicit none
private

! ==== public interfaces =====================================================
!public :: bisect    ! finds a position of point in array of bounds
!public :: lin_int   ! linear interpolation
public :: ludcmp, lubksb ! LU decomposition and back substitution
!public :: tridiag   ! tri-diagonal system solver
public :: nearest   ! nearest point search

!public :: horiz_remap_type
!public :: horiz_remap_del
!public :: horiz_remap_print

! Matrix inversions
public :: m33inv, m44inv
! incomplete gamma
public :: gammaL, gammaU
!public :: alngam, alnorm

public :: rank_descending ! rank the input array in descending order

public :: numerics_init
! ==== end of public interfaces ==============================================


!     Linear interpolation.
!interface lin_int
!   module procedure lin_int0
!   module procedure lin_int1
!   module procedure lin_int2
!   module procedure lin_int1m
!   module procedure lin_int2m
!end interface

interface nearest
   module procedure nearest1D, nearest2D
end interface

logical :: module_is_initialized =.FALSE.
! module constants
character(len=*), parameter :: &
     mod_name = 'land_numerics_mod', &
     version  = '$Id: land_numerics.F90,v 1.1.2.3 2012/06/19 18:34:54 pjp Exp $', &
     tagname  = '$Name: no_fms_b_pjp $'

! ==== public type ===========================================================
! this data structure describes the horizontal remapping: that is, the operation 
! of copying the data from the source points to the destination points. The source
! points are not necessarily on the same PE as destination points.
type :: horiz_remap_type
   integer :: n = 0 ! number of points that need remapping on this PE
   integer, pointer :: &
       dst_i(:)=>NULL(), & ! x-indices of destination points
       dst_j(:)=>NULL()    ! y-indices of destination points
   integer, pointer :: &
       src_i(:)=>NULL(), & ! x-indices of source points
       src_j(:)=>NULL(), & ! y-indices of source points
       src_p(:)=>NULL()    ! processor number of source points
   ! data distribution map: for each processor pair that communicate 
   ! (unidirectionally), an entry in the srcPE and dstPE arrays holds their 
   ! numbers. This map is the same on each of the PEs that participate in
   ! remapping.
   integer :: mapSize = 0
   integer, pointer :: &     ! for each index:
       srcPE(:) => NULL(), & ! PE that provides the data
       dstPE(:) => NULL()    ! PE that requests and then uses the data
end type horiz_remap_type

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
! Initializes the numerics module.
subroutine numerics_init()

  module_is_initialized =.TRUE. 
  write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)

end subroutine numerics_init


! ============================================================================
!  Finds a position of point in array of bounds. Returns i, such that x is
!  between xx(i) and xx(i+1).
!  Usage:
!     value=bisect( xx, x1, periodic )
function bisect(xx, x1, periodic)
  real, intent(in)              :: xx(:)     ! array of boundaries
  real, intent(in)              :: x1        ! point to locate
  logical, intent(in), optional :: periodic  ! if present and true, the data
                                             ! domain is assumed to be periodic
  ! ---- result type ---------------------------------------------------
  integer bisect

  ! ---- local vars ----------------------------------------------------
  real    :: x              ! duplicate of input value
  integer :: low, high, mid
  integer :: n              ! size of the input array
  logical :: ascending      ! if true, the coordinates are in ascending order

  n = size(xx)
  x = x1

  ! bring the point inside bounds of the period
  if (present(periodic).AND.periodic) then
     __ASSERT__(xx(n)-xx(1)/=0,"periodic bisect: period equal to zero")
     x = modulo(x-min(xx(1),xx(n)),abs(xx(n)-xx(1)))+min(xx(1),xx(n))
  endif

  ! find the coordinates
  if (x >= xx(1).and.x<=xx(n)) then
     low = 1; high = n
     ascending = xx(n) > xx(1)
     do while (high-low > 1)
        mid = (low+high)/2
        if (ascending.eqv.xx(mid) <= x) then
           low = mid
        else
           high = mid
        endif
     enddo
     bisect = low
  else
     bisect = -1
  endif

end function bisect


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int0(data, xx, x, res)
  real, intent(in) :: data(:)    ! data to interpolate
  real, intent(in) :: xx(:)      ! coordinates of the data points
  real, intent(in) :: x          ! coordinates to interpolate to
  real, intent(inout) :: res     ! result of interpolation

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(i1)*f1+data(i2)*f2

end subroutine lin_int0


! ==============================================================================
!     Linearly interpolates 1-D data.
subroutine lin_int1(data, xx, x, res)

  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coordinates of the data points
  real, intent(in) :: x            ! coordinates to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1

  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! update the result
  res = data(:,i1)*f1+data(:,i2)*f2

end subroutine lin_int1


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2(data, tt, t, res)

  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  res = data(:,:,i1)*f1+data(:,:,i2)*f2

end subroutine lin_int2


!     Linearly interpolates 1-D data.
subroutine lin_int1m(data, xx, x, res, mask)
  real, intent(in) :: data(:,:)    ! data to interpolate
  real, intent(in) :: xx(:)        ! coordinates of the data points
  real, intent(in) :: x            ! coordinates to interpolate to
  real, intent(inout) :: res(:)    ! result of interpolation
  logical, intent(in) :: mask(:)   ! valid data mask

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(xx,x)
  i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(xx),"Coordinate is out of range")

  f2 = (x-xx(i1))/(xx(i2)-xx(i1))
  f1 = 1.0-f2

  ! finally, update the result
  where (mask) 
     res = data(:,i1)*f1+data(:,i2)*f2
  endwhere

end subroutine lin_int1m


! ==============================================================================
!     Interpolates prescribed over time.
subroutine lin_int2m(data, tt, t, res, mask)
  real, intent(in) :: data(:,:,:)  ! data to interpolate
  real, intent(in) :: tt(:)        ! time moments corresponding to data points
  real, intent(in) :: t            ! time to interpolate to
  real, intent(inout) :: res(:,:)  ! result
  logical, intent(in) :: mask(:,:) ! interpolation mask

  ! ---- local vars ----------------------------------------------------------
  integer :: i1, i2
  real    :: f1, f2

  ! find where is our time point and calculate weights
  i1 = bisect(tt,t); i2 = i1+1
  __ASSERT__(i1>0.AND.i1<size(tt),"Coordinate is out of range")

  f2 = (t-tt(i1))/(tt(i2)-tt(i1))
  f1 = 1-f2

  ! update the result
  where (mask) 
     res = data(:,:,i1)*f1+data(:,:,i2)*f2
  endwhere
end subroutine lin_int2m


! ==============================================================================
! given a matrix a(n,n) replaces it by LU decomposition of a row-wise permutation
! of itself indx(n) is an output that records the permutation. This routine is
! used in combination with lubksb to solve linear equations or invert a matrix
! example:
!    call ludcmp(a,indx)
!    call lubksb(a,indx,b1)
!    call lubksb(a,indx,b2)
subroutine ludcmp(a,indx,status)
  real,    intent(inout) :: a(:,:) ! matrix that gets replaced by its LU decomposition
  integer, intent(out)   :: indx(:) ! index of row permutations affected by partial pivoting
  integer, intent(out), optional :: status

  integer, parameter :: TINY = 1.0e-20
  integer :: n ! size of the matrix 
  integer :: i,j,k,imax
  real    :: aamax,dum,sum
  real    :: vv(size(a,1)) ! implicit scaling for each row

  n = size(a,1)
  if(present(status))status = 0

  ! find largest element in each row and calculate scaling 
  do i = 1,n
     aamax = 0.0
     do j = 1,n
        if(abs(a(i,j)) > aamax)aamax = abs(a(i,j))
     enddo
     if(.not.(aamax /= 0.0)) then 
        if(present(status))then
           status = -1; aamax = TINY
        else
           call error_mesg('ludcmp','Matrix is singular', FATAL)
        endif
     endif
     vv(i) = 1.0/aamax
  enddo

  ! loop over the columns of Crout's method
  do j=1,n
     do i = 1,j-1
        sum = a(i,j)
        do k = 1,i-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
     enddo
     aamax = 0.0 ! initialize the search for the largest pivot element
     do i = j,n
        sum = a(i,j)
        do k = 1,j-1
           sum = sum-a(i,k)*a(k,j)
        enddo
        a(i,j) = sum
        dum = vv(i)*abs(sum) ! figure of merit for the pivot
        if (dum >= aamax) then ! is it better than the best so far?
           imax = i
           aamax = dum
        endif
     enddo
     if (j /= imax) then ! do we need to interchange rows?
        ! Yes, do so
        do k=1,n
           dum=a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        enddo
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     ! if the pivot element is zero, then the matrix is singular (at least to the
     ! precision of the algorithm). For some applications on singular matrices, it 
     ! is desirable to substitute TINY for zero
     if(a(j,j)==0.0) a(j,j) = TINY

     if (j/=n)then
        ! Finally, divide by the pivot element
        dum = 1.0/a(j,j)
        do i = j+1,n
           a(i,j) = a(i,j)*dum
        enddo
     endif
  enddo ! loop over the columns
end subroutine ludcmp


! ==============================================================================
! given a LU decomposition of matrix a(n,n), permutation vector indx, and right-
! hand side b, solves the set of linear equations A*X = B 
subroutine lubksb(a,indx,b)
  real,    intent(in)    :: a(:,:)  ! LU-decomposed matrix
  integer, intent(in)    :: indx(:) ! permutation vector, as returned by the ludcmp
  real,    intent(inout) :: b(:)    ! right-hand side vector, returns solution

  integer :: i,ii,j,ll,n
  real    :: sum

  n = size(a,1)
  ii = 0
  do i = 1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii/=0) then
        do j = ii,i-1
           sum = sum - a(i,j)*b(j)
        enddo
     else if (sum /= 0.0) then
        ii = i
     endif
     b(i) = sum
  enddo
  do i = n, 1, -1
     sum = b(i)
     do j = i+1,n
        sum = sum-a(i,j)*b(j)
     enddo
     b(i) = sum/a(i,i)
  enddo
end subroutine lubksb


! ============================================================================
! given values of the tri-diagonal matrix coefficients, computes a solution
subroutine tridiag(a,b,c,r,u)
  real, intent(in)  :: a(:),b(:),c(:),r(:)
  real, intent(out) :: u(:)

  integer :: j
  real :: bet, gam(size(a))
  
  ! check that the sizes are the same
  if(size(a)/=size(b).or.size(a)/=size(c).or.size(a)/=size(r)) &
       call error_mesg('tridiag','sizes of input arrays are not equal',FATAL)
  if(size(u)<size(a)) &
       call error_mesg('tridiag','size of the result is insufficient',FATAL)
  ! check that a(1)==0 and c(N)==0
  if(a(1)/=0.or.c(size(a))/=0) &
       call error_mesg('tridiag','a(1) and c(N) must be equal to 0',FATAL)
  ! decomposition and forward substitution
  bet = b(1)
  u(1) = r(1)/bet
  do j = 2,size(a)
     gam(j) = c(j-1)/bet
     bet = b(j)-a(j)*gam(j)
     if(bet==0) &
          call error_mesg('tridiag','system is ill-defined',FATAL)
     u(j) = (r(j)-a(j)*u(j-1))/bet
  enddo
  ! backward substitution
  do j = size(a)-1,1,-1
     u(j) = u(j)-gam(j+1)*u(j+1)
  enddo
end subroutine tridiag

! ============================================================================
! finds nearest point that is not masked out in input data
! NOTE: implemented in very naive and inefficient way
subroutine nearest1D(mask, lon, lat, plon, plat, iout, jout, dist)
  logical, intent(in) :: mask(:,:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:)     ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:)     ! latitudes of input grid central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: iout, jout ! indices of nearest valid (unmasked) point
  real, optional, intent(out):: dist ! distance to the point 

  ! ---- local constants
  character(*),parameter :: mod_name='nearest1D'
  ! ---- local vars
  integer :: i,j
  real    :: r,r1

  __ASSERT__(size(mask,1)==size(lon),'sizes of "mask" and "lon" are inconsistent')
  __ASSERT__(size(mask,2)==size(lat),'sizes of "mask" and "lat" are inconsistent')
  
  r = HUGE(r)  ! some value larger than any possible distance

  do j = 1, size(mask,2)
  do i = 1, size(mask,1)
     if (.not.mask(i,j)) cycle
     r1 = distance(plon,plat,lon(i),lat(j))
     if ( r1 < r ) then
        iout = i
        jout = j
        r = r1
     endif
  enddo
  enddo
  if (present(dist)) dist = r
end subroutine nearest1D

! ============================================================================
! finds nearest point that is not masked out in input data
! this version works with 2D lon and lat fields
subroutine nearest2D(mask, lon, lat, plon, plat, iout, jout, dist)
  logical, intent(in) :: mask(:,:)  ! mask of valid input points (.true. if valid point)
  real,    intent(in) :: lon(:,:)   ! longitudes of input grid central points, radian
  real,    intent(in) :: lat(:,:)   ! latitudes of input grid central points, radian
  real,    intent(in) :: plon, plat ! coordinates of destination point, radian
  integer, intent(out):: iout, jout ! indices of nearest valid (unmasked) point
  real, optional, intent(out):: dist! distance to the point 

  ! ---- local constants
  character(*),parameter :: mod_name='nearest2D'
  ! ---- local vars 
  integer :: i,j
  real    :: r,r1

  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lon)),'shapes of "mask" and "lon" are different')
  __ASSERT__(ALL(SHAPE(mask)==SHAPE(lat)),'shapes of "mask" and "lat" are different')

  r = HUGE(r)  ! some value larger than any possible distance

  do j = 1, size(mask,2)
  do i = 1, size(mask,1)
     if (.not.mask(i,j)) cycle
     r1 = distance(plon,plat,lon(i,j),lat(i,j))
     if ( r1 < r ) then
        iout = i
        jout = j
        r = r1
     endif
  enddo
  enddo
  if (present(dist)) dist=r
end subroutine nearest2D

! ============================================================================
! private functions that calculates the distance between two points given their 
! coordinates
function distance(lon1, lat1, lon2, lat2) ; real distance
  ! calculates distance between points on unit square
  real, intent(in) :: lon1,lat1,lon2,lat2
  
  real :: x1,y1,z1, x2,y2,z2
  real :: dlon
  dlon = (lon2-lon1)
  
  z1 = sin(lat1) ;  z2 = sin(lat2)
  y1 = 0.0       ;  y2 = cos(lat2)*sin(dlon)
  x1 = cos(lat1) ;  x2 = cos(lat2)*cos(dlon)
  
  ! distance = acos(x1*x2 + z1*z2)
  distance = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
end function distance


! ============================================================================
! deallocate memory associated with the 
subroutine horiz_remap_del(map)
   type(horiz_remap_type), intent(inout) :: map
#define __DEALLOC__(x)\
if (associated(x)) then; deallocate(x); x=>NULL(); endif
   __DEALLOC__(map%dst_i)
   __DEALLOC__(map%dst_j)
   __DEALLOC__(map%src_i)
   __DEALLOC__(map%src_j)
   __DEALLOC__(map%src_p)
   map%n=0
      
   __DEALLOC__(map%srcPE)
   __DEALLOC__(map%dstPE)
   map%mapSize=0   
#undef __DEALLOC__
end subroutine

! ============================================================================
! prints remapping information
subroutine horiz_remap_print(map, prefix)
   type(horiz_remap_type), intent(in) :: map
   character(*), intent(in) :: prefix

   integer :: k

   do k = 1, map%n
      write(*,100) prefix,&
         map%src_i(k),map%src_j(k),map%src_p(k),&
         map%dst_i(k),map%dst_j(k)
   enddo
100 format(a,'(I:',i4.4,' J:',i4.4,' PE:',i4.4,') -> (I:',i4.4,' J:',i4.4,')')
end subroutine

! ======================================================================
! ranks array x in descending order: on return, idx() contains indices
! of elements of array x in descending order of x values
subroutine rank_descending(x,idx)
   real,    intent(in)  :: x(:)
   integer, intent(out) :: idx(:)

   integer :: i,n
   integer, allocatable :: t(:)
   
   n = size(x)
   do i = 1,n
      idx(i) = i
   enddo
   
   allocate(t((n+1)/2))
   call mergerank(x,idx,n,t)
   deallocate(t)
end subroutine 

! =====================================================================
! based on:
! http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
subroutine merge(x,a,na,b,nb,c,nc)
   integer, intent(in) :: na,nb,nc ! Normal usage: NA+NB = NC
   real, intent(in)       :: x(*)
   integer, intent(in)    :: a(na)    ! B overlays C(NA+1:NC)
   integer, intent(in)    :: b(nb)
   integer, intent(inout) :: c(nc)
 
   integer :: i,j,k
 
   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) >= x(b(j))) then
         c(k) = a(i) ; i = i+1
      else
         c(k) = b(j) ; j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= na)
      c(k) = a(i) ; i = i + 1 ; k = k + 1
   enddo
end subroutine merge
 
recursive subroutine mergerank(x,a,n,t)
  integer, intent(in) :: n
  real,    intent(in) :: x(*)
  integer, dimension(n), intent(inout) :: a
  integer, dimension((n+1)/2), intent (out) :: t

  integer :: na,nb
  integer :: v

  if (n < 2) return
  if (n == 2) then
     if ( x(a(1)) < x(a(2)) ) then
        v = a(1) ; a(1) = a(2) ; a(2) = v
     endif
     return
  endif      
  na=(n+1)/2
  nb=n-na

  call mergerank(x,a,na,t)
  call mergerank(x,a(na+1),nb,t)

  if (x(a(na)) < x(a(na+1))) then
     t(1:na)=a(1:na)
     call merge(x,t,na,a(na+1),nb,a,n)
  endif
end subroutine mergerank

! ==============================================================================
! Reports error, including file name and line.
subroutine my_error(mod_name, message, mode, file, line)

  character(len=*), intent(in) :: mod_name
  character(len=*), intent(in) :: message
  integer,          intent(in) :: mode
  character(len=*), intent(in), optional :: file
  integer,          intent(in), optional :: line

  ! ---- local vars ----------------------------------------------------------
  character(len=512) :: mesg
  if(present(file)) then ! assume that file and line are either both present or not
  write(mesg,'("File ",a," Line ",i4.4," :: ",a)')&
       file, line, trim(message)
  else
    mesg = trim(message)
  endif
  call error_mesg(mod_name, mesg, mode)
end subroutine

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!  http://www.davidgsimpson.com/software.html
!
!  Date:         February 6, 2009
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

END SUBROUTINE M33INV

!***********************************************************************************************************************************
!  M44INV  -  Compute the inverse of a 4x4 matrix.
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!  http://www.davidgsimpson.com/software.html
!
!  Date:         February 6, 2009
!
!  A       = input 4x4 matrix to be inverted
!  AINV    = output 4x4 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

SUBROUTINE M44INV (A, AINV, OK_FLAG)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
    DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: AINV
    LOGICAL, INTENT(OUT) :: OK_FLAG

    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
    DOUBLE PRECISION :: DET
    DOUBLE PRECISION, DIMENSION(4,4) :: COFACTOR


    DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
           A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
           A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
           A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
           A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

    IF (ABS(DET) .LE. EPS) THEN
       AINV = 0.0D0
       OK_FLAG = .FALSE.
       RETURN
    END IF

    COFACTOR(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    COFACTOR(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    COFACTOR(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    COFACTOR(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    COFACTOR(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    COFACTOR(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    COFACTOR(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    COFACTOR(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    COFACTOR(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    COFACTOR(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    COFACTOR(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    COFACTOR(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    COFACTOR(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    COFACTOR(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    COFACTOR(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    COFACTOR(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    AINV = TRANSPOSE(COFACTOR) / DET

    OK_FLAG = .TRUE.

    RETURN

END SUBROUTINE M44INV

! Gamma functions from 
! http://people.sc.fsu.edu/~jburkardt/f_src/asa239/asa239.html

function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real XVALUE, the argument of the Gamma function.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real :: alngam
  real, parameter :: alr2pi = 0.918938533204673D+00
  integer :: ifault
  real, dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real, dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real, dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real, dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real :: x, x1, x2, xvalue, y
  real, parameter :: xlge = 5.10D+05
  real, parameter :: xlgst = 1.0D+30

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end function alngam

function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real, parameter :: a1 = 5.75885480458D+00
  real, parameter :: a2 = 2.62433121679D+00
  real, parameter :: a3 = 5.92885724438D+00
  real :: alnorm
  real, parameter :: b1 = -29.8213557807D+00
  real, parameter :: b2 = 48.6959930692D+00
  real, parameter :: c1 = -0.000000038052D+00
  real, parameter :: c2 = 0.000398064794D+00
  real, parameter :: c3 = -0.151679116635D+00
  real, parameter :: c4 = 4.8385912808D+00
  real, parameter :: c5 = 0.742380924027D+00
  real, parameter :: c6 = 3.99019417011D+00
  real, parameter :: con = 1.28D+00
  real, parameter :: d1 = 1.00000615302D+00
  real, parameter :: d2 = 1.98615381364D+00
  real, parameter :: d3 = 5.29330324926D+00
  real, parameter :: d4 = -15.1508972451D+00
  real, parameter :: d5 = 30.789933034D+00
  real, parameter :: ltone = 7.0D+00
  real, parameter :: p = 0.398942280444D+00
  real, parameter :: q = 0.39990348504D+00
  real, parameter :: r = 0.398942280385D+00
  logical :: up, upper
  real, parameter :: utzero = 18.66D+00
  real :: x, y, z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 & 
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
  
end function alnorm


function gammaL ( x, p )

!*****************************************************************************80
!
!! gammaL computes the Lower Incomplete Gamma Integral
!
!  Auxiliary functions:
!
!    ALOGAM = logarithm of the gamma function, 
!    ALNORM = algorithm AS66
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by B Shea.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    B Shea,
!    Algorithm AS 239:
!    Chi-squared and Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X < 0 or P <= 0.
!
!    Output, real gammaL, the value of the lower incomplete 
!    Gamma integral.
!
  implicit none

  real :: a, an, arg, b, c, gammaL
  real :: p, pn1, pn2, pn3, pn4, pn5, pn6, rn, x
  real, parameter :: elimit = - 88.0D+00
  real, parameter :: oflo = 1.0D+37
  real, parameter :: plimit = 1000.0D+00
  real, parameter :: tol = 1.0D-14
  real, parameter :: xbig = 1.0D+08
  integer :: ifault
  logical upper

  gammaL = 0.0D+00
!
!  Check the input.
!
  if ( x < 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( p <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0

  if ( x == 0.0D+00 ) then
    gammaL = 0.0D+00
    return
  end if
!
!  If P is large, use a normal approximation.
!
  if ( plimit < p ) then

    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
    + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

    upper = .false.
    gammaL = alnorm ( pn1, upper )
    return

  end if
!
!  If X is large set gammaL = 1.
!
  if ( xbig < x ) then
    gammaL = 1.0D+00
    return
  end if
!
!  Use Pearson's series expansion.
!  (Note that P is not large enough to force overflow in ALOGAM).
!  No need to test IFAULT on exit since P > 0.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
    c = 1.0D+00
    gammaL = 1.0D+00
    a = p

    do

      a = a + 1.0D+00
      c = c * x / a
      gammaL = gammaL + c

      if ( c <= tol ) then
        exit
      end if

    end do

    arg = arg + log ( gammaL )

    if ( elimit <= arg ) then
      gammaL = exp ( arg )
    else
      gammaL = 0.0D+00
    end if
!
!  Use a continued fraction expansion.
!
  else 

    arg = p * log ( x ) - x - alngam ( p, ifault )
    a = 1.0D+00 - p
    b = a + x + 1.0D+00
    c = 0.0D+00
    pn1 = 1.0D+00
    pn2 = x
    pn3 = x + 1.0D+00
    pn4 = x * b
    gammaL = pn3 / pn4

    do

      a = a + 1.0D+00
      b = b + 2.0D+00
      c = c + 1.0D+00
      an = a * c
      pn5 = b * pn3 - an * pn1
      pn6 = b * pn4 - an * pn2

      if ( pn6 /= 0.0D+00 ) then

        rn = pn5 / pn6

        if ( abs ( gammaL - rn ) <= min ( tol, tol * rn ) ) then
          exit
        end if

        gammaL = rn

      end if

      pn1 = pn3
      pn2 = pn4
      pn3 = pn5
      pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
      if ( oflo <= abs ( pn5 ) ) then
        pn1 = pn1 / oflo
        pn2 = pn2 / oflo
        pn3 = pn3 / oflo
        pn4 = pn4 / oflo
      end if

    end do

    arg = arg + log ( gammaL )

    if ( elimit <= arg ) then
      gammaL = 1.0D+00 - exp ( arg )
    else
      gammaL = 1.0D+00
    end if

  end if

  return
end function gammaL

function gammaU ( x, p )

  real :: x, p, gammaU

  gammaU = 1.0D+00 - gammaL ( x, p )

  return
end function gammaU

end module land_numerics_mod
