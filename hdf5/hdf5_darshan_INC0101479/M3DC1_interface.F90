!> \file
!> \brief load RMP file from M3DC1 code
!> \author
!> Lei Shi
!> \date
!> Oct 11, 2016 -- Mar 23, 2017

!######################################################################
!> \brief Module for reading plasma configuration file prepared by M3D-C1 code,
!> specifically with magnetic islands created by RMP.
!
!> The module contains:
!
!> constants for
!> the physical unit conversion coefficients from m3dc1 output file
!> quantities to CGS quantities
!
!> Public Subroutines: 
!
!> read_m3dc1_dims, read_m3dc1_profs, check_Jacobian,
!> load_m3dc1_alpha,
!######################################################################
module m3dc1

use precision, only: lk
use global_parameters, only: r0, b0, pi, fielddir, gtcout, mtoroidal
use equilibrium, only: lst, lsp, lszeta, ndim, psiw, ped, spdpsi, spdtheta, &
  spdzeta, qpsi, gpsi, rpsi, bsp, xsp, zsp, gsp, jsp, rd, cpsi
use spline_function, only: construct_spline1d, construct_spline2d, &
  construct_spline3d
use smooth_function, only: smooth_1d

implicit none
private
public read_m3dc1_dims, read_m3dc1_eq, check_Jacobian, load_m3dc1_alpha

! The physical unit conversion coefficients from m3dc1 output file
! quantities to CGS quantities
  real(lk), parameter :: m3dc1_length_unit=100.0_lk, m3dc1_time_unit=1.0_lk, &
       m3dc1_magnetic_unit=10000.0_lk

contains

!#####################################################################

! Subroutine read_m3dc1_dims
!
!> \brief read in the dimensions of psi and theta grids
!
!> Read in dimention of psi and theta grid in M3DC1 file
!> and set the lsp and lst parameters.
!
!> @param[in] m3dfile
!> A string contains the path to the M3DC1 file
!
subroutine read_m3dc1_dims(m3dfile)

  include 'netcdf.inc'

  ! netcdf file name should be passing in
  character(len=*), intent(in):: m3dfile

  ! some preset dimension names
  character(len=*), parameter :: psi_dimname = 'npsi', theta_dimname = 'mpol'

  ! netcdf ids
  integer :: ncid, psi_dimid, theta_dimid


  ! read in lsp and lst from the netcdf file

  ! open the netcdf file
  call check_nf( nf_open(m3dfile, nf_nowrite, ncid))

  ! get the dimension ids
  call check_nf( nf_inq_dimid(ncid, psi_dimname, psi_dimid))
  call check_nf( nf_inq_dimid(ncid, theta_dimname, theta_dimid))

  ! get the dimension lengths
  call check_nf( nf_inq_dimlen(ncid, psi_dimid, lsp))
  call check_nf( nf_inq_dimlen(ncid, theta_dimid, lst))

  ! close the netcdf file
  call check_nf( nf_close(ncid))

  ! assign the global variables lsp and lst
  ! Note that in m3dc1 output files, inner most psi is at the first stepsize
  ! delta_psi, and the out most psi is the one right inside the separatrix, psi_ed-delta_psi.
  ! So, to include the magnetic axis psi=0 and the Last Closed Flux Surface (LCFS)
  ! psi=psi_ed, we need to add two value points, one at axis, one at LCFS.
  lsp = lsp+2

  ! For theta, m3dc1 output file range from 0 to 2pi-delta_theta, to include the periodic
  ! data point at 2pi, we add one data point.
  lst = lst+1

  ! write out the read in dimensions
  write(gtcout, *), "*** Read in M3DC1 dimensions: lsp = ", lsp, ", lst = ", lst
end subroutine read_m3dc1_dims
!##############################################################################################


! Subroutine read_m3dc1_eq
!
! Description:
!> read in the axisymmetric equilibrium quantities
!>
!> @param[in] m3dfile
!> A string contains the path to the M3DC1 file
!>
!> @note
!> The subroutine assigns values to the following global variables:
!>  - Spline constants:
!>     -# `spdpsi`
!>  - New radial coordinate spline coefficient
!>     -# `rpsi` (outer mid-plane R-R0 as function of psi_p)
!>  - 2D spline coefficients:
!>     -# Total B field : `bsp`,
!>     -# R, Z coordinates : `xsp`, `zsp`,
!>     -# Jacobian : `gsp`,
!>     -# Toroidal current: `cpsi` (1D), `rd` (essentially 1D, coeffs in dtheta are all 0), duplicated with cpsi
!>  - 1D spline coefficients:
!>     -# Safety factor : `qpsi`,
!>     -# Poloidal current : `gpsi`,
!>     -# Minor radius : `rpsi`,
!>     -# Toroidal current: `cpsi`

subroutine read_m3dc1_eq(m3dfile)
  include 'netcdf.inc'

  ! netcdf file name should be passing in
  character(len=*), intent(in):: m3dfile

  ! netcdf variable names
  character(len=*), parameter :: q_name='q', g_name='F', I_name='current', &
       psi_name='psi', psin_name='psi_norm', r_name='rpath', z_name='zpath', &
       bp_name='Bp'
  ! netcdf ids
  integer :: ncid, q_id, g_id, I_id, psi_id, psin_id, r_id, z_id, bp_id

  ! Note that M3DC1 data is stored as single precision arrays
  ! temporary container for 1D variables
  real(4), dimension(:), allocatable :: temp1d
  ! temporary container for 2D variables
  real(4), dimension(:,:), allocatable :: temp2d

  ! loop index
  integer :: i, j
  ! temporary psi information storage
  ! we need to use these values to calculate the correct psiw since poloidal flux in M3DC1 may have different definition on axis and on wall
  real(lk) :: psi_0, psi_1, psin_0, psin_1

  ! The vacuum permeability constant in SI
  real(lk), parameter:: mu0 = 1.2566370614D-6

  ! open m3dc1 netcdf file
  call check_nf( nf_open(m3dfile, nf_nowrite, ncid))

  ! get the variable dimensions
  call check_nf( nf_inq_varid(ncid, q_name, q_id))
  call check_nf( nf_inq_varid(ncid, g_name, g_id))
  call check_nf( nf_inq_varid(ncid, I_name, I_id))
  call check_nf( nf_inq_varid(ncid, psi_name, psi_id))
  call check_nf( nf_inq_varid(ncid, psin_name, psin_id))
  call check_nf( nf_inq_varid(ncid, r_name, r_id))
  call check_nf( nf_inq_varid(ncid, z_name, z_id))
  call check_nf( nf_inq_varid(ncid, bp_name, bp_id))

  ! read in the variables

  ! 1D variables
  allocate (temp1d(lsp-2))

  ! Poloidal flux in physical unit (weber * meter**2): psi0, psi1
  call check_nf( nf_get_var(ncid, psi_id, temp1d))
  ! store the first and last values of psi
  psi_0 = temp1d(1)
  psi_1 = temp1d(lsp-2)

  ! Normalized poloidal flux in range [0, 1]: psin0, psin1
  call check_nf( nf_get_var(ncid, psin_id, temp1d))
  ! store the first and last values of psi
  psin_0 = temp1d(1)
  psin_1 = temp1d(lsp-2)

  ! The last closed flux surface ped
  ! Calculate ped in terms of psi_0, psi_1, psin_0, psin_1
  ! (psi_0-psi_1)/ped = (psin_0 - psin_1)
  ped = (psi_0-psi_1)/(psin_0-psin_1)
  ! convert unit
  ped = ped * m3dc1_length_unit**2 * m3dc1_magnetic_unit
  ! set psiw to be the same as ped
  psiw = ped

  ! Safety factor: qpsi
  call check_nf( nf_get_var(ncid, q_id, temp1d))

  ! fill in the values of spline coefficients
  qpsi(1, 2:lsp-1)=temp1d(:)

  ! on-axis and at psiw values need to be extrapolated
  ! both sides use linear extrapolation with derivative calculated by '431' formula
  call linear_extrapolation(qpsi(1,:), 1, -1)
  call linear_extrapolation(qpsi(1,:), lsp, 1)

  ! Covariant toroidal magnetic field: gpsi
  call check_nf( nf_get_var(ncid, g_id, temp1d))

  ! fill in the values of spline coefficients
  ! on-axis value set to be the same as the first point
  gpsi(1, 2:lsp-1)=temp1d(:)
  ! extrapolated on both sides
  call linear_extrapolation(gpsi(1,:), 1, -1)
  call linear_extrapolation(gpsi(1,:), lsp, 1)
  ! convert unit
  gpsi(1, :) = gpsi(1, :) * m3dc1_magnetic_unit * m3dc1_length_unit

  ! Covariant poloidal magnetic field: cpsi
  call check_nf( nf_get_var(ncid, I_id, temp1d))
  ! fill in the values of spline coefficients
  ! on-axis I value should be zero since no toroidal current enclosed in the "zero width" flux surface
  cpsi(1, 1) = 0
  cpsi(1, 2:lsp-1)=temp1d(:)
  ! extrapolate on outside
  call linear_extrapolation(cpsi(1,:), lsp, 1)

  ! The read in quantity is toroidal current, times 1/2pi and magnetic permeability mu0 to obtain the
  ! correct poloidal magnetic field quantity
  cpsi(1, :) = cpsi(1, :)*mu0/(2*pi)
  ! convert the unit
  cpsi(1, :) = cpsi(1, :)* m3dc1_length_unit * m3dc1_magnetic_unit

  ! 2D variables
  ! Note that in M3DC1 output file, 2D quantities are stored so that the fastest changing index is in theta
  ! So, in Fortran, the first dimension should be theta, the second dimension is psi
  allocate (temp2d(lst-1, lsp-2))

  ! Major Radius : R
  call check_nf( nf_get_var(ncid, r_id, temp2d))
  ! Be careful that in GTC spline coefficients, psi is in front of theta

  do i=2, lsp-1
     do j=1, lst-1
        xsp(1, i, j) = temp2d(j, i-1)
     enddo
     ! periodic in theta
     xsp(1, i, lst) = temp2d(1, i-1)
  enddo
  ! on-axis value take the average of the inner most flux surface
  xsp(1, 1, 1)=0
  do j=1, lst-1
     xsp(1, 1, 1) = xsp(1,1,1)+ temp2d(j, 1)
  enddo
  xsp(1, 1, 1) = xsp(1, 1, 1)/(lst-1)
  xsp(1, 1, :) = xsp(1, 1, 1)

  ! Extrapolate to obtain the values at LCFS.
  ! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
  do j=1, lst
     call linear_extrapolation(xsp(1,:,j), lsp, 1)
  enddo
  ! convert from M3DC1 unit to CGS unit
  xsp(1,:,:) = xsp(1,:,:) * m3dc1_length_unit
  ! reset the magnetic axis R0 value
  r0 = xsp(1, 1, 1)
  ! reset the magnetic axis B0 value
  b0 = abs(gpsi(1, 1)/r0)

  ! Vertical Coordinate: Z
  call check_nf( nf_get_var(ncid, z_id, temp2d))
  ! Be careful that in GTC spline coefficients, psi is in front of theta
  do i=2, lsp-1
     do j=1, lst-1
        zsp(1, i, j) = temp2d(j, i-1)
     enddo
     ! periodic in theta
     zsp(1, i, lst) = temp2d(1, i-1)
  enddo
  ! Since theta=0 line is defined as horizontal, on-axis Z value should equal Z values at theta=0
  zsp(1, 1, :)=zsp(1, 2, 1)

  ! Z value at separatrix needs to be extrapolated from inside.
  ! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
  do j=1, lst
     call linear_extrapolation(zsp(1,:,j), lsp, 1)
  enddo
  ! convert the unit
     zsp(1,:,:) = zsp(1,:,:) * m3dc1_length_unit

  ! Poloidal magnetic field strength: Bp
  ! We'll temporarily store this quantity in bsp
  call check_nf( nf_get_var(ncid, bp_id, temp2d))
  ! Be careful that in GTC spline coefficients, psi is in front of theta
  do i=2, lsp-1
     do j=1, lst-1
        bsp(1, i, j) = temp2d(j, i-1)
     enddo
     ! periodic in theta
     bsp(1, i, lst) = temp2d(1, i-1)
  enddo
  ! on-axis Bp should be zero
  bsp(1, 1, :)=0

  ! Bp value at psiw needs to be extrapolated from inside.
  ! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
  do j=1, lst
     call linear_extrapolation(bsp(1,:,j), lsp, 1)
  enddo

  ! convert the unit
  bsp(1, :, :) = bsp(1, :, :)* m3dc1_magnetic_unit

  ! Calculate the total B based on B_p, gpsi and xsp
  do i=1, lsp
     do j=1, lst
        bsp(1,i,j) = sqrt(bsp(1,i,j)**2 + (gpsi(1,i)/xsp(1,i,j))**2)
     enddo
  enddo

  ! close netcdf file
  call check_nf (nf_close(ncid))

  ! convert the coordinates
  call coordinate_conversion

  ! Normalize the quantities into GTC unit
  ! GTC units:
  ! Length : R0
  ! Magnetic field : B0
  ! Time : 1/omega_cp
  call normalize_profs

  ! Additional smooth the periodic direction for noise control
  ! Uncomment the following section for additional smoothing
  !do i=1, lsp
  !   call smooth_1d(xsp(1,i,:), lst, boundary_condition=1, smooth_type=0)
  !   call smooth_1d(zsp(1,i,:), lst, boundary_condition=1, smooth_type=0)
  !   call smooth_1d(bsp(1,i,:), lst, boundary_condition=1, smooth_type=0)
  !enddo

  ! generate the spline coefficients
  call create_splines

  ! populate rd array using the calculated cpsi
  ! cpsi is actually used in GTC calculation, rd is just for equilibrium output
  rd(1:3,:,1) = cpsi(:,:)
  do i=1,3
     do j=1,lsp
        rd(i,j,:) = rd(i,j,1)
     enddo
  enddo
  rd(4:9,:,:)=0

  ! check the self-consistency of the equilibrium, write out warning when error is too large
  ! Jacobian splines will be created here.
  ! Jacobians are checked at the end of eqdata subroutine. No need to check here.
  ! call check_m3dc1

  ! write out a message
  write(gtcout, *) '***M3DC1 equilibrium profiles successfully read.'
  write(gtcout, *) 'q_axis = ', qpsi(1,2),'q_LCFS = ', qpsi(1,lsp-1)
  write(gtcout, *) 'R0=', r0, 'R_min/R0=', xsp(1, lsp, lst/2), 'R_max/R0=', xsp(1, lsp, 1)
  write(gtcout, *) 'Z0/R0=', zsp(1,1,1), 'Z_min(max)/R0=', minval(zsp(1,:,:)), 'Z_max(min)/R0=', maxval(zsp(1,:,:))
  write(gtcout, *) 'B0=', b0, 'B_out/B0=', bsp(1, lsp, 1), 'B_in/B0=', bsp(1, lsp, lst/2)
  if (fielddir == 0) then
     write(gtcout, *) 'fielddir=', fielddir, 'I:out, B:out'
  else if (fielddir ==1) then
     write(gtcout, *) 'fielddir=', fielddir, 'I:out, B:in'
  else if (fielddir ==2) then
     write(gtcout, *) 'fielddir=', fielddir, 'I:in, B:in'
  else
     write(gtcout, *) 'fielddir=', fielddir, 'I:in, B:out'
  endif
  write(gtcout, *) '**********************************************'
end subroutine read_m3dc1_eq
!####################################################################################

! Subroutine normalize_profs
!> Purpose: Normalize the equilibrium quantities into GTC units
!>
!> GTC units:
!> - Length : R0 (major radius of magnetic axis)
!> - Magnetic field : B0 (on-axis magnetic field)
!> - Time : 1/omega_cp (invert cyclotron frequency of proton)
subroutine normalize_profs
  ! some useful derived units
  real(lk) :: flux_unit

  ! normalize lengths
  xsp(:,:,:) = xsp(:,:,:)/r0
  zsp(:,:,:) = zsp(:,:,:)/r0

  ! normalize magnetic field
  bsp(:,:,:) = bsp(:,:,:)/b0

  ! normalize covariant field quantities
  cpsi(1,:) = cpsi(1,:)/(r0*b0)
  gpsi(1,:) = gpsi(1,:)/(r0*b0)

  ! normalize flux
  flux_unit = r0*r0*b0
  psiw = psiw/flux_unit
  ped = ped/flux_unit

end subroutine normalize_profs
!#####################################################################################

! Subroutine coordinate_conversion
!> Purpose: convert all the quantities in M3DC1 coordinates into GTC coordinates
!>
!> See the documentation (Doc/DeveloperNotes/M3DC1_loader_notes.pdf) for detailed
!> convention information

subroutine coordinate_conversion

  integer:: I_sign, B_sign, j
  real(lk):: temp(lsp)

  ! The signs of two quantities are required to determine the Ip and B_T directions
  ! For simplicity, we just use the signs of toroidal current and the g(psi)
  if (cpsi(1, lsp)>0) then
     I_sign = 1
  else
     I_sign = -1
  endif

  if (gpsi(1, 1)>0) then
     B_sign = 1
  else
     B_sign = -1
  endif
  ! Check if sign of qpsi is compatible with I_sign and B_sign
  ! sign of q should always equal I_sign*B_sign
  if (I_sign+B_sign == 0 .and. qpsi(1,1)>0) then
     write(gtcout, *) 'M3DC1 EQUILIBRIUM ERROR: Signs of toroidal current and toroidal B &
          field don''t agree with sign of q. Check the equilibrium file.'
     write(gtcout, *) 'cpsi(1,lsp) = ', cpsi(1,lsp), ', I_sign = ', I_sign, ', gpsi(1,1)= ', gpsi(1,1), &
          ', B_sign = ', B_sign
     stop
  elseif (I_sign+B_sign /= 0 .and. qpsi(1,1)<0) then
     write(gtcout, *) 'M3DC1 EQUILIBRIUM ERROR: Signs of toroidal current and toroidal B &
               field don''t agree with sign of q. Check the equilibrium file.'
      write(gtcout, *) 'cpsi(1,lsp) = ', cpsi(1,lsp), ', I_sign = ', I_sign, ', gpsi(1,1)= ', gpsi(1,1), &
          ', B_sign = ', B_sign
     stop
  endif

  gpsi(1,:) = I_sign * gpsi(1,:)
  cpsi(1,:) = abs(cpsi(1,:))
  psiw = abs(psiw)
  ped = abs(ped)

  ! if I_sign is negative, then M3DC1 theta and GTC theta are opposite
  ! All 2D quantities need to be flipped in theta,
  ! i.e. f_gtc(theta_gtc) = f_m3d(2*pi - theta_gtc)
  ! Need a buffer for switching the two
  if (I_sign<0) then
     do j=1,lst/2  ! Switch the first half with the second half
        temp(:) = xsp(1,:,j)
        xsp(1,:,j) = xsp(1,:,lst-j+1)
        xsp(1,:,lst-j+1) = temp(:)

        temp(:) = zsp(1,:,j)
        zsp(1,:,j) = zsp(1,:,lst-j+1)
        zsp(1,:,lst-j+1) = temp(:)

        temp(:) = bsp(1,:,j)
        bsp(1,:,j) = bsp(1,:,lst-j+1)
        bsp(1,:,lst-j+1) = temp(:)
     enddo
  endif

  ! set fielddir value to corresponding cases, for constructing field line
  ! aligned GTC mesh
  ! 4 cases: fielddir =
  ! 0 : I out, B out
  ! 1 : I out, B in
  ! 2 : I in, B in
  ! 3 : I in, B out
  if (I_sign < 0) then
     if (B_sign < 0) then
        fielddir = 0
     else
        fielddir = 1
     endif
  else
     if (B_sign < 0) then
        fielddir = 3
     else
        fielddir = 2
     endif
  endif

end subroutine coordinate_conversion
!############################################################################

! Subroutine create_splines
!> Purpose: generate the spline coefficients for all the read in quantities

subroutine create_splines
  integer:: i,j

  ! global parameter spdpsi will be used in the future
  spdpsi = psiw/(lsp-1)

  ! Create the new minor radius coordinate rpsi
  ! rpsi is simply the radial distance along the outer mid-plane
  ! note that xsp is already normalized to R0, so rpsi is automatically normalized
  do i=1, lsp
     rpsi(1, i) = xsp(1,i,1)-1
  enddo

  call construct_spline1d(lsp, spdpsi, rpsi,1)
  call construct_spline1d(lsp,spdpsi,qpsi,0)
  call construct_spline1d(lsp,spdpsi,gpsi,0)
  call construct_spline1d(lsp,spdpsi,cpsi,0)

  ! global paramter spdtheta will be used in the future
  spdtheta = 2*pi/(lst-1)

  call construct_spline2d(lsp, lst, spdpsi, spdtheta, xsp, 1, 2)
  call construct_spline2d(lsp, lst, spdpsi, spdtheta, zsp, 1, 2)
  call construct_spline2d(lsp, lst, spdpsi, spdtheta, bsp, 1, 2)

end subroutine create_splines
!##########################################################################################

! Subroutine check_m3dc1
!> Purpose: check the consistence of the read in quantities
!>
!> Method:
!>
!> 1. Check Jacobian
!>    \f[
!>      1/\mathcal{J} = \nabla \psi \times \nabla \theta \cdot \nabla
!>      \zeta
!>    \f]
!>    \f[
!>      1/\mathcal{J} = \frac{\vec{B} \cdot \vec{B}}{gq+I}
!>    \f]
!>
!>    The first formula can be evaluated from \f$R(\psi,\theta)\f$ and \f$Z(\psi,\theta)\f$.

!> @todo
!> 2. Check poloidal flux
!>    Poloidal flux is a read-in quantity. On the other hand, it can be calculated by integrating
!>    B_pol along major radius on mid-plane.

subroutine check_m3dc1
  call check_Jacobian
  !call check_poloidal_flux
end subroutine check_m3dc1

!> See check_m3dc1 for details.
subroutine check_Jacobian

  integer:: j
  real(lk):: error(lsp, lst), maxerr
  ! calculate Jacobian using J^-1 = (Grad psi) X (Grad theta) . Grad(zeta) formula
  ! Note that J = R*|(dR_dpsi*dZ_dtheta - dR_dtheta*dZ_dpsi)|

  jsp(1,:,:) = abs((xsp(2,:,:)*zsp(4,:,:)-xsp(4,:,:)*zsp(2,:,:))*xsp(1,:,:))
  do j=1,lst
     gsp(1,:,j) = (gpsi(1,:)*qpsi(1,:) + cpsi(1,:))/(bsp(1,:,j)**2)
  enddo

  if (mod(lst,2)==0) then
     call construct_spline2d(lsp, lst, spdpsi, spdtheta, gsp, 1, 2)
     call construct_spline2d(lsp, lst, spdpsi, spdtheta, jsp, 1, 2)
  else
     call construct_spline2d(lsp, lst, spdpsi, spdtheta, gsp, 1, 0)
     call construct_spline2d(lsp, lst, spdpsi, spdtheta, jsp, 1, 0)
  endif

  !error = abs(jacobian_metric-jacobian_boozer)/jacobian_boozer
  !maxerr = maxval(error)
  !if ( maxval(error) > 1e-4 ) then
  !   print *, 'Max relative Jacobian error exceeds 1e-4, maxerr= ', maxerr
  !endif

end subroutine check_Jacobian

!> @todo
!> Implement the poloidal flux check.
!
!> see check_m3dc1 for details.
subroutine check_poloidal_flux
! Checking the self-consistency of the poloidal field and the poloidal flux
! NOT IMPLEMENTED YET
end subroutine check_poloidal_flux
!############################################################

!##########################################################################################
!> subroutine to check the status of netcdf calls
subroutine check_nf(status)
  include 'netcdf.inc'
  integer, intent ( in) :: status

  if(status /= nf_noerr) then
     write(gtcout,*) 'NetCDF ERROR:', nf_strerror(status)
     stop "Stopped"
  end if
end subroutine check_nf

!#######################################################################
! Some used short subroutines
!####################################################################
! Description:
!> Linear extrapolation using the '341' formula to estimate the derivative at boundary
!>
!> @param f
!> array that needs to be extrapolated at point n
!
!> @param n
!> the index on which the value is needed
!
!> @param dir
!> extrapolation direction
!>   - dir==1 means forward, using n-1, n-2, n-3 to extrapolate n
!>   - dir==-1 means backward, using n+1, n+2, n+3 to extrapolate n

!> @note
!> This formula is only valid for equally space function,
!> n must be greater than 3.
subroutine linear_extrapolation(f, n, dir)
  real(lk), intent(inout) :: f(:)
  integer, intent(in) :: n, dir
  if (dir==1) then
     f(n) = 2.5_lk*f(n-1) - 2.0_lk*f(n-2) + 0.5_lk*f(n-3)
  else
     f(n) = 2.5_lk*f(n+1) - 2.0_lk*f(n+2) + 0.5_lk*f(n+3)
  endif
end subroutine linear_extrapolation

!###########################################################################
!### RESONANT MAGNETIC PERTURBATION (RMP) LOADING PART #########

! Written by Lei Shi, Feb. 2, 2017

! The M3DC1 file provides the perturbed magnetic field with
! a parallel component of the perturbed vector potential.

! GTC uses the following convention
! delta_B = curl (alpha * B_0)

! where alpha is a scalar function of (psi, theta, zeta)

! M3D-C1 calculates the Fourier components of alpha in
! theta and zeta direction, alpha_mn, based on the radial component
! of delta_B calculated from M3DC1 solution.

! Detailed information about the M3DC1 calculation of alpha_mn
! can be found in Doc/Developer_Notes/M3DC1_Coordinates_and_RMP_islands.pdf

! Each M3DC1 file contains single n number data
! multiple n requires multiple M3DC1 input files
!##########################################################################

! Description
!> Main alpha loading subroutine. Read in multiple n number files
!
!> @param[in] m3dfile_base
!> Base file name, the n number will be added in the tail to form the full filename
!> like m3dfile_base\\"01.nc" and so on
!
!> @param[out] alpha_sp
!> spline coefficients array for alpha
!
!> @param[in] rescale
!> constant scaling factor to control the magnitude of the RMP fields
!
!> @callgraph
subroutine load_m3dc1_alpha(m3dfile_base, alpha_sp, rescale)
  ! Base file name, the n number will be added in the tail to form the full
  ! filename
  character(len=*), intent(in):: m3dfile_base
  ! alpha 3D spline coefficients
  real(lk), intent(out):: alpha_sp(27, lsp, lst, lszeta)
  ! artificially rescale the perturbation level
  real(lk), intent(in), optional:: rescale
  ! alpha coefficients, sum over m, and sum over n arrays
  real(lk) alpha_mn_re(ndim-1,lsp,lst-1), alpha_mn_im(ndim-1,lsp,lst-1), &
    alpha_n_re(ndim-1,lsp,lst-1), alpha_n_im(ndim-1,lsp,lst-1), &
    alpha(lsp, lst, lszeta), alpha_mn_re_tmp(lsp,lst-1), &
    alpha_mn_im_tmp(lsp,lst-1)
  ! local scale factor to store the optional rescale argument
  real(lk) local_rescale
  ! read-in poloidal and toroidal mode numbers
  integer:: m(lst-1), n(ndim-1)
  ! loop index and other  temporary variables
  integer:: i,j,k
  real(lk) r0_inv
  character(len=50):: fullfilename, n_string

  ! set the default perturbation level to be the read-in
  if(present(rescale)) then
    local_rescale = rescale
  else
    local_rescale = 1.0_lk
  endif

  ! loop through M3DC1 files and load the alpha coefficients
  ! Note that ndim==1 corresponds to 2D equilibrium only
  do i=1,ndim-1
    ! create the file name
    write(n_string, '(i2.2)') i
    fullfilename = trim(m3dfile_base)//trim(n_string)//'.nc'
    ! load single n file
    call load_singlen_alpha(fullfilename, alpha_mn_re_tmp, alpha_mn_im_tmp, m, n(i))
    do k=1, lst-1
      do j=1, lsp
        alpha_mn_re(i,j,k) = alpha_mn_re_tmp(j,k)
        alpha_mn_im(i,j,k) = alpha_mn_im_tmp(j,k)
      enddo
    enddo
  enddo

  ! Generate alpha values on equilibrium mesh based on the Fourier
  ! transformation
  call poloidal_harmonic_sum(alpha_mn_re, alpha_mn_im, alpha_n_re, alpha_n_im, &
                             m)
  call toroidal_harmonic_sum(alpha, alpha_n_re, alpha_n_im, n)

  ! Convert unit
  ! WARNING: Gauss to Tesla conversion is hardcoded. Double check the
  ! convention in your m3dc1 input file. Contact the provider of that file
  ! if you are in doubt.
  alpha = alpha * m3dc1_length_unit * 0.0001_lk

  ! Normalization
  r0_inv = 1.0_lk/r0
  ! array multiplication is faster than division
  alpha = alpha * r0_inv

  ! create 3D spline, rescale the perturbation by the chosen factor
  alpha_sp(1,:,:,:) = alpha * local_rescale
  call construct_spline3d(lsp, lst, lszeta, spdpsi, spdtheta, spdzeta, alpha_sp, &
                          1, 2, 2)
  write(gtcout,*) 'M3DC1 RMP Island loaded successfully. Values below are in &
    GTC units'
  write(gtcout,*) 'toroidal harmonics=', n, "Max alpha=", maxval(alpha), "Min alpha=", &
    minval(alpha), "rescale factor=", local_rescale

end subroutine load_m3dc1_alpha

! Description:
!> load single n alpha from one M3DC1 file
!
!> @param[in] m3dfile
!> file name containing the desired toroidal n data
!
!> @param[out] alpha_mn_re
!> real part of the read-in alpha Fourier components
!
!> @param[out] alpha_mn_im
!> imaginary part of the read-in alpha Fourier components
!
!> @param[out] m
!> poloidal mode numbers
!
!> @param[out] n
!> toroidal mode number
subroutine load_singlen_alpha(m3dfile, alpha_mn_re, alpha_mn_im, m, n)
  include 'netcdf.inc'
  character(len=*), intent(in):: m3dfile
  ! variable names in M3DC1 netcdf file.
  ! m contains the harmonic mode numbers in theta
  ! ntol is the n number chosen in the M3DC1 calculation, and is an attribute
  ! of the netcdf file
  character(len=*), parameter :: alpha_mn_re_name='alpha_real', &
    alpha_mn_im_name='alpha_imag', m_name='m', n_name='ntor'
  integer :: ncid, alpha_mn_re_id, alpha_mn_im_id, m_id

  ! temporary array for the read in data, note that the order of the indexes
  ! are exchanged in the M3DC1 file. The theta direction is the fastest
  ! changing index. Also, the float type must be single precision.
  real(4) temp2d(lst-1, lsp-2)
  integer(2) temp1d(lst-1)
  ! loop index
  integer:: i, j
  ! toroidal mode number
  integer, intent(out):: n
  ! poloidal mode numbers
  integer, intent(out):: m(lst-1)
  ! real and imaginary part of the Fourier components of alpha
  ! The psi_n=0 and psi_n=1 values will be set to zero
  real(lk), intent(out):: alpha_mn_re(lsp, lst-1), alpha_mn_im(lsp, lst-1)

  ! open m3dc1 netcdf file
  call check_nf( nf_open(m3dfile, nf_nowrite, ncid))

  ! get the ntor attribute
  call check_nf( nf_get_att(ncid, nf_global, n_name, n))
  ! get the variable ids
  call check_nf( nf_inq_varid(ncid, alpha_mn_re_name, alpha_mn_re_id))
  call check_nf( nf_inq_varid(ncid, alpha_mn_im_name, alpha_mn_im_id))
  call check_nf( nf_inq_varid(ncid, m_name, m_id))

  ! get the variables
  ! alpha_mn real
  call check_nf( nf_get_var(ncid, alpha_mn_re_id, temp2d))
  do i=2, lsp-1
    do j=1, lst-1
      alpha_mn_re(i,j) = temp2d(j, i-1)
    enddo
  enddo
  alpha_mn_re(1,:)=0
  alpha_mn_re(lsp,:)=0

  ! alpha_mn imag
  call check_nf( nf_get_var(ncid, alpha_mn_im_id, temp2d))
  do i=2, lsp-1
    do j=1, lst-1
      alpha_mn_im(i,j) = temp2d(j, i-1)
    enddo
  enddo
  alpha_mn_im(1,:)=0
  alpha_mn_im(lsp,:)=0

  ! m array
  call check_nf( nf_get_var(ncid, m_id, temp1d))
  m(:) = temp1d(:)
  
  ! close file
  call check_nf( nf_close(ncid))

end subroutine load_singlen_alpha


! Description:
!> 3D data creation from m and n harmonics
!> \f[
!>   \alpha(\psi,\theta,\zeta)= \sum_{n,m} \alpha_{mn}(\psi)
!>                              \textrm{exp}(\textrm{i}n\zeta - \textrm{i}m\theta)
!> \f]
!> First sum over m
!
!> @param[in] alphamn_re
!> real part of alpha_mn
!
!> @param[in] alphamn_im
!> imaginary part of alpha_mn
!
!> @param[out] alphan_re
!> real part after summation over m
!
!> @param[out] alphan_im
!> imaginary part after summation over m
!
!> @param[in] m
!> poloidal mode numbers
subroutine poloidal_harmonic_sum(alphamn_re, alphamn_im, alphan_re, alphan_im, m)
  integer, intent(in):: m(lst-1)
  real(lk), intent(in):: alphamn_re(ndim-1,lsp,lst-1), alphamn_im(ndim-1,lsp,lst-1)
  real(lk), intent(out):: alphan_re(ndim-1,lsp,lst-1), alphan_im(ndim-1,lsp,lst-1)
  integer i,j,k,l
  real(lk) theta, thetam, cos_thetam, sin_thetam

  alphan_re = 0.0_lk
  alphan_im = 0.0_lk
  do i=1, lst-1
     theta = real(i-1,lk)*spdtheta
     do j=1, lst-1
        thetam = theta*real(m(j),lk)
        cos_thetam = cos(thetam)
        sin_thetam = sin(thetam)
        do k=1, lsp
           do l=1, ndim-1
              alphan_re(l, k, i) = alphan_re(l,k,i)+ alphamn_re(l,k,j)*cos_thetam + &
                   alphamn_im(l,k,j)*sin_thetam
              alphan_im(l, k, i) = alphan_im(l,k,i)- alphamn_re(l,k,j)*sin_thetam + &
                   alphamn_im(l,k,j)*cos_thetam
           enddo
        enddo
     enddo
  enddo
end subroutine poloidal_harmonic_sum

! Description:
!> Sum over n. See poloidal_harmonic_sum for summation formula
!
!> @param[in] alphan_re
!> real part of alpha after sum over m
!
!> @param[in] alphan_im
!> imaginary part of alpha after sum over m
!
!> @param[out] alpha
!> RMP alpha as a function of psi,theta,zeta
!
!> @param[in] n
!> toroidal mode numbers

subroutine toroidal_harmonic_sum(alpha, alphan_re, alphan_im, n)
  integer, intent(in):: n(ndim-1)
  real(lk), intent(in):: alphan_re(ndim-1, lsp, lst-1), alphan_im(ndim-1, lsp, lst-1)
  real(lk), intent(out):: alpha(lsp, lst, lszeta)

  integer i,j,k,l
  real(lk) dzeta, zeta, zetan, cos_zetan, sin_zetan
  alpha = 0.0_lk
  do i=1, lszeta
     zeta=spdzeta*real(i-1,lk)
     do j=1, ndim-1
        zetan = zeta*real(n(j),lk)
        cos_zetan = cos(zetan)
        sin_zetan = sin(zetan)
        do k=1, lst-1
           do l=1, lsp
              alpha(l,k,i) = alpha(l,k,i) + alphan_re(j,l,k)*cos_zetan - alphan_im(j,l,k)*sin_zetan
           enddo
        enddo
     enddo
  enddo
  ! periodicity in theta
  alpha(:,lst,1:mtoroidal) = alpha(:,1,1:mtoroidal)
  ! periodicity in zeta
  alpha(:,:,mtoroidal+1) = alpha(:,:,1)
end subroutine

end module m3dc1
