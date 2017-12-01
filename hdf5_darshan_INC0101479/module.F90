! This file is part of GTC version 3 
! GTC version 3 is released under the 3-Clause BSD license:

! Copyright (c) 2002,2010,2016, GTC Team (team leader: Zhihong Lin, zhihongl@uci.edu)
! All rights reserved.

! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.

! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.

! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission.
! ==============================================================================

!> \file 
!> Common data modules and general preload modules

!> general preload modules: mpi, omp_lib, openacc, cudafor
module system_env
  use mpi
  use omp_lib
#ifdef _OPENACC
  use openacc
  use cudafor
#endif
end module system_env

!> float precision control + preload basic modules
module precision
  use system_env, only : MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_REAL,&
                         MPI_COMPLEX
  !> system double precision float length
  integer, parameter :: doubleprec=kind(1.0d0)
  !> system single precision float length
  integer, parameter :: singleprec=kind(1.0e0)
  !> system default precision float length
  integer, parameter :: defaultprec=kind(0.0)

#ifdef DOUBLE_PRECISION
  !> float precession, lk stands for Length of real Kind
  integer, parameter :: lk=doubleprec
  integer, parameter :: mpi_Rsize=MPI_DOUBLE_PRECISION
  integer, parameter ::  mpi_Csize=MPI_DOUBLE_COMPLEX
#else
  !> float precession, lk stands for Length of real Kind
  integer, parameter :: lk=singleprec
  integer, parameter :: mpi_Rsize=MPI_REAL
  integer, parameter :: mpi_Csize=MPI_COMPLEX
#endif
  !> machine accuracy for current chosen float precision
  real(lk), parameter:: machineEpsilon=10.0_lk*epsilon(1.0_lk)
  save
end module precision

module global_parameters
  use precision, only : lk
  !> gtc.out IO unit number
  integer,parameter :: gtcout=11
  !> number of (n,m) pairs accepted from gtc.in, and diagnosed
  integer,parameter :: num_modes=8

  !> \defgroup control_parameters ''control parameters''
  !> @{

  !> 0/1 switch
  integer :: track_particles
  !> 3D data output 0/1 switch
  integer :: ndata3d
  !> total perturbation field grids per toroidal plane
  integer ::  mgrid
  !> largest index of perturbation field psi grid. Note psi grid range is [0,mpsi]
  integer ::  mpsi
  !> perturbation field theta grid number on iflux surface. theta grid on other flux surfaces are rescaled according local ion temperature.
  integer ::  mthetamax
  !> perturbation field toroidal cross-section number
  integer ::  mtoroidal
  !> UNKNOWN
  integer ::  iodiag
  !> UNKNOWN
  integer ::  iodata1d
  !> number of included particle species
  integer ::  nspecies
  !> current time step
  integer ::  istep
  !> time history diagnostic frequency (in time steps)
  integer ::  ndiag
  !> total snapshot number (spatial diagnostic and restart points)
  integer ::  msnap
  !> total simulation ion time step
  integer ::  mstep
  !> accumulated diagnostic steps read from restart history
  integer ::  mstepall
  !> UNKNOWN
  integer ::  izonal
  !> radial boundary treatment control
  integer ::  nbound
  !> assymetric radial boundary control
  integer ::  nboundR
  !> run ID. 0 for fresh run, >1 for restart runs
  integer ::  irun
  !> Runge-Kutta step counter
  integer ::  irk
  !> diagnostic flag, 0 means it's diagnostic time
  integer ::  idiag
  !> thermal electron sub-steps within 1st ion RK step. Note that for second ion RK step, 2*ncyclee electron sub-steps are used.
  integer ::  ncyclee
  !> sub-cycle for fast electrons, similar to ncyclee
  integer ::  ncyclefe
  !> theta and zeta grid number used for diagnostics and mode
  !> filtering
  integer ::  mtdiag
  !> total electron recursive loop number in hybrid model
  integer ::  nhybrid
  !> current electron recursive loop
  integer ::  ihybrid
  !> ES/EM control, 0/1 switch
  integer ::  magnetic

  !> Linear/NL control, 0/1 switch
  integer ::  nonlinear

  !> mode filtering control.
  !> value   | meaning
  !> ------- | -------
  !> 0       | keep all modes
  !> 1       | select n-mode in setup.F90
  !> 2       | select n & m mode
  !> >2      | select n&m mode and \f$ k_\parallel \ll k_\perp \f$
  integer :: nfilter
  !> equilibrium input type
  !> value | meaning
  !> ----- | --------
  !> 0     | analytic
  !> 1     | EFIT
  !> 2     | VMEC
  !> 3     | M3DC1
  integer :: numereq
  !> UNKNOWN
  integer ::  neop
  !> UNKNOWN
  integer ::  neot
  !> UNKNOWN
  integer ::  neoz
  !> UNKNOWN
  integer ::  eqcurrent
  !> toroidal magnetic field and current relation, viewing from
  !> zeta=0 plane
  !> value | meaning
  !> ----- | -------
  !> 0     | B out, I out
  !> 1     | B in, I out
  !> 2     | B in, I in
  !> 3     | B out, I in
  integer ::  fielddir
  !> total number of information slots for each thermal ion marker
  integer ::  nparami
  !> information slots number for fast ion marker
  integer ::  nparamf
  !> information slots number for fast electron marker
  integer :: nparamfe
  !> information slots number for thermal electron marker
  integer :: nparame
  !> UNKNOWN
  integer :: mpsilow
  !> UNKNOWN
  integer :: mgridlow
  !> UNKNOWN
  integer :: mpsihigh
  !> UNKNOWN
  integer :: mgridhigh
  !> boundary condition formula
  !> value | meaning
  !> ----- | -------
  !> 0     | fixed zero boundary
  !> 1     | linear inner boundary
  integer :: bcond
  !> analytic equilibrium field model
  !> value | meaning
  !> ----- | -------
  !> 0     | s-alpha like (cyclone)
  !> 1     | first order (in r/R_0) model with parallel current
  integer :: fieldmodel
  !> UNKNOWN
  integer :: ilaplacian
  !> UNKNOWN
  integer :: antenna
  !> UNKNOWN
  integer :: tearingmode
  !> UNKNOWN
  integer :: myth0
  !> toroidal mode numbers of the chosen modes
  integer :: n_modes(num_modes)
  !> poloidal mode numbers of the chosen modes
  integer :: m_modes(num_modes)
  !> analytic island option, 0/1 switch
  integer :: island
  !> UNKNOWN
  integer :: fem
  !> UNKNOWN
  integer :: mgrid_fem
  !> UNKNOWN
  integer :: trilength
  !> number of parallel and perpendicular smoothing times
  integer :: ismooth
  !> input equilibrium values standard
  !> value | meaning
  !> ----- | -------
  !> 0     | on-axis values
  !> 1     | on iflux values
  integer :: inorm
  !> UNKNOWN
  integer :: irestore
  !> equilibrium rotation control
  !> value | meaning
  !> ----- | -------
  !> 0     | no rotation, no equilibrium Er
  !> 1     | Er only
  !> 2     | Rotation and Er satisfy force balance
  integer :: irotation
  !> UNKNOWN
  integer :: idiagonal
  !> output for Synthetic Diagnostics Platform, 0/1 switch
  integer :: sdp_output
  !> numerical island existence, 0/1 switch
  integer :: eq_island
  !> partial torus period
  integer :: toroidaln

  !> parallel non-linearity
  real(lk) :: paranl
  !> lower simulation boundary in psi 
  real(lk) :: psi0
  !> upper simulation boundary in psi
  real(lk) :: psi1
  !> UNKNOWN 
  real(lk) :: rg0
  !> ratio of a circle's circumference to its diameter, \f$ \pi \f$
  real(lk) :: pi
  !> \f$ 2\pi \f$ 
  real(lk) :: pi2
  !> \f$ 1/2\pi \f$
  real(lk) :: pi2_inv
  !> simulation time step in unit \f$ R_0/C_s \f$, where \f$ R_0 \f$ is the
  !> major radius, and \f$ C_s \equiv \sqrt{T_e/m_i} \f$ is the ion acoustic
  !> speed
  real(lk) :: tstep
  !> GTC length unit, equals \f$ R_0 \f$
  real(lk) :: ulength
  !> GTC time unit, equals \f$ 1/\omega_{ci} \f$, 
  !> where \f$ \omega_{ci} \equiv eB_0/m_i \f$ is the ion cyclotron frequency
  real(lk) :: utime
  !> ion gyro radius \f$ \rho_0 \equiv v_{thi}/\omega_{ci} \f$, 
  !> where \f$ v_{thi} \equiv \sqrt{T_i/m_i} \f$ is the ion thermal speed
  real(lk) :: rho0
  !> UNKNOWN 
  real(lk) :: maxwell(100001)
  !> major radius in centi-meter
  real(lk) :: r0
  !> on-axis magnetic field strength in Gauss
  real(lk) :: b0
  !> input electron temperature in eV
  real(lk) :: etemp0
  !> input electron density in 1/cm^3 
  real(lk) :: eden0
  !> UNKNOWN
  real(lk) :: hypr1
  !> UNKNOWN
  real(lk) :: hypr2
  !> UNKNOWN
  real(lk) :: omega_antenna
  !> UNKNOWN
  real(lk) :: eta
#ifdef _FRC
  !> partial torus toroidal angle range
  real(16) torbound
#else
  real(lk) torbound
#endif
  !> @}

#ifndef GPU_UM
  !$acc declare create(nparami,nparame,nparamf,nparamfe)
#endif
!> \defgroup decomposition_parameters ''MPI toroidal and particle decomposion''
!> @{

  !> current MPI process ID
  integer :: mype
  !> total MPI processes number
  integer :: numberpe
  !> MPI processes number within a single toroidal domain
  integer :: npartdom
  !> toroidal communicator, for MPI communications __between__ toroidal domains with
  !> the same particle domain rank
  integer :: toroidal_comm
  !> particle domain communicator, for communications between particle domains __within__ the same
  !> toroidal domain
  integer :: partd_comm
  !> size of the current particle domain world, should equal npartdom
  integer :: nproc_partd
  !> current rank in particle domain world, should equal
  !> particle_domain_location
  integer :: myrank_partd
  !> size of the current toroidal domain world, should equal mtoroidal
  integer :: nproc_toroidal
  !> current rank in toroidal domain world, should equal
  !> toroidal_domain_location
  integer :: myrank_toroidal
  !> toroidal rank of the MPI corresponds to myrank_toroidal-1, periodicity is
  !> considered
  integer :: left_pe
  !> toroidal rank corresponds to myrank_toroidal+1
  integer :: right_pe
  !> identifier for current PE's toroidal domain
  integer :: toroidal_domain_location
  !> identifier for current PE's particle domain
  integer :: particle_domain_location

#ifdef _OPENMP
  !> total OpenMP threads number per MPI process
  integer nthreads
#endif
!> @}

!> \defgroup restart_parameters ''restart parameters''
!> @{

! XY rolling restart

  !> restart rolling index
  integer :: irest
  !> FileExit IO unit number
  integer :: FileExit
  !> path to restart folder 1
  character(len=10) :: restart_dir1
  !> path to restart folder 2
  character(len=10) :: restart_dir2
!> @}
! keep changes made by subroutines
  save
end module global_parameters

!> equilibrium field evaluation related quantities
!> equilibrium fields are represented with quadratic spline functions on an
!> equilibrium mesh.
!> \sa spline_function
module equilibrium
  use precision,only: lk
  !> spline grid points in psi
  integer:: lsp
  !> spline grid points in theta (including the duplicate at 0 and 2pi)
  integer:: lst
  !> total spline grid points in zeta (including duplicates at 0 and 2pi)
  integer:: lszeta
  !> number of spline sections in zeta per toroidal section of simulation planes
  integer:: nzsp_sec
  !> number of toroidal modes used in 3D field (VMEC or M3DC1)
  integer:: ndim
#ifdef _TOROIDAL3D
  !> spline array dimension for 3D equilibrium 
  integer,parameter :: spdim=27 
#else
  !> spline array diemnsion for 2D equilibrium
  integer,parameter :: spdim=9 
#endif
  !> \f$ \psi_p \f$ at machine wall
  real(lk) psiw
  !> \f$ \psi_p \f$ at separatrix
  real(lk) ped
  !> \f$ \Delta \psi_p \f$ used for spline
  real(lk) spdpsi
  !> \f$ \Delta \theta \f$ used for spline
  real(lk) spdtheta
  !> \f$ \Delta \zeta \f$ used for spline
  real(lk) spdzeta
  !> \f$ \Delta rg \f$ used for spline
  real(lk) spdrg
  !> \deprecated
  !> \f$ \Delta \psi_t \f$ used for spline
  real(lk) spdtor
  !> \f$ 1/\Delta \psi_p \f$
  real(lk) spdpsi_inv
  !> \f$ 1/\Delta \theta \f$
  real(lk) spdtheta_inv
  !> \f$ 1/\Delta \zeta \f$
  real(lk) spdzeta_inv
  !> \f$ 1/\Delta rg \f$
  real(lk) spdrg_inv
  !> \f$ 1/\Delta \psi_t \f$
  real(lk) spdtor_inv

  !> \defgroup 1d_spline_arrays ''1D spline coefficients''
  !> @{

  !> UNKNOWN profile
  real(lk),dimension(:),allocatable :: stpp
  !> UNKNOWN profile
  real(lk),dimension(:),allocatable :: mipp
  !> UNKNOWN profile
  real(lk),dimension(:),allocatable :: mapp
  !> equilibrium radial E field on perturbation grids
  real(lk),dimension(:),allocatable :: mesher
  !> safety factor spline coefficients
  real(lk),dimension(:,:),allocatable :: qpsi
  !> covariant toroidal B field spline coefficients
  real(lk),dimension(:,:),allocatable :: gpsi
  !> pressure spline coefficients
  real(lk),dimension(:,:),allocatable :: ppsi
  !> minor radius spline coefficients
  real(lk),dimension(:,:),allocatable :: rpsi
  !> toroidal flux spline coefficients
  real(lk),dimension(:,:),allocatable :: torpsi
  !> electron density spline coefficients
  real(lk),dimension(:,:),allocatable :: nepp
  !> electron temperature spline coefficients
  real(lk),dimension(:,:),allocatable :: tepp
  !> ion temperature spline coefficients
  real(lk),dimension(:,:),allocatable :: tipp
  !> ion density spline coefficients
  real(lk),dimension(:,:),allocatable :: nipp
  !> fast ion temperature spline coefficients
  real(lk),dimension(:,:),allocatable :: tfpp
  !> fast ion density spline coefficients
  real(lk),dimension(:,:),allocatable :: nfpp
  !> fast electron temperature spline coefficients
  real(lk),dimension(:,:),allocatable :: tfepp
  !> fast electron density spline coefficients
  real(lk),dimension(:,:),allocatable :: nfepp
  !> effective charge profile spline coefficients
  real(lk),dimension(:,:),allocatable :: zepp
  !> plasma rotation profile spline coefficients
  real(lk),dimension(:,:),allocatable :: ropp
  !> radial electric field spline coefficients
  real(lk),dimension(:,:),allocatable :: erpp
  !> cos function spline coefficients
  real(lk),dimension(:,:),allocatable :: spcos
  !> sin function spline coefficients
  real(lk),dimension(:,:),allocatable :: spsin
  !> rg(psi_p) spline coefficients
  real(lk),dimension(:,:),allocatable :: rgpsi
  !> psi_p(rg) spline coefficients
  real(lk),dimension(:,:),allocatable :: psirg
  !> psi_p(psi_t) spline coefficients
  real(lk),dimension(:,:),allocatable :: psitor
  !> covariant poloidal B field spline coefficients
  real(lk),dimension(:,:),allocatable :: cpsi
  !> X-Z coordinates for all perturbation grids
  real(lk),dimension(:,:),allocatable :: xygrid
  !> @}

  !> \defgroup 2d_spline_arrays ''2D spline coefficients''
  !> @{

  !> B field magnitude \f$ B(\psi_p,\theta) \f$ spline
  real(lk),dimension(:,:,:),allocatable :: bsp
  !> Major radius \f$ X(\psi_p,\theta) \f$ spline
  real(lk),dimension(:,:,:),allocatable :: xsp
  !> vertical coordinate \f$ Z(\psi_p,\theta) \f$ spline
  real(lk),dimension(:,:,:),allocatable :: zsp
  !> Jacobian obtained from Boozer expression
  !> \f[
  !> \mathcal{J} = (gq + I)/B^2
  !> \f]
  !> \sa m3dc1::check_Jacobian
  real(lk),dimension(:,:,:),allocatable :: gsp
  !> Jacobian obtained from metric tensor
  !> \f[
  !> \mathcal{J} = \frac{\partial (X,Z,\Phi)}{\partial (\psi_p,\theta,\zeta)}
  !> \f]
  !> \sa m3dc1::check_Jacobian
  real(lk),dimension(:,:,:),allocatable :: jsp
  !> difference between \f$ \zeta \f$ and \f$ \Phi \f$
  !> \f[
  !> \zeta = \Phi + \nu(\psi_p, \theta, \zeta)
  !> \f]
  !> Only avaible for VMEC
  real(lk),dimension(:,:,:),allocatable :: fsp
  !> \deprecated
  !> toroidal plasma current, proportional to cpsi
  real(lk),dimension(:,:,:),allocatable :: rd
  !> \deprecated
  !> same meaning as fsp, but set to 0
  real(lk),dimension(:,:,:),allocatable :: nu
  !> \deprecated
  !> covariant radial B field component
  !> set to zero
  real(lk),dimension(:,:,:),allocatable :: dl
  !> \deprecated
  !> UNKNOWN quantity
  real(lk),dimension(:,:,:),allocatable :: ha
  !> \deprecated
  !> UNKNOWN quantity
  real(lk),dimension(:,:,:),allocatable :: hb
  !> @}

  !> \defgroup 3d_spline_arrays ''3D spline coefficients''
  !> @{

  !> spline coefficients for equilibrium islands \f$ \alpha(\psi_p,\theta,\zeta) \f$
  !> \sa m3dc1::load_m3dc1_alpha
  real(lk),dimension(:,:,:,:),allocatable :: alphasp
  !> additional scaling factor on read-in alpha
  real(lk) alpha_scale
  !> @}
#ifndef GPU_UM
  !$acc declare create(qpsi,gpsi,ropp,erpp,rgpsi,cpsi)
  !$acc declare create(bsp,xsp,mesher,alphasp)
#endif
  save
end module equilibrium

module magnetic_island
  use precision, only: lk
  integer,parameter :: l=1 !island number
  integer wi,wo,ires,qs
  integer,dimension(:),allocatable :: isn,ism

  real(lk),dimension(:,:),allocatable :: hmesh_total,hmesh_perturb,hangle,alphaIs !helical flux , alpha_zero, ksi in Eq(15), alpha in Eq(11)
  !> gradient of alphaIs
  real(lk),dimension(:,:,:),allocatable :: gradalphaIs 
#ifndef GPU_UM
  !$acc declare create(gradalphaIs,alphaIs)
#endif
  save
end module magnetic_island

module cylindricalRZ
  use precision, only: lk
  integer  lsr,lsz,nbbbs,limitr
  real(lk) rdim,zdim,rleft,psi_sep,spdR,spdZ
  real(lk),dimension(:),allocatable :: fpol,rbbbs,zbbbs,rlim,zlim
  real(lk),dimension(:,:),allocatable :: psirz
  real(lk),dimension(:,:,:),allocatable :: psirz_eq,psieq_p
  save
end module cylindricalRZ


!> particle related arrays delaration.
!>
!> All the particle arrays are defined on current MPI
module particle_array
  use precision, only: lk
  !> particle diagnostics: # of quantities per species in history.out
  integer,parameter :: mpdiag=10
  !> particle diagnostics: # of quantities per species in data1d.out
  integer,parameter :: mpdata1d=3

  !> \defgroup electron_array ''Thermal Electron Arrays''
  !> @{

  !> electron particle number
  integer :: me
  !> UNKNOWN
  integer :: me1
  !> electron particle array size 
  integer :: memax
  !> electron gyro-averaging points
  integer :: ngyroe
  !> electron load scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no kinetic particle
  !> 1     | uniform profiles
  !> 2     | real profiles, marker uniform in space, Maxwellian in velocity
  !> 3     | marker proportional to real distribution
  integer :: eload
  !> electron pitch angle filtering scheme
  !> value | meaning
  !> ----- | -------
  !> 1     | load trapped particles only
  !> not 1 | load all particles
  integer :: etrap
  !> electron collision scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no collision
  !> >0    | with collision
  integer :: ecoll
  !> electron charge
  real(lk) :: qelectron
  !> electron mass
  real(lk) :: aelectron
  !> electron beta
  real(lk) :: betae
  !> electron-electron collison time
  real(lk) :: tauee
  !> electron-ion collision time
  real(lk) :: tauei
  !> electron time diagnostic array
  !> \n See also: diagnosis, PushParticle
  real(lk),dimension(mpdiag) :: diagelectron
  !> electron location array, 
  !> index of the lower theta grid on inner flux surface
  integer,dimension(:,:),allocatable :: jtelectron0
  !> electron location array, 
  !> index of the lower theta grid on outer flux surface
  integer,dimension(:,:),allocatable :: jtelectron1
  !> electron location array, 
  !> linear interpolation weight in zeta
  real(lk),dimension(:),allocatable :: wzelectron
  !> electron temperature profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshte
  !> electron density profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshne
  !> electron \f$ \kappa_{ne} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{ne} \equiv -\frac{1}{n_e} \frac{d n_e}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapane
  !> electron \f$ \kappa_{te} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{te} \equiv -\frac{1}{T_e} \frac{d T_e}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapate
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: dtndpsi
  !> UNKNOWN 
  real(lk),dimension(:),allocatable :: pmarke
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonale
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalce
  real(lk),dimension(:),allocatable :: markere
  real(lk),dimension(:),allocatable :: rdteme
  real(lk),dimension(:),allocatable :: pfluxe
  real(lk),dimension(:),allocatable :: markeret
  !> electron particle information array. \n
  !> The first index indicates the different information for each marker. \n
  !> Their meanings are:
  !> index | meaning
  !> ----- | ----------
  !> 1     | psi
  !> 2     | theta
  !> 3     | zeta
  !> 4     | rho_para
  !> 5     | weight
  !> 6     | sqrt(mu)
  !> 7     | f/g
  !> 8,9   | particle tracking labels
  real(lk),dimension(:,:),allocatable :: zelectron
  !> temporary container for initial particle informations before each
  !> ion Runge-Kutta push. See zelectron for the format.
  real(lk),dimension(:,:),allocatable :: zelectron0
  !> temporary container for initial particle informations before each
  !> electron sub Runge-Kutta push. See zelectron for the format.
  real(lk),dimension(:,:),allocatable :: zelectron1
  !> psi weight array used for perturbed field evaluation at particle location.
  !> weight for outer flux surface.
  !> \n Shape: (ngyroe,me)
  real(lk),dimension(:,:),allocatable :: wpelectron
  !> theta weight on upper theta grid point on inner flux surface.
  !> \n Shape: (ngyroe,me)
  real(lk),dimension(:,:),allocatable :: wtelectron0
  !> theta weight on upper theta grid point on outer flux surface.
  !> \n Shape: (ngyroe,me)
  real(lk),dimension(:,:),allocatable :: wtelectron1
  !> electron number density on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that density is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: densitye
  !> electron parallel flow field on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that flowe is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: flowe
  !> 1D diagnostic data for electron
  real(lk),dimension(:,:),allocatable :: data1de
  !> parallel non-adiabatic electron pressure
  real(lk),dimension(:,:),allocatable :: pressureepara
  !> perpendicular non-adiabatic electron pressure
  real(lk),dimension(:,:),allocatable :: pressureeperp
  !> @}
#ifndef GPU_UM
  !$acc declare create(diagelectron,jtelectron0,jtelectron1)
  !$acc declare create(wzelectron,meshte,meshne,kapane,kapate,rdteme)
  !$acc declare create(meshte,meshne,kapane,kapate)
  !$acc declare create(zelectron,zelectron0,zelectron1,wpelectron,wtelectron0)
  !$acc declare create(wtelectron1,densitye,flowe,data1de,phit,dnet,pressureepara,pressureeperp)
#endif
  !> \f$ \frac{\partial \Phi_{eff}}{\partial t}\f$
  !> \n See hybrid_model documents for details
  real(lk),dimension(:,:),allocatable :: phit
  !> UNKNOWN
  real(lk),dimension(:,:,:),allocatable :: phisave
  !> \f$ \frac{\partial n_e}{\partial t}\f$
  !> \n See hybrid_model documents for details
  real(lk),dimension(:,:),allocatable :: dnet
  !> UNKNOWN
  real(lk),dimension(:,:,:),allocatable :: dnesave
  !> UNKNOWN
  real(lk) tfracn
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: tfrac

  !> \defgroup ion_arrays ''Thermal Ion Arrays''
  !> @{
  ! themal ion (main ion)
  
  !> ion particle number
  integer :: mi
  !> ion particle array size 
  integer :: mimax
  !> ion gyro-averaging points
  integer :: ngyroi
  !> ion load scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no kinetic particle
  !> 1     | uniform profiles
  !> 2     | real profiles, marker uniform in space, Maxwellian in velocity
  !> 3     | marker proportional to real distribution
  integer :: iload
  !> ion collision scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no collision
  !> >0    | with collision
  integer :: icoll
  !> ion charge
  real(lk) :: qion
  !> ion mass
  real(lk) :: aion
  !> ion-ion collision time
  real(lk) :: tauii
  !> ion time diagnostic array
  !> \n See also: diagnosis, PushParticle
  real(lk),dimension(mpdiag) :: diagion
  !> ion location array, 
  !> index of the lower theta grid on inner flux surface
  integer,dimension(:,:),allocatable :: jtion0
  !> index of the lower theta grid on outer flux surface
  integer,dimension(:,:),allocatable :: jtion1
  !> linear interpolation weight in zeta
  real(lk),dimension(:),allocatable :: wzion
  !> ion temperature profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshti
  !> ion density profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshni
  !> ion \f$ \kappa_{ni} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{ni} \equiv -\frac{1}{n_i} \frac{d n_i}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapani
  !> ion \f$ \kappa_{ti} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{ti} \equiv -\frac{1}{T_i} \frac{d T_i}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapati
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: jacobianpsi
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pmarki
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonali
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalci
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markeri
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: rdtemi
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pfluxi
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markerit
  !> ion particle information array. \n
  !> The first index indicates the different information for each marker. \n
  !> Their meanings are:
  !> index | meaning
  !> ----- | ----------
  !> 1     | psi
  !> 2     | theta
  !> 3     | zeta
  !> 4     | rho_para
  !> 5     | weight
  !> 6     | sqrt(mu)
  !> 7     | f/g
  !> 8,9   | particle tracking labels
  real(lk),dimension(:,:),allocatable :: zion
  !> temporary container for initial particle informations before each
  !> ion Runge-Kutta push. See zion for the format.
  real(lk),dimension(:,:),allocatable :: zion0
  !> psi weight array used for perturbed field evaluation at particle location.
  !> weight for outer flux surface.
  !> \n Shape: (ngyroi,mi)
  real(lk),dimension(:,:),allocatable :: wpion
  !> theta weight on upper theta grid point on inner flux surface.
  !> \n Shape: (ngyroi,mi)
  real(lk),dimension(:,:),allocatable :: wtion0
  !> theta weight on upper theta grid point on outer flux surface.
  !> \n Shape: (ngyroi,mi)
  real(lk),dimension(:,:),allocatable :: wtion1
  !> ion number density on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that density is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: densityi
  !> ion parallel flow field on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that flowi is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: flowi
  !> 1D diagnostic data for ion 
  real(lk),dimension(:,:),allocatable :: data1di
  !> @}
#ifndef GPU_UM
  !$acc declare create(diagion,jtion0,jtion1)
  !$acc declare create(wzion,meshti,meshni,kapani,kapati,rdtemi)
  !$acc declare create(zion,zion0,wpion,wtion0,wtion1,densityi,flowi,data1di)
#endif

  !> \defgroup fast_ion_arrays ''Fast Ion Arrays''
  !> @{
  ! fast ion
  
  !> fast ion particle number
  integer :: mf
  !> fast ion particle array size 
  integer :: mfmax
  !> fast ion gyro-averaging points
  integer :: ngyrof
  !> fast ion load scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no kinetic particle
  !> 1     | uniform profiles
  !> 2     | real profiles, marker uniform in space, Maxwellian in velocity
  !> 3     | marker proportional to real distribution
  integer :: fload
  !> fast ion charge
  real(lk) :: qfast
  !> fast ion mass
  real(lk) :: afast
  !> fast ion slowing down distribution parameters
  !> UNKNOWN
  real(lk) :: sd_v0
  !> UNKNOWN
  real(lk) :: sd_vc
  !> UNKNOWN
  real(lk) :: sd_l0
  !> UNKNOWN
  real(lk) :: sd_widthInv
  !> fast ion time diagnostic array
  !> \n See also: diagnosis, PushParticle
  real(lk),dimension(mpdiag) :: diagfast
  !> fast ion location array, 
  !> index of the lower theta grid on inner flux surface
  integer,dimension(:,:),allocatable :: jtfast0
  !> index of the lower theta grid on outer flux surface
  integer,dimension(:,:),allocatable :: jtfast1
  !> linear interpolation weight in zeta
  real(lk),dimension(:),allocatable :: wzfast
  !> fast ion temperature profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshtf
  !> fast ion density profile on simulation psi grid
  real(lk),dimension(:),allocatable :: meshnf
  !> fast ion \f$ \kappa_{nf} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{nf} \equiv -\frac{1}{n_f} \frac{d n_f}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapanf
  !> fast ion \f$ \kappa_{tf} \f$ profile on simulation psi grid
  !> \f[
  !> \kappa_{tf} \equiv -\frac{1}{T_f} \frac{d T_f}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapatf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pmarkf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalcf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markerf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: rdtemf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pfluxf
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markerft
  !> fast ion particle information array. \n
  !> The first index indicates the different information for each marker. \n
  !> Their meanings are:
  !> index | meaning
  !> ----- | ----------
  !> 1     | psi
  !> 2     | theta
  !> 3     | zeta
  !> 4     | rho_para
  !> 5     | weight
  !> 6     | sqrt(mu)
  !> 7     | f/g
  !> 8,9   | particle tracking labels
  real(lk),dimension(:,:),allocatable :: zfast
  !> temporary container for initial particle informations before each
  !> ion Runge-Kutta push. See zfast for the format.
  real(lk),dimension(:,:),allocatable :: zfast0
  !> psi weight array used for perturbed field evaluation at particle location.
  !> weight for outer flux surface.
  !> \n Shape: (ngyrof,mf)
  real(lk),dimension(:,:),allocatable :: wpfast
  !> theta weight on upper theta grid point on inner flux surface.
  !> \n Shape: (ngyrof,mf)
  real(lk),dimension(:,:),allocatable :: wtfast0
  !> theta weight on upper theta grid point on outer flux surface.
  !> \n Shape: (ngyrof,mf)
  real(lk),dimension(:,:),allocatable :: wtfast1
  !> fast ion number density on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that density is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: densityf
  !> fast ion parallel flow field on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that flow is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: flowf
  !> 1D diagnostic data for fast ion
  real(lk),dimension(:,:),allocatable :: data1df
  !> @}
#ifndef GPU_UM
  !$acc declare create(diagfast,jtfast0,jtfast1)
  !$acc declare create(wzfast,meshtf,meshnf,kapanf,kapatf,rdtemf)
  !$acc declare create(zfast,zfast0,wpfast,wtfast0,wtfast1,densityf,flowf,data1df)
#endif

  !> \defgroup fast_electron_arrays ''Fast Electron Arrays''
  !> @{
!fast electron

  !> fast electron particle number
  integer :: mfe
  !> UNKNOWN
  integer :: mfe1
  !> fast electron particle array size 
  integer :: mfemax
  !> fast electron gyro-averaging points
  integer :: ngyrofe
  !> fast electron load scheme
  !> value | meaning
  !> ----- | -------
  !> 0     | no kinetic particle
  !> 1     | uniform profiles
  !> 2     | real profiles, marker uniform in space, Maxwellian in velocity
  !> 3     | marker proportional to real distribution
  integer :: feload
  integer :: fetrap
  !> fast electron charge
  real(lk) :: qfaste
  !> fast electron mass
  real(lk) :: afaste
  !> fast electron time diagnostic array
  !> \n See also: diagnosis, PushParticle
  real(lk),dimension(mpdiag) :: diagfaste
  !> fast electron locatelectron array, 
  !> index of the lower theta grid on inner flux surface
  integer,dimension(:,:),allocatable :: jtfaste0
  !> index of the lower theta grid on outer flux surface
  integer,dimension(:,:),allocatable :: jtfaste1
  !> linear interpolation weight in zeta
  real(lk),dimension(:),allocatable :: wzfaste
  !> fast electron temperature profile on simulatelectron psi grid
  real(lk),dimension(:),allocatable :: meshtfe
  !> fast electron density profile on simulatelectron psi grid
  real(lk),dimension(:),allocatable :: meshnfe
  !> fast electron \f$ \kappa_{nfe} \f$ profile on simulatelectron psi grid
  !> \f[
  !> \kappa_{nfe} \equiv -\frac{1}{n_{fe}} \frac{d n_{fe}}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapanfe
  !> fast electron \f$ \kappa_{tfe} \f$ profile on simulatelectron psi grid
  !> \f[
  !> \kappa_{tfe} \equiv -\frac{1}{T_{fe}} \frac{d T_{fe}}{d\psi_p}
  !> \f]
  real(lk),dimension(:),allocatable :: kapatfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pmarkfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: zonalcfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markerfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: rdtemfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: pfluxfe
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: markerfet
  !> fast electron particle informatelectron array. \n
  !> The first index indicates the different information for each marker. \n
  !> Their meanings are:
  !> index | meaning
  !> ----- | ----------
  !> 1     | psi
  !> 2     | theta
  !> 3     | zeta
  !> 4     | rho_para
  !> 5     | weight
  !> 6     | sqrt(mu)
  !> 7     | f/g
  !> 8,9   | particle tracking labels
  real(lk),dimension(:,:),allocatable :: zfaste
  !> temporary container for initial particle informations before each
  !> ion Runge-Kutta push. See zfaste for the format.
  real(lk),dimension(:,:),allocatable :: zfaste0
  !> temporary container for initial particle informations before each
  !> electron sub Runge-Kutta push. See zfaste for the format.
  real(lk),dimension(:,:),allocatable :: zfaste1
  !> psi weight array used for perturbed field evaluation at particle location.
  !> weight for outer flux surface.
  !> \n Shape: (ngyrofe,mfe)
  real(lk),dimension(:,:),allocatable :: wpfaste
  !> theta weight on upper theta grid point on inner flux surface.
  !> \n Shape: (ngyrofe,mfe)
  real(lk),dimension(:,:),allocatable :: wtfaste0
  !> theta weight on upper theta grid point on outer flux surface.
  !> \n Shape: (ngyrofe,mfe)
  real(lk),dimension(:,:),allocatable :: wtfaste1
  !> fast electron number density on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that density is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: densityfe
  !> fast electron parallel flow field on grid. \n
  !> Shape: (0:1,mgrid)
  !> \n first dimension stands for the two zeta planes. Note that flow is
  !> only fully collected on zeta1 plane for field solving. zeta0 plane only
  !> temporarily stores contributions from current MPI.
  real(lk),dimension(:,:),allocatable :: flowfe
  !> 1D diagnostic data for fast ion
  real(lk),dimension(:,:),allocatable :: data1dfe
  !> UNKNOWN
  real(lk) :: fetfracn
  !> UNKNOWN
  real(lk),dimension(:),allocatable :: fetfrac
  !> @}
#ifndef GPU_UM
  !$acc declare create(diagfaste,jtfaste0,jtfaste1)
  !$acc declare create(wzfaste,meshtfe,meshnfe,kapanfe,kapatfe,rdtemfe)
  !$acc declare create(zfaste,zfaste0,zfaste1,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,data1dfe)
#endif
  save
end module particle_array

module particle_tracking
  use precision, only: lk
  real(lk),dimension(:,:),allocatable :: ptrackedi,ptrackede,ptrackedf
  integer,dimension(3) :: ntrackp
  save
end module particle_tracking

module field_array
  use precision, only: lk
! PIC global fieldline following mesh
  !> The distance between two neibour radial grids for \f$ \theta=0 \f$
  real(lk) deltar
  !> The angle distance between two neibour poloidal plane 
  real(lk) deltaz
  !> The toroidal angle of the poloidal plane for current MPI
  real(lk) zeta1
  !> The toroidal angle of the poloidal plane for previous MPI
  real(lk) zeta0
  !> The number of the grid index shift for field aligned mesh on each poloidal plane
  integer,dimension(:),allocatable :: itran
  !> The grid index array for each poloidal plane (for finite difference solver)
  integer,dimension(:),allocatable :: igrid
  !> The poloidal grid number for each flux surfaces
  integer,dimension(:),allocatable ::  mtheta
  !> UNKNOW 
  integer,dimension(:),allocatable ::  nshift
  !> The grid index array for each poloidal plane (for finite element solver)
  integer,dimension(:),allocatable ::  igrid_fem
  !> \f$\delta\psi\f$ between each two neigour flux surfaces
  real(lk),dimension(:),allocatable :: deltap
  !> \f$\delta\theta\f$ between two neibour poloidal grids
  real(lk),dimension(:),allocatable ::  deltat 
  !> safety factor calculate by using the field aligned mesh
  real(lk),dimension(:),allocatable ::  qtinv
  !> realistic safety factor
  real(lk),dimension(:),allocatable ::  qmesh
  !> magnetic field on the simulation grid position
  real(lk),dimension(:),allocatable ::  bmesh
  !> jacobian value on the simulation grid position
  real(lk),dimension(:),allocatable ::  jmesh
  !> poloidal \psi value on the simulation grid position
  real(lk),dimension(:),allocatable ::  psimesh
  !> toroidal \psi value on the simulation grid position
  real(lk),dimension(:),allocatable ::  tormesh
  !> the effective charge on the simulation grid (for collision operator)
  real(lk),dimension(:),allocatable ::  meshze
#ifndef GPU_UM
  !$acc declare create(igrid,mtheta)
  !$acc declare create(deltap,deltat,qtinv,psimesh)
#endif
  !> flux surface averged <g^{\psi\psi} (for zonal solver)
  real(lk),dimension(:),allocatable :: gpsi200
  !> flux surface averaged b*b
  real(lk),dimension(:),allocatable :: b2m00
  !> UNKNOWN (output diagnosis related)
  real(lk),dimension(:),allocatable :: wtshift
  !> metric tensor element \f$ g^{\psi\psi} \f$ on the simulation grid position
  real(lk),dimension(:),allocatable :: gupp
  !> metric tensor element \f$ g^{\psi\theta} \f$ on the simulation grid position
  real(lk),dimension(:),allocatable :: gupt
  !> metric tensor element g^{\theta\theta} on the simulation grid position 
  real(lk),dimension(:),allocatable :: gutt
  !> metric tensor element g^{\zeta\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: guzz
  !> metric tensor element g^{\psi\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gupz
  !> metric tensor element g^{\theta\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gutz
  !> ion gyro radius value on the simulation grid position for the gyro average
  real(lk),dimension(:),allocatable :: rhom
  !> metric tensor element g_{\psi\psi} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdpp
  !> metric tensor element g_{\psi\theta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdpt
  !> metric tensor element g_{\theta\theta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdtt
  !> metric tensor element g_{\zeta\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdzz
  !> metric tensor element g_{\psi\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdpz
  !> metric tensor element g_{\theta\zeta} on the simulation grid position
  real(lk),dimension(:),allocatable :: gdtz
  !> toroidal mode number for the semi-spectral solver of FRC in magnetic coordinates
  real(lk),dimension(:),allocatable :: spectraln

! fields on mesh: phi, apara, fluidne, fluidue, zonal and gradient
  !> zonal electrostatic potential
  real(lk),dimension(:),allocatable :: phi00
  !> zonal radial electic field
  real(lk),dimension(:),allocatable :: phip00
  !> zonal parallel vector potential
  real(lk),dimension(:),allocatable :: apara00
  !> zonal parallel vector potential from Zhixuan's new term:\f$ \delta B_{\perp}\cdot\nabla\phi \f$
  real(lk),dimension(:),allocatable :: apara00nl
  !> initial value of apara00nl for rk2 advance algorithm
  real(lk),dimension(:),allocatable :: apara00nl0
  !> zonal electron fluid density
  real(lk),dimension(:),allocatable :: fluidne00
  !> external high k modes damping term in electron continuity equation
  real(lk),dimension(:),allocatable :: d4fluidne
  !> laplacian of \f$\delta A_{||} \f$
  real(lk),dimension(:),allocatable :: d2apara
  !> electrostatic potential
  real(lk),dimension(:,:),allocatable :: phi
  !> parallel vector potential for time advance in Ohm's law
  real(lk),dimension(:,:),allocatable :: apara
  !> fluid electron density for time advance in electron continuity equation
  real(lk),dimension(:,:),allocatable :: fluidne
  !> electron perturbed parallel velocity calculated by inverting Ampere's law, \f$n_0\delta u_e/B_0\f$
  real(lk),dimension(:,:),allocatable :: fluidue
  !> initial value of apara for rk2 time advance algorithm
  real(lk),dimension(:,:),allocatable :: apara0
  !> initial value of fluidne for rk2 time advance algorithm 
  real(lk),dimension(:,:),allocatable :: fluidne0
  !> perturbed poloidal flux \f$\delta\psi\f$ in Eq.(15) of Holod 2009pop
  real(lk),dimension(:,:),allocatable :: deltapsi
  !> initial value of deltapsi for rk2 time advance algrithm(Eq. (18) in Holod 2009pop)
  real(lk),dimension(:,:),allocatable :: deltapsi0
  !> Smoothed deltapsi for particle pusher
  real(lk),dimension(:,:),allocatable :: sdeltapsi
  !> Smoothed apara for particle pusher (both zonal and nonzonal components)
  real(lk),dimension(:,:),allocatable :: sapara
  !> Smoothed fluidne for solving Poisson's equation
  real(lk),dimension(:,:),allocatable :: sfluidne
  !> Smoothed apara for calculating laplacian \f$\nabla^2\delta A_{||}\f$  (only nonzonal component)
  real(lk),dimension(:,:),allocatable :: sdelapara
  !> MHD displacement for diagnosis and comparision with NOVA-K(Zhixuan2015pop)
  real(lk),dimension(:,:),allocatable :: MHDxi_mag
  !> save phi adiabatic to be corrected in hybrid loop
  real(lk),dimension(:,:),allocatable :: phi_zero
  !> gradient of the electrostatic potential
  real(lk),dimension(:,:,:),allocatable :: gradphi
  !> gradient of parallel vector potential
  real(lk),dimension(:,:,:),allocatable :: gradapara
  !> gradient of fluidue(note fluidue is \f$n_0\delta u_e/B_0\f$)
  real(lk),dimension(:,:,:),allocatable :: gradue
  !> gradient of electron fluid density(fluidne)
  real(lk),dimension(:,:,:),allocatable :: gradne
  !> gradient of internal antenna fields with zero boundary condition
  real(lk),dimension(:,:,:),allocatable :: gradext
  !> gradient of smoothed perturbed poloidal flux: sdeltapsi
  real(lk),dimension(:,:,:),allocatable :: gradpsi
  !> gradient of \f$\phi_{eff}\f$
  real(lk),dimension(:,:,:),allocatable :: gradphieff
  !> gradient of high order kinetic electron density: densitye
  real(lk),dimension(:,:,:),allocatable :: gradkne
  !> gradient of parallel perturbed electron pressure
  real(lk),dimension(:,:,:),allocatable :: gradpepara
  !> gradient of perpendicular perturbed electron pressure
  real(lk),dimension(:,:,:),allocatable :: gradpeperp
  !> output by Zhixuan(Zhixuan2015pop Eq.(A3))
  real(lk),dimension(:,:,:),allocatable :: gradgaugef
  !> output by Zhixuan
  real(lk),dimension(:,:,:),allocatable :: MHDxiperp
#ifndef GPU_UM
  !$acc declare create(phip00,sapara)
  !$acc declare create(gradphi,gradapara,gradpsi,gradphieff)
#endif

! external field
  !> external antenna frequency 
  real(lk),dimension(:),allocatable :: omega
  !> external antenna potential field
  real(lk),dimension(:,:),allocatable :: phiext
  !> laplacian of phiext for non-zero boundary condition
  real(lk),dimension(:,:),allocatable :: dn_ext

! diagnostics and filtering
  !> diagnosis flux surface index
  integer iflux
  !> number of diagnosis modes
  integer modes
  !> FRC solver methd control(0:radial fluxtube,1:radial global)
  integer solvermethod
  !> diagnosis toroidal mode value 
  integer,dimension(:),allocatable :: nmodes
  !> diagnosis poloidal mode value
  integer,dimension(:),allocatable :: mmodes
  !> selected toroidal mode number for FRC spectral solver in magnetic coordinates
  integer nmode

! radial interpolation
  !> lower poloidal grid index on flux surfaces in laplacian and smooth
  integer,dimension(:,:),allocatable :: jtp1
  !> UNKNOWN
  integer,dimension(:,:),allocatable :: jtp2
  !> the weight of upper poloidal grid index on flux surfaces in laplacian and smooth
  real(lk),dimension(:,:),allocatable :: wtp1
  !> UNKNOWN
  real(lk),dimension(:,:),allocatable :: wtp2

! laplacian coefficients
  !> the grid number used for calculating tokamak laplacian on each points
  integer mindexlap
  !> the grid number used for calculating FRC laplacian on each points
  integer mindexlap2
  !> the grid number used for calculating tokamak laplacian by fem
  integer mindex_fem
  !> UNKNOWN(Yong2015pop)
  integer,dimension(:),allocatable :: nindexlap
  integer,dimension(:),allocatable :: nindexlap2
  integer,dimension(:),allocatable :: nindex_fem
  integer,dimension(:,:),allocatable :: indexplap
  integer,dimension(:,:),allocatable :: indexlap2
  integer,dimension(:,:),allocatable :: indexp_fem
  integer,dimension(:,:),allocatable :: trilist
  real(lk),dimension(:,:),allocatable :: lapmat
  real(lk),dimension(:,:),allocatable :: lapmat2
  real(lk),dimension(:,:),allocatable :: lmat
  real(lk),dimension(:,:),allocatable :: dmat

! theta of max and min b-field for particle boundary cond.
  !> maximum b field on inner and outer simulation flux surfaces
  real(lk) maxbfield(0:1)
  !> mininum b field on inner and outer simulation flux surfaces
  real(lk) minbfield(0:1)
  !> theta index corresponding to the mininum b field on inner and outer simulation flux surfaces
  real(lk) thetabmin(0:1)
  !> theta index corresponding to the maximum b field on inner and outer simulation flux surfaces
  real(lk) thetabmax(0:1)
  !> theta index array corresponding to b field bins in increasing theta direction from thetabmin to thetabmax
  integer,dimension(:,:),allocatable :: thetaupp
  !> theta index array corresponding to b field bins in decreasing theta direction from thetabmin to thetabmax
  integer,dimension(:,:),allocatable :: thetadwn

! gyro averaging
  real(lk),dimension(:,:),allocatable :: pgyro
  real(lk),dimension(:,:),allocatable :: tgyro
  real(lk),dimension(:,:),allocatable :: pgyro2
  real(lk),dimension(:,:),allocatable :: tgyro2
#ifndef GPU_UM
  !$acc declare create(pgyro,tgyro,pgyro2,tgyro2)
#endif
  save
end module field_array

module petsc_array
  use precision, only: lk
  integer newcomm,nproc_newcomm,myrank_newcomm
  integer,dimension(:),allocatable :: userp,users,luserp,lusers,luserp2,lusers2
  real(lk),dimension(:),allocatable :: usera,userb,userx,lusera,luserb,luserx,lusera2,luserb2,luserx2
  save
end module petsc_array

module data_type
  integer,parameter :: bp_char=0,bp_short=1,&
                bp_int=2,bp_long=3,bp_longlong=4,bp_float=5,&
                bp_double=6,bp_longdouble=7,bp_pointer=8,&
                bp_string=9,bp_uchar=50,bp_ushort=51,bp_uint=52,&
                bp_ulong=53,bp_ulonglong=54
  save
end module data_type

module LAMYRIDGEeq
  use precision, only: lk
  integer numr,numz,nrzgrid
  integer,dimension(:),allocatable :: igrid_frc,nz_frc

  real(lk) zmax,psi0_frc,psi1_frc,Zsimulation0,Zsimulation1,min_r,&
           lz,delta_Z,deltapsi_frc
  real(lk),dimension(:),allocatable :: r_Eq,z_Eq,psiEq0,meshni_cy,meshne_cy,&
                                       meshti_cy,meshte_cy,psi_frc,psifrc_tmp,&
                                       deltar_frc,dpsifrcdr,dpsifrcdz
  real(lk),dimension(:,:),allocatable ::RZ_reg
  real(lk),dimension(:,:,:),allocatable :: B0mag_Eq,B0R_Eq,B0Z_Eq,lrne_Eq,lrni_Eq,lrti_Eq,lrte_Eq
  save
end module LAMYRIDGEeq

module interfaces
  interface
    subroutine push(species_name,icycle,irksub,ihybrid)
      implicit none
      character(*),intent(in) :: species_name
      integer,intent(in),optional :: icycle,irksub,ihybrid
    end subroutine push

    subroutine axisExtrapolate(farray)
      use precision, only: lk
      use global_parameters,only: mgrid

      real(lk),dimension(mgrid) :: farray
    end subroutine axisExtrapolate

  end interface
end module interfaces
