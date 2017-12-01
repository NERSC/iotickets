# 1 "PushParticle.F90"
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



# 26


!> \brief Unified particle push subroutine.
!
!> This subroutine is called through the wrapper subroutine \c push in file \c push.F90.\n
!
!> The particle information array \c zpart is updated based on the particle
!! species, Runge-Kutta step, push type ( \c hybrid or \c gyrokinetic ),
!> subcycling setting, and plasma field evaluated at each particle location.\n
!>\n
!>
!> Essentially, 5 particle quantities are pushed:
!>  - psi ( \c zpart(1) ): poloidal flux, effective radial coordinate
!>  - theta ( \c zpart(2) ): poloidal angle
!>  - zeta ( \c zpart(3) ): toroidal angle
!>  - rho_para ( \c zpart(4) ): v_para/omega_c, effective parallel gyroradius
!>  - weight ( \c zpart(5) ): distribution function carried by the marker
!>  .
!> \n
!>
!> A push consists of the following parts:
!>  -# set time step based on specific push type
!>    - First ion RK ( \c irk==1 ): half ion step
!>    - Second ion RK ( \c irk==2 ): full ion step
!>    - First electron RK ( \c irksub==1 ): half electron step
!>    - Second electron RK ( \c irksub==2 ): full electron step
!>  -# keep temporary record of initial values if in RK 1
!>  -# evaluate perturbed fields using the chosen gyro-averaging scheme
!>  -# evaluate equilibrium fields at particle GC location
!>  -# calculate the time derivatives of the 5 quantities
!>  -# update the 5 quantities
!>  -# deal with out of boundary particles and restore profiles
!>  -# calculate time diagnostics
!>  .
!> \n
!>
!> @param[in,out] zpart
!> particle information array, stands for zion,zelectron,zfast, or zfaste
!
!> @param[in,out] zpart0
!> initial values of zpart, set in RK 1, used in RK 2
!
!> @param[in] wppart
!> linear interpolation weight in psi for gyro averaging locations, calculated
!> in subroutine locateParticle
!
!> @param[in] wtpart0
!> linear interpolation weight for upper theta grid on inner flux surface
!
!> @param[in] wtpart1
!> linear interpolation weight for upper theta grid on outer flux surface
!
!> @param[out] data1d
!> 1D diagnostic data
!
!> @param[in] jtpart0
!> index of the lower theta grid on inner flux surface
!
!> @param[in] jtpart1
!> index of the lower theta grid on outer flux surface
!
!> @param[out] diagpart
!> time diagnostic data
!
!> @param[in] wzpart
!> linear interpolation weight in zeta
!
!> @param[in] rdtem
!> UNKNOWN
!
!> @param[in] kapan
!> density gradient term:\f$ \frac{1}{n}\frac{dn}{d\psi} \f$
!
!> @param[in] kapat
!> temperature gradient term:\f$ \frac{1}{T}\frac{dT}{d\psi} \f$
!
!> @param[in] meshn
!> equilibrium density on simulation grid
!
!> @param[in] mesht
!> equilibrium temperature on simulation grid
!
!> @param[in] pflux
!> UNKNOWN
!
!> @param[in] qpart
!> particle charge
!
!> @param[in] apart
!> particle mass
!
!> @param[in] pload
!> particle load scheme, see iload for details.
!
!> @param[in] ngyro
!> gyro-averaging scheme, normally either 1 or 4
!
!> @param[in] mpmax
!> maximum particle number on single MPI, determines array size
!
!> @param[in] mp
!> actual particle number on current MPI
!
!> @param[in] nparam
!> number of particle information entries, different kinds of particles may
!> have different number of information entries in zpart array.
!
!> @param[in] ncyclep
!> total subcycling number
!
!> @param[in] icycle
!> current subcycle count
!
!> @param[in] irksub
!> current RK step in a subcycle
!
!> @param[in,out] mp1
!> UNKNOWN
!
!> @param[in,out] zpart1
!> UNKNOWN
!
!> @param[in] ihybrid
!> current hybrid loop number
!
!> @param[in] phit
!> time derivative of electic potential used in hybrid model:
!> \f$ \frac{\partial \phi}{\partial t} \f$
!
!> @param[in] dnt
!> time derivative of density used in hybrid model:
!> \f$ \frac{\partial n}{\partial t}\f$
subroutine hybridPushParticle(zpart,zpart0,wppart,wtpart0,wtpart1,&
    data1d,jtpart0,jtpart1,diagpart,wzpart,rdtem,kapan,kapat,meshn,mesht,&
    pflux,qpart,apart,pload,ngyro,mpmax,mp,nparam,ncyclep,icycle,irksub,mp1,&
    zpart1,ihybrid,phit,dnt)
  use system_env
  use precision, only: lk, mpi_Rsize
  use global_parameters,only: nonlinear,irk,tstep,magnetic,rg0,eqcurrent,&
    paranl,psi1,pi,pi2,psi0,idiag,mpsi,mpsilow,island,gtcout,mype,istep,&
    torbound,rho0,track_particles,nhybrid,numereq,mthetamax,irestore
  use particle_array,only: mpdiag,mpdata1d,meshne,meshte,kapane,betae,&
    sd_vc,sd_v0,sd_widthInv,sd_l0
  use field_array,only: deltar,deltap,deltat,deltaz,gradphi,gradapara,&
    gradphieff,sapara,psimesh,nmodes,nmode,zeta0,phip00,solvermethod,gradpsi,&
    maxbfield,minbfield,thetabmin,thetabmax,thetaupp,thetadwn
  use equilibrium,only:&
    lsp,spdpsi_inv,spdpsi,lst,spdtheta_inv,spdtheta,lszeta,spdzeta,spdzeta_inv,&
    nzsp_sec,qpsi,gpsi,cpsi,ropp,erpp,rgpsi,bsp,xsp,spdim,mesher,alphasp,ndim
  use magnetic_island,only: gradalphaIs,alphaIs
  implicit none

!declaration of the dummy arguments
  integer pload,ngyro,mpmax,mp,nparam
  integer,optional :: ncyclep,icycle,irksub,ihybrid,mp1
  integer,dimension(:,:) :: jtpart0,jtpart1
  real(lk) qpart,apart
  real(lk),dimension(:) :: wzpart
  real(lk),dimension(0:) :: kapan,kapat,meshn,mesht,pflux,rdtem
  real(lk),dimension(mpdiag) :: diagpart
  real(lk),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
  real(lk),dimension(0:,:) :: data1d
  real(lk),dimension(:,:),optional :: zpart1
  real(lk),dimension(0:,:),optional :: phit,dnt

!declaration of the local variables
  integer i,m,ii,ierror,igyro,ij,irsmooth,icount,iir(mpmax),isp,jst,ksz,bind,nb
  real(lk) rdum,tdum,dtime,q,q_inv,rg,b,g,gp,ri,rip,dbdp,dbdt,dbdz,dedb,deni,&
    gqi,deltaf,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,&
    wdot,wdrift,delr,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,tp_inv,epara,wpara,&
    wdrive,perturb,psitmp,thetatmp,zetatmp,cmratio,cinv,dapdp,dapdt,dapdz,dpx,dp2,dtx,&
    dt2,dzx,dx(spdim),lam,dlamdp,dlamdt,dlamdz,plampt,&
    kappan,kappat,b_inv,phipEq,dzonal,&
    angmom,fullf,dp1,upara0,dupara0dp,kappav,dedb0,energy0,majorr,dmajorrdp,&
    gyrodum,paxis,bdum,vtot,dist_term,eperp,&
    dtem(0:mpsi),dden(0:mpsi),ddum(0:mpsi),vdrtmp(0:mpsi),dmark(0:mpsi),&
    diagtmp(mpdiag),data1dtmp(0:mpsi,mpdata1d),delp(mpsi),deltab(0:1),&
    pdum,cost,sint,xdot,ydot,xdum,ydum,&
! analytic island terms
    wisland,lami,dlamidp,dlamidt,dlamidz,i1,i2,i3,i4,wigc(4,mpmax),&
! equilibrium island terms
! equilibrium island spline delta in zeta temporarily different to the
! geometry spline delta zeta. And the delta array is always 3D,thus 27
! elements.
    dzxisland,dxisland(27),lami_eq,dlamidp_eq,dlamidt_eq,dlamidz_eq,vap,&
    

    dphieffdp,dphieffdt,dphieffdz,kappas,dupara0dt,kappavt,wdrift0,wpara0,&
    qp,gp2,rip2,vet,wdriftind,delz,np_inv,dphiinddp,dphiinddt,dphiinddz,&
    pptpt,dpsidp,dpsidt,dpsidz,e1,e2,e3,e4,b1,b2,b3,b4,b5,b6,b7,b8,&
    b9,b10,delt(0:mpsi),wpgc(4,mpmax),wbgc(10,mpmax)
# 220



!$acc declare create(iir)
!$acc declare create(delp,wpgc,vdrtmp,wbgc)
!$acc declare create(dtem,dmark,dden)
!$acc declare create(wigc)

  logical subcycle
# 231

  real(lk) kperprhoi,zetagyravg

  subcycle=.false.
  if(present(ncyclep))then
    if(ncyclep>0) subcycle=.true.
  endif

!paxis=0.5*(8.0*rho0)**2
  paxis=0.0
  if(psi0<1.0e-8)paxis=psimesh(1)
  delr=1.0/deltar

  delt=1.0/deltat !!no 2pi/deltar in version3
  delz=1.0/deltaz

  cmratio=qpart/apart
  cinv=1.0/qpart
  perturb=real(nonlinear)
  vdrtmp=0.0
  delp=1.0/deltap


!$acc update device(delp)


! for out-of-boundary particles
  deltab(0)=(maxbfield(0)-minbfield(0))/real(mthetamax/2 -1)
  deltab(1)=(maxbfield(1)-minbfield(1))/real(mthetamax/2 -1)


  if(.not. subcycle)then
    if(irk==1)then
! 1st step of Runge-Kutta method
      dtime=0.5*tstep
# 268

!$omp parallel do private(m)

      do m=1,mp
        zpart0(1:nparam,m)=zpart(1:nparam,m)
      enddo
!$acc end parallel

! 2nd step of Runge-Kutta method
    else
      dtime=tstep
      if(nonlinear==1)vdrtmp=pflux
    endif
  else
    if(irksub==1)then
! 1st step of Runge-Kutta method
      dtime=0.25*tstep/real(ncyclep)
      if(icycle==1)then
        if(.not. (irk/=1 .or. present(ihybrid) .and. ihybrid/=1))then
          mp1=mp

# 291

!$omp parallel do private(m)

          do m=1,mp
            zpart1(1:nparam,m)=zpart(1:nparam,m)
          enddo
!$acc end parallel
        else
          mp=mp1
# 302

!$omp parallel do private(m)

          do m=1,mp
            zpart(1:nparam,m)=zpart1(1:nparam,m)
          enddo
!$acc end parallel
        endif
      endif
# 313

!$omp parallel do private(m)

      do m=1,mp
        zpart0(1:nparam,m)=zpart(1:nparam,m)
      enddo
!$acc end parallel

! 2nd step of Runge-Kutta method
    else
      dtime=0.5*tstep/real(ncyclep)
      if(.not. (nonlinear/=1 .or. irk/=2 .or. present(ihybrid) .and.&
        ihybrid/=nhybrid))vdrtmp=pflux
    endif
  endif


!$acc update device(vdrtmp)


! gather e_field using ngyro-point gyro-averaging
  gyrodum=1.0/real(ngyro)

! island modification of lambda and lambda gradient, represented by wigc(1,m) and wigc(2,m) wigc(3,m) wigc(4,m)
  if(island==1)then
# 340

!$omp parallel do private(m,igyro,i1,i2,i3,i4,wz1,wz0,wp1,wp0,wt00,wt10,wt01,wt11,ij)

    do m=1,mp
      i1=0.0
      i2=0.0
      i3=0.0
      i4=0.0
      wz1=wzpart(m)                !weight for upper toroidal grid
      wz0=1.0-wz1                 !weight for lower toroidal grid

      do igyro=1,ngyro
        wp1=wppart(igyro,m)       !outer flux surface
        wp0=1.0-wp1              !inner flux surface

        wt10=wp0*wtpart0(igyro,m) !upper poloidal grid on inner flux surface
        wt00=wp0-wt10            !lower poloidal grid on inner flux surface

        wt11=wp1*wtpart1(igyro,m) !upper poloidal grid on outer flux surface
        wt01=wp1-wt11            !lower poloidal grid on outer flux surface

        ij=jtpart0(igyro,m)       !lower poloidal grid on inner flux surface
        i1=i1+wt00*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt00*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt00*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt00*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))


        ij=ij+1                  !upper poloidal grid on inner flux surface
        i1=i1+wt10*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt10*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt10*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt10*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))

        ij=jtpart1(igyro,m)       !lower poloidal grid on outer flux surface
        i1=i1+wt01*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt01*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt01*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt01*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))

        ij=ij+1                  !upper poloidal grid on outer flux surface
        i1=i1+wt11*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt11*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt11*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt11*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))
      enddo

      wigc(1,m)=gyrodum*i1
      wigc(2,m)=gyrodum*i2
      wigc(3,m)=gyrodum*i3
      wigc(4,m)=gyrodum*i4
    enddo
!$acc end parallel
  endif

! electrostatic fluctuations
  if(magnetic==0)then
# 399


!$omp parallel do private(m,igyro,e1,e2,e3,wz1,wz0,wp0,wp1,wt00,wt10,wt01,wt11,ij,e4)
# 404


    do m=1,mp
      e1=0.0
      e2=0.0
      e3=0.0

      e4=0.0


      wz1=wzpart(m)                !weight for upper toroidal grid
      wz0=1.0-wz1                 !weight for lower toroidal grid

      do igyro=1,ngyro
        wp1=wppart(igyro,m)       !outer flux surface
        wp0=1.0-wp1              !inner flux surface

        wt10=wp0*wtpart0(igyro,m) !upper poloidal grid on inner flux surface
        wt00=wp0-wt10            !lower poloidal grid on inner flux surface

        wt11=wp1*wtpart1(igyro,m) !upper poloidal grid on outer flux surface
        wt01=wp1-wt11            !lower poloidal grid on outer flux surface

        ij=jtpart0(igyro,m)       !lower poloidal grid on inner flux surface
        e1=e1+wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt00*(wz0*phit(0,ij)+wz1*phit(1,ij))


        ij=ij+1                  !upper poloidal grid on inner flux surface
        e1=e1+wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt10*(wz0*phit(0,ij)+wz1*phit(1,ij))


        ij=jtpart1(igyro,m)       !lower poloidal grid on outer flux surface
        e1=e1+wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt01*(wz0*phit(0,ij)+wz1*phit(1,ij))


        ij=ij+1                  !upper poloidal grid on outer flux surface
        e1=e1+wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt11*(wz0*phit(0,ij)+wz1*phit(1,ij))

      enddo

      wpgc(1,m)=gyrodum*e1
      wpgc(2,m)=gyrodum*e2
      wpgc(3,m)=gyrodum*e3

      wpgc(4,m)=gyrodum*e4


      wbgc(:,m)=0.0
    enddo
!$acc end parallel
  else
! electromagnetic fields
# 474


!$omp parallel do private(m,igyro,e1,e2,e3,b1,b2,b3,b4,b5,wz1,wz0,wp0,wp1,wt00,wt10,wt01,wt11,ij)&
!$omp& private(e4,b6,b7,b8,b9,b10)
# 480


    do m=1,mp
!$acc cache(readonly:gradphi,dnt,gradapara,gradpsi,sapara)
      e1=0.0
      e2=0.0
      e3=0.0

      e4=0.0

      b1=0.0
      b2=0.0
      b3=0.0
      b4=0.0
      b5=0.0

      b6=0.0
      b7=0.0
      b8=0.0
      b9=0.0
      b10=0.0


      wz1=wzpart(m)                !weight for upper toroidal grid
      wz0=1.0-wz1                 !weight for lower toroidal grid
      do igyro=1,ngyro
        wp1=wppart(igyro,m)       !outer flux surface
        wp0=1.0-wp1              !inner flux surface

        wt10=wp0*wtpart0(igyro,m) !upper poloidal grid on inner flux surface
        wt00=wp0-wt10            !lower poloidal grid on inner flux surface

        wt11=wp1*wtpart1(igyro,m) !upper poloidal grid on outer flux surface
        wt01=wp1-wt11            !lower poloidal grid on outer flux surface

        ij=jtpart0(igyro,m)       !lower poloidal grid on inner flux surface
        e1=e1+wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt00*(wz0*dnt(0,ij)+wz1*dnt(1,ij))

        b1=b1+wt00*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
        b2=b2+wt00*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
        b3=b3+wt00*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
        b4=b4+wt00*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
        b5=b5+wt00*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

        b6=b6+wt00*(wz0*gradphieff(1,0,ij)+wz1*gradphieff(1,1,ij))
        b7=b7+wt00*(wz0*gradphieff(2,0,ij)+wz1*gradphieff(2,1,ij))
        b8=b8+wt00*(wz0*gradpsi(1,0,ij)+wz1*gradpsi(1,1,ij))
        b9=b9+wt00*(wz0*gradpsi(2,0,ij)+wz1*gradpsi(2,1,ij))
        b10=b10+wt00*(wz0*gradpsi(3,0,ij)+wz1*gradpsi(3,1,ij))


        ij=ij+1                  !upper poloidal grid on inner flux surface
        e1=e1+wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt10*(wz0*dnt(0,ij)+wz1*dnt(1,ij))

        b1=b1+wt10*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
        b2=b2+wt10*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
        b3=b3+wt10*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
        b4=b4+wt10*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
        b5=b5+wt10*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

        b6=b6+wt10*(wz0*gradphieff(1,0,ij)+wz1*gradphieff(1,1,ij))
        b7=b7+wt10*(wz0*gradphieff(2,0,ij)+wz1*gradphieff(2,1,ij))
        b8=b8+wt10*(wz0*gradpsi(1,0,ij)+wz1*gradpsi(1,1,ij))
        b9=b9+wt10*(wz0*gradpsi(2,0,ij)+wz1*gradpsi(2,1,ij))
        b10=b10+wt10*(wz0*gradpsi(3,0,ij)+wz1*gradpsi(3,1,ij))


        ij=jtpart1(igyro,m)       !lower poloidal grid on outer flux surface
        e1=e1+wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt01*(wz0*dnt(0,ij)+wz1*dnt(1,ij))

        b1=b1+wt01*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
        b2=b2+wt01*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
        b3=b3+wt01*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
        b4=b4+wt01*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
        b5=b5+wt01*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

        b6=b6+wt01*(wz0*gradphieff(1,0,ij)+wz1*gradphieff(1,1,ij))
        b7=b7+wt01*(wz0*gradphieff(2,0,ij)+wz1*gradphieff(2,1,ij))
        b8=b8+wt01*(wz0*gradpsi(1,0,ij)+wz1*gradpsi(1,1,ij))
        b9=b9+wt01*(wz0*gradpsi(2,0,ij)+wz1*gradpsi(2,1,ij))
        b10=b10+wt01*(wz0*gradpsi(3,0,ij)+wz1*gradpsi(3,1,ij))


        ij=ij+1                  !upper poloidal grid on outer flux surface
        e1=e1+wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
        e2=e2+wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
        e3=e3+wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

        e4=e4+wt11*(wz0*dnt(0,ij)+wz1*dnt(1,ij))

        b1=b1+wt11*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
        b2=b2+wt11*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
        b3=b3+wt11*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
        b4=b4+wt11*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
        b5=b5+wt11*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

        b6=b6+wt11*(wz0*gradphieff(1,0,ij)+wz1*gradphieff(1,1,ij))
        b7=b7+wt11*(wz0*gradphieff(2,0,ij)+wz1*gradphieff(2,1,ij))
        b8=b8+wt11*(wz0*gradpsi(1,0,ij)+wz1*gradpsi(1,1,ij))
        b9=b9+wt11*(wz0*gradpsi(2,0,ij)+wz1*gradpsi(2,1,ij))
        b10=b10+wt11*(wz0*gradpsi(3,0,ij)+wz1*gradpsi(3,1,ij))

      enddo

      wpgc(1,m)=gyrodum*e1
      wpgc(2,m)=gyrodum*e2
      wpgc(3,m)=gyrodum*e3

      wpgc(4,m)=gyrodum*e4

      wbgc(1,m)=gyrodum*b1
      wbgc(2,m)=gyrodum*b2
      wbgc(3,m)=gyrodum*b3
      wbgc(4,m)=gyrodum*b4
      wbgc(5,m)=gyrodum*b5

      wbgc(6,m)=gyrodum*b6
      wbgc(7,m)=gyrodum*b7
      wbgc(8,m)=gyrodum*b8
      wbgc(9,m)=gyrodum*b9
      wbgc(10,m)=gyrodum*b10

    enddo
!$acc end parallel
  endif

! update GC position. Assuming psi0>spdpsi; will be relaxed later
# 621

!$omp parallel do private(m,psitmp,thetatmp,isp,dpx,dp1,dp2,jst,dtx,dt2,dzx,dx,&
!$omp& dxisland,rg,q,g,&
!$omp& gp,ri,rip,b,dbdp,dbdt,dbdz,ii,dedb,deni,gqi,upara,energy,dptdp,dptdt,&
!$omp& dptdz,dapdp,dapdt,dapdz,dlamidt,dlamidz,dlamidp_eq,dlamidt_eq,dlamidz_eq,&
!$omp& epara,vdr,wdrive,wpara,wdrift,wdot,rdot,&
!$omp& pdot,tdot,zdot,kappa,lam,dlamdp,dlamdt,dlamdz,plampt,kappan,&
!$omp& kappat,tp_inv,upara0,dupara0dp,kappav,b_inv,majorr,dmajorrdp,dedb0,&
!$omp& energy0,phipEq,kperprhoi,zetagyravg,vtot,dist_term,eperp,dzonal)&

!$omp& private(kappavt,vet,wdriftind,gp2,rip2,pptpt,wp0,wp1,dupara0dt,&
!$omp& wdrift0,wpara0,dphieffdp,dphieffdt,dphieffdz,qp,kappas,&
!$omp& q_inv,np_inv,dphiinddp,dphiinddt,dphiinddz,dpsidp,dpsidt,dpsidz)
# 636


  do m=1,mp
    psitmp=zpart(1,m)
    isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
! For 1D profiles, spline is regular near the axis
    dpx=psitmp-spdpsi*real(isp-1)
    dp2=dpx*dpx
! poloidal spline is regular
    thetatmp=zpart(2,m)
    jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
    dtx=thetatmp-spdtheta*real(jst-1)
    dt2=dtx*dtx
! zeta spline for 3D equilibrium is temporarily in its own format
!> @todo
!> update the 3D VMEC equilibrium to use the same subsectioned
!> zeta spline scheme, as in M3DC1 alpha spline.
    zetatmp = zpart(3,m)
    dzx=zetatmp-zeta0

    dx(1)=1.0
    dx(2)=dpx
    dx(3)=dp2

! 1D spline in psi
    q=    qpsi(1,isp)   +qpsi(2,isp)*dx(2)   +qpsi(3,isp)*dx(3)
    g=    gpsi(1,isp)   +gpsi(2,isp)*dx(2)   +gpsi(3,isp)*dx(3)
    gp=                  gpsi(2,isp)       +gpsi(3,isp)*dx(2)*2.0
    ri=   cpsi(1,isp)   +cpsi(2,isp)*dx(2)   +cpsi(3,isp)*dx(3)
    rip=                 cpsi(2,isp)       +cpsi(3,isp)*dx(2)*2.0
    upara0=ropp(1,isp)  +ropp(2,isp)*dx(2)   +ropp(3,isp)*dx(3)
    dupara0dp=           ropp(2,isp)       +ropp(3,isp)*dx(2)*2.0
!    er=   erpp(1,isp)   +erpp(2,isp)*dx(2)   +erpp(3,isp)*dx(3)

! To make the 2nd derivatives of I and g continuous,
! currently we use linear interpolations as a temporary solution.
! In the future, higher order spline functions of I and g need to be used
    gp2=2.0*(dpx*spdpsi_inv*(gpsi(3,isp+1)-gpsi(3,isp))+gpsi(3,isp))
    rip2=2.0*(dpx*spdpsi_inv*(cpsi(3,isp+1)-cpsi(3,isp))+cpsi(3,isp))


! For coordinates splines, first cell needs to be specially treated.
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
    if(isp==1)dpx=sqrt(dpx)
    dp2=dpx*dpx

    dx(1)=1.0
    dx(2)=dpx
    dx(3)=dp2
    dx(4:6)=dx(1:3)*dtx
    dx(7:9)=dx(1:3)*dt2
# 690

! Evaluate rg with special dpsi.
    rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dx(2)  +rgpsi(3,isp)*dx(3)
! Evaluate b, majorr, dbdp, dbdt, dbdz
! spline in (psi, theta, zeta)
    b=0.0
    majorr=0.0
    do ii = 1, spdim
      b=b +bsp(ii,isp,jst)*dx(ii)
      majorr=majorr +xsp(ii,isp,jst)*dx(ii)
    enddo
! zeta-derivative
    dbdz=0.0
# 709

! theta-derivative
    dx(4:6)=dx(1:3)*1.0
    dx(7:9)=dx(1:3)*dtx*2.0
! L.Shi 03/14/2017:
! It is essential to set dx(1:3)=0 to correctly calculate
! theta derivatives in 3D case.
    dx(1:3)=0.0
# 720

    dbdt=0.0
    do ii = 4, spdim
      dbdt=dbdt +bsp(ii,isp,jst)*dx(ii)
    enddo

! psi-derivative
    dx(1)=0.0
    dx(2)=1.0
    if(isp==1)dx(2)=0.5/dpx
    dx(3)=dpx*2.0
    if(isp==1)dx(3)=1.0
    dx(4:6)=dx(1:3)*dtx
    dx(7:9)=dx(1:3)*dt2
# 737

    dbdp=0.0
    dmajorrdp=0.0
    do ii = 2, spdim
      dbdp=dbdp +bsp(ii,isp,jst)*dx(ii)
      dmajorrdp=dmajorrdp +xsp(ii,isp,jst)*dx(ii)
    enddo

    b_inv=1.0/b

    zetagyravg=1.0

! Equilibrium island quantities. These calculations are hard coded inline to accelerate the execution.
! Created by L. Shi, 03/22/17
! Equilibrium island terms are by default set to 0
    lami_eq = 0.0_lk
    dlamidp_eq = 0.0_lk
    dlamidt_eq = 0.0_lk
    dlamidz_eq = 0.0_lk
! We calculate equilibrium island when m3dc1 islands are included
    if(numereq==3 .and. ndim>1) then
! equilibrium island is expressed as a 3D spline
! We first calculate the spline indexes. The psi and theta indexesare the same
! as calculated from other equilibrium splines.
      dxisland(1)=1
      dxisland(2)=dpx
      dxisland(3)=dp2
      dxisland(4:6)=dxisland(1:3)*dtx
      dxisland(7:9)=dxisland(1:3)*dt2
! dzeta needs to be calculated again because of sub-section used in equilibrium
! island spline.
      ksz=max(1,min(nzsp_sec,ceiling((zetatmp-zeta0)*spdzeta_inv)))
      dzxisland = zetatmp - (zeta0 + spdzeta*(ksz-1))
      dxisland(10:18)=dxisland(1:9)*dzxisland
      dxisland(19:27)=dxisland(10:18)*dzxisland
! Now we evaluate alpha_eq using spline coefficients
      do ii=1,27
        lami_eq = lami_eq+alphasp(ii,isp,jst,ksz)*dxisland(ii)
      enddo
! psi derivative needs special care when isp==1. But we assume here that the
! perturbed magnetic field is zero at the axis, and ignore this part.
!> @todo
!> isp==1 case needs to be taken care correctly when alpha_eq is not sufficiently
!> small near the axis.
      dxisland(1)=0
      dxisland(2)=1
      dxisland(3)=2*dpx
      dxisland(4:6)=dxisland(1:3)*dtx
      dxisland(7:9)=dxisland(1:3)*dt2
      dxisland(10:18)=dxisland(1:9)*dzxisland
      dxisland(19:27)=dxisland(10:18)*dzxisland
      do ii=2,27
        dlamidp_eq = dlamidp_eq+alphasp(ii,isp,jst,ksz)*dxisland(ii)
      enddo
! theta derivative can be easily evaluated
      dxisland(1:3)=0.0_lk
      dxisland(4)=1
      dxisland(5)=dpx
      dxisland(6)=dp2
      dxisland(7:9)=dxisland(4:6)*dtx*2
      dxisland(10:18)=dxisland(1:9)*dzxisland
      dxisland(19:27)=dxisland(10:18)*dzxisland
      do ii=4,27
        dlamidt_eq = dlamidt_eq+alphasp(ii,isp,jst,ksz)*dxisland(ii)
      enddo
! zeta derivative can be easily evaluated
      dxisland(1:9)=0.0_lk
      dxisland(10)=1
      dxisland(11)=dpx
      dxisland(12)=dp2
      dxisland(13:15)=dxisland(10:12)*dtx
      dxisland(16:18)=dxisland(13:15)*dtx
      dxisland(19:27)=dxisland(10:18)*dzxisland*2
      do ii=10,27
        dlamidz_eq = dlamidz_eq+alphasp(ii,isp,jst,ksz)*dxisland(ii)
      enddo
    endif
# 827

    q_inv=1.0/q


! perturbed electric field
    dptdp=wpgc(1,m)
    dptdt=wpgc(2,m)
    dptdz=wpgc(3,m)-wpgc(2,m)*q_inv

    pptpt=wpgc(4,m)

! effective potential
    dphieffdp=wbgc(6,m)
    dphieffdt=wbgc(7,m)
    dphieffdz=wbgc(4,m)-wbgc(7,m)*q_inv

! induced potential
    dphiinddp=dphieffdp-dptdp
    dphiinddt=dphieffdt-dptdt
    dphiinddz=dphieffdz-dptdz

! perturbed magnetic flux
    dpsidp=wbgc(8,m)
    dpsidt=wbgc(9,m)
    dpsidz=wbgc(10,m)-wbgc(9,m)*q_inv


! perturbed vector potential
    dapdp=wbgc(1,m)
    dapdt=wbgc(2,m)
    dapdz=wbgc(3,m)-wbgc(2,m)*q_inv

! analytic island perturbed magnetic field
    lami=0.0_lk
    dlamidp=0.0_lk
    dlamidt=0.0_lk
    dlamidz=0.0_lk
    if(island==1)then
      lami=wigc(4,m)
      dlamidp=wigc(1,m)
      dlamidt=wigc(2,m)
! Kaisheng & Lei: Add dalpha/dtheta correction to evaluation of
! dalpha/dzeta, since gradient evaluates parallel derivative instead of
! zeta derivative.
      dlamidz=wigc(3,m)-wigc(2,m)*q_inv
    endif

! analytic and equilibrium island fields are summed up
    lami = lami + lami_eq
    dlamidp = dlamidp + dlamidp_eq
    dlamidt = dlamidt + dlamidt_eq
    dlamidz = dlamidz + dlamidz_eq

! self-consistent lambda term is calculated from delta_Apara
! lambda=apara*b_inv
    lam=wbgc(5,m)*b_inv
    dlamdp=(dapdp-lam*dbdp)*b_inv
    dlamdt=(dapdt-lam*dbdt)*b_inv
    dlamdz=(dapdz-lam*dbdz)*b_inv

! we use lami_eq and other _eq variables to temporarily store the
! self-consistent lambda terms. For later, when determining the non-linear
! run conditions, these terms are still needed.
    lami_eq = lam
    dlamidp_eq = dlamdp
    dlamidt_eq = dlamdt
    dlamidz_eq = dlamdz

! Now, we sum all the lambda terms for weight calculation
! analytic and equilibrium islands contributions are added to lambda
! the island terms are properly set to zero if turned off, so no need
! for conditional control
    lam = lam + lami + lami_eq
    dlamdp = dlamdp + dlamidp + dlamidp_eq
    dlamdt = dlamdt + dlamidt + dlamidt_eq
    dlamdz = dlamdz + dlamidz + dlamidz_eq

    ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
    dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface

    phipEq=-(dp1*mesher(ii-1)+(1.0-dp1)*mesher(ii)) ! mesher = -dphi/dpsi
    kappan=dp1*kapan(ii-1)+(1.0-dp1)*kapan(ii)
    kappat=dp1*kapat(ii-1)+(1.0-dp1)*kapat(ii)
    tp_inv=1.0/(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))

    np_inv=1.0/(dp1*meshn(ii-1)+(1.0-dp1)*meshn(ii))

    dzonal=dp1*phip00(ii-1)+(1.0-dp1)*phip00(ii)
    deni=1.0/(g*q + ri + (zpart(4,m)+lam)*(g*rip-ri*gp))
    gqi=1.0/(g*q+ri)

    upara0=upara0*majorr
    dupara0dp=dupara0dp*majorr+upara0*dmajorrdp


! kappas=(partial S/partial psi/S), S defined in equilibrium current paper
    kappas=(g*rip2-ri*gp2)/(g*rip-ri*gp)-(gp*q+g*qp+rip)*gqi

    dupara0dt=0.0
    if(eqcurrent==1)then
      upara0=upara0-(2.0_lk*rho0*rho0/betae)*b*np_inv*(g*rip-ri*gp)*gqi
      dupara0dp=dupara0dp-(2.0_lk*rho0*rho0/betae)*np_inv*(g*rip-ri*gp)*gqi*(dbdp+b*(kappan+kappas))
      dupara0dt=-(2.0_lk*rho0*rho0/betae)*dbdt*np_inv*(g*rip-ri*gp)*gqi
    endif

    upara=zpart(4,m)*b*cmratio
    dedb=zpart(4,m)*zpart(4,m)*b*cmratio+cinv*zpart(6,m)*zpart(6,m)

    dedb0=-cinv*upara*upara0*apart*b_inv
# 938

    energy=0.5*apart*upara*upara+zpart(6,m)*zpart(6,m)*b
    energy0=0.5*apart*(upara-upara0)*(upara-upara0)+zpart(6,m)*zpart(6,m)*b

! Additional term for weight equation with slowing distribution
    if(pload==11)then !slowing down
      vtot=min(sqrt(2.0_lk*energy/apart),sd_v0)
      eperp=2*zpart(6,m)*zpart(6,m)/(apart*vtot**2)
      dist_term=-2*qpart*eperp**2*sd_widthInv**2/(zpart(6,m)*zpart(6,m))*(eperp-sd_l0)+3.0*qpart*vtot/(vtot**3+sd_vc**3)/apart
    else ! local maxwellian
        dist_term=qpart*tp_inv
    endif

    kappa=kappan+kappat*(energy0*tp_inv-1.5)
    kappav=-(upara-upara0)*dupara0dp*apart*tp_inv

    kappavt=-(upara-upara0)*dupara0dt*apart*tp_inv


! dne/n0=e phi/t for electrostatic case
! parallel electric field epara, & dlambda/dt=-epara-gradphi(3,m)
    if (magnetic==0) then
# 962

      epara=-wpgc(3,m)*q*b*deni

! static island doesn't have time derivative, so no extra term
! due to island here.
!> @todo
!> additional terms for non-static island.
      plampt=0.0
    else
      epara=-wbgc(4,m)*b*q*deni
      plampt=-(epara+wpgc(3,m)*q*b*deni)*b_inv
    endif

! ExB drift in radial direction for w-dot and flux diagnostics
    vdr=(ri*dptdz-g*dptdt)*deni
! vap: the magnetic field line perturbation induced radial drift velocity.
! Only used for diagnosis. Effect of vap in weight equation and particle
! trajectory is included in wisland and lam terms in the calculation of
! wdot and pdot.
    vap=-b*(ri*dlamdz-g*dlamdt)*upara*deni ! upara0 effect is kept

    vet=g*dptdp*deni
    wdrive=vdr*(kappa+kappav)+vet*kappavt

    if (magnetic==0) then
      wpara=pptpt*qpart*tp_inv
      wdrift=(g*dedb*dbdt - ri*dedb*dbdz)*deni*dzonal*qpart*tp_inv+(g*dptdt-ri*dptdz)*deni*(dzonal+phipEq)*qpart*tp_inv
      wdriftind=0.0
    else
      wpara=-pptpt*np_inv+(kappan-kappa)*dphiinddt*q_inv
      wdriftind=-(g*dbdt*dphiinddp-g*dbdp*dphiinddt+ri*dbdp*dphiinddz-ri*dbdz*dphiinddp)*deni*dedb*tp_inv*qpart-&
      (g*dbdt*dpsidp-ri*dbdz*dpsidp-g*dbdp*dpsidt+ri*dbdp*dpsidz)*deni*dedb*(kappa-qpart*tp_inv*phipEq)+&
      (g*dbdt-ri*dbdz)*dzonal*deni*dedb
      if(eqcurrent==1)then
        wdriftind=wdriftind+apart*upara*upara*tp_inv*deni*(gp*dphiinddt-rip*dphiinddz)+&
          +apart*upara*upara*deni*kappa*cinv*(gp*dpsidt-rip*dpsidz)
      endif
      wdrift=(g*(dphieffdt*tp_inv*qpart+kappa*dpsidt-qpart*phipEq*tp_inv*dpsidt)-ri*(dphieffdz*tp_inv*qpart+kappa*dpsidz-qpart*phipEq*tp_inv*dpsidz))*deni*(dzonal+phipEq)
    endif

    wdrift0=(g*dbdp*dptdt-g*dbdt*dptdp-ri*dbdp*dptdz+ri*dbdz*dptdp)*gqi*tp_inv*apart*upara*upara0*b_inv
    if(eqcurrent==1)then
      wdrift0=wdrift0-apart*upara*upara0*tp_inv*deni*(gp*dptdt-rip*dptdz)
    endif
    wpara0=-upara0*qpart*tp_inv*epara

    wdot=(zpart(7,m)-paranl*zpart(5,m))*(wdrive+wpara+wdrift+wdriftind+wdrift0+wpara0)

# 1036


! Particle orbit calculation
! only self-consistent fields are controlled by the non-linear flag
    dptdp=dptdp*perturb + phipEq
    dptdt=dptdt*perturb
    dptdz=dptdz*perturb
! lami_eq and other _eq variables are storing the self-consistent lambda
! terms
    dlamdp=dlamidp_eq*perturb + dlamidp
    dlamdt=dlamidt_eq*perturb + dlamidt
    dlamdz=dlamidz_eq*perturb + dlamidz
    lam=lami_eq*perturb + lami
    plampt=plampt*perturb

! particle velocity
! rho_para dot partially checked by L.Shi, 03/22/2017
! dB/dzeta term due to non-axisymmetric equilibrium not checked.
! paranl is added to necessary terms
    rdot = ((gp*(zpart(4,m)+lam)+g*dlamdp-1.0)*(dedb*dbdt+paranl*dptdt)-&
      (q+rip*(zpart(4,m)+lam)+ri*dlamdp)*(dedb*dbdz+paranl*dptdz)+&
      (ri*dlamdz-g*dlamdt)*(dedb*dbdp+paranl*(dptdp-phipEq)+phipEq))*deni-paranl*plampt
! psi dot partially checked by L.Shi, 03/22/2017.
! dB/dzeta term due to non-axisymmetric equilibrium not checked.
! vdrtmp term should be turned off by default
!pdot = (-g*dedb*dbdt +ri*dedb*dbdz-g*dptdt+ri*dptdz+upara*b*(g*dlamdt-ri*dlamdz))*deni-vdrtmp(ii)
    pdot = (-g*dedb*dbdt +ri*dedb*dbdz-g*dptdt+ri*dptdz+upara*b*(g*dlamdt-ri*dlamdz))*deni
! theta dot checked by L.Shi, 03/22/2017
    tdot = (upara*b*(1.0-gp*(zpart(4,m)+lam)-g*dlamdp)+g*(dedb*dbdp+dptdp))*deni
! zeta dot checked by L.Shi, 03/22/2017
    zdot = (upara*b*(q+rip*(zpart(4,m)+lam)+ri*dlamdp)-ri*(dedb*dbdp+dptdp))*deni

! update particle position
    if(zpart0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
      pdum=sqrt(zpart0(1,m))
      cost=pdum*cos(zpart0(2,m))
      sint=pdum*sin(zpart0(2,m))
      pdum=1.0/zpart(1,m)
      xdot   = 0.5*pdot*cost*pdum-sint*tdot
      ydot   = 0.5*pdot*sint*pdum+cost*tdot
      xdum   = cost + dtime*xdot
      ydum   = sint + dtime*ydot
      zpart(1,m) = max(1.0e-8_lk*psi1,xdum*xdum+ydum*ydum)
      zpart(2,m) = sign(1.0_lk,ydum)*acos(max(-1.0_lk,min(1.0_lk,xdum/sqrt(zpart(1,m)))))
    else
      zpart(1,m) = max(1.0e-8_lk*psi1,zpart0(1,m)+dtime*pdot)
      zpart(2,m) = zpart0(2,m)+dtime*tdot
    endif

    zpart(3,m) = zpart0(3,m)+dtime*zdot
    zpart(4,m) = zpart0(4,m)+dtime*rdot
    zpart(5,m) = zpart0(5,m)+dtime*wdot

    zpart(2,m)=modulo(zpart(2,m),pi2)
    zpart(3,m)=modulo(zpart(3,m),torbound)

! store GC information for flux measurements
    wpgc(1,m)=vdr
    wpgc(2,m)=energy
    wpgc(3,m)=upara*majorr
! store magnetic flutter velocity in island case. L.Shi,04/17/2017
! use wbgc(1,m) as temporary storage since wbgc array is not used after
! this line
    wbgc(1,m)=vap
  enddo
!$acc end parallel

  if(pload>99) then
!$acc parallel loop
     do m=1,mp
        zpart(5,m)=1.0 ! for full-f simulation
     enddo
!$acc end parallel
  endif

  if(subcycle .and. irksub==2 .or. (.not. subcycle) .and. irk==2)then

! out of boundary particle
    nb=mthetamax/2 !B array index limit, see bminmax in setup.F90
# 1118

!$omp parallel do private(m,tdum,isp,jst,dx,dtx,bdum,bind,ii)

    do m=1,mp
      if(zpart(1,m) > psi1)then
        zpart(1,m)=zpart0(1,m)
        zpart(3,m)=zpart0(3,m)
        zpart(4,m)=zpart0(4,m)
        zpart(5,m)=zpart0(5,m)
        if(numereq==0)then
          zpart(2,m)=2.0*pi-zpart0(2,m) ! analytic equilibrium assumes up-down symmetry
        else
          tdum=zpart0(2,m)
          isp=max(1,min(lsp-1,ceiling(psi1*spdpsi_inv)))
          jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
          dx(1)=1.0
          dx(2)=psi1-spdpsi*real(isp-1)
          dx(3)=dx(2)*dx(2)
          dtx=tdum-spdtheta*real(jst-1)
          dx(4:6)=dx(1:3)*dtx
          dx(7:9)=dx(4:6)*dtx
! b-field at leaving point
          bdum=0.0
          do ii = 1, 9
             bdum=bdum +bsp(ii,isp,jst)*dx(ii)
          enddo
          bind=floor((bdum-minbfield(1))/deltab(1))+1 !index of b-field bin
! bdum could be smaller than minbfield or larger than maxbfield due
! to low resolution in theta grids
! we need to keep ``bind`` in the correct range
! Note that the max index in thetaupp and thetadwn is mthetamax/2
          bind = min(nb, max(1,bind))
          tdum=mod(pi2+tdum-thetabmin(1),pi2) ! poloidal angle w.r.t. min b-field angle
          if(tdum<thetabmax(1))then ! return in the  opposite half-plane
             zpart(2,m)=mod(deltat(mpsi)*(real(thetaupp(1,bind))+0.5),pi2)
          else
             zpart(2,m)=mod(deltat(mpsi)*(real(thetadwn(1,bind))+0.5),pi2)
          endif
        endif
      elseif(zpart(1,m) < psi0)then
        zpart(1,m)=zpart0(1,m)
        zpart(3,m)=zpart0(3,m)
        zpart(4,m)=zpart0(4,m)
        zpart(5,m)=zpart0(5,m)
        if(numereq==0)then
          zpart(2,m)=2.0*pi-zpart0(2,m) ! analytic equilibrium assumes up-down symmetry
        else
          tdum=zpart0(2,m)
          isp=max(1,min(lsp-1,ceiling(psi0*spdpsi_inv)))
          jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
          dx(1)=1.0
          dx(2)=psi0-spdpsi*real(isp-1)
          if(isp==1)dx(2)=sqrt(dx(2))
          dx(3)=dx(2)*dx(2)
          dtx=tdum-spdtheta*real(jst-1)
          dx(4:6)=dx(1:3)*dtx
          dx(7:9)=dx(4:6)*dtx
! b-field at leaving point
          bdum=0.0
          do ii = 1, 9
            bdum=bdum +bsp(ii,isp,jst)*dx(ii)
          enddo
          bind=floor((bdum-minbfield(0))/deltab(0))+1
! same as the outer boundary case, we need to bound ``bind`` in the
! correct range
          bind = min(nb, max(1, bind))
          tdum=mod(pi2+tdum-thetabmin(0),pi2)
          if(tdum<thetabmax(0))then ! return in the  opposite half-plane
             zpart(2,m)=mod(deltat(0)*(real(thetaupp(0,bind))+0.5),pi2)
          else
             zpart(2,m)=mod(deltat(0)*(real(thetadwn(0,bind))+0.5),pi2)
          endif
        endif
      endif
    enddo
!$acc end parallel


    if(irestore==1)then
! Restore temperature profile when running a nonlinear calculation

      if((nonlinear==1 .and. pload==1 .and. irk==2).and. (.not. subcycle .or. subcycle .and. &
        icycle==2*ncyclep) .and. ihybrid==nhybrid)then
# 1204


# 1208

!$omp parallel do private(m,psitmp,isp,dpx,rg)

        do m=1,mp
          psitmp=zpart(1,m)
          isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
          dpx=psitmp-spdpsi*real(isp-1)
          if(isp==1)dpx=sqrt(dpx)
          rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dpx*dpx
          iir(m)=max(0,min(mpsi,int((rg-rg0)*delr+0.5)))
        enddo
!$acc end parallel

!$acc parallel loop
        do ii=0,mpsi
          dtem(ii)=0.0
          dden(ii)=0.0
          dmark(ii)=0.0
        enddo
!$acc end parallel

!Todo: why cannot use reduction with openacc
# 1232

!$omp parallel do private(m,ii,fullf) reduction(+:dtem,dmark,dden)

        do m=1,mp
          ii=iir(m)
          fullf=zpart(7,m)
          deltaf=fullf*zpart(5,m)
!$acc atomic update
          dtem(ii)=dtem(ii)+wpgc(2,m)*deltaf
!$acc atomic update
          dmark(ii)=dmark(ii)+wpgc(1,m)*fullf
!$acc atomic update
          dden(ii)=dden(ii)+fullf
        enddo
!$acc end parallel


!$acc update host(dden,dmark,dtem)

        icount=mpsi+1
        ddum=0.0
        call MPI_ALLREDUCE(dtem,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dtem=ddum
        call MPI_ALLREDUCE(dmark,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dmark=ddum
        call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dden=ddum

        irsmooth=int(sqrt(real(mpsi)))
        dmark=dmark/max(1.0_lk,dden) !radial marker flux
        do i=1,irsmooth
          rdum=dmark(1)
          tdum=dmark(mpsi-1)
          dmark(1:mpsi-1)=0.5*dmark(1:mpsi-1)+0.25*(dmark(0:mpsi-2)+dmark(2:mpsi))
          dmark(0)=0.5*(dmark(0)+rdum)
          dmark(mpsi)=0.5*(dmark(mpsi)+tdum)
        enddo
        tdum=0.1
        pflux=(1.0-tdum)*pflux+tdum*dmark

! remove small scale temperature perturbation
        irsmooth=mpsi
        dtem=dtem/(mesht*max(1.0_lk,dden))
        do i=1,irsmooth
          rdum=dtem(1)
          tdum=dtem(mpsi-1)
          dtem(1:mpsi-1)=0.5*dtem(1:mpsi-1)+0.25*(dtem(0:mpsi-2)+dtem(2:mpsi))
          dtem(0)=0.5*(dtem(0)+rdum)
          dtem(mpsi)=0.5*(dtem(mpsi)+tdum)
        enddo
        tdum=0.01
        rdtem=(1.0-tdum)*rdtem+tdum*dtem


!$acc update device(rdtem)

# 1290

!$omp parallel do private(m,ii)

        do m=1,mp
          ii=iir(m)
          zpart(5,m)=zpart(5,m)-(wpgc(2,m)/mesht(ii)-1.5)*rdtem(ii)
        enddo
!$acc end parallel
      endif
    endif
  endif

  if(idiag==0 .and. (.not. subcycle .or. subcycle .and. icycle*irksub==1))then
! fluxes diagnose at irk=1
# 1306

!$omp parallel do private(m)

    do m=1,mpdiag
      diagpart(m)=0.0_lk
    enddo
!$acc end parallel

# 1316

!$omp parallel do private(m)

    do m=0,mpsi
      dden(m)=0.0_lk
    enddo
!$acc end parallel

# 1326

!$omp parallel do private(i,m)

    do i=1,mpdata1d
      do m=0,mpsi
        data1d(m,i)=0.0_lk
      enddo
    enddo
!$acc end parallel

!$acc parallel loop gang vector private(dx(spdim))
    do m=1,mp
      psitmp=zpart(1,m)
      isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
      dpx=psitmp-spdpsi*real(isp-1)
      if(isp==1)dpx=sqrt(dpx)
      dp2=dpx*dpx

! poloidal spline is regular
      thetatmp=zpart(2,m)
      jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
      dtx=thetatmp-spdtheta*real(jst-1)
      dt2=dtx*dtx
      dzx=zpart(3,m)-zeta0

!spline in psi, theta, and zeta
      dx(1)=1.0
      dx(2)=dpx
      dx(3)=dp2
      dx(4:6)=dx(1:3)*dtx
      dx(7:9)=dx(1:3)*dt2
# 1360

      majorr=0.0
!$acc loop seq
      do ii = 1, spdim
        majorr=majorr+xsp(ii,isp,jst)*dx(ii)
      enddo

      rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
      upara0=ropp(1,isp)  +ropp(2,isp)*dpx   +ropp(3,isp)*dp2
      upara0=upara0*majorr
!! r= rpsi(1,isp)  +rpsi(2,isp)*dpx  +rpsi(3,isp)*dp2
!! q=    qpsi(1,isp)   +qpsi(2,isp)*dpx   +qpsi(3,isp)*dp2

      ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
      dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface

# 1386

      tp_inv=1.0/(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii)+0.5*apart*upara0*upara0)

      fullf=zpart(7,m)
      deltaf=fullf*zpart0(5,m)
      energy=wpgc(2,m)*tp_inv-1.5
      vdr=wpgc(1,m)!!!*q/r
      angmom=wpgc(3,m)
      vap=wbgc(1,m) ! magnetic flutter velocity

! radial profile of particle and energy flux
!$acc atomic update
      dden(ii-1)=dden(ii-1)+fullf*dp1
!$acc atomic update
      dden(ii)  =dden(ii)+  fullf*(1.0-dp1)

! Use full radial velocity including EXB and magnetic fluttervelocities
! to calculate the total particle flux.
! L. Shi, 04/17/2017
!$acc atomic update
      data1d(ii-1,1)=data1d(ii-1,1)+(vdr+vap)*deltaf*dp1
!$acc atomic update
      data1d(ii,  1)=data1d(ii,  1)+(vdr+vap)*deltaf*(1.0-dp1)

! use full radial velocity to calculate the total energy and toroidal
! momentum fluxes. L. Shi, 05/02/2017
!$acc atomic update
      data1d(ii-1,2)=data1d(ii-1,2)+(vdr+vap)*deltaf*energy*dp1
!$acc atomic update
      data1d(ii,  2)=data1d(ii,  2)+(vdr+vap)*deltaf*energy*(1.0-dp1)
! radial profiles of momentum flux
!$acc atomic update
      data1d(ii-1,3)=data1d(ii-1,3)+(vdr+vap)*deltaf*angmom*dp1
!$acc atomic update
      data1d(ii,3)=data1d(ii,3)+(vdr+vap)*deltaf*angmom*(1.0-dp1)

!!! volume averaged: density,entropy,flow,energy,fluxes of particle,momentum,heat
!$acc atomic update
      diagpart(1)=diagpart(1)+deltaf
!$acc atomic update
      diagpart(2)=diagpart(2)+deltaf*deltaf
!$acc atomic update
      diagpart(3)=diagpart(3)+angmom
!$acc atomic update
      diagpart(4)=diagpart(4)+angmom*deltaf
!$acc atomic update
      diagpart(5)=diagpart(5)+energy
!$acc atomic update
      diagpart(6)=diagpart(6)+energy*deltaf
! use EXB and magnetic flutter velocities to calculate total particle
! flux. L. Shi 04/17/2017
!$acc atomic update
      diagpart(7)=diagpart(7)+(vdr+vap)*deltaf
! use total radial velocities to calculate toroidal momentum flux and
! energy flux. L. Shi 05/02/2017
!$acc atomic update
      diagpart(8)=diagpart(8)+(vdr+vap)*angmom*deltaf
!$acc atomic update
      diagpart(9)=diagpart(9)+(vdr+vap)*energy*deltaf
    enddo
!$acc end parallel

!$acc update host(dden,data1d,diagpart)

    diagpart(10)=real(mp)

! sum over all MPI processes
    icount=mpdiag
    diagtmp=0.0
    call MPI_REDUCE(diagpart,diagtmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    diagpart=diagtmp

    icount=(mpsi+1)*mpdata1d
    data1dtmp=0.0
    call MPI_ALLREDUCE(data1d,data1dtmp,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

    icount=mpsi+1
    call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

! radial profile data normalized by marker #
!$omp parallel do private(i,m)
    do i=1,mpdata1d
      do m=0,mpsi
        data1d(m,i)=data1dtmp(m,i)/max(1.0,ddum(m))
      enddo
    enddo
  endif
end subroutine hybridPushParticle
