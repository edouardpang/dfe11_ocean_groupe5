      SUBROUTINE init
!
#include "parameter.h"
#include "common.h"

      CHARACTER*255 cdfparameters /''/
      CHARACTER*255 clbuff
      CHARACTER*7 car7
      REAL omega, f0, beta, theta, rday, pi, oceanMeanDepth, tau, ra, h00
      REAL x, y, fx, fy, fxp, fyp
      REAL aux, xlh, zaux0, zaux1, zaux2
      INTEGER i, j
      INTEGER icpp1, icpp2, ii0, ij0, id0, it
      INTEGER ipi2, ipi2p1, ipi4
      INTEGER ijlim
!
      NAMELIST /physics/ omega, g, theta, beta, oceanMeanDepth, islip, xnu, xnuh, xlamda, tau, itau_stop
      NAMELIST /space_axis/ dx, dy, ii0, ij0
      NAMELIST /time_axis/ dt, nit000, nitend, nwrite
!
      COMMON /opaparameters/ cdfparameters
!
      fx(x,y)  = dx * ( x - 0.5 - FLOAT(ii0) )
      fy(x,y)  = dy * ( y - 0.5 - FLOAT(ij0) )
      fxp(x,y) = dx
      fyp(x,y) = dy
!
! -----------------------
! START OF USER SECTION 1
! -----------------------
!

!
! 0. Variables initializations
! ----------------------------
!
! useful parameters
      pi = acos(-1.)
      rday = 86400.
      ra = 6400000.
!
! time step and horizontal resolution
! Mean depth, gravity and coriolis
!
#if defined dynWarmPool
! equateur
      theta = 0.
#else
! 45 degré Nord
      theta = 45. * pi / 180.
#endif
!
#if defined realWorld
!
      dx = 50000.
      dy = 50000.
      dt = 5000.
! 
!     g = 10. 
      g = 0.01
!
      omega = 2 * pi / rday
      f0 = 2. * omega * sin(theta)
      beta = 2. * omega / ra * cos(theta)
!
      oceanMeanDepth = 1000.
!
#else
!
      dx = 0.02
      dy = 0.02
      dt = 0.02
!
      g = 0.01
! pour avoir f=1
      omega = 0.5
      f0 = 2. * omega
      beta = 1./10.
! 
      oceanMeanDepth = 1.
!
#endif
!
! Parameter for the time loop and numerical schemes
!
      nit000 = 1
      nitend = 500
      nwrite = 100
      gamm = 0.01
!
! format de ndate0 yyyymmdd
! nbisex =  0  'no leap year'
! nbisex = +1  'with leap year'
! nbisex = -1  'perpetual year' (360 days in 12 months of 30 days)
!
      ndate0 = 19700101
      nbisex = -1
      raajj = 360.
      raamo = 12.
      rjjhh = 24.
      rhhmm = 60.
      rmmss = 60.
! Tipical value for the wind 
      tau = 0.0005
      itau_stop = nitend + 1 ! stops the wind after itau_stop
! islip = 1 pour glissement, islip = 0 pour frottement latéral
      islip = 1
! Coefficient for linear drag at the bottom
      xlamda = 0.0
! Arbitrary value for viscosity which ensure numerical stability
      xnu = dx*dy/800./dt
      xnuh = 0.
!
! Read namelists
! 
      OPEN(10,file='namelist.txt')
      READ(10,nml=physics)
      READ(10,nml=time_axis)
      READ(10,nml=space_axis)
      CLOSE(10)
      WRITE(*,nml=physics)
      WRITE(*,nml=time_axis)
      WRITE(*,nml=space_axis)
!
! Now dt is final
!
      rdt = dt
      dt2 = 2.*dt
!
#if defined dynRestart
!
      clbuff = 'Restart  at time step '
      WRITE(car7,'(I7)')nit000
      CALL concattxt(clbuff,car7)
      CALL procerr('info','main',clbuff)
!
#endif
!
!
! 1. Domain initializations
! -------------------------
! 
! Define the horizontal grid including x or y periodicity
!
#include "inigrid.h"
!
! Needed only for the netcdf file 
!
      gdept(1)=1.
      gdepw(1)=0.
      glamt4 = glamt
      glamu4 = glamu
      glamv4 = glamv
      glamf4 = glamf
      gphit4 = gphit
      gphiu4 = gphiu
      gphiv4 = gphiv
      gphif4 = gphif
      gdept4 = gdept
      gdepw4 = gdepw
!
! Define de geometry of the basin and its depth
!
#include "inidepth.h"
!
      DO jk=1,jpk
        IF (iiperio.EQ.1) THEN
            prof( 1 ,:,jk) = prof(jpim1,:,jk)
            prof(jpi,:,jk) = prof(  2  ,:,jk)
        ENDIF
        IF (ijperio.EQ.1) THEN
            prof(:, 1 ,jk) = prof(:,jpjm1,jk)
            prof(:,jpj,jk) = prof(:,  2  ,jk)
        ENDIF
      ENDDO
!
!
! Define relaxation coefficient
!
#include "inirelax.h"
!
! Define the Coriolis parameter
!
#include "inicoriolis.h"
!
! ---------------------
! END OF USER SECTION 1
! ---------------------
!
      IF (iiperio.EQ.1) THEN
          f( 1 ,:) = f(jpim1,:)
          f(jpi,:) = f(  2  ,:)
      ENDIF
      IF (ijperio.EQ.1) THEN
          f(:, 1 ) = f(:,jpjm1)
          f(:,jpj) = f(:,  2  )
      ENDIF
!
      DO jj = 1, jpj
        DO ji = 1, jpim1
          profu(ji,jj,:) = (prof(ji+1,jj,:)+prof(ji,jj,:))/2.
        ENDDO
      ENDDO
      profu(jpi,:,:) = profu(jpim1,:,:)
      DO jj = 1, jpjm1
        DO ji = 1, jpi
          profv(ji,jj,:) = (prof(ji,jj+1,:)+prof(ji,jj,:))/2.
        ENDDO
      ENDDO
      profv(:,jpj,:) = profv(:,jpjm1,:)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          proff(ji,jj,:) = ( prof(ji,jj,:)+prof(ji+1,jj+1,:) + prof(ji+1,jj,:)+prof(ji,jj+1,:) )/4.
        ENDDO
      ENDDO
      DO ji = 1, jpim1
        proff(ji,jpj,:) = proff(ji,jpjm1,:)
      ENDDO
      DO jj = 1, jpjm1
        proff(jpi,jj,:) = proff(jpim1,jj,:)
      ENDDO
      proff(jpi,jpj,:) = proff(jpi,jpjm1,:)
!
      ub = 0.
      un = 0.
      ua = 0.
      vb = 0.
      vn = 0.
      va = 0.
      hb = 0.
      hn = 0.
      ha = 0.
      rotn = 0.
      rotb = 0.
      potvor(:,:) = f(:,:)/proff(:,:,1)
#if defined diaTrends
      utrd = 0.
      vtrd = 0.
      htrd = 0.
#endif
!
! Build masks for Dirichlet boundary conditions
!
      hmask = 0.
      umask = 0.
      vmask = 0.
      fmask = 0.
!
      DO jk = 1,jpk
        WHERE ( mbathy .NE. 0) 
            hmask(:,:,jk) = 1.
        END WHERE
      ENDDO
!
      DO jj = 1,jpjm1
        DO ji = 1,jpim1
          umask(ji,jj,:) = amin1(hmask(ji,jj,:),hmask(ji+1,jj,:))
          vmask(ji,jj,:) = amin1(hmask(ji,jj,:),hmask(ji,jj+1,:))
          fmask(ji,jj,:) = amin1(hmask(ji,jj,:),hmask(ji+1,jj+1,:),hmask(ji+1,jj,:),hmask(ji,jj+1,:))
        ENDDO 
      ENDDO
!
      IF (iiperio.EQ.1) THEN
          hmask( 1 ,:,:) = hmask(jpim1,:,:)
          hmask(jpi,:,:) = hmask(  2  ,:,:)
          umask( 1 ,:,:) = umask(jpim1,:,:)
          umask(jpi,:,:) = umask(  2  ,:,:)
          vmask( 1 ,:,:) = vmask(jpim1,:,:)
          vmask(jpi,:,:) = vmask(  2  ,:,:)
          fmask( 1 ,:,:) = fmask(jpim1,:,:)
          fmask(jpi,:,:) = fmask(  2  ,:,:)
          profu( 1 ,:,:) = profu(jpim1,:,:)
          profu(jpi,:,:) = profu(  2  ,:,:)
          profv( 1 ,:,:) = profv(jpim1,:,:)
          profv(jpi,:,:) = profv(  2  ,:,:)
          proff( 1 ,:,:) = proff(jpim1,:,:)
          proff(jpi,:,:) = proff(  2  ,:,:)
      ENDIF
      IF (ijperio.EQ.1) THEN
          hmask(:, 1 ,:) = hmask(:,jpjm1,:)
          hmask(:,jpj,:) = hmask(:,  2  ,:)
          umask(:, 1 ,:) = umask(:,jpjm1,:)
          umask(:,jpj,:) = umask(:,  2  ,:)
          vmask(:, 1 ,:) = vmask(:,jpjm1,:)
          vmask(:,jpj,:) = vmask(:,  2  ,:)
          fmask(:, 1 ,:) = fmask(:,jpjm1,:)
          fmask(:,jpj,:) = fmask(:,  2  ,:)
          profu(:, 1 ,:) = profu(:,jpjm1,:)
          profu(:,jpj,:) = profu(:,  2  ,:)
          profv(:, 1 ,:) = profv(:,jpjm1,:)
          profv(:,jpj,:) = profv(:,  2  ,:)
          proff(:, 1 ,:) = proff(:,jpjm1,:)
          proff(:,jpj,:) = proff(:,  2  ,:)
      ENDIF
! 
! Save the first parameters in nc file
!
      CALL cdfparammsg(cdfparameters,'g',g)
      CALL cdfparammsg(cdfparameters,'omega',omega)
      CALL cdfparammsg(cdfparameters,'beta',beta)
      CALL cdfparammsg(cdfparameters,'oceanMeanDepth',oceanMeanDepth)
      CALL cdfparammsg(cdfparameters,'xlamda',xlamda)
      CALL cdfparammsg(cdfparameters,'xnu',xnu)
!
!
! 2. Initialize initial values for dynamics
! -----------------------------------------
!
! -----------------------
! START OF USER SECTION 2
! -----------------------
!
#include "inidyn.h"
!
! ---------------------
! END OF USER SECTION 2
! ---------------------
!
#if ! defined dynOBC
      hb = hb * hmask
      ub = ub * umask
      vb = vb * vmask
#endif
!
      IF (iiperio.EQ.1) THEN
          hb( 1 ,:,ik) = hb(jpim1,:,ik)
          hb(jpi,:,ik) = hb(  2  ,:,ik)
          ub( 1 ,:,ik) = ub(jpim1,:,ik)
          ub(jpi,:,ik) = ub(  2  ,:,ik)
      ENDIF
      IF (ijperio.EQ.1) THEN
          hb(:, 1 ,ik) = hb(:,jpjm1,ik)
          hb(:,jpj,ik) = hb(:,  2  ,ik)
          ub(:, 1 ,ik) = ub(:,jpjm1,ik)
          ub(:,jpj,ik) = ub(:,  2  ,ik)
      ENDIF
!
      hn = hb
      un = ub
      vn = vb
!
! Restart
!
#if defined dynRestart
!
      READ(51)un,ub,vn,vb,hn,hb,rotb,rotn,it
      IF( it+1 .NE. nit000) THEN
          WRITE(*,*)'Restart mais mauvais pas de temps nit000/lu dans restart = ',nit000,'/',it+1
          WRITE(*,*)'STOP'
          CALL EXIT(60)
      ENDIF
      CLOSE(51)
!
#endif
!
      RETURN
!
      END
