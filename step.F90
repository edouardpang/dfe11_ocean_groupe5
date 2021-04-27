      SUBROUTINE step(kt)
!
#include "parameter.h"
#include "common.h"
!
      REAL hmin, hmax
      REAL zua, zuad, zva, zvad, zha
      REAL zx, zxx, zy, zyy
!
!
! 0. Initialize and loop over layers
! ----------------------------------
!
      gu = 0.
      gv = 0.
!
      DO 1000 jk=1,jpk
!
!
! 1. Compute total height on the whole mesh
! -----------------------------------------
!
! compute total water column height at T-point
#if defined nonlinear
        htot(:,:) = hn(:,:,jk) + prof(:,:,jk)
#else
        htot(:,:) = prof(:,:,jk)
#endif
! compute total water column height at U-point
        DO jj = 1, jpj
          DO ji = 1, jpim1
#if defined nonlinear
            htotu(ji,jj) = (hn(ji+1,jj,jk)+hn(ji,jj,jk))/2. + profu(ji,jj,jk)
#else
            htotu(ji,jj) = profu(ji,jj,jk)
#endif
          ENDDO
        ENDDO
        htotu(jpi,:) = htotu(jpim1,:)
! compute total water column height at V-point
        DO jj = 1, jpjm1
          DO ji = 1, jpi
#if defined nonlinear
            htotv(ji,jj) = (hn(ji,jj+1,jk)+hn(ji,jj,jk))/2. + profv(ji,jj,jk)
#else
            htotv(ji,jj) = profv(ji,jj,jk)
#endif
          ENDDO
        ENDDO
        htotv(:,jpj) = htotv(:,jpjm1)
! compute total water column height at F-point
        DO jj = 1, jpj
          DO ji = 1, jpi
#if defined nonlinear
            htotf(ji,jj) = (hn(ji+1,jj+1,jk)+hn(ji,jj,jk)+hn(ji+1,jj,jk)+hn(ji,jj+1,jk))/4. + proff(ji,jj,jk)
#else
            htotf(ji,jj) = proff(ji,jj,jk)
#endif
          ENDDO
        ENDDO
        htotf(:,jpj) = htotf(:,jpjm1)
        htotf(jpi,:) = htotf(jpim1,:)
        htotf(jpi,jpj) = htotf(jpi,jpjm1)

!
! compute potential vorticity
!
#if defined tensorial
        zku(:,:,1) = un(:,:,jk) * e1u(:,:)        
        zkv(:,:,1) = vn(:,:,jk) * e2v(:,:)        
#endif
!
        DO jj = 1,jpjm1
          DO ji = 1,jpim1
#if defined tensorial
            rotn(ji,jj,jk) = ( zkv(ji+1,jj,1)-zkv(ji,jj,1) -zku(ji,jj+1,1)+zku(ji,jj,1) ) &
                                   / e1f(ji,jj) / e2f(ji,jj) 
#else
            rotn(ji,jj,jk) = + ( vn(ji+1,jj,jk)-vn(ji,jj,jk) )/dx - ( un(ji,jj+1,jk)-un(ji,jj,jk) )/dy
#endif
          ENDDO
        ENDDO
! Boundary condition
        IF( islip .EQ. 1) THEN
! slip (i.e.neumann) boundary condition (condition de glissement)
            rotn(:,:,jk) = rotn(:,:,jk) * fmask(:,:,jk)
        ELSE
! no slip boundary condition (condition de frottement u=v=0)
            rotn(:,:,jk) = rotn(:,:,jk) * ( 2. - fmask(:,:,jk) )
        ENDIF
!
        IF (iiperio.EQ.1) THEN
            rotn( 1 ,:,jk) = rotn(jpim1,:,jk)
            rotn(jpi,:,jk) = rotn(  2  ,:,jk)
        ENDIF
        IF (ijperio.EQ.1) THEN
            rotn(:, 1 ,jk) = rotn(:,jpjm1,jk)
            rotn(:,jpj,jk) = rotn(:,  2  ,jk)
        ENDIF
!
#if defined nonlinear
        potvor(:,:) = f(:,:) + rotn(:,:,jk)
#else
        potvor = f
#endif
        potvor = potvor / htotf
!
! compute kinetic energy and bernouilli function
!
        DO jj = 2,jpj
          DO ji = 2, jpi
#if defined nonlinear
            zku(ji,jj,jk) =0.25*( un(ji-1, jj ,jk)*un(ji-1, jj ,jk) &
                                + un( ji , jj ,jk)*un( ji , jj ,jk) &
                                + vn( ji ,jj-1,jk)*vn( ji ,jj-1,jk) &
                                + vn( ji , jj ,jk)*vn( ji , jj ,jk)  )
            zku(ji,jj,jk) = zku(ji,jj,jk) + g * hn(ji,jj,jk)
#else
            zku(ji,jj,jk) =  + g * hn(ji,jj,jk)
#endif
          ENDDO
        ENDDO
!
!
! 3. u-momentum equation
! ----------------------
!
! compute v-momentum flux
!
#if defined tensorial
        vflux(:,:) = vn(:,:,jk) * e1v(:,:) * htotv(:,:)
#else
        vflux(:,:) = vn(:,:,jk) * htotv(:,:)
#endif
!
! store first in ua the temporal derivative 
! 
        DO jj = 2, jpj 
          DO ji = 1, jpim1
#if defined tensorial
            zx = e1u(ji,jj)
            zxx = e1u(ji,jj)
#else
            zx = 1.
            zxx = dx
#endif
!
#if ! defined enstrophy
            zua = + ( potvor(ji,jj)   * ( vflux(ji+1,jj)+vflux(ji ,jj)    ) &
                    + potvor(ji,jj-1) * ( vflux(ji,jj-1)+vflux(ji+1,jj-1) ) ) / zx / 4.
#else
            zua = + ( potvor(ji,jj)+potvor(ji,jj-1) ) * &
                    ( vflux( ji ,jj)+vflux( ji ,jj-1)+vflux(ji+1,jj)+vflux(ji+1,jj-1) ) / zx / 8.
#endif
            zuad = ( zku(ji+1,jj,jk)-zku(ji,jj,jk) ) / zxx
            ua(ji,jj,jk) = zua - zuad
#if defined diaTrends
            zha   = g * ( hn(ji+1,jj,jk)-hn(ji,jj,jk) ) / zxx
            utrd(ji,jj,jk,1) = zua - zuad + zha
            utrd(ji,jj,jk,2) = - zha
#endif
          ENDDO
        ENDDO
!
! Add Wind stress and linear drag
!
      IF ( kt .EQ. itau_stop) THEN
          taux = 0.
          tauy = 0.
      ENDIF
!
      IF(jk.EQ.1)THEN
#if defined diaTrends
          utrd(:,:,jk,3) = taux(:,:) / htotu(:,:)
          ua(:,:,jk) = ua(:,:,jk) + utrd(:,:,jk,3)
#else
          ua(:,:,jk) = ua(:,:,jk) + taux(:,:) / htotu(:,:)
#endif
      ENDIF
      IF(jk.EQ.jpk)THEN
#if defined diaTrends
          utrd(:,:,jk,4) = - xlamda * ub(:,:,jk) / htotu(:,:)
          ua(:,:,jk) = ua(:,:,jk) + utrd(:,:,jk,4)
#else
          ua(:,:,jk) = ua(:,:,jk) - xlamda * ub(:,:,jk) / htotu(:,:)
#endif
      ENDIF
!
! Add a relaxation
!
#if defined diaTrends
        utrd(:,:,jk,5) = - relaxu(:,:) * ub(:,:,jk)
        ua(:,:,jk) = ua(:,:,jk) + utrd(:,:,jk,5)
#else
        ua(:,:,jk) = ua(:,:,jk) - relaxu(:,:) * ub(:,:,jk)
#endif
!
      DO jj = 1, jpj
          if(kt.eq.5000)write(*,*)jj,va(7,jj,1),ua(7,jj,1)
      ENDDO
!
!
! 4. v-momentum equation
! ----------------------
!
! compute u-momentum flux
!
#if defined tensorial
        uflux(:,:) = un(:,:,jk) * e2u(:,:) * htotu(:,:)
#else
        uflux(:,:) = un(:,:,jk) * htotu(:,:)
#endif
!
! store first in va the temporal derivative 
! 
        DO jj = 1, jpjm1 
          DO ji = 2, jpi
#if defined tensorial
            zy = e2v(ji,jj)
            zyy = e2v(ji,jj)
#else
            zy = 1.
            zyy = dy
#endif
!
#if ! defined enstrophy
            zva = - ( potvor(ji,jj)  * ( uflux(ji,jj+1)+uflux(ji,jj) ) &
                    + potvor(ji-1,jj)* ( uflux(ji-1,jj+1)+uflux(ji-1,jj) ) ) / zy / 4.
#else
            zva = - ( potvor(ji-1,jj)+potvor(ji,jj) ) * &
                   ( uflux( ji ,jj)+uflux( ji ,jj+1)+uflux(ji-1,jj)+uflux(ji-1,jj+1) ) / zy / 8.
#endif
            va(ji,jj,jk) = zva - ( zku(ji,jj+1,jk)-zku(ji,jj,jk) ) / zyy
            zvad = ( zku(ji,jj+1,jk)-zku(ji,jj,jk) ) / zyy
            va(ji,jj,jk) = zva - zvad
#if defined diaTrends
            zha   = g * ( hn(ji,jj+1,jk)-hn(ji,jj,jk) ) / zyy
            vtrd(ji,jj,jk,1) = zva - zvad + zha
            vtrd(ji,jj,jk,2) = - zha
#endif
          ENDDO
        ENDDO 
!
! Add Wind stress and linear drag
!
      IF(jk.EQ.1)THEN
#if defined diaTrends
          vtrd(:,:,jk,3) = tauy(:,:) / htotv(:,:)
          va(:,:,jk) = va(:,:,jk) + vtrd(:,:,jk,3)
#else
          va(:,:,jk) = va(:,:,jk) + tauy(:,:) / htotv(:,:)
#endif
      ENDIF
      IF(jk.EQ.jpk)THEN
#if defined diaTrends
          vtrd(:,:,jk,4) = - xlamda * vb(:,:,jk) / htotv(:,:)
          va(:,:,jk) = va(:,:,jk) + vtrd(:,:,jk,4)
#else
          va(:,:,jk) = va(:,:,jk) - xlamda * vb(:,:,jk) / htotv(:,:)
#endif
      ENDIF
!
! Add a relaxation
!
#if defined diaTrends
        vtrd(:,:,jk,5) = - relaxu(:,:) * vb(:,:,jk)
        va(:,:,jk) = va(:,:,jk) + vtrd(:,:,jk,5)
#else
        va(:,:,jk) = va(:,:,jk) - relaxu(:,:) * vb(:,:,jk)
#endif
!
!
! 4. add horizontal diffusion and compute vertical average
! --------------------------------------------------------
!
        CALL qdmdiff(kt)
!
!       gu(:,:)=gu(:,:)+htotu(:,:)*ua(:,:,jk) / (profu(:,:,jk)+eps)
!       gv(:,:)=gu(:,:)+htotv(:,:)*va(:,:,jk) / (profv(:,:,jk)+eps)
!
! 5. continuity equation
! ----------------------
!
! store first in ha the temporal derivative 
! 
!       DO jj = 1, jpjm1
!         DO ji = 1, jpim1
! if an upstream scheme is necessary
!           IF ( un(ji,jj,jk) .ge. 0. ) THEN
!               zku(ji,jj,jk) = un(ji,jj,jk) * htot(ji,jj)
!             ELSE
!               zku(ji,jj,jk) = un(ji,jj,jk) * htot(ji+1,jj)
!           ENDIF
!           IF ( vn(ji,jj,jk) .ge. 0. ) THEN
!               zkv(ji,jj,jk) = vn(ji,jj,jk) * htot(ji,jj)
!             ELSE
!               zkv(ji,jj,jk) = vn(ji,jj,jk) * htot(ji,jj+1)
!           ENDIF
!         ENDDO
!       ENDDO
!       DO jj = 2,jpjm1
!         DO ji = 2, jpim1
!           ha(ji,jj,jk) = - ( zku(ji,jj,1)-zku(ji-1,jj,1) )/dx &
!                          - ( zkv(ji,jj,1)-zkv(ji,jj-1,1) )/dy
!         ENDDO
!       ENDDO 
!
        DO jj = 2,jpjm1
          DO ji = 2, jpim1
#if defined tensorial
            ha(ji,jj,jk) = - ( uflux(ji,jj)-uflux(ji-1,jj) + vflux(ji,jj)-vflux(ji,jj-1) ) &
                               / e1t(ji,jj) / e2t(ji,jj)
#else
            ha(ji,jj,jk) = - ( uflux(ji,jj)-uflux(ji-1,jj) ) / dx &
                           - ( vflux(ji,jj)-vflux(ji,jj-1) ) / dy
#endif
#if defined diaTrends
            htrd(ji,jj,jk,1) = ha(ji,jj,jk) 
#endif
          ENDDO
        ENDDO 
!
! Add a relaxation
!
        ha(:,:,jk) = ha(:,:,jk) - relaxh(:,:) * ( hb(:,:,jk) - hinit(:,:) )
#if defined diaTrends
        htrd(:,:,jk,2) = - relaxh(:,:) * ( hb(:,:,jk) - hinit(:,:) )
#endif
! 
        if(xnuh .NE. 0.)CALL hdiff(kt)
!
 1000 CONTINUE
!
#if defined multilayer
!
! 5. Compute stream function
! --------------------------
!
      DO jk=1,jpk
        ua(:,:,jk) = ua(:,:,jk) - gu(:,:)
        va(:,:,jk) = va(:,:,jk) - gv(:,:)
      ENDDO
!
      CALL psibar(kt)
!
#endif
!
! 6. Update from t to t+dt with leapfrog scheme
! ---------------------------------------------
!
      DO 2000 jk=1,jpk
!
#if ! defined dynRestart
        IF (kt.EQ.nit000) THEN
            zdt = dt
          ELSE
            zdt = dt2
        ENDIF
#else
        zdt = dt2
#endif
!
        ua(:,:,jk) = ub(:,:,jk) + zdt * ua(:,:,jk)
        va(:,:,jk) = vb(:,:,jk) + zdt * va(:,:,jk)
        ha(:,:,jk) = hb(:,:,jk) + zdt * ha(:,:,jk)
! 
!
! 7. Boundary conditions 
! ---------------------- 
!
        ua(:,:,jk) = ua(:,:,jk) * umask(:,:,jk)
        va(:,:,jk) = va(:,:,jk) * vmask(:,:,jk)
        ha(:,:,jk) = ha(:,:,jk) * hmask(:,:,jk)
#if defined dynOBC
        ha(1,:,jk) = h0(:)
        ua(1,:,jk) = u0(:)
!       va(1,:,jk) = v0(:)
        ha(2,:,jk) = h0(:)
        ua(2,:,jk) = u0(:)
!       va(2,:,jk) = v0(:)
#endif
!
        IF (iiperio.EQ.1) THEN
            ha( 1 ,:,jk) = ha(jpim1,:,jk)
            ha(jpi,:,jk) = ha(  2  ,:,jk)
            ua( 1 ,:,jk) = ua(jpim1,:,jk)
            ua(jpi,:,jk) = ua(  2  ,:,jk)
            va( 1 ,:,jk) = va(jpim1,:,jk)
            va(jpi,:,jk) = va(  2  ,:,jk)
        ENDIF
        IF (ijperio.EQ.1) THEN
            ha(:, 1 ,jk) = ha(:,jpjm1,jk)
            ha(:,jpj,jk) = ha(:,  2  ,jk)
            ua(:, 1 ,jk) = ua(:,jpjm1,jk)
            ua(:,jpj,jk) = ua(:,  2  ,jk)
            va(:, 1 ,jk) = va(:,jpjm1,jk)
            va(:,jpj,jk) = va(:,  2  ,jk)
        ENDIF
!
#if defined diaTrends
! on place diag avant le swap des tableaux et Asselin
        CALL diag(kt) ! attention diag est dans la boucle jk. Valable uniquement si jk=1
#endif
!  
!
! 8. Filter with Asselin and swap arrays
! --------------------------------------
!
#if ! defined dynRestart
        IF (kt.EQ.nit000) THEN
            un(:,:,jk) = ua(:,:,jk)
            ua(:,:,jk) = 0.
            vn(:,:,jk) = va(:,:,jk)
            va(:,:,jk) = 0.
            hn(:,:,jk) = ha(:,:,jk)
            ha(:,:,jk) = 0.
          ELSE
#endif
! Asselin
            un(:,:,jk) = un(:,:,jk) + gamm * (ua(:,:,jk)+ub(:,:,jk)-2.*un(:,:,jk))
            vn(:,:,jk) = vn(:,:,jk) + gamm * (va(:,:,jk)+vb(:,:,jk)-2.*vn(:,:,jk))
            hn(:,:,jk) = hn(:,:,jk) + gamm * (ha(:,:,jk)+hb(:,:,jk)-2.*hn(:,:,jk))
! Swap
            ub(:,:,jk) = un(:,:,jk) 
            un(:,:,jk) = ua(:,:,jk) 
            ua(:,:,jk) = 0.
            vb(:,:,jk) = vn(:,:,jk) 
            vn(:,:,jk) = va(:,:,jk) 
            va(:,:,jk) = 0.
            hb(:,:,jk) = hn(:,:,jk) 
            hn(:,:,jk) = ha(:,:,jk) 
            ha(:,:,jk) = 0.
            rotb(:,:,jk) = rotn(:,:,jk)
#if ! defined dynRestart
        ENDIF
#endif
!
! End of jk loop
!
 2000 CONTINUE
!
      RETURN
!
      END


      SUBROUTINE hdiff(kt)
!
#include "parameter.h"
#include "common.h"
!
      REAL zha
!
! Loop over layers
! ----------------
!
      DO 1000 jk=1,jpk
!
!
! Fick law with constant viscosity
!
        DO jj = 1,jpjm1
          DO ji = 1,jpim1
#if defined tensorial
            zku(ji,jj,jk) = xnu * (hb(ji+1,jj,jk)-hb(ji,jj,jk)) * e2u(ji,jj) / e1u(ji,jj)
            zkv(ji,jj,jk) = xnu * (hb(ji,jj+1,jk)-hb(ji,jj,jk)) * e1v(ji,jj) / e2v(ji,jj)
#else
            zku(ji,jj,jk) = xnu * (hb(ji+1,jj,jk)-hb(ji,jj,jk))/dx
            zkv(ji,jj,jk) = xnu * (hb(ji,jj+1,jk)-hb(ji,jj,jk))/dy
#endif
          ENDDO
        ENDDO
!
! boundary condition (pas de flux)
!
        zku(:,:,jk) = zku(:,:,jk) * umask(:,:,jk)
        zkv(:,:,jk) = zkv(:,:,jk) * vmask(:,:,jk)
        IF (iiperio.EQ.1) THEN
            zku( 1 ,:,jk) = zku(jpim1,:,jk)
            zku(jpi,:,jk) = zku(  2  ,:,jk)
            zkv( 1 ,:,jk) = zkv(jpim1,:,jk)
            zkv(jpi,:,jk) = zkv(  2  ,:,jk)
        ENDIF
        IF (ijperio.EQ.1) THEN
            zku(:, 1 ,jk) = zku(:,jpjm1,jk)
            zku(:,jpj,jk) = zku(:,  2  ,jk)
            zkv(:, 1 ,jk) = zkv(:,jpjm1,jk)
            zkv(:,jpj,jk) = zkv(:,  2  ,jk)
        ENDIF
!
!
! second derivative
!
        DO jj = 2,jpj
          DO ji = 2,jpi
#if defined tensorial
!           ha(ji,jj,jk) = ha(ji,jj,jk) + &
            zha          = ( zku(ji,jj,jk)-zku(ji-1,jj,jk) + zkv(ji,jj,jk)-zkv(ji,jj-1,jk) ) &
                           / e1t(ji,jj) / e2t(ji,jj)
#else
!           ha(ji,jj,jk) = ha(ji,jj,jk) + &
            zha          = ( (zku(ji,jj,jk)-zku(ji-1,jj,jk))/dx &
                           + (zkv(ji,jj,jk)-zkv(ji,jj-1,jk))/dy )
#endif
            ha(ji,jj,jk) = ha(ji,jj,jk) + zha
#if defined diaTrends
            htrd(ji,jj,jk,3) = zha
#endif
          ENDDO
        ENDDO
!
! End of jk loop
!
 1000 CONTINUE
!
!
      RETURN
!
      END


      SUBROUTINE qdmdiff(kt)

#include "parameter.h"
#include "common.h"
!
      REAL uta, vta, hmin, hmax
!
!
! Loop over layers
! ----------------
!
      zkf = rotb
#if defined tensorial
      DO jk=1,jpk
        zku(:,:,jk) = ub(:,:,jk) * e2u(:,:)
        zkv(:,:,jk) = vb(:,:,jk) * e1v(:,:)
      ENDDO
#else
      zku = ub
      zkv = vb
#endif
!
      DO 1000 jk=1,jpk
!
! Compute curl and div
!
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
#if defined tensorial
            zkh(ji,jj,jk) = ( zku( ji ,jj,jk)-zku(ji-1,jj,jk) + zkv(ji, jj ,jk)-zkv(ji,jj-1,jk) )/ e1t(ji,jj) / e2t(ji,jj) ! div
#else
            zkh(ji,jj,jk) = (zku( ji ,jj,jk)-zku(ji-1,jj,jk))/dx + (zkv(ji, jj ,jk)-zkv(ji,jj-1,jk))/dy ! div
#endif
          ENDDO
        ENDDO
!
! Fick law with constant viscosity
!
        zkf(:,:,jk) = xnu * zkf(:,:,jk)
        zkh(:,:,jk) = xnu * zkh(:,:,jk) * hmask(:,:,jk)
!
! zku is already computed on boundaries and periodicity done
!
        IF (iiperio.EQ.1) THEN
            zkh( 1 ,:,jk) = zkh(jpim1,:,jk)
            zkh(jpi,:,jk) = zkh(  2  ,:,jk)
        ENDIF
        IF (ijperio.EQ.1) THEN
            zkh(:, 1 ,jk) = zkh(:,jpjm1,jk)
            zkh(:,jpj,jk) = zkh(:,  2  ,jk)
        ENDIF
!
! second derivative
!
        DO jj = 2,jpj
          DO ji = 2,jpi
#if defined tensorial
!           ua(ji,jj,jk) = ua(ji,jj,jk) &
!                           - (zkf( ji ,jj,jk)-zkf(ji,jj-1,jk)) / e2u(ji,jj) &
!                           + (zkh(ji+1,jj,jk)-zkh(ji, jj ,jk)) / e1u(ji,jj)
!           va(ji,jj,jk) = va(ji,jj,jk) &
!                           + (zkf(ji, jj ,jk)-zkf(ji-1,jj,jk)) / e1v(ji,jj) &
!                           + (zkh(ji,jj+1,jk)-zkh( ji ,jj,jk)) / e2v(ji,jj)
            uta =           - (zkf( ji ,jj,jk)-zkf(ji,jj-1,jk)) / e2u(ji,jj) &
                            + (zkh(ji+1,jj,jk)-zkh(ji, jj ,jk)) / e1u(ji,jj)
            ua(ji,jj,jk) = ua(ji,jj,jk) + uta
            vta =           + (zkf(ji, jj ,jk)-zkf(ji-1,jj,jk)) / e1v(ji,jj) &
                            + (zkh(ji,jj+1,jk)-zkh( ji ,jj,jk)) / e2v(ji,jj)
            va(ji,jj,jk) = va(ji,jj,jk) + vta
#else
            uta =           - (zkf( ji ,jj,jk)-zkf(ji,jj-1,jk)) / dy &
                            + (zkh(ji+1,jj,jk)-zkh(ji, jj ,jk)) / dx
            ua(ji,jj,jk) = ua(ji,jj,jk) + uta
            vta =           + (zkf(ji, jj ,jk)-zkf(ji-1,jj,jk)) / dx &
                            + (zkh(ji,jj+1,jk)-zkh( ji ,jj,jk)) / dy
            va(ji,jj,jk) = va(ji,jj,jk) + vta
#endif
#if defined diaTrends
            utrd(ji,jj,jk,6) = uta
            vtrd(ji,jj,jk,6) = vta
#endif
          ENDDO
        ENDDO
!
! End of jk loop
!
 1000 CONTINUE
!
      RETURN
!
      END

      SUBROUTINE psibar
!
#include "parameter.h"
#include "common.h"
!
!
! Loop over layers
! ----------------
!
!
      DO jj=1,jpjm1
        DO ji=1,jpim1
          cur(ji,jj)=-(gv(ji+1,jj)-gv(ji,jj))/dx &
                     +(gu(ji,jj+1)-gu(ji,jj))/dy
        ENDDO
      ENDDO 
!
      RETURN
!
      END

      SUBROUTINE minmax(tab, tabm, tmin,tmax)
#include "parameter.h"
      REAL tab(jpi,jpj), tabm(jpi,jpj)
      REAL tmin,tmax
      INTEGER ji, jj
      tmax = -999.
      tmin = 999.
      DO jj = 2,jpjm1
        DO ji = 2,jpim1
!         IF(tabm(ji,jj) .GE. 0.5) THEN
              tmax = amax1(tmax,tab(ji,jj))
              tmin = amin1(tmin,tab(ji,jj))
!         ENDIF
        ENDDO
      ENDDO
!
      END
