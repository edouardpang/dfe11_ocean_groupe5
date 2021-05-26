      SUBROUTINE diag(kt)
!
! Attention : diag est ecrit pour jk=1. Pas valable en multicouche
!
#include "parameter.h"
#include "common.h"
!
#if defined diaTrends
!
      REAL*4 zeptrd(3), zectrd(6), zep, zec, zens
      INTEGER j4
      REAL, save :: zepn1 = 0.   
      REAL, save :: zepn2 = 0.   
!
        htrd(:,:,:,1) = htrd(:,:,:,1)*hmask(:,:,:)
        htrd(:,:,:,2) = htrd(:,:,:,2)*hmask(:,:,:)
        htrd(:,:,:,3) = htrd(:,:,:,3)*hmask(:,:,:)
        DO j4 = 1, 6
          utrd(:,:,:,j4) = utrd(:,:,:,j4)*umask(:,:,:)
          vtrd(:,:,:,j4) = vtrd(:,:,:,j4)*vmask(:,:,:)
        ENDDO
!
! diagnose integral momentum and thickness
! TO BE DONE
!
! diagnose potential and kinetic energy and transfers
!
        zectrd = 0.
        zeptrd = 0.
! Potential energy
        zkh(:,:,1) = ha(:,:,1) * ha(:,:,1)
        zep = g * SUM(zkh(:,:,1)) / 2. 
! Kinetic 
        zku(:,:,1) = ua(:,:,1) * ua(:,:,1) * htotu(:,:)
        zkv(:,:,1) = va(:,:,1) * va(:,:,1) * htotv(:,:)
        zec = ( SUM(zku(:,:,1)) + SUM(zkv(:,:,1)) ) / 2.
! Enstrophy
        zkf(:,:,1) = rotn(:,:,1) * rotn(:,:,1)
        zens = SUM(zkf(:,:,1)) / 2.
! Transfer potential->kinetic
        zkh(:,:,1) = htrd(:,:,1,1) * hn(:,:,1)
        zeptrd(1) = g * SUM(zkh(:,:,1)) ! Transfer potential->kinetic
! PE input/relax
        zkh(:,:,1) = htrd(:,:,1,2) * hn(:,:,1)
        zeptrd(2) = SUM(zkh(:,:,1))
! PE dissipation
        zkh(:,:,1) = htrd(:,:,1,3) * hn(:,:,1)
        zeptrd(3) = SUM(htrd(:,:,1,3))
! KE trend for coriolis and non-linear terms (=0 if energy conserving)
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,1)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,1)
        zectrd(1) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
! Transfer kinetic->potential
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,2)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,2)
        zectrd(2) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
! KE input from wind
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,3)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,3)
        zectrd(3) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
! KE dissipation due to bottom friction
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,4)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,4)
        zectrd(4) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
! KE dissipation due to damping
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,5)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,5)
        zectrd(5) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
! KE dissipation due to lateral diffusion
        zku(:,:,1) = un(:,:,1) * htotu(:,:) * utrd(:,:,1,6)
        zkv(:,:,1) = vn(:,:,1) * htotv(:,:) * vtrd(:,:,1,6)
        zectrd(6) = SUM(zku(:,:,1)) + SUM(zkv(:,:,1))
!
        write(90,*)kt, zep, zec, zens, zeptrd, zectrd
        zepn2 = zepn1
        zepn1 = zep
!
#endif
!
      RETURN
!
      END
