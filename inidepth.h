!
! Define a closed rectangular basin
!
      mbathy = 0
      mbathy(2:jpim1,2:jpjm1) = 1
!
!
! Exemples de geometrie simples : decommenter pour activer
!
!
!
! ajout d un mur en epi de la cote Sud 
!
!     iim0 = 10
!     ijm0 = 11
!     mbathy(1:iim0,ijm0)=0.
!
! ajout d une ile en iil0,ijl0 avec un point de part et d autre
!
!     iil0 = 60
!     ijl0 = jpj/2
!     mbathy(iil0,ijl0+1)=0.
!     mbathy(iil0,ijl0)=0.
!     mbathy(iil0,ijl0-1)=0.
!
      IF (iiperio.EQ.1) THEN
          mbathy( 1 ,:) = mbathy(jpim1,:)
          mbathy(jpi,:) = mbathy(  2  ,:)
      ENDIF
      IF (ijperio.EQ.1) THEN
          mbathy(:, 1 ) = mbathy(:,jpjm1)
          mbathy(:,jpj) = mbathy(:,  2  )
      ENDIF
!
!
! Define the depth
!
!
      prof(:,:,ik) = oceanMeanDepth
!
! Cas d un seuil
!
#if defined sill
!
      ijtopo=60
      prof(1:itopo    ,:,ik) = oceanMeanDepth
      prof(itopo+1:jpi,:,ik) = oceanMeanDepth * 0.5
!
#endif
!
! Cas d une montagne dans un canal zonal
!
#if defined dynZonalFlow
!
      ijtopo=50
      DO ji = 1, jpi
        prof(ji,:,ik) = oceanMeanDepth * (1-0.4*exp(-(ji-ijtopo)**2/0.01))
      ENDDO 
!
#endif
!
