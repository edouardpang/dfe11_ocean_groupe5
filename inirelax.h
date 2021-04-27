      relaxh = 0.
      relaxu = 0.
      hinit = 0.
! 
! Exemple de tableau pour la relaxation
! dh/dt = ... - relaxh * (hb - hinit)
! du/dt = ... - relaxu * ub
! dv/dt = ... - relaxu * vb
! 
! Jouer avec des fonction tanh ou des fonction lineaires pour avoir des tableaux continus
! Exemple :
!     DO jj = 1,jpj
!       DO ji = jpi-15,jpi
!         relaxh(ji,jj) = tanh(float(ji-jpi+9))/100./dt2 * tanh(float(jj-5))
!         relaxu(ji,jj) = relaxh(ji,jj) / 10.
!       ENDDO 
!     ENDDO
!
! Ci-dessous une zone de rappel au N et au S symetrique (pour le cas Nino par exemple)
!
!     ijlim = 10
!     DO jj = 2,ijlim
!       DO ji = 2,jpim1
!         relaxh(ji,jj) = 1./300./dt2 * float(ijlim - jj)/float(ijlim - 2)
!         relaxu(ji,jj) = relaxh(ji,jj) / 10.
!         relaxh(ji,jpj - jj + 1) = 1./300./dt2 * float(ijlim - jj)/float(ijlim - 2)
!         relaxu(ji,jpj - jj + 1) = relaxh(ji,jpj - jj + 1) / 10.
!       ENDDO 
!     ENDDO
