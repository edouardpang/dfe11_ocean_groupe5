!
!
! 2. Initialize initial values for dynamics
! -----------------------------------------
!
! -----------------------
! START OF USER SECTION 2
! -----------------------
!
      icpp1 = 0 ! pour compter les clefs cpp pour la dynamique
      icpp2 = 0 ! 
!
#if defined dynFreeConfig
!
      STOP 'Define your initial state'
      icpp1 = icpp1 + 1
!
#endif
!
!
! Initialize the wind stress
!
      taux(:,:) = 0.0
      tauy(:,:) = 0.0
#if defined windStress
!
! Vent pour le Gulf Stream simple gyre (decommenter les 3 lignes)
!
!     aux = - 2. * pi / (gphif(1,jpjm1) - gphif(1,1))
!     taux(:,:) = 0.078608 * sin(0.042513 * (gphif(:,:) - gphif(:,1)))
!     tauy(ji,jj) = 0.0
!
! Alizés avec minimum à l'équateur et composante méridienne = dérivé de la composante zonale ! A vérifier
!
      aux = 0.05 * (gphif(1,jpjm1) - gphif(1,1))
      xlh = 800000.
      DO jj = 1,jpj/2
        DO ji = 1,jpi
          zaux1 = gphiu(ji,jj) - gphif(1,jpj/2) + xlh
          zaux2 = exp( - zaux1 **2 / aux / aux )
          taux(ji,jj) = - tau * zaux2
          tauy(ji,jj) = - tau * zaux2 * (-2.) * zaux1 /aux 
        ENDDO
      ENDDO
      DO jj = jpj/2+1,jpj
        DO ji = 1,jpi
          zaux1 = gphiu(ji,jj) - gphif(1,jpj/2) - xlh
          zaux2 = exp( - zaux1 **2 / aux / aux )
          taux(ji,jj) = - tau * zaux2
          tauy(ji,jj) = - tau * zaux2 * (-2.) * zaux1 /aux
        ENDDO
      ENDDO
!
! Vent zonal gaussien centré au milieu du bassin
!
!     zaux0 = 0.20 * (gphif(1,jpjm1) - gphif(1,1))
!     DO jj = 1,jpj
!       DO ji = 1,jpi
!         zaux1 = gphiu(ji,jj) - gphif(1,jpj/2)
!         zaux2 = exp( - (zaux1/zaux0) **2 )
!         taux(ji,jj) = - tau * zaux2
!       ENDDO
!     ENDDO
!
      CALL cdfparammsg(cdfparameters,'tau',tau)
      CALL cdfparammsg(cdfparameters,'itau_stop',FLOAT(itau_stop))
      icpp2 = 1
#endif
!
! Initialize a step
!
#if defined dynGill
!
      ipi2 = jpi / 2
      ipi4 = jpi / 4
      ipi2p1 = jpi / 2 + 1
      h00 = 0.01
      hb = -h00 ! on fixe la hauteur partout et ci-dessous, h est augmenté à h00 dans une région donnée
!
! 1 front au centre du bassin, les hautes pressions à gauche
!
      hb(2:ipi2,:,ik) = h00 
!
! 2 fronts de meme amplitude de part et d'autre du milieu du bassin
!
!     hb(ipi2-ipi4:ipi2+ipi4,:,ik) = h00 
!
      CALL cdfparammsg(cdfparameters,'h00',h00)
      icpp1 = icpp1 + 1
!
#endif
!
! Initialize a warm pool
!
#if defined dynWarmPool
!
#if ! defined realWorld
      WRITE(*,*)'La clef dynWarmPool demande que la clef realWorld soit active.'
      WRITE(*,*)'STOP'
      STOP
#endif
      ij0 = jpj/2
      id0 = 10
      h00 = 0.01
! Comme ci-dessus, mais ici on impose une stricte symétrie par rapport à l'équateur 
! La 'warmpool' se trouve à l'Ouest du bassin
      hb(1:2*id0,ij0-id0:ij0+id0+1,ik) = h00
      CALL cdfparammsg(cdfparameters,'h00',h00)
      icpp1 = icpp1 + 1
!
#endif
!
! Initialize the velocity profile of ex 4.6
!
#if defined dynOBC
!
      IF (iiperio .NE. 0) THEN
          WRITE(*,*)'La clef obc demande que iiperio = 0, mais iiperio = ',iiperio
          WRITE(*,*)'STOP'
          STOP
      ENDIF
!
      il0 = 11
      xu0 = 0.4
      xl0 = float(il0)
      u0 = 0.
      v0 = 0.
      h0 = 0.
      zku(1,:,1) = 0.
!
      DO jj=2,il0
        u0(jj) = xu0/(xl0-2.)*(xl0-jj)
      ENDDO
      DO jj=jpj,2,-1
        h0(jj) = h0(jj+1)+f(1,jj)/g*(u0(jj+1)+u0(jj))/2.*dy
      ENDDO
      DO jj=3,jpjm1
        zku(1,jj,1) = (u0(jj+1)+u0(jj-1)-2.*u0(jj))/dy/dy
      ENDDO
      DO jj=2,jpjm1
        v0(jj+1) = -v0(jj) -2./f(1,jj)*xnu*zku(1,jj,1)
      ENDDO
      do jj=1,jpj
        print '(1x,i2,1x,4(E10.4,1X))',jj,U0(jj),v0(jj),zku(1,jj,1),h0(jj) 
      enddo
!
      hb(1,:,jk) = h0(:)
      ub(1,:,jk) = u0(:)
      hb(2,:,jk) = h0(:)
      ub(2,:,jk) = u0(:)
      icpp1 = icpp1 + 1
!
#endif
!
! Initialize a kelvin wave
!
#if defined dynKelvinWave 
!
      IF (iiperio .NE. 1) THEN
          WRITE(*,*)'La clef kelvinWave demande que iiperio = 1, mais iiperio = ',iiperio
          WRITE(*,*)'STOP'
          STOP
      ENDIF
!
      h00 = 0.01 ! amplitude de l'onde
!
      zx = (jpi-2) * dx ! La longueur d'onde est celle du canal. Ajuster si nécessaire
      xk = 2.*acos(-1.) / zx ! nombre d'onde
      ros = sqrt(g*prof(2,2,1))/f(1,1)
!
      DO jj = 1, jpj
        DO ji = 1, jpi
          xi  = (ji-1) * dx
          xp5 = (ji-0.5) * dx
          hb(ji,jj,ik) = h00 * exp (-(jj-2)*dy/ros)*sin(xi*xk)
          ub(ji,jj,ik) = h00 * exp (-(jj-2)*dy/ros)*sin(xp5*xk)/ros/f(1,1)*g ! equilibre geostrophique
        ENDDO
      ENDDO
      icpp1 = icpp1 + 1
!
#endif
!
! Initialize a rossby wave
!
#if defined dynRossbyWave
!
      IF (iiperio .NE. 1) THEN
          WRITE(*,*)'La clef rossbyWave demande que iiperio = 1, mais iiperio = ',iiperio
          WRITE(*,*)'STOP'
          STOP
      ENDIF
!
      zx = (jpi-2) * dx ! La longueur d'onde selon x est celle du canal. Ajuster si nécessaire
      zy = (jpj-2) * dy
!
      h00 = 0.01 ! amplitude de l'onde
      u00 = 0.04 ! on rajoute un courant moyen 
!
      xk = 2.*acos(-1.) / zx ! nombre d'onde
      xl = acos(-1.)/ zy ! Il faut que v soit nul aux frontières N et S
!
      DO jj = 2, jpj
        DO ji = 2, jpi
          hb(ji,jj,ik) =  cos(glamt(ji,jj)*xk)*sin(gphit(ji,jj)*xl)*h00 - (f(ji,jj)+f(ji,jj-1))/2./g*u00*gphit(ji,jj)
          ub(ji,jj,ik) = -cos(glamu(ji,jj)*xk)*cos(gphiu(ji,jj)*xl)*h00 * xl * g / (f(ji,jj-1)+f(ji,jj)) * 2. +u00
          vb(ji,jj,ik) = -sin(glamv(ji,jj)*xk)*sin(gphiv(ji,jj)*xl)*h00 * xk * g / (f(ji-1,jj)+f(ji,jj)) * 2.
        ENDDO
      ENDDO
      CALL cdfparammsg(cdfparameters,'h00',h00)
      CALL cdfparammsg(cdfparameters,'u00',u00)
      icpp1 = icpp1 + 1
!
#endif
!
! Initialize a constant transport zonal flow
!
#if defined dynZonalFlow
!
      IF (iiperio .NE. 1) THEN
          WRITE(*,*)'La clef zonalFlow demande que iiperio = 1, mais iiperio = ',iiperio
          WRITE(*,*)'STOP'
          STOP
      ENDIF
!
      h00 = 0.01
      u00 = 1.
!
      DO jj = 2, jpj
        DO ji = 2, jpi
          ub(ji,jj,ik) = u00 * profu(1,jj,ik)/profu(ji,jj,ik)
          hb(ji,jj,ik) = -(f(ji,jj)+f(ji,jj-1))/ 2./g*ub(ji,jj,ik)*gphit(ji,jj)
        ENDDO
      ENDDO
      CALL cdfparammsg(cdfparameters,'h00',h00)
      CALL cdfparammsg(cdfparameters,'u00',u00)
      icpp1 = icpp1 + 1
!
#endif
!
      IF (icpp1 .NE. 1 .AND. (icpp1 .EQ. 0 .AND. icpp2 .EQ. 0) ) THEN
          WRITE(*,*)'Il faut initialiser la dynamique, et icpp1/icpp2 = ',icpp1,'/',icpp2
          WRITE(*,*)'STOP'
          STOP
      ENDIF
!
! ---------------------
! END OF USER SECTION 2
! ---------------------
!