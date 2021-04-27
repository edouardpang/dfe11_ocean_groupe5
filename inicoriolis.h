!
      DO jj = 1, jpj
        DO ji = 1, jpi
#if defined betaPlane
          f(ji,jj) = f0 + beta * (gphif(ji,jj)-gphif(ji,jpj/2))
#else
          f(ji,jj) = f0
#endif
        ENDDO 
      ENDDO
!
      CALL cdfparammsg(cdfparameters,'theta',theta)
      CALL cdfparammsg(cdfparameters,'beta',beta)
