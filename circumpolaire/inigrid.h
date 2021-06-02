! iiperio = 1 means x periodicity
! jiperio = 1 means y preiodicity
      iiperio = 1
      ijperio = 0
!
! Change ii0 and ij0 to move the left/bottom corner of the grid
! ii0=-1,ij0=-1 permet d avoir comme origine 
! le point de vorticité d indice (1,1)
      ii0 = +1
      ij0 = +1
! ij0 = jpj/2 pour avoir un bassin symetrique de part et d autre de
! l equateur AVEC JPJ PAIR !! Les points H sont à -dy/2 et +dy/2 de 
! l équateur 
#if defined dynWarmPool
      IF (MOD(jpj,2) .NE. 0) THEN
          WRITE(*,*)'La clef cpp dynWarmPool demande que jpj soit pair.'
          WRITE(*,*)'STOP'
          STOP
      ENDIF
      ij0 = jpj/2
#endif
#
      DO jj = 1, jpj
        DO ji = 1, jpi
          glamt(ji,jj) = fx( FLOAT(ji)    ,FLOAT(jj)     )
          glamu(ji,jj) = fx( FLOAT(ji)+0.5,FLOAT(jj)     )
          glamv(ji,jj) = fx( FLOAT(ji)    ,FLOAT(jj)+0.5 )
          glamf(ji,jj) = fx( FLOAT(ji)+0.5,FLOAT(jj)+0.5 )
        ENDDO 
      ENDDO 
!
      DO jj = 1, jpj
        DO ji = 1, jpi
          gphit(ji,jj) = fy( FLOAT(ji)    ,FLOAT(jj)     )
          gphiu(ji,jj) = fy( FLOAT(ji)+0.5,FLOAT(jj)     )
          gphiv(ji,jj) = fy( FLOAT(ji)    ,FLOAT(jj)+0.5 )
          gphif(ji,jj) = fy( FLOAT(ji)+0.5,FLOAT(jj)+0.5 )
        ENDDO 
      ENDDO 
!
      DO jj = 1, jpj
        DO ji = 1, jpi
          e1t = fxp( FLOAT(ji)    ,FLOAT(jj)     )
          e1u = fxp( FLOAT(ji)+0.5,FLOAT(jj)     )
          e1v = fxp( FLOAT(ji)    ,FLOAT(jj)+0.5 )
          e1f = fxp( FLOAT(ji)+0.5,FLOAT(jj)+0.5 )
        ENDDO
      ENDDO
!
      DO jj = 1, jpj
        DO ji = 1, jpi
          e2t = fyp( FLOAT(ji)    ,FLOAT(jj)     )
          e2u = fyp( FLOAT(ji)+0.5,FLOAT(jj)     )
          e2v = fyp( FLOAT(ji)    ,FLOAT(jj)+0.5 )
          e2f = fyp( FLOAT(ji)+0.5,FLOAT(jj)+0.5 )
        ENDDO
      ENDDO
