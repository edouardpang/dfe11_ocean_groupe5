CDIR$ LIST
      SUBROUTINE onc_putstdvar(cdncfil)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE onc_putstdvar
CCC                       *********************
CCC file:%M% (%F%) version:%I% %R% %S% modified at %H% %T%
CCC
CCC  PURPOSE :
CCC  ---------
CCC     
CCC
CC   METHOD :	
CC   -------
CC      
CC        
CC        
CC   USAGE :	
CC   -------
CC      
CC        
CC
CC   INPUT :                    
CC   -----
CC
CC   OUTPUT :		
CC   ------
CC
CC   EXTERNAL :                 
CC   -------- 			 
CC				 
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : creation at 24-10-97 15:17:21 by thierry philipovitch
CC----------------------------------------------------------------------
CC      IMPLICIT NONE
CC parameters and commons
CC ======================
CDIR$ NOLIST
# include "parameter.h"
# include "common.h"
CDIR$ LIST
      CHARACTER*(*) cdncfil
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      REAL bsfn, tn, sn, wn
      REAL zbuff(JPI,JPJ,JPK)
#if defined key_sigma || defined key_s_coord
      REAL zbufv(JPI,JPJ,JPK)
      REAL zvv(JPI,JPJ,JPK)
#endif
CCC---------------------------------------------------------------------
CCC  package:cdfpack module:onc_putstdvar.F, CETIIS 97
CCC---------------------------------------------------------------------
      CALL onc_putvarflt(cdncfil,'UN',un,zbuff)
      CALL onc_putvarflt(cdncfil,'VN',vn,zbuff)
#if defined key_sigma || defined key_s_coord
      zbuff = 0.
      zbufv = 0.
      DO jk=2,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            zbuff(ji,jj,jk)=(un(ji,jj,jk-1)+un(ji,jj,jk))/4.
            zbufv(ji,jj,jk)=(vn(ji,jj,jk-1)+vn(ji,jj,jk))/4.
          ENDDO
        ENDDO
      ENDDO
C
      DO jk=1,jpk
        DO jj=1,jpjm1
          DO ji=1,jpim1
            zsigi = (gdepu(ji,jj,jk)-gdepu(ji-1,jj,jk))/e1t(ji,jj)
            zsigj = (gdepu(ji,jj,jk)-gdepu(ji,jj-1,jk))/e2t(ji,jj)
            zvv(ji,jj,jk)=wn(ji,jj,jk)
     $            +(zbuff(ji,jj,jk)+zbuff(ji-1,jj,jk))*zsigi
     $            +(zbufv(ji,jj,jk)+zbufv(ji,jj-1,jk))*zsigj
          ENDDO
        ENDDO
      ENDDO
      CALL onc_putvarflt(cdncfil,'WN',zvv,zbuff)
#else
      CALL onc_putvarflt(cdncfil,'WN',wn,zbuff)
#endif
      CALL onc_putvarflt(cdncfil,'TN',tn,zbuff)
      CALL onc_putvarflt(cdncfil,'SN',sn,zbuff)
      CALL onc_putvarflt(cdncfil,'ROTN',rotn,zbuff)
      CALL onc_putvarflt(cdncfil,'BSFN',bsfn,zbuff)
C
      END
