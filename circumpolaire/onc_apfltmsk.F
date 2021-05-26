CDIR$ LIST
      SUBROUTINE onc_apfltmsk(pvalue,cdmask,pmissing,kndim,klevel,pbuff)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE onc_apfltmsk
CCC                       ********************
CCC file:%M% (%F%) version:%I% %R% %S% modified at %H% %T%
CCC
CCC  PURPOSE : to apply a mask onto a variable.
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
CC      original : creation at 04-02-98 09:43:38 by thierry philipovitch
CC----------------------------------------------------------------------
CC      IMPLICIT NONE
CC parameters and commons
CC ======================
CDIR$ NOLIST
# include "parameter.h"
# include "common.h"
CDIR$ LIST
      CHARACTER*5   cdmask
      REAL          pvalue(*),pmissing,pbuff(*)
      INTEGER       kndim,klevel
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER*128 clbuff
      INTEGER       istart,ilevel,iend
CC----------------------------------------------------------------------
CC intrinsic functions
CC ===================
CC----------------------------------------------------------------------
CC external functions
CC ==================
      CHARACTER*10  cfitostr
      LOGICAL       lfnullstr
      EXTERNAL      cfitostr,lfnullstr
      CHARACTER*128 cfflttostr
      EXTERNAL      cfflttostr
CC----------------------------------------------------------------------
CC statement functions
CC ===================
CC----------------------------------------------------------------------
CC DATA statements
CC ===============
CCC---------------------------------------------------------------------
CCC  package:cdfpack3 module:onc_apfltmsk.F, CETIIS 97
CCC---------------------------------------------------------------------

      IF (lfnullstr(cdmask)) THEN
          clbuff = 'Le nom du masque est vide.'
          CALL procerr('bad_arg','onc_apfltmsk',clbuff)
      ENDIF

# if defined(opaCDF_debug) 
      clbuff = 'Application du masque '
      CALL concattxt(clbuff,cdmask)
      CALL concattxt(clbuff,' avec pour missing value:')
      CALL concattxt(clbuff,cfflttostr(pmissing))
      CALL concattxt(clbuff,' au niveau ')
      CALL concattxt(clbuff,cfitostr(klevel))
      CALL procerr('debug','onc_apfltmsk',clbuff)
# endif

      IF(kndim .EQ. 2) THEN
          IF (klevel .LE. 0) THEN
              CALL procerr('bab_arg','onc_apfltmsk',
     $        'Vous devez donner un niveau dans le cas dim=2')
          ELSE IF (klevel .GT. JPK) THEN
              clbuff = 'Le niveau '
              CALL concattxt(clbuff,cfitostr(klevel))
              CALL concattxt(clbuff,
     $        ' donne est plus grd que JPK=')
              CALL concattxt(clbuff,cfitostr(JPK))
              CALL procerr('bad_arg','onc_apfltmsk',clbuff)
          ELSE
              ilevel = klevel
              istart = 1 + JPI * JPJ * (klevel - 1)
              iend   = JPI * JPJ
          ENDIF
      ELSE
          IF(kndim .NE. 3) THEN
              clbuff = 'Le nombre de dimensions: '
              CALL concattxt(clbuff,cfitostr(kndim))
              CALL concattxt(clbuff,' doit etre 2 ou 3')
              CALL procerr('bad_arg','onc_apfltmsk',clbuff)
          ENDIF
          ilevel = 1
          istart = 1
          iend   = JPI * JPJ * JPK
      ENDIF
      IF (cdmask .EQ. 'TMASK') THEN
          CALL applyfltmask(pbuff,pvalue(istart),
     $         TMASK(1,1,ilevel),pmissing,iend)
      ELSE IF (cdmask .EQ. 'UMASK') THEN
          CALL applyfltmask(pbuff,pvalue(istart),
     $         UMASK(1,1,ilevel),pmissing,iend)
      ELSE IF (cdmask .EQ. 'VMASK') THEN
          CALL applyfltmask(pbuff,pvalue(istart),
     $         VMASK(1,1,ilevel),pmissing,iend)
      ELSE IF (cdmask .EQ. 'FMASK') THEN
          CALL applyfltmask(pbuff,pvalue(istart),
     $         FMASK(1,1,ilevel),pmissing,iend)
      ELSE
          clbuff = 'Le masque donne: '
          CALL concattxt(clbuff,cdmask)
          CALL concattxt(clbuff,' est inconnu.')
          CALL procerr('bad_arg','onc_apfltmsk',clbuff)
      ENDIF
      END
      
