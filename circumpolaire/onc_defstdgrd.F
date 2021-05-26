CDIR$ LIST
      SUBROUTINE onc_defstdgrd(cdncfil)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE onc_defstdgrd
CCC                       **********************
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
CC      original : creation at 03-02-98 16:42:32 by thierry philipovitch
CC      modif.   : coordonnees sigma 25-03-99 by laurent mortier
CC----------------------------------------------------------------------
CC      IMPLICIT NONE
CC parameters and commons
CC ======================
CDIR$ NOLIST
# include "parameter.h"
# include "common.h"
# include "cdfpack.h"
CDIR$ LIST
      CHARACTER*(*) cdncfil
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER*_LG_DIM_NAME_ cldimn(4)
      INTEGER   incnum,incid,ivarid
      CHARACTER*80 clbuff
CC----------------------------------------------------------------------
CC external functions
CC ==================
      INTEGER  kfcdfaddvar
      EXTERNAL kfcdfaddvar
CCC---------------------------------------------------------------------
CCC  package:cdfpack3 module:onc_defstdgrd.F, CETIIS 97
CCC---------------------------------------------------------------------
      CALL onc_chkinit('onc_defstdgrd')
      IF (ngeoref .EQ. 2) THEN
C
C 1.informations about the nc file:
C ---------------------------------
C
          CALL onc_setncf(cdncfil,incnum,incid)
C
C 4.Axis definition:
C ------------------
C
          cldimn(1) = 'JPI'
          cldimn(2) = 'JPJ'
          ivarid = kfcdfaddvar(incid,'GLAMT',4004,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GLAMU',4005,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GLAMV',4006,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GLAMF',4007,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GPHIT',4008,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GPHIU',4009,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GPHIV',4010,'real',2,cldimn)
          ivarid = kfcdfaddvar(incid,'GPHIF',4011,'real',2,cldimn)
          CALL cdfputfltval(incid,'GLAMT',-1,GLAMT4,lflag_sync)
          CALL cdfputfltval(incid,'GLAMU',-1,GLAMU4,lflag_sync)
          CALL cdfputfltval(incid,'GLAMV',-1,GLAMV4,lflag_sync)
          CALL cdfputfltval(incid,'GLAMF',-1,GLAMF4,lflag_sync)
          CALL cdfputfltval(incid,'GPHIT',-1,GPHIT4,lflag_sync)
          CALL cdfputfltval(incid,'GPHIU',-1,GPHIU4,lflag_sync)
          CALL cdfputfltval(incid,'GPHIV',-1,GPHIV4,lflag_sync)
          CALL cdfputfltval(incid,'GPHIF',-1,GPHIF4,lflag_sync)
#if defined key_sigma || defined key_s_coord
          cldimn(1) = 'JPI'
          cldimn(2) = 'JPJ'
          cldimn(3) = 'JPK'
          ivarid = kfcdfaddvar(incid,'GDEPT',4012,'real',3,cldimn)
          ivarid = kfcdfaddvar(incid,'GDEPW',4013,'real',3,cldimn)
          ivarid = kfcdfaddvar(incid,'GDEPU',4111,'real',3,cldimn)
          ivarid = kfcdfaddvar(incid,'GDEPV',4112,'real',3,cldimn)
          ivarid = kfcdfaddvar(incid,'GDEPF',4113,'real',3,cldimn)
          CALL cdfputfltval(incid,'GDEPT',-1,GDEPT4,lflag_sync)
          CALL cdfputfltval(incid,'GDEPW',-1,GDEPW4,lflag_sync)
          CALL cdfputfltval(incid,'GDEPU',-1,GDEPU4,lflag_sync)
          CALL cdfputfltval(incid,'GDEPV',-1,GDEPV4,lflag_sync)
          CALL cdfputfltval(incid,'GDEPF',-1,GDEPF4,lflag_sync)
#else
          cldimn(1) = 'JPK'
          ivarid = kfcdfaddvar(incid,'GDEPT',4012,'real',1,cldimn)
          ivarid = kfcdfaddvar(incid,'GDEPW',4013,'real',1,cldimn)
          CALL cdfputfltval(incid,'GDEPT',-1,GDEPT4,lflag_sync)
          CALL cdfputfltval(incid,'GDEPW',-1,GDEPW4,lflag_sync)
#endif
          cldimn(1) = 'JPI'
          cldimn(2) = 'JPJ'
          ivarid = kfcdfaddvar(incid,'MBATHY',4003,'integer',2,cldimn)
          CALL cdfputintval(incid,'MBATHY',-1,mbathy,lflag_sync)

C
C 3. information of user:
C -----------------------
C
          clbuff='Les grilles sont ajoutees au fichier '
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('info','onc_defstdgrd',clbuff)
C
      ENDIF
      END
