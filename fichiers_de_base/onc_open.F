CDIR$ LIST
      SUBROUTINE onc_open(kgeoref)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE onc_open
CCC                       ****************
CCC  file:%M% (%F%) version:%I% %R% %S% modified at %H% %T%
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
CC      26-10-97 11:00:51 by thierry philipovitch: creation 
CC      28-01-98 by thierry philipovitch: addition of georeferencing flag
CC      06-04-98 by thierry philipovitch: addition of synchro. flag
CC----------------------------------------------------------------------
CC    IMPLICIT NONE
CC parameters and commons
CC ======================
CDIR$ NOLIST
# include "parameter.h"
# include "common.h"
# include "netcdf.h"
# include "cdfpack.h"
# if defined(key_cpp_save)
# include "cppkeys.h"
# endif
      INTEGER kgeoref
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER jvar,jncfile
# if defined(key_cpp_save)
      INTEGER jkey
# endif
      INTEGER      itest
      INTEGER      ndxncref
      CHARACTER*20 clenvmnt
      CHARACTER*80 clbuff
CC----------------------------------------------------------------------
CC external functions
CC ===================
      INTEGER    kfgetenvvar
      LOGICAL    lfucstrncmp
CCC---------------------------------------------------------------------
CCC  package:cdfpack module:onc_open.F, CETIIS 97
CCC---------------------------------------------------------------------
C
      IF (ncdfinit .EQ. jpncinitflg) THEN
          clbuff='L''initialisation a deja ete faite...'
          CALL procerr('warning','onc_open',clbuff)
          RETURN
      ENDIF
      IF ((kgeoref .LT. 0) .OR. (kgeoref .GT. 2)) THEN
          clbuff='L''argument de onc_open doit etre 0,1,2'
          CALL procerr('bad_arg','onc_open',clbuff)
      ENDIF
C     testing the version of the netCDF library used:
C     ----------------------------------------------
      clbuff='Version de la librairie netcdf: '
      CALL concattxt(clbuff,nf_inq_libvers())
      CALL procerr('info','onc_open',clbuff)
C
C     getting from OPA some parameters:
C     ---------------------------------
C
      IF (rdt .LE. 0.) THEN
          CALL procerr('bad_arg','onc_open'
     $         ,'Le pas de temps en seconde est <= 0')
      ENDIF
      nresx       = jpi
      nresy       = jpj
      nresz       = jpk
      fstpsec     = rdt
      nopadate0   = ndate0
      nstpbegin   = nit000
      nstpend     = nitend
      ntypyear    = nbisex
      IF (ntypyear .EQ. -1) THEN
         gdaypyear = 360.
         gmthpyear = 12.
      ELSE
         gdaypyear   = RAAJJ
         gmthpyear   = RAAMO
      ENDIF
      ghrpday     = RJJHH
      gmnphour    = RHHMM
      gsecpmin    = RMMSS
# if defined(key_cpp_save)
      CALL setupkeys()
      nbcppkeys   = jpnbkey
      IF (nbcppkeys .GT. 0) THEN
          copakeys = ccppkeys(1)
          DO jkey = 2 , nbcppkeys
            CALL concattxt(copakeys,':')
            CALL concattxt(copakeys,ccppkeys(jkey))
          ENDDO
      ELSE
          copakeys = 'none'
      ENDIF
# else
C Get cpp keys from the shell
      call getenv("OPA_CPP_KEYS",copakeys)
      nbcppkeys = 1
C
# endif
C
      ncid_grid   = -1
      cdfgridfile = ' '
      nbncfile    = 0
      ndxncref    = -1
      ngeoref     = kgeoref
C     setting defaut int/real missing values:
C
      DO jncfile=1,jpnbncmax
        ccurdate(jncfile) = '00-XXX-0000:00:00:00'
        ncurstep(jncfile) = -1
        ncidfile(jncfile) = -1
        nbphysrec(jncfile)= -1
        cdncfile(jncfile) = ' '
        DO jvar=1,jpnbvarmax
          nphysrecid(jncfile,jvar) = -1
          nphysrec(jncfile,jvar)   = -1
        ENDDO
      ENDDO
C
      CALL iofindcfgfile(cpepicdeff,cepicfile)
      clbuff='Fichier EPIC utilise: '
      CALL concattxt(clbuff,cepicfile)
      CALL procerr('info','onc_open',clbuff)
C
C     SYNCHRONIZATION MODE:
C     ---------------------
C
C     default value:
C     Getting some environment variable if any:
      itest = kfgetenvvar('OPACDF_SYNC',clenvmnt)
      IF (itest .GT. 0) THEN
         lflag_sync = lfucstrncmp(clenvmnt,'yes',3)
         IF (lflag_sync) THEN
            clenvmnt = 'yes'
         ELSE
            clenvmnt = 'no'
         ENDIF
      ELSE
         clenvmnt = 'no'
         lflag_sync = .FALSE.
      ENDIF
      clbuff = 'SYNCHRONISATION AUTOMATIQUE DES SORTIES NETCDF: '
      CALL concattxt(clbuff,clenvmnt)
      CALL procerr('info','onc_open',clbuff)
      ncdfinit    = jpncinitflg
      END
