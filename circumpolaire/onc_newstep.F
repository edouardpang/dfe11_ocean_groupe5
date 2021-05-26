CDIR$ LIST
      SUBROUTINE onc_newstep(cdncfil,kstep,ccurrent)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE onc_newstep
CCC                       *******************
CCC file:%M% (%F%) version:%I% %R% %S% modified at %H% %T%
CCC
CCC  PURPOSE :
CCC  ---------
CCC   Write a new record (in time and date variable) in the netcdf file
CCC   containing the physical variables. Check if all record variable
CCC   was written.
CC   METHOD :
CC   -------
CC
CC   USAGE :
CC   -------
CC   CALL onc_newstep(cdncfil,kstep)
CC        
CC
CC   INPUT :
CC   -----
CC   kstep [INTEGER]: the new record time step.
CC   OUTPUT :
CC   ------
CC
CC   EXTERNAL :
CC   -------- 
CC
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : creation at 20-11-97 16:32:50 by thierry philipovitch
CC      02-02-98, by T. Philipovitch: test: kt must be <= to NITEND
CC----------------------------------------------------------------------
      IMPLICIT NONE
CC parameters and commons
CC ======================
CDIR$ NOLIST
# include "netcdf.h"
# include "cdfpack.h"
CDIR$ LIST
      INTEGER kstep
      CHARACTER*(*) cdncfil
      CHARACTER*4 ccurrent
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER                 itimid,idateid,jrec
      CHARACTER*256           clbuff
      CHARACTER*20            clbuf
      CHARACTER               cldatab(jpdatlg)
      INTEGER                 inctest,jc,istart(2),icount(2)
      INTEGER                 jvaridx,idatidx,itimidx
      LOGICAL                 llfound
      REAL                    zsecond,zskip
      CHARACTER*_LG_VAR_NAME_ clvarname
      INTEGER                 incnum,incid
      LOGICAL                 llfisrtrec
CC----------------------------------------------------------------------
CC external functions
CC ==================
      LOGICAL       lfvardefined
      EXTERNAL      lfvardefined
      CHARACTER*10  cfitostr
      EXTERNAL      cfitostr
      CHARACTER*128 cfflttostr
      EXTERNAL      cfflttostr
CCC---------------------------------------------------------------------
CCC  package:cdfpack module:onc_newstep.F, CETIIS 97
CCC---------------------------------------------------------------------
C
C
C 0.getting the ID of the netCDF file:
C --------------------------------------
      CALL onc_chkinit('onc_newstep')
      CALL onc_setncf(cdncfil,incnum,incid)
      llfisrtrec = ncurstep(incnum) .EQ. -1

# if defined(opaCDF_debug)
C     divers messages pour debuggage:
      clbuff='Debut de l''enregistrement du  pas de temps '
      CALL concattxt(clbuff,cfitostr(kstep))
      CALL concattxt(clbuff,'dans le fichier ')
      CALL concattxt(clbuff,cdncfil)
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Numero interne du fichier: '
      CALL concattxt(clbuff,cfitostr(incnum))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='ID de ce fichier ='
      CALL concattxt(clbuff,cfitostr(incid))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Pas de temps actuel: '
      CALL concattxt(clbuff,cfitostr(ncurstep(incnum)))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Nombre de variables dependant du temps: '
      CALL concattxt(clbuff,cfitostr(nbphysrec(incnum)))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='ID de la variable de reference: '
      CALL concattxt(clbuff,
     $     cfitostr(nphysrecid(incnum,1)))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Nom de la variable de reference: '
      inctest = nf_inq_varname(incid,
     $     nphysrecid(incnum,1),clvarname)
      CALL concattxt(clbuff,clvarname)
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Nombre d''enregistrement de la variable de reference: '
      CALL concattxt(clbuff,cfitostr(nphysrec(incnum,1)))
      CALL procerr('debug','onc_newstep',clbuff)

      clbuff='Nombre d enregistrement des variables dependant du temps:'
      CALL procerr('debug','onc_newstep',clbuff)

      DO jrec = 1,nbphysrec(incnum)
         clbuff='  -->'
         inctest=nf_inq_varname(incid,nphysrecid(incnum,jrec),clvarname)
         CALL concattxt(clbuff,clvarname)
         CALL concattxt(clbuff,'(record=')
         CALL concattxt(clbuff,cfitostr(jrec))
         CALL concattxt(clbuff,',ID=')
         CALL concattxt(clbuff,cfitostr(nphysrecid(incnum,jrec)))
         CALL concattxt(clbuff,'): ')
         CALL concattxt(clbuff,cfitostr(nphysrec(incnum,jrec)))
         CALL procerr('debug','onc_newstep',clbuff)
      ENDDO

      IF (llfisrtrec) THEN
         CALL procerr('debug','onc_newstep',
     $        'C''est le premier enregistrement...')
      ENDIF
      
      IF (lstdout(incnum)) THEN
         CALL procerr('debug','onc_newstep',
     $   'Les variables standards seront copiees automatiquement...')
      ELSE
         CALL procerr('debug','onc_newstep',
     $        'Les variables standards ne seront pas copiees...')
      ENDIF
# endif

C 1.checking arguments and misc. checks:
C --------------------------------------
C     1.1 checking given step value:
      IF (kstep .LT. (nstpbegin-1)) THEN
          clbuff = 'Le pas de temps '
          CALL concattxt(clbuff,cfitostr(kstep))
          CALL concattxt(clbuff,' pour le fichier ')
          CALL concattxt(clbuff,cdncfil)
          CALL concattxt(clbuff,
     $         ' doit etre plus grand que le pas de temps=')
          CALL concattxt(clbuff,cfitostr(nstpbegin))
          CALL procerr('fatal_error_misc','onc_newstep',clbuff)
      ENDIF
      IF (kstep .EQ. (nstpbegin-1)) THEN
          clbuff = 'Le pas de temps '
          CALL concattxt(clbuff,cfitostr(kstep))
          CALL concattxt(clbuff,' pour le fichier ')
          CALL concattxt(clbuff,cdncfil)
          CALL concattxt(clbuff,
     $         ' doit etre plus grand que le pas de temps=')
          CALL concattxt(clbuff,cfitostr(nstpbegin))
          CALL procerr('warning','onc_newstep',clbuff)
      ENDIF
      IF (kstep .GT. nstpend) THEN
          clbuff = 'Le pas de temps '
          CALL concattxt(clbuff,cfitostr(kstep))
          CALL concattxt(clbuff,cfitostr(kstep))
          CALL concattxt(clbuff,' pour le fichier ')
          CALL concattxt(clbuff,cdncfil)
          CALL concattxt(clbuff,
     $         ' est plus grd que le dernier pdt ')
          CALL concattxt(clbuff,cfitostr(nstpend))
          CALL procerr('cdf_warning','onc_newstep',clbuff)
      ENDIF
C     ----------------------------------------------------------------------
C     debut du premier pas de temps:
C     CONVENTIONS:
C     1. nopadate0 est la date du debut du run (origine absolue)
C     2. Le pas de temps initial est le nombre de pdt debut la date d origine.
C     ----------------------------------------------------------------------

      IF (llfisrtrec) THEN
          CALL opadate2str_2(nopadate0,0,0,0,ccurdate(incnum))
          ncurstep(incnum) = kstep
          zsecond = fstpsec * REAL(kstep)
          zskip   = zsecond

# if defined(opaCDF_debug)
          clbuff='(premier pdt) '
          CALL concattxt(clbuff,
     $         'Date de reference dans OPA en numerique: ')
          CALL concattxt(clbuff,cfitostr(nopadate0))
          clbuff='(premier pdt) '
          CALL concattxt(clbuff,
     $         'Date de reference dans OPA en caracteres: ')
          CALL concattxt(clbuff,ccurdate(incnum))
          CALL procerr('debug','onc_newstep',clbuff)
# endif

C
C     Les autres pas de temps:
C
      ELSE
          IF (kstep .LE. ncurstep(incnum)) THEN
              clbuff = 'Le pas de temps '
              CALL concattxt(clbuff,cfitostr(kstep))
              CALL concattxt(clbuff,' dans le fichier ')
              CALL concattxt(clbuff,cdncfil)
              CALL concattxt(clbuff,' est deja enregistre.')
              CALL concattxt(clbuff,' Le dernier pas de temps etait ')
              CALL concattxt(clbuff,cfitostr(ncurstep(incnum)))
              CALL procerr('bad_arg','onc_newstep',clbuff)
          ENDIF
          zsecond = fstpsec * REAL(kstep)
          zskip   = fstpsec * REAL(kstep - ncurstep(incnum))
          ncurstep(incnum) = kstep
      ENDIF

# if defined(opaCDF_debug)
          clbuff='Nouveau pas de temps en secondes: '
          CALL concattxt(clbuff,cfflttostr(zsecond))
          CALL procerr('debug','onc_newstep',clbuff)

          clbuff='Saut par rapport au dernier pas de temps: '
          CALL concattxt(clbuff,cfflttostr(zskip))
          CALL procerr('debug','onc_newstep',clbuff)
# endif

C     1.2 checking if time and date was defined:
      IF (.NOT. lfvardefined(incid,'time',itimid)) THEN
          clbuff = 'La variable '
          CALL concattxt(clbuff,'time')
          CALL concattxt(clbuff,'n''est pas definie dans ')
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('cdf_fatal_misc','onc_newstep',clbuff)
      ENDIF
C
      IF (.NOT. lfvardefined(incid,'date',idateid)) THEN
          clbuff = 'La variable '
          CALL concattxt(clbuff,'date')
          CALL concattxt(clbuff,'n''est pas definie dans ')
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('cdf_fatal_misc','onc_newstep',clbuff)
      ENDIF

C     Recherche des numeros internes du temps et de la date:
      llfound = .FALSE.
      DO jvaridx=1,nbphysrec(incnum)
        llfound = ( nphysrecid(incnum,jvaridx) .EQ. itimid )
        IF (llfound) GOTO 10
      ENDDO
 10   IF (llfound) THEN
          itimidx = jvaridx
      ELSE
          clbuff = 'La variable '
          CALL concattxt(clbuff,'time')
          CALL concattxt(clbuff,'n''est pas enregistree pour ')
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('cdf_fatal_misc','onc_newstep',clbuff)
      ENDIF
C
      llfound = .FALSE.
      DO jvaridx=1,nbphysrec(incnum)
        llfound = ( nphysrecid(incnum,jvaridx) .EQ. idateid )
        IF (llfound) GOTO 20
      ENDDO
 20   IF (llfound) THEN
          idatidx = jvaridx
      ELSE
          clbuff = 'La variable '
          CALL concattxt(clbuff,'date')
          CALL concattxt(clbuff,'n''est pas enregistree pour ')
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('cdf_fatal_misc','onc_newstep',clbuff)
      ENDIF

C     1.3 checking if reference record was defined:
      IF (nbphysrec(incnum) .LT. 0) THEN
          clbuff='La variable de reference n''est pas definie dans '
          CALL concattxt(clbuff,cdncfil)
          CALL procerr('cdf_fatal_misc','onc_newstep',clbuff)
      ELSE IF (llfisrtrec) THEN
          inctest = nf_inq_varname(incid
     $         ,nphysrecid(incnum,1),clvarname)
          CALL cdfcheckmes(inctest,'onc_newstep','nf_inq_varname')
          clbuff='La variable de reference dans '
          CALL concattxt(clbuff,cdncfil)
          CALL concattxt(clbuff,' est ')
          CALL concattxt(clbuff,clvarname)
          CALL procerr('info','onc_newstep',clbuff)
      ENDIF
C     1.4 checking if all records were done:
      DO jrec = 1,nbphysrec(incnum)
        IF (nphysrec(incnum,jrec) .LT. 
     $       nphysrec(incnum,1)) THEN
           inctest = nf_inq_varname(incid
     $          ,nphysrecid(incnum,jrec),clvarname)
           CALL cdfcheckmes(inctest,'onc_newstep','nf_inq_varname')
           clbuff = 'La variable ' 
           CALL concattxt(clbuff,clvarname)
           CALL concattxt(clbuff,
     $     ' n''a pas ete enregistree (nrec=')
           CALL concattxt(clbuff,cfitostr(nphysrec(incnum,jrec)))
           CALL concattxt(clbuff,' ref.rec=')
           CALL concattxt(clbuff
     $          ,cfitostr(nphysrec(incnum,1)))
           CALL concattxt(clbuff,')')
           CALL concattxt(clbuff,' dans le fichier ')
           CALL concattxt(clbuff,cdncfil)
           CALL procerr('warning','onc_newstep',clbuff)
        ENDIF
      ENDDO
C
C 2.Validation of the record:
C ---------------------------
C
C     2.0 incrementation of the date:
      CALL addtim2date(ccurdate(incnum),zskip,ntypyear,ccurdate(incnum))

# if defined(opaCDF_debug)
      clbuff='L''incrementation de la date donne:'
      CALL concattxt(clbuff,ccurdate(incnum))
      CALL procerr('debug','onc_newstep',clbuff)
# endif

      DO jc=1,jpdatlg
        cldatab(jc) = ccurdate(incnum)(jc:jc)
      ENDDO
C     2.1 writing the date:
      IF (llfisrtrec) THEN
          nphysrec(incnum,idatidx) = 1
      ELSE
          nphysrec(incnum,idatidx) = nphysrec(incnum,idatidx) + 1
      ENDIF
      istart(1) = 1
      istart(2) = nphysrec(incnum,idatidx)
      icount(1) = jpdatlg
      icount(2) = 1

# if defined(opaCDF_debug)
      clbuff='Enregistrement de la variable DATE:'
      CALL procerr('debug','onc_newstep',clbuff)
      clbuff='  istart=['
      CALL concattxt(clbuff,cfitostr(istart(1)))
      CALL concattxt(clbuff,',')
      CALL concattxt(clbuff,cfitostr(istart(2)))
      CALL concattxt(clbuff,']')
      CALL procerr('debug','onc_newstep',clbuff)
      clbuff='  icount=['
      CALL concattxt(clbuff,cfitostr(icount(1)))
      CALL concattxt(clbuff,',')
      CALL concattxt(clbuff,cfitostr(icount(2)))
      CALL concattxt(clbuff,']')
      CALL procerr('debug','onc_newstep',clbuff)
# endif

      inctest = nf_put_vara_text(incid,idateid,istart,icount,cldatab)
      CALL cdfcheckmes(inctest,'onc_newstep','nf_put_vara_text')
C     2.2 writing time in second:
      IF (llfisrtrec) THEN
          nphysrec(incnum,itimidx) = 1
      ELSE
          nphysrec(incnum,itimidx) = nphysrec(incnum,itimidx) + 1
      ENDIF
      istart(1) = nphysrec(incnum,itimidx)

# if defined(opaCDF_debug)
      clbuff='Enregistrement de la variable TIME:'
      CALL procerr('debug','onc_newstep',clbuff)
      clbuff='  istart='
      CALL concattxt(clbuff,cfitostr(istart(1)))
      CALL procerr('debug','onc_newstep',clbuff)
      clbuff='  valeur (en sec.)='
      CALL concattxt(clbuff,cfflttostr(zsecond))
      CALL procerr('debug','onc_newstep',clbuff)
# endif

      inctest   = nf_put_var1_real(incid,itimid,istart,zsecond)
      CALL cdfcheckmes(inctest,'onc_newstep','nf_put_var1_real')

C
C 2. Synchronization:
C -------------------
      IF (lflag_sync) THEN
         inctest = nf_sync(incid)
         CALL cdfcheckmes(inctest,'onc_newstep','nf_sync')
      ENDIF
      clbuff = 'Le fichier '
      CALL concattxt(clbuff,cdncfil)
      CALL concattxt(clbuff,
     $     ' est pret pour recevoir les enregistrements')
      CALL concattxt(clbuff,' du pas de temps: ')
      CALL concattxt(clbuff,cfitostr(ncurstep(incnum)))
      CALL concattxt(clbuff,' et de la date:')
      CALL concattxt(clbuff,ccurdate(incnum))
      WRITE(clbuf,'(A13,A4,A3)')'onc_newstep (',ccurrent,' s)' 
      !CALL procerr('info','onc_newstep',clbuff)
      CALL procerr('info',clbuf,clbuff)
C
C 4. writing standart values:
C ---------------------------
C
      IF (lstdout(incnum)) THEN
         CALL onc_putstdvar(cdncfil)
      ENDIF
      END
