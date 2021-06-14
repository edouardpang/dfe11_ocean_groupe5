      PROGRAM shallow_water 
!
#include "parameter.h"
#include "common.h"
!
      CHARACTER*32 cddims(4)
      CHARACTER(len=255) netcdf_file,diag_file,restart_file_out,restart_file_in
      CHARACTER(len=255) clbuff
      CHARACTER(len=4)ccurrent
      EXTERNAL fcost,fcos1
      REAL hmin, hmax
      REAL current, current0
!
!
! 1. Initializations
! ------------------
!
#if defined dynRestart
!
      CALL getenv("RESTART_FILE_IN", restart_file_in)
      OPEN(51,file = restart_file_in, form = 'unformatted')
      clbuff = 'Restart with file '
      CALL concattxt(clbuff,restart_file_in)
      CALL procerr('info','main',clbuff)
!
#endif
!
      CALL init
!
!
! 2. Open output files
! --------------------
!
      CALL getenv("DIAG_OUTPUT_FILE", diag_file)
      CALL getenv("RESTART_FILE_OUT", restart_file_out)
      OPEN(unit = 90, file = TRIM(diag_file))
!
! Open NECDF output
!
      CALL getenv("NETCDF_OUTPUT_FILE", netcdf_file)
      netcdf_file = TRIM(netcdf_file)
      CALL onc_open(2)
      CALL onc_opnncfile(netcdf_file)
      CALL onc_addtime(netcdf_file)

      CALL onc_mknewdim(netcdf_file,'JPI',jpi)
      CALL onc_mknewdim(netcdf_file,'JPJ',jpj)
      CALL onc_mknewdim(netcdf_file,'JPK',jpk) 
      CALL onc_defstdgrd(netcdf_file)
      cddims(1) = 'JPI'
      cddims(2) = 'JPJ'
      cddims(3) = 'JPK'
      cddims(4) = 'time'
      CALL onc_addvar(netcdf_file,'UN',3042,'FLOAT',4,cddims,'GLAMU','GPHIU','GDEPT',' ','no',-1)
      CALL onc_addvar(netcdf_file,'VN',3043,'FLOAT',4,cddims,'GLAMV','GPHIV','GDEPT',' ','no',-1)
      CALL onc_addvar(netcdf_file,'HN',3060,'FLOAT',4,cddims,'GLAMT','GPHIT','GDEPT',' ','no',-1)
      CALL onc_addvar(netcdf_file,'PV',3061,'FLOAT',4,cddims,'GLAMF','GPHIF','GDEPT',' ','no',-1)
      CALL onc_addvar(netcdf_file,'ROTN',6007,'FLOAT',4,cddims,'GLAMF','GPHIF','GDEPT',' ','no',-1)
      xmiss=-999999.
      CALL onc_misvarflt(netcdf_file,'UN',xmiss)
      CALL onc_misvarflt(netcdf_file,'VN',xmiss)
      CALL onc_misvarflt(netcdf_file,'HN',xmiss)
      CALL onc_misvarflt(netcdf_file,'PV',xmiss)
      CALL onc_misvarflt(netcdf_file,'ROTN',xmiss)
      DO jk = 1,jpk
        DO jj = 1,jpj
          DO ji = 1,jpi
            hcdf(ji,jj,jk)=hn(ji,jj,jk)
            ucdf(ji,jj,jk)=un(ji,jj,jk)
            vcdf(ji,jj,jk)=vn(ji,jj,jk)
            pcdf(ji,jj,jk)=potvor(ji,jj) 
            rcdf(ji,jj,jk)=rotn(ji,jj,jk) 
            IF(hmask(ji,jj,jk).EQ.0)hcdf(ji,jj,jk) = xmiss 
            IF(umask(ji,jj,jk).EQ.0)ucdf(ji,jj,jk) = xmiss
            IF(vmask(ji,jj,jk).EQ.0)vcdf(ji,jj,jk) = xmiss
            IF(fmask(ji,jj,jk).EQ.0)pcdf(ji,jj,jk) = xmiss
            IF(fmask(ji,jj,jk).EQ.0)rcdf(ji,jj,jk) = xmiss
          ENDDO
        ENDDO
      ENDDO
!
! Output initial state
!
          

      CALL onc_newstep(netcdf_file,nit000,'0')
      CALL onc_putvarflt(netcdf_file,'UN',ucdf,wkk)
      CALL onc_putvarflt(netcdf_file,'VN',vcdf,wkk)
      CALL onc_putvarflt(netcdf_file,'HN',hcdf,wkk)
      CALL onc_putvarflt(netcdf_file,'PV',pcdf,wkk)
      CALL onc_putvarflt(netcdf_file,'ROTN',rcdf,wkk)
!
!
! 3. Temporal loop
! ----------------
!
      call cpu_time(current0)
!
      DO jt = nit000,nitend
!
! Advance one time step
!
        CALL step(jt)
          

!
! Output to file
!
        IF (mod(jt,nwrite).EQ.0) THEN
            call cpu_time(current)
            current = current - current0
            WRITE(ccurrent,'(i4)')int(current)
!
            DO jk = 1,jpk
              DO jj = 1,jpj
               DO ji = 1,jpi
                 hcdf(ji,jj,jk)=hn(ji,jj,jk)
                 ucdf(ji,jj,jk)=un(ji,jj,jk)
                 vcdf(ji,jj,jk)=vn(ji,jj,jk)
                 pcdf(ji,jj,jk)=potvor(ji,jj) 
                 rcdf(ji,jj,jk)=rotn(ji,jj,jk) 
                 IF(hmask(ji,jj,jk).EQ.0)hcdf(ji,jj,jk) = xmiss 
                 IF(umask(ji,jj,jk).EQ.0)ucdf(ji,jj,jk) = xmiss
                 IF(vmask(ji,jj,jk).EQ.0)vcdf(ji,jj,jk) = xmiss
                 IF(fmask(ji,jj,jk).EQ.0)pcdf(ji,jj,jk) = xmiss
                 IF(fmask(ji,jj,jk).EQ.0)rcdf(ji,jj,jk) = xmiss
               ENDDO
             ENDDO
           ENDDO
           CALL onc_newstep(netcdf_file,jt,ccurrent)
           CALL onc_putvarflt(netcdf_file,'UN',ucdf,wkk)
           CALL onc_putvarflt(netcdf_file,'VN',vcdf,wkk)
           CALL onc_putvarflt(netcdf_file,'HN',hcdf,wkk)
           CALL onc_putvarflt(netcdf_file,'PV',pcdf,wkk)
           CALL onc_putvarflt(netcdf_file,'ROTN',rcdf,wkk)
        ENDIF
!
        CALL minmax(un(1,1,ik), umask(1,1,ik), hmin,hmax)
        hmax = AMAX1(ABS(hmax),ABS(hmin))
        IF ( hmax .GT. 1000.) THEN
            CALL onc_close()
            PRINT *, 'Le modele diverge au pas de temps ',jt
            CALL EXIT(50)
        ENDIF
!
      ENDDO
!
!
! 5. Close files and additionnal diagnostics
! ------------------------------------------
!
      OPEN(unit = 50, file = TRIM(restart_file_out), form='UNFORMATTED')
      write(50)un,ub,vn,vb,hn,hb,rotb,rotn,jt-1
      CLOSE(90)
      CLOSE(50)
      CALL onc_close()
!
      END
