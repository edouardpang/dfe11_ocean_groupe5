!
! Attention  !!
! common.h est utilisé par des routines en F77 et F90
! Il faut garder la syntaxe F77 des declarations des tableaux
! On peut dépasser 72 caractères par lignes avec la bonne directive de compilation (-ffixed-line-length-132)
! Mais on garde les 6 espaces en début de lignes
!
      REAL glamt(jpi,jpj),glamu(jpi,jpj),glamv(jpi,jpj),glamf(jpi,jpj)
      REAL gphit(jpi,jpj),gphiu(jpi,jpj),gphiv(jpi,jpj),gphif(jpi,jpj)
      REAL gdept(jpk),gdepw(jpk)
      INTEGER mbathy(jpi,jpj)
      COMMON /grid_array/glamt,glamu,glamv,glamf,gphit,gphiu,gphiv,gphif,gdept,gdepw,mbathy
!
      REAL e1t(jpi,jpj),e1u(jpi,jpj),e1v(jpi,jpj),e1f(jpi,jpj)
      REAL e2t(jpi,jpj),e2u(jpi,jpj),e2v(jpi,jpj),e2f(jpi,jpj)
      COMMON /array1/ e1t,e1u,e1v,e1f, e2t,e2u,e2v,e2f
!
      REAL f(jpi,jpj) 
      REAL prof(jpi,jpj,jpk)
      REAL profu(jpi,jpj,jpk),profv(jpi,jpj,jpk),proff(jpi,jpj,jpk)
      COMMON /array5/ prof,profu,profv,proff,f
!
      REAL htot(jpi,jpj),htotu(jpi,jpj),htotv(jpi,jpj),htotf(jpi,jpj)
      COMMON /array4/ htot,htotu,htotv,htotf
!
      REAL hmask(jpi,jpj,jpk),umask(jpi,jpj,jpk),vmask(jpi,jpj,jpk)
      REAL fmask(jpi,jpj,jpk),tmask(jpi,jpj,jpk)
      COMMON /array3/ hmask,umask,vmask,tmask,fmask
!
      REAL hb(jpi,jpj,jpk),hn(jpi,jpj,jpk),ha(jpi,jpj,jpk)
      REAL ub(jpi,jpj,jpk),un(jpi,jpj,jpk),ua(jpi,jpj,jpk)
      REAL vb(jpi,jpj,jpk),vn(jpi,jpj,jpk),va(jpi,jpj,jpk) 
      COMMON /array2/ hb,hn,ha,ub,un,ua,vb,vn,va
! 
! Local arrays used by step and diag (not in common block)
      REAL zkh(jpi,jpj,jpk),zku(jpi,jpj,jpk),zkv(jpi,jpj,jpk)
      REAL zkc(jpi,jpj,jpk)
      REAL zkf(jpi,jpj,jpk)
!
! REAL*4 array for writing netcdf file (lfloat in ibcdfpack are REAL*4, ...)
! These declarations are not overwritten by -fdefault-real-8 compilation directive
      REAL*4 hcdf(jpi,jpj,jpk),ucdf(jpi,jpj,jpk),vcdf(jpi,jpj,jpk)
      REAL*4 pcdf(jpi,jpj,jpk),rcdf(jpi,jpj,jpk),wkk(jpi,jpj,jpk)
      REAL*4 glamt4(jpi,jpj),glamu4(jpi,jpj),glamv4(jpi,jpj),glamf4(jpi,jpj)
      REAL*4 gphit4(jpi,jpj),gphiu4(jpi,jpj),gphiv4(jpi,jpj),gphif4(jpi,jpj)
      REAL*4 gdept4(jpk),gdepw4(jpk)
      COMMON /grid_array4/glamt4,glamu4,glamv4,glamf4,gphit4,gphiu4,gphiv4,gphif4,gdept4,gdepw4
!
      REAL uflux(jpi,jpj),vflux(jpi,jpj)
      REAL potvor(jpi,jpj),coriou(jpi,jpj),coriov(jpi,jpj)
      REAL rotn(jpi,jpj,jpk),rotb(jpi,jpj,jpk)
      REAL taux(jpi,jpj),tauy(jpi,jpj)
      REAL gu(jpi,jpj),gv(jpi,jpj),cur(jpi,jpj)
      COMMON /array6/ taux,tauy,gu,gv,cur, potvor,rotn,rotb
!
!
#if diaTrends
      REAL utrd(jpi,jpj,jpk,6), vtrd(jpi,jpj,jpk,6), htrd(jpi,jpj,jpk,3)
      COMMON /array8/ utrd, vtrd, htrd
#endif
! 
      REAL g,dt,dt2,zdt,gamm
      REAL xnu,xnuh,dx,dy, xlamda
      INTEGER iiperio,ijperio, islip, itau_stop
      COMMON /shallow_var/dx,dy,dt,dt2,gamm,xnu,xnuh,g,xlamda,iiperio,ijperio,islip,itau_stop
!
      REAL relaxh(jpi,jpj), relaxu(jpi,jpj), hinit(jpi,jpj)
      COMMON /relaxation/relaxh, relaxu, hinit 
!
      INTEGER il0
      REAL h0(jpj)
      REAL xu0,xl0,u0(jpj),v0(jpj)
      COMMON /initial/h0,u0,v0,xu0,xl0,il0
!
      INTEGER ndate0,nit000,nitend,nbisex,nwrite
      REAL*4 raajj,raamo,rjjhh,rhhmm,rmmss
      REAL*4 xmiss,rdt
      COMMON /for_opacdf/rdt,ndate0,nit000,nitend,nbisex,nwrite,raajj,raamo,rjjhh,rhhmm,rmmss,xmiss
      INTEGER ji,jj,jk,jt,kt
!
      REAL ff(jpi-2,2001)
      REAL f1(jpi-2,2001)
      REAL xl,xk,xd9,jcost,ros
      COMMON /for_imsl/ff,f1,xl,xk,xd9,jcost,ros
!
      REAL ctb(jpi,jpj,jpk,jptra),ctn(jpi,jpj,jpk,jptra),cta(jpi,jpj,jpk,jptra)
      COMMON /shallow_traceur/ctb,ctn,cta
  
