c     ******************************************
c     O significado de cada variável dinâmica está 
c     listado no subprograma alocadin.for
c     *******************************************

      !>Stores the frame structure dynamic variables

      MODULE BIG_STUFF
      integer*4,allocatable:: IRN(:),ICN(:),IW(:),IKEEP(:,:),IW1(:)
     $     ,ic(:,:),ko(:),ko1(:),inct(:,:),kn(:),im(:,:),itip(:,:)
     $     ,indx(:),in(:,:),irna(:),icna(:),glgl(:,:) ,itr(:,:),ngl(:)
     $     ,iglv(:,:),notl(:),ipar_nduplo(:,:),iseq_nosdp(:,:)
     $     ,iconect_nosdp(:,:),iresistcrit(:)

      REAL*8,ALLOCATABLE::A(:),w27(:),p0(:),p0i(:),rkl(:,:,:),Aag(:)
     $     ,re(:),ri(:),rc(:),resp(:),ra(:),f(:),aux(:),vv(:),p(:),df(:)
     $     ,rkls(:,:),fi(:),dfi(:,:),dp(:),v(:),h(:,:),vs(:),as(:)
     $     ,as1(:), rhook(:,:,:),sx(:,:),sy(:,:),sxy(:,:),ro(:),ram(:)
     $     ,rs(:),qs(:), qsi(:),w(:),vint(:),vtl(:,:),rg(:),qsi3(:)
     $     ,w3(:),reacao(:), rl1(:),rl2(:),xs1(:),xs2(:),pl0(:,:,:),pl(:
     $     ,:,:),qsl(:,:,:), drl1(:),drl2(:),rh(:),rf(:),rp(:),RHS(:),
     $     drex(:,:),drey(:,:),drez(:,:),drexy(:,:),drexz(:,:),dreyz(:
     $     ,:), p1(:,:,:),p2(:,:,:),p3(:,:,:),daf(:,:,:,:),dch(:,:,:,:),
     $     aurxy(:,:),aur1xy(:,:),aurxz(:,:),aur1xz(:,:),rsl(:,:,:)
     $     ,auryz(:,:),aur1yz(:,:),valor(:,:),fno(:),rmp(:,:,:),rmg(:,:
     $     ,:), fri(:),frib(:),rmfi(:,:),fd(:,:),brd(:,:),yb(:),hb(:,:)
     $     ,rib(:) ,qx(:,:),qy(:,:),qn(:,:),qt(:,:),rmz(:,:),rvy(:,:)
     $     ,rnx(:,:),rlambdalag(:),vlambda(:),EsfResist(:,:),vs_aux(:),
     $     vs_ant(:)

c     *************************************************************
c     para plasticidade e temperatura temos as seguintes variáveis
c     seus significados também estão no subprograma alocadina
c     *************************************************************
      REAL*8,ALLOCATABLE::prsi0(:,:,:,:),ep(:,:,:,:,:),eeq(:,:,:,:)
     $     ,defl(:,:),rhdg(:,:),rhvg(:,:),temc(:),temb(:),dtemc(:)
     $     ,dtemb(:) ,calt(:)

      
      CHARACTER*32 entrada      !nome do arquivo de entrada
      
      END MODULE BIG_STUFF

      !>Stores the frame structure static variables

      MODULE BLOCLAGATBIG

      INTEGER*4 nnos,nel,nnr,nnc,n,ig1,ig2,nfc,ipt ,k2,li,lj,ig3
     $     ,neltr,kkk,nzexp,nzexp1,kkka,nprint,nnosdp,nnosmt
      INTEGER*4 nelcs,ngs,jlado,idv,ig,ia,it,ih1 ,igf,nelcar,nelnco
     $     ,iposic,nnoel

      INTEGER*4 npt,ngf,ipc,j,noel,ng,ng3,ngf1,not
      REAL*8 rnorma,tol,rmrd,rvrd,rnrd,rmsd,rnsd,rvsd

      REAL*8  t,t1,t2,ba,bb,bc,bd,be,bf,bg,bh,bri,bj,bk ,ca,cb,cc,cd,ce
     $     ,cf,cg,crh,cri,cj,ck
      REAL*8  dxg1,dxg2,dxp1,dxp2,dyg1,dyg2
      REAL*8  dyp1,dyp2,rjac,rjac0,salfa,calfa,tf1,tf2,tfa,tfb
      REAL*8  pi,teta,rjac2,a1
      
      REAL*8 re1,re2,ren,rex,rey,rexy,rez,rexz,reyz

      REAL*8 rr,alfa,bea,gam,alfa1,alfa2,beta1,beta2 ,gama1 ,gama2
     $     ,dgaxydx,dgaxydy,dgaxzdx,dgaxzdz,dgayzdy, dgayzdz ,d2gxydxdy
     $     ,d2gxydydx,d2gxydxdx,d2gxydydy, d2gxzdxdz,d2gxzdzdx
     $     ,d2gxzdxdx,d2gxydzdz, d2gyzdydz,d2gyzdzdy,d2gyzdydy,d2gyzdzdz


      REAL*8 velo,rbn,rgn,dt,xc,yc,xp,yp,rimp


      REAL*8 dpl(3,2),xsi1,xsi2,xsi3,qi3
     

      REAL*8 pxsi(3),peta(3),pdel(3),a0(2,2),af(2,2),
     $     a0inv(2,2),aaux(3,3),ch(2,2),up(3),vp(3),wp(3),aup(3),
     $     aaux1(3,3)
      REAL*8  rco1,rco2,aup2(3),rco3,drco1x,drco2x,drco1y,drco2y,d2ue(6
     $     ,6),aval(2),acp(2,2),dcp(3,3),drco3x,drco3y,rjac1,marca(17)
     $     ,nordem(17)

c      common /bloco23/ due(6),daf1(3,3),d2af(3,3),
c     #aup1(3),rlex,rley,rlez,rlexy,aup12(3)
c     #,rlexz,rleyz,ham(7,3),wh(7)

      REAL*8  due(6),daf1(3,3),d2ch(3,3,6,6),d2af(3,3),aup1(3),rlex
     $     ,rley,rlez,rlexy,aup12(3),ham(12,3),wh(12) ,rlexz,rleyz,pp1(3
     $     ,6,6),pp2(3,6,6),pp3(3,6,6),d2rex(6,6), d2rey(6,6),d2rez(6,6)
     $     ,d2rexy(6,6),d2rexz(6,6),d2reyz(6,6)



      INTEGER*4 ICNTL(30),INFO(20),LIW,LA,IFLAG,I,NSTEPS,MAXFRT
      REAL*8 CNTL(5),OPS,kac
      
      INTEGER*4 n27,nz

 
c     blocos para a plasticidade e temperatura
      INTEGER*4  np,ip,iaret,nt,ntem
      REAL*8 s(6),rl(6),rhret(6),rt(6),ceret(6,6),cp(6,6),dpret(6,6)
     $     ,rpret(6,6),s0,tpos,ds(6),rps(6),rph(6),rpt(6),ee(6),de(6,6)
     $     ,et(6),rpe(6),et1(6),rhd(2),rhv(2)

      END MODULE BLOCLAGATBIG
