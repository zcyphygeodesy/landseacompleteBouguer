      subroutine LTerAllBLH(BLH,dtm,sfh,nlat,nlon,hd,dr,GRS,ter)
      !dr-积分半径m，简化公式，严密性不变
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj,astat(5)
	real*8::dr,dtm(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst(7),ter(7)
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3),BLH2(3),XYZ2(3)
	real*8 CGrdPntD2,zch,L0,L1,L2,K1,K2,P1,P2,V1,V2,G1,G2,T1,T2,tfi,ff1,ff2
	real*8::GRS(6),ge,tp,dh,hh,tmp,fx
	real*8 rlat,rlon,rlat1,rlon1,sin2f,cos2f,sina,cosa,u!,cos2,cos3
!-----------------------------------------------------------------
      ge=6.67428d-11;tp=2.67d3!万有引力常数!地球平均半径!地形密度
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;gp=ge*tp
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      zch=CGrdPntD2(BLH(2),BLH(1),dtm,nlat,nlon,hd)
      ter=0.d0;mdr=r0*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/r0/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      rlat=rln(2)*RAD;rlon=rln(3)*RAD;hh=BLH(3)-BLH0(3);if(hh<1.d0)hh=0.d0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=sfh(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          if(L0>dr)goto 9101
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          if(L1<mdr)then!计算奇异积分
             call LTerAllsgn(BLH,dtm,sfh,nlat,nlon,hd,i,j,4,GRS,rst)
             ter=ter+rst; goto 9101 
          endif
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          sin2f=dsqrt((1.d0-tt)/2.d0);u=2.d0*sin2f*dsqrt((1.d0+tt)/2.d0)
          cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/u
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/u
          fx=dlog((L1+r1-rr*tt)/(L1-r1+rr*tt))/2.d0
          P1=(r1+3.d0*rr*tt)*L1/2.d0+rr**2/2.d0*(3.d0*tt**2-1.d0)*fx
          K1=(r1**2+3.d0*rr**2)*tt+rr*r1*(1.d0-6.d0*tt**2)
          K1=K1/L1-rr*(3.d0*tt**2-1.d0)*fx
          V1=(rr-r1*tt)/L1/u
          G1=r1/rr**2/L1+r1*(rr-r1*tt)/rr/L1**3
          BLH2=BLH1;BLH2(3)=sfh(i,j)+dtm(i,j)-zch;dh=dtm(i,j)-zch
          call BLH_XYZ(GRS,BLH2,XYZ2)
          L2=dsqrt((XYZ2(1)-XYZ(1))**2+(XYZ2(2)-XYZ(2))**2+(XYZ2(3)-XYZ(3))**2)
          r2=dsqrt(XYZ2(1)**2+XYZ2(2)**2+XYZ2(3)**2)
          fx=dlog((L2+r2-rr*tt)/(L2-r2+rr*tt))/2.d0
          P2=(r2+3.d0*rr*tt)*L2/2.d0+rr**2/2.d0*(3.d0*tt**2-1.d0)*fx
          K2=(r2**2+3.d0*rr**2)*tt+rr*r2*(1.d0-6.d0*tt**2)
          K2=K2/L2-rr*(3.d0*tt**2-1.d0)*fx
          V2=(rr-r2*tt)/L2/u
          G2=r2/rr**2/L2+r2*(rr-r2*tt)/rr/L2**3
          ter(1)=ter(1)+(P2-P1)*gp*ds/gr/r1**2
          ter(2)=ter(2)-(K2-K1)*gp*ds/r1**2
          if(dabs(u)>1.d-10)then
            ter(3)=ter(3)+(V2-V1)*cosa*gp*ds/gr/rr
            ter(4)=ter(4)+(V2-V1)*sina*gp*ds/gr/rr
          endif
          ter(5)=ter(5)+(G2-G1)*gp*ds
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
!--------------------------------------------------------------------------------
      subroutine LTerAllsgn(BLH,dtm,sfh,nlat,nlon,hd,i0,j0,m,GRS,rst)
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
!-------------------------------------------------------------
      implicit none
	integer::m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,dtm(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rv,lon,lat,rst(7)
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3),BLH2(3),XYZ2(3)
	real*8 CGrdPntD2,zch,zch1,L0,L1,L2,P1,P2,K1,K2,V1,V2,G1,G2,T1,T2,tfi,ff1,ff2
	real*8::GRS(6),ge,tp,dh,hh,fx
	real*8 rlat,rlon,rlat1,rlon1,sin2f,cos2f,sina,cosa,u,cos2,cos3
!-----------------------------------------------------------------
      ge=6.67428d-11;tp=2.67d3!万有引力常数!地球平均半径!地形密度
      rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0;gp=ge*tp
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      zch=CGrdPntD2(BLH(2),BLH(1),dtm,nlat,nlon,hd)
      rst=0.d0;mdr=r0*rv*RAD*dcos(rln(2)*RAD)/dble(m)/2.d0  !奇异点判断
      rlat=rln(2)*RAD;rlon=rln(3)*RAD;hh=BLH(3)-BLH0(3)
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),sfh,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          zch1=CGrdPntD2(BLH1(2),BLH1(1),dtm,nlat,nlon,hd)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=rv**2*RAD**2*dcos(rln1(2)*RAD)*r1**2
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          if(L1<mdr)L1=mdr
          sin2f=dsqrt((1.d0-tt)/2.d0);u=2.d0*sin2f*dsqrt((1.d0+tt)/2.d0)
          if(dabs(u)<1.d-18)u=1.d-18
          cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/u
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/u
          fx=dlog((L1+r1-rr*tt)/(L1-r1+rr*tt))/2.d0
          P1=(r1+3.d0*rr*tt)*L1/2.d0+rr**2/2.d0*(3.d0*tt**2-1.d0)*fx
          K1=(r1**2+3.d0*rr**2)*tt+rr*r1*(1.d0-6.d0*tt**2)
          K1=K1/L1-rr*(3.d0*tt**2-1.d0)*fx
          V1=(rr-r1*tt)/L1/u
          G1=r1/rr**2/L1+r1*(rr-r1*tt)/rr/L1**3
          BLH2=BLH1;BLH2(3)=BLH1(3)+zch1-zch
          call BLH_XYZ(GRS,BLH2,XYZ2)
          L2=dsqrt((XYZ2(1)-XYZ(1))**2+(XYZ2(2)-XYZ(2))**2+(XYZ2(3)-XYZ(3))**2)
          r2=dsqrt(XYZ2(1)**2+XYZ2(2)**2+XYZ2(3)**2)
          if(L2<mdr)L2=mdr;dh=zch1-zch
          fx=dlog((L2+r2-rr*tt)/(L2-r2+rr*tt))/2.d0
          P2=(r2+3.d0*rr*tt)*L2/2.d0+rr**2/2.d0*(3.d0*tt**2-1.d0)*fx
          K2=(r2**2+3.d0*rr**2)*tt+rr*r2*(1.d0-6.d0*tt**2)
          K2=K2/L2-rr*(3.d0*tt**2-1.d0)*fx
          V2=(rr-r2*tt)/L2/u
          G2=r2/rr**2/L2+r2*(rr-r2*tt)/rr/L2**3
          if(dabs(P2-P1)>1.d-10)rst(1)=rst(1)+(P2-P1)*gp*ds/gr/r1**2
          if(dabs(K2-K1)>1.d-10)rst(2)=rst(2)-(K2-K1)*gp*ds/r1**2
          if(dabs(G2-G1)>1.d-10)rst(5)=rst(5)+(G2-G1)*gp*ds
	  enddo
	enddo
9002	return
      end
