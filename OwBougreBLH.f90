      subroutine OwBougreBLH(BLH,sea,sfh,nlat,nlon,hd,dr,GRS,ter)
      !按严密积分公式计算地球外空间BLH点的扰动重力海水完全布格影响
      !dr-积分半径km
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj,astat(5)
	real*8::dr,sea(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst(7),ter(7)
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3),BLH2(3),XYZ2(3)
	real*8 CGrdPntD2,zch,L0,L1,L2,K1,K2,P1,P2,fx
	real*8::GRS(6),ge,re,tp,wp
!-----------------------------------------------------------------
      ge=6.67428d-11;tp=2.67d3;wp=1.03d3!地形密度海水密度
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;gp=ge*(tp-wp)
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      ter=0.d0;mdr=r0*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/r0/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
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
             call OwBougresgn(BLH,sea,sfh,nlat,nlon,hd,i,j,4,GRS,rst)
             ter=ter+rst; goto 9101 
          endif
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          P1=dlog(r1-rr*tt+L1);K1=-r1/rr/L1
          BLH2=BLH1;BLH2(3)=sfh(i,j)+sea(i,j)!-zch
          call BLH_XYZ(GRS,BLH2,XYZ2)
          L2=dsqrt((XYZ2(1)-XYZ(1))**2+(XYZ2(2)-XYZ(2))**2+(XYZ2(3)-XYZ(3))**2)
          r2=dsqrt(XYZ2(1)**2+XYZ2(2)**2+XYZ2(3)**2)
          P2=dlog(r2-rr*tt+L2);K2=-r2/rr/L2
          ter(1)=ter(1)+(P2-P1)*gp*ds/gr
          ter(2)=ter(2)-(K2-K1)*gp*ds
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
!--------------------------------------------------------------------------------
      subroutine OwBougresgn(BLH,sea,sfh,nlat,nlon,hd,i0,j0,m,GRS,rst)
      !细化核函数，计算BLH点的扰动重力海水完全布格影响的奇异积分
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
      implicit none
	integer::m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,sea(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst(7),rv,lon,lat
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3),BLH2(3),XYZ2(3)
	real*8 CGrdPntD2,zch,zch1,L0,L1,L2,P1,P2,K1,K2,fx,dpth
	real*8::GRS(6),ge,re,tp,wp
!-----------------------------------------------------------------
      ge=6.67428d-11;tp=2.67d3;wp=1.03d3!地形密度海水密度
      rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0;gp=ge*(tp-wp)
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      rst=0.d0;mdr=r0*rv*RAD*dcos(rln(2)*RAD)/dble(m)/2.d0  !奇异点判断
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),sfh,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          dpth=CGrdPntD2(BLH1(2),BLH1(1),sea,nlat,nlon,hd)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=rv**2*RAD**2*dcos(rln1(2)*RAD)*r1**2
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          if(L1<mdr)L1=mdr
          P1=dlog(r1-rr*tt+L1);K1=-r1/rr/L1
          BLH2=BLH1;BLH2(3)=BLH1(3)+dpth
          call BLH_XYZ(GRS,BLH2,XYZ2)
          L2=dsqrt((XYZ2(1)-XYZ(1))**2+(XYZ2(2)-XYZ(2))**2+(XYZ2(3)-XYZ(3))**2)
          r2=dsqrt(XYZ2(1)**2+XYZ2(2)**2+XYZ2(3)**2)
          if(L2<mdr)L2=mdr
          P2=dlog(r2-rr*tt+L2);K2=-r2/rr/L2
          if(dabs(P2-P1)>1.d-10)rst(1)=rst(1)+(P2-P1)*gp*ds/gr
          if(dabs(K2-K1)>1.d-10)rst(2)=rst(2)-(K2-K1)*gp*ds
	  enddo
	enddo
9002	return
      end
