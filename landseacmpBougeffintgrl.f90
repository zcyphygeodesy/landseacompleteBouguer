!  landseacmpBougeffintgrl.f90 
!
!  FUNCTIONS:
!  landseacmpBougeffintgrl - Entry point of console application.
!
!****************************************************************************

      program landseacmpBougeffintgrl
      !calculate the land-sea unified complete Bouguer effect on gravity using numerical integral
      implicit none
 	character*800::line,line0
      integer knd,row,nk,sn,len,astat(8)
	integer i,j,nlon,nlat,kk
	real*8::hd(6),hd1(6),rec(800),CGrdPntD2
	real*8::GRS(6),BLH(3),XYZ(3),pi,RAD,tmp(280000)
	real*8::dr1,dr2,rr,zch,r0,ge,tp,seag,lndg,sph,rst,ter(7)
	real*8,allocatable::lnd(:,:),sea(:,:),sfh(:,:)
	integer::status=0
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5) = 1.d0/298.2564619427d0
   	pi=datan(1.d0)*4.d0;RAD=pi/180.d0;ge=6.67428d-11; tp=2.67d3
	dr1=9.d4!dr1 - 局部地形影响积分半径 the integral radius (m) for local terrain effect
	dr2=30.d4!dr2 - 海水完全布格影响积分半径 the integral radius (m) for sear-water complete Bouguer effect
      open(unit=8,file="dtm5m.dat",status="old",iostat=status)
      if(status/=0)goto 902  !the land-sea terrain model
      read(8,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/real(nlon)
	hd(6)=(hd(4)-hd(3))/real(nlat)
 	allocate(lnd(nlat,nlon), stat=astat(1))
 	allocate(sea(nlat,nlon), stat=astat(2))
 	allocate(sfh(nlat,nlon), stat=astat(3))!The land-sea surface ellipsoidal height grid
	if (sum(astat(1:3)) /= 0) then
          close(8);goto 902
      endif
      lnd=0.d0;sea=0.d0
 	do i=1,nlat
	   read(8,*,end=903)(tmp(j),j=1,nlon)
         do j=1,nlon
           if(tmp(j)>=0.d0)lnd(i,j)=tmp(j)
           if(tmp(j)<0.d0)sea(i,j)=tmp(j)
         enddo
      enddo
903   close(8)
      open(unit=10,file="dbmhgt5m.dat",status="old",iostat=status)
      if(status/=0)goto 904  
      read(10,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd1(1:6)=rec(1:6)
      if(sum(hd1-hd)>1.d-5)then  !格网规格不同 The grid specifications are different
         close(10);goto 904
      endif
 	do i=1,nlat
	   read(10,*,end=905)(sfh(i,j),j=1,nlon)
      enddo
905   close(10)
      open(unit=8,file="surfhgt.txt",status="old",iostat=status)
      if(status/=0)goto 904
      open(unit=10,file="reslt.txt",status="replace")
      read(8,'(a)') line  !读取头文件
      write(10,101)trim(line);kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line0
         call PickRecord(line0,len,rec,sn)
         if(sn<4)goto 906
         BLH(2)=rec(2);BLH(1)=rec(3)
         BLH(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!land-sea surface ellipsoidal
         call BLH_XYZ(GRS,BLH,XYZ)
         r0=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
         zch=CGrdPntD2(BLH(2),BLH(1),lnd,nlat,nlon,hd)
         BLH(3)=rec(4)!ellipsoidal height of calculation point
         call BLH_XYZ(GRS,BLH,XYZ)
         rr=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
         !计算局部地形影响
         !calculate the local terrain effect
         call LTerAllBLH(BLH,lnd,sfh,nlat,nlon,hd,dr1,GRS,ter)
         lndg=ter(2)*1.d5;kk=kk+1
         !计算球壳布格影响
         !calculate the the spherical shell Bouguer effect
	   sph=4.d0*pi*ge*tp*(r0**2*zch+r0*zch**2+zch**3/3.d0)/rr**2*1.d5
         !计算海水完全布格
         !calculate the sear-water complete Bouguer effect
         call OwBougreBLH(BLH,sea,sfh,nlat,nlon,hd,dr2,GRS,ter)
         seag=ter(2)*1.d5;rst=lndg+sph+seag
         write(10,101)trim(line0),lndg,sph,seag,rst
         if(kk/250*250==kk)write(*, '(a,i9)'), '    Calculated point number: ',kk
       enddo
906   close(8)
      close(10)
904   deallocate(sea,lnd,sfh)
902   continue
101   format(a,40F12.4)
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
