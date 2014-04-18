c      SUBROUTINE READPT (nflpar,partcl,InOut,nptot0,Nptot0N,
c     *                   x,y,z,Timeb,px,py,pz,Tbp,dNbp,Chb,Ebp,
c     *                   Chsb,Pbp,Rpfl,Qbp,ifluid,
c     *                 ApN,AtN,EkinN,trl,tfreg,pibint,interEoS,
c     *                 nroN,ieos,edefr0,ehdfr0,RpN,RtN,hxN,hyN,hzN,dtN,
c     *                 uN,guN,V0N,g0N,RoN,RoiN,NML,jFreez,timefrz,
c     *                 Npr,Ntg1,Nptotin,Qp,Qt,dvofr,alf0,Kdump,
c     *                 gamh,beta,strm,iextra,gamq,betaq,iextraq,
c     *                 Qbstop,spectat,Lpercell,dxdz,dxdt,Lnucl,
c     *                 hxmin,hymin,hzmin,dtmin,hxnew,hynew,hznew,dtnew,
c     *                 Npmin,Npmax,Nptot)
      subroutine readpt
      parameter   (Npi0=120000,Nptot0=163000)
      dimension   x(nptot0),y(nptot0),z(nptot0),Timeb(nptot0),
     *            px(nptot0),py(nptot0),pz(nptot0),
     *            Tbp(nptot0),dNbp(nptot0),Chb(nptot0),Ebp(nptot0),
     *            Chsb(nptot0),Pbp(nptot0),Rpfl(nptot0),
     *            Qbp(nptot0),ifluid(nptot0)
      logical     tfreg
      integer     pibint
      character*50 partcl
      partcl = 'Au30mix_b6_Npi2_i3_Bps.dat'
*-------------------------------------------------------------*
       IF (InOut.eq.0) THEN
        OPEN (unit=nflpar,file=partcl,
     *        status='unknown',form='unformatted')
        READ (Nflpar) ApN,AtN,EkinN,trl,tfreg,pibint,interEoS,
     *                nroN,Lpercell,dxdz,dxdt,Lnucl,ieos,edefr0,ehdfr0,
     *                RpN,RtN,hxN,hyN,hzN,dtN,hxmin,hymin,hzmin,dtmin,
     *                hxnew,hynew,hznew,dtnew,
     *                uN,guN,V0N,g0N,RoN,RoiN,NML,jFreez,timefrz,
     *                Npr,Ntg1,Nptotin,Nptot0N,Qp,Qt,dvofr,alf0,Kdump,
     *                gamh,beta,strm,iextra,gamq,betaq,iextraq,
     *                Qbstop,spectat,Npmin,Npmax,Nptot
	  DO i=1,Nptot
         read (Nflpar) ifluid(i),x(i),y(i),z(i),Timeb(i),
     *        px(i),py(i),pz(i),Tbp(i),Chb(i),dNbp(i),Ebp(i),
     *        Chsb(i),Pbp(i),Qbp(i),Rpfl(i)
	  END DO
       ELSE IF (InOut.eq.1) THEN
        OPEN (unit=nflpar,file=partcl,
     *        status='unknown',form='formatted')
        READ (Nflpar,*) ApN,AtN,EkinN,trl,tfreg,pibint,interEoS,
     *                nroN,Lpercell,dxdz,dxdt,Lnucl,ieos,edefr0,ehdfr0,
     *                RpN,RtN,hxN,hyN,hzN,dtN,hxmin,hymin,hzmin,dtmin,
     *                hxnew,hynew,hznew,dtnew,
     *                uN,guN,V0N,g0N,RoN,RoiN,NML,jFreez,timefrz,
     *                Npr,Ntg1,Nptotin,Nptot0N,Qp,Qt,dvofr,alf0,Kdump,
     *                gamh,beta,strm,iextra,gamq,betaq,iextraq,
     *                Qbstop,spectat,Npmin,Npmax,Nptot
	  DO i=1,Nptot
         read (Nflpar,*) ifluid(i),x(i),y(i),z(i),Timeb(i),
     *        px(i),py(i),pz(i),Tbp(i),Chb(i),dNbp(i),Ebp(i),
     *        Chsb(i),Pbp(i),Qbp(i),Rpfl(i)
	  END DO
       END IF
       CLOSE (Nflpar)
       RETURN
      END
