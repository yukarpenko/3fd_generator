      SUBROUTINE READFRBL (Nfpisr,ffpisr,InOut,Npi0,
     *           xpi,ypi,zpi,Timepi,pxpi,pypi,pzpi,eppi,Etpip,Tpip,Vpip,
     *           ApP,AtP,EkinP,trlP,tfregP,pibintP,interEoSp,
     *           nroP,ieosP,edefr0P,ehdfr0P,RpP,RtP,hxP,hyP,hzP,dtP,
     *           uP,guP,V0P,g0P,Npipc,taupif,dvofr,alf0,Kdump,RoP,RoiP,
     *           NML,jFreezP,timefrzP,Lpercell,dxdz,dxdt,Lnucl,
     *           gamh,betaP,strmP,iextraP,QbstopP,spectatP,
     *           gamqP,betaqP,iextraqP,hxmin,hymin,hzmin,dtmin,
     *           hxnew,hynew,hznew,dtnew,Npmin,Npmax,Jpi)
     *           
      dimension   xpi(Npi0), ypi(Npi0), zpi(Npi0),Timepi(Npi0),
     *            pxpi(Npi0),pypi(Npi0),pzpi(Npi0),eppi(Npi0),
     *            Etpip(Npi0),Tpip(Npi0),Vpip(Npi0)
      logical     tfregP
      integer     pibintP
      character*50 ffpisr
*-------------------------------------------------------------*
        IF (InOut.eq.0) THEN
         OPEN (Nfpisr,file=ffpisr,status='unknown',
     *         form='unformatted')
         READ (Nfpisr) ApP,AtP,EkinP,trlP,tfregP,pibintP,interEoSp,
     *        nroP,Lpercell,dxdz,dxdt,Lnucl,ieosP,edefr0P,ehdfr0P,
     *        RpP,RtP,hxP,hyP,hzP,dtP,uP,guP,V0P,g0P,Npipc,taupif,dvofr,
     *        alf0,Kdump,RoP,RoiP,NML,jFreezP,timefrzP,
     *        hxmin,hymin,hzmin,dtmin,hxnew,hynew,hznew,dtnew,
     *        gamh,betaP,strmP,iextraP,QbstopP,spectatP,
     *        gamqP,betaqP,iextraqP,Npmin,Npmax,Jpi
	   DO j=1,Jpi
          read (Nfpisr) xpi(j),ypi(j),zpi(j),Timepi(j),
     *                 pxpi(j),pypi(j),pzpi(j),eppi(j),Vpip(j),
     *                 Tpip(j),Etpip(j)
	   END DO
        ELSE IF (InOut.eq.1) THEN
         OPEN (Nfpisr,file=ffpisr,status='unknown',
     *         form='formatted')
         READ (Nfpisr,*) ApP,AtP,EkinP,trlP,tfregP,pibintP,interEoSp,
     *        nroP,Lpercell,dxdz,dxdt,Lnucl,ieosP,edefr0P,ehdfr0P,
     *        RpP,RtP,hxP,hyP,hzP,dtP,uP,guP,V0P,g0P,Npipc,taupif,dvofr,
     *        alf0,Kdump,RoP,RoiP,NML,jFreezP,timefrzP,
     *        hxmin,hymin,hzmin,dtmin,hxnew,hynew,hznew,dtnew,
     *        gamh,betaP,strmP,iextraP,QbstopP,spectatP,
     *        gamqP,betaqP,iextraqP,Npmin,Npmax,Jpi
	   DO j=1,Jpi
          read (Nfpisr,*) xpi(j),ypi(j),zpi(j),Timepi(j),
     *                 pxpi(j),pypi(j),pzpi(j),eppi(j),Vpip(j),
     *                 Tpip(j),Etpip(j)
	   END DO
        END IF
        CLOSE (Nfpisr)
       RETURN
      END