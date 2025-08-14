!******************************************************************************
   subroutine SOILFLUXES(i, j, nzg, freedrain, dtll, slz, dz, soiltxt, smoiwtd, transp, transpdeep &
                         , smoi, wtd, rech, deeprech, precip, pet_s, et_s, runoff, flux, fdepth, qlat &
                         , qlatflux, qrf, qrfcorrect, flood, icefactor &
                         , smoieq, o18, precipo18, tempsfc, qlato18, transpo18)!pppendepth)
      implicit none

      integer :: nzg, freedrain, nsoil, nsoil1, k, iwtd, kwtd, i, j
      real, dimension(nzg + 1) :: slz
      real, dimension(nzg) :: dz
      real, dimension(nzg) ::vctr2, vctr4, vctr5, vctr6
      real, dimension(nzg) :: transp, smoi, kfmid, diffmid &
                              , aa, bb, cc, rr, smoieq, o18, o18dz
      real, dimension(nzg) :: smoiold, o18ratio
      real*8, dimension(nzg + 1) :: vt3di, o18flux, gravflux, capflux
      real, dimension(nzg + 1) :: flux
      real, dimension(0:nzg + 1) :: qlatflux
      integer, dimension(2) :: soiltxt
      integer*1, dimension(nzg) :: icefactor
      real :: precip, runoff, pet_s, et_s, transpdeep, pppendepth, precipo18, tempsfc, qgwo18
      real :: wgpmid, kfup, kfdw, hydcon, newwgp, smoiwtd, rech, deeprech, wtd, deeptemp &
              , fracliqwtd, wmid, wtdold, dzup, vt3dbdw, vt3dcdw, dtll, smoibot, icefac, ddw, dup &
              , smoisat, psisat, smoicp, fdepth, qlat, qlato18, qrf, qrfcorrect, qgw, flood, transpo18
      real*8 :: alpha, o18evap, dsmoi, transptot, qlatlayer, qrflayer, o18frac, o18out, o18tot, fluxdiff

      do k = 1, nzg
         vctr2(k) = 1./dz(k)
         vctr4(k) = 0.5*(slz(k) + slz(k + 1))
      end do
      do k = 2, nzg
         vctr5(k) = vctr4(k) - vctr4(k - 1)
         vctr6(k) = 1./vctr5(k)
      end do

      kfmid = 0.
      diffmid = 0.

      vt3di = 0.
      gravflux = 0.
      capflux = 0.

      rech = 0.
      runoff = 0.

!o18 ratio
      o18ratio = o18/smoi
      smoiold = smoi

      qgw = qlat - qrf

!top boundary condition, infiltration + potential et from soil
      nsoil = soiltxt(1)
      smoicp = soilcp(nsoil)
      if (smoi(nzg) .le. smoicp) pet_s = 0.
      vt3di(nzg + 1) = (-precip + pet_s)*1.e-3 - flood

      if (-vt3di(nzg + 1) .gt. slcons(nsoil)*dtll) then
         runoff = -vt3di(nzg + 1) - slcons(nsoil)*dtll
         vt3di(nzg + 1) = -slcons(nsoil)*dtll
      end if

!         smoisat = slmsts(nsoil)
!         dsmoi = max((smoisat-smoi(nzg))*dz(nzg)+transp(nzg),0.)
!         if(-vt3di(nzg+1).gt.dsmoi)then
!             runoff = -vt3di(nzg+1)-dsmoi
!             vt3di(nzg+1)=-dsmoi
!         endif

      if (freedrain .eq. 0) then
         do k = 1, nzg
            if (wtd .lt. slz(k)) exit
         end do
         iwtd = k
      else
         iwtd = 0
      end if

!k=max(iwtd-1,1)
!qlatflux(k)=qlatflux(k)+qgw

!         do k = 2,nzg

      do k = max(iwtd - 1, 2), nzg

!gmmdiffusivity and conductivity at the interlayers

         wgpmid = smoi(k) + (smoi(k) - smoi(k - 1))*(slz(k) - vctr4(k))*vctr6(k)

         if (slz(k) .lt. -0.30) then
            nsoil = soiltxt(1)
         else
            nsoil = soiltxt(2)
         end if

         hydcon = slcons(nsoil)*max(min(exp((slz(k) + 1.5)/fdepth), 1.), 0.1)
         smoisat = slmsts(nsoil)*max(min(exp((slz(k) + 1.5)/fdepth), 1.), 0.1)
         psisat = slpots(nsoil)*min(max(exp(-(slz(k) + 1.5)/fdepth), 1.), 10.)

         wgpmid = min(wgpmid, smoisat)
!            icefac=fracliq(k)** (2. * slbs(nsoil) + 3.)
!            icefac=1.
         if (icefactor(k) .eq. 0) then
            icefac = 1.
         else
            icefac = 0.
         end if

         kfmid(k) = icefac*hydcon &
                    *(wgpmid/smoisat)**(2.*slbs(nsoil) + 3.)
         diffmid(k) = -icefac*(hydcon*psisat*slbs(nsoil)/smoisat) &
                      *(wgpmid/smoisat)**(slbs(nsoil) + 2.)

!write(6,*)k,diffmid(k),kfdw,kfup,ddw,dup,smoi(k),smoi(k-1)

      end do

!calculate tridiagonal matrix elements

!       do k=2,nzg-1
!         do k=max(iwtd-2,2),nzg
      do k = max(iwtd, 3), nzg

         aa(k) = diffmid(k)*vctr6(k)
         cc(k) = diffmid(k + 1)*vctr6(k + 1)
         bb(k) = -(aa(k) + cc(k) + dz(k)/dtll)
         rr(k) = -smoi(k)*dz(k)/dtll - kfmid(k + 1) + kfmid(k) + transp(k)/dtll
!            if(k.eq.iwtd-1)rr(k) = rr(k) - qgw/dtll

      end do

!boundary conditions

!top boundary

      if (iwtd - 1 .eq. nzg) then
         aa(nzg) = 0.
         cc(nzg) = 0.
         bb(nzg) = -dz(nzg)/dtll
            rr(nzg) = vt3di(nzg+1)/dtll -smoi(nzg)*dz(nzg)/dtll + transp(k)/dtll + min( kfmid(nzg)+ diffmid(nzg)*vctr6(nzg)*(smoi(nzg)-smoi(nzg-1)) , 0. )
      else
         aa(nzg) = diffmid(nzg)*vctr6(nzg)
         cc(nzg) = 0.
         bb(nzg) = -aa(nzg) - dz(nzg)/dtll
         rr(nzg) = vt3di(nzg + 1)/dtll - smoi(nzg)*dz(nzg)/dtll + kfmid(nzg) + transp(nzg)/dtll
!            if(iwtd-1.eq.nzg)rr(nzg) = rr(nzg) - qgw/dtll
      end if

!now bottom boundary condition

      IF (freedrain .ne. 1) then

         if (iwtd .le. 2) then
            aa(1) = 0.
            cc(1) = diffmid(2)*vctr6(2)
            bb(1) = -(cc(1) + dz(1)/dtll)
            rr(1) = -smoi(1)*dz(1)/dtll - kfmid(2) + transp(1)/dtll
!            if(iwtd.le.2)rr(1) = rr(1) - qgw/dtll

            k = 2
            aa(k) = diffmid(k)*vctr6(k)
            cc(k) = diffmid(k + 1)*vctr6(k + 1)
            bb(k) = -(aa(k) + cc(k) + dz(k)/dtll)
            rr(k) = -smoi(k)*dz(k)/dtll - kfmid(k + 1) + kfmid(k) + transp(k)/dtll
         else
            do k = 1, iwtd - 3
               aa(k) = 0.
               cc(k) = 0.
               bb(k) = 1.
               rr(k) = smoi(k)
            end do

            k = iwtd - 1 ! layer where the water table is
            aa(k) = 0.
            cc(k) = diffmid(k + 1)*vctr6(k + 1)
            bb(k) = -(cc(k) + dz(k)/dtll)
       rr(k) = -smoi(k)*dz(k)/dtll - kfmid(k + 1) + transp(k)/dtll + min(kfmid(k) + diffmid(k)*vctr6(k)*(smoi(k) - smoi(k - 1)), 0.)

            k = iwtd - 2
            aa(k) = 0.
            cc(k) = 0.
            bb(k) = -dz(k)/dtll
            rr(k) = -smoi(k)*dz(k)/dtll + max(-kfmid(k + 1) - diffmid(k + 1)*vctr6(k + 1)*(smoi(k + 1) - smoi(k)), 0.)
         end if

      ELSE

!gmmgravitational drainage at the bottom
         nsoil = soiltxt(1)

         hydcon = slcons(nsoil)*max(min(exp((slz(1) + 1.5)/fdepth), 1.), 0.1)
         smoisat = slmsts(nsoil)*max(min(exp((slz(1) + 1.5)/fdepth), 1.), 0.1)

         kfmid(1) = hydcon &
                    *(smoi(1)/smoisat)**(2.*slbs(nsoil) + 3.)

         aa(1) = 0.
         cc(1) = diffmid(2)*vctr6(2)
         bb(1) = -(cc(1) + dz(1)/dtll)
         rr(1) = -smoi(1)*dz(1)/dtll - kfmid(2) + kfmid(1) + transp(1)/dtll

      END IF

!solve tridiagonal system and update smoi

      call tridag(aa, bb, cc, rr, smoi, nzg)

!calculate the fluxes

      do k = max(iwtd, 3), nzg
         gravflux(k) = -kfmid(k)*dtll
         capflux(k) = -aa(k)*(smoi(k) - smoi(k - 1))*dtll
         vt3di(k) = capflux(k) + gravflux(k)
!                if(k.le.iwtd-1)vt3di(k)=max(vt3di(k),0.)
      end do
      if (iwtd .le. 2) then
         capflux(1) = 0.
         gravflux(1) = 0.
         vt3di(1) = 0.
         gravflux(2) = -kfmid(2)*dtll
         capflux(2) = -aa(2)*(smoi(2) - smoi(1))*dtll
         vt3di(2) = capflux(2) + gravflux(2)
      else
         do k = 1, iwtd - 2
            capflux(k) = 0.
            gravflux(k) = 0.
            vt3di(k) = 0.
         end do

         k = iwtd - 1
         gravflux(k) = -kfmid(k)*dtll
         capflux(k) = -diffmid(k)*vctr6(k)*(smoiold(k) - smoiold(k - 1))*dtll
         if (capflux(k) .gt. -gravflux(k)) then
            vt3di(k) = capflux(k) + gravflux(k)
         elseif (capflux(k) .gt. 0.) then
            gravflux(k) = -capflux(k)
            vt3di(k) = 0.
         else
            capflux(k) = 0.
            gravflux(k) = 0.
            vt3di(k) = 0.
         end if

      end if

      if (freedrain .eq. 0) then
         vt3di(1) = 0.
      else
         vt3di(1) = -kfmid(1)*dtll
      end if

      smoi = smoiold
!recalculate soil moisture
      do k = 1, nzg
         smoiold(k) = smoiold(k) + (vt3di(k) - vt3di(k + 1) - transp(k))*vctr2(k)
!          if(k.eq.iwtd-1)then
!                  smoiold(k) = smoiold(k) + qgw * vctr2(k)
!          endif
      end do

!if(i.eq.25.and.j.eq.30)write(6,*)'fluxes antes',(vt3di(k),k=1,nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'capillarity',(-aa(k)*(smoi(k)-smoi(k-1)),k=2,nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'drainage',(-kfmid(k)*dtll,k=2,nzg)

! now check that soil moisture values are within bounds (slmsts and soilcp)
! if not, correct fluxes

      do k = 1, nzg

         if (slz(k) .lt. -0.30) then
            nsoil = soiltxt(1)
         else
            nsoil = soiltxt(2)
         end if

         smoisat = slmsts(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)

!if(i.eq.69.and.j.eq.34)write(6,*)'soilflux',k,smoi(k),smoiold(k),smoisat,vt3di(k+1),transp(k),iwtd-1,qgw,qlat,-qrf

         if (smoiold(k) .gt. smoisat) then
            dsmoi = max((smoiold(k) - smoisat)*dz(k), 0.)
            if (k .lt. nzg) then
               smoiold(k + 1) = smoiold(k + 1) + dsmoi*vctr2(k + 1)
               vt3di(k + 1) = vt3di(k + 1) + dsmoi
            else
               vt3di(k + 1) = vt3di(k + 1) + dsmoi
               runoff = runoff + dsmoi
            end if
            smoiold(k) = smoisat
            if (capflux(k + 1) .lt. 0) then
               gravflux(k + 1) = gravflux(k + 1) + capflux(k + 1)
               capflux(k + 1) = 0.
            end if
            gravflux(k + 1) = gravflux(k + 1) + dsmoi
            if (gravflux(k + 1) .gt. 0.) then
               capflux(k + 1) = capflux(k + 1) + gravflux(k + 1)
               gravflux(k + 1) = 0.
            end if
         end if
!if(i.eq.69.and.j.eq.34)write(6,*)'soilflux 2',k,smoiold(k),vt3di(k+1),dsmoi,iwtd,wtd
      end do

!        do k=nzg,1

!            if(slz(k).lt.-0.30)then
!                  nsoil=soiltxt(1)
!            else
!                  nsoil=soiltxt(2)
!            endif
!            smoicp = soilcp(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth),1.),0.1)

!             if(smoiold(k).lt.smoicp)then
!                 dsmoi=max((smoicp-smoiold(k))*dz(k),0.)
!                 if(vt3di(k).lt.0.)then
!                      vt3di(k)=vt3di(k)+dsmoi
!                      dsmoi=max(vt3di(k),0.)
!                      vt3di(k)=min(vt3di(k),0.)
!                 endif
!                 smoiold(k+1)=smoiold(k+1)-dsmoi*vctr2(k+1)
!                 vt3di(k+1)=vt3di(k+1)-dsmoi
!                 smoiold(k)=smoicp
!             endif
!        enddo

      k = nzg

      if (slz(k) .lt. -0.30) then
         nsoil = soiltxt(1)
      else
         nsoil = soiltxt(2)
      end if

      smoicp = soilcp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)

      if (smoiold(k) .lt. smoicp) then
!first reduce soil evaporation from PET
         dsmoi = max((smoicp - smoiold(k))*dz(k), 0.)
         if (vt3di(k + 1) .gt. dsmoi) then
            et_s = max(0., pet_s - dsmoi*1.e3)
            smoiold(k) = smoicp
            vt3di(k + 1) = vt3di(k + 1) - dsmoi
         else
            et_s = max(0., pet_s - max(vt3di(k + 1), 0.)*1.e3)
            vt3di(k + 1) = min(vt3di(k + 1), 0.)
            smoiold(k) = smoiold(k) + max(vt3di(k + 1), 0.)/dz(k)
!take water from below
            dsmoi = max((smoicp - smoiold(k))*dz(k), 0.)
            smoiold(k - 1) = smoiold(k - 1) - dsmoi*vctr2(k - 1)
            vt3di(k) = vt3di(k) + dsmoi
            smoiold(k) = smoicp
         end if
      else
         et_s = pet_s
      end if

!then go down all the way to the bottom
      do k = nzg - 1, 1, -1
         if (slz(k) .lt. -0.30) then
            nsoil = soiltxt(1)
         else
            nsoil = soiltxt(2)
         end if

         smoicp = soilcp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
         if (smoiold(k) .lt. smoicp) then
            !take water from below
            dsmoi = max((smoicp - smoiold(k))*dz(k), 0.)
            if (k .gt. 1) smoiold(k - 1) = smoiold(k - 1) - dsmoi*vctr2(k - 1)
            vt3di(k) = vt3di(k) + dsmoi
            smoiold(k) = smoicp
         end if
      end do

      if (vt3di(1) .gt. 0.) then
         qrfcorrect = -min(vt3di(1), max(qrf, 0.))
!       write(6,*)'too much qrf',i,j,qrf,qrfcorrect,vt3di(1),qlat,wtd
      else
         qrfcorrect = 0.
      end if

!save rain penetration depth

      flux = flux + vt3di

!     do k=1,nzg
!       if(vt3di(k).lt.-1.e-6)then
!            if(pppendepth.gt.slz(k))pppendepth=slz(k)
!            exit
!       endif
!     enddo
      IF (freedrain .eq. 1) then

!accumulate gravitational drainage
         rech = vt3di(1)
!smoiwtd is now the bucket of water at the bottom, to save the water and put it later into the rivers
         smoiwtd = smoiwtd - vt3di(1)

      END IF

!now o18 tracer

      vt3di(1) = 0.
      capflux(1) = 0.
      gravflux(1) = 0.

      transpo18 = 0.

      do k = 2, nzg
         if (capflux(k) .lt. 0) then
            gravflux(k) = gravflux(k) + capflux(k)
            capflux(k) = 0.
         end if
         fluxdiff = vt3di(k) - (gravflux(k) + capflux(k))
         if (fluxdiff .gt. 0.) then
            capflux(k) = capflux(k) + fluxdiff
         else
            gravflux(k) = gravflux(k) + fluxdiff
         end if

      end do

!calculate tridiagonal matrix elements

      do k = 1, nzg - 1

         bb(k) = 1.+0.5*vctr2(k)*transp(k)/smoiold(k)
         rr(k) = o18(k) - 0.5*vctr2(k)*transp(k)*o18ratio(k)

         bb(k) = bb(k) - 0.5*vctr2(k)*gravflux(k)/smoiold(k)
         rr(k) = rr(k) + 0.5*vctr2(k)*gravflux(k)*o18ratio(k)
         aa(k) = -0.5*vctr2(k)*capflux(k)/smoiold(k - 1)
         rr(k) = rr(k) + 0.5*vctr2(k)*capflux(k)*o18ratio(k - 1)

         cc(k) = 0.5*vctr2(k)*gravflux(k + 1)/smoiold(k + 1)
         rr(k) = rr(k) - 0.5*vctr2(k)*gravflux(k + 1)*o18ratio(k + 1)
         bb(k) = bb(k) + 0.5*vctr2(k)*capflux(k + 1)/smoiold(k)
         rr(k) = rr(k) - 0.5*vctr2(k)*capflux(k + 1)*o18ratio(k)

      end do

!top boundary

      bb(nzg) = 1.+0.5*vctr2(nzg)*transp(nzg)/smoiold(nzg)
      rr(nzg) = o18(nzg) - 0.5*vctr2(nzg)*transp(nzg)*o18ratio(nzg)

      bb(nzg) = bb(nzg) - 0.5*vctr2(nzg)*gravflux(nzg)/smoiold(nzg)
      rr(nzg) = rr(nzg) + 0.5*vctr2(nzg)*gravflux(nzg)*o18ratio(nzg)
      aa(nzg) = -0.5*vctr2(nzg)*capflux(nzg)/smoiold(nzg - 1)
      rr(nzg) = rr(nzg) + 0.5*vctr2(nzg)*capflux(nzg)*o18ratio(nzg - 1)

      alpha = 1./(exp(1137./tempsfc**2.-0.4156/tempsfc - 0.0020667))

      bb(nzg) = bb(nzg) + 0.5*vctr2(nzg)*(alpha*et_s*1.e-3)/smoiold(nzg)
      rr(nzg) = rr(nzg) - 0.5*vctr2(nzg)*(alpha*et_s*1.e-3)*o18ratio(nzg)

!         if(runoff.lt.precip *1.e-3 + flood )then
!              rr(nzg) = rr(nzg) + precipo18 * max(precip *1.e-3 + flood  - runoff ,0.) * vctr2(nzg)
      if (vt3di(nzg + 1) - et_s*1.e-3 .lt. 0.) then
         rr(nzg) = rr(nzg) + precipo18*max(et_s*1.e-3 - vt3di(nzg + 1), 0.)*vctr2(nzg)
      else
!              bb(nzg) = bb(nzg) + 0.5 * vctr2(nzg) * (runoff - precip *1.e-3 - flood ) / smoiold(nzg)
!              rr(nzg) = rr(nzg) - 0.5 * vctr2(nzg) * (runoff - precip *1.e-3 - flood ) *  o18ratio(nzg)
         bb(nzg) = bb(nzg) + 0.5*vctr2(nzg)*(vt3di(nzg + 1) - et_s*1.e-3)/smoiold(nzg)
         rr(nzg) = rr(nzg) - 0.5*vctr2(nzg)*(vt3di(nzg + 1) - et_s*1.e-3)*o18ratio(nzg)
      end if

!solve tridiagonal system and update smoi
!k=40
      o18dz = o18

      call tridag(aa, bb, cc, rr, o18, nzg)

!            o18 = o18dz / dz

      do k = 1, nzg
         transpo18 = transpo18 + 0.5*(o18(k)/smoiold(k) + o18ratio(k))*transp(k)
          if(o18(k).lt.0.)write(6,*)'O18 less than zero!!!',o18dz(k),o18(k),i,j,k,iwtd-1,vt3di(k),capflux(k),gravflux(k),vt3di(k+1),capflux(k+1),gravflux(k+1),o18ratio(k)*transp(k)*vctr2(k)
  if (o18(k) .gt. smoiold(k)) write (6, *) 'O18 greater than smoi!!!', o18dz(k), o18(k), vt3di(k), vt3di(k + 1), smoiold(k), i, j, k
      end do

!update wtd
      call updateshallowwtd(i, j, nzg, freedrain, slz, dz, soiltxt, smoieq, smoiwtd, smoiold, wtd, rech, fdepth)
      do k = 1, nzg
         if (wtd .lt. slz(k)) exit
      end do
      iwtd = k

      kwtd = max(iwtd - 1, 1)
!now lateral flow

!if(i.eq.20.and.j.eq.82)write(6,*)'now qlat',qlat

      qlatlayer = qgw
      if (qrf .gt. 0) then
         qgwo18 = qlato18 - o18ratio(kwtd)*qrf
      else
         qgwo18 = qlato18 - precipo18*qrf
      end if

!if(i.eq.20.and.j.eq.82)write(6,*)'now qlat',qlat,kwtd,o18(kwtd),qlato18/dz(k)

      if (qgw .gt. 0.) then
         do k = max(kwtd - 1, 1), nzg

            if (slz(k) .lt. -0.30) then
               nsoil = soiltxt(1)
            else
               nsoil = soiltxt(2)
            end if

            smoisat = slmsts(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
            dsmoi = (smoisat - smoiold(k))*dz(k)
            if (qlatlayer .le. dsmoi) then
               smoiold(k) = smoiold(k) + qlatlayer*vctr2(k)
               o18frac = qlatlayer/qgw
               o18frac = min(max(o18frac, 0.), 1.)
               o18(k) = o18(k) + o18frac*qgwo18*vctr2(k)
!if(i.eq.20.and.j.eq.82)write(6,*)'now qlat +',qlatlayer,k,o18(k)
               qlatflux(k) = qlatflux(k) + qlatlayer
               exit
            elseif (k .eq. kwtd - 1) then
               smoiold(k) = smoisat
               o18frac = dsmoi/qgw
               o18frac = min(max(o18frac, 0.), 1.)
               o18(k) = o18(k) + o18frac*qgwo18*vctr2(k)

               qlatlayer = qlatlayer - dsmoi
               qlatflux(k) = qlatflux(k) + dsmoi
!if(i.eq.100.and.j.eq.35)write(6,*)'now qlat +',qlatlayer,k,o18(k),dsmoi,( 1.e6 * (o18(k)/smoiold(k)) / 2005.2 - 1. ) * 1.e3

            else
               smoiold(k) = smoiold(k) + qlatlayer*vctr2(k)
               o18frac = qlatlayer/qgw
               o18frac = min(max(o18frac, 0.), 1.)
               o18(k) = o18(k) + o18frac*qgwo18*vctr2(k)
               o18ratio(k) = o18(k)/smoiold(k)
               !after mixing, take out what is over saturation
               smoiold(k) = smoisat
               o18(k) = smoisat*o18ratio(k)

               qlatlayer = qlatlayer - dsmoi
               qgwo18 = o18ratio(k)*qlatlayer
               qgw = qlatlayer
               qlatflux(k) = qlatflux(k) + dsmoi

!if(i.eq.100.and.j.eq.35)write(6,*)'now qlat +',qlatlayer,k,o18(k),dsmoi,( 1.e6 * (o18(k)/smoiold(k)) / 2005.2 - 1. ) * 1.e3
               if (k .eq. nzg) then
                  runoff = runoff + qlatlayer
                  exit
               end if
            end if
         end do
      elseif (qgw .lt. 0.) then
         do k = kwtd, 1, -1
            dsmoi = max((smoiold(k) - smoieq(k))*dz(k), 0.)
            if (k .eq. 1) then
               nsoil = soiltxt(2)
               smoicp = soilcp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
               dsmoi = max((smoiold(k) - smoicp)*dz(k), 0.)
            end if
            if (-qlatlayer .le. dsmoi) then
               smoiold(k) = smoiold(k) + qlatlayer*vctr2(k)
               o18frac = qlatlayer/qgw
               o18frac = min(max(o18frac, 0.), 1.)
               o18out = -o18frac*qgwo18
               if (o18(k)*dz(k) .lt. o18out .and. k .gt. 1) then !take part from the layer below
                  o18tot = o18(k)*dz(k) + o18(k - 1)*dz(k - 1)
                  o18frac = o18(k)*dz(k)/o18tot
                  o18(k) = o18(k) - o18frac*o18out*vctr2(k)
                  o18(k - 1) = o18(k - 1) - (1.-o18frac)*o18out*vctr2(k - 1)
               else
                  o18(k) = o18(k) - o18out*vctr2(k)
               end if
               qlatflux(k) = qlatflux(k) + qlatlayer
!if(i.eq.20.and.j.eq.82)write(6,*)'now qlat -',qlatlayer,k,o18(k),dsmoi,o18frac
               exit
            else
               qlatlayer = qlatlayer + dsmoi
               smoiold(k) = smoiold(k) - dsmoi*vctr2(k)

               o18frac = -dsmoi/qgw
               o18frac = min(max(o18frac, 0.), 1.)
               o18(k) = o18(k) + o18frac*qgwo18*vctr2(k)

               qlatflux(k) = qlatflux(k) + dsmoi
!if(i.eq.20.and.j.eq.82)write(6,*)'now qlat -',qlatlayer,k,o18(k),dsmoi,o18frac
               if (k .eq. 1) then
                  qrfcorrect = qrfcorrect - qlatlayer
                  exit
               end if
            end if
         end do
      end if

!          if(iwtd.eq.1)o18(1) = o18(1) + o18ratio (1) * qgw * vctr2(1)

      do k = 1, nzg
!      if(i.eq.100.and.j.eq.35.and.k.gt.35)write(6,*)'mirar o18 despues',i,j,k,o18(k),o18(k)/smoiold(k) &
!          , ( 1.e6 * (o18(k)/smoiold(k)) / 2005.2 - 1. ) * 1.e3
         if (o18(k) .lt. 0.) write (6, *) 'O18 less than zero!!!', o18(k), i, j, k, iwtd - 1, qlat, qlato18, qrf
         if (o18(k) .gt. smoiold(k)) write (6, *) 'O18 greater than smoi!!!', o18(k), smoiold(k), i, j, k, iwtd - 1, qlato18, qrf
      end do

      o18 = max(o18, 0.)

!          o18ratio = o18 / smoiold

      smoi = smoiold

   end subroutine soilfluxes
