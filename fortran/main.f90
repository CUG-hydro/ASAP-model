!now run the model
! if(freedrain.eq.0.and.hour.eq.nint(deltat/3600.)) then
call LATERAL(n2, n3, is, ie, js, je, nzg, soiltxt, wtd, qlat, fdepth, topo, landmask, deltatwtd, area, lats, dxy, slz &
             , o18, smoi, qlato18, qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum)
qslat = 0.
qlatsum = qlatsum + qlat*1.e3
! endif

if (riverswitch .eq. 1) then
   call GW2RIVER(n2, n3, is, ie, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth &
                 , riverwidth, riverlength, area, fdepth, qrf)
end if

t3 = mpi_wtime()

call ROOTDEPTH(freedrain,is,ie,js,je,nzg,slz,dz,deltat,landmask,veg,hveg,soiltxt,wind,temp,qair,press,netrad,rshort &
               , lai, precip, qsrun, smoi, smoieq, smoiwtd, wtd, waterdeficit, watext, watextdeep, rech, deeprech &
               , et_s, et_i, et_c, intercepstore, ppacum, pppendepth, pppendepthold &
               , qlat*deltat/deltatwtd, qslat, qsprings, inactivedays, maxinactivedays, fieldcp, fdepth, steps, floodheight &
               , qrf, delsfcwat, icefactor &
        ,wtdflux(is,js,daypast),et_s_daily(is,js,daypast),et_c_daily(is,js,daypast),transptop(is,js,daypast),infilk(is,js,daypast) &
               , infilflux, infilfluxday, infilcounter, hour, o18, o18ratiopp, tempsfc, qlato18*deltat/deltatwtd, transpo18 &
               , upflux)

wtdmax = max(wtdmax, wtd)
wtdmin = min(wtdmin, wtd)

qrfsum = qrfsum + qrf*1.e3
qsrunsum = qsrunsum + qsrun*1.e3
delsfcwatsum = delsfcwatsum - delsfcwat*1.e3

if (mod(daysfromstart - 1, 5) .eq. 0 .and. hour .eq. 0) then
   wtd_daily(:, :, daypast) = wtd(:, :)
   smoi_daily(is:ie, js:je, daypast) = smoi(nzg, is:ie, js:je)
   daysforoutput(daypast) = daysfromstart - 1
   daypast = daypast + 1 ! advance counter
end if

!this is not needed now, all layers are treaded the same
!       if(freedrain.eq.0.and.hour.eq.3) call WTABLE(n2,n3,is,ie,js,je,nzg,slz,dz,area,soiltxt,wtd,deeprech,rech,qslat,fdepth &
!                                     ,topo,landmask,deltatwtd,smoi,smoieq,smoiwtd,qsprings)

if (riverswitch .eq. 1) then

call FLOODING(n2, n3, is, ie, js, je, deltat, fd, bfd, topoflood, area, riverwidth, riverlength, riverdepth, floodheight, delsfcwat)

   do iter_river = 1, niter_river

      call RIVERS_KW_FLOOD(n2, n3, is, ie, js, je, deltat, dtlr, fd, bfd, riverflow, qsrun, qrf, delsfcwat &
                           , slope, riverdepth, riverwidth, riverlength, maxdepth, area, riverarea, floodarea, riverchannel &
                           , riverflowmean, floodheight, topoflood)

   end do

   qsrun = 0.
   delsfcwat = 0.

end if
