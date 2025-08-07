program driver

   use module_parallel
   use module_forcings
   use module_rootdepth
   use module_io
   use module_wtable
   use module_initial

   implicit none

!integer, parameter :: n2=100,n3=100,nzg=40,restart=1,freedrain=0,riverswitch=1,nvar_out=15,writepar=0
   integer, parameter :: n2 = 7320, n3 = 8520, nzg = 40, restart = 1, freedrain = 0, riverswitch = 1, nvar_out = 15, writepar = 1
!integer, parameter :: n2=500,n3=500,nzg=32,restart=0,freedrain=0,nvar_out=13 !change this when writing in parallel
!integer, parameter :: n2=1100,n3=300,nzg=32,restart=1,freedrain=0,nvar_out=14,writepar=0 !change this when writing in parallel
!integer, parameter :: n2=400,n3=400,nzg=32,restart=1,freedrain=0,nvar_out=14,writepar=0
!integer, parameter :: n2=100,n3=100,nzg=40,restart=0,freedrain=0,nvar_out=12,writepar=0
!real, parameter :: deltat=1.*3600.,deltatwtd=24.*3600.,deltatriver=5.*60.,dxy=1./120.,steps=1.
   real, parameter :: deltat = 1.*3600., deltatwtd = 1.*3600., deltatriver = 5.*60., dxy = 1./120., steps = 1.
   integer, parameter :: maxinactivedays = 365*24*nint(steps)

   integer :: i, j, js, je, is, ie, k, nmax_x, nmax_y, icount, n, nday, niter_river, iter_river, monthhours

   real, dimension(nzg + 1) :: slz
   real, dimension(nzg) :: dz
   real :: inpair(2), outpair(2)
   integer, allocatable, dimension(:, :, :) :: soiltxt, inactivedays
   integer*1, allocatable, dimension(:, :, :) :: icefactor
   integer, allocatable, dimension(:, :) :: landmask, fd, bfd
   real, allocatable, dimension(:, :, :) :: smoi, watext, smoieq, infilflux, infilfluxday, o18, smoimax, smoimin, upflux
real, allocatable, dimension(:,:) :: topo,topoera5,area,lats,lons,wind,temp,qair,netrad,rshort,rlon,press,precip,snow,snow_past,snow_fut,rain,lai &
                                        , lai_past, lai_fut, o18ratiopp, o18_past, o18_fut, tempsfc &
                   , press_past, press_fut, temp_past, temp_fut, qair_past, qair_fut, wind_past, wind_fut, tempsfc_past, tempsfc_fut
real, allocatable, dimension(:,:) :: veg,hveg,wtd,smoiwtd,fdepth,rech,deeprech,qsrun,qsrunsum,qslat,qlat,qlato18,qlatsum,qlatinsum,qlatoutsum,qsprings &
          , et_s, et_i, et_c, intercepstore, ppacum, waterdeficit, watextdeep, pppendepth, transpo18, transpo18ratio, wtdmax, wtdmin
   real, allocatable, dimension(:, :) :: qlatino18sum, qlatouto18sum, qlatino18ratio, qlatouto18ratio
   integer*1, allocatable, dimension(:, :) :: pppendepthold
   real, allocatable, dimension(:, :) :: riverflow, qrf, qrfsum, delsfcwat, delsfcwatsum &
   , slope, riverdepth, riverwidth, riverlength, maxdepth, riverflowmean, floodheight, topoflood, riverarea, floodarea, riverchannel
   real, allocatable, dimension(:, :, :) :: wtdflux, et_s_daily, et_c_daily, transptop, wtd_daily, smoi_daily
   integer*1, allocatable, dimension(:, :, :) :: infilk
   integer*2, allocatable, dimension(:, :, :) :: infilcounter
   real, dimension(nzg, nstyp) :: fieldcp

   real, allocatable, dimension(:, :, :) :: varoutput, varoutput2, varoutput3, varhis, varhis2, varhis3
   real, allocatable, dimension(:, :, :, :) :: varforcing, varforcing2, varforcing5
   real, allocatable, dimension(:, :, :) :: varforcing4
   real, allocatable, dimension(:, :) :: varforcing3

   character*200 filename
integer :: hour,day,month,year,dayforc,monthforc,yearforc,lastday,daylai,monthlai,yearlai,jday,jdaylaipast,jdaylaifut,mdays,daypast,monthpast,yearpast,daysfromstart &
              , montho18, jdayo18past, jdayo18fut

   integer :: daysforoutput(8)
   integer :: hoursn, daysn, monthsn, yearsn
   real :: tfact, dtlr, swlat, swlon, midpoint
   real*8 :: t1, t2, t3, t4, t5
   integer :: mondays(12)
   data mondays/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
   character*3 monthname(12)
   data monthname/'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

   character*100 date
   integer :: irec, istart, request, requestsnow, req(nzg), req2(nzg), req3(nvar_out), rc, nsoil &
              , istartforcing, istartforcing2, istartforcing3, istartforcing4, istartforcing5
   integer :: istarthis, reqhis(nzg), reqhis2(0:nzg + 1), reqhis3(3:7)
   integer, allocatable, dimension(:) :: reqforcing, reqforcing2, reqforcing3, reqforcing4, reqforcing5

   character(len=300) :: filesoil, filetopo, filef, filewtd, fileveg, filehveg, filesmoieq, filerivers, fileo18, pathoutput

   swlat = -55.9958333333333
   swlon = -92.9958333333333

   filesoil = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/sa_soil_modified.dat'
   filef = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/fdepth_sa.nc'
   filetopo = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/RIVERPARAMETERS/riverparameters_definitive_yama.dat'
   filewtd = '/mnt/netapp2/Store_uscfmgmm/GLOBALWTD/sciencepaper/ncfileshaibin_newrun/S_America_model_wtd_v2.nc'
   fileveg = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/vegsam_modis.nc'
   filehveg = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/hveg_sam.nc'
   filesmoieq = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/smoieq_SA_new.nc'
   filerivers = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/RIVERPARAMETERS/riverparameters_definitive_yama.dat'
   fileo18 = '/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/modelwithERA5/isoratio_sa.nc'
   pathoutput = '/mnt/lustre/hsm/nlsas/notape/home/usc/fm/gfnlmeteo/GONZALO/O18RUNS/SAMERICA/outputparera5/'
!pathoutput='/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/outputparera5/'
!gmm intialize mpi stuff

   call MPI_INIT(ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *, 'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if

   call INITIALIZEDOMAIN(n2, n3, nzg, filetopo)

!Y dimension for each node

   NMAX_x = nend_x(pid) - nini_x(pid) + 1
   NMAX_y = nend_y(pid) - nini_y(pid) + 1
   is = nini_x(pid)
   ie = nend_x(pid)
   js = nini_y(pid)
   je = nend_y(pid)

!gmm allocate variables

   allocate (landmask(is:ie, js:je))
   allocate (veg(is:ie, js:je))
   allocate (hveg(is:ie, js:je))
   allocate (soiltxt(2, is:ie, js:je))
   allocate (smoi(nzg, is:ie, js:je))
   allocate (o18(nzg, is:ie, js:je))
   allocate (watext(nzg, is:ie, js:je))
   allocate (smoieq(nzg, is:ie, js:je))
!  allocate(smoimax(nzg,is:ie,js:je))
!  allocate(smoimin(nzg,is:ie,js:je))

   allocate (wtd(is:ie, js:je))
   allocate (smoiwtd(is:ie, js:je))
   allocate (fdepth(is:ie, js:je))
   allocate (rech(is:ie, js:je))
   allocate (deeprech(is:ie, js:je))
   allocate (qsrun(is:ie, js:je))
   allocate (qsrunsum(is:ie, js:je))
   allocate (qslat(is:ie, js:je))
   allocate (qlat(is:ie, js:je))
   allocate (qlato18(is:ie, js:je))
   allocate (qlatsum(is:ie, js:je))
   allocate (qlatinsum(is:ie, js:je))
   allocate (qlatoutsum(is:ie, js:je))
   allocate (qlatino18sum(is:ie, js:je))
   allocate (qlatouto18sum(is:ie, js:je))
   allocate (qlatino18ratio(is:ie, js:je))
   allocate (qlatouto18ratio(is:ie, js:je))
   allocate (qsprings(is:ie, js:je))
   allocate (et_s(is:ie, js:je))
   allocate (et_i(is:ie, js:je))
   allocate (et_c(is:ie, js:je))
   allocate (transpo18(is:ie, js:je))
   allocate (transpo18ratio(is:ie, js:je))
   allocate (intercepstore(is:ie, js:je))
   allocate (ppacum(is:ie, js:je))
   allocate (waterdeficit(is:ie, js:je))
   allocate (watextdeep(is:ie, js:je))
   allocate (pppendepth(is:ie, js:je))
   allocate (pppendepthold(is:ie, js:je))
   allocate (inactivedays(0:nzg + 1, is:ie, js:je))
   allocate (wtdmax(is:ie, js:je))
   allocate (wtdmin(is:ie, js:je))

   allocate (topo(is:ie, js:je))
   allocate (area(is:ie, js:je))
   allocate (lats(is:ie, js:je))
   allocate (lons(is:ie, js:je))

   allocate (riverflow(is:ie, js:je))
   allocate (qrf(is:ie, js:je))
   allocate (qrfsum(is:ie, js:je))
   allocate (delsfcwat(is:ie, js:je))
   allocate (delsfcwatsum(is:ie, js:je))
   allocate (slope(is:ie, js:je))
   allocate (riverdepth(is:ie, js:je))
   allocate (riverwidth(is:ie, js:je))
   allocate (riverlength(is:ie, js:je))
   allocate (maxdepth(is:ie, js:je))
   allocate (riverflowmean(is:ie, js:je))
   allocate (floodheight(is:ie, js:je))
   allocate (topoflood(is:ie, js:je))
   allocate (riverarea(is:ie, js:je))
   allocate (floodarea(is:ie, js:je))
   allocate (riverchannel(is:ie, js:je))
   allocate (fd(is:ie, js:je))
   allocate (bfd(is:ie, js:je))

   allocate (wind(is:ie, js:je))
   allocate (temp(is:ie, js:je))
   allocate (tempsfc(is:ie, js:je))
   allocate (qair(is:ie, js:je))
   allocate (netrad(is:ie, js:je))
   allocate (rshort(is:ie, js:je))
   allocate (press(is:ie, js:je))
   allocate (rain(is:ie, js:je))
   allocate (snow(is:ie, js:je))
   allocate (snow_past(is:ie, js:je))
   allocate (snow_fut(is:ie, js:je))
   allocate (precip(is:ie, js:je))
   allocate (lai(is:ie, js:je))
   allocate (lai_past(is:ie, js:je))
   allocate (lai_fut(is:ie, js:je))
   allocate (o18ratiopp(is:ie, js:je))
   allocate (o18_past(is:ie, js:je))
   allocate (o18_fut(is:ie, js:je))
   allocate (press_past(is:ie, js:je))
   allocate (press_fut(is:ie, js:je))
   allocate (temp_past(is:ie, js:je))
   allocate (temp_fut(is:ie, js:je))
   allocate (qair_past(is:ie, js:je))
   allocate (qair_fut(is:ie, js:je))
   allocate (wind_past(is:ie, js:je))
   allocate (wind_fut(is:ie, js:je))
   allocate (tempsfc_past(is:ie, js:je))
   allocate (tempsfc_fut(is:ie, js:je))
   allocate (icefactor(is:ie, js:je, nzg - 14:nzg))
   allocate (topoera5(is:ie, js:je))

   allocate (infilflux(nzg, is:ie, js:je))
   allocate (infilfluxday(nzg, is:ie, js:je))
   allocate (infilcounter(nzg, is:ie, js:je))
   allocate (upflux(nzg, is:ie, js:je))

   allocate (wtdflux(is:ie, js:je, 8))
   allocate (et_s_daily(is:ie, js:je, 8))
   allocate (et_c_daily(is:ie, js:je, 8))
   allocate (transptop(is:ie, js:je, 8))
   allocate (infilk(is:ie, js:je, 8))
   allocate (wtd_daily(is:ie, js:je, 8))
   allocate (smoi_daily(is:ie, js:je, 8))

   if (writepar .eq. 0) then
      allocate (varoutput(is:ie, js:je, nzg))
      allocate (varoutput2(is:ie, js:je, nzg))
      if (freedrain .eq. 0) then
         allocate (varoutput3(is:ie, js:je, nvar_out))
      else
         allocate (varoutput3(is:ie, js:je, 5))
      end if

      allocate (varhis(is:ie, js:je, nzg))
      allocate (varhis2(is:ie, js:je, 0:nzg + 1))
      allocate (varhis3(is:ie, js:je, 3:7))
   end if

   allocate (reqforcing(numtasks - 2))
   allocate (varforcing(dimerax, dimeray, 0:23, 11))
!allocate(reqforcing2(numtasks-2))
!allocate(varforcing2(xdim,ydim,0:23,3))
!allocate(reqforcing3(numtasks-2))
!allocate(varforcing3(dimerax,dimeray))
   allocate (reqforcing4(numtasks - 2))
   allocate (varforcing4(dimeralandx, dimeralandy, 0:23))
!allocate(reqforcing5(numtasks-2))
!allocate(varforcing5(xdim,ydim,0:23,4))

!timestep for rivers
   niter_river = max(1, nint(deltat/deltatriver + 0.4))
!        niter_river = 1
   dtlr = deltat/float(niter_river)

   if (pid .eq. 0) write (6, *) 'niter_river', niter_river, dtlr

!counter to tell the output routine that it is first time

   istart = 1
   istarthis = 1
   istartforcing = 1
   istartforcing2 = 1
   istartforcing3 = 1
   istartforcing4 = 1
   istartforcing5 = 1
   irec = 0

!initialize some soil parameters

   call init_soil_param(fieldcp, nzg)

   if (pid .eq. 0) write (6, *) 'slwilt', slwilt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!read in fixed fields

   if (pid .eq. 0) write (6, *) 'reading soiltxt,topo and landmask'

   call READINITIAL(n2, n3, is, ie, js, je, soiltxt, topo, fdepth, landmask, filesoil, filetopo, filef)

!soiltxt=13

   call READVEG(n2, n3, is, ie, js, je, veg, fileveg)

   call READHVEG(n2, n3, is, ie, js, je, hveg, filehveg)

   call READLATLON(n2, n3, is, ie, js, je, lats, lons, area, dxy, swlat, swlon)

   if (pid .eq. 1) write (6, *) 'SW corner of the domain', lats(is, js), lons(is, js)

   call READTOPOERA5(n2, n3, is, ie, js, je, landmask, lats, lons, topoera5)

!initialize variables

   if (riverswitch .eq. 1) then
      call READFLOWDIRECTION(n2, n3, is, ie, js, je, fd, bfd, filerivers)
      call READRIVERPARAMETERS(n2, n3, is, ie, js, je, riverlength, filerivers, 2)
      call READRIVERPARAMETERS(n2, n3, is, ie, js, je, riverwidth, filerivers, 4)
      call READRIVERPARAMETERS(n2, n3, is, ie, js, je, slope, filerivers, 5)
      call READRIVERPARAMETERS(n2, n3, is, ie, js, je, topoflood, filerivers, 9)
      call READRIVERPARAMETERS(n2, n3, is, ie, js, je, maxdepth, filerivers, 11)
      riverarea = riverwidth*riverlength
      floodarea = max(area - riverarea, 0.)
      riverchannel = maxdepth*riverarea
   end if

!        open(56,file='testinput.dat' &
!            ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)

!                          write(56,rec=1) ((topoflood(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=2) ((float(soiltxt(1,i,j)),i=1,n2),j=1,n3)
!                          write(56,rec=3) ((float(soiltxt(2,i,j)),i=1,n2),j=1,n3)
!                          write(56,rec=4) ((fdepth(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=5) ((veg(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=6) ((hveg(i,j),i=1,n2),j=1,n3)
!           close(56)
!           write(6,*)'done with test'

!stop

   if (pid .eq. 0) write (6, *) 'now to initialize soil depth and initial data'

   call INITIALIZESOILDEPTH(nzg, slz, dz)

   if (pid .eq. 0) write (6, *) 'soil layers', (slz(k), k=1, nzg + 1), (dz(k), k=1, nzg)

   if (pid .eq. 0) write (6, *) 'soil nodes', (-0.5*(slz(k) + slz(k + 1)), k=1, nzg)

   pppendepthold = nzg + 1
   wtdflux = 0.
   et_s_daily = 0.
   et_c_daily = 0.
   transptop = 0.
   infilk = nzg + 1
   infilflux = 0.
   infilfluxday = 0.
   infilcounter = 0
   upflux = 0.
!       o18 = 2005.20 * 1.e-6 * smoi
!       smoimax=0.
!       smoimin=1.

   daypast = 1 !5day period number for output

   if (restart .eq. 1) then
!                filename='/mnt/gluster/distributed/home/usc/fm/gmm/SAMERICA/outputSA/history_wt_01aug2013.nc'
!                filename='/mnt/netapp2/Store_uscfmgmm/GLUSTER/SAMERICA/outputSAinitial/history_wt_01jan2014.nc'
!               filename='/mnt/netapp2/Store_uscfmgmm/GLUSTER/SAMERICA/outputera5SAini/history_wt_01jan2014.nc'
!                filename='/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/outputera5SA/history_wt_01jan2014.nc'
                filename='/mnt/lustre/hsm/nlsas/notape/home/usc/fm/gfnlmeteo/GONZALO/O18RUNS/SAMERICA/outputera5SA/historyfull/history_wt_01jan2023.nc'
      call READHISTORYNC(n2, n3, is, ie, js, je, nzg, smoi, smoiwtd, intercepstore, wtd, inactivedays, filename)
      if (freedrain .eq. 0) then
!                 call READSMOIEQ(n2,n3,nzg,js,je,smoieq,filesmoieq)
         call EQSOILMOISTUREtheor(is, ie, js, je, nzg, slz, dz, soiltxt, landmask, fdepth, smoieq)
      end if
      if (riverswitch .eq. 1) then
         call READHISTORYVARNC(n2, n3, is, ie, js, je, riverflow, 'RIVERFLOW', filename)
         call READHISTORYVARNC(n2, n3, is, ie, js, je, riverdepth, 'RIVERDEPTH', filename)
         call READHISTORYVARNC(n2, n3, is, ie, js, je, floodheight, 'FLOODHEIGHT', filename)
      end if

!careful with this when starting from the very beginning
      call READHISTORYVAR3DNC(n2, n3, is, ie, js, je, nzg, o18, 'O18', filename)

!                  call READHISTORYBYTEVARNC(n2,n3,is,ie,js,je,pppendepthold,'PPPENDEPTHOLD',filename)

!                  call READHISTORYVARNC(n2,n3,is,ie,js,je,wtdflux(1,js,1),'WTDFLUX',filename)
!                  call READHISTORYVARNC(n2,n3,is,ie,js,je,et_s_daily(1,js,1),'ET_S_DAILY',filename)
!                  call READHISTORYVARNC(n2,n3,is,ie,js,je,et_c_daily(1,js,1),'ET_C_DAILY',filename)
!                  call READHISTORYVARNC(n2,n3,is,ie,js,je,transptop(1,js,1),'TRANSPTOP',filename)
!                  call READHISTORYBYTEVARNC(n2,n3,is,ie,js,je,infilk(1,js,1),'INFILK',filename)

      where (veg .lt. 0.5) landmask = 0

!!!CHANGE THIS
      goto 222
      inactivedays = maxinactivedays + 1
      inactivedays(nzg + 1, is:ie, js:je) = 0

      riverdepth = max(riverdepth, 0.)
      where (riverwidth .lt. 1.00001) riverdepth = max(min(riverdepth, 5.), floodheight)
      where (riverwidth .lt. 0.5) riverdepth = floodheight
      where (riverwidth .lt. 0.5) riverflow = 0.

!comment this out when not starting from the very beginning!!!initilize O18 ratio from climatology
      call READO18CLIM(n2, n3, is, ie, js, je, 13, fileo18, o18_past)

      do j = js, je
      do i = is, ie
      do k = 1, nzg
         if (landmask(i, j) .eq. 0) then
            if (slz(k) .lt. -0.3) then
               nsoil = soiltxt(1, i, j)
            else
               nsoil = soiltxt(2, i, j)
            end if
            midpoint = 0.5*(slz(k) + slz(k + 1))
            smoi(k, i, j) = slmsts(nsoil)*max(min(exp((midpoint + 1.5)/fdepth(i, j)), 1.), 0.1)
!!                smoimax(k,i,j)=smoi(k,i,j)
!!                smoimin(k,i,j)=smoi(k,i,j)
         end if

         o18(k, i, j) = o18_past(i, j)*smoi(k, i, j)
      end do
      end do
      end do

222   continue
!!!!

   else

      if (freedrain .eq. 0) then
         call READWTDNC(n2, n3, is, ie, js, je, wtd, filewtd)
         where (wtd .lt. -1.e5) wtd = 0.
         wtd = max(wtd, slz(1))
!                 call READSMOIEQ(n2,n3,nzg,js,je,smoieq,filesmoieq)
!                 call EQSOILMOISTURE(is,ie,js,je,nzg,slz,dz,deltat,soiltxt,landmask,smoieq)
         call EQSOILMOISTUREtheor(is, ie, js, je, nzg, slz, dz, soiltxt, landmask, fdepth, smoieq)
      end if

      where (veg .lt. 0.5) landmask = 0
      where (landmask .eq. 0) wtd = 0.
      call INITIALIZE(n2, n3, is, ie, js, je, nzg, freedrain, slz, dz, soiltxt, wtd, smoi, smoieq &
                      , fdepth, topo, landmask, deltat, area)
      intercepstore = 0.

      inactivedays = maxinactivedays + 1
      inactivedays(nzg + 1, is:ie, js:je) = 0

!initialize O18 ratio from climatology
      call READO18CLIM(n2, n3, is, ie, js, je, 13, fileo18, o18_past)
      do k = 1, nzg
         o18(k, :, :) = o18_past(:, :)*smoi(k, :, :)
      end do
!!!
      floodheight = 0.

      if (riverswitch .eq. 1) then
         call READRIVERPARAMETERS(n2, n3, is, ie, js, je, riverflow, filerivers, 6)
         call READRIVERPARAMETERS(n2, n3, is, ie, js, je, riverdepth, filerivers, 3)
      end if
   end if

   if (pid .eq. 1) write (6, *) 'smoi', smoi(5, 5, 5), landmask(5, 5), smoieq(5, 5, 5)

   delsfcwat = 0.
   delsfcwatsum = 0.
   qrf = 0.
   qrfsum = 0.
   qsrun = 0.
   qsrunsum = 0.
   qsprings = 0.
   rech = 0.
   deeprech = 0.
   qslat = 0.
   qlat = 0.
   qlato18 = 0.
   qlatsum = 0.
   qlatinsum = 0.
   qlatoutsum = 0.
   qlatino18sum = 0.
   qlatouto18sum = 0.
   et_s = 0.
   et_i = 0.
   et_c = 0.
   transpo18 = 0.
   ppacum = 0.
   watext = 0.
   watextdeep = 0.
   waterdeficit = 0.
   pppendepth = 0.
   wind_past = 0.
   wind_fut = 0.
   temp_past = 0.
   temp_fut = 0.
   tempsfc_past = 0.
   tempsfc_fut = 0.
   press_past = 0.
   press_fut = 0.
   qair_past = 0.
   qair_fut = 0.
   rain = 0.
   snow = 0.
   snow_past = 0.
   snow_fut = 0.
   rshort = 0.
   netrad = 0.
   wtdmax = -1000.
   wtdmin = 0.

!initial time

   hour = 0
   day = 1
   month = 1
   year = 2023

   daysfromstart = daynumber(day, month, year)

   daylai = 1
   monthlai = 1

!read past temp,wind,qair,pres

!       icount = (day-1)*8 + hour/3+1
   icount = (day - 1)*24 + hour + 1

   write (filename, '(i4.4,a1,i2.2,a3)') year, '-', month, '.nc'

   if (pid .eq. 0) write (6, *) 'reading now forcings first ', filename

   call READFORCINGS(is,ie,js,je,filename,icount,hour,topo,lats,lons,wind_past,temp_past,press_past,qair_past,tempsfc_past,landmask&
                     , topoera5, varforcing, request, istartforcing)

!       call READFORCINGSACC(is,ie,js,je,filename,icount,hour,lats,lons,rain,rshort,netrad,landmask &
!                        ,varforcing,reqforcing,istartforcing)

   daysn = day
   hoursn = hour
   yearsn = year
   monthsn = month

!          icount = (daysn-1)*8 + hoursn/3 +1

   write (filename, '(i4.4,a1,i2.2,a3)') yearsn, '-', monthsn, '.nc'

   call READFORCINGSSNOW(is, ie, js, je, filename, icount, hoursn, lats, lons, snow_fut, landmask &
                         , varforcing4, requestsnow, istartforcing4)

!read lai. the first time of each file is always 1 jan

!       call READLAI(is,ie,js,je,lats,lons,year,monthlai,daylai,lai)

!read lai

   jday = julday(month, day, year)

   jdaylaipast = jday - mod(jday - 1, 8)
   yearlai = year

   write (filename, '(i4.4,a13,i4.4,a1,i3.3,a3)') yearlai, '/SAhires_30s_', yearlai, '_', jdaylaipast, '.nc'
   filename = '/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

   call READLAICHINA(n2, n3, is, ie, js, je, filename, lai_past)

   jdaylaifut = jdaylaipast + 8 !this works because initial time is always the first day of a month

   write (filename, '(i4.4,a13,i4.4,a1,i3.3,a3)') yearlai, '/SAhires_30s_', yearlai, '_', jdaylaifut, '.nc'
   filename = '/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

   call READLAICHINA(n2, n3, is, ie, js, je, filename, lai_fut)

   tfact = float(jday - jdaylaipast)/float(jdaylaifut - jdaylaipast)
   lai = lai_past + (lai_fut - lai_past)*tfact

   if (pid .eq. 0) write (6, *) 'read lai', jday, jdaylaipast, jdaylaifut

!read isotopic ratio of o18 in precipitation

   jday = julday(month, day, year)
   mdays = mondays(month)
   if (month .eq. 2 .and. mod(year, 4) .eq. 0.) mdays = 29

   if (day .gt. mdays/2) then
      montho18 = month
      call READO18CLIM(n2, n3, is, ie, js, je, montho18, fileo18, o18_past)
      jdayo18past = julday(montho18, mdays/2, year)

      montho18 = month + 1
      if (montho18 .eq. 13) then
         montho18 = 1
         jdayo18fut = julday(12, 31, year) + julday(montho18, mondays(1)/2, year + 1)
      elseif (montho18 .eq. 2 .and. mod(year, 4) .eq. 0.) then
         jdayo18fut = julday(montho18, 29/2, year)
      else
         jdayo18fut = julday(montho18, mondays(montho18)/2, year)
      end if

      call READO18CLIM(n2, n3, is, ie, js, je, montho18, fileo18, o18_fut)

   else

      montho18 = month
      call READO18CLIM(n2, n3, is, ie, js, je, montho18, fileo18, o18_fut)
      jdayo18fut = julday(montho18, mdays/2, year)

      montho18 = month - 1
      if (montho18 .eq. 0) then
         montho18 = 12
         jdayo18past = -17
      elseif (montho18 .eq. 2 .and. mod(year, 4) .eq. 0.) then
         jdayo18past = julday(montho18, 29/2, year)
      else
         jdayo18past = julday(montho18, mondays(montho18)/2, year)
      end if

      call READO18CLIM(n2, n3, is, ie, js, je, montho18, fileo18, o18_past)

   end if

   tfact = float(jday - jdayo18past)/float(jdayo18fut - jdayo18past)
   o18ratiopp = o18_past + (o18_fut - o18_past)*tfact

   if (pid .eq. 0) write (6, *) 'read o18', jday, jdayo18past, jdayo18fut

!write initial state
   if (restart .eq. 0) then
      if (freedrain .eq. 0) then

         write (date, '(i4,a1,i2.2,a1,i2.2,a9)') year, '-', month, '-', day, '_00:00:00'

         if (writepar .eq. 1) then

            write (filename, '(a13,i2.2,a3,i4.4,a1,i3.3,a3)') 'rootdaily_wt_', day, monthname(month), year, '_', pid, '.nc'
            filename = pathoutput(1:len_trim(pathoutput))//filename
            call writeoutputnc_par(n2, n3, is, ie, js, je, nzg, dz, smoi, waterdeficit, watext, watextdeep &
                                   , wtd, smoiwtd, qsrunsum, rech, qsprings, qlatsum, et_s, et_i, et_c, ppacum, pppendepth &
                                   , riverflow, qrfsum, delsfcwatsum &
                                   , filename, irec, istart, req, req2, req3, date, nvar_out)
         else

            write (filename, '(a13,i2.2,a3,i4.4,a3)') 'rootdaily_wt_', day, monthname(month), year, '.nc'
            filename = '/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/output/'//filename

            call writeoutputnc(n2, n3, is, ie, js, je, nzg, dz, smoi, waterdeficit, watext, watextdeep &
                               , wtd, smoiwtd, qsrunsum, rech, qsprings, qlatsum, et_s, et_i, et_c, ppacum, pppendepth &
                               , riverflow, qrfsum, delsfcwatsum &
                               , lai, filename, irec, istart, varoutput, varoutput2, varoutput3, req, req2, req3, date, nvar_out)
         end if
      else
!                         write(filename,'(a14,i2.2,a4)')'smoiwatext_fd_',month-1,'.dat'
!                         if(year.eq.1982)filename='smoiwatext_fd_12.dat'
         write (filename, '(a13,i2.2,a3,i4.4,a6)') 'rootdaily_fd_', day, monthname(month), year, '.grads'
         filename = 'outputfd/'//filename

         call writeoutputfd(n2, n3, is, ie, js, je, nzg, smoi, waterdeficit, watext &
                            , qsrunsum, rech, et_c, ppacum &
                            , filename, irec, istart, req, req2, req3)
      end if
   end if

!!!!!! do something if restarting on day day different than jan 1

   DO WHILE (year .ne. 2024)
!DO WHILE(month.ne.2)

!read icefactor before updating time
      icount = (day - 1)*24 + hour + 1

      write (filename, '(i4.4,a1,i2.2,a3)') year, '-', month, '.nc'

!       icount = (julday(month,day,year)-1)*8 + hour/3+1
!       write(filename,'(i4.4)')year
      call READFORCINGSSOILT(is, ie, js, je, nzg, filename, icount, hour, lats, lons, icefactor, landmask, topo, topoera5 &
                             , varforcing, reqforcing, istartforcing)

!advance time to read forcing
      monthpast = month
      yearpast = year
      hour = hour + 1
      if (hour .ge. 24) then
         hour = hour - 24
         day = day + 1
         lastday = mondays(month)
         if (month .eq. 2 .and. mod(year, 4) .eq. 0) lastday = 29
         if (day .eq. lastday + 1) then
            day = 1
            month = month + 1
            if (month .eq. 13) then
               month = 1
               year = year + 1
            end if
         end if
         daysfromstart = daynumber(day, month, year)
      end if

      if (pid .eq. 0) write (6, *) 'Hour, day, month, year', hour, day, month, year, lastday, daysfromstart

      t1 = mpi_wtime()

      icount = (day - 1)*24 + hour + 1

      write (filename, '(i4.4,a1,i2.2,a3)') year, '-', month, '.nc'

      call READFORCINGS(is,ie,js,je,filename,icount,hour,topo,lats,lons,wind_fut,temp_fut,press_fut,qair_fut,tempsfc_fut,landmask &
                        , topoera5, varforcing, request, istartforcing)

      press = 0.5*(press_past + press_fut)
      wind = 0.5*(wind_past + wind_fut)
      temp = 0.5*(temp_past + temp_fut)
      qair = 0.5*(qair_past + qair_fut)
      tempsfc = 0.5*(tempsfc_past + tempsfc_fut)

!move forcings to past
      press_past = press_fut
      wind_past = wind_fut
      temp_past = temp_fut
      qair_past = qair_fut
      tempsfc_past = tempsfc_fut

      call READFORCINGSACC(is, ie, js, je, filename, icount, hour, lats, lons, rain, rshort, netrad, landmask &
                           , varforcing, reqforcing, istartforcing)

!read snowmelt for the next 3h

!      hoursn = hour + 2 ! so, add 2 + 1 = 3h from present time

!     if(mod(hoursn,3).eq.0)then

!      daysn = day
!      monthsn = month
!      yearsn = year

!      if(hoursn.ge.24)then
!         hoursn=hoursn-24
!         daysn=day+1
!            lastday=mondays(month)
!            if(month.eq.2.and.mod(year,4).eq.0)lastday=29
!            if(daysn.eq.lastday+1)then
!               daysn=1
!               monthsn=month+1
!               if(monthsn.eq.13)then
!                    monthsn=1
!                    yearsn=year+1
!               endif
!             endif
!       endif

!          icount = (daysn-1)*8 + hoursn/3 +1

!now snow comes every hour
      yearsn = year
      monthsn = month
      hoursn = hour

      write (filename, '(i4.4,a1,i2.2,a3)') yearsn, '-', monthsn, '.nc'

      call READFORCINGSSNOW(is, ie, js, je, filename, icount, hoursn, lats, lons, snow_fut, landmask &
                            , varforcing4, requestsnow, istartforcing4)

!          snow = ( snow_fut - snow_past ) * (deltat/(3.*3600.))
      snow = (snow_fut - snow_past)*(deltat/3600.)

      if (hoursn .eq. 0) then
         snow_past = 0.
      else
         snow_past = snow_fut
      end if

!       endif

      precip = rain + snow

!        open(56,file='testinput.dat' &
!           ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)

!                          write(56,rec=1) ((wind(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=2) ((temp(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=3) ((press(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=4) ((qair(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=5) ((rshort(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=6) ((netrad(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=7) ((rain(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=8) ((topoera5(i,j),i=1,n2),j=1,n3)
!           close(56)
!stop
!           write(6,*)'done with test'

!write(6,*)press(100,100),wind(100,100),temp(100,100),qair(100,100),rshort(100,100),netrad(100,100),lai(100,100),veg(100,100)

      t2 = mpi_wtime()

!now run the model
!       if(freedrain.eq.0.and.hour.eq.nint(deltat/3600.)) then
      call LATERAL(n2, n3, is, ie, js, je, nzg, soiltxt, wtd, qlat, fdepth, topo, landmask, deltatwtd, area, lats, dxy, slz &
                   , o18, smoi, qlato18, qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum)
      qslat = 0.
      qlatsum = qlatsum + qlat*1.e3
!       endif

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

      t4 = mpi_wtime()

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

      t5 = mpi_wtime()

      if (pid .eq. numtasks/2) write (6, '(a12,5f7.3)') 'CPU(sec)', T5 - T1, t2 - t1, t3 - t2, t4 - t3, t5 - t4

      inpair(1) = T5 - T1
      inpair(2) = float(pid)

      call MPI_REDUCE(inpair, outpair, 1, MPI_2REAL, MPI_MAXLOC, numtasks - 1, MPI_COMM_WORLD, ierr)

      if (pid .eq. numtasks - 1) write (6, '(a27,f7.3,f7.0)') 'max CPU total step, pid', outpair(1), outpair(2)
!if(pid.eq.1)write(6,'(a12,5f7.3)')'pid 1 CPU(sec)',T5-T1,t2-t1,t3-t2,t4-t3,t5-t4

!output

!        if(year.ge.1981)then
!         if(day.eq.2)then

      if (day .eq. 1 .and. hour .eq. 0) then

         riverflowmean = riverflowmean/(float(lastday)*86400.)
         delsfcwatsum = delsfcwatsum/float(lastday)
         qrfsum = qrfsum/float(lastday)
         qsrunsum = qsrunsum/float(lastday)
         rech = rech/float(lastday)
         qsprings = qsprings/float(lastday)
         qlatsum = qlatsum/float(lastday)
         where (qlatinsum .gt. 0.)
            qlatino18ratio = qlatino18sum/qlatinsum
         elsewhere
            qlatino18ratio = 0.
         end where
         where (qlatoutsum .gt. 0.)
            qlatouto18ratio = qlatouto18sum/qlatoutsum
         elsewhere
            qlatouto18ratio = 0.
         end where
         qlatinsum = qlatinsum/float(lastday)
         qlatoutsum = qlatoutsum/float(lastday)
         where (et_c .gt. 0.)
            transpo18ratio = transpo18/et_c
         elsewhere
            transpo18ratio = 0.
         end where
         et_s = et_s/float(lastday)
         et_i = et_i/float(lastday)
         et_c = et_c/float(lastday)
         ppacum = ppacum/float(lastday)
         waterdeficit = waterdeficit/float(lastday)
         watext = watext/float(lastday)
         watextdeep = watextdeep/float(lastday)
         infilflux = infilflux/float(lastday)
         upflux = upflux/float(lastday)

         monthhours = lastday*24

!              irec=0
!          endif
!              if(month.gt.1.or.year.eq.1983)then
         if (freedrain .eq. 0) then

            write (date, '(i4,a1,i2.2,a1,i2.2,a9)') year, '-', month, '-', day, '_00:00:00'

            if (writepar .eq. 1) then

               write (filename, '(a13,i2.2,a3,i4.4,a1,i3.3,a3)') 'rootdaily_wt_', day, monthname(month), year, '_', pid, '.nc'
               filename = pathoutput(1:len_trim(pathoutput))//filename

               call writeoutputnc_par(n2, n3, is, ie, js, je, nzg, dz, smoi, waterdeficit, watext, watextdeep &
                                      , wtd, smoiwtd, qsrunsum, rech, qsprings, qlatsum, et_s, et_i, et_c, ppacum, pppendepth &
                                      , riverflowmean, qrfsum, delsfcwatsum &
                                      !                                   ,wind,temp,qair,press,netrad,rshort,pet,lai  &
                                      , filename, irec, istart, req, req2, req3, date, nvar_out)

!                   write(filename,'(a12,a3,i4.4,a1,i3.3,a3)')'dailyoutput_',monthname(monthpast),yearpast,'_',pid,'.nc'
!                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputparera5/'//filename

!                         call writeoutputnc_daily_par(n2,n3,is,ie,js,je,nzg,daypast-1 &
!                                             ,wtdflux,et_s_daily,et_c_daily,transptop,infilk,smoi_daily,wtd_daily &
!                                             ,daysforoutput,filename )

               write (filename, '(a12,a3,i4.4,a1,i3.3,a3)') 'infiloutput_', monthname(monthpast), yearpast, '_', pid, '.nc'
               filename = pathoutput(1:len_trim(pathoutput))//filename

               call writeoutputnc_infil_par(n2, n3, is, ie, js, je, nzg, lastday &
                                            , infilflux, o18, smoi, transpo18ratio, upflux, wtdmax, wtdmin &
                                            , qlatinsum, qlatoutsum, qlatino18ratio, qlatouto18ratio &
                                            , filename)

            else

               write (filename, '(a13,i2.2,a3,i4.4,a3)') 'rootdaily_wt_', day, monthname(month), year, '.nc'
               filename = '/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/output/'//filename

               call writeoutputnc(n2, n3, is, ie, js, je, nzg, dz, smoi, waterdeficit, watext, watextdeep &
                                  , wtd, smoiwtd, qsrunsum, rech, qsprings, qlatsum, et_s, et_i, et_c, ppacum, pppendepth &
                                  , riverflowmean, qrfsum, delsfcwatsum &
                                  , lai, filename, irec, istart, varoutput, varoutput2, varoutput3, req, req2, req3, date, nvar_out)

!                   write(filename,'(a12,a3,i4.4,a3)')'dailyoutput_',monthname(monthpast),yearpast,'.nc'
!                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/output/'//filename

!                         call writeoutputnc_daily_par(n2,n3,is,ie,js,je,nzg,daypast-1 &
!                                             ,wtdflux,et_s_daily,et_c_daily,transptop,infilk,smoi_daily,wtd_daily &
!                                             ,daysforoutput,filename )

               write (filename, '(a12,a3,i4.4,a3)') 'infiloutput_', monthname(monthpast), yearpast, '.nc'
               filename = '/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/output/'//filename

               call writeoutputnc_infil_par(n2, n3, is, ie, js, je, nzg, lastday &
                                            , infilflux, o18, smoi, transpo18ratio, upflux, wtdmax, wtdmin &
                                            , qlatinsum, qlatoutsum, qlatino18ratio, qlatouto18ratio &
                                            , filename)

            end if

         else
!                         write(filename,'(a14,i2.2,a4)')'smoiwatext_fd_',month-1,'.dat'
!                         if(year.eq.1982)filename='smoiwatext_fd_12.dat'
            write (filename, '(a13,i2.2,a3,i4.4,a6)') 'rootdaily_fd_', day, monthname(month), year, '.grads'
            filename = 'outputfd/'//filename

            call writeoutputfd(n2, n3, is, ie, js, je, nzg, smoi, waterdeficit, watext &
                               , qsrunsum, rech, et_c, ppacum &
                               , filename, irec, istart, req, req2, req3)
         end if
!               endif

         riverflowmean = 0.
         delsfcwatsum = 0.
         qrfsum = 0.
         qsrunsum = 0.
         rech = 0.
         qsprings = 0.
         qlatsum = 0.
         qlatinsum = 0.
         qlatoutsum = 0.
         qlatino18sum = 0.
         qlatouto18sum = 0.
         transpo18 = 0.
         et_s = 0.
         et_i = 0.
         et_c = 0.
         ppacum = 0.
         waterdeficit = 0.
         watext = 0.
         watextdeep = 0.
         pppendepth = 0.

         wtdflux(:, :, 1:daypast - 1) = 0.
         et_s_daily(:, :, 1:daypast - 1) = 0.
         et_c_daily(:, :, 1:daypast - 1) = 0.
         transptop(:, :, 1:daypast - 1) = 0.
         infilk(:, :, 1:daypast - 1) = nzg + 1

         infilflux = 0.
         infilcounter = 0
         upflux = 0.
         wtdmax = -1000.
         wtdmin = 0.

      end if
!        endif
!       endif

!history

      if (day .eq. 1 .and. hour .eq. 0) then

         write (date, '(i4,a1,i2.2,a1,i2.2,a9)') year, '-', month, '-', day, '_00:00:00'

         if (writepar .eq. 1) then

            write (filename, '(a11,i2.2,a3,i4.4,a1,i3.3,a3)') 'history_wt_', day, monthname(month), year, '_', pid, '.nc'
            filename = pathoutput(1:len_trim(pathoutput))//filename
            if (pid .eq. 1) write (6, *) 'writing history file', filename, date, daypast

            call writehistorync_par(n2, n3, is, ie, js, je, nzg, smoi, intercepstore, wtd, inactivedays &
                                    , riverflow, riverdepth, floodheight &
                                    , o18 &
                                    !                                          ,pppendepthold &
                                    !                                          ,wtdflux(1,js,daypast),et_s_daily(1,js,daypast),et_c_daily(1,js,daypast),transptop(1,js,daypast),infilk(1,js,daypast) &
                                    !                                           ,riverflow,netrad,floodheight &
                                    , filename, date)
         else

            write (filename, '(a11,i2.2,a3,i4.4,a3)') 'history_wt_', day, monthname(month), year, '.nc'
            filename = '/mnt/lustre/scratch/nlsas/home/usc/fm/gmm/SAMERICA/output/'//filename
            if (pid .eq. 1) write (6, *) 'writing history file', filename, date

            call writehistorync(n2, n3, is, ie, js, je, nzg, smoi, smoiwtd, intercepstore, wtd, inactivedays &
                                , riverflow, riverdepth, floodheight &
                                , filename, date, istarthis, varhis, varhis2, varhis3, reqhis, reqhis2, reqhis3)

         end if

!               wtdflux(:,:,1)=wtdflux(:,:,daypast)
!               et_s_daily(:,:,1)=et_s_daily(:,:,daypast)
!               et_c_daily(:,:,1)=et_c_daily(:,:,daypast)
!               transptop(:,:,1)=transptop(:,:,daypast)
!               infilk(:,:,1)=infilk(:,:,daypast)

!               wtdflux(:,:,daypast)=0.
!               et_s_daily(:,:,daypast)=0.
!               et_c_daily(:,:,daypast)=0.
!               transptop(:,:,daypast)=0.
!               infilk(:,:,daypast)=nzg+1

         daypast = 1  !bring back counter to 1

      end if

!now LAI from MODIS
!       if(mod( julday(month,day,year) -1 , 4 ) .eq. 0 .and.hour.eq.0)then
!                 if(pid.eq.0)write(6,*)'time to read LAI',julday(month,day,year)
!                 call READLAI(is,ie,js,je,lats,lons,year,month,day,lai)
!       endif

      if (hour .eq. 0) then

         jday = julday(month, day, year)
         if (jday .eq. 1) jdaylaifut = 1

         if (jday .eq. jdaylaifut) then

            jdaylaipast = jdaylaifut
            jdaylaifut = jday + 8
            yearlai = year

            lastday = 365
            if (mod(year, 4) .eq. 0) lastday = 366

            if (jdaylaifut .gt. lastday) then
               jdaylaifut = lastday + 1
               yearlai = year + 1
               write (filename, '(i4.4,a13,i4.4,a1,i3.3,a3)') yearlai, '/SAhires_30s_', yearlai, '_', 1, '.nc'
            else
               write (filename, '(i4.4,a13,i4.4,a1,i3.3,a3)') yearlai, '/SAhires_30s_', yearlai, '_', jdaylaifut, '.nc'
            end if

            filename = '/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

            lai_past = lai_fut
            call READLAICHINA(n2, n3, is, ie, js, je, filename, lai_fut)

            if (pid .eq. 0) write (6, *) 'jdaylaipast,jdaylaifut,yearlai', jdaylaipast, jdaylaifut, yearlai

         end if

         tfact = float(jday - jdaylaipast)/float(jdaylaifut - jdaylaipast)
         lai = lai_past + (lai_fut - lai_past)*tfact

!now O18 ratio in precipitation

         mdays = mondays(month)
         if (month .eq. 2 .and. mod(year, 4) .eq. 0.) mdays = 29

         if (day .eq. mdays/2 + 1) then

            o18_past = o18_fut
            jdayo18past = jdayo18fut

            montho18 = month + 1
            if (montho18 .eq. 13) then
               montho18 = 1
               jdayo18fut = julday(12, 31, year) + julday(montho18, mondays(1)/2, year + 1)
            elseif (montho18 .eq. 2 .and. mod(year, 4) .eq. 0.) then
               jdayo18fut = julday(montho18, 29/2, year)
            else
               jdayo18fut = julday(montho18, mondays(montho18)/2, year)
            end if

            call READO18CLIM(n2, n3, is, ie, js, je, montho18, fileo18, o18_fut)

            write (6, *) 'read o18', jday, jdayo18past, jdayo18fut

         elseif (jday .le. 15) then
            jdayo18past = -17
            jdayo18fut = 15
         end if

         tfact = float(jday - jdayo18past)/float(jdayo18fut - jdayo18past)
         o18ratiopp = o18_past + (o18_fut - o18_past)*tfact

      end if

   END DO

   call MPI_FINALIZE(ierr)

end
