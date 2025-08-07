MODULE module_io

use module_parallel

implicit none

CONTAINS

!******************************************************************************************************
subroutine READINITIAL(n2,n3,is,ie,js,je,soiltxt,topo,fdepth,landmask,filesoil,filetopo,filef)
use netcdf

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
integer, dimension(2,is:ie,js:je) :: soiltxt
integer, dimension(is:ie,js:je) :: landmask
real, dimension(is:ie,js:je) :: var,topo,fdepth

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request
integer :: ncid,varid,statusnc

character (len = *) :: filesoil,filetopo,filef

real*4, allocatable :: varreadbig(:,:)
real, allocatable :: varread(:,:)

!gmmread in soil textures

do k=1,2

if(pid.eq.0)then

write(6,*)'ready to read soil textures ',k

     iun=42
      open(iun,file=filesoil &
!      open(iun,file='/home/usc/fm/gmm/mnt/store/GLOBALWTD/sciencepaper/sa_soil.dat'
!      &
      ,form='unformatted',convert='little_endian',access='direct',recl=4*n2big*n3big)

     irec=0

     allocate(varreadbig(n2big,n3big))
     allocate(varread(n2,n3))

write(6,*)'reading soil textures ',k

     irec=irec+1
     irec=1  !for now
     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

!      statusnc = nf90_open(filesoil, 0, ncid)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      statusnc = nf90_inq_varid(ncid,'STXT',varid)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      statusnc = nf90_get_var(ncid,varid,varreadbig)!,start,count)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      statusnc = nf90_close(ncid)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)



     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

write(6,*)'sending soil textures to nodes',k

         do n=1,numtasks-2
         call MPI_isend(varread(1,1),1,arraysection(n),n,50+k,MPI_COMM_WORLD,req(n),ierr)
         enddo

     if(numtasks.eq.1)then

       do j=js,je
            do i=is,ie
               soiltxt(k,i,j)=nint(varread(i,j))
            enddo
          enddo
      else
           call MPI_waitall(numtasks-2,req,stats,ierr)
      endif

      deallocate(varread,varreadbig)
      close(iun)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,50+k,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      soiltxt(k,:,:)=nint(var(:,:))

endif

enddo


where(soiltxt.eq.0)soiltxt=5

!gmmread in topo

if(pid.eq.0)then

write(6,*)'sending topo to nodes'

      iun=52
!      open(iun,file='/home/usc/fm/gmm/mnt/store/GLOBALWTD/sciencepaper/sa_topo4.dat' &
      open(iun,file=filetopo&
!      ,form='unformatted',convert='little_endian',access='direct',recl=4*n2big*n3big)
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2big*n3big)

!      irec=1
     irec=9

     allocate(varreadbig(n2big,n3big))
     allocate(varread(n2,n3))

     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,23,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            topo(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(topo(is:ie,js:je),1,domblock(pid),0,23,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif


      landmask=1
      where(topo.lt.-1.e5)landmask=0
      where(topo.lt.-1.e5)topo=0.

!gmmread in fdepth

if(pid.eq.0)then

write(6,*)'sending fdepth to nodes'

!      iun=62
!      open(iun,file=filef&
!      ,form='unformatted',convert='little_endian',access='direct',recl=4*n2big*n3big)

!      irec=1

     allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

!     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

      statusnc = nf90_open(filef, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'F',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)!,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

     
     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

     where(varread.le.0.)varread=100.

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,33,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            fdepth(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(fdepth(is:ie,js:je),1,domblock(pid),0,33,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif

where(fdepth.lt.1.e-6)fdepth=100.

end subroutine readinitial

!     ******************************************************************
subroutine READLATLON(n2,n3,is,ie,js,je,lats,lons,area,dx,swlat,swlon)

integer :: n2,n3,is,ie,js,je,i,j,ii,jj
real, dimension(is:ie,js:je):: lats,lons,area
real :: swlat,swlon,xswlat,xswlon,dx,dy,rr,pii,deg2rad,xn,xs

!dx=1./120.

!xswlat=-55.9958333333333 + float(ns-1)*dx
!xswlon=-92.9958333333333 + float(nw-1)*dx !+ 360.
xswlat=swlat + float(ns-1)*dx
xswlon=swlon + float(nw-1)*dx !+ 360.
if(xswlon.gt.180.)xswlon=xswlon-360.


      do j=js,je
         do i=is,ie
             lats(i,j)=xswlat + float(j-1)*dx
         enddo
      enddo



      do j=js,je
         do i=is,ie
             lons(i,j)=xswlon + float(i-1)*dx
             if(lons(i,j).gt.180.)lons(i,j)=lons(i,j)-360.
         enddo
      enddo


rr=6370997.
pii=4.*atan(1.)
deg2rad=pii/180.
dy=rr*deg2rad*dx

     do j=js,je
        do i=is,ie
                  xn=(lats(i,j) + 0.5*dx)*deg2rad
                  xs=(lats(i,j) - 0.5*dx)*deg2rad
                  area(i,j)=(sin(xn)-sin(xs))*dy*rr
        enddo
      enddo



end subroutine readlatlon

!******************************************************************************************************
subroutine READWTD(n2,n3,is,ie,js,je,wtd)

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: wtd

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

if(pid.eq.0)then

write(6,*)'sending wtd to nodes'

      iun=2
      open(iun,file='/mnt/EMC/Store_uscfm/uscfmgmm/GLOBALWTD/sciencepaper/ncfileshaibin_newrun/sa_ewtd4.dat' &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2big*n3big)

      irec=1

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n2,1:n3)=min(-varreadbig(nw:ne,ns:nn),0.)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,53,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            wtd(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(wtd(is:ie,js:je),1,domblock(pid),0,53,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif


end subroutine readwtd

!******************************************************************************************************
subroutine READWTDNC(n2,n3,is,ie,js,je,wtd,filewtd)
use netcdf

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: wtd

character (len = *) :: filewtd

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
!integer :: start(3),count(3)

if(pid.eq.0)then

write(6,*)'sending wtd to nodes'

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))


!     start=(/1,1,1/)
!     count=(/n2big,n3big,1/)

      statusnc = nf90_open(filewtd, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'WTD',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)!,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

     varread(1:n2,1:n3)=min(-varreadbig(nw:ne,ns:nn),0.)

!write(6,*)'mirar veg',varread(100,100),varreadbig(nw+100,n3big*3+100)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,53,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            wtd(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(wtd(is:ie,js:je),1,domblock(pid),0,53,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif


end subroutine readwtdnc


!******************************************************************************************************
subroutine READVEG(n2,n3,is,ie,js,je,veg,fileveg)
use netcdf

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: veg
character (len = *) :: fileveg

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
!integer :: start(3),count(3)

if(pid.eq.0)then

write(6,*)'sending veg to nodes'

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))


!     start=(/1,1,1/)
!     count=(/n2big,n3big,1/)

      statusnc = nf90_open(fileveg, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'VEG',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)!,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!     varread(1:n2,1:n3)=varreadbig(nw:ne,(ns+n3big*3):(nn+n3big*3))
     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

!write(6,*)'mirar veg',varread(100,100),varreadbig(nw+100,n3big*3+100)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,83,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            veg(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(veg(is:ie,js:je),1,domblock(pid),0,83,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif


end subroutine readveg

!******************************************************************************************************
subroutine READHVEG(n2,n3,is,ie,js,je,hveg,filehveg)
use netcdf

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: hveg
character (len = *) :: filehveg

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
!integer :: start(3),count(3)

if(pid.eq.0)then

write(6,*)'sending hveg to nodes'

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))


!     start=(/1,1,1/)
!     count=(/n2big,n3big,1/)

      statusnc = nf90_open(filehveg, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'HVEG',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)!,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

     varreadbig=max(varreadbig,0.1)
!     varread(1:n2,1:n3)=varreadbig(nw:ne,(ns+n3big*3):(nn+n3big*3))
     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

!write(6,*)'mirar hveg',varread(100,100),varreadbig(nw+100-1,n3big*3+100)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,83,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            hveg(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(hveg(is:ie,js:je),1,domblock(pid),0,83,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif


end subroutine readhveg

!******************************************************************************************************
subroutine READSMOIEQ(n2,n3,nzg,is,ie,js,je,smoieq,filesmoieq)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(nzg,is:ie,js:je) :: smoieq
real, dimension(is:ie,js:je) :: var
character (len = *) :: filesmoieq

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)


DO k=1,nzg

if(pid.eq.0)then

write(6,*)'sending smoieq to nodes',k

!      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))



      statusnc = nf90_open(filesmoieq, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'EQSMOI',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      start=(/1,1,k/)
!      count=(/n2big,n3big,1/)
      count=(/n2,n3,1/)
!      statusnc = nf90_get_var(ncid,varid,varreadbig,start,count)
      statusnc = nf90_get_var(ncid,varid,varread,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!     varread(1:n2,1:n3)=varreadbig(nw:ne,(ns+n3big*3):(nn+n3big*3))
!     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,89,MPI_COMM_WORLD,req(n),ierr)
  enddo



   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            smoieq(k,i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif


deallocate(varread)!,varreadbig)



elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,89,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      smoieq(k,:,:)=var(:,:)

endif

ENDDO

end subroutine readsmoieq
!******************************************************************************************************
subroutine READFLOWDIRECTION(n2,n3,is,ie,js,je,fd,bfd,filerivers)

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
integer, dimension(is:ie,js:je) :: fd,bfd
real, dimension(is:ie,js:je) :: var

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)
character (len = *) :: filerivers


!flow direction
if(pid.eq.0)then

write(6,*)'sending fd to nodes'

      iun=2
      open(iun,file=filerivers &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2big*n3big)

      irec=1

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,53,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            fd(i,j)=nint(varread(i,j))
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(var(is,js),1,domblock(pid),0,53,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
      fd=nint(var)
endif


!now back flow direction
if(pid.eq.0)then

write(6,*)'sending bfd to nodes'

      iun=2
      open(iun,file=filerivers &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2big*n3big)

      irec=10

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,53,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            bfd(i,j)=nint(varread(i,j))
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(var(is,js),1,domblock(pid),0,53,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
      bfd=nint(var)
endif



end subroutine readflowdirection

!******************************************************************************************************
subroutine READRIVERPARAMETERS(n2,n3,is,ie,js,je,var,filerivers,irec)

integer :: n2,n3,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: var

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)
character (len = *) :: filerivers


!flow direction
if(pid.eq.0)then

write(6,*)'sending river par  to nodes',irec

      iun=2
      open(iun,file=filerivers &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2big*n3big)

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

     read(iun,rec=irec) ((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,53,MPI_COMM_WORLD,req(n),ierr)
  enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            var(i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

close(iun)

elseif(pid.lt.numtasks-1)then
      call MPI_irecv(var(is,js),1,domblock(pid),0,53,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)
endif

end subroutine readriverparameters
!     ******************************************************************

            subroutine handle_err(statusnc)
      use netcdf
      integer statusnc
      if (statusnc .ne. nf90_noerr) then
!        print*, nf_strerror(statusnc)
        stop 'Stopped'
      endif
      end subroutine handle_err

!     ******************************************************************

subroutine READHISTORY(n2,n3,is,ie,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,filename)
use netcdf

integer :: n2,n3,is,ie,js,je,nzg,k,irec,i,j,n
real, dimension(nzg,is:ie,js:je) :: smoi
real, dimension(is:ie,js:je) :: var,intercepstore,wtd,smoiwtd
real, allocatable :: varread(:,:)
character*80 :: filename

integer :: tag,request,req(nzg,numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2)


irec=1

if(pid.eq.0)then

      open(33,file=filename &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)

allocate(varread(n2,n3))


  DO k=1,nzg

      read(33,rec=irec) ((varread(i,j),i=1,n2),j=1,n3)
      irec=irec+1

   do n=1,numtasks-2
       tag=k
       call MPI_isend(varread(1,1),1,arraysection(n),n,tag,MPI_COMM_WORLD,req(k,n),ierr)
   enddo


   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            smoi(k,i,j)=varread(i,j)
         enddo
      enddo

   endif

  ENDDO

    call MPI_waitall(nzg*(numtasks-2),req,stats,ierr)

deallocate(varread)

close(33)

elseif(pid.lt.numtasks-1)then

  DO k=1,nzg

      tag=k
      call MPI_irecv(var(is,js),1,domblock(pid),0,tag,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      smoi(k,:,:)=var(:,:)

  ENDDO

endif



end subroutine readhistory

!******************************************************************************************************
subroutine READHISTORYNC(n2,n3,is,ie,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,inactivedays,filename)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(nzg,is:ie,js:je) :: smoi
integer, dimension(0:nzg+1,is:ie,js:je) :: inactivedays
real, dimension(is:ie,js:je) :: smoiwtd,wtd,intercepstore
real, dimension(is:ie,js:je) :: var
integer, dimension(is:ie,js:je) :: varint

character*200 filename

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)
integer*2, allocatable :: varreadint(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)


DO k=1,nzg

if(pid.eq.0)then

write(6,*)'sending smoi to nodes',k

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))


      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'SMOI',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      start=(/1,1,k/)
      count=(/n2big,n3big,1/)
      statusnc = nf90_get_var(ncid,varid,varreadbig,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,89,MPI_COMM_WORLD,req(n),ierr)
  enddo



   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            smoi(k,i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif


deallocate(varread,varreadbig)



elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,89,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      smoi(k,:,:)=var(:,:)

endif

ENDDO

!goto 1000

DO k=0,nzg+1

if(pid.eq.0)then

write(6,*)'sending inactivedasy to nodes',k

      allocate(varreadint(n2big,n3big))
      allocate(varread(n2,n3))


      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'INACTIVEDAYS',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      start=(/1,1,k+1/)
      count=(/n2big,n3big,1/)
      statusnc = nf90_get_var(ncid,varid,varreadint,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=float(varreadint(nw:ne,ns:nn))


  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,899,MPI_COMM_WORLD,req(n),ierr)
  enddo



   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            inactivedays(k,i,j)=nint(varread(i,j))
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif


deallocate(varreadint)
deallocate(varread)



elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,899,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      inactivedays(k,:,:)=nint(var(:,:))

endif

ENDDO

!1000 continue


if(pid.eq.0)then

write(6,*)'sending wtd to nodes'

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'WTD',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,91,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            wtd(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,91,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      wtd(:,:)=var(:,:)

endif

if(pid.eq.0)then

write(6,*)'sending intercepstore to nodes'

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'INTERCEPSTORE',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,92,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            intercepstore(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,92,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      intercepstore(:,:)=var(:,:)

endif


end subroutine readhistorync
!******************************************************************************************************
subroutine READHISTORYVARNC(n2,n3,is,ie,js,je,var,varname,filename)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: var

character (len = *) :: filename,varname

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)


if(pid.eq.0)then

write(6,*)'sending to nodes',varname(1:len_trim(varname))

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,varname(1:len_trim(varname)),varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,92,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            var(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,92,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

endif


end subroutine readhistoryvarnc

!******************************************************************************************************
subroutine READHISTORYBYTEVARNC(n2,n3,is,ie,js,je,var,varname,filename)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n
integer*1, dimension(is:ie,js:je) :: var

character (len = *) :: filename,varname

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

integer*1, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)


if(pid.eq.0)then

write(6,*)'sending to nodes',varname(1:len_trim(varname))

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,varname(1:len_trim(varname)),varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(nini_x(n):nend_x(n),nini_y(n):nend_y(n)),1,domblockbyte(n),n,97,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            var(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblockbyte(pid),0,97,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

endif


end subroutine readhistorybytevarnc
!******************************************************************************************************
subroutine READHISTORYVAR3DNC(n2,n3,is,ie,js,je,nzg,var3d,varname,filename)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n
real, dimension(is:ie,js:je) :: var
real, dimension(nzg,is:ie,js:je) :: var3d

character (len = *) :: filename,varname

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable :: varread(:,:),varreadbig(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)

DO k=1,nzg

if(pid.eq.0)then

write(6,*)'sending o18 to nodes',k

      allocate(varreadbig(n2big,n3big))
      allocate(varread(n2,n3))


      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,varname(1:len_trim(varname)),varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      start=(/1,1,k/)
      count=(/n2big,n3big,1/)
      statusnc = nf90_get_var(ncid,varid,varreadbig,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns:nn)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,189,MPI_COMM_WORLD,req(n),ierr)
  enddo



   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            var3d(k,i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif


deallocate(varread,varreadbig)



elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,189,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      var3d(k,:,:)=var(:,:)

endif

ENDDO

end subroutine readhistoryvar3dnc
!******************************************************************************************************
subroutine READHISTORYNCblock(n2,n3,is,ie,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,filename)
use netcdf

integer :: n2,n3,nzg,is,ie,js,je,i,j,irec,iun,k,n,n3big2,ns2,nn2
real, dimension(nzg,is:ie,js:je) :: smoi
real, dimension(is:ie,js:je) :: smoiwtd,wtd,intercepstore
real, dimension(is:ie,js:je) :: var

character*200 filename

integer :: tag,req(numtasks-2),stats(MPI_STATUS_SIZE,numtasks-2),request

real, allocatable ::  varreadbig(:,:),varread(:,:)

integer :: ncid,varid,statusnc
integer :: start(3),count(3)

!n2big2=
!ns2=
!nn2=


DO k=1,nzg

if(pid.eq.0)then

write(6,*)'sending smoi to nodes',k

      allocate(varreadbig(n2big,n3big2))
      allocate(varread(n2,n3))


      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'SMOI',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      start=(/1,1,k/)
      count=(/n2,n3,1/)
      statusnc = nf90_get_var(ncid,varid,varreadbig,start,count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns2:nn2)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,89,MPI_COMM_WORLD,req(n),ierr)
  enddo



   if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
            smoi(k,i,j)=varread(i,j)
         enddo
      enddo

   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif


deallocate(varread,varreadbig)



elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,89,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      smoi(k,:,:)=var(:,:)

endif

ENDDO

if(pid.eq.0)then

write(6,*)'sending smoiwtd to nodes'

      allocate(varreadbig(n2big,n3big2))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'SMOIWTD',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns2:nn2)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,90,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            smoiwtd(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,90,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      smoiwtd(:,:)=var(:,:)

endif

if(pid.eq.0)then

write(6,*)'sending wtd to nodes'

      allocate(varreadbig(n2big,n3big2))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'WTD',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns2:nn2)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,91,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            wtd(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,91,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      wtd(:,:)=var(:,:)

endif

if(pid.eq.0)then

write(6,*)'sending intercepstore to nodes'

      allocate(varreadbig(n2big,n3big2))
      allocate(varread(n2,n3))

      statusnc = nf90_open(filename, 0, ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_inq_varid(ncid,'INTERCEPSTORE',varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid,varid,varreadbig)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_close(ncid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      varread(1:n2,1:n3)=varreadbig(nw:ne,ns2:nn2)

  do n=1,numtasks-2
      call MPI_isend(varread(1,1),1,arraysection(n),n,92,MPI_COMM_WORLD,req(n),ierr)
  enddo

   if(numtasks.eq.1)then
      do j=js,je
         do i=is,ie
            intercepstore(i,j)=varread(i,j)
         enddo
      enddo
   else
         call MPI_waitall(numtasks-2,req,stats,ierr)
   endif

deallocate(varread,varreadbig)

elseif(pid.lt.numtasks-1)then

      call MPI_irecv(var(is,js),1,domblock(pid),0,92,MPI_COMM_WORLD,request,ierr)
      call MPI_wait(request,status,ierr)

      intercepstore(:,:)=var(:,:)

endif


end subroutine readhistoryncblock

!     ******************************************************************
subroutine WRITEOUTPUTFD(n2,n3,is,ie,js,je,nzg,smoi,waterdeficit,watext,qsrun,rech,pet,ppacum &
                      ,filename,irec,istart,req,req2,req3)

integer :: n2,n3,is,ie,js,je,nzg,k,irec,n,nss_x,nee_x,nss_y,nee_y,istart,i,j
real, dimension(nzg,is:ie,js:je) :: smoi,watext
real, dimension(is:ie,js:je) :: qsrun,rech,pet,ppacum,waterdeficit
real, dimension(is:ie,js:je,nzg) :: var,var2
real, dimension(is:ie,js:je,5) :: var3
real, allocatable :: varout(:,:)
character (len = *) :: filename
         
integer :: tag,request(numtasks-2),req(nzg),req2(nzg),req3(8),stats(MPI_STATUS_SIZE,numtasks-2)

             
IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

    open(33,file=filename & 
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)

  allocate(varout(n2,n3))
    
  DO k=1,nzg
    
    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    tag=k
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
               varout(i,j)=smoi(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

      irec=irec+1
      write(33,rec=irec) ((varout(i,j),i=1,n2),j=1,n3)
write(6,*)'written smoi, layer',k

  ENDDO

  deallocate(varout)

  close(33)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)call MPI_wait(req,status,ierr)

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

   do k=1,nzg

    var(:,:,k)=smoi(k,:,:)
    tag=k
    call MPI_isend(var(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req(k),ierr)

   enddo

ENDIF


!!!!!!!!now ncount

IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

    open(33,file=filename &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)


  allocate(varout(n2,n3))

  DO k=1,nzg

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1
    tag=k+nzg
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
!               varout(i,j)=float(ncount(k,i,j))
               varout(i,j)=watext(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

      irec=irec+1
      write(33,rec=irec) ((varout(i,j),i=1,n2),j=1,n3)

write(6,*)'written ncount, layer',k

  ENDDO

  deallocate(varout)

  close(33)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)call MPI_wait(req2,status,ierr)

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

   do k=1,nzg

!    var2(:,:,k)=float(ncount(k,:,:))
    var2(:,:,k)=watext(k,:,:)
    tag=k+nzg
    call MPI_isend(var2(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req2(k),ierr)

   enddo

ENDIF

!!!!!now wtd,smoiwtd,qsrun,qsprings,rech,qlat

IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

    open(33,file=filename &
      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)


  allocate(varout(n2,n3))

  DO k=1,5

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    tag=k+2*nzg
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

     select case(k)
        case(1)
           varout=qsrun
        case(2)
           varout=rech
        case(3)
           varout=pet
        case(4)
           varout=ppacum
        case(5)
           varout=waterdeficit
      end select


  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

      irec=irec+1
      write(33,rec=irec) ((varout(i,j),i=1,n2),j=1,n3)

write(6,*)'written other variable',k

  ENDDO

  deallocate(varout)

  close(33)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)call MPI_wait(req2,status,ierr)
   istart=0

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    var3(:,:,1)=qsrun(:,:)
    var3(:,:,2)=rech(:,:)
    var3(:,:,3)=pet(:,:)
    var3(:,:,4)=ppacum(:,:)
    var3(:,:,5)=waterdeficit(:,:)

   do k=1,5
    tag=k+2*nzg
    call MPI_isend(var3(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req3(k),ierr)

   enddo

ENDIF


end subroutine writeoutputfd

!     ******************************************************************
subroutine WRITEOUTPUTNC(n2,n3,is,ie,js,je,nzg,dz &
                         ,smoi,waterdeficit,watext,watextdeep,wtd,smoiwtd,qsrun,rech,qsprings &
                         ,qlat,et_s,et_i,et_c,ppacum,pppendepth,riverflowmean,qrf,delsfcwat,lai &
                         ,filename,irec,istart,var,var2,var3,req,req2,req3,date,nvar_out)

use netcdf
integer, parameter :: DateStrLen_len=19,NVAR_total=17
integer :: n2,n3,is,ie,js,je,nzg,nvar_out,k,irec,n,nss_x,nee_x,nss_y,nee_y,istart,i,j
real, dimension(nzg) :: dz
real, dimension(nzg,is:ie,js:je) :: smoi,watext
real, dimension(is:ie,js:je) :: wtd,smoiwtd,qsrun,rech,qsprings,qlat,et_s,et_i,et_c,ppacum &
                            ,waterdeficit,watextdeep,pppendepth,riverflowmean,qrf,delsfcwat,lai
real, dimension(is:ie,js:je,nzg) :: var,var2
real, dimension(is:ie,js:je,nvar_out) :: var3
real, allocatable :: varout(:,:)
integer*2, allocatable :: varoutint(:,:)
character (len = *) :: filename

integer :: tag,request(numtasks-2),req(nzg),req2(nzg),req3(nvar_out),stats(MPI_STATUS_SIZE,numtasks-2)

character*100 :: date
! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,DateStrLen_dim
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      double precision, dimension(NVAR_total) :: scale_factor, add_offset
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
     'soil moisture','SMOI','m3 m-3'& 
     ,'root activity','ROOT','mm'&
     ,'water table depth','WTD','m' &
     ,'surface runoff','QSRUN','mm' &
     ,'grouwater surface sip','QSPRING','mm'&
     ,'recharge','RECH','mm' &
     ,'lateral flow','QSLAT','mm' &
     ,'soil evaporation','ETS','mm' &
     ,'canopy transpiration','ETC','mm' &
     ,'evaporation from interception','ETI','mm' &
     ,'water deficit','WATDEF','mm'&
     ,'pp penetration depth','PPDEPTH','m' &
     ,'precipitation','PP','mm' &
     ,'LAI','LAI','' &
     ,'riverflow','RIVERFLOW','m3 s-1' &
     ,'baseflow','QRF','mm' &
     ,'infiltration from flooding','DELSFCWAT','mm'/

      real :: bound(2,NVAR_total)  !max and min of data
      data bound/&
          0.492,   0.0 &
         ,65.  ,   0. &
         ,0.   ,-1000. &
         ,650. ,   0. &
         ,65.  ,   0. &
         ,1000. ,-1000.    &
         ,100.  , -100. &
         ,65.   ,   0. & 
         ,65.   ,   0. &
         ,65.   ,   0. &
         ,65.   ,   0. &
         ,0.    , -1000. &
         ,650. ,   0. &
         ,100. ,  0.  &
        ,130000. ,  0. &
         ,650. ,   0. &
         ,6500. ,   0. /


if(pid.eq.numtasks-1.or.numtasks.eq.1)then

write(6,*)'now to open output file',filename

! enter define mode
      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
write(6,*)'now to define dimensions of filename',filename
! define dimensions
    iret = nf90_def_dim( ncid, 'time', 1, time_dim )
    call handle_err(iret)
    iret = nf90_def_dim(ncid, 'DateStrLen', DateStrLen_len, DateStrLen_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim', n3 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', n2 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)


write(6,*)'now defining variables'
!time
          type=NF90_CHAR
          VAR_rank=1
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = DateStrLen_dim
          iret = nf90_def_var(ncid, 'time', type,VAR_dims,time_id)
          call handle_err(iret)
          deallocate(VAR_dims)
!now the variables
      DO nvar=1,NVAR_total

       if(nvar.eq.1.or.nvar.eq.2)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim
          type = NF90_SHORT
        else
          VAR_rank=2
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
            select case (nvar)
                case(3,6,7,15,16)
                     type = NF90_REAL
                case default
                     type = NF90_SHORT
            end select
        endif

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar), 'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar), 'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)
       
          if(nvar.ne.3.and.nvar.ne.6.and.nvar.ne.7.and.nvar.ne.15.and.nvar.ne.16)then
          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**16-1.d+0)
          iret = nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)
          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**15. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)
          endif

       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)

write(6,*)'writing out date',date(1:DateStrLen_len)

!put time first
             iret = nf90_put_var(ncid, time_id, date(1:DateStrLen_len) )
             call handle_err(iret)

write(6,*)'date written'

endif

!now writeout soil moisture'
IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varout(n2,n3))
  allocate(varoutint(n2,n3))

  DO k=1,nzg

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    tag=k
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo
    


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
             varout(i,j) = smoi(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)


  endif

             nvar=1 !soil moisture
             varoutint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))

             count = (/ n2, n3,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint,start,count)
             call handle_err(iret)


write(6,*)'written smoi, layer',k

  ENDDO

  deallocate(varout,varoutint)


ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
                 do k=1,nzg
                   call MPI_wait(req(k),status,ierr)
                 enddo
   endif
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

   do k=1,nzg

    var(:,:,k)=smoi(k,:,:)
    tag=k
    call MPI_isend(var(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req(k),ierr)

   enddo


ENDIF

!!!!now root activity
IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varout(n2,n3))
  allocate(varoutint(n2,n3))

  DO k=1,nzg

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    tag=k+nzg
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
               varout(i,j)=watext(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

             nvar=2 !root activity

             varoutint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))

             count = (/ n2, n3,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint,start,count )
             call handle_err(iret)


write(6,*)'written root activity, layer',k

  ENDDO

  deallocate(varout,varoutint)


ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
             do k=1,nzg
               call MPI_wait(req2(k),status,ierr)
             enddo
   endif

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

   do k=1,nzg

    var2(:,:,k)=watext(k,:,:)
    tag=k+nzg
    call MPI_isend(var2(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req2(k),ierr)

   enddo

ENDIF

!!!!!now wtd,smoiwtd,qsrun,qsprings,rech,qlat

IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varout(n2,n3))
  allocate(varoutint(n2,n3))

  DO k=1,nvar_out

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1
    tag=k+2*nzg
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

     select case(k)
        case(1)
           varout=wtd
        case(2)
           varout=qsrun
        case(3)
           varout=qsprings
        case(4)
           varout=rech
        case(5)
           varout=qlat
        case(6)
           varout=et_s
        case(7)
           varout=et_c
        case(8)
           varout=et_i
        case(9)
           varout=waterdeficit
        case(10)
           varout=pppendepth
        case(11)
           varout=ppacum
        case(12)
           varout=lai
        case(13)
           varout=riverflowmean
        case(14)
           varout=qrf
        case(15)
           varout=delsfcwat
      end select


  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

             nvar=k+2 
           if(nvar.ne.3.and.nvar.ne.6.and.nvar.ne.7.and.nvar.ne.15.and.nvar.ne.16)then
             varoutint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint )
             call handle_err(iret)
           else
             iret = nf90_put_var(ncid, VAR_id(nvar), varout )
             call handle_err(iret)
           endif

write(6,*)'written other variable',k,VAR_tags(2,k+2)(1:len_trim(VAR_tags(2,k+2)))

  ENDDO

  deallocate(varout,varoutint)

     iret = nf90_close(ncid)
      call handle_err(iret)


ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
          do k=1,nvar_out
            call MPI_wait(req3(k),status,ierr)
          enddo

   endif
   istart=0

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    var3(:,:,1)=wtd(:,:)
    var3(:,:,2)=qsrun(:,:)
    var3(:,:,3)=qsprings(:,:)
    var3(:,:,4)=rech(:,:)
    var3(:,:,5)=qlat(:,:)
    var3(:,:,6)=et_s(:,:)
    var3(:,:,7)=et_c(:,:)
    var3(:,:,8)=et_i(:,:)
    var3(:,:,9)=waterdeficit(:,:)
    var3(:,:,10)=pppendepth(:,:)
    var3(:,:,11)=ppacum(:,:)
    var3(:,:,12)=lai(:,:)
    var3(:,:,13)=riverflowmean(:,:)
    var3(:,:,14)=qrf(:,:)
    var3(:,:,15)=delsfcwat(:,:)

   do k=1,nvar_out
    tag=k+2*nzg
    call MPI_isend(var3(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req3(k),ierr)
   enddo

ENDIF


end subroutine writeoutputnc

!     ******************************************************************
subroutine WRITEOUTPUTNC_par(n2,n3,is,ie,js,je,nzg,dz &
                         ,smoi,waterdeficit,watext,watextdeep,wtd,smoiwtd,qsrun,rech,qsprings &
                         ,qlat,et_s,et_i,et_c,ppacum,pppendepth,riverflowmean,qrf,delsfcwat &
                         ,filename,irec,istart,req,req2,req3,date,nvar_out)

use netcdf
integer, parameter :: DateStrLen_len=19,NVAR_total=16
integer :: n2,n3,is,ie,js,je,nzg,nvar_out,k,irec,n,nss_x,nee_x,nss_y,nee_y,istart,i,j
real, dimension(nzg) :: dz
real, dimension(nzg,is:ie,js:je) :: smoi,watext,var4
real, dimension(is:ie,js:je) :: wtd,smoiwtd,qsrun,rech,qsprings,qlat,et_s,et_i,et_c &
                             ,ppacum,waterdeficit,watextdeep,pppendepth,riverflowmean,qrf,delsfcwat
real, allocatable :: varout(:,:)
integer*2, allocatable :: varoutint(:,:)
integer*1, allocatable :: varoutint1(:,:)
integer*2 :: missvalue
character (len = *) :: filename

integer :: nmode,tag,request(numtasks-2),req(nzg),req2(nzg),req3(nvar_out),stats(MPI_STATUS_SIZE,numtasks-2)

character*100 :: date
! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,DateStrLen_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      double precision, dimension(NVAR_total) :: scale_factor, add_offset
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
     'root activity','ROOT','mm'&
     ,'soil moisture','SMOI','m3 m-3'&
     ,'water table depth','WTD','m' &
     ,'surface runoff','QSRUN','mm' &
     ,'grouwater surface sip','QSPRING','mm'&
     ,'recharge','RECH','mm' &
     ,'lateral flow','QSLAT','mm' &
     ,'soil evaporation','ETS','mm' &
     ,'canopy transpiration','ETC','mm' &
     ,'evaporation from interception','ETI','mm' &
     ,'water deficit','WATDEF','mm'&
     ,'pp penetration depth','PPDEPTH','m' &
     ,'precipitation','PP','mm' &
     ,'riverflow','RIVERFLOW','m3 s-1' &
     ,'baseflow','QRF','mm' &
     ,'infiltration from flooding','DELSFCWAT','mm'/

      real :: bound(2,NVAR_total)  !max and min of data
      data bound/&
          65.  ,   0. &
         ,0.492,   0. &
         ,0.   ,-1000. &
         ,650. ,   0. &
         ,650.  ,   0. &
         ,55. ,-600.    &
         ,327.  , -327. &
         ,65.   ,   0. &
         ,65.   ,   0. &
         ,65.   ,   0. &
         ,65.   ,   0. &
         ,0.    ,  -1000. &
         ,650. ,   0. &
         ,1300000., 0. &
         ,327. ,   -327. &
         ,6500.,      0. /


if(pid.gt.0.and.pid.lt.numtasks-1)then

    n=pid
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

          if(nini_x(n).eq.1)nss_x=1
          if(nini_y(n).eq.1)nss_y=1
          if(nend_x(n).eq.n2)nee_x=n2
          if(nend_y(n).eq.n3)nee_y=n3


    if(numtasks.eq.1)then
          nss_x=nini_x(n)
          nee_x=nend_x(n)
          nss_y=nini_y(n)
          nee_y=nend_y(n)
    endif

!write(6,*)'dims in output',pid,nss_x,nee_x,nss_y,nee_y,nini(pid),nend(pid)
!create file

      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
!write(6,*)'now to define dimensions of filename',filename,pid
! define dimensions
    iret = nf90_def_dim( ncid, 'time', 1, time_dim )
    call handle_err(iret)
    iret = nf90_def_dim(ncid, 'DateStrLen', DateStrLen_len, DateStrLen_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim', nee_y-nss_y+1 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', nee_x-nss_x+1 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_start', nss_y , y_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_end', nee_y , y_dim_end)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_start', nss_x , x_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_end', nee_x , x_dim_end)
    call handle_err(iret)


!time
          type=NF90_CHAR
          VAR_rank=1
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = DateStrLen_dim
          iret = nf90_def_var(ncid, 'time', type,VAR_dims,time_id)
          call handle_err(iret)
          deallocate(VAR_dims)
!now the variables
      DO nvar=1,NVAR_total

       if(nvar.eq.1)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim
          type = NF90_SHORT
        else
          VAR_rank=2
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
            select case (nvar)
                case(14)
                     type = NF90_REAL
                case default
                    type = NF90_SHORT
            end select
        endif

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar),'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar), 'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)


!          if(nvar.eq.1)then
!          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**8-1.d+0)
!          iret = nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
!          call handle_err(iret)
!          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**7. * scale_factor(nvar)
!          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
!          call handle_err(iret)

          if(nvar.ne.14)then

          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**16-1.d+0)
          iret =nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)

          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**15. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)

          endif

       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)

!put time first
             iret = nf90_put_var(ncid, time_id, date(1:DateStrLen_len) )
             call handle_err(iret)

write(6,*)'date written',pid


  allocate(varoutint(is:ie,js:je))


  DO k=1,nzg

             nvar=1 !root activity

!          do j=nss,nee
!           do i=is,ie
!           if(watext(k,i,j).ge.1.e-6)then
             varoutint(:,:)=nint((dble(watext(k,:,:)) - add_offset(nvar)) / scale_factor(nvar))
!             varoutint(i,j)=nint((dble(watext(k,i,j)) - add_offset(nvar)) / scale_factor(nvar))
!           else
!             varoutint(i,j)=0
!           endif
!           enddo
!          enddo

             count = (/ nee_x-nss_x+1, nee_y-nss_y+1,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y),start,count )
             call handle_err(iret)

  ENDDO

!  deallocate(varoutint1)
!write(6,*)'written root activity, layer',k

!  allocate(varoutint(is:ie,js:je))
  allocate(varout(is:ie,js:je))

   do k=1,nzg
       var4(k,:,:)=smoi(k,:,:)*dz(k)
   enddo

  DO k=1,nvar_out

     select case(k)
        case(1)
           varout=sum(var4,dim=1)/sum(dz)
        case(2)
           varout=wtd
        case(3)
           varout=qsrun
        case(4)
           varout=qsprings
        case(5)
           varout=rech
        case(6)
           varout=qlat
        case(7)
           varout=et_s
        case(8)
           varout=et_c
        case(9)
           varout=et_i
        case(10)
           varout=waterdeficit
        case(11)
           varout=pppendepth
        case(12)
           varout=ppacum
        case(13)
           varout=riverflowmean
        case(14)
           varout=qrf
        case(15)
           varout=delsfcwat
      end select


             nvar=k+1
           if(nvar.ne.14)then
!            where(varout.gt.1.e-6)
             varoutint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))
!            elsewhere
!             varoutint=0
!            endwhere
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)
           else
             iret = nf90_put_var(ncid, VAR_id(nvar), varout(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)
           endif

   ENDDO

  deallocate(varout,varoutint)

     iret = nf90_close(ncid)
      call handle_err(iret)


endif

end subroutine writeoutputnc_par
!     ******************************************************************
subroutine WRITEHISTORYNC_par(n2,n3,is,ie,js,je,nzg,smoi,intercepstore,wtd,inactivedays &
                              ,riverflow,riverdepth,floodheight &
                              ,o18 &
                              ,filename,date)

use netcdf
integer, parameter :: DateStrLen_len=19,NVAR_total=8
integer :: n2,n3,is,ie,js,je,nzg,k,n,nss_x,nee_x,nss_y,nee_y,i,j
real, dimension(nzg,is:ie,js:je) :: smoi,o18
integer, dimension(0:nzg+1,is:ie,js:je) :: inactivedays
real, dimension(is:ie,js:je) :: wtd,intercepstore,riverflow,riverdepth,floodheight
real, allocatable :: varout(:,:),varout3d(:,:,:)
integer*2, allocatable :: varoutint(:,:,:)
character (len = *) :: filename
integer :: stats(MPI_STATUS_SIZE,numtasks-2)

character*100 :: date
! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,zext_dim,DateStrLen_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
     'soil moisture','SMOI','m3 m-3'&
     ,'inactive days for roots','INACTIVEDAYS','days' &
     ,'water table depth','WTD','m' &
     ,'intercepted water','INTERCEPSTORE','mm' &
     ,'riverflow','RIVERFLOW','m3 s-1' &
     ,'riverdepth','RIVERDEPTH','m' &
     ,'surface water depth','FLOODHEIGHT','m' &
     ,'O18','O18','m3 m-3' /

if(pid.gt.0.and.pid.lt.numtasks-1)then

    n=pid

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

          if(nini_x(n).eq.1)nss_x=1
          if(nini_y(n).eq.1)nss_y=1
          if(nend_x(n).eq.n2)nee_x=n2
          if(nend_y(n).eq.n3)nee_y=n3


    if(numtasks.eq.1)then
          nss_x=nini_x(n)
          nee_x=nend_x(n)
          nss_y=nini_y(n)
          nee_y=nend_y(n)
    endif

!write(6,*)'dims in history',pid,nss_x,nee_x,nss_y,nee_y,nini(pid),nend(pid)
!create file

      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
if(pid.eq.1)write(6,*)'now to define dimensions of filename',filename,pid
! define dimensions
    iret = nf90_def_dim( ncid, 'time', 1, time_dim )
    call handle_err(iret)
    iret = nf90_def_dim(ncid, 'DateStrLen', DateStrLen_len, DateStrLen_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim', nee_y-nss_y+1 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', nee_x-nss_x+1 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zextdim', nzg+2 , zext_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_start', nss_y , y_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_end', nee_y , y_dim_end)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_start', nss_x , x_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_end', nee_x , x_dim_end)
    call handle_err(iret)


if(pid.eq.1)write(6,*)'finished with dimensions',filename,pid

!time
          type=NF90_CHAR
          VAR_rank=1
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = DateStrLen_dim
          iret = nf90_def_var(ncid, 'time', type,VAR_dims,time_id)
          call handle_err(iret)
          deallocate(VAR_dims)
!now the variables
      DO nvar=1,NVAR_total

       if(nvar.eq.1.or.nvar.eq.8)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim
          type = NF90_REAL
       elseif(nvar.eq.2)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = zext_dim
          type = NF90_SHORT
        else
          VAR_rank=2
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          type = NF90_REAL
        endif

!if(pid.eq.1)write(6,*)'defining variable',filename,pid,nvar

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid,VAR_id(nvar),'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar),'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)

!if(pid.eq.1)write(6,*)'finished with variable',filename,pid,nvar
       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)

!if(pid.eq.1)write(6,*)'finished with all variables',filename,pid

!put time first
             iret = nf90_put_var(ncid, time_id, date(1:DateStrLen_len) )
             call handle_err(iret)

write(6,*)'date written history',pid

  allocate(varout3d(nss_x:nee_x,nss_y:nee_y,nzg))

  do k=1,nzg
     do j=nss_y,nee_y
       do i=nss_x,nee_x
             varout3d(i,j,k)=smoi(k,i,j)
          enddo
       enddo
  enddo

             nvar=1 !soil moisture
             iret = nf90_put_var(ncid, VAR_id(nvar),varout3d)
             call handle_err(iret)

  do k=1,nzg
     do j=nss_y,nee_y
       do i=nss_x,nee_x
             varout3d(i,j,k)=o18(k,i,j)
          enddo
       enddo
  enddo

             nvar=8 !o18 ratio
             iret = nf90_put_var(ncid, VAR_id(nvar),varout3d)
             call handle_err(iret)

  deallocate(varout3d)


  allocate(varoutint(nss_x:nee_x,nss_y:nee_y,0:nzg+1))

  do k=0,nzg+1
     do j=nss_y,nee_y
       do i=nss_x,nee_x
             varoutint(i,j,k)=inactivedays(k,i,j)
          enddo
       enddo
  enddo

             nvar=2 !inactive days

             iret = nf90_put_var(ncid,VAR_id(nvar),varoutint)
             call handle_err(iret)

  deallocate(varoutint)

  allocate(varout(is:ie,js:je))

  DO k=3,7

     select case(k)
        case(3)
           varout=wtd
        case(4)
           varout=intercepstore
        case(5)
           varout=riverflow
        case(6)
           varout=riverdepth
        case(7)
           varout=floodheight
     end select

        nvar=k

            iret = nf90_put_var(ncid, VAR_id(nvar), varout(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)

   ENDDO

  deallocate(varout)

     iret = nf90_close(ncid)
      call handle_err(iret)


endif

end subroutine writehistorync_par

!     ******************************************************************

subroutine WRITEHISTORYNC(n2,n3,is,ie,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,inactivedays &
                         ,riverflow,riverdepth,floodheight,filename,date,istart,var,var2,var3,req,req2,req3)

use netcdf
integer, parameter :: DateStrLen_len=19,NVAR_total=7
integer :: n2,n3,is,ie,js,je,nzg,k,n,nss_x,nee_x,nss_y,nee_y,i,j,istart
real, dimension(nzg,is:ie,js:je) :: smoi
integer, dimension(0:nzg+1,is:ie,js:je) :: inactivedays
real, dimension(is:ie,js:je) :: wtd,smoiwtd,intercepstore,riverflow,riverdepth,floodheight
real, dimension(is:ie,js:je,nzg) :: var
real, dimension(is:ie,js:je,0:nzg+1) :: var2
real, dimension(is:ie,js:je,3:NVAR_total) :: var3
real, allocatable :: varout(:,:)
integer*2, allocatable :: varoutint(:,:)
character (len = *) :: filename
integer :: stats(MPI_STATUS_SIZE,numtasks-2),tag,request(numtasks-2),req(nzg),req2(0:nzg+1),req3(3:NVAR_total)

character*100 :: date
! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,zext_dim,DateStrLen_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
     'soil moisture','SMOI','m3 m-3'&
     ,'inactive days for roots','INACTIVEDAYS','days' &
     ,'water table depth','WTD','m' &
     ,'intercepted water','INTERCEPSTORE','mm' &
     ,'riverflow','RIVERFLOW','m3 s-1' &
     ,'riverdepth','RIVERDEPTH','m' &
     ,'surface water depth','FLOODHEIGHT','m'/


if(pid.eq.numtasks-1.or.numtasks.eq.1)then

!write(6,*)'dims in history',pid,nss_x,nee_x,nss_y,nee_y,nini(pid),nend(pid)
!create file

      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
!write(6,*)'now to define dimensions of filename',filename,pid
! define dimensions
    iret = nf90_def_dim( ncid, 'time', 1, time_dim )
    call handle_err(iret)
    iret = nf90_def_dim(ncid, 'DateStrLen', DateStrLen_len, DateStrLen_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim', n3 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', n2 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zextdim', nzg+2 , zext_dim)
    call handle_err(iret)


!time
          type=NF90_CHAR
          VAR_rank=1
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = DateStrLen_dim
          iret = nf90_def_var(ncid, 'time', type,VAR_dims,time_id)
          call handle_err(iret)
          deallocate(VAR_dims)
!now the variables
      DO nvar=1,NVAR_total

       if(nvar.eq.1)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim
          type = NF90_REAL
       elseif(nvar.eq.2)then
          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = zext_dim
          type = NF90_SHORT
        else
          VAR_rank=2
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          type = NF90_REAL
        endif

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid,VAR_id(nvar),'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar),'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)

       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)

!put time first
             iret = nf90_put_var(ncid, time_id, date(1:DateStrLen_len) )
             call handle_err(iret)

write(6,*)'date written'

endif

IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varout(n2,n3))

  DO k=1,nzg

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1    
    tag=k
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
             varout(i,j) = smoi(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif



             nvar=1 !soil moisture

             count = (/ n2, n3,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid,VAR_id(nvar),varout,start,count)
             call handle_err(iret)


!write(6,*)'written smoi, layer',k

  ENDDO

  deallocate(varout)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
                 do k=1,nzg
                   call MPI_wait(req(k),status,ierr)
                 enddo
   endif
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1   

   do k=1,nzg

    var(:,:,k)=smoi(k,:,:)
    tag=k
    call MPI_isend(var(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req(k),ierr)

   enddo

ENDIF

!now inactivedays

IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varoutint(n2,n3))
  allocate(varout(n2,n3))

  DO k=0,nzg+1

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1    
    tag=k
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo


  if(numtasks.eq.1)then

      do j=js,je
         do i=is,ie
             varoutint(i,j) = inactivedays(k,i,j)
         enddo
      enddo

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

      varoutint=nint(varout)

  endif



             nvar=2 !inactivedays

             count = (/ n2, n3,  1 /)
             start = (/ 1, 1, k+1 /)
             iret = nf90_put_var(ncid,VAR_id(nvar),varoutint,start,count)
             call handle_err(iret)


!write(6,*)'written smoi, layer',k

  ENDDO

  deallocate(varoutint)
  deallocate(varout)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
                 do k=0,nzg+1
                   call MPI_wait(req2(k),status,ierr)
                 enddo
   endif
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1   

   do k=0,nzg+1

    var2(:,:,k)=float(inactivedays(k,:,:))
    tag=k
    call MPI_isend(var2(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req2(k),ierr)

   enddo

ENDIF



IF(pid.eq.numtasks-1.or.numtasks.eq.1)then

  allocate(varout(n2,n3))

  DO k=3,NVAR_total

    do n=1,numtasks-2
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1    
    tag=k+nzg
    call MPI_irecv(varout(nss_x:nee_x,nss_y:nee_y),1,domblocksmall(n),n,tag,MPI_COMM_WORLD,request(n),ierr)
    enddo

  if(numtasks.eq.1)then

     select case(k)
        case(3)
           varout=wtd
        case(4)
           varout=intercepstore
        case(5)
           varout=riverflow
        case(6)
           varout=riverdepth
        case(7)
           varout=floodheight
     end select

  else
      call MPI_waitall(numtasks-2,request,stats,ierr)

  endif

        nvar=k

            iret = nf90_put_var(ncid, VAR_id(nvar), varout )
             call handle_err(iret)

   ENDDO

  deallocate(varout)

     iret = nf90_close(ncid)
      call handle_err(iret)

ELSEIF(pid.gt.0)then

!wait for previous output to be completed
   if(istart.eq.0.and.numtasks.gt.1)then
          do k=3,NVAR_total
            call MPI_wait(req3(k),status,ierr)
          enddo

   endif
   istart=0

          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

    var3(:,:,3)=wtd(:,:)
    var3(:,:,4)=intercepstore(:,:)
    var3(:,:,5)=riverflow(:,:)
    var3(:,:,6)=riverdepth(:,:)
    var3(:,:,7)=floodheight(:,:)


   do k=3,NVAR_total
    tag=k+nzg
    call MPI_isend(var3(nss_x:nee_x,nss_y:nee_y,k),1,domblocksmall(pid),numtasks-1,tag,MPI_COMM_WORLD,req3(k),ierr)
   enddo

ENDIF


end subroutine writehistorync
!     ******************************************************************
subroutine WRITEEQSMOINC(n2,n3,is,ie,js,je,nzg,smoieq)

use netcdf
integer :: n2,n3,is,ie,js,je,nzg,k,n,nss_x,nee_x,nss_y,nee_y,i,j,istart
real, dimension(nzg,is:ie,js:je) :: smoieq
real, allocatable :: varout(:,:)
integer :: ncid,x_dim,y_dim,z_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end,iret,VAR_id,type
integer :: VAR_dims(3),start(3),count(3)
character*200 :: filename
integer :: stats(MPI_STATUS_SIZE,numtasks-2),tag,request(numtasks-2),req(nzg)


if(pid.gt.0.and.pid.lt.numtasks-1)then

    n=pid
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1    

    if(numtasks.eq.1)then
          nss_x=nini_x(n)
          nee_x=nend_x(n)
          nss_y=nini_y(n)
          nee_y=nend_y(n)
    endif

write(filename,'(a14,i3.3,a3)')'smoieq_SA_new_',pid,'.nc'
filename='/mnt/EMC/Store_uscfm/uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/outputpareqsmoi/'//filename

!create file
      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
write(6,*)'now to define dimensions of filename',filename
! define dimensions

    iret =  nf90_def_dim( ncid, 'ydim', nee_y-nss_y+1 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', nee_x-nss_x+1 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_start', nss_y , y_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_end', nee_y , y_dim_end)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_start', nss_x , x_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_end', nee_x , x_dim_end)    
    call handle_err(iret)


          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim
          type = NF90_REAL

          iret = nf90_def_var(ncid, 'EQSMOI' &
                     , type, VAR_dims, VAR_id,deflate_level = 4)
          call handle_err(iret)

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)



  allocate(varout(is:ie,js:je))

   DO k=1,nzg

    varout(:,:)=smoieq(k,:,:)


             count = (/ nee_x-nss_x+1, nee_y-nss_y+1,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid,VAR_id,varout(nss_x:nee_x,nss_y:nee_y),start,count)
             call handle_err(iret)


!write(6,*)'written smoi, layer',k

  ENDDO

  deallocate(varout)

     iret = nf90_close(ncid)
      call handle_err(iret)

endif

end subroutine writeeqsmoinc

!     ******************************************************************
subroutine WRITEOUTPUTNC_DAILY_par(n2,n3,is,ie,js,je,nzg,ntim &
                        ,wtdflux,et_s_daily,et_c_daily,transptop,infilk,topsmoi,wtd_daily &
                        ,daysforoutput,filename ) 
use netcdf

integer, parameter :: NVAR_total=7
integer :: n2,n3,is,ie,js,je,nzg,k,n,nss_x,nee_x,nss_y,nee_y,istart,i,j,ntim
real, dimension(is:ie,js:je,8) :: wtdflux,et_s_daily,et_c_daily,transptop,transpfrtop,topsmoi,varout,wtd_daily
integer*1, dimension(is:ie,js:je,8) :: infilk,varoutint
integer*2, dimension(is:ie,js:je,8) :: varint
integer :: daysforoutput(8)
character (len = *) :: filename

! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,DateStrLen_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      double precision, dimension(NVAR_total) :: scale_factor, add_offset
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
      '0-0.1m soil moisture','SMOI','mm' &
     ,'soil evaporation','ETS','mm/day' &
     ,'total transpiration','ETC','mm/day' &
     ,'root uptake from 0-0.1m','TOPTRANS','mm/day' &
     ,'water table depth','WTD','m' &
     ,'groundwater recharge','WTDFLUX','mm/day' &
     ,'layer number of infiltration depth','INFILK','' /

      real :: bound(2,NVAR_total)  !max and min of data
      data bound/&
          50. , 0.    &
         ,50.   , 0.    &
         ,50.   ,   0. &
         ,50.   ,   0. &
         ,0      ,-1000. &
         ,24.5   ,  0. &
         ,254.   ,   0. /


if((pid.gt.0.and.pid.lt.numtasks-1).or.(pid.eq.0.and.numtasks.eq.1))then

    n=pid
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

          if(nini_x(n).eq.1)nss_x=1
          if(nini_y(n).eq.1)nss_y=1
          if(nend_x(n).eq.n2)nee_x=n2
          if(nend_y(n).eq.n3)nee_y=n3

    if(numtasks.eq.1)then
          nss_x=nini_x(n)
          nee_x=nend_x(n)
          nss_y=nini_y(n)
          nee_y=nend_y(n)
    endif

!write(6,*)'dims in output',pid,nss_x,nee_x,nss_y,nee_y,nini(pid),nend(pid)
!create file

      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
!write(6,*)'now to define dimensions of filename',filename,pid
! define dimensions
    iret = nf90_def_dim( ncid, 'time', ntim, time_dim )
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim', nee_y-nss_y+1 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', nee_x-nss_x+1 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_start', nss_y , y_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_end', nee_y , y_dim_end)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_start', nss_x , x_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_end', nee_x , x_dim_end)
    call handle_err(iret)


!time
          type=NF90_INT
          VAR_rank=1
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = time_dim
          iret = nf90_def_var(ncid, 'time', type,VAR_dims,time_id)
          call handle_err(iret)
          deallocate(VAR_dims)

          iret = nf90_put_att(ncid, time_id, 'units','days since 01-01-2003' )
          call handle_err(iret)

!now the variables
      DO nvar=1,NVAR_total

          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = time_dim

          if(nvar.ge.6)then
                    type = NF90_BYTE
          else
                    type = NF90_SHORT
          endif

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar),'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar), 'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)


          if(nvar.lt.6)then

          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**16-1.d+0)
          iret =nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)

          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**15. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)

          elseif(nvar.eq.6)then

          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**8-1.d+0)
          iret =nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)

          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**7. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)

          endif

        

       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)


!write time variable

               iret = nf90_put_var(ncid, time_id, daysforoutput(1:ntim))
               call handle_err(iret)

          do nvar=1,NVAR_total

             select case(nvar)
                 case(1)
                   varout = topsmoi
                 case(2)
                   varout = min(et_s_daily/5.,bound(1,nvar)) 
                 case(3)
                   varout = min(et_c_daily/5.,bound(1,nvar))
                 case(4)
                   varout = min(transptop/5.,bound(1,nvar))
                 case(5)
                   varout = wtd_daily
                 case(6)
                   varout = max(min(wtdflux/5.,bound(1,nvar)),bound(2,nvar))   
             end select

             if(nvar.le.5)then
               varint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))

               iret = nf90_put_var(ncid, VAR_id(nvar), varint(nss_x:nee_x,nss_y:nee_y,1:ntim))
               call handle_err(iret)

             elseif(nvar.eq.6)then
               varoutint=nint((dble(varout) - add_offset(nvar)) / scale_factor(nvar))

               iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y,1:ntim))
               call handle_err(iret)

             else
               iret = nf90_put_var(ncid, VAR_id(nvar), infilk(nss_x:nee_x,nss_y:nee_y,1:ntim))
               call handle_err(iret)
             endif

            
           enddo

     iret = nf90_close(ncid)
      call handle_err(iret)

endif

end subroutine writeoutputnc_daily_par
!     ******************************************************************
subroutine WRITEOUTPUTNC_INFIL_par(n2,n3,is,ie,js,je,nzg,monthhours &
                        ,infilflux,o18,smoi,transpo18ratio &
                        ,upflux,wtdmax,wtdmin,qlatin,qlatout,qlatino18ratio,qlatouto18ratio &
                        ,filename )
use netcdf

integer, parameter :: NVAR_total=10
integer :: n2,n3,is,ie,js,je,nzg,k,n,nss_x,nee_x,nss_y,nee_y,istart,i,j,monthhours
real, dimension(nzg,is:ie,js:je) :: infilflux,o18,smoi,upflux
real, dimension(is:ie,js:je) :: transpo18ratio,wtdmax,wtdmin,qlatin,qlatout,qlatino18ratio,qlatouto18ratio
!integer*2, dimension(nzg,is:ie,js:je) :: infilcounter
integer*2, dimension(is:ie,js:je) :: varoutint
integer*1, dimension(is:ie,js:je) :: varoutbyte
real, dimension(is:ie,js:je) :: var

character (len = *) :: filename

! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids

      integer time_dim
      integer x_dim,y_dim,z_dim,DateStrLen_dim,x_dim_start,x_dim_end,y_dim_start,y_dim_end
      integer lon_id,lat_id,time_id,varid
      integer :: len,ln,nvar,itime
      character*80 outfile

      integer type
      double precision, dimension(NVAR_total) :: scale_factor, add_offset
      integer :: start(3), count(3)
! variable id
      integer  Times_id
      integer VAR_id(100)
! rank (number of dimensions) for each variable
      integer ::VAR_rank
! variable shapes
      integer, allocatable :: VAR_dims(:)
! data variables
!descriptions
      character*500 :: VAR_tags(3,NVAR_total),text,descript
      data VAR_tags/ &
      'downward flux at the top of the layer','INFILFLUX','mm/day' &
      ,'wtd max','WTDMAX','m' &
      ,'wtd min','WTDMIN','m' &
      ,'lateral flow in','QLATIN','mm/day' &
      ,'lateral flow out','QLATOUT','mm/day' &
      ,'O18 isotopic ratio','O18','' &
      ,'upward flux at the top of the layer','UPFLUX','mm/day' &
      ,'O18 ratio in transpiration','TRANSPO18','' &
      ,'O18 ratio in lateral flow in','QLATINO18','' &
      ,'O18 ratio in lateral flow out','QLATOUTO18','' /

!     ,'fraction of time with infiltration','INFILFR','' /
!     ,'number of days with infiltration','INFILFR','days' /

      real :: bound(2,NVAR_total)  !max and min of data
      data bound/&
          0.  ,  -65. &
         ,0.   , -1000. &
         ,0.   , -1000. &
         ,327. , 0.  &
         ,327. , 0.  &
         ,50.   ,  -50. &
         ,50.,   0. & 
         ,50.   ,  -50. &
         ,50.   ,  -50. &
         ,50.   ,  -50. /


if((pid.gt.0.and.pid.lt.numtasks-1).or.(pid.eq.0.and.numtasks.eq.1))then

    n=pid
          nss_x=nini_x(n)+1
          nee_x=nend_x(n)-1
          nss_y=nini_y(n)+1
          nee_y=nend_y(n)-1

          if(nini_x(n).eq.1)nss_x=1
          if(nini_y(n).eq.1)nss_y=1
          if(nend_x(n).eq.n2)nee_x=n2
          if(nend_y(n).eq.n3)nee_y=n3

    if(numtasks.eq.1)then
          nss_x=nini_x(n)
          nee_x=nend_x(n)
          nss_y=nini_y(n)
          nee_y=nend_y(n)
    endif


!write(6,*)'dims in output',pid,nss_x,nee_x,nss_y,nee_y,nini(pid),nend(pid)
!create file

      iret = nf90_create(filename, NF90_NETCDF4, ncid)
      call handle_err(iret)
!write(6,*)'now to define dimensions of filename',filename,pid
! define dimensions
    iret =  nf90_def_dim( ncid, 'ydim', nee_y-nss_y+1 , y_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim', nee_x-nss_x+1 , x_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'zdim', nzg , z_dim)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_start', nss_y , y_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'ydim_end', nee_y , y_dim_end)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_start', nss_x , x_dim_start)
    call handle_err(iret)
    iret =  nf90_def_dim( ncid, 'xdim_end', nee_x , x_dim_end)
    call handle_err(iret)


!now the variables
      DO nvar=1,NVAR_total

          select case(nvar)
              case(1,6,7)

          VAR_rank=3
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim
          VAR_dims(3) = z_dim

              case(2,3,4,5,8,9,10)
          VAR_rank=2
          allocate(VAR_dims(VAR_rank))
          VAR_dims(1) = x_dim
          VAR_dims(2) = y_dim

          end select

          if(nvar.le.6)then
                    type = NF90_SHORT
          else
                    type = NF90_BYTE
          endif

          ln=len_trim(VAR_tags(2,nvar))
          iret = nf90_def_var(ncid, VAR_tags(2,nvar)(1:ln) &
                     , type, VAR_dims, VAR_id(nvar),deflate_level = 4)
          call handle_err(iret)
          deallocate(VAR_dims)

          ln=len_trim(VAR_tags(1,nvar))
          iret = nf90_put_att(ncid,VAR_id(nvar),'description',VAR_tags(1,nvar)(1:ln))
          call handle_err(iret)

          ln=len_trim(VAR_tags(3,nvar))
          iret = nf90_put_att(ncid, VAR_id(nvar), 'units',VAR_tags(3,nvar)(1:ln))
          call handle_err(iret)


          if(nvar.le.6)then

          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**16-1.d+0)
          iret =nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)

          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**15. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)

          else

          scale_factor(nvar) = dble(bound(1,nvar)-bound(2,nvar))/(2.d+0**8-1.d+0)
          iret =nf90_put_att(ncid,VAR_id(nvar),'scale_factor',scale_factor(nvar))
          call handle_err(iret)

          add_offset(nvar) = dble(bound(2,nvar)) + 2.d+0**7. * scale_factor(nvar)
          iret = nf90_put_att(ncid, VAR_id(nvar),'add_offset',add_offset(nvar))
          call handle_err(iret)

          iret = nf90_put_att(ncid,VAR_id(nvar),'missing_value',-128)

          endif



       ENDDO

! leave define mode
      iret = nf90_enddef(ncid)
      call handle_err(iret)

nvar=1

            infilflux = max(min(infilflux,bound(1,nvar)),bound(2,nvar))     

  DO k=1,nzg

             varoutint(:,:)=nint((dble(infilflux(k,:,:)) - add_offset(nvar)) / scale_factor(nvar))

             count = (/ nee_x-nss_x+1, nee_y-nss_y+1,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y),start,count )
             call handle_err(iret)

  ENDDO

do nvar=2,5
          select case(nvar)
             case(2)
                var=wtdmax
             case(3)
                var=wtdmin
             case(4)
                var=qlatin
             case(5)
                var=qlatout
          end select

             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutint(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))

             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)

enddo

nvar=6

  DO k=1,nzg
             var(:,:) = ( 1.e6 * (o18(k,:,:)/smoi(k,:,:)) / 2005.2 - 1. ) * 1.e3
             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutint(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))

             count = (/ nee_x-nss_x+1, nee_y-nss_y+1,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutint(nss_x:nee_x,nss_y:nee_y),start,count )
             call handle_err(iret)

  ENDDO

do nvar=7,7

  DO k=1,nzg
!          select case(nvar)
!             case(6)
!             var(:,:) = ( 1.e6 * (o18(k,:,:)/smoi(k,:,:)) / 2005.2 - 1. ) * 1.e3
!             case(7)
             var(:,:) = upflux(k,:,:)
!           end select
             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutbyte(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))
!             varoutbyte(:,:)=infilcounter(k,:,:)

             count = (/ nee_x-nss_x+1, nee_y-nss_y+1,  1 /)
             start = (/ 1, 1, k /)
             iret = nf90_put_var(ncid, VAR_id(nvar), varoutbyte(nss_x:nee_x,nss_y:nee_y),start,count )
             call handle_err(iret)

  ENDDO

enddo

nvar=8

             var(:,:) = ( 1.e6 * transpo18ratio(:,:) / 2005.2 - 1. ) * 1.e3
             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutbyte(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))

             iret = nf90_put_var(ncid, VAR_id(nvar), varoutbyte(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)

nvar=9

             var(:,:) = ( 1.e6 * qlatino18ratio(:,:) / 2005.2 - 1. ) * 1.e3
             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutbyte(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))

             iret = nf90_put_var(ncid, VAR_id(nvar), varoutbyte(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)

nvar=10

             var(:,:) = ( 1.e6 * qlatouto18ratio(:,:) / 2005.2 - 1. ) * 1.e3
             var = max(min(var,bound(1,nvar)),bound(2,nvar))
             varoutbyte(:,:)=nint((dble(var(:,:)) - add_offset(nvar)) / scale_factor(nvar))

             iret = nf90_put_var(ncid, VAR_id(nvar), varoutbyte(nss_x:nee_x,nss_y:nee_y) )
             call handle_err(iret)


     iret = nf90_close(ncid)
      call handle_err(iret)

endif

end subroutine writeoutputnc_infil_par

!     ******************************************************************
END MODULE MODULE_IO
