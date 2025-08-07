MODULE module_parallel
   implicit none
   include 'mpif.h'

!maximumn number of processors
   integer, parameter :: npmax = 1000

   integer, save, dimension(0:npmax) :: domblock, domblock2dint, domblock3d, domblock3dint &
                                        , domblocksmall, domblocksmallint, domblock3dsmall &
                                        , nini_y, nend_y, nini_x, nend_x, n_right, n_left, n_up, n_down &
                                        , rcountblock, rcountblocksmall, disp, domblockbyte &
                                        , leftrighthalo, updownhalo, arraysection, arraysectionint &
                                        , borderupdown, borderleftright
   integer, save :: columntype, columntype2
   integer, save :: pid, numtasks
   integer :: status(MPI_STATUS_SIZE), ierr

!integer, parameter :: n2big=7320,n3big=8520,nw=2501,ne=2600,ns=401,nn=500
!integer, parameter :: n2big=7320,n3big=8520,nw=2901,ne=3000,ns=6301,nn=6400
!integer, parameter :: n2big=7320,n3big=8520,nw=4151,ne=4250,ns=2701,nn=2800
!integer, parameter :: n2big=7320,n3big=8520,nw=3031,ne=3130,ns=6681,nn=6780
!integer, parameter :: n2big=7320,n3big=8520,nw=4401,ne=4500,ns=4481,nn=4580
!integer, parameter :: n2big=7320,n3big=8520,nw=4401,ne=4500,ns=5151,nn=5250
!manaosinteger, parameter :: n2big=7320,n3big=8520,nw=4001,ne=4500,ns=6391,nn=6890
!pantanalinteger, parameter :: n2big=7320,n3big=8520,nw=4201,ne=4500,ns=4501,nn=4800
!testinteger, parameter :: n2big=7320,n3big=8520,nw=4601,ne=4700,ns=6301,nn=6400
!integer, parameter :: n2big=7320,n3big=8520,nw=2851,ne=3250,ns=3251,nn=3650
!argentinainteger, parameter :: n2big=7320,n3big=8520,nw=2651,ne=2750,ns=2801,nn=2900
!integer, parameter :: n2big=7320,n3big=8520,nw=2950,ne=3051,ns=2801,nn=2900
!integer, parameter :: n2big=7320,n3big=8520,nw=3901,ne=4800,ns=4101,nn=5000
!integer, parameter :: n2big=7320,n3big=8520,nw=4001,ne=4100,ns=3251,nn=3350
   integer, parameter :: n2big = 7320, n3big = 8520, nw = 1, ne = n2big, ns = 1, nn = n3big

CONTAINS

   SUBROUTINE INITIALIZEDOMAIN(n2, n3, nzg, filetopo)
      integer :: n2, n3, nzg
      integer :: n, nmax_y, nmax_x
      integer :: tasktype
      integer :: request
      integer, allocatable :: req(:), stats(:, :)
      character(len=300) :: filetopo
      integer, dimension(2) :: array_of_sizes, array_of_subsizes, array_of_starts

      call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

      if (numtasks .eq. 1) then
         nini_x(0) = 1
         nend_x(0) = n2
         nini_y(0) = 1
         nend_y(0) = n3
         return
      end if

      rcountblocksmall(0) = n2
      rcountblock(0) = 0
      disp(0) = 0

      rcountblocksmall(numtasks - 1) = n2
      rcountblock(numtasks - 1) = 0
      disp(numtasks - 1) = 0

      call MPI_TYPE_CONTIGUOUS(numtasks, MPI_INTEGER, tasktype, ierr)
      call MPI_Type_commit(tasktype, ierr)

      if (pid .eq. 0) then

!task 0 and task numtasks-1 are not used for normal calculations

         call dividedomain(n2, n3, nini_x, nend_x, nini_y, nend_y, n_right, n_left, n_up, n_down, numtasks, filetopo)
!   call dividedomainold(n2,n3,nini_y,filetopo)

         do n = 1, numtasks - 1
            call MPI_send(nini_x(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(nend_x(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(nini_y(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(nend_y(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(n_right(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(n_left(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(n_up(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

         do n = 1, numtasks - 1
            call MPI_send(n_down(0), 1, tasktype, n, 1, MPI_COMM_WORLD, ierr)
         end do

      else
         call MPI_recv(nini_x(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(nend_x(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(nini_y(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(nend_y(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(n_right(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(n_left(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(n_up(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(n_down(0), 1, tasktype, 0, 1, MPI_COMM_WORLD, status, ierr)
      end if

      call MPI_TYPE_FREE(tasktype, ierr)

!do n=2,numtasks-2
!  nend_y(n-1)=nini_y(n)+1
!enddo

!gmmdeclare pieces to be send and received

      do n = 1, numtasks - 2

         nmax_x = nend_x(n) - nini_x(n) + 1
         nmax_y = nend_y(n) - nini_y(n) + 1

         if (pid .eq. 0) write (6, *) nini_x(n), nend_x(n), nmax_x, n_left(n), n_right(n), n
         if (pid .eq. 0) write (6, *) nini_y(n), nend_y(n), nmax_y, n_up(n), n_down(n), n

         rcountblocksmall(n) = (nmax_x - 2)*(nmax_y - 2)
!  if(n.eq.1.or.n.eq.numtasks-2)rcountblocksmall(n)=n2*(nmax_y-1)
         rcountblock(n) = nmax_x*nmax_y

         disp(n) = rcountblocksmall(n - 1) + disp(n - 1)

         call MPI_TYPE_CONTIGUOUS(nmax_x*nmax_y, MPI_REAL, domblock(n), ierr)
         call MPI_Type_commit(domblock(n), ierr)
         call MPI_TYPE_CONTIGUOUS(nmax_x*nmax_y, MPI_INTEGER, domblock2dint(n), ierr)
         call MPI_Type_commit(domblock2dint(n), ierr)
         call MPI_TYPE_CONTIGUOUS(2*nmax_x*nmax_y, MPI_INTEGER, domblock3dint(n), ierr)
         call MPI_Type_commit(domblock3dint(n), ierr)
         call MPI_TYPE_CONTIGUOUS(rcountblocksmall(n), MPI_REAL, domblocksmall(n), ierr)
         call MPI_Type_commit(domblocksmall(n), ierr)
         call MPI_TYPE_CONTIGUOUS(rcountblocksmall(n), MPI_INTEGER, domblocksmallint(n), ierr)
         call MPI_Type_commit(domblocksmallint(n), ierr)
         call MPI_TYPE_CONTIGUOUS(nmax_x*nmax_y, MPI_BYTE, domblockbyte(n), ierr)
         call MPI_Type_commit(domblockbyte(n), ierr)

         call MPI_TYPE_CONTIGUOUS(nmax_x, MPI_REAL, borderupdown(n), ierr)
         call MPI_Type_commit(borderupdown(n), ierr)

         call MPI_TYPE_CONTIGUOUS(nmax_y, MPI_REAL, borderleftright(n), ierr)
         call MPI_Type_commit(borderleftright(n), ierr)

         call MPI_TYPE_VECTOR(nmax_x, 1, 1, MPI_REAL, updownhalo(n), ierr)
         call MPI_Type_commit(updownhalo(n), ierr)

         call MPI_TYPE_VECTOR(nmax_y, 1, nmax_x, MPI_REAL, leftrighthalo(n), ierr)
         call MPI_Type_commit(leftrighthalo(n), ierr)

         array_of_sizes(1) = n2
         array_of_sizes(2) = n3
         array_of_subsizes(1) = nmax_x
         array_of_subsizes(2) = nmax_y
         array_of_starts(1) = nini_x(n) - 1
         array_of_starts(2) = nini_y(n) - 1
   call MPI_TYPE_CREATE_SUBARRAY(2,array_of_sizes,array_of_subsizes,array_of_starts,MPI_ORDER_FORTRAN,MPI_REAL,arraysection(n),ierr)
         call MPI_Type_commit(arraysection(n), ierr)
call MPI_TYPE_CREATE_SUBARRAY(2,array_of_sizes,array_of_subsizes,array_of_starts,MPI_ORDER_FORTRAN,MPI_INTEGER,arraysectionint(n),ierr)
         call MPI_Type_commit(arraysectionint(n), ierr)

      end do

      call MPI_Type_CONTIGUOUS(n2, MPI_REAL, columntype, ierr)
      call MPI_Type_commit(columntype, ierr)

      call MPI_Type_CONTIGUOUS(2*n2, MPI_REAL, columntype2, ierr)
      call MPI_Type_commit(columntype2, ierr)

   end subroutine initializedomain
!**********************************************************************

   subroutine dividedomainold(n1, n2, ini, filetopo)
      integer :: n1, n2
      integer :: ini(0:numtasks - 1)
      real, dimension(n2big, n3big) :: varreadbig
      real, dimension(n1, n2) :: varread
      integer, dimension(n2) :: ncells
      integer :: ntotal, ncount, i, j, n
      real :: nperpid
      character(len=*) :: filetopo

      write (6, *) 'reading soil data'
!read in topo data to use as mask

      open (21, file=filetopo &
            , access='direct', convert='big_endian', recl=4*n2big*n3big)

      read (21, rec=9) ((varreadbig(i, j), i=1, n2big), j=1, n3big)

      varread(1:n1, 1:n2) = varreadbig(nw:ne, ns:nn)

      close (21)

!varread=0.
      write (6, *) 'test soil data', varread(2, 2), varread(n1/2, n2/2)

      ntotal = count(varread >= -1.E+05)

      nperpid = float(ntotal)/float(numtasks - 2)

      write (6, *) 'total number of land cells', ntotal

      ncells = count(varread >= -1.E+05, 1)

      ncount = 0
      ini(1) = 1
      n = 2
      do j = 1, n2
         ncount = ncount + ncells(j)
         if (ncount .ge. nint(float(n - 1)*nperpid)) then
!              write(6,*)ncount,ncount-ncells(j-1),j,n
            ini(n) = j - 1
!              ncount=ncells(j)
            n = n + 1
         end if
         if (n .eq. numtasks - 1) exit
      end do

!      do n=0,numtasks-1
!         write(6,*)nini_y(n),n
!      enddo

!test

!      ncount=0
!      n=0
!      do j=nini_y(10),nini(11)-1
!        n=ncells(j)+n
!        do i=1,n1
!         if(varread(i,j).gt.0.5)ncount=ncount+1
!        enddo
!      enddo

!     write(6,*)ncount,ntotal/numtasks,n

   end subroutine dividedomainold
!     ******************************************************************
   subroutine dividedomain(n2, n3, is, ie, js, je, n_right, n_left, n_up, n_down, nprocmax, filetopo)
      integer :: n2, n3, n, nprocmax, nproc, nx_proc, ny_proc, iter, tag, totproc
      integer, dimension(0:npmax) :: is, ie, js, je, n_right, n_left, n_up, n_down
      integer, allocatable, dimension(:) :: i_ini, i_end, j_ini, j_end, procmask
      real, dimension(n2big, n3big) :: varreadbig
      real, dimension(n2, n3) :: varread
      integer, dimension(n2, n3) :: landmask, procid
      integer :: i, j, kk
      character(len=*) :: filetopo

      write (6, *) 'reading soil data'
!read in topo data to use as mask

      open (21, file=filetopo &
            , access='direct', convert='big_endian', recl=4*n2big*n3big)

      read (21, rec=9) ((varreadbig(i, j), i=1, n2big), j=1, n3big)

      varread(1:n2, 1:n3) = varreadbig(nw:ne, ns:nn)

      close (21)

      landmask = 0
      where (varread >= -1.E+05) landmask = 1

!now start the division from the max number of processors
!keep the same proportion in x and y as the dimensions
!processor 0 and numtasks-1 are reserved for input output

      nproc = nprocmax - 2
      iter = 0

      DO

         iter = iter + 1

         ny_proc = sqrt(float(nproc*n3/n2))
         nx_proc = nproc/ny_proc

         nproc = nx_proc*ny_proc

         write (6, *) 'nx_proc,ny_proc,nproc', nx_proc, ny_proc, nproc, nprocmax

         allocate (i_ini(nproc))
         allocate (i_end(nproc))
         allocate (j_ini(nproc))
         allocate (j_end(nproc))
         allocate (procmask(nproc))

!number patches from bottom left corner, right and up, and get dimensions, for now no haloes
         n = 0
         do j = 1, ny_proc
            do i = 1, nx_proc
               n = n + 1
               j_ini(n) = n3/ny_proc*(j - 1) + 1
               j_end(n) = j_ini(n) + n3/ny_proc - 1
               if (j .eq. ny_proc) j_end(n) = n3
               i_ini(n) = n2/nx_proc*(i - 1) + 1
               i_end(n) = i_ini(n) + n2/nx_proc - 1
               if (i .eq. nx_proc) i_end(n) = n2
            end do
         end do

!now check how many patches are entirely ocean, including the haloes (procmask=0)
         procmask = 0

         DO n = 1, nproc
         do j = max(j_ini(n) - 1, 1), min(j_end(n) + 1, n3)
            if (procmask(n) .eq. 1) exit
            do i = max(i_ini(n) - 1, 1), min(i_end(n) + 1, n2)
               if (landmask(i, j) .eq. 1) then
                  procmask(n) = 1
                  exit
               end if
            end do
         end do
         END DO

         write (6, *) 'number of land patches', sum(procmask)

!now number patches, -1 is for ocean patches
         DO n = 1, nproc
            do j = j_ini(n), j_end(n)
               do i = i_ini(n), i_end(n)
                  if (procmask(n) .eq. 0) then
                     procid(i, j) = -1
                  else
                     procid(i, j) = n
                  end if
               end do
            end do
         END DO

         totproc = sum(procmask)
         deallocate (i_ini, i_end, j_ini, j_end, procmask)
         if (totproc .gt. nprocmax - 2) exit

!if there are less land patches than processors available, try making them smaller until possible
!make the same division increasing the number of processors one at a time
         nproc = nprocmax - 2 + iter

         write (6, *) "nproc=", nproc
      END DO

      write (6, *) "nproc in final iteration=", totproc

!recalculate going back one iteration, before the max number of available processors is exceeded
      nproc = nproc - 1

      ny_proc = sqrt(float(nproc*n3/n2))
      nx_proc = nproc/ny_proc

      nproc = nx_proc*ny_proc

      write (6, *) nx_proc, ny_proc, nproc

      allocate (i_ini(nproc))
      allocate (i_end(nproc))
      allocate (j_ini(nproc))
      allocate (j_end(nproc))
      allocate (procmask(nproc))

      n = 0
      do j = 1, ny_proc
         do i = 1, nx_proc
            n = n + 1
            j_ini(n) = n3/ny_proc*(j - 1) + 1
            j_end(n) = j_ini(n) + n3/ny_proc - 1
            if (j .eq. ny_proc) j_end(n) = n3
            i_ini(n) = n2/nx_proc*(i - 1) + 1
            i_end(n) = i_ini(n) + n2/nx_proc - 1
            if (i .eq. nx_proc) i_end(n) = n2
         end do
      end do

      procmask = 0

      DO n = 1, nproc
      do j = max(j_ini(n) - 1, 1), min(j_end(n) + 1, n3)
         if (procmask(n) .eq. 1) exit
         do i = max(i_ini(n) - 1, 1), min(i_end(n) + 1, n2)
            if (landmask(i, j) .eq. 1) then
               procmask(n) = 1
               exit
            end if
         end do
      end do
      END DO

      write (6, *) sum(procmask)

      totproc = sum(procmask)

      write (6, *) "final nproc=", totproc

!now renumber the processors skipping the ocean ones and get the final dimensions with the haloes

      is = n2
      ie = 1
      js = n3
      je = 1
      n_right = -1
      n_left = -1
      n_up = -1
      n_down = -1

      procid = -1
!now renumber land patches and get dimensions

      kk = 0
      do n = 1, nproc
         if (procmask(n) .eq. 0) then
            tag = -1
         else
            kk = kk + 1
            tag = kk
            is(tag) = i_ini(n)
            ie(tag) = i_end(n)
            js(tag) = j_ini(n)
            je(tag) = j_end(n)
         end if

         do j = j_ini(n), j_end(n)
            do i = i_ini(n), i_end(n)
               procid(i, j) = tag
            end do
         end do
      end do

      do n = 1, totproc
         if (je(n) + 1 .gt. n3) then
            n_up(n) = -1
         else
            n_up(n) = procid(is(n), je(n) + 1)
         end if
         if (js(n) - 1 .lt. 1) then
            n_down(n) = -1
         else
            n_down(n) = procid(is(n), js(n) - 1)
         end if

         if (ie(n) + 1 .gt. n2) then
            n_right(n) = -1
         else
            n_right(n) = procid(ie(n) + 1, je(n))
         end if
         if (is(n) - 1 .lt. 1) then
            n_left(n) = -1
         else
            n_left(n) = procid(is(n) - 1, js(n))
         end if

      end do

      do n = 1, totproc
         is(n) = max(is(n) - 1, 1)
         ie(n) = min(ie(n) + 1, n2)
         js(n) = max(js(n) - 1, 1)
         je(n) = min(je(n) + 1, n3)
      end do

      deallocate (i_ini, i_end, j_ini, j_end, procmask)

   end subroutine dividedomain
!     ******************************************************************
   subroutine SENDBORDERS(n2, js, je, wtd, reqsu, reqsd, reqru, reqrd)
      integer :: n2, js, je, reqsu, reqsd, reqru, reqrd
      real, dimension(n2, js:je):: wtd

      if (pid .eq. 1) then
         call MPI_isend(wtd(1, je - 1), 1, columntype, 2, 200, MPI_COMM_WORLD, reqsu, ierr)
         call MPI_irecv(wtd(1, je), 1, columntype, 2, 201, MPI_COMM_WORLD, reqru, ierr)

      elseif (pid .eq. numtasks - 2) then
         call MPI_isend(wtd(1, js + 1), 1, columntype, pid - 1, 201, MPI_COMM_WORLD, reqsd, ierr)
         call MPI_irecv(wtd(1, js), 1, columntype, pid - 1, 200, MPI_COMM_WORLD, reqrd, ierr)

      elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
         call MPI_isend(wtd(1, je - 1), 1, columntype, pid + 1, 200, MPI_COMM_WORLD, reqsu, ierr)
         call MPI_isend(wtd(1, js + 1), 1, columntype, pid - 1, 201, MPI_COMM_WORLD, reqsd, ierr)
         call MPI_irecv(wtd(1, js), 1, columntype, pid - 1, 200, MPI_COMM_WORLD, reqrd, ierr)
         call MPI_irecv(wtd(1, je), 1, columntype, pid + 1, 201, MPI_COMM_WORLD, reqru, ierr)

      end if

   end subroutine sendborders
!     ******************************************************************
   subroutine SENDBORDERSFLOOD(n2, js, je, wtd, reqsu, reqsd, reqru, reqrd)
      integer :: n2, js, je, reqsu, reqsd, reqru, reqrd
      real, dimension(n2, js:je):: wtd

      if (pid .eq. 1) then
         call MPI_isend(wtd(1, je), 1, columntype, 2, 200, MPI_COMM_WORLD, reqsu, ierr)
         call MPI_irecv(wtd(1, je - 1), 1, columntype, 2, 201, MPI_COMM_WORLD, reqru, ierr)

      elseif (pid .eq. numtasks - 2) then
         call MPI_isend(wtd(1, js), 1, columntype, pid - 1, 201, MPI_COMM_WORLD, reqsd, ierr)
         call MPI_irecv(wtd(1, js + 1), 1, columntype, pid - 1, 200, MPI_COMM_WORLD, reqrd, ierr)

      elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
         call MPI_isend(wtd(1, je), 1, columntype, pid + 1, 200, MPI_COMM_WORLD, reqsu, ierr)
         call MPI_isend(wtd(1, js), 1, columntype, pid - 1, 201, MPI_COMM_WORLD, reqsd, ierr)
         call MPI_irecv(wtd(1, js + 1), 1, columntype, pid - 1, 200, MPI_COMM_WORLD, reqrd, ierr)
         call MPI_irecv(wtd(1, je - 1), 1, columntype, pid + 1, 201, MPI_COMM_WORLD, reqru, ierr)

      end if

   end subroutine sendbordersflood

!     ******************************************************************
   subroutine SENDBORDERS4blocking(is, ie, js, je, wtd)
      integer :: is, ie, js, je
      real, dimension(is:ie, js:je):: wtd

!first up and down then left and right, to ensure halo corners are passed
      if (n_down(pid) .gt. 0) then
         call MPI_send(wtd(is, js + 1), 1, updownhalo(pid), n_down(pid), 200, MPI_COMM_WORLD, ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_send(wtd(is, je - 1), 1, updownhalo(pid), n_up(pid), 201, MPI_COMM_WORLD, ierr)
      end if

      if (n_down(pid) .gt. 0) then
         call MPI_recv(wtd(is, js), 1, updownhalo(pid), n_down(pid), 201, MPI_COMM_WORLD, status, ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_recv(wtd(is, je), 1, updownhalo(pid), n_up(pid), 200, MPI_COMM_WORLD, status, ierr)
      end if

      if (n_left(pid) .gt. 0) then
         call MPI_send(wtd(is + 1, js), 1, leftrighthalo(pid), n_left(pid), 202, MPI_COMM_WORLD, ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_send(wtd(ie - 1, js), 1, leftrighthalo(pid), n_right(pid), 203, MPI_COMM_WORLD, ierr)
      end if

      if (n_left(pid) .gt. 0) then
         call MPI_recv(wtd(is, js), 1, leftrighthalo(pid), n_left(pid), 203, MPI_COMM_WORLD, status, ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_recv(wtd(ie, js), 1, leftrighthalo(pid), n_right(pid), 202, MPI_COMM_WORLD, status, ierr)
      end if

   end subroutine sendborders4blocking

!     ******************************************************************
   subroutine SENDBORDERS4(is, ie, js, je, wtd)
      integer :: is, ie, js, je
      real, dimension(is:ie, js:je):: wtd
      integer, dimension(2) :: requ, reqd, reql, reqr
      integer :: status2(MPI_STATUS_SIZE, 2)

!first up and down then left and right, to ensure halo corners are passed

      if (n_down(pid) .gt. 0) then
         call MPI_irecv(wtd(is, js), 1, updownhalo(pid), n_down(pid), 201, MPI_COMM_WORLD, reqd(1), ierr)
         call MPI_isend(wtd(is, js + 1), 1, updownhalo(pid), n_down(pid), 200, MPI_COMM_WORLD, reqd(2), ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_irecv(wtd(is, je), 1, updownhalo(pid), n_up(pid), 200, MPI_COMM_WORLD, requ(1), ierr)
         call MPI_isend(wtd(is, je - 1), 1, updownhalo(pid), n_up(pid), 201, MPI_COMM_WORLD, requ(2), ierr)
      end if

      if (n_down(pid) .gt. 0) then
         call MPI_waitall(2, reqd, status2, ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_waitall(2, requ, status2, ierr)
      end if

      if (n_left(pid) .gt. 0) then
         call MPI_irecv(wtd(is, js), 1, leftrighthalo(pid), n_left(pid), 203, MPI_COMM_WORLD, reql(1), ierr)
         call MPI_isend(wtd(is + 1, js), 1, leftrighthalo(pid), n_left(pid), 202, MPI_COMM_WORLD, reql(2), ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_irecv(wtd(ie, js), 1, leftrighthalo(pid), n_right(pid), 202, MPI_COMM_WORLD, reqr(1), ierr)
         call MPI_isend(wtd(ie - 1, js), 1, leftrighthalo(pid), n_right(pid), 203, MPI_COMM_WORLD, reqr(2), ierr)
      end if

      if (n_left(pid) .gt. 0) then
         call MPI_waitall(2, reql, status2, ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_waitall(2, reqr, status2, ierr)
      end if

   end subroutine sendborders4

!     ******************************************************************
   subroutine SENDBORDERSFLOOD4(is, ie, js, je, var, borderu, borderd, borderl, borderr)
      integer :: is, ie, js, je
      real, dimension(is:ie, js:je):: var
      real, dimension(is:ie) :: borderu, borderd
      real, dimension(js:je) :: borderl, borderr
      integer, dimension(2) :: requ, reqd, reql, reqr
      integer :: status2(MPI_STATUS_SIZE, 2)

      borderu = 0.
      borderd = 0.
      borderl = 0.
      borderr = 0.

!first up and down then left and right, to ensure halo corners are passed

      if (n_down(pid) .gt. 0) then
         call MPI_irecv(borderd(is), 1, borderupdown(pid), n_down(pid), 201, MPI_COMM_WORLD, reqd(1), ierr)
         call MPI_isend(var(is, js), 1, updownhalo(pid), n_down(pid), 200, MPI_COMM_WORLD, reqd(2), ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_irecv(borderu(is), 1, borderupdown(pid), n_up(pid), 200, MPI_COMM_WORLD, requ(1), ierr)
         call MPI_isend(var(is, je), 1, updownhalo(pid), n_up(pid), 201, MPI_COMM_WORLD, requ(2), ierr)
      end if

      if (n_down(pid) .gt. 0) then
         call MPI_waitall(2, reqd, status2, ierr)
      end if
      if (n_up(pid) .gt. 0) then
         call MPI_waitall(2, requ, status2, ierr)
      end if

!to pass the corners from the corner domains use the corners of the halo, which really have no use
      var(is, je) = borderu(is)
      var(ie, je) = borderu(ie)
      var(is, js) = borderd(is)
      var(ie, js) = borderd(ie)

      if (n_left(pid) .gt. 0) then
         call MPI_irecv(borderl(js), 1, borderleftright(pid), n_left(pid), 203, MPI_COMM_WORLD, reql(1), ierr)
         call MPI_isend(var(is, js), 1, leftrighthalo(pid), n_left(pid), 202, MPI_COMM_WORLD, reql(2), ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_irecv(borderr(js), 1, borderleftright(pid), n_right(pid), 202, MPI_COMM_WORLD, reqr(1), ierr)
         call MPI_isend(var(ie, js), 1, leftrighthalo(pid), n_right(pid), 203, MPI_COMM_WORLD, reqr(2), ierr)
      end if

      if (n_left(pid) .gt. 0) then
         call MPI_waitall(2, reql, status2, ierr)
      end if
      if (n_right(pid) .gt. 0) then
         call MPI_waitall(2, reqr, status2, ierr)
      end if

   end subroutine sendbordersflood4

!     ******************************************************************

END MODULE MODULE_PARALLEL
