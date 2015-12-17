!--------------------------------------------------------------------------------
!                  read files pmss files from gadget runs
!
!                  (nx,ny,nz) - number of boxes in each direction
!                  rhalo  = width of buffer in mpch around each region
module setarrs

integer*4  ::    nx, ny, nz, nbuffer                               ! number of domains
real*4     ::     xl,xr,yl,yr,zl,zr,dbuffer   ! boundaries of domain
                     ! ------------------------- cosmology ----------------------
real*4     ::  aexpn, om0,oml0,hubble        ! cosmology
real*4     ::  box,massone,fraction                   !  current simulation
integer*4  ::  jstep
integer*4  ::  np,ntot                          ! number of particles

                     ! -------------------------  particles ----------------------

real*4,       allocatable,   dimension(:) ::             &
                                      xpp,  ypp,  zpp!,   &      !   coords
!                                      wxp,  wyp,  wzp
integer*8       ::  icount

contains
!--------------------------------------------------------------------------------
!                          define configuration of files
!                          read control info from the first file

subroutine readconfig
  character*80 :: name
  logical      :: ext,found

  write(*,'(a,$)')' enter snapshot number ='
  read(*,*)jstep
                         ! read header information from the first existing
  found = .false.        ! default false, if any file found becomes true
  do jfile = 1, 512!?
     write(name,'(a,i3.3,a,i4.4,a)') 'pmss.snap_',jstep,'.',jfile,'.dat'
     inquire(file=name,exist =ext)
     if(ext)then
        open(1,file=name,form='unformatted')
        read(1) aexpn,om0,oml0,hubble,box,massone     !global variables all, but anx-z
        read(1) node,anx,any,anz,dbuffer,anbuffer
        read(1) xl,xr,yl,yr,zl,zr
        write(*,'(a,f9.4,5x,a,2i5)')' aexpn =',aexpn
        write(*,'(a,f9.4,5x,a,2i5)')' om0   =',om0
        write(*,'(a,f9.4,5x,a,2i5)')' h     =',hubble
        write(*,'(a,f9.2,5x,a,2i5)')' box   =',box
        write(*,'(a,1p,g13.5)')     ' mass_1=',massone
        nx = anx; ny = any; nz = anz;
        nbuffer = anbuffer
        write(*,'(a,3i5)')          ' domain=',nx,ny,nz
        write(*,'(a,f8.2,a,i10)')   ' buffer=',dbuffer,'(mpch)  nbuffer=',nbuffer
        found = .true.
        close(1)
        exit
     end if
  enddo

   if(.not.found)stop ' no files for this snapshot are found'
  write(name,'(a,i4.4,a)') 'part.',jstep,'.dat'
  open(2,file=name,form='formatted')
  write(2,'(a,f9.4,5x,a,f9.2)')'# aexpn=',aexpn,' box=',box

end subroutine readconfig



end module setarrs
!--------------------------------------------------------------------------------
!
subroutine readfiles
  use setarrs
  logical :: inside

  call readconfig                     ! set parameters, user enters snapshot number, reads single header
                                      ! reads global variables if finds file from snapshot
  nbuffer = 512*512;
  allocate(xpp(nbuffer),ypp(nbuffer),zpp(nbuffer))
!  !allocate(wxp(nbuffer),wyp(nbuffer),wzp(nbuffer))

  icount = 0       ! total selected particles

  nx = 8
  ny = 8
  nz = 8

   do k0=1,nz
   do j0=1,ny
   do i0=1,nx
             call readpmss(i0,j0,k0) !read file header for existing files, changing some globals
   enddo ! k
   enddo ! j
   enddo ! i

   write(*,*) ' total number of particles =',icount
   close(2)
   stop
end subroutine readfiles

!---------------------------------------------------------
!                                        read particles
subroutine readpmss(i0,j0,k0)
!--------------------------------------------------------------
  use setarrs
  character*80 :: name
  logical      :: ext,found, inside
  real*8       :: xxl,xxr,yyl,yyr,zzl,zzr,db

  node  = i0+(j0-1)*nx +(k0-1)*nx*ny
  write(name,'(a,i3.3,a,i4.4,a)') 'pmss.snap_',jstep,'.',node,'.dat'

  inquire(file=trim(name),exist =ext)
  if(ext)then
     open(1,file=trim(name),form='unformatted')
     read(1) aexpn,om0,oml0,hubble,box,massone
     read(1) node1,anx1,any1,anz1,dbuffer1,anbuffer1
     read(1) xl,xr,yl,yr,zl,zr
     read(1) ntot
     if(node1 /= node)stop '-- error in header --'

     write(*,'(5x,a,i5,3x,a,3(2f8.2,2x),3x,a,i10)') 'node= ',node, &
              'limits= ',xl,xr,yl,yr,zl,zr,' ntotal in node=',ntot

     nin = 0
     db  = dbuffer

     xxl = xl+db ; xxr = xr-db
     yyl = yl+db ; yyr = yr-db
     zzl = zl+db ; zzr = zr-db

     do
        read(1,iostat=ierr)nrecord         ! number of particles in this record

        if(ierr/=0)exit
        nin = nin + nrecord
        read(1)     (xpp(i),ypp(i),zpp(i),i=1,nrecord)
!        read(1)     (wxp(i),wyp(i),wzp(i),i=1,nrecord)
          do i=1,nrecord
             inside = (xxl <= xpp(i)).and. (xpp(i)<xxr) .and. &
                      (yyl <= ypp(i)).and. (ypp(i)<yyr) .and. &
                      (zzl <= zpp(i)).and. (zpp(i)<zzr)

             if(inside)then
                   icount = icount +1
!                   write(2,'(3f12.6,1p,4g14.6)') xpp(i),ypp(i),zpp(i),wxp(i),wyp(i),wzp(i)
                   write(2,'(3f12.6)') xpp(i),ypp(i),zpp(i)
              endif ! inside
          enddo
     enddo
     close(1)

     if(nin /= ntot)stop ' -- wrong number of particles --'

  endif
end subroutine readpmss

subroutine teststuff

  write(*,*) "FORTRAN BITCHES"
end subroutine teststuff
