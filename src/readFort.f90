!--------------------------------------------------------------------------------
!                  read files pmss files from gadget runs
!
!                  (nx,ny,nz) - number of boxes in each direction
!                  rhalo  = width of buffer in mpch around each region

!jstep is snapshot number
!filestart is something like ../data/Part. and is used to write out files
!filestartlength is the length of filestart, needed to get the format for the name variable right. Was having issues opening without this
function readpmss( jstep, filestart, filestartlength )

  integer*4,    intent(in) :: jstep, filestartlength
  character*60, intent(in) :: filestart

  integer*4 :: nx = 8
  integer*4 :: ny = 8
  integer*4 :: nz = 8
  integer*4 :: i0,j0, k0, totparticles

  character*80:: name3, name2
  character*20 :: formatstr

  !Creates write statement for output file names, using length of filestart
  if (filestartlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',filestartlength,',i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',filestartlength,',i4.4,a)'
  endif

call srand(1165313)

  !Opens file to write particles to
  write(name2,formatstr) trim(filestart),jStep,'.DAT'
  open(2,file=name2,form='formatted')

  !Opens header file
  write(name3,formatstr) trim(filestart),jStep,'.header.dat'
  open(3,file=name3,form='formatted')

  totParticles = 0

  do k0=1,nz
  do j0=1,ny
  do i0=1,nx
            totparticles = totparticles + readfile(i0,j0,k0,nx,ny,nz,jstep,totparticles,filestart,filestartlength-5)
  enddo ! k
  enddo ! j
  enddo ! i

  write(*,*) ' Total Number of particles =',totparticles
  close(2)
  write(*,'(2a)') 'Wrote file: ',name2

  write(3,'(a8,i12)') 'N_part  ',totparticles
  close(3)
  write(*,'(2a)') 'Wrote file: ',name3


  readpmss = totparticles
  return

end function readpmss

!---------------------------------------------------------
!                                        read particles
function readfile(i0,j0,k0,nx,ny,nz,jstep,totparticles,dirpath,dirpathlength)
!--------------------------------------------------------------

  integer*4,    intent(in) :: i0, j0, k0, jstep, totparticles, dirpathlength
  integer*4,    intent(in) :: nx, ny, nz
  character*60, intent(in) :: dirpath


  integer*4  ::  nNbuffer                            ! number of domains
  integer*4  ::  np, ntot                            ! Number of particles
  integer*4  ::   i, ierr, icount = 0
  integer*8  :: nin, node

  real*8     :: xxl,xxr,yyl,yyr,zzl,zzr,    dB
  real*4     ::  xl,xr,   yl,yr,     zl,zr, dbuffer  ! boundaries of domain
  real*4     ::  aexpn,     om0,      oml0, hubble   ! cosmology
  real*4     ::    box, massone, fraction           !  current simulation

  real*4,       allocatable,   dimension(:) :: xpp,  ypp,  zpp

  character*100 :: name
  character*30  :: formatstr
  logical       :: ext,found, inside


  nbuffer = 1024*1024!???????!512*512; !=2^18=2^9*2^9
  nin = 0

  allocate(xpp(nbuffer),ypp(nbuffer),zpp(nbuffer))
  node  = i0+(j0-1)*nx +(k0-1)*nx*ny


  !Creates write statement for file name, using length of dirpath
  if (dirpathlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',dirpathlength,',a,i3.3,a,i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',dirpathlength,',a,i3.3,a,i4.4,a)'
  endif

  !File to attempt to open, if it exists, read data
  write(name,formatstr) dirpath,'PMss.snap_',jstep,'.',node,'.DAT'

  inquire(file=trim(name),exist =ext)
  if(ext)then

    write(*,*) " Reading file: ",name

     !Reads header info from file
     open(1,file=trim(Name),form='unformatted')
     read(1) aexpn,om0,oml0,hubble,box,massone
     read(1) node1,anx1,any1,anz1,dbuffer1,anbuffer1
     read(1) xl,xr,yl,yr,zl,zr
     read(1) ntot
     if(node1 /= node)stop '-- Error in Header --'
     node    = node1
     dbuffer = dbuffer1
     nbuffer = anbuffer1

     write(*,'(5x,a,i5,3x,a,3(2f8.2,2x),3x,a,i10)') 'node= ',node, &
              'limits= ',xl,xr,yl,yr,zl,zr,' Ntotal in node=',ntot

    !First time through, write data to header file
    if (totparticles.eq.0) then
      write(3,'(a8,f12.6)') 'a       ', aexpn
      write(3,'(a8,f12.6)') 'Om_m    ',om0
      write(3,'(a8,f12.6)') 'Om_l    ',oml0
      write(3,'(a8,f12.6)') 'h       ',hubble
      write(3,'(a8,1p,e12.5)') 'M_part  ',massone
      write(3,'(a8,i12)'  ) 'node    ',node
      write(3,'(a8,f12.6)') 'Buffer  ',dbuffer
    endif

     !Buffer size from file
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

             inside = (xxl <= xpp(i)) .and. (xpp(i)<xxr) .and. ( 1e-8 < xpp(i) ) .and. &
                      (yyl <= ypp(i)) .and. (ypp(i)<yyr) .and. ( 1e-8 < ypp(i) ) .and. &
                      (zzl <= zpp(i)) .and. (zpp(i)<zzr) .and. ( 1e-8 < zpp(i) )

if(inside.and.(rand()<0.00001))then
                   icount = icount +1
                   write(2,'(3f12.6)') xpp(i),ypp(i),zpp(i)
              endif ! inside

          enddo
     enddo
     close(1)

     if(nin /= ntot)stop ' -- Wrong number of particles --'

  endif

  readfile = icount
  return

end function readfile