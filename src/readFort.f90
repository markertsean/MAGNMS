!--------------------------------------------------------------------------------
!                  read files pmss files from gadget runs
!
!                  (nx,ny,nz) - number of boxes in each direction
!                  rhalo  = width of buffer in mpch around each region
module setarrs

  integer*4, parameter                           :: MaxDomains =2048
  integer*4                                       :: nx, ny, nz, nbuffer  ! Number of domains
  integer*4                                       :: jstep                ! For loops
  integer*4                                       :: ntot, nrecord        ! Number of things
  integer*4                                       :: dirpathlength
  integer*4                                       :: PMssFirst, PMssLast  ! First and last values for PMss Files

  real*4                                        :: xl=0, xr=0, yl=0, yr=0, zl=0, zr=0, dbuffer  ! Boundaries
  real*4                                        :: aexpn, om0, oml0, hubble, ovdens ! Cosmology
  real*4                                        :: box, massone                     ! Simulation
  real*4,                Dimension(MaxDomains) :: xLeft,xRight, yLeft,yRight,zLeft,zRight
  real*4,  ALLOCATABLE, DIMENSION(:)          :: Xp, Yp, Zp


  character*30 :: formatstr
  character*60 :: dirpath

contains

subroutine readconfig

  character*80 :: name
  logical      :: ext, found
  real*4       :: qx, qy, qz
  integer*4    :: fileIndex

  !Creates write statement for file name, using length of dirpath
  if (dirpathlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',dirpathlength-5,',a,i3.3,a,i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',dirpathlength-5,',a,i3.3,a,i4.4,a)'
  endif


  found = .false.

  do jfile = PMssFirst, PMssLast

     !File to attempt to open, if it exists, read data
     write(name,formatstr) dirpath,'PMss.snap_',jstep,'.',jfile,'.DAT'

     Inquire(file=name,exist =ext)

     If(ext)Then

        fileIndex = jfile

        open(1,file=name,form='unformatted')
        read(1) aexpn,om0,oml0,hubble,box,massone
        read(1) node,nx,ny,nz,dbuffer,nbuffer
        read(1) xl,xr,yl,yr,zl,zr
        read(1) ntot
        read(1) nrecord
        close(1)

        write(*,*) ' '
        write(*,'(a,a)'            ) '     Found PMSS:  ',name
        write(*,'(a,f9.4,5x,a,2i5)') '         aexpn   =',aexpn
        write(*,'(a,f9.4,5x,a,2i5)') '         Om0     =',Om0
        write(*,'(a,f9.4,5x,a,2i5)') '         h       =',hubble
        write(*,'(a,f9.2,5x,a,2i5)') '         Box     =',Box
        write(*,'(a,1p,g13.5)'     ) '         Mass_1  =',MassOne
        write(*,'(a,3i5)'          ) '         Domain  =',nx,ny,nz
        write(*,'(a,f8.2)'         ) '         Buffer  =',dBuffer
        write(*,'(a,i10)'          ) '         Nbuffer =',nBuffer

        found = .true.

        exit

     end If ! File exists

  enddo ! File loop

  if( .not.found ) stop ' No files for this snapshot found'


  ! The size of a node on a side
  qx = Box / nx
  qy = Box / ny
  qz = Box / nz

  ! All of the boundaries for the boxes
  do i=1, nx
     xLeft(i) = ( i - 1 ) * qx
    xRight(i) =   i       * qx
  end do

  do i=1, ny
     yLeft(i) = ( i - 1 ) * qy
    yRight(i) =   i       * qy
  end do

  do i=1, nz
     zLeft(i) = ( i - 1 ) * qz
    zRight(i) =   i       * qz
  end do


  xl =  xLeft( fileIndex )
  xr = xRight( fileIndex )
  yl =  yLeft( fileIndex )
  yr = yRight( fileIndex )
  zl =  zLeft( fileIndex )
  zr = zRight( fileIndex )


end subroutine readconfig

end module setarrs



!jstep is snapshot number
!filestart is something like ../data/Part. and is used to write out files
!filestartlength is the length of filestart, needed to get the format for the name variable right. Was having issues opening without this
integer*8 function readpmss( initjstep, filestart, filestartlength, firstPMss, lastPMss )

  use setarrs

  integer*4,    intent(in) :: initjstep, filestartlength, firstPMss, lastPMss
  character*60, intent(in) :: filestart

  integer*4 :: i0,j0, k0
  integer*8 :: totparticles

  character*80 :: particle_name, header_name


  ! Save input variables to common block, configure common block variables
  jstep         = initjstep
  dirpath       = filestart
  dirpathlength = filestartlength
  PMssFirst     = firstPMss
  PMssLast      =  lastPMss

  ! Creates write statement for output file names, using length of filestart
  if (filestartlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',filestartlength,',i4.4,a,i4.4,a,i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',filestartlength,',i4.4,a,i4.4,a,i4.4,a)'
  endif

  ! Generates file names
  write( particle_name, formatstr) trim(filestart),jStep,'.',PMssFirst,'.',PMssLast,'.DAT'
  write(   header_name, formatstr) trim(filestart),jStep,'.',PMssFirst,'.',PMssLast,'.header.dat'


  call readconfig   ! Configures variables about the simulation into our common block




 ! Opens files
  open( 2, file=particle_name, form='formatted' )
  open( 3, file=  header_name, form='formatted' )

  totParticles = 0




  do k0=1,nz
  do j0=1,ny
  do i0=1,nx
            totparticles = totparticles + readfile( i0, j0, k0, totparticles)
  enddo ! i
  enddo ! j
  enddo ! k





  write(*,*) ' Total Number of particles =',totparticles
  close(2)
  write(*,'(2a)') '  Wrote file: ', particle_name

  write(3,'(a8,f12.6)'   ) 'a       ', aexpn
  write(3,'(a8,f12.6)'   ) 'Om_m    ', om0
  write(3,'(a8,f12.6)'   ) 'Om_l    ', oml0
  write(3,'(a8,f12.6)'   ) 'h       ', hubble
  write(3,'(a8,1p,e12.5)') 'M_part  ', massone
  write(3,'(a8,i12)'     ) 'node    ', node
  write(3,'(a8,f12.6)'   ) 'Buffer  ', dbuffer
  write(3,'(a8,f12.6)'   ) 'xl      ', xl
  write(3,'(a8,f12.6)'   ) 'xr      ', xr
  write(3,'(a8,f12.6)'   ) 'yl      ', yl
  write(3,'(a8,f12.6)'   ) 'yr      ', yr
  write(3,'(a8,f12.6)'   ) 'zl      ', zl
  write(3,'(a8,f12.6)'   ) 'zr      ', zr
  write(3,'(a8,i12)'     ) 'N_part  ', totparticles
  close(3)
  write(*,'(2a)') '  Wrote file: ', header_name

  write(*,*) ' '

  readpmss = totparticles
  return

end function readpmss



!---------------------------------------------------------
!                                        read particles
function readfile( i0, j0, k0, totparticles )
!--------------------------------------------------------------

  use setarrs

  integer*4,    intent(in) :: i0, j0, k0
  integer*8,    intent(in) :: totparticles

  integer*4     ::  nNbuffer                            ! number of domains
  integer*4     ::  i, ierr, icount= 0
  integer*8     ::  nin, node, jcnt

  real*8        :: xxl, xxr, yyl, yyr, zzl, zzr, dB

  character*100 :: name
  logical       :: ext,found, inside, validNode


  real*4,    allocatable, dimension(:) :: Xb, Yb, Zb, Vxb, Vyb, Vzb
  integer*8, allocatable, dimension(:) :: Ib


  node   = i0+(j0-1)*nx +(k0-1)*nx*ny
  nin    = 0
  icount = 0

  validNode = ( node.le.PMssLast).and.( node.ge.PMssfirst)

  !File to attempt to open, if it exists, read data
  write(name,formatstr) dirpath,'PMss.snap_',jstep,'.',node,'.DAT'


  inquire(file=trim(name),exist =ext)
  if(ext.and.validNode)then

     write(*,*) "  Reading file: ",name

     ! Reads header info from file
     open(1,file=trim(Name),form='unformatted')
     read(1) aexpn,om0,oml0,hubble,box,massone
     read(1) node1,anx1,any1,anz1,dbuffer1,anbuffer1
     read(1) txl,txr,tyl,tyr,tzl,tzr
     read(1) ntot
     if(node1 /= node)stop '-- Error in Header --'


     node    = node1
     dbuffer = dbuffer1
     nbuffer = anbuffer1


     xxl =  xLeft( i0 )
     xxr = xRight( i0 )
     yyl =  yLeft( j0 )
     yyr = yRight( j0 )
     zzl =  zLeft( k0 )
     zzr = zRight( k0 )

     xl =  min( xl, xxl )
     yl =  min( yl, yyl )
     zl =  min( zl, zzl )
     xr =  max( xr, xxr )
     yr =  max( yr, yyr )
     zr =  max( zr, zzr )

     write(*,'(5x,a,i5,3x,a,3(2f8.2,2x),3x,a,i10)') 'node= ',node1, &
              'limits= ',xxl,xxr,yyl,yyr,zzl,zzr,' Ntotal in node=',ntot


     do
        read(1,iostat=ierr)nrecord         ! number of particles in this record
        if(ierr/=0)exit

        allocate(  Xb(nrecord),  Yb(nrecord),  Zb(nrecord) )
        allocate( Vxb(nrecord), Vyb(nrecord), Vzb(nrecord) )
        allocate(  Ib(Nrecord))


        nin  = nin  + nrecord


        read(1)(  Xb(m),  Yb(m),  Zb(m), &
                 Vxb(m), Vyb(m), Vzb(m), Ib(m), m=1, Nrecord)


          do i=1,nrecord

            inside = (xxl <= xb(i)) .and. (xb(i)<xxr) .and. &
                     (yyl <= yb(i)) .and. (yb(i)<yyr) .and. &
                     (zzl <= zb(i)) .and. (zb(i)<zzr)


            if( inside ) then
              icount = icount + 1

              write(2,'(3f12.6)') Xb(i), Yb(i), Zb(i)

            endif ! inside

          enddo ! nrecord loop

        deallocate(  xb,  yb,  zb)
        deallocate( vxb, vyb, vzb)
        deallocate(  Ib )

     enddo      ! read record loop
     close(1)

     if(nin /= ntot)stop ' -- Wrong number of particles --'


  endif ! File exists

  readfile = icount
  return

end function readfile
