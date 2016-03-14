!--------------------------------------------------------------------------------
!                  read files pmss files from gadget runs
!
!                  (nx,ny,nz) - number of boxes in each direction
!                  rhalo  = width of buffer in mpch around each region
module setarrs

integer*4 :: nx, ny, nz, nbuffer              ! Number of domains
integer*4 :: istep, jstep                     ! For loops
integer*4 :: ntot, iseed                      ! Number of things
integer*4 :: dirpathlength


real*4    :: xl, xr, yl, yr, zl, zr, dbuffer  ! Boundaries
real*4    :: aexpn, om0, oml0, hubble, ovdens ! Cosmology
real*4    :: box, massone                     ! Simulation


real*4, allocatable, dimension(:) :: xpp, ypp, zpp, &
                                       wxp, wyp, wzp


character*30 :: formatstr
character*60 :: dirpath


contains

subroutine readconfig

  character*80 :: name
  logical      :: ext, found


  !Creates write statement for file name, using length of dirpath
  if (dirpathlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',dirpathlength-5,',a,i3.3,a,i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',dirpathlength-5,',a,i3.3,a,i4.4,a)'
  endif


  found = .false.

  do jfile = 1, 1024

     !File to attempt to open, if it exists, read data
     write(name,formatstr) dirpath,'PMss.snap_',jstep,'.',jfile,'.DAT'

     Inquire(file=name,exist =ext)

     If(ext)Then

        open(1,file=name,form='unformatted')
        read(1) aexpn,om0,oml0,hubble,box,massone
        read(1) node,anx,any,anz,dbuffer,anbuffer
        read(1) xl,xr,yl,yr,zl,zr

        write(*,*) ' '
        write(*,'(a,a)'            ) '   Found PMSS: ',name
        write(*,'(a,f9.4,5x,a,2i5)') '     aexpn  =',aexpn
        write(*,'(a,f9.4,5x,a,2i5)') '     Om0    =',Om0
        write(*,'(a,f9.4,5x,a,2i5)') '     h      =',hubble
        write(*,'(a,f9.2,5x,a,2i5)') '     Box    =',Box
        write(*,'(a,1p,g13.5)')      '     Mass_1 =',MassOne
        nx = aNx; ny = aNy; nz = aNz;
        nBuffer = anBuffer
        write(*,'(a,3i5)')           '     Domain =',nx,ny,nz
        write(*,'(a,f8.2,a,i10)')    '     Buffer =',dBuffer,'(Mpch)  Nbuffer=',nBuffer
        found = .true.
        close(1)

        write(3,'(a8,f12.6)'   ) 'a       ', aexpn
        write(3,'(a8,f12.6)'   ) 'Om_m    ', om0
        write(3,'(a8,f12.6)'   ) 'Om_l    ', oml0
        write(3,'(a8,f12.6)'   ) 'h       ', hubble
        write(3,'(a8,1p,e12.5)') 'M_part  ', massone
        write(3,'(a8,i12)'     ) 'node    ', node
        write(3,'(a8,f12.6)'   ) 'Buffer  ', dbuffer

        exit

     end If ! File exists

  enddo ! File loop



  if( .not.found ) stop ' No files for this snapshot found'


end subroutine readconfig

end module setarrs



!jstep is snapshot number
!filestart is something like ../data/Part. and is used to write out files
!filestartlength is the length of filestart, needed to get the format for the name variable right. Was having issues opening without this
integer*8 function readpmss( initjstep, filestart, filestartlength )

  use setarrs

  integer*4,    intent(in) :: initjstep, filestartlength
  character*60, intent(in) :: filestart

  integer*4 :: i0,j0, k0
  integer*8 :: totparticles

  character*80 :: particle_name, header_name




  ! Save input variables to common block, configure common block variables
  jstep         = initjstep
  dirpath       = filestart
  dirpathlength = filestartlength

  ! Creates write statement for output file names, using length of filestart
  if (filestartlength.ge.10) then
    write(formatstr,'(a,i2,a)') '(a',filestartlength,',i4.4,a)'
  else
    write(formatstr,'(a,i1,a)') '(a',filestartlength,',i4.4,a)'
  endif

  ! Generates file names
  write( particle_name, formatstr) trim(filestart),jStep,'.DAT'
  write(   header_name, formatstr) trim(filestart),jStep,'.header.dat'


  call readconfig   ! Configures variables about the simulation into our common block



nz=8
ny=8
nx=8



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
  integer*8     ::  nin, node

  real*8        :: xxl,xxr,yyl,yyr,zzl,zzr,    dB
  real*4        :: txl,txr,tyl,tyr,tzl,tzr              ! temp boundaries on positions

  character*100 :: name
  logical       :: ext,found, inside


  node  = i0+(j0-1)*nx +(k0-1)*nx*ny

  !File to attempt to open, if it exists, read data
  write(name,formatstr) dirpath,'PMss.snap_',jstep,'.',node,'.DAT'
nbuffer=1024*1024

  allocate(xpp(nbuffer),ypp(nbuffer),zpp(nbuffer))
  allocate(wxp(nbuffer),wyp(nbuffer),wzp(nbuffer))

  nin    = 0
  icount = 0

  inquire(file=trim(name),exist =ext)
  if(ext)then

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



     write(*,'(5x,a,i5,3x,a,3(2f8.2,2x),3x,a,i10)') 'node= ',node, &
              'limits= ',txl,txr,tyl,tyr,tzl,tzr,' Ntotal in node=',ntot


    ! First time through, use the boundaries of the box
    if (totparticles.eq.0) then
      xl = txl + dbuffer
      yl = tyl + dbuffer
      zl = tzl + dbuffer
      xr = txr - dbuffer
      yr = tyr - dbuffer
      zr = tzr - dbuffer
    else

      ! Double check the boundaries with every read in
      if ( txl < xl- dbuffer ) then
      xl = txl + dbuffer
      endif
      if ( tyl < yl- dbuffer ) then
      yl = tyl + dbuffer
      endif
      if ( tzl < zl- dbuffer ) then
      zl = tzl + dbuffer
      endif

      if ( txr > xr+ dbuffer ) then
      xr = txr - dbuffer
      endif
      if ( tyr > yr+ dbuffer ) then
      yr = tyr - dbuffer
      endif
      if ( tzr > zr+ dbuffer ) then
      zr = tzr - dbuffer
      endif

    endif


     ! Buffer size from file, need to use the one for each node
     db  = dbuffer
     xxl = xl+db ; xxr = xr-db
     yyl = yl+db ; yyr = yr-db
     zzl = zl+db ; zzr = zr-db


     do
        read(1,iostat=ierr)nrecord         ! number of particles in this record
        if(ierr/=0)exit

        nin = nin + nrecord
        read(1)     (xpp(i),ypp(i),zpp(i),i=1,nrecord)
        read(1)     (wxp(i),wyp(i),wzp(i),i=1,nrecord)

          do i=1,nrecord

             inside = (xxl <= xpp(i)) .and. (xpp(i)<xxr) .and. ( 1e-8 < xpp(i) ) .and. &
                      (yyl <= ypp(i)) .and. (ypp(i)<yyr) .and. ( 1e-8 < ypp(i) ) .and. &
                      (zzl <= zpp(i)) .and. (zpp(i)<zzr) .and. ( 1e-8 < zpp(i) )

!if(inside.and.(rand()<0.001))then
                   icount = icount +1
                   write(2,'(3f12.6)') xpp(i),ypp(i),zpp(i)
!endif ! inside

          enddo
     enddo
     close(1)

     if(nin /= ntot)stop ' -- Wrong number of particles --'

  endif

  deallocate( xpp, ypp, zpp )
  deallocate( wxp, wyp, wzp )


  readfile = icount
  return

end function readfile
