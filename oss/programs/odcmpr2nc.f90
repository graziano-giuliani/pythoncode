program odcmpr2nc
  use netcdf
  implicit none
  integer :: narg
  integer , parameter :: stderr = 0
  integer , parameter :: stdout = 6
  character(len=256) :: basef
  character(len=256) :: infile1 , infile2
  character(len=256) :: outfile
  integer :: ierr , iomode
  integer :: ncid , ncstat , ichdimid , icwvnvar , imolidfixvar , imolidvar
  integer :: iprefvar , itmptabvar , iwvptabvar , icoefvar
  integer :: iichmapvar , inchvar , iselsvar , ivwvnvar , imolsvar
  integer :: ikfixvar , idkh2ovar , ikh2ovar , ikvarvar
  integer :: imolidfixdimid , imoliddimid , ilevdimid , intmpoddimid
  integer :: inlayoddimid , inmoltabdimid , infdimid , inchmaxdimid
  integer :: imxmolsdimid
  integer , parameter :: iunit1 = 101, iunit2 = 102
  integer :: uid_sel , uid_od
  integer :: nchan , nf_sel , nchmax , nfsmp , nmolx , nn
  integer :: nlayod , ntmpod , nlev , ismp , i , j , k , n
  integer :: mxmols
  integer(2) :: nmolfix , nmol , nmols
  integer(2) , allocatable , dimension(:) :: molidfix , molid
  real , dimension(:) , allocatable :: cwvn
  real , dimension(:) , allocatable :: pref
  real(8) , dimension(:) , allocatable :: vwvn
  real , dimension(:,:) , allocatable :: tmptab , wvptab , coef
  real , dimension(:,:,:) , allocatable :: kfix , dkh2o , kh2o
  real , dimension(:,:,:,:) , allocatable :: kvar
  integer , dimension(:,:) , allocatable :: ichmap
  integer , dimension(:) , allocatable :: isels
  integer(2) , dimension(:) , allocatable :: nch
  integer(2) , dimension(:,:) , allocatable :: imols
  character(len=100) :: instr_info
  integer , dimension(4) :: dimids
  real(8) :: v1 , v2

  call purpose

  narg = command_argument_count( )
  if ( narg < 1 ) then
    write(stderr,*) 'Not enough arguments. Need basis file name'
    call usage( )
  end if

  call get_command_argument(1,basef)
  write(infile1,'(a,a)') trim(basef), '.od.cmpr'
  write(infile2,'(a,a)') trim(basef), '.sel.cmpr'
  write(outfile,'(a,a)') trim(basef), '.nc'

  write(stdout,*) 'Opening input file OD  : ',trim(infile1)
  write(stdout,*) 'Opening input file SEL : ',trim(infile2)

  open(iunit1,file=infile1,status='old',form='unformatted', &
       action='read',convert='BIG_ENDIAN')
  open(iunit2,file=infile2,status='old',form='unformatted', &
       action='read',convert='BIG_ENDIAN')

  read(iunit2,iostat=ierr) uid_sel
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading uid_sel from SEL file'
    stop
  end if
  read(iunit1,iostat=ierr) uid_od
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading uid_od from OD file'
    stop
  end if
  read(iunit2,iostat=ierr) ! Skip a record, in file this is integer 0 
  read(iunit2,iostat=ierr) instr_info
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading instr_info from SEL file'
    stop
  end if
  read(iunit1,iostat=ierr) v1 , v2 , nfsmp
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading v1 v2 nfsmp from SEL file'
    stop
  end if
  read(iunit2,iostat=ierr) nchan , nf_sel , nchmax
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading nchan nf_sel nchmax from SEL file'
    stop
  end if
  if ( nf_sel /= nfsmp ) then
    write(stderr,*) 'ERROR: the OD and SEL file do not match!'
    write(stderr,*) 'SEL nf_sel   : ',nf_sel
    write(stderr,*) 'OD  nfsmp : ',nfsmp
    stop
  end if
  allocate(cWvn(nchan),nch(nf_sel),coef(nchmax,nf_sel),   &
           isels(nf_sel),ichmap(nchmax,nf_sel),vwvn(nf_sel), &
           stat=ierr)
  coef = -9999.0
  ichmap = -1
  if ( ierr /= 0 ) then
    write(stderr,*) 'Allocation error!'
    stop
  end if
  read(iunit2,iostat=ierr) cwvn
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading cwvn from SEL file'
    stop
  end if
  read(iunit1,iostat=ierr) nmolfix , nmol
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading nmolfix nmol from OD file'
    stop
  end if
  mxmols = nmolfix+nmol
  allocate(molidfix(nmolfix),molid(nmol),stat=ierr)
  if ( ierr /= 0 ) then
    write(stderr,*) 'Allocation error!'
    stop
  end if
  read(iunit1,iostat=ierr) molidfix , molid
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading molidfix molid from OD file'
    stop
  end if
  read(iunit1,iostat=ierr) nlayod , ntmpod
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading nlayod ntmpod from OD file'
    stop
  end if
  nlev = nlayod+1
  allocate(pref(nlev),tmptab(ntmpod,nlayod),wvptab(nmol+2,nlayod), &
           imols(mxmols,nf_sel),kfix(nlayod,ntmpod,nf_sel),     &
           dkh2o(nlayod,ntmpod,nf_sel),                         &
           kh2o(nlayod,ntmpod,nf_sel),                          &
           kvar(mxmols,nlayod,ntmpod,nf_sel),  stat=ierr)
  if ( ierr /= 0 ) then
    write(stderr,*) 'Allocation error!'
    stop
  end if
  imols = -1
  kfix = -9999.0
  dkh2o = -9999.0
  kh2o = -9999.0
  kvar = -9999.0
  read(iunit1,iostat=ierr) pref , tmptab , wvptab
  if ( ierr /= 0 ) then
    write(stderr,*) 'Error reading pref tmptab wvptab from OD file'
    stop
  end if
  do ismp = 1 , nf_sel
    read(iunit2,iostat=ierr) isels(ismp) , nch(ismp)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Error reading isels nch from SEL file at record ',ismp
      stop
    end if
    read(iunit2,iostat=ierr) coef(1:nch(ismp),ismp), &
                             ichmap(1:nch(ismp),ismp)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Error reading coef ichmap from SEL '// &
                      'file at record ',ismp
      stop
    end if
  end do
  do ismp = 1 , nf_sel
    read(iunit1,iostat=ierr) vwvn(ismp) , nmols
    if ( ierr /= 0 ) then
      write(stderr,*) 'Error reading vwvn nmols from OD '// &
                      'file at record ',ismp
      stop
    end if
    read(iunit1,iostat=ierr) imols(1:nmols,ismp)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Error reading imols from OD file at record ',ismp
      stop
    end if
    read(iunit1,iostat=ierr) kfix(:,:,ismp),dkh2o(:,:,ismp),&
                             kh2o(:,:,ismp),kvar(1:nmols-1,:,:,ismp)
    if ( ierr /= 0 ) then
      write(stderr,*) 'Error reading kfix dkh2o kh2o kvar from '//&
                      'OD file at record ',ismp
      stop
    end if
  end do

  write(stdout,*) 'Opening output file    : ',trim(outfile)
  iomode = ior(ior(nf90_classic_model,nf90_clobber),nf90_netcdf4)
  ncstat = nf90_create(outfile,iomode,ncid)
  call check_ncerror('Cannot create output file '//trim(outfile))

  ncstat = nf90_put_att(ncid,nf90_global,'id_sel',uid_sel)
  call check_ncerror('Cannot add global attribute id_sel')
  ncstat = nf90_put_att(ncid,nf90_global,'id_od',uid_od)
  call check_ncerror('Cannot add global attribute id_od')
  ncstat = nf90_put_att(ncid,nf90_global,'info',trim(instr_info))
  call check_ncerror('Cannot add global attribute info')
  ncstat = nf90_put_att(ncid,nf90_global,'v1',v1)
  call check_ncerror('Cannot add global attribute v1')
  ncstat = nf90_put_att(ncid,nf90_global,'v2',v2)
  call check_ncerror('Cannot add global attribute v2')

  ncstat = nf90_def_dim(ncid,'nchan',nchan,ichdimid)
  call check_ncerror('Cannot create nchan dimension')
  ncstat = nf90_def_dim(ncid,'nf',nf_sel,infdimid)
  call check_ncerror('Cannot create nf dimension')
  ncstat = nf90_def_dim(ncid,'nchmax',nchmax,inchmaxdimid)
  call check_ncerror('Cannot create nchmax dimension')
  nmolx = nmolfix
  ncstat = nf90_def_dim(ncid,'nmolfix',nmolx,imolidfixdimid)
  call check_ncerror('Cannot create nmolfix dimension')
  nmolx = nmol
  ncstat = nf90_def_dim(ncid,'nmol',nmolx,imoliddimid)
  call check_ncerror('Cannot create nmol dimension')
  ncstat = nf90_def_dim(ncid,'nlev',nlev,ilevdimid)
  call check_ncerror('Cannot create nlev dimension')
  ncstat = nf90_def_dim(ncid,'ntmpod',ntmpod,intmpoddimid)
  call check_ncerror('Cannot create ntmpod dimension')
  ncstat = nf90_def_dim(ncid,'nlayod',nlayod,inlayoddimid)
  call check_ncerror('Cannot create nlayod dimension')
  ncstat = nf90_def_dim(ncid,'nmoltab',nmol+2,inmoltabdimid)
  call check_ncerror('Cannot create nmoltab dimension')
  ncstat = nf90_def_dim(ncid,'mxmols',mxmols,imxmolsdimid)
  call check_ncerror('Cannot create mxmols dimension')

  dimids(1) = ichdimid
  ncstat = nf90_def_var(ncid,'cWvn',nf90_float,dimids(1:1),icwvnvar)
  call check_ncerror('Cannot create var cWvn')
  dimids(1) = imolidfixdimid
  ncstat = nf90_def_var(ncid,'molidfix',nf90_short,dimids(1:1),imolidfixvar)
  call check_ncerror('Cannot create var molidfix')
  dimids(1) = imoliddimid
  ncstat = nf90_def_var(ncid,'molid',nf90_short,dimids(1:1),imolidvar)
  call check_ncerror('Cannot create var molid')
  dimids(1) = ilevdimid
  ncstat = nf90_def_var(ncid,'pref',nf90_float,dimids(1:1),iprefvar)
  call check_ncerror('Cannot create var pref')
  dimids(1) = intmpoddimid
  dimids(2) = inlayoddimid
  ncstat = nf90_def_var(ncid,'tmptab',nf90_float,dimids(1:2),itmptabvar)
  call check_ncerror('Cannot create var tmptab')
  dimids(1) = inmoltabdimid
  dimids(2) = inlayoddimid
  ncstat = nf90_def_var(ncid,'wvptab',nf90_float,dimids(1:2),iwvptabvar)
  call check_ncerror('Cannot create var wvptab')
  dimids(1) = inchmaxdimid
  dimids(2) = infdimid
  ncstat = nf90_def_var(ncid,'coef',nf90_float,dimids(1:2),icoefvar)
  call check_ncerror('Cannot create var coef')
  ncstat = nf90_put_att(ncid,icoefvar,'_FillValue',-9999.0)
  call check_ncerror('Cannot add attribute _FillValue for coef variable')
  dimids(1) = inchmaxdimid
  dimids(2) = infdimid
  ncstat = nf90_def_var(ncid,'ichmap',nf90_int,dimids(1:2),iichmapvar)
  call check_ncerror('Cannot create var ichmap')
  ncstat = nf90_put_att(ncid,iichmapvar,'_FillValue',-1)
  call check_ncerror('Cannot add attribute _FillValue for ichmap variable')
  dimids(1) = infdimid
  ncstat = nf90_def_var(ncid,'nch',nf90_short,dimids(1:1),inchvar)
  call check_ncerror('Cannot create var nch')
  ncstat = nf90_def_var(ncid,'isels',nf90_int,dimids(1:1),iselsvar)
  call check_ncerror('Cannot create var isels')
  ncstat = nf90_def_var(ncid,'vwvn',nf90_double,dimids(1:1),ivwvnvar)
  call check_ncerror('Cannot create var vwvn')
  dimids(1) = imxmolsdimid
  dimids(2) = infdimid
  ncstat = nf90_def_var(ncid,'imols',nf90_short,dimids(1:2),imolsvar)
  call check_ncerror('Cannot create var imols')
  ncstat = nf90_put_att(ncid,imolsvar,'_FillValue',-1_2)
  call check_ncerror('Cannot add attribute _FillValue for imols variable')
  dimids(1) = inlayoddimid
  dimids(2) = intmpoddimid
  dimids(3) = infdimid
  ncstat = nf90_def_var(ncid,'kfix',nf90_real,dimids(1:3),ikfixvar)
  call check_ncerror('Cannot create var kfix')
  ncstat = nf90_put_att(ncid,ikfixvar,'_FillValue',-9999.0)
  call check_ncerror('Cannot add attribute _FillValue for kfix variable')
  ncstat = nf90_def_var_deflate(ncid,ikfixvar,1,1,9)
  call check_ncerror('Cannot set deflate level to var kfix')
  ncstat = nf90_def_var(ncid,'dkh2o',nf90_real,dimids(1:3),idkh2ovar)
  call check_ncerror('Cannot create var dkh2o')
  ncstat = nf90_put_att(ncid,idkh2ovar,'_FillValue',-9999.0)
  call check_ncerror('Cannot add attribute _FillValue for dkh2o variable')
  ncstat = nf90_def_var_deflate(ncid,idkh2ovar,1,1,9)
  call check_ncerror('Cannot set deflate level to var dkh2o')
  ncstat = nf90_def_var(ncid,'kh2o',nf90_real,dimids(1:3),ikh2ovar)
  call check_ncerror('Cannot create var kh2o')
  ncstat = nf90_put_att(ncid,ikh2ovar,'_FillValue',-9999.0)
  call check_ncerror('Cannot add attribute _FillValue for kh2o variable')
  ncstat = nf90_def_var_deflate(ncid,ikh2ovar,1,1,9)
  call check_ncerror('Cannot set deflate level to var kh2o')
  dimids(1) = imxmolsdimid
  dimids(2) = inlayoddimid
  dimids(3) = intmpoddimid
  dimids(4) = infdimid
  ncstat = nf90_def_var(ncid,'kvar',nf90_real,dimids(1:4),ikvarvar)
  call check_ncerror('Cannot create var kvar')
  ncstat = nf90_put_att(ncid,ikvarvar,'_FillValue',-9999.0)
  call check_ncerror('Cannot add attribute _FillValue for kvar variable')
  ncstat = nf90_def_var_deflate(ncid,ikvarvar,1,1,9)
  call check_ncerror('Cannot set deflate level to var kvar')

  ncstat = nf90_enddef(ncid)
  call check_ncerror('Cannot finalize output file '//trim(outfile))

  ncstat = nf90_put_var(ncid,icwvnvar,cwvn)
  call check_ncerror('Cannot write variable cWvn')
  ncstat = nf90_put_var(ncid,imolidfixvar,molidfix)
  call check_ncerror('Cannot write variable imolidfix')
  ncstat = nf90_put_var(ncid,imolidvar,molid)
  call check_ncerror('Cannot write variable imolid')
  ncstat = nf90_put_var(ncid,iprefvar,pref)
  call check_ncerror('Cannot write variable pref')
  ncstat = nf90_put_var(ncid,itmptabvar,tmptab)
  call check_ncerror('Cannot write variable tmptab')
  ncstat = nf90_put_var(ncid,iwvptabvar,wvptab)
  call check_ncerror('Cannot write variable wvptab')
  ncstat = nf90_put_var(ncid,icoefvar,coef)
  call check_ncerror('Cannot write variable coef')
  ncstat = nf90_put_var(ncid,iichmapvar,ichmap)
  call check_ncerror('Cannot write variable ichmap')
  ncstat = nf90_put_var(ncid,inchvar,nch)
  call check_ncerror('Cannot write variable nch')
  ncstat = nf90_put_var(ncid,iselsvar,isels)
  call check_ncerror('Cannot write variable isels')
  ncstat = nf90_put_var(ncid,ivwvnvar,vwvn)
  call check_ncerror('Cannot write variable vwvn')
  ncstat = nf90_put_var(ncid,imolsvar,imols)
  call check_ncerror('Cannot write variable imols')
  ncstat = nf90_put_var(ncid,ikfixvar,kfix)
  call check_ncerror('Cannot write variable kfix')
  ncstat = nf90_put_var(ncid,idkh2ovar,dkh2o)
  call check_ncerror('Cannot write variable dkh2o')
  ncstat = nf90_put_var(ncid,ikh2ovar,kh2o)
  call check_ncerror('Cannot write variable kh2o')
  ncstat = nf90_put_var(ncid,ikvarvar,kvar)
  call check_ncerror('Cannot write variable kvar')

  ncstat = nf90_close(ncid)
  call check_ncerror('Cannot close output file '//trim(outfile))

  close(iunit1)
  close(iunit2)

  contains

    subroutine purpose( )
      implicit none
      write(stderr,*) 'Converts HITRAN pre-computed Look Up Tables in netCDF'
    end subroutine purpose

    subroutine usage( )
      implicit none
      character(len=256) :: myname
      call get_command_argument(0,myname)
      write(stderr,*) 'Usage: '
      write(stderr,*) '       ',trim(myname),' basis'
      write(stderr,*)
      write(stderr,*) 'Example: ',trim(myname),' leo.airs.0.05'
      write(stderr,*)
      stop
    end subroutine usage

    subroutine check_ncerror(message)
      implicit none
      character(len=*) , intent(in) :: message
      if ( ncstat /= nf90_noerr ) then
        write(stderr,*) message
        write(stderr,*) nf90_strerror(ncstat)
        stop
      end if
    end subroutine check_ncerror

end program odcmpr2nc
