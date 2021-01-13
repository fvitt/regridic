program main
  use netcdf
  use esmf
  use datetime_mod

  implicit none

  integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: r4 = selected_real_kind( 6) ! 4 byte real

  character(len=512) :: infile='NONE', outfile='NONE', imeshfile='NONE', omeshfile='NONE'
  integer :: ncid, out_ncid

  integer :: varid, ndims, nvars
  integer :: vid, did, len, id, len1, len2
  integer, allocatable :: dimids(:)
  character(len=64) :: name
  integer :: coldimid, timedimid, levdimid, ilevdimid

  type(ESMF_Mesh) :: imesh, omesh
  integer :: rc
  integer                 :: spatialDim
  integer                 :: numOwnedElements
  real(ESMF_KIND_R8), allocatable :: ownedElemCoords(:)
  real(r8), allocatable :: imeshlats(:), imeshlons(:)
  real(r8), allocatable :: omeshlats(:), omeshlons(:)
  integer :: ncols, n, ncolso, xtype, natts, t, ntimes, nlevs, nilevs, k
  real(r8), allocatable :: ilats(:), ilons(:)

  type(ESMF_Field) :: ifield, ofield
  type(ESMF_ArraySpec) :: arrayspec
  type(ESMF_RouteHandle) :: routeHandle

  integer, parameter :: nchars = 8
  real(r8) :: x
  integer :: xi
  character(len=nchars) :: curdate
  character(len=nchars) :: curtime
  integer, allocatable :: int_arr(:)
  real(r8), allocatable :: dbl_arr(:)
  real(r8), allocatable :: dbl_arr2d(:,:)
  real(r4), allocatable :: input_arr(:)
  real(r4), allocatable :: output_arr(:)

  character(len=*), parameter :: nml_file = 'options_nml'
  integer :: unitn

  namelist /options/ infile, outfile, imeshfile, omeshfile

  write(*,*) 'START REGRID ...'

  open(newunit=unitn, file=trim(nml_file), status='old')
  read(unit=unitn, nml=options)
  close(unitn)

  call ESMF_Initialize()

  write(*,*) 'infile: '//trim(infile)
  call check( nf90_open(infile, nf90_nowrite, ncid) )

  call check( nf90_inquire(ncid, nDimensions=ndims, nVariables=nvars, unlimitedDimid=timedimid) )
  allocate( dimids(ndims) )
  call check( nf90_inq_dimid(ncid, 'ncol', coldimid) )
  call check( nf90_inq_dimid(ncid, 'lev',  levdimid) )
  call check( nf90_inq_dimid(ncid, 'ilev', ilevdimid) )

  call check( nf90_inquire_dimension(ncid, timedimid, len=ntimes))
  call check( nf90_inquire_dimension(ncid, ilevdimid, len=nilevs))
  call check( nf90_inquire_dimension(ncid, levdimid, len=nlevs))
  call check( nf90_inquire_dimension(ncid, coldimid, len=ncols))

  call check( nf90_create(outfile, NF90_64BIT_OFFSET, out_ncid) )

  allocate( input_arr(ncols) )

  allocate( ilats(ncols), ilons(ncols) )
  call check( nf90_inq_varid(ncid, 'lat', varid) )
  call check( nf90_get_var( ncid, varid, ilats ) )
  call check( nf90_inq_varid(ncid, 'lon', varid) )
  call check( nf90_get_var( ncid, varid, ilons ) )

  imesh = ESMF_MeshCreate(imeshfile, ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  omesh = ESMF_MeshCreate(omeshfile, ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! obtain mesh lats and lons
  call ESMF_MeshGet(imesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !write(*,*) 'spatialDim:',spatialDim
  !write(*,*) 'numOwnedElements:',numOwnedElements
  allocate( imeshlats(numOwnedElements), imeshlons(numOwnedElements) )
  allocate(ownedElemCoords(spatialDim*numOwnedElements))

  call ESMF_MeshGet(imesh, ownedElemCoords=ownedElemCoords)

  do n = 1,numOwnedElements
     imeshlons(n) = ownedElemCoords(2*n-1)
     imeshlats(n) = ownedElemCoords(2*n)
     if (abs(ilons(n)-imeshlons(n))>1.e-6) then
        stop "input lons do not match"
     endif
     if (abs(ilats(n)-imeshlats(n))>1.e-6) then
        stop "input lats do not match"
     endif
  end do
  deallocate(ownedElemCoords)

  ! obtain mesh lats and lons
  call ESMF_MeshGet(omesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !write(*,*) 'spatialDim:',spatialDim
  !write(*,*) 'numOwnedElements:',numOwnedElements
  ncolso = numOwnedElements
  allocate( output_arr(ncolso) )

  allocate( omeshlats(numOwnedElements), omeshlons(numOwnedElements) )
  allocate(ownedElemCoords(spatialDim*numOwnedElements))

  call ESMF_MeshGet(omesh, ownedElemCoords=ownedElemCoords)

  do n = 1,numOwnedElements
     omeshlons(n) = ownedElemCoords(2*n-1)
     omeshlats(n) = ownedElemCoords(2*n)
  end do
  deallocate(ownedElemCoords)

  ! 2D fields
  call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ifield = ESMF_FieldCreate(imesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ofield = ESMF_FieldCreate(omesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldRegridStore( srcField=ifield, dstField=ofield, routeHandle=routeHandle, &
                              polemethod=ESMF_POLEMETHOD_ALLAVG, &
                              regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! define dims
  do did = 1,ndims
     call check( nf90_inquire_dimension( ncid, did, name=name, len=len ) )
     if ( name=='ncol' ) then
        len=ncolso
     elseif ( did==timedimid ) then
        len=nf90_unlimited
     endif
     call check( nf90_def_dim( out_ncid, name, len, id ))
  enddo

  ! define vars and copy attributes
  do vid = 1,nvars
     call check( nf90_inquire_variable(ncid, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid, vid, dimids=dimids(:ndims)) )

     call check( nf90_def_var(out_ncid, name, xtype, dimids(:ndims), id ) )

     do n = 1, natts
        call check( nf90_inq_attname( ncid, vid, n, name ) )
        call check( nf90_copy_att( ncid, vid, name, out_ncid, id) )
     enddo

  enddo
  call datetime(curdate,curtime)
  call getenv('USER',name)

  ! global attributes
  call check( nf90_put_att( out_ncid, nf90_global, 'created_by', &
              trim(name)//' using regridIC ESMF tool '//curdate//' '//curtime) )
  call check( nf90_put_att( out_ncid, nf90_global, 'input_data', trim(infile) ) )
  call check( nf90_put_att( out_ncid, nf90_global, 'src_grid', trim(imeshfile) ) )
  call check( nf90_put_att( out_ncid, nf90_global, 'dst_grid', trim(omeshfile) ) )


  ! end of define phase
  call check( nf90_enddef( out_ncid ) )


  do vid = 1,nvars
     call check( nf90_inquire_variable(ncid, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid, vid, dimids=dimids(:ndims)) )

     if (name == 'lat') then
        call check( nf90_put_var( out_ncid, vid, omeshlats ) )
     elseif (name == 'lon') then
        call check( nf90_put_var( out_ncid, vid, omeshlons ) )
     elseif (name == 'date_written') then
        do t=1,ntimes
           call check( nf90_put_var( out_ncid, vid, curdate, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif (name == 'time_written') then
        do t=1,ntimes
           call check( nf90_put_var( out_ncid, vid, curtime, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif ( .not. any(dimids(:ndims) == coldimid)) then
        write(*,*) 'copy '//trim(name)
        ! simple copy ...
        if (xtype==nf90_char) then
        else if (xtype==nf90_int) then
           if (ndims==0) then
              call check( nf90_get_var( ncid, vid, xi ) )
              call check( nf90_put_var( out_ncid, vid, xi ) )
           elseif (ndims==1) then
              call check( nf90_inquire_dimension( ncid, dimids(1), len=len ) )
              allocate(int_arr(len))
              call check( nf90_get_var( ncid, vid, int_arr ) )
                 call check( nf90_put_var( out_ncid, vid, int_arr ) )
                 deallocate(int_arr)
              endif
           else if (xtype==nf90_float) then
        else if (xtype==nf90_double) then
           if (ndims==0) then
              call check( nf90_get_var( ncid, vid, x ) )
              call check( nf90_put_var( out_ncid, vid, x ) )
           elseif (ndims==1) then
              call check( nf90_inquire_dimension( ncid, dimids(1), len=len ) )
              allocate(dbl_arr(len))
              call check( nf90_get_var( ncid, vid, dbl_arr ) )
              call check( nf90_put_var( out_ncid, vid, dbl_arr ) )
              deallocate(dbl_arr)
           elseif (ndims==2) then
              call check( nf90_inquire_dimension( ncid, dimids(1), len=len1 ) )
              call check( nf90_inquire_dimension( ncid, dimids(2), len=len2 ) )
              allocate(dbl_arr2d(len1,len2))
              call check( nf90_get_var( ncid, vid, dbl_arr2d ) )
              call check( nf90_put_var( out_ncid, vid, dbl_arr2d ) )
              deallocate(dbl_arr2d)
           endif

        else
           write(*,*) "type not recognized for copy for var "//trim(name)
           stop "type not recognized for copy"
        end if

     else
        write(*,*) 'regrid '//trim(name)
        ! regrid column dependent varaibles
        if ( any(dimids(:ndims) == levdimid) ) then
           do t=1,ntimes
              do k=1,nlevs

                 call check( nf90_get_var( ncid, vid, input_arr, start=(/1,k,t/), count=(/ncols,1,1/) ) )
                 call regrid( input_arr, output_arr )
                 call check( nf90_put_var( out_ncid, vid, output_arr, start=(/1,k,t/), count=(/ncolso,1,1/) ) )

              end do
           end do
        elseif ( any(dimids(:ndims) == ilevdimid) ) then
           do t=1,ntimes
              do k=1,nilevs

                 call check( nf90_get_var( ncid, vid, input_arr, start=(/1,k,t/), count=(/ncols,1,1/) ) )
                 call regrid( input_arr, output_arr )
                 call check( nf90_put_var( out_ncid, vid, output_arr, start=(/1,k,t/), count=(/ncolso,1,1/) ) )

              end do
           end do
        else
           do t=1,ntimes

              call check( nf90_get_var( ncid, vid, input_arr, start=(/1,t/), count=(/ncols,1/) ) )
              call regrid( input_arr, output_arr )
              call check( nf90_put_var( out_ncid, vid, output_arr, start=(/1,t/), count=(/ncolso,1/) ) )

           end do
        endif
     endif


  end do

  call check( nf90_close( ncid ) )
  call check( nf90_close( out_ncid ) )


  call ESMF_FieldRegridRelease(routeHandle, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldDestroy( ifield, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldDestroy( ofield, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_Finalize()
  deallocate(dimids)
  deallocate( imeshlats, imeshlons )
  deallocate( omeshlats, omeshlons )
  deallocate( ilats, ilons )

  write(*,*) 'outfile: '//trim(outfile)

  write(*,*) 'END REGRID.'

contains

  subroutine regrid( indata, outdata )

    real(r4), intent(in) :: indata(:)
    real(r4), intent(out) :: outdata(:)

    real(ESMF_KIND_R8), pointer :: fptr1d(:)
    integer :: i
    integer :: lbnd(1), ubnd(1)

    call ESMF_FieldGet( ifield, localDe=0, farrayPtr=fptr1d, &
         computationalLBound=lbnd, computationalUBound=ubnd, rc=rc )
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    fptr1d = -huge(1._r8)
    do i = lbnd(1), ubnd(1)
       fptr1d(i) = indata(i)
    end do

    call ESMF_FieldRegrid( ifield, ofield, routeHandle, rc=rc )
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_FieldGet(ofield, localDe=0, farrayPtr=fptr1d, rc=rc)
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    outdata(:) = fptr1d(:)

  end subroutine regrid

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program main
