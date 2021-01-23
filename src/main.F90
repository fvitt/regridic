program main
  use netcdf
  use esmf
  use datetime_mod

  implicit none

  integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: r4 = selected_real_kind( 6) ! 4 byte real

  character(len=512) :: infile='NONE', outfile='NONE', imeshfile='NONE', imeshfile_d='NONE', omeshfile='NONE'
  integer :: ncid_in, ncid_out

  integer :: varid, ndims, nvars, nvars_out
  integer :: vid, vid_in, did, len, id, len1, len2
  integer, allocatable :: dimids(:)
  integer, allocatable :: dimids_out(:)
  character(len=64) :: name, dimname
  integer :: timedimid, levdimid, ilevdimid
  integer :: did_ncol, did_ncol_d
  integer :: did_ncol_out

  type(ESMF_Mesh) :: imesh, imesh_d, omesh
  integer :: rc

  real(r8), pointer :: omeshlats(:), omeshlons(:)
  real(r8), pointer :: omeshlats_d(:), omeshlons_d(:)
  integer :: ncols, ncols_d, n, ncolso, xtype, natts, t, ntimes, nlevs, nilevs, k
  integer :: num_cols_in
  real(r8), allocatable :: ilats(:), ilons(:)
  real(r8), allocatable :: ilats_d(:), ilons_d(:)

  type(ESMF_ArraySpec) :: arrayspec

  type(ESMF_Field), target :: ifield, ofield
  type(ESMF_Field), target :: ifield_d
  type(ESMF_Field), pointer :: ifld_ptr
  type(ESMF_RouteHandle), target :: routeHandle
  type(ESMF_RouteHandle), target :: routeHandle_d
  type(ESMF_RouteHandle), pointer :: rh_ptr

  integer, parameter :: nchars = 8
  real(r8) :: x
  integer :: xi
  character(len=nchars) :: curdate
  character(len=nchars) :: curtime
  integer, allocatable :: int_arr(:)
  real(r8), allocatable :: dbl_arr(:)
  real(r8), allocatable :: dbl_arr2d(:,:)
  real(r4), pointer :: input_arr(:)
  real(r4), pointer :: output_arr(:)
  real(r4), pointer :: input_arr_d(:)
  real(r4), pointer :: output_arr_d(:)
  real(r4), pointer :: input_arr_ptr(:)

  character(len=*), parameter :: nml_file = 'options_nml'
  integer :: unitn
  logical :: has_ncol_in, has_ncol_d_in

  did_ncol = -huge(1)
  did_ncol_d = -huge(1)

  namelist /options/ infile, outfile, imeshfile, imeshfile_d, omeshfile

  write(*,*) 'START REGRID ...'

  open(newunit=unitn, file=trim(nml_file), status='old')
  read(unit=unitn, nml=options)
  close(unitn)

  call ESMF_Initialize()

  write(*,*) 'infile: '//trim(infile)
  call check( nf90_open(infile, nf90_nowrite, ncid_in) )

  rc = nf90_inq_dimid(ncid_in, 'ncol',   did_ncol)
  has_ncol_in = rc == nf90_noerr

  rc = nf90_inq_dimid(ncid_in, 'ncol_d', did_ncol_d)
  has_ncol_d_in = rc == nf90_noerr

  if ((.not.has_ncol_in).and.(.not.has_ncol_d_in)) then
     write(*,*) 'input file must have ncol or ncol_d dimension'
     stop "ERROR: invalid input data file"
  endif

  call check( nf90_inquire(ncid_in, nDimensions=ndims, nVariables=nvars, unlimitedDimid=timedimid) )
  allocate( dimids(ndims) )
  call check( nf90_inq_dimid(ncid_in, 'lev',  levdimid) )
  call check( nf90_inq_dimid(ncid_in, 'ilev', ilevdimid) )

  call check( nf90_inquire_dimension(ncid_in, timedimid, len=ntimes))
  call check( nf90_inquire_dimension(ncid_in, ilevdimid, len=nilevs))
  call check( nf90_inquire_dimension(ncid_in, levdimid, len=nlevs))

  call check( nf90_create(outfile, NF90_64BIT_OFFSET, ncid_out) )

  if (has_ncol_in) then
     call check( nf90_inquire_dimension(ncid_in, did_ncol, len=ncols))
     allocate( ilats(ncols), ilons(ncols) )
     call check( nf90_inq_varid(ncid_in, 'lat', varid) )
     call check( nf90_get_var( ncid_in, varid, ilats ) )
     call check( nf90_inq_varid(ncid_in, 'lon', varid) )
     call check( nf90_get_var( ncid_in, varid, ilons ) )
     call create_mesh( imeshfile, imesh, 'NCOL', lons_in=ilons, lats_in=ilats )
  end if

  if (has_ncol_d_in) then
     call check( nf90_inquire_dimension(ncid_in, did_ncol_d, len=ncols_d))
     allocate( ilats_d(ncols_d), ilons_d(ncols_d) )
     call check( nf90_inq_varid(ncid_in, 'lat_d', varid) )
     call check( nf90_get_var( ncid_in, varid, ilats_d ) )
     call check( nf90_inq_varid(ncid_in, 'lon_d', varid) )
     call check( nf90_get_var( ncid_in, varid, ilons_d ) )
     call create_mesh( imeshfile_d, imesh_d, 'NCOL_D', lons_in=ilons_d, lats_in=ilats_d )
  end if

  call create_mesh( omeshfile, omesh, 'OUTPUT', ncols_out=ncolso,  lons_out=omeshlons,lats_out=omeshlats )

  allocate( input_arr(ncols) )
  if (has_ncol_d_in) then
     allocate( input_arr_d(ncols_d) )
  endif

  allocate( output_arr(ncolso) )

  ! 2D fields
  call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (has_ncol_in) then
     ifield = ESMF_FieldCreate(imesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  if (has_ncol_d_in) then
     ifield_d = ESMF_FieldCreate(imesh_d, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  ofield = ESMF_FieldCreate(omesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (has_ncol_in) then
     call ESMF_FieldRegridStore( srcField=ifield, dstField=ofield, routeHandle=routeHandle, &
          polemethod=ESMF_POLEMETHOD_ALLAVG, &
          regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  if (has_ncol_d_in) then
     call ESMF_FieldRegridStore( srcField=ifield_d, dstField=ofield, routeHandle=routeHandle_d, &
          polemethod=ESMF_POLEMETHOD_ALLAVG, &
          regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  allocate(dimids_out(ndims))

  ! define dims
  do did = 1,ndims
     call check( nf90_inquire_dimension( ncid_in, did, name=name, len=len ) )

     if ( name=='ncol_d' .and. has_ncol_in ) then
        cycle
     else if( name=='ncol_d' .and. .not.has_ncol_in ) then
        name = 'ncol'
     end if

     if ( name=='ncol' ) then
        len=ncolso
     elseif ( did==timedimid ) then
        len=nf90_unlimited
     endif

!     ndims_out = ndims_out + 1
     call check( nf90_def_dim( ncid_out, name, len, id ) )
  enddo

  ! define vars and copy attributes
  do vid = 1,nvars
     call check( nf90_inquire_variable(ncid_in, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid_in, vid, dimids=dimids(:ndims) ))

     if (name == 'area_d') then
        cycle
     endif
     if (has_ncol_in) then
        if (name == 'lat_d' .or. name == 'lon_d') then
           cycle
        endif
     else
        if (name == 'lat_d') then
           name = 'lat'
        endif
        if (name == 'lon_d') then
           name = 'lon'
        endif
     endif

     do did = 1,ndims
        call check( nf90_inquire_dimension(ncid_in, dimids(did), dimname ) )
        if (dimname == 'ncol_d') dimname = 'ncol'
        call check( nf90_inq_dimid(ncid_out, dimname, dimids_out(did)) )
     end do

     if ( xtype==nf90_double .and. any(dimids(:ndims)==timedimid) .and. &
         (any(dimids(:ndims)==did_ncol).or.any(dimids(:ndims)==did_ncol_d)) ) then
        xtype = nf90_float
     endif

     call check( nf90_def_var(ncid_out, name, xtype, dimids_out(:ndims), id ) )

     do n = 1, natts
        call check( nf90_inq_attname( ncid_in, vid, n, name ) )
        call check( nf90_copy_att( ncid_in, vid, trim(name), ncid_out, id) )
     enddo

  enddo

  call datetime(curdate,curtime)
  call getenv('USER',name)

  ! global attributes
  call check( nf90_put_att( ncid_out, nf90_global, 'created_by', &
              trim(name)//' using regridIC ESMF tool '//curdate//' '//curtime) )
  call check( nf90_put_att( ncid_out, nf90_global, 'input_data', trim(infile) ) )
  call check( nf90_put_att( ncid_out, nf90_global, 'src_grid', trim(imeshfile) ) )
  if (has_ncol_d_in) then
     call check( nf90_put_att( ncid_out, nf90_global, 'src_d_grid', trim(imeshfile_d) ) )
  endif
  call check( nf90_put_att( ncid_out, nf90_global, 'dst_grid', trim(omeshfile) ) )

  ! end of define phase
  call check( nf90_enddef( ncid_out ) )
  call check( nf90_inq_dimid(ncid_out, 'ncol', did_ncol_out) )

  call check( nf90_inquire(ncid_out, nVariables=nvars_out) )

  do vid = 1,nvars_out
     call check( nf90_inquire_variable(ncid_out, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid_out, vid, dimids=dimids_out(:ndims)) )

     if (name == 'lat') then
        call check( nf90_put_var( ncid_out, vid, omeshlats ) )
     elseif (name == 'lon') then
        call check( nf90_put_var( ncid_out, vid, omeshlons ) )
     elseif (name == 'date_written') then
        do t=1,ntimes
           call check( nf90_put_var( ncid_out, vid, curdate, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif (name == 'time_written') then
        do t=1,ntimes
           call check( nf90_put_var( ncid_out, vid, curtime, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif ( .not. any(dimids_out(:ndims) == did_ncol_out)) then
        write(*,*) 'copy '//trim(name)
        call check( nf90_inq_varid(ncid_in, name, vid_in) )
        call check( nf90_inquire_variable(ncid_in, vid_in, dimids=dimids(:ndims)) )
        ! simple copy ...
        if (xtype==nf90_char) then
        else if (xtype==nf90_int) then
           if (ndims==0) then
              call check( nf90_get_var( ncid_in, vid_in, xi ) )
              call check( nf90_put_var( ncid_out, vid, xi ) )
           elseif (ndims==1) then
              call check( nf90_inquire_dimension( ncid_in, dimids(1), len=len ) )
              allocate(int_arr(len))
              call check( nf90_get_var( ncid_in, vid_in, int_arr ) )
              call check( nf90_put_var( ncid_out, vid, int_arr ) )
              deallocate(int_arr)
           endif
        else if (xtype==nf90_float) then
        else if (xtype==nf90_double) then
           if (ndims==0) then
              call check( nf90_get_var( ncid_in, vid_in, x ) )
              call check( nf90_put_var( ncid_out, vid, x ) )
           elseif (ndims==1) then
              call check( nf90_inquire_dimension( ncid_in, dimids(1), len=len ) )
              allocate(dbl_arr(len))
              call check( nf90_get_var( ncid_in, vid_in, dbl_arr ) )
              call check( nf90_put_var( ncid_out, vid, dbl_arr ) )
              deallocate(dbl_arr)
           elseif (ndims==2) then
              call check( nf90_inquire_dimension( ncid_in, dimids(1), len=len1 ) )
              call check( nf90_inquire_dimension( ncid_in, dimids(2), len=len2 ) )
              allocate(dbl_arr2d(len1,len2))
              call check( nf90_get_var( ncid_in, vid_in, dbl_arr2d ) )
              call check( nf90_put_var( ncid_out, vid, dbl_arr2d ) )
              deallocate(dbl_arr2d)
           endif

        else
           write(*,*) "type not recognized for copy for var "//trim(name)
           stop "type not recognized for copy"
        end if

     else
        call check( nf90_inq_varid(ncid_in, name, vid_in) )
        call check( nf90_inquire_variable(ncid_in, vid_in, dimids=dimids(:ndims)) )
        write(*,*) 'regrid '//trim(name)

        if ( any(dimids(:ndims) == did_ncol_d) ) then
           rh_ptr => routehandle_d
           ifld_ptr => ifield_d
           input_arr_ptr => input_arr_d
           num_cols_in = ncols_d
        else
           rh_ptr => routehandle
           ifld_ptr => ifield
           input_arr_ptr => input_arr
           num_cols_in = ncols
        endif

        ! regrid column dependent varaibles
        if ( any(dimids(:ndims) == levdimid) ) then
           do t=1,ntimes
              do k=1,nlevs

                 call check( nf90_get_var( ncid_in, vid_in, input_arr_ptr, start=(/1,k,t/), count=(/num_cols_in,1,1/) ) )
                 call regrid( rh_ptr, ifld_ptr, ofield, input_arr_ptr, output_arr )
                 call check( nf90_put_var( ncid_out, vid, output_arr, start=(/1,k,t/), count=(/ncolso,1,1/) ) )

              end do
           end do
        elseif ( any(dimids(:ndims) == ilevdimid) ) then
           do t=1,ntimes
              do k=1,nilevs

                 call check( nf90_get_var( ncid_in, vid_in, input_arr_ptr, start=(/1,k,t/), count=(/num_cols_in,1,1/) ) )
                 call regrid( rh_ptr, ifld_ptr, ofield, input_arr_ptr, output_arr )
                 call check( nf90_put_var( ncid_out, vid, output_arr, start=(/1,k,t/), count=(/ncolso,1,1/) ) )

              end do
           end do
        else
           do t=1,ntimes

              call check( nf90_get_var( ncid_in, vid_in, input_arr_ptr, start=(/1,t/), count=(/num_cols_in,1/) ) )
              call regrid( rh_ptr, ifld_ptr, ofield, input_arr_ptr, output_arr )
              call check( nf90_put_var( ncid_out, vid, output_arr, start=(/1,t/), count=(/ncolso,1/) ) )

           end do
        endif
     endif
     nullify(rh_ptr)
     nullify(ifld_ptr)
     nullify(input_arr_ptr)

  end do

  call check( nf90_close( ncid_in ) )
  call check( nf90_close( ncid_out ) )

  nullify(rh_ptr)

  if (has_ncol_in) then
     call ESMF_FieldRegridRelease(routeHandle, rc=rc)
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
     call ESMF_FieldDestroy( ifield, rc=rc )
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif
  if (has_ncol_d_in) then
     call ESMF_FieldRegridRelease(routeHandle_d, rc=rc)
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
     call ESMF_FieldDestroy( ifield_d, rc=rc )
     if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  call ESMF_FieldDestroy( ofield, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_Finalize()
  deallocate(dimids)
  deallocate(dimids_out)
  deallocate( omeshlats, omeshlons )
  deallocate( input_arr )
  deallocate( output_arr )
  if (associated(input_arr_d)) nullify(input_arr_d)
  if (allocated(ilats)) deallocate( ilats )
  if (allocated(ilons)) deallocate( ilons )
  if (allocated(ilats_d)) deallocate( ilats_d )
  if (allocated(ilons_d)) deallocate( ilons_d )

  write(*,*) 'outfile: '//trim(outfile)

  write(*,*) 'END REGRID.'

contains

  subroutine create_mesh( meshfile, esmfmesh, errstr, lons_in, lats_in, ncols_out, lons_out, lats_out)
    character(len=*) :: meshfile
    type(ESMF_Mesh), intent(out) :: esmfmesh
    character(len=*) :: errstr
    real(r8), optional, intent(in) :: lons_in(:)
    real(r8), optional, intent(in) :: lats_in(:)
    integer, optional, intent(out) :: ncols_out
    real(r8), optional, pointer :: lons_out(:)
    real(r8), optional, pointer :: lats_out(:)

    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(ESMF_KIND_R8), allocatable :: ownedElemCoords(:)

    real(r8), pointer :: meshlats(:), meshlons(:)

    esmfmesh = ESMF_MeshCreate(meshfile, ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! obtain mesh lats and lons
    call ESMF_MeshGet(esmfmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    allocate(meshlats(numOwnedElements), meshlons(numOwnedElements))

    allocate(ownedElemCoords(spatialDim*numOwnedElements))

    call ESMF_MeshGet(esmfmesh, ownedElemCoords=ownedElemCoords)

    do n = 1,numOwnedElements
       meshlons(n) = ownedElemCoords(2*n-1)
       meshlats(n) = ownedElemCoords(2*n)
       if (present(lons_in)) then
          if (abs(lons_in(n)-meshlons(n))>1.e-6_r8) then
          if (abs(360._r8-abs(lons_in(n)-meshlons(n)))>1.e-6_r8) then
             write(*,*) trim(errstr)//' : meshfile: '//trim(meshfile)
             write(*,*) trim(errstr)//' : lons_in(n), meshlons(n) : ',lons_in(n), meshlons(n)
             stop "create_mesh: lons do not match -- "
          endif
          endif
       endif
       if (present(lats_in)) then
          if (abs(lats_in(n)-meshlats(n))>1.e-6_r8) then
             write(*,*) trim(errstr)//' : meshfile: '//trim(meshfile)
             write(*,*) trim(errstr)//' : lat_in(n), meshlats(n) : ',lats_in(n), meshlats(n)
             stop "create_mesh: lats do not match"
          endif
       endif
    end do

    if (present(lons_out)) then
       lons_out => meshlons
    else
       deallocate(meshlons)
    endif

    if (present(lats_out)) then
       lats_out => meshlats
    else
       deallocate(meshlats)
    endif

    if (present(ncols_out)) then
       ncols_out = numOwnedElements
    endif

    deallocate(ownedElemCoords)

  end subroutine create_mesh


  subroutine regrid( myroutehandle, infield, outfield, indata, outdata )

    type(ESMF_RouteHandle), intent(inout) :: myrouteHandle
    type(ESMF_Field), intent(inout) :: infield
    type(ESMF_Field), intent(inout) :: outfield
    real(r4), intent(in) :: indata(:)
    real(r4), intent(out) :: outdata(:)

    real(ESMF_KIND_R8), pointer :: fptr1d(:)
    integer :: i
    integer :: lbnd(1), ubnd(1)

    call ESMF_FieldGet( infield, localDe=0, farrayPtr=fptr1d, &
         computationalLBound=lbnd, computationalUBound=ubnd, rc=rc )
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    fptr1d = -huge(1._r8)
    do i = lbnd(1), ubnd(1)
       fptr1d(i) = indata(i)
    end do

    call ESMF_FieldRegrid( infield, outfield, myrouteHandle, rc=rc )
    if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_FieldGet(outfield, localDe=0, farrayPtr=fptr1d, rc=rc)
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
