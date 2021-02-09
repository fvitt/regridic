program main
  use netcdf
  use esmf
  use datetime_mod

  implicit none

  integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: r4 = selected_real_kind( 6) ! 4 byte real

  integer, parameter :: MAX_NPFLDS=100

  character(len=512) :: infile='NONE', outfile='NONE', imeshfile_p='NONE', imeshfile_d='NONE', omeshfile_p='NONE', omeshfile_d='NONE'
  character(len=32) :: phys_flds(MAX_NPFLDS) = ' '
  integer :: ncid_in, ncid_out

  integer :: varid, ndims, nvars, nvars_out
  integer :: vid, vid_in, did, len, id, len1, len2
  integer, allocatable :: dimids_in(:)
  integer, allocatable :: dimids_out(:)
  character(len=64) :: name, dimname
  integer :: timedimid, levdimid, ilevdimid
  integer :: did_ncol_p_in, did_ncol_d_in
  integer :: did_ncol_p_out, did_ncol_d_out

  type(ESMF_Mesh) :: imesh_p, imesh_d, omesh_p, omesh_d
  integer :: rc

  integer :: ncols_p_in, ncols_d_in, ncols_p_out, ncols_d_out
  integer :: xtype, natts, t, n, ntimes, nlevs, nilevs, numlevs, k
  integer :: num_cols_in, num_cols_out
  real(r8), pointer :: ilats_p(:), ilons_p(:)
  real(r8), pointer :: ilats_d(:), ilons_d(:)
  real(r8), pointer :: olats_p(:), olons_p(:)
  real(r8), pointer :: olats_d(:), olons_d(:)

  type(ESMF_ArraySpec) :: arrayspec

  type(ESMF_Field), target :: ifield_p, ofield_p
  type(ESMF_Field), target :: ifield_d, ofield_d
  type(ESMF_Field), pointer :: ifld_ptr
  type(ESMF_Field), pointer :: ofld_ptr
  type(ESMF_RouteHandle), target :: routeHandle_pp
  type(ESMF_RouteHandle), target :: routeHandle_dd
  type(ESMF_RouteHandle), target :: routeHandle_dp
  type(ESMF_RouteHandle), target :: routeHandle_pd
  type(ESMF_RouteHandle), pointer :: rh_ptr

  integer, parameter :: nchars = 8
  real(r8) :: x
  integer :: xi
  character(len=nchars) :: curdate
  character(len=nchars) :: curtime
  integer, allocatable :: int_arr(:)
  real(r8), allocatable :: dbl_arr(:)
  real(r8), allocatable :: dbl_arr2d(:,:)
  real(r4), pointer :: input_arr_p(:)
  real(r4), pointer :: output_arr_p(:)
  real(r4), pointer :: input_arr_d(:)
  real(r4), pointer :: output_arr_d(:)
  real(r4), pointer :: input_arr_ptr(:)
  real(r4), pointer :: output_arr_ptr(:)

  character(len=*), parameter :: nml_file = 'options_nml'
  integer :: unitn
  logical :: has_ncol_p_in, has_ncol_d_in
  logical :: ouput_phys_grid
  integer :: nphysflds = 0

  did_ncol_p_in = -huge(1)
  did_ncol_d_in = -huge(1)
  did_ncol_p_out = -huge(1)
  did_ncol_d_out = -huge(1)

  namelist /options/ infile, outfile, imeshfile_p, imeshfile_d, omeshfile_p, omeshfile_d, phys_flds

  write(*,*) 'START REGRID ...'
  open(newunit=unitn, file=trim(nml_file), status='old')
  read(unit=unitn, nml=options)
  close(unitn)

  call ESMF_Initialize()

  nphysflds = count( len_trim(phys_flds(:))>0 )
  ouput_phys_grid = nphysflds>0

  write(*,*) 'infile: '//trim(infile)
  call check( nf90_open(infile, nf90_nowrite, ncid_in) )

  rc = nf90_inq_dimid(ncid_in, 'ncol', did_ncol_p_in)
  has_ncol_p_in = rc == nf90_noerr

  rc = nf90_inq_dimid(ncid_in, 'ncol_d', did_ncol_d_in)
  has_ncol_d_in = rc == nf90_noerr

  if ((.not.has_ncol_p_in).and.(.not.has_ncol_d_in)) then
     write(*,*) 'input file must have ncol or ncol_d dimension'
     stop "ERROR: invalid input data file"
  endif

  call check( nf90_inquire(ncid_in, nDimensions=ndims, nVariables=nvars, unlimitedDimid=timedimid) )
  call check( nf90_inq_dimid(ncid_in, 'lev',  levdimid) )
  call check( nf90_inq_dimid(ncid_in, 'ilev', ilevdimid) )

  call check( nf90_inquire_dimension(ncid_in, timedimid, len=ntimes))
  call check( nf90_inquire_dimension(ncid_in, ilevdimid, len=nilevs))
  call check( nf90_inquire_dimension(ncid_in, levdimid, len=nlevs))

  call check( nf90_create(outfile, NF90_64BIT_OFFSET, ncid_out) )

  if (has_ncol_p_in) then
     call check( nf90_inquire_dimension(ncid_in, did_ncol_p_in, len=ncols_p_in))
     allocate( ilats_p(ncols_p_in), ilons_p(ncols_p_in) )
     call check( nf90_inq_varid(ncid_in, 'lat', varid) )
     call check( nf90_get_var( ncid_in, varid, ilats_p ) )
     call check( nf90_inq_varid(ncid_in, 'lon', varid) )
     call check( nf90_get_var( ncid_in, varid, ilons_p ) )
     call create_mesh( imeshfile_p, imesh_p, 'INPUT_P', lons_in=ilons_p, lats_in=ilats_p )
  else
     call create_mesh( imeshfile_p, imesh_p, 'INPUT_P', ncols_out=ncols_p_in, lons_out=ilons_p, lats_out=ilats_p )
  end if

  if (has_ncol_d_in) then
     call check( nf90_inquire_dimension(ncid_in, did_ncol_d_in, len=ncols_d_in))
     allocate( ilats_d(ncols_d_in), ilons_d(ncols_d_in) )
     call check( nf90_inq_varid(ncid_in, 'lat_d', varid) )
     call check( nf90_get_var( ncid_in, varid, ilats_d ) )
     call check( nf90_inq_varid(ncid_in, 'lon_d', varid) )
     call check( nf90_get_var( ncid_in, varid, ilons_d ) )
     call create_mesh( imeshfile_d, imesh_d, 'NPUT_D', lons_in=ilons_d, lats_in=ilats_d )
  else
     call create_mesh( imeshfile_d, imesh_d, 'INPUT_D', ncols_out=ncols_d_in, lons_out=ilons_d, lats_out=ilats_d )
  end if

  call create_mesh( omeshfile_p, omesh_p, 'OUTPUT_P', ncols_out=ncols_p_out, lons_out=olons_p, lats_out=olats_p )
  call create_mesh( omeshfile_d, omesh_d, 'OUTPUT_D', ncols_out=ncols_d_out, lons_out=olons_d, lats_out=olats_d )
  allocate( input_arr_p(ncols_p_in) )
  allocate( input_arr_d(ncols_d_in) )

  allocate( output_arr_p(ncols_p_out) )
  allocate( output_arr_d(ncols_d_out) )

  ! 2D fields
  call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ifield_p = ESMF_FieldCreate(imesh_p, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ifield_d = ESMF_FieldCreate(imesh_d, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ofield_p = ESMF_FieldCreate(omesh_p, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ofield_d = ESMF_FieldCreate(omesh_d, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! p->p
  call ESMF_FieldRegridStore( srcField=ifield_p, dstField=ofield_p, routeHandle=routeHandle_pp, &
       polemethod=ESMF_POLEMETHOD_ALLAVG, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! p->d
  call ESMF_FieldRegridStore( srcField=ifield_p, dstField=ofield_d, routeHandle=routeHandle_pd, &
       polemethod=ESMF_POLEMETHOD_ALLAVG, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! d->p
  call ESMF_FieldRegridStore( srcField=ifield_d, dstField=ofield_p, routeHandle=routeHandle_dp, &
       polemethod=ESMF_POLEMETHOD_ALLAVG, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! d->d
  call ESMF_FieldRegridStore( srcField=ifield_d, dstField=ofield_d, routeHandle=routeHandle_dd, &
       polemethod=ESMF_POLEMETHOD_ALLAVG, regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  allocate( dimids_in(ndims) )
  allocate(dimids_out(ndims))

  ! define dims
  do did = 1,ndims
     call check( nf90_inquire_dimension( ncid_in, did, name=name, len=len ) )

     if ( name=='ncol' ) then
        len=ncols_p_out
     elseif ( name=='ncol_d' ) then
        len=ncols_d_out
     elseif ( did==timedimid ) then
        len=nf90_unlimited
     endif

     call check( nf90_def_dim( ncid_out, name, len, id ) )
  enddo

  ! make sure the output file has the phys and dyn column dimensions
  rc = nf90_inq_dimid(ncid_out, 'ncol', did_ncol_p_out)
  if (rc /= nf90_noerr) then
     call check( nf90_def_dim( ncid_out, 'ncol', ncols_p_out, did_ncol_p_out ) )
  end if
  rc = nf90_inq_dimid(ncid_out, 'ncol_d', did_ncol_d_out)
  if (rc /= nf90_noerr) then
     call check( nf90_def_dim( ncid_out, 'ncol_d', ncols_d_out, did_ncol_d_out ) )
  end if

  ! define vars and copy attributes
  do vid = 1,nvars
     call check( nf90_inquire_variable(ncid_in, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid_in, vid, dimids=dimids_in(:ndims) ))

     if (name == 'area_d' .or. name == 'lat' .or. name == 'lon' .or. name == 'lat_d' .or. name == 'lon_d') then
        cycle
     endif

     if ( xtype==nf90_double .and. any(dimids_in(:ndims)==timedimid) .and. &
         (any(dimids_in(:ndims)==did_ncol_p_in).or.any(dimids_in(:ndims)==did_ncol_d_in)) ) then
        xtype = nf90_float
     endif

     do did = 1,ndims
        call check( nf90_inquire_dimension(ncid_in, dimids_in(did), dimname ) )
        call check( nf90_inq_dimid(ncid_out, dimname, dimids_out(did)) )
     end do

     if ( any(dimids_in(:ndims)==timedimid) ) then
        if (any(phys_flds(:)==trim(name))) then
           do did=1,ndims
              if (dimids_out(did)==did_ncol_d_out) dimids_out(did)=did_ncol_p_out
           end do
        else
           do did=1,ndims
              if (dimids_out(did)==did_ncol_p_out) dimids_out(did)=did_ncol_d_out
           end do
        endif
     endif

     call check( nf90_def_var(ncid_out, name, xtype, dimids_out(:ndims), id ) )

     do n = 1, natts
        call check( nf90_inq_attname( ncid_in, vid, n, name ) )
        call check( nf90_copy_att( ncid_in, vid, trim(name), ncid_out, id) )
     enddo

  enddo
  ! make sure the output file has the phys and dyn column coordinate varialbles
  rc = nf90_inq_varid(ncid_out, 'lon', vid)
  if (rc /= nf90_noerr) then
     call check( nf90_inq_dimid(ncid_out, 'ncol', did) )
     call check( nf90_def_var(ncid_out, 'lon', nf90_double, (/did/), id ) )
     call check( nf90_put_att(ncid_out, id, "long_name", "column longitude"))
     call check( nf90_put_att(ncid_out, id, "units", "degrees east"))
  endif
  rc = nf90_inq_varid(ncid_out, 'lat', vid)
  if (rc /= nf90_noerr) then
     call check( nf90_inq_dimid(ncid_out, 'ncol', did) )
     call check( nf90_def_var(ncid_out, 'lat', nf90_double, (/did/), id ) )
     call check( nf90_put_att(ncid_out, id, "long_name", "column latitude"))
     call check( nf90_put_att(ncid_out, id, "units", "degrees north"))
  endif
  rc = nf90_inq_varid(ncid_out, 'lon_d', vid)
  if (rc /= nf90_noerr) then
     call check( nf90_inq_dimid(ncid_out, 'ncol_d', did) )
     call check( nf90_def_var(ncid_out, 'lon_d', nf90_double, (/did/), id ) )
     call check( nf90_put_att(ncid_out, id, "long_name", "column longitude"))
     call check( nf90_put_att(ncid_out, id, "units", "degrees east"))
  endif
  rc = nf90_inq_varid(ncid_out, 'lat_d', vid)
  if (rc /= nf90_noerr) then
     call check( nf90_inq_dimid(ncid_out, 'ncol_d', did) )
     call check( nf90_def_var(ncid_out, 'lat_d', nf90_double, (/did/), id ) )
     call check( nf90_put_att(ncid_out, id, "long_name", "column latitude"))
     call check( nf90_put_att(ncid_out, id, "units", "degrees north"))
  endif

  call datetime(curdate,curtime)
  call getenv('USER',name)

  ! global attributes
  call check( nf90_put_att( ncid_out, nf90_global, 'created_by', &
       trim(name)//' using regridIC ESMF tool '//curdate//' '//curtime) )
  call check( nf90_put_att( ncid_out, nf90_global, 'input_data', trim(infile) ) )
  call check( nf90_put_att( ncid_out, nf90_global, 'src_p_grid', trim(imeshfile_p) ) )
  call check( nf90_put_att( ncid_out, nf90_global, 'src_d_grid', trim(imeshfile_d) ) )
  call check( nf90_put_att( ncid_out, nf90_global, 'dst_p_grid', trim(omeshfile_p) ) )
  call check( nf90_put_att( ncid_out, nf90_global, 'dst_d_grid', trim(omeshfile_d) ) )

  ! end of define phase
  call check( nf90_enddef( ncid_out ) )

  ! set the phys and dyn column coordinate varialbles in the output file
  rc = nf90_inq_varid(ncid_out, 'lon', vid)
  if (rc == nf90_noerr) then
     call check( nf90_put_var( ncid_out, vid, olons_p ) )
  endif
  rc = nf90_inq_varid(ncid_out, 'lat', vid)
  if (rc == nf90_noerr) then
     call check( nf90_put_var( ncid_out, vid, olats_p ) )
  endif
  rc = nf90_inq_varid(ncid_out, 'lon_d', vid)
  if (rc == nf90_noerr) then
     call check( nf90_put_var( ncid_out, vid, olons_d ) )
  endif
  rc = nf90_inq_varid(ncid_out, 'lat_d', vid)
  if (rc == nf90_noerr) then
     call check( nf90_put_var( ncid_out, vid, olats_d ) )
  endif

  call check( nf90_inquire(ncid_out, nVariables=nvars_out) )

  do vid = 1,nvars_out
     call check( nf90_inquire_variable(ncid_out, vid, name=name, xtype=xtype, ndims=ndims, natts=natts ))
     call check( nf90_inquire_variable(ncid_out, vid, dimids=dimids_out(:ndims)) )

     if (name == 'lat' .or. name == 'lon' .or. name == 'lat_d' .or. name == 'lon_d') then
        cycle
     elseif (name == 'date_written') then
        do t=1,ntimes
           call check( nf90_put_var( ncid_out, vid, curdate, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif (name == 'time_written') then
        do t=1,ntimes
           call check( nf90_put_var( ncid_out, vid, curtime, start=(/1,t/), count=(/nchars,1/) ) )
        end do
     elseif (.not.(any(dimids_out(:ndims)==did_ncol_p_out).or.any(dimids_out(:ndims)==did_ncol_d_out))) then

        write(*,*) 'copy '//trim(name)
        call check( nf90_inq_varid(ncid_in, name, vid_in) )
        call check( nf90_inquire_variable(ncid_in, vid_in, dimids=dimids_in(:ndims)) )
        ! simple copy ...
        if (xtype==nf90_char) then
        else if (xtype==nf90_int) then
           if (ndims==0) then
              call check( nf90_get_var( ncid_in, vid_in, xi ) )
              call check( nf90_put_var( ncid_out, vid, xi ) )
           elseif (ndims==1) then
              call check( nf90_inquire_dimension( ncid_in, dimids_in(1), len=len ) )
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
              call check( nf90_inquire_dimension( ncid_in, dimids_in(1), len=len ) )
              allocate(dbl_arr(len))
              call check( nf90_get_var( ncid_in, vid_in, dbl_arr ) )
              call check( nf90_put_var( ncid_out, vid, dbl_arr ) )
              deallocate(dbl_arr)
           elseif (ndims==2) then
              call check( nf90_inquire_dimension( ncid_in, dimids_in(1), len=len1 ) )
              call check( nf90_inquire_dimension( ncid_in, dimids_in(2), len=len2 ) )
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
        call check( nf90_inquire_variable(ncid_in, vid_in, dimids=dimids_in(:ndims)) )
        write(*,*) 'regrid '//trim(name)

        if ( any(dimids_in(:ndims) == did_ncol_p_in) .and. any(dimids_out(:ndims) == did_ncol_p_out) ) then
           ifld_ptr => ifield_p
           ofld_ptr => ofield_p
           rh_ptr => routehandle_pp
           input_arr_ptr => input_arr_p
           output_arr_ptr => output_arr_p
           num_cols_in = ncols_p_in
           num_cols_out = ncols_p_out
        else if ( any(dimids_in(:ndims) == did_ncol_p_in) .and. any(dimids_out(:ndims) == did_ncol_d_out) ) then
           ifld_ptr => ifield_p
           ofld_ptr => ofield_d
           rh_ptr => routehandle_pd
           input_arr_ptr => input_arr_p
           output_arr_ptr => output_arr_d
           num_cols_in = ncols_p_in
           num_cols_out = ncols_d_out
        else if ( any(dimids_in(:ndims) == did_ncol_d_in) .and. any(dimids_out(:ndims) == did_ncol_p_out) ) then
           ifld_ptr => ifield_d
           ofld_ptr => ofield_p
           rh_ptr => routehandle_dp
           input_arr_ptr => input_arr_d
           output_arr_ptr => output_arr_p
           num_cols_in = ncols_d_in
           num_cols_out = ncols_p_out
        else if ( any(dimids_in(:ndims) == did_ncol_d_in) .and. any(dimids_out(:ndims) == did_ncol_d_out) ) then
           ifld_ptr => ifield_d
           ofld_ptr => ofield_d
           rh_ptr => routehandle_dd
           input_arr_ptr => input_arr_d
           output_arr_ptr => output_arr_d
           num_cols_in = ncols_d_in
           num_cols_out = ncols_d_out
        endif

        numlevs = 0
        if  ( any(dimids_in(:ndims) == levdimid) ) then
           numlevs = nlevs
        elseif ( any(dimids_in(:ndims) == ilevdimid) ) then
           numlevs = nilevs
        endif

        ! regrid column dependent varaibles

        if (numlevs>0) then
           do t=1,ntimes
              do k=1,numlevs

                 call check( nf90_get_var( ncid_in, vid_in, input_arr_ptr, start=(/1,k,t/), count=(/num_cols_in,1,1/) ) )
                 call regrid( rh_ptr, ifld_ptr, ofld_ptr, input_arr_ptr, output_arr_ptr )
                 call check( nf90_put_var( ncid_out, vid, output_arr_ptr, start=(/1,k,t/), count=(/num_cols_out,1,1/) ) )

              end do
           end do
        else
           do t=1,ntimes

              call check( nf90_get_var( ncid_in, vid_in, input_arr_ptr, start=(/1,t/), count=(/num_cols_in,1/) ) )
              call regrid( rh_ptr, ifld_ptr, ofld_ptr, input_arr_ptr, output_arr_ptr )
              call check( nf90_put_var( ncid_out, vid, output_arr_ptr, start=(/1,t/), count=(/num_cols_out,1/) ) )

           end do
        endif

     endif
     nullify(rh_ptr)
     nullify(ifld_ptr)
     nullify(ofld_ptr)
     nullify(input_arr_ptr)
     nullify(output_arr_ptr)

  end do

  call check( nf90_close( ncid_in ) )
  call check( nf90_close( ncid_out ) )

  call ESMF_FieldRegridRelease(routeHandle_pp, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldRegridRelease(routeHandle_dp, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldRegridRelease(routeHandle_pd, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldRegridRelease(routeHandle_dd, rc=rc)
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldDestroy( ifield_p, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldDestroy( ifield_d, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldDestroy( ofield_p, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_FieldDestroy( ofield_d, rc=rc )
  if (rc/=ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_Finalize()
  deallocate(dimids_in)
  deallocate(dimids_out)

  deallocate( ilats_p, ilons_p )
  deallocate( ilats_d, ilons_d )
  nullify( ilats_p, ilons_p )
  nullify( ilats_d, ilons_d )

  deallocate( olats_p, olons_p )
  deallocate( olats_d, olons_d )
  nullify( olats_p, olons_p )
  nullify( olats_d, olons_d )

  deallocate( input_arr_p )
  deallocate( input_arr_d )
  deallocate( output_arr_p )
  deallocate( output_arr_d )

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
