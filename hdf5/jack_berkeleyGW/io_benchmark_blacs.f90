! Compile with:
! ftn io_benchmark_blacs.f90 -L/opt/cray/hdf5-parallel/1.8.8/pgi/119/lib -lhdf5_fortran -lhdf5 -lz -o io_bench_blacs.x

program io_benchmark

  use hdf5
  implicit none

  include 'mpif.h'

  real(kind(1D0)), allocatable :: myarray(:,:,:)

  real(kind(1.0d0)) :: starttime, currenttime, mygbytes

  CHARACTER(len=32) :: arg
  integer :: rank,error,mype,npes,i,j,k,nwrites

  integer :: nmtx, nprow, npcol, myrow, mycol, nbl, npr, npc, irow, icol
  CHARACTER(len=120) :: ddir

  CHARACTER(len=120) :: indfname
  CHARACTER(len=120) :: colfname
  CHARACTER(len=6) :: ind1 = "ind.h5", col1 = "col.h5"

  CALL getarg(1, arg)
  READ(arg,*) nmtx
  CALL getarg(2, arg)
  READ(arg,*) nwrites
  CALL getarg(3, ddir)
  !READ(arg,*) ddir
  !PRINT *,"ddir is ", ddir

  indfname = trim(ddir) // ind1
  colfname = trim(ddir) // col1
  indfname = trim(indfname)
  colfname = trim(colfname)

  !PRINT *, indfname
  !PRINT *, ddir

  !call EXIT(0) 

  call MPI_INIT(error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, error)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mype, error)

  call layout_scalapack(nmtx, nbl, npes, mype, nprow, npcol, myrow, mycol)

  if (mype .eq. 0) then
    WRITE (6,*) 'Running on npes: ',npes
    WRITE (6,*) 'nmtx: ',nmtx
    WRITE (6,*) 'nwrites: ',nwrites
    WRITE (6,*) 'nbl: ',nbl
    WRITE (6,*) 'nprow: ',nprow
    WRITE (6,*) 'npcol: ',npcol
  endif
  call flush(6)

  npr = 0
  npc = 0

  do j = 1, nmtx
    irow=MOD(INT(((j-1)/nbl)+1D-6),nprow)
    if (irow .eq. myrow) then
      npr = npr+1
    endif
  enddo

  do k = 1, nmtx
    icol=MOD(INT(((k-1)/nbl)+1D-6),npcol)
    if(icol .eq. mycol) then
      npc = npc+1
    endif
  enddo

  mygbytes = 1D0*nmtx*1D0*nmtx*1D0*nwrites*8D0/1024D0/1024D0/1024D0

  !write(6,*) mype, nmtx, nwrites, npr, npc, nbl, nprow, npcol, myrow, mycol

  allocate(myarray(npr,npc,1))
  myarray(:,:,:) = mype

  ! HDF5 Collective Test
  call write_hdf5(nmtx,nwrites,npes,mype,0,mygbytes,myarray,nbl,npcol,nprow,myrow,mycol,npr,npc,colfname)

  ! HDF5 Independent Test
  call write_hdf5(nmtx,nwrites,npes,mype,1,mygbytes,myarray,nbl,npcol,nprow,myrow,mycol,npr,npc,indfname)

  ! File Per Proc Test
  call write_fpp(nmtx,nwrites,npes,mype,mygbytes,myarray,npr,npc)

  deallocate(myarray)

  call MPI_FINALIZE(error)

end program

subroutine timget(wall)
  real(kind(1.0d0)) :: wall

  integer :: values(8)

  call date_and_time(VALUES=values)
  wall=((values(3)*24.0d0+values(5))*60.0d0 &
    +values(6))*60.0d0+values(7)+values(8)*1.0d-3

  return
end subroutine timget

subroutine write_hdf5(nmtx,nwrites,npes,mype,iflag,mygbytes,data,nbl,npcol,nprow,myrow,mycol,npr,npc,fname)

  use hdf5
  implicit none

  include 'mpif.h'

  integer :: rank,error,mype,npes,i,iflag,nwrites
  integer :: nmtx,comm,info,ijk,nbl,npcol,nprow,myrow,mycol,npr,npc
  CHARACTER(len=120) :: fname
  !write(6,*) "fname:",fname
  real(kind(1D0)) :: data(nmtx,nmtx,1)
  real(kind(1.0d0)) :: starttime, endtime, tottime, mygbytes

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: data_id       ! Dataset identifier
  integer(HID_T) :: plist_id
  integer(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: memspace     ! Dataspace identifier in file
  integer(HID_T) :: infospace
  integer(HSIZE_T) :: count(3), countm(3), offset(3), offsetm(3), stride(3), block(3)

  CHARACTER(len=32) :: arg

  tottime=0D0

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  call h5open_f(error)
  !PRINT *, fname
  !write(6,*) "in hdf", nmtx

  if (mype .eq. 0) then
    PRINT *, fname
    if (iflag .eq. 0) call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
    if (iflag .eq. 1) call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

    rank = 3
    count(1) = nmtx
    count(2) = nmtx
    count(3) = nwrites

    call h5screate_simple_f(rank, count, filespace, error)
    call h5dcreate_f(file_id, 'data', H5T_NATIVE_DOUBLE, filespace, &
                      data_id, error)

    call h5dclose_f(data_id, error)
    call h5sclose_f(filespace, error)
    call h5fclose_f(file_id, error)
  endif

  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, comm, info, error)
  if (iflag .eq. 0) call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  if (iflag .eq. 1) call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  call h5pclose_f(plist_id,error)

  call h5dopen_f(file_id, 'data', data_id, error)

  !write(6,*) mype, "Waiting at barrier"

  call MPI_BARRIER(MPI_COMM_WORLD,error)

  !write(6,*) mype, "past barrier"

  do ijk=1,nwrites

    rank = 3

    countm(1) = npr - mod(npr,nbl)
    countm(2) = npc - mod(npc,nbl)
    countm(3) = 1

    offsetm(:) = 0

    count(1) = npr/nbl
    count(2) = npc/nbl
    count(3) = 1

    offset(1) = (myrow) * nbl
    offset(2) = (mycol) * nbl
    offset(3) = ijk-1

    stride(1) = nprow * nbl
    stride(2) = npcol * nbl
    stride(3) = 1

    block(1) = nbl
    block(2) = nbl
    block(3) = 1

    call h5screate_simple_f(rank, countm, memspace, error)
    call h5dget_space_f(data_id,filespace,error)

    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (iflag .eq. 0) then
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,error)
    call timget(starttime)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
        xfer_prp = plist_id)
    !call MPI_BARRIER(MPI_COMM_WORLD,error)
    call timget(endtime)

    tottime=tottime+endtime-starttime

    call h5pclose_f(plist_id, error)
    call h5sclose_f(memspace, error)
    call h5sclose_f(filespace, error)

  enddo

  if (mype .eq. 0) then
    if (iflag .eq. 0) write(6,*) "HDF5 Collective time: ", tottime, mygbytes/tottime, "GB/s"
    if (iflag .eq. 1) write(6,*) "HDF5 Independent time: ", tottime, mygbytes/tottime, "GB/s"
  endif

  call h5dclose_f(data_id, error)
  call h5fclose_f(file_id, error)

end subroutine write_hdf5

subroutine write_fpp(nmtx,nwrites,npes,mype,mygbytes,data,npr,npc)

  implicit none

  include 'mpif.h'

  integer :: nmtx,mype,npes,error,ijk,nwrites,npr,npc

  real(kind(1D0)) :: data(npr,npc,1)
  real(kind(1.0d0)) :: starttime, endtime,mygbytes

  call MPI_BARRIER(MPI_COMM_WORLD,error)
  call timget(starttime)
  do ijk = 1, nwrites
    write(2000+mype) data(:,:,1)
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  call timget(endtime)
 
  if (mype .eq. 0) write(6,*) "FPP time: ", endtime-starttime, mygbytes/(endtime-starttime), "GB/s"

end subroutine write_fpp

subroutine layout_scalapack(matsize, nbl, nproc, mype, nprow, npcol, myrow, mycol)
    integer, intent(in) :: matsize
    integer, intent(out) :: nbl
    integer, intent(in) :: nproc,mype
    integer, intent(out) :: nprow, npcol, myrow, mycol

    integer :: i,j

!------------------
! Find processor grid

    nprow = int(sqrt(dble(nproc) + 1.0d-6))
    !PRINT *,'nprow is ',nprow
    do i = nprow, 1, -1
      if(mod(nproc, i) .eq. 0) exit
    enddo
    
    nprow = i
    npcol = nproc/nprow

!-------------------
! Now for the block size

    nbl = min(1024, matsize/(max(nprow, npcol)))

!-------------------
! Ensure nonzero

    nbl = max(nbl, 1)

!-------------------
! Ensure nonzero

    nbl = max(nbl, 1)

    icount = 0
    do i = 1, nprow
      do j = 1, npcol
        if (icount .eq. mype) then
          myrow=i-1
          mycol=j-1
        endif
        icount = icount + 1
      enddo
    enddo

    return
end subroutine layout_scalapack
