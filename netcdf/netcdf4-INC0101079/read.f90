program simple_xy_rd 
use netcdf 
implicit none 

integer :: ncid, varid,ierr 

ierr= nf90_open('//project/projectdirs/mpccc/fbench/greenland_4km_epsg3413_bheatflxPos_c161227.nc', NF90_NOWRITE, ncid) 
write(*,*) ierr,ncid
!ierr= nf90_open('greenland_4km_epsg3413_bheatflxPos_c161227.nc', NF90_NOWRITE, ncid) 
!write(*,*) ierr,ncid 
end 
