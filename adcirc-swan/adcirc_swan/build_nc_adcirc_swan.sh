#!bin/bash
# this is sh file assumes you have raw adcirc source already opened
# it also assumes this sh file is within the adcirc source
# open cmplrflags.mk
#	find the section with stampede
#	modify  NETCDFHOME to :=/opt/apps/intel17/netcdf/4.3.3.1/x86_64
#	save the file

#unzip adcirc_v51.zip
cd adcirc_v51
cd work
chmod +x config.guess
module load intel/17.0.4 netcdf python3 hdf5 zlib
make adcirc compiler=intel MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
make padcirc compiler=intel MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
cd ../swan
perl ./platform.pl # it's good practice to check to make sure it's using ifort
make punswan
make clobber
cd ../work
make padcswan compiler=intel MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
make adcprep compiler=intel SWAN=enable MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
make hstime compiler=intel MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
make aswip compiler=intel MACHINENAME=stampede NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable




