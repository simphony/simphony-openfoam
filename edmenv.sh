export PATH=/usr/lib64/openmpi/bin/:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib/:$LD_LIBRARY_PATH
# openfoam bashrc script has the annoying property of removing our current path, 
# leaving us with nothing of our custom installations
oldpath=$PATH
. /opt/OpenFOAM-2.3.0/etc/bashrc
export PATH=$oldpath:$PATH
