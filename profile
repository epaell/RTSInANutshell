# Sample .profile for SuSE Linux
# rewritten by Christian Steinruecken <cstein@suse.de>
#
# This file is read each time a login shell is started.
# All other interactive shells will only read .bashrc; this is particularly
# important for language settings, see below.

test -z "$PROFILEREAD" && . /etc/profile || true

# Most applications support several languages for their output.
# To make use of this feature, simply uncomment one of the lines below or
# add your own one (see /usr/share/locale/locale.alias for more codes)
# This overwrites the system default set in /etc/sysconfig/language
# in the variable RC_LANG.
#
#export LANG=de_DE.UTF-8	# uncomment this line for German output
#export LANG=fr_FR.UTF-8	# uncomment this line for French output
#export LANG=es_ES.UTF-8	# uncomment this line for Spanish output


# Some people don't like fortune. If you uncomment the following lines,
# you will have a fortune each time you log in ;-)

#if [ -x /usr/bin/fortune ] ; then
#    echo
#    /usr/bin/fortune
#    echo
#fi

if [[ $LOADEDMODULES == *"PrgEnv-cray"* ]] 
then
    module swap PrgEnv-cray PrgEnv-gnu
fi

module load python
module load pawseytools

if [[ "$PAWSEY_CLUSTER" == "galaxy" ]]
then
	module load lapack
	module load astropy
	module unload gcc
	# use gcc 4.8 for CUDA/RTS, use version 4.9 for boost
	module load gcc/4.8.2
	module load cray-libsci
	module load linecache2
	module load traceback2
	module load unittest2
	module load h5py
	module load numpy/1.10.1
	module load python-dateutil/2.4.2
	module load pytz/2015.7
	module load matplotlib/1.5.0
	module load pyephem
	module load psycopg2
	module load healpy/1.9.1
	module load cudatoolkit
	module load gsl
	module load boost
	module load cmake
	module load cfitsio
	module load wcslib
	module load scipy
	module load readline
	module use /group/mwasci/software/cle52up04/modulefiles
	module load casacore/2.3.0
elif [[ "$PAWSEY_CLUSTER" == "magnus" ]]
	module load lapack
	module load astropy
	module unload gcc
	# use gcc 4.8 for CUDA/RTS, use version 4.9 for boost
	module load gcc/4.8.2
	module load cray-libsci
	module load linecache2
	module load traceback2
	module load unittest2
	module load h5py
	module load numpy/1.10.1
	module load python-dateutil/2.4.2
	module load pytz/2015.7
	module load matplotlib/1.5.0
	module load pyephem
	module load psycopg2
	module load healpy/1.9.1
	module load cudatoolkit
	module load gsl
	module load boost
	module load cmake
	module load cfitsio
	module load wcslib
	module load scipy
	module load readline
	module use /group/mwasci/software/cle52up04/modulefiles
	module load casacore/2.3.0
else
	module load astropy
	module load matplotlib
	module load numpy
	module load setuptools
fi

