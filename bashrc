#
export HISTSIZE=10000
test -s ~/.alias && . ~/.alias || true
#
export MWA_OPS_DIR=/group/mwaops
export MWA_SCI_DIR=/group/mwasci
#
# switch to GCC/GNU enironment
module switch PrgEnv-cray PrgEnv-gnu
module unload gcc
# use gcc 4.8 for CUDA/RTS, use version 4.9 for boost
module load gcc/4.8.2
module load cray-libsci
module load cmake
module load scipy
module load lapack
module load cudatoolkit
module load astropy
module load cfitsio
module load boost
module load casacore
module load ephem
module load readline
module load psycopg2
module load gsl
#
module load matplotlib
# set path to CFITSIO library
export CFITSIO_DIR=/ivec/cle52/galaxy/devel/PrgEnv-gnu/5.2.25/cfitsio/3370
source ~/MIRRC.sh
# set path to RTS binary
export RTSDIR=$MWA_OPS_DIR/CODE/RTS
export RTSBIN=$RTSDIR/bin
# set path to path where global sky models can be found
export RTS_CAT_PATH="$HOME/catalog/"
# set up library paths
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MWA_OPS_DIR/CODE/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${MIRLIB}"
# set up binary paths
export PATH="${PATH}:${MWA_OPS_DIR}/CODE/bin:${RTSBIN}:${MIRBIN}"
export PATH="${PATH}:${MWA_SCI_DIR}/code/MWA_Tools/scripts"
# set up Python paths
export PYTHONPATH=${PYTHONPATH}:$MWA_SCI_DIR/code/MWA_Tools/
export PYTHONPATH=${PYTHONPATH}:$MWA_SCI_DIR/code/MWA_Tools/scripts
export PYTHONPATH=${PYTHONPATH}:$MWA_SCI_DIR/code/MWA_Tools/mwapy
export PYTHONPATH=${PYTHONPATH}:$MWA_SCI_DIR/code/MWA_Tools/configs
export PYTHONPATH=${PYTHONPATH}:$HOME/bin

# env vars for psql client
export PGPASSWORD=BowTie
export PGHOST=mwa-metadata01.pawsey.org.au
export PGDATABASE=mwa
export PGUSER=mwa
