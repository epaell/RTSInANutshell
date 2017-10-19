#
export HISTSIZE=10000
test -s ~/.alias && . ~/.alias || true
#
export MWA_OPS_DIR=/group/mwaops
export MWA_SCI_DIR=/group/mwasci
# set path to CFITSIO library
export CFITSIO_DIR=/ivec/cle52/galaxy/devel/PrgEnv-gnu/5.2.25/cfitsio/3370
source ~/MIRRC.sh
# set path to RTS binary
export RTSDIR=$MWA_OPS_DIR/CODE/RTS
export RTSBIN=$RTSDIR/bin
# set path to path where global sky models can be found
export RTS_CAT_PATH=$HOME/catalog/
# set up library paths
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${MIRLIB}"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MWA_SCI_DIR/code/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MWA_OPS_DIR/CODE/lib"
# set up binary paths
export PATH="${PATH}:${MWA_SCI_DIR}/code/MWA_Tools/scripts"
export PATH="${PATH}:${MWA_SCI_DIR}/code/bin"
export PATH="${PATH}:$RTSBIN"
export PATH="${PATH}:${MIRBIN}"
# set up Python paths
export PYTHONPATH=${PYTHONPATH}:$MWA_SCI_DIR/code/lib/python2.6/site-packages/
export PYTHONPATH=${PYTHONPATH}:$HOME/bin
