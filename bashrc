#
export HISTSIZE=10000
test -s ~/.alias && . ~/.alias || true
#
# set path to path where global sky models can be found
export RTS_CAT_PATH=$HOME/catalog/
# set up Python paths
export PYTHONPATH=${PYTHONPATH}:$HOME/bin
