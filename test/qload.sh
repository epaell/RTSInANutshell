#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --partition=copyq
#SBATCH --clusters=zeus
#SBATCH --account=mwasci
#SBATCH --export=NONE
#
module load mpifileutils
#
rm -f load.log
mpirun -n 1 ./load.py