#PBS -j oe
#PBS -V
#PBS -A feb-684-ac

cd $PBS_O_WORKDIR

Rscript --verbose cummerbund.R \
&> cummerbund.sh.log
