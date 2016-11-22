#!/bin/bash
#files=($1)

outdir=$HOME/tmhbxx3/scratch/XMLhttp/QuickXMLs

mkdir -p $outdir
for i in `seq 0 29`;
do
script=$(mktemp)
    cat > $script <<EOF
#!/bin/bash
#PBS -N GSMXMLQcuik$i
#PBS -o $outdir/$i.out
#PBS -e $outdir/$i.err
#PBS -m e
#PBS -M bxia@houstonmethodist.org
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

module load python/2.7.11
cd /home/tmhbxx3/scratch/XMLhttp/
python ./xmlHttpGSMDownloader.py $i

EOF
qsub $script
done
