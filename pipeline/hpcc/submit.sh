logdir=hpcc/log
mkdir -p $logdir

QSUB="qsub -l ncpus={threads} "
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem}"
QSUB="$QSUB -o $logdir -e $logdir"

snakemake                               \
    -j 3000                             \
    --cluster-config hpcc/cluster.yaml  \
    --js hpcc/jobscript.sh              \
    --rerun-incomplete                  \
    --keep-going                        \
    --latency-wait 10                   \
    --use-conda                         \
    --cluster "$QSUB" $@
