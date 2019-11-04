#!/bin/bash -l
#SBATCH -t 00:20:00
#SBATCH -J Build_SC
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --array=0-80                                                            
#SBATCH --no-requeue                                                            
#SBATCH --mem-per-cpu=16000                                                     
#SBATCH --constraint=hsw

bulkStart=3200 # set to intial bulk
bulkEnd=4000 # set to final bulk

# How many jobs? SLURM_ARRAY_TASK_COUNT does not work on all systems            
# so calculate job count (or set it manually to match the array                 
# argument given to sbatch).
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 ))

# find job array index                                                          
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

bulkEndC=$(( $bulkEnd + 1 )) # Need to iterate to 1 past final bulk
totalBulk=$(( $bulkEndC - $bulkStart )) # Total bulk to calculate
increment=$(( $totalBulk / $jobcount )) # amount of bulk per job (rounded down)

# Calculate remainder
remainder=$(( $totalBulk - $jobcount * $increment ))

start=$(( $bulkStart + $index * $increment ))
end=$(( $start + $increment ))

# Remainder frames are divided out evenly among tasks
if [ $index -lt $remainder ];
then
    start=$(( $start + $index ))
    end=$(( $end + $index + 1 ))
else
    start=$(( $start + $remainder ))
    end=$(( $end + $remainder ))
fi;

# Ensure final job gets correct last bulk
if [ $SLURM_ARRAY_TASK_ID -eq $SLURM_ARRAY_TASK_MAX ];
then
    echo Verifying final bulk: $end $bulkEndC
    end=$bulkEndC
fi;

module load mayavi2
export OMP_NUM_THREADS=24
export PYTHONPATH=$PYTHONPATH:~/appl_taito/analysator
export PTNONINTERACTIVE=1 # Turns off x-windowing and Mayavi2 loading
python array_SC.py $start $end
echo Job $SLURM_ARRAY_TASK_ID complete.

