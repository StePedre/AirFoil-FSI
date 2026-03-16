#!/bin.bash             # use bash as command interpreter
#$ -cwd                 # currentWorkingDirectory
#$ -N pimpleSim         # jobName
#$ -j y                 # merges output and errors
#$ -S /bin/bash         # scripting language
#$ -l h_rt=12:00:0      # jobDuration hh:mm:ss
#$ -q all.q            # queueName
#$ -pe mpi 32           # cpuNumber
#---------------------------------------------------------

### LOAD THE OPENFOAM ENVIRONMENT
module use module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv

#---------------------------------------------------------

### EXECUTE COMMANDS

rm -r log.*
rm -r processor*
rm -r postProcessing

decomposePar >& log.decomposeParSim
mpirun --hostfile machinefile.$JOB_ID potentialFoam -parallel -writephi >& log.potentialFoam
mpirun --hostfile machinefile.$JOB_ID pimpleFoam -parallel >& log.pimpleFoam
reconstructPar -latestTime >& log.reconstructParSim

echo End Parallel Run
