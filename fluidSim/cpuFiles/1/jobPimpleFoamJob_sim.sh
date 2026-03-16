#!/bin.bash             # use bash as command interpreter
#$ -cwd                 # currentWorkingDirectory
#$ -N pimpleSim         # jobName
#$ -j y                 # merges output and errors
#$ -S /bin/bash         # scripting language
#$ -l h_rt=12:00:0      # jobDuration hh:mm:ss
#$ -q dida.q            # queueName
#$ -pe mpi 1           # cpuNumber
#---------------------------------------------------------

### LOAD THE OPENFOAM ENVIRONMENT
module use module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv

#---------------------------------------------------------

### EXECUTE COMMANDS

rm -r log.*
rm -r processor*
rm -r postProcessing

potentialFoam -writephi >& log.potentialFoam
pimpleFoam >& log.pimpleFoam

echo End Run
