#!/bin.bash             # use bash as command interpreter
#$ -cwd                 # currentWorkingDirectory
#$ -N pimpleMesh     # jobName
#$ -j y                 # merges output and errors
#$ -S /bin/bash         # scripting language
#$ -l h_rt=3:00:00      # jobDuration hh:mm:ss
#$ -q dida.q             # queueName
#$ -pe mpi 1            # cpuNumber
#---------------------------------------------------------

### LOAD THE OPENFOAM ENVIRONMENT
module use module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv

#---------------------------------------------------------

### EXECUTE COMMANDS

rm -r processor*
rm log.*

blockMesh >& log.blockMesh
surfaceFeatureExtract >& log.featureExtract
snappyHexMesh -overwrite >& log.snappyHexMesh
checkMesh >& log.checkMesh

echo End Run
