#!/bin.bash             # use bash as command interpreter
#$ -cwd                 # currentWorkingDirectory
#$ -N pimpleFull        # jobName
#$ -j y                 # merges output and errors
#$ -S /bin/bash         # scripting language
#$ -l h_rt=12:00:00      # jobDuration hh:mm:ss
#$ -q dida.q             # queueName
#$ -pe mpi 16           # cpuNumber
#---------------------------------------------------------

### LOAD THE OPENFOAM ENVIRONMENT
module use module use /software/spack/spack/share/spack/modules/linux-rocky8-sandybridge/
module load openfoam/2306-gcc-13.2.0-tnytlfv

#---------------------------------------------------------

### EXECUTE COMMANDS

rm -rf [1-9]*
rm -rf 0.*
rm log.*
rm -r processor*
rm -r postProcessing

blockMesh >& log.blockMesh
surfaceFeatureExtract >& log.featureExtract
decomposePar >& log.decomposeParMesh
mpirun --hostfile machinefile.$JOB_ID snappyHexMesh -parallel -overwrite >& log.snappyHexMesh
reconstructParMesh -constant >& log.reconstructParMesh
checkMesh >& log.checkMesh

rm -r processor*

decomposePar >& log.decomposeParSim
mpirun --hostfile machinefile.$JOB_ID pimpleFoam -parallel >& log.pimpleFoam
reconstructPar -constant >& log.reconstructParSim

echo End Parallel Run
