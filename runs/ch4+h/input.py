$training_data
DataTrain             train.dat
DataTest              test.dat
NumWrite              2
Units                 10 12
RangeParameter        100.
RefEnergy             15.
NumRanges             3
EnergyRanges          300. 150. 50.
$end

$functional_form
NumAtoms              6
Symbols               C H H H H  H
AtomGroups            1 2 2 2 2  2
ReadBasis             F
FactorOrder           2
TotalOrder            2
IMode                 1
NumChannels           1
FragmentGroups        1 1 1 1 2  2
$end

$fortran_execution
UseCL                 T
CommandLine           ../../src/pip.x < input > output
$end

# Can also send a submission script to job queue
CommandLine           sbatch pippy_submit.sh
