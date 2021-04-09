$training_data
DataTrain             train.dat
DataTest              test.dat
NumWrite              -1
RangeParameter        650.
RefEnergy             -650.
NumRanges             4
EnergyRanges          12000. 4000. 0. -100.
$end

$functional_form
NumAtoms              5
Symbols               O O H  N N
AtomGroups            1 1 2  3 3
ReadBasis             F
FactorOrder           3
TotalOrder            3
IMode                 1
NumChannels           1
FragmentGroups        1 1 1  2 2
$end

$fortran_execution
UseCL                 T
CommandLine           ../../src/pip.x < input > output
$end

# Can also send a submission script to job queue
CommandLine           sbatch pippy_submit.sh
