
i=100

while [ $i -lt 200 ]
do
  echo -e "\n-------------------------------------------------Generation Size $i--------------------------------------------------------------------------------\n"


  # echo -e "\n-------------------------------------------------Parallel Version using MPI of TSP-------------------------------------------------\n"
  mpicc -fopenmp -o test Parallel_project_with_collectiveComm.c -lm
  mpirun -n 4 ./test wi29.txt $i


  # echo -e "\n-------------------------------------------------Parallel Version using Both openMP and MPI of TSP-------------------------------------------------\n"
  mpicc -fopenmp -o test Parallel_project_without_collectiveComm.c -lm
  mpirun -n 4 ./test wi29.txt $i




   i=`expr $i + 100`
   echo -e "\n-------------------------------------------------End --------------------------------------------------------------------------------\n"

done









#  run on grid
#  ssh sfarhat@access.grid5000.fr
# scp wi29.txt sfarhat@access.grid5000.fr:grenoble/FirstTest/wi29.txt
#  oarsub -I -l nodes=1,walltime=03:00
#
