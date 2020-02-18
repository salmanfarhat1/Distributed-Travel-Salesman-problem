
i=100

while [ $i -lt 1000 ]
do
  echo -e "\n-------------------------------------------------Generation Size $i--------------------------------------------------------------------------------\n"


  # echo -e "\n-------------------------------------------------Parallel Version using MPI of TSP-------------------------------------------------\n"
  mpicc -fopenmp -o test1 Parallel_project_with_collectiveComm.c -lm
  mpirun --mca pml ob1 --mca btl ^openib -n 80 -hostfile machinefile ./test1 wi29.txt $i


  # echo -e "\n-------------------------------------------------Parallel Version using Both openMP and MPI of TSP-------------------------------------------------\n"
  mpicc -fopenmp -o test2 Parallel_project_without_collectiveComm.c -lm
  mpirun --mca pml ob1 --mca btl ^openib -n 80 -hostfile machinefile ./test2 wi29.txt $i




  # echo -e "\n-------------------------------------------------Sequential Version of TSP-------------------------------------------------\n"
  # gcc -o serial Serial_TSP.c -lm
  # ./serial wi29.txt $i


   i=`expr $i + 100`
   echo -e "\n-------------------------------------------------End --------------------------------------------------------------------------------\n"

done









#  run on grid
#  ssh sfarhat@access.grid5000.fr
# scp wi29.txt sfarhat@access.grid5000.fr:grenoble/FirstTest/wi29.txt
#  oarsub -I -l nodes=1,walltime=03:00
#
