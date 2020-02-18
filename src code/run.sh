
i=100

while [ $i -lt 200 ]
do
  echo -e "\n-------------------------------------------------Generation Size $i--------------------------------------------------------------------------------\n"
  # echo -e "\n-------------------------------------------------Parallel Version using Both openMP and MPI of TSP-------------------------------------------------\n"
  mpicc -fopenmp -o test Parallel_project_with_openMP.c -lm
  mpirun -n 4 ./test wi29.txt $i

  # echo -e "\n-------------------------------------------------Parallel Version using MPI of TSP-------------------------------------------------\n"
  mpicc -o test Parallel_project_without_openMP.c -lm
  mpirun -n 4 ./test wi29.txt $i

  # echo -e "\n-------------------------------------------------Sequential Version of TSP-------------------------------------------------\n"
  gcc -o serial Serial_TSP.c -lm
  ./serial wi29.txt $i


   i=`expr $i + 100`
   echo -e "\n-------------------------------------------------End --------------------------------------------------------------------------------\n"

done


killall test



#  run on grid
#  ssh sfarhat@access.grid5000.fr
# scp wi29.txt sfarhat@access.grid5000.fr:grenoble/FirstTest/wi29.txt
# scp *  sfarhat@access.grid5000.fr:grenoble/project/
#  oarsub -I -l nodes=1,walltime=03:00
#
