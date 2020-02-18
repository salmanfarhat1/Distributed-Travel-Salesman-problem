#include <stdio.h>
#include <stddef.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>



// Check whether there is a logger for C
// Or define own logger with var args to printf
#define DEBUG 0

int CHROMOSOME_LENGTH; // cities number

int POPULATION_SIZE; // nb of chromosomes
int MAX_NB_GENERATIONS;

int SUB_POPULATION_SIZE_PER_PROCESS;
int NB_OF_GENERATION_PER_PROCESS;

char cities_file[255];

typedef struct {
    int id;
    float x;
    float y;
} City;
City *cities;

float **distances_matrix;

typedef struct Chromosome{
    // genes represent a sequence of cities, i.e. a path
    int genes[29];
    float fitness;
} Chromosome;

void parse_arguments(int argc, char **argv);
void init_distances_matrix();
void fill_sample_chromosome(Chromosome *ptr_chromosome);
void print_chromosome(Chromosome *ptr_chromosome);
void fill_randomly_the_chromosome(Chromosome *chrom);
void print_the_received_population(Chromosome *population );
int getRandomNumber();
void print_population(Chromosome *population );
void swap_chromosomes( Chromosome *pop , int src , int dest);
void sort_population(Chromosome *population , int full_population);
void selection(Chromosome *pop);
int get_random_index_of_chrom();
void print_fitness(Chromosome *pop , int size);
void crossoverV2(Chromosome *pop);
void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro);
float percentage_of_difference(Chromosome chro1 , Chromosome chro2);
int if_exist(Chromosome *chrom , int x);
void mutation(Chromosome *pop );
void calculate_population_fitness(Chromosome *population);
void fix_reduntancy(Chromosome *population);
int check_if_same_Chom(Chromosome c1 , Chromosome c2);

void main(int argc, char **argv) {
    // TODO read parameters from configuration file
    Chromosome * population;
    Chromosome * sub_population;
    POPULATION_SIZE = 1000 ;
    NB_OF_GENERATION_PER_PROCESS =  atoi(argv[2]);
    SUB_POPULATION_SIZE_PER_PROCESS = 250;
    parse_arguments(argc, argv);
    init_distances_matrix();
    population = (Chromosome *)malloc(POPULATION_SIZE*sizeof(Chromosome));
    sub_population = (Chromosome *)malloc(SUB_POPULATION_SIZE_PER_PROCESS*sizeof(Chromosome));


    int numtasks, myID;
    MPI_Status status;


    MPI_Datatype newChroType;
    MPI_Datatype type[2] = { MPI_INT, MPI_FLOAT};
    int blocklen[2] = { CHROMOSOME_LENGTH, 1 };
    MPI_Aint disp[2];

    MPI_Init(&argc, &argv);
    disp[0] = offsetof(struct Chromosome , genes);
    disp[1] = offsetof(struct Chromosome , fitness);
    MPI_Type_create_struct(2, blocklen, disp, type, &newChroType);
    MPI_Type_commit(&newChroType);

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);

    // create new datatype for send and recv

    int *global;

    if(myID == 0 )
       global = (int *)malloc (8 * sizeof(int));

    int *local ;
    local= (int *)malloc (2 * sizeof(int));

    local[0] = myID;
    local[1] = myID+4;

    for(int i = 0 ; i < SUB_POPULATION_SIZE_PER_PROCESS  ; i++ ){
      fill_randomly_the_chromosome(&sub_population[i]);
    }
    sort_population(sub_population , SUB_POPULATION_SIZE_PER_PROCESS);

    // start measure time after initializing
    double start_time;
    double end_time;
    if(myID == 0)
       start_time = MPI_Wtime();
    int i = 0 ;
    while(i != NB_OF_GENERATION_PER_PROCESS)
    {
      i++;
      selection(sub_population);
      crossoverV2(sub_population);
      mutation(sub_population);
      calculate_population_fitness(sub_population);
      sort_population(sub_population ,  SUB_POPULATION_SIZE_PER_PROCESS);
    }
    // printf("\n\nrank = %d my sub fitness array is : " , myID );
    // print_fitness(sub_population , SUB_POPULATION_SIZE_PER_PROCESS);
    int c=0;

    // if(myID == 0 )
    // {
    //   for(int i = 0 ; i < 1000000000 ; i++){
    //     c++;
    //   }
    //   printf("C %d " , c);
    // }

    MPI_Gather( sub_population, SUB_POPULATION_SIZE_PER_PROCESS, newChroType, population, SUB_POPULATION_SIZE_PER_PROCESS, newChroType, 0, MPI_COMM_WORLD);

    if(myID == 0)
    {
      sort_population(population , POPULATION_SIZE);
      end_time = MPI_Wtime();
      printf("the best path founded is :");
      print_chromosome(&population[0]);
      // print_fitness(population , POPULATION_SIZE);
      // printf("\nNumber of generations: %d\npopulation size per process: %d\npopulation size of master: (since 4 nodes we r runing 4*%d): %d" ,NB_OF_GENERATION_PER_PROCESS , SUB_POPULATION_SIZE_PER_PROCESS , SUB_POPULATION_SIZE_PER_PROCESS , POPULATION_SIZE );
      printf("\nparallel Version using MPI: time measured on my personal pc: %f\n" , (end_time - start_time));
    }
    MPI_Finalize();
}

void print_population(Chromosome *population ){
  for(int i = 0 ; i < SUB_POPULATION_SIZE_PER_PROCESS ; i++){
    print_chromosome(&population[i]);
  }
}
void print_the_received_population(Chromosome *population ){
  for(int i = 0 ; i < POPULATION_SIZE ; i++){
    print_chromosome(&population[i]);
  }
}

void parse_arguments(int argc, char **argv) {
    // TODO can be made more elegant, namely by using getopt.h
    if (argc < 2) {
        printf("Must submit TSP locations test file name as argument.\n");
        exit(0);
    } else {
        strcpy(cities_file, argv[1]);

    }
}

void read_cities_from_file() {
    int nb_cities;
    FILE *fp;
    int i, line;
    char buffer[1024];
    float x, y;
    fp = fopen(cities_file, "r");
    for (i = 0; i < 4; i++) {
        fgets(buffer, 1024, fp);
    }
    if (fscanf(fp, "DIMENSION : %d", &nb_cities) == 0) {
        printf("Illegal TSP locations file format. Expecting the DIMENSION at line 5.\n");
        exit(0);
    }
    CHROMOSOME_LENGTH = nb_cities;
    for (i = 0; i < 2; i++) {
        fgets(buffer, 1024, fp);
    }
    cities = (City *) malloc(sizeof (City) * nb_cities);
    rewind(fp);
    for (i = 0; i < 7; i++) {
        fgets(buffer, 1024, fp);
    }
    while (fscanf(fp, "%d %f %f", &line, &x, &y) > 0 && line <= nb_cities) {
        cities[line - 1].id = line;
        cities[line - 1].x = x;
        cities[line - 1].y = y;
    }

    fclose(fp);
}

float get_distance(City city1, City city2) {
    return sqrt(pow(city1.x - city2.x, 2) + pow(city1.y - city2.y, 2));
}

void print_distances_matrix() {
    printf("\n\t");
    for (int i = 0; i < CHROMOSOME_LENGTH; ++i) {
        printf("\n");
        printf("|%d|\t", i);
        for (int j = 0; j < CHROMOSOME_LENGTH; j++) {
            printf("(%d,%d) = %.4f\t",i ,j, distances_matrix[i][j]);
        }
    }
}

void init_distances_matrix() {
    read_cities_from_file();
    distances_matrix = malloc(sizeof (float *) * CHROMOSOME_LENGTH);
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        distances_matrix[i] = calloc(CHROMOSOME_LENGTH, sizeof (float));
    }
    for (int i = 0; i < CHROMOSOME_LENGTH - 1; i++) {
        for (int j = i + 1; j < CHROMOSOME_LENGTH; j++) {
            float distance = get_distance(cities[i], cities[j]);
            distances_matrix[i][j] = distances_matrix[j][i] = distance;
        }
    }
    free(cities);
}

void calculate_fitness(Chromosome *ptr_chromosome ){
  float fitness= 0;
  int i= 0;
  for ( i = 0; i < CHROMOSOME_LENGTH-1; i++) {
      fitness += distances_matrix[ptr_chromosome->genes[i] -1][ptr_chromosome->genes[i+1] -1];
  }
  fitness += distances_matrix[ptr_chromosome->genes[i]-1 ][0];
  ptr_chromosome-> fitness= fitness;
}

void calculate_population_fitness(Chromosome *population){
  for(int i = 0 ; i < SUB_POPULATION_SIZE_PER_PROCESS ; i++){
    calculate_fitness(&population[i]);
  }
}

void fill_randomly_the_chromosome(Chromosome *chrom){
  int array[CHROMOSOME_LENGTH];
  // chrom->genes = malloc(CHROMOSOME_LENGTH * sizeof (int));
  for(int i = 0 ; i < CHROMOSOME_LENGTH; i++){
    array[i] = i+1;
  }
  int nbRand , tmp;
  for(int i = 0; i < CHROMOSOME_LENGTH ; i++){
    nbRand = getRandomNumber()%(CHROMOSOME_LENGTH -i);
    tmp = array[nbRand];
    array[nbRand] = array[CHROMOSOME_LENGTH - i - 1];
    array[CHROMOSOME_LENGTH - i - 1] = tmp;
    chrom->genes[i] = tmp;
  }
  calculate_fitness(chrom);
}

int getRandomNumber(){
	int seed = (unsigned)(time(NULL)+rand());
	srand(seed);
	return rand()%CHROMOSOME_LENGTH;
}

void getRandomChromosome(Chromosome *chro ){
	int array[CHROMOSOME_LENGTH] ;
	int temp = 0,nb = 0 ;
	for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++ ){
		array[i] = i+1;
	}
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		nb = getRandomNumber()%(CHROMOSOME_LENGTH - i);
		temp = array[nb];
		array[nb] = array[CHROMOSOME_LENGTH - i -1];
		array[CHROMOSOME_LENGTH - i -1 ] = temp;
		chro -> genes[i] = temp;
	}
	calculate_fitness(chro );
}

void print_chromosome(Chromosome *ptr_chromosome) {

    printf("\nFitness = %f\t, Genes = ", ptr_chromosome->fitness);
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        printf("%d_", ptr_chromosome->genes[i]);
    }
    printf("\n");
}

void sort_population(Chromosome *population , int full_population){
  for(int i = 0 ; i < full_population ; i++){
    for(int j = i+1 ; j < full_population ; j++){
      if(population[i].fitness >population[j].fitness)
        swap_chromosomes(population , i, j);
    }
  }
}

void swap_chromosomes( Chromosome *pop , int src , int dest){
	Chromosome chrom;
	chrom  = pop[src];
	pop[src] = pop[dest];
	pop[dest] = chrom;
}

void fix_reduntancy(Chromosome *population){
  Chromosome prev_chrome =population[0];
  for(int i = 1 ; i < SUB_POPULATION_SIZE_PER_PROCESS - 2  ; i++){
    if(check_if_same_Chom(prev_chrome , population[i]) == 1){
      fill_randomly_the_chromosome(&population[i]);
    }
    else{
      prev_chrome = population[i];
    }
  }
}

int check_if_same_Chom(Chromosome c1 , Chromosome c2){
  for(int i = 0 ; i < CHROMOSOME_LENGTH ; i++ ){
    if(c1.genes[i] != c2.genes[i])
      return 0;
  }
  return 1;
}

int get_random_index_of_chrom(){
  int seed = (unsigned)(time(NULL)+rand());
  srand(seed);
  return rand()%SUB_POPULATION_SIZE_PER_PROCESS;
}

void selection(Chromosome *pop ){

  int n = (40*SUB_POPULATION_SIZE_PER_PROCESS)/100;
  int randNb;
  for(int i = 0 ; i < (10*SUB_POPULATION_SIZE_PER_PROCESS)/100 ; i++ ){
    randNb =(SUB_POPULATION_SIZE_PER_PROCESS/2) + get_random_index_of_chrom()%(SUB_POPULATION_SIZE_PER_PROCESS/2);
    swap_chromosomes(pop ,n+i ,randNb );
  }
}



void print_fitness(Chromosome *population , int size ){
  printf("\n--------------------------------------------------------------------------------------Fitness---------------------------------------------------------------------------\n");
  for(int i = 0 ; i < size ; i++){
    printf("%.3f - " ,population[i].fitness );
  }
  printf("\n--------------------------------------------------------------------------------------END---------------------------------------------------------------------------\n");
}

float percentage_of_difference(Chromosome chro1 , Chromosome chro2){
	float sum = 0;
	for(int i = 0 ; i < CHROMOSOME_LENGTH ;i++){
		if(chro1.genes[i] != chro2.genes[i]){
			sum++;
		}
	}
	//printf("%f" , (sum*100)/CHROMOSOME_LENGTH);
	return (sum*100)/CHROMOSOME_LENGTH;
}

int if_exist(Chromosome *chrom , int x){
	//printf("chrom : %d , %d %d \n" ,chrom.genes[0] , chrom.genes[1] , x );
	for(int i =0 ; i < CHROMOSOME_LENGTH ; i++){
		if(x == chrom->genes[i]){
			return 1;
		}
	}
	return 0;
}

void crossoverV2(Chromosome *pop){
	int j=0 , nb=0;
	for(int i = 0 ; i <( SUB_POPULATION_SIZE_PER_PROCESS/2)  ; i++){
		do{
			nb= getRandomNumber()%(SUB_POPULATION_SIZE_PER_PROCESS/2);
		}while(nb == i && percentage_of_difference(pop[i] , pop[nb]) < 70);
		create_ChildV2(pop[i] , pop[nb] , &pop[(SUB_POPULATION_SIZE_PER_PROCESS/2) +i]);
	}
	//create_ChildV2(pop[(POPULATION_SIZE/2)-1] , pop[0] , &pop[POPULATION_SIZE-1]);
}

void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro){
	int n=0 , i=0,z=1;
	n = getRandomNumber()%(CHROMOSOME_LENGTH );

	for(i=0 ; i < CHROMOSOME_LENGTH ;i++){
		Chro->genes[i] =0 ;
	}

	for( i = n ; i < n+((CHROMOSOME_LENGTH*30)/100);i++){
		z=i%CHROMOSOME_LENGTH;
		Chro->genes[z]=p.genes[z];
	}
	int c=0;
	i=(z+1)%CHROMOSOME_LENGTH;
	while( i!=z){
		c = c%CHROMOSOME_LENGTH;
		//printf("\nin loop %d\n",i);
		if(if_exist(Chro , m.genes[c]) != 1){
			Chro->genes[i] = m.genes[c];
		}
		else{
			if(Chro->genes[i] == 0){
				while(if_exist(Chro , m.genes[c]) == 1){
					c++;
				}
				Chro->genes[i] = m.genes[c];
			}
		}
		c++;
		i++;
		i=i%CHROMOSOME_LENGTH;
	}
	//rintf("z %d",z);
}

void mutation(Chromosome *pop){
	int i,j,k,z;
	for( z =(50  *SUB_POPULATION_SIZE_PER_PROCESS/100) ; z < SUB_POPULATION_SIZE_PER_PROCESS - 1 ; z++){
		i = getRandomNumber()%(CHROMOSOME_LENGTH );
		j = getRandomNumber()%(CHROMOSOME_LENGTH );
		// k = getRandomNumber()%(POPULATION_SIZE -(20*POPULATION_SIZE/100));
		int temp = pop[z].genes[j];
    pop[z].genes[j] = pop[z].genes[i];
    pop[z].genes[i] = temp;
		// pop[(20*POPULATION_SIZE/100)+k].genes[j] = pop[(20*POPULATION_SIZE/100)+k].genes[i];
		// pop[(20*POPULATION_SIZE/100)+k].genes[i] = temp;
	}
}
