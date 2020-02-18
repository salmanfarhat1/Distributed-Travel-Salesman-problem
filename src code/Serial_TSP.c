#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <math.h>

// Check whether there is a logger for C
// Or define own logger with var args to printf
#define DEBUG 0

int CHROMOSOME_LENGTH; // cities number
int POPULATION_SIZE; // nb of chromosomes
int MAX_NB_GENERATIONS;

char cities_file[255];

typedef struct {
    int id;
    float x;
    float y;
} City;
City *cities;

float **distances_matrix;

typedef struct {
    // genes represent a sequence of cities, i.e. a path
    int * genes;
    float fitness;
} Chromosome;

Chromosome * population;


void parse_arguments(int argc, char **argv);
void init_distances_matrix();
void fill_sample_chromosome(Chromosome *ptr_chromosome);
void print_chromosome(Chromosome *ptr_chromosome);
void fill_randomly_the_chromosome(Chromosome *chrom);
int getRandomNumber();
void print_population(Chromosome *population );
void swap_chromosomes( Chromosome *pop , int src , int dest);
void sort_population(Chromosome *population);
void selection(Chromosome *pop);
int get_random_index_of_chrom();
void print_fitness();
void crossoverV2(Chromosome *pop);
void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro);
float percentage_of_difference(Chromosome chro1 , Chromosome chro2);
int if_exist(Chromosome *chrom , int x);
void mutation(Chromosome *pop);
void calculate_population_fitness(Chromosome *population);


void main(int argc, char **argv) {
    // TODO read parameters from configuration file
    POPULATION_SIZE = 1000 ;
    MAX_NB_GENERATIONS =  atoi(argv[2]);
    parse_arguments(argc, argv);
    init_distances_matrix();
    population = (Chromosome *)malloc(POPULATION_SIZE*sizeof(Chromosome));
    clock_t begin = clock();

    for(int i = 0 ; i < POPULATION_SIZE  ; i++ ){
      fill_randomly_the_chromosome(&population[i]);
    }
    sort_population(population);
    // print_fitness();

    int  i = 0;
    while(i < MAX_NB_GENERATIONS){
      selection(population);
      crossoverV2(population);
      mutation(population);
      calculate_population_fitness(population);
      sort_population(population);
      i++;
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("the best path founded is :");
    print_chromosome(&population[0]);
    // print_fitness();
    // printf("\nNumber of generations: %d\npopulation size: %d" ,MAX_NB_GENERATIONS , POPULATION_SIZE );
    printf("\n seq Version 1: time measured on my personal pc: %f\n" , time_spent);

    // print_population(population);
}

void print_population(Chromosome *population ){
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
  for(int i = 0 ; i < POPULATION_SIZE ; i++){
    calculate_fitness(&population[i]);
  }
}

void fill_randomly_the_chromosome(Chromosome *chrom){
  int array[CHROMOSOME_LENGTH];
  chrom->genes = malloc(CHROMOSOME_LENGTH * sizeof (int));
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

void sort_population(Chromosome *population){
  for(int i = 0 ; i < POPULATION_SIZE ; i++){
    for(int j = i+1 ; j < POPULATION_SIZE ; j++){
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

int get_random_index_of_chrom(){
  int seed = (unsigned)(time(NULL)+rand());
  srand(seed);
  return rand()%POPULATION_SIZE;
}

void selection(Chromosome *pop){

  int n = (40*POPULATION_SIZE)/100;
  int randNb;
  for(int i = 0 ; i < (10*POPULATION_SIZE)/100 ; i++ ){
    randNb =(POPULATION_SIZE/2) + get_random_index_of_chrom()%(POPULATION_SIZE/2);
    swap_chromosomes(population ,n+i ,randNb );
  }
}

void print_fitness(){
  printf("\n--------------------------------------------------------------------------------------Fitness---------------------------------------------------------------------------\n");
  for(int i = 0 ; i < POPULATION_SIZE ; i++){
    printf("%.3f - " ,population[i].fitness );
  }
  printf("\n--------------------------------------------------------------------------------------Fitness---------------------------------------------------------------------------------\n");
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
	for(int i = 0 ; i <( POPULATION_SIZE/2)  ; i++){
		do{
			nb= getRandomNumber()%(POPULATION_SIZE/2);
		}while(nb == i && percentage_of_difference(pop[i] , pop[nb]) < 70);
		create_ChildV2(pop[i] , pop[nb] , &pop[(POPULATION_SIZE/2) +i]);
	}
	//create_ChildV2(pop[(POPULATION_SIZE/2)-1] , pop[0] , &pop[POPULATION_SIZE-1]);
}

void create_ChildV2(Chromosome p , Chromosome m , Chromosome *Chro){
	//printf("Here");
	int n=0 , i=0,z=1;
	n = getRandomNumber()%(CHROMOSOME_LENGTH );
	//printf("{%d , %d}", n , n+((CHROMOSOME_LENGTH*30)/100));

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
	int i,j,k;
	for(int z =0 ; z < 5 ; z++){
		i = getRandomNumber()%(CHROMOSOME_LENGTH );
		j = getRandomNumber()%(CHROMOSOME_LENGTH );
		k = getRandomNumber()%(POPULATION_SIZE -(20*POPULATION_SIZE/100));
		int temp = pop[(20*POPULATION_SIZE/100)+k].genes[j];
		pop[(20*POPULATION_SIZE/100)+k].genes[j] = pop[(20*POPULATION_SIZE/100)+k].genes[i];
		pop[(20*POPULATION_SIZE/100)+k].genes[i] = temp;
	}
}
