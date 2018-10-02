import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.*;
import java.io.File;
import java.io.IOException;

public class RDE implements ContestSubmission
{

    // PRE-PROGRAMMED STUFF

	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public RDE()
	{
		rnd_ = new Random();
	}
	
	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;
		
		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

        String eval_type = props.getProperty("Evaluation");

        System.out.println(eval_type);
		// Do sth with property values, e.g. specify relevant settings of your algorithm
        if(isMultimodal){
            // Do sth
        }else{
            // Do sth else
        }
    }

    // MAIN ALGORITH PARAMETERS

    // Static for assignment
    public int PHENOTYPE_DIM = 10; // each phenotype has 10 dimensions in our assignment
    public int DIM_LOWER_BOUND = -5; // each dimension ranges from [-5, 5]
    public int DIM_UPPER_BOUND = 5;

    // Changeable params
	public int POP_SIZE = 29;
	public double SCALING_FACTOR = 0.015165767714657719;
	public double CROSSOVER_PROB = 0.23955957811705852;

    // Params for DE operators (different versions of algorithm)
	public int NR_PERTURBATION_VECTORS = 3;
	public String BASE_VECTOR = "best";
	public String CROSSOVER_SCHEME = "exp";


    // HELPER FUNCTIONS FOR MAIN
    private double[][] init_population(int pop_size, int phenotype_dim){
        // function to intialize a random population of size pop_size

        double[][] pop = new double[pop_size][phenotype_dim];

        for(int i = 0; i < pop_size;i++){
            // init phenotype
            double[] pheno = new double[phenotype_dim];

            for(int j = 0; j < phenotype_dim; j++){
                // initialize values randomly
                pheno[j] = DIM_LOWER_BOUND + (DIM_UPPER_BOUND - DIM_LOWER_BOUND) * rnd_.nextDouble();
                }
            pop[i] = pheno;
        }
        return pop;
    }
    
    private double[] eval_pop(double[][] pop){
        // function to evaluate each phenotype in a population
    	double[] fitness_values = new double[POP_SIZE];
    	
    	for (int i = 0; i < POP_SIZE; i++) {
    		fitness_values[i] = (double) evaluation_.evaluate(pop[i]);
		}

        return fitness_values;
    }


    //Get the fittest individual
   public double getFittest(double[] fitness_values) {
       double maxFit = Double.MIN_VALUE;
       int maxFitIndex = 0;
       for (int i = 0; i < POP_SIZE; i++) {
           if (maxFit <= fitness_values[i]) {
               maxFit = fitness_values[i];
               maxFitIndex = i;
           }
       }
       return maxFit;
   }

   //Get least fit individual
   public double getLeastFittest(double[] fitness_values) {
       double minFitVal = Double.MAX_VALUE;
       int minFitIndex = 0;
       for (int i = 0; i < POP_SIZE; i++) {
           if (minFitVal >= fitness_values[i]) {
               minFitVal = fitness_values[i];
               minFitIndex = i;
           }

       }
       return minFitVal;
   }

   // Diversity
   public double getDiversity(double[] fitness_values){
       double diversity;
       double fittest;
       double least_fittest;

       // For now: max fitness -  min fitness
       fittest = getFittest(fitness_values);
       least_fittest = getLeastFittest(fitness_values);
       diversity = fittest-least_fittest;
       return diversity;

   }


    private double[][] get_mutant_vector(double[][] pop){ // Efi
        // function to create a new mutant population based on pop
        // result dims should be [POP_SIZE][PHENOTYPE_DIM]

        // create placeholder for population
        double[][] mutants = new double[POP_SIZE][PHENOTYPE_DIM];

        // pick 3 invividuals from the population  P[i1],P[i2],P[i3]
        // where i1!= i2 != i3
        int a,b,c;

        for(int j = 0; j < POP_SIZE; j++){
            a = rnd_.nextInt(POP_SIZE);
            
            do{
                b = rnd_.nextInt(POP_SIZE);
            }while(b==a);
            do{
                c = rnd_.nextInt(POP_SIZE);
            }while(c == a || c == b);

            // create three agent individuals
            double[] individual1 = pop[a];
            double[] individual2 = pop[b];
            double[] individual3 = pop[c];

            // mutation process
            // create difference vector based on NR_PERTURBATION_VECTORS

            for(int n = 0; n < PHENOTYPE_DIM; n++){
                mutants[j][n] = (individual1[n] + SCALING_FACTOR
                        * (individual2[n] - individual3[n]));
            }
        }
        // Sample base vector based on BASE_VECTOR 
        return mutants;
    }

    private double[][] get_trial_vector(double[][] pop, double[][] mutants){ // Hans
        // function to create the trial vectors T
        // apply crossover based on CROSSOVER_SCHEME
        // return T
	    double[][] T = new double[POP_SIZE][PHENOTYPE_DIM];
	    for(int i = 0; i < POP_SIZE; i++) {
            T[i] = crossover(pop[i], mutants[i]);
	    }

	    return T;
    }

    private double[] crossover(double[] parent1, double[] parent2){
        // function to perform crossover between two parents, return new child
    	double[] child = new double[PHENOTYPE_DIM];
    	int R = rnd_.nextInt(PHENOTYPE_DIM); //index on which at least crossover is applied
    	
    	for (int i = 0; i < PHENOTYPE_DIM; i++) {
    		double r_i = rnd_.nextDouble();
    		if(r_i < CROSSOVER_PROB || i == R){
    			child[i] = parent2[i];
    		}else{
    			child[i] = parent1[i];
    		}
        }
        return child;
    }

    private List<? extends Object> survival_selection(double[][] parents, double[] parents_fitness, double[][] children){ // Niels
         // function to select new surivors based on fitness
        double[][] survivors = new double[POP_SIZE][PHENOTYPE_DIM];
        double[] survivor_fitness = new double[POP_SIZE];
        // evaluate children
        double[] child_fitness = eval_pop(children);
        // compare each parent to child and save fittest in new population
        for(int i = 0; i < parents_fitness.length; i++){
            if(parents_fitness[i] > child_fitness[i]){
                survivors[i] = parents[i].clone();
                survivor_fitness[i] = parents_fitness[i];
            }else{
                survivors[i] = children[i].clone();
                survivor_fitness[i] = child_fitness[i];
            }
        }
        // return list of objects, to be converted back outside function
        return Arrays.asList(survivors, survivor_fitness);
    }

    // MAIN FUNCTION
	public void run()
	{
        int evals = POP_SIZE;
        // init population
        double[][] pop = init_population(POP_SIZE, PHENOTYPE_DIM);
        double diversity;
        // calculate fitness
        double[] fitness_scores = eval_pop(pop);
        while(evals<evaluations_limit_){

            //decreasing scaling factor
            //for (int i = 0; i < PHENOTYPE_DIM; i++) {
            //SCALING_FACTOR = SCALING_FACTOR[i]*0.999;
            //}

            diversity = getDiversity(fitness_scores);
            //System.out.println(String.valueOf(diversity));
            // Select parents
            // Apply crossover / mutation operators
        	double[][] mutant_pop =  get_mutant_vector(pop);
        	double[][] trial_pop = get_trial_vector(pop, mutant_pop);
        	// Selection
        	List res = survival_selection(pop, fitness_scores, trial_pop);
        	pop = (double[][])res.get(0);
        	fitness_scores = (double[])res.get(1);
        	evals += POP_SIZE;

        }
	}
	
	public void main(String[] argv){
    	run();
    }
}