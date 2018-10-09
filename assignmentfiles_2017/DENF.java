import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.*;
import java.io.File;
import java.io.IOException;

public class DENF implements ContestSubmission
{

    // PRE-PROGRAMMED STUFF

	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public DENF()
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
	public int POP_SIZE = 241;
	public double SCALING_FACTOR = 0.15860938844837258;
	public double CROSSOVER_PROB = 0.37555635842296986;

    // Params for DE operators (different versions of algorithm)
	public int NR_PERTURBATION_VECTORS = 2;
	public String BASE_VECTOR = "rand";
	public String CROSSOVER_SCHEME = "bin";


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
   private double new_mutation_factor(List<Integer> candidates, double[] fitness_values,int evals) {
       // Function to calculate the self-adaptive scaling factor
        if (evals + POP_SIZE <= (int)((double)evaluations_limit_*0.9)){
            // use static scaling factor for first 90% of runs
            return SCALING_FACTOR;
        }else{ // compute NF in last 10% of runs

            double group1 = 0;
            double group2 = 0;

            // compute scores for both sides of difference vector
            int nr_candidates = candidates.size();
            for(int i = 0; i < nr_candidates/2;i++){
                group1 += fitness_values[candidates.get(i)];
                group2 += fitness_values[candidates.get(nr_candidates/2+i)];
            }

            // compute NF
            double NF = ( Math.abs(group1 - group2) ) /
                    ( Math.abs(group1) + Math.abs(group2) );

            return NF;
        }
    }


    private double[][] get_mutant_vector(double[][] pop, double[] fitness_values, int evals){ // Efi
        // function to create a new mutant population based on pop
        // result dims should be [POP_SIZE][PHENOTYPE_DIM]
        // create placeholder for population
    	double[][] mutants = new double[POP_SIZE][PHENOTYPE_DIM];
        
        // Only initialize best once
        // Initialization of variables
        int indexIndividual1 = 0;
        double[] individual1 = pop[indexIndividual1];
        if(BASE_VECTOR == "best") {
        	int best = 0;
        	for (int i = 1; i < POP_SIZE; i++){
        	    if ( fitness_values[i] > fitness_values[best] ) {
        	  	  best = i;
        	    }
        	}
        	indexIndividual1 = best;
        	individual1 = pop[indexIndividual1];
        }

    	// MAIN LOOP: create all mutants 
        for(int j = 0; j < POP_SIZE; j++){
        	
        	// Random base vectors
        	if(BASE_VECTOR == "rand") {
                indexIndividual1 = rnd_.nextInt(POP_SIZE);  
    	        individual1 = pop[indexIndividual1];
            }
	        // perturbation candidates
	        Set<Integer> candidates = new HashSet<Integer>();
	        for(int i = 0; i < (2*NR_PERTURBATION_VECTORS); i++) {
	        	int setLength = candidates.size();
	        	do{
	                int randomCandidateIndex = rnd_.nextInt(POP_SIZE);
	                if(randomCandidateIndex != indexIndividual1){
	                	candidates.add(randomCandidateIndex);
	                }
	            }while(candidates.size()==setLength);
	        }

            
	       
	        List<Integer> randomCandidateList = new ArrayList<Integer>(candidates);
			double NF = new_mutation_factor(randomCandidateList, fitness_values, evals);
            System.out.println(NF);
            double[] difference = pop[randomCandidateList.get(0)].clone();
	        for(int n = 0; n < PHENOTYPE_DIM; n++){
	        	for(int i=1; i<randomCandidateList.size(); i++){
	    			if(i < (randomCandidateList.size()/2)){
	    				difference[n] += pop[randomCandidateList.get(i)][n];
	    			}else{
	    				difference[n] -= pop[randomCandidateList.get(i)][n];
	    			}
	    		}
                double d = individual1[n] + SCALING_FACTOR * difference[n];
                    if (d < -5.0){
                        d = (-1.0*d)%5;
                        d = -1.0*d;
                    }
                    if (d > 5){
                        d = d % 5;
                    }
                    mutants[j][n] = d;
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
    	
    	if(CROSSOVER_SCHEME=="bin") {
	    	for (int i = 0; i < PHENOTYPE_DIM; i++) {
	    		double r_i = rnd_.nextDouble();
	    		if(r_i < CROSSOVER_PROB || i == R){
	    			child[i] = parent2[i];
	    		}else{
	    			child[i] = parent1[i];
	    		}
	        }
    	} else if(CROSSOVER_SCHEME=="exp") {
    		int i = 0;
    		while(rnd_.nextDouble()<CROSSOVER_PROB && i < PHENOTYPE_DIM) {
    			child[i] = parent2[i];
    			i++;
    		}
    		for (int j = i; i < PHENOTYPE_DIM; i++) {
    			child[j] = parent1[j];
    		}
    		child[R] = parent2[R];
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
            // Select parents
            // Apply crossover / mutation operators
        	double[][] mutant_pop =  get_mutant_vector(pop, fitness_scores, evals);
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