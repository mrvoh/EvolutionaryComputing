import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.*;
import java.io.File;
import java.io.IOException;

public class DifferentialEvolution implements ContestSubmission
{

    // PRE-PROGRAMMED STUFF

	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
	
	public DifferentialEvolution()
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
    public int POP_SIZE = 100; // mu
    public double SCALING_FACTOR = 0.5; //F
    public double CROSSOVER_PROB = 0.5; // Cr
    public double CROSSOVER_RATE = 0.1; // Pc

    // Params for DE operators (different versions of algorithm)
    public int NR_PERTURBATION_VECTORS = 1;
    public String BASE_VECTOR = "random";
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
    private double[] init_fitness_scores(int pop_size){
        // function to initialize an array of fitnesses with only zeros

        double[] fitness = new double[pop_size];// java initializes an array with zeros
        return fitness;
    }

    private double[] eval_pop(double[][] pop){ //Ronald
        // function to evaluate each phenotype in a population

        double[] dummy = new double[POP_SIZE];
        return  dummy;
    }

    private double[][] get_mutant_vector(double[][] pop){ // Efi
        // function to create a new mutant population based on pop
        // result dims should be [POP_SIZE][PHENOTYPE_DIM]

        // Sample base vector based on BASE_VECTOR 

        // create difference vector based on NR_PERTURBATION_VECTORS

        // create final population
        double[][] dummy = new double[POP_SIZE][PHENOTYPE_DIM];
        return dummy;
    }

    private double[][] get_trial_vector(double[][] pop, double[][] mutants){ // Hans
        // function to create the trial vector T

        // apply crossover based on CROSSOVER_SCHEME


        // return T
        double[][] dummy = new double[POP_SIZE][PHENOTYPE_DIM];
        return dummy;
    }

    private double[] crossover(double[] parent1, double[] parent2){ // Ronald
        // function to perform crossover between two parents, return new child

        double[] dummy = new double[PHENOTYPE_DIM];
        return dummy;
    }

    private List<Object> survival_selection(double[][] parents, double[] parents_fitness, double[][] children){ // Niels
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
        

        int evals = 0;
        // init population
        // calculate fitness
        while(evals<evaluations_limit_){
            // Select parents
            // Apply crossover / mutation operators
            double child[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            // Check fitness of unknown fuction
            Double fitness = (Double) evaluation_.evaluate(child);
            evals++;
            // Select survivors
            // EXAMPLE USAGE OF survival_selection()
            // List res = survival_selection(pop, fitness, pop);
            // double[][] survivors = (double[][])res.get(0);
            // double[] surv_fitness = (double[])res.get(1);
        }

	}
	
	public void main(String[] argv){
    	run();
    }
}
