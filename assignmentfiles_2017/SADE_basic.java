import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.*;
import java.io.File;
import java.io.IOException;

public class SADE_basic implements ContestSubmission
{

    // PRE-PROGRAMMED STUFF

    Random rnd_;
    ContestEvaluation evaluation_;
    private int evaluations_limit_;

    public SADE_basic()
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
	public int POP_SIZE = 313;
	public double SCALING_FACTOR = 0.6170124706444252;
    //public double[] SCALING_FACTOR_MULTI = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
	public double CROSSOVER_PROB = 0.9377236940855606;


    // Params for DE operators (different versions of algorithm)
	public int NR_PERTURBATION_VECTORS = 1;
	public String BASE_VECTOR = "rand";
	public String CROSSOVER_SCHEME = "exp";
    //public String SCALING_FACTOR_SCHEME = "multi"; // multi or 1d


    // HELPER FUNCTIONS FOR MAIN
    private double[][] init_population(int pop_size, int phenotype_dim, double l, double u){
        // function to intialize a random population of size pop_size

        double[][] pop = new double[pop_size][phenotype_dim];

        for(int i = 0; i < pop_size;i++){
            // init phenotype
            double[] pheno = new double[phenotype_dim];

            for(int j = 0; j < phenotype_dim; j++){
                // initialize values randomly
                pheno[j] = l + (u - l) * rnd_.nextDouble();
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
    public int getFittest(double[] fitness_values) {
        double maxFit = Double.MIN_VALUE;
        int maxFitIndex = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            if (maxFit <= fitness_values[i]) {
                maxFit = fitness_values[i];
                maxFitIndex = i;
            }
        }
        return maxFitIndex;
    }

    //Get least fit individual
    public int getLeastFittest(double[] fitness_values) {
        double minFitVal = Double.MAX_VALUE;
        int minFitIndex = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            if (minFitVal >= fitness_values[i]) {
                minFitVal = fitness_values[i];
                minFitIndex = i;
            }

        }
        return minFitIndex;
    }


    public double getDiversity(double[][] pop, double[] fitness_values){

        // Computes the standard deviation (sd) in every dimension
        // diversity is equal to the sum of the sds

        double diversity = 0.0;

        double standardDeviation=0.0;
        double sum = 0.0;

        for (int i = 0; i < PHENOTYPE_DIM; i++) {
            sum=0.0;
            standardDeviation=0.0;

            for (int n = 0; n < POP_SIZE; n++) {
                sum += pop[n][i];
            }
            double mean = sum/POP_SIZE;

            for(int n = 0; n < POP_SIZE; n++) {
                standardDeviation += Math.pow(pop[n][i] - mean, 2);
            }
            standardDeviation = standardDeviation/POP_SIZE;
        }
        diversity += standardDeviation;
        return diversity/PHENOTYPE_DIM;

    }
    private double new_mutation_factor(List<Integer> candidates, double[] fitness_values,int evals) {
        // Function to calculate the self-adaptive scaling factor
        if (evals + POP_SIZE <= (int)((double)evaluations_limit_*10.0)){
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

    private double[] get_mutant_f(double[][] pop, int base_indiv, List<Integer> candidates, double[][] scaling_factors, double[] fitness_values, int evals){

        // Initialize F and mutant scaling factor (mf)
        double NF = new_mutation_factor(candidates, fitness_values, evals);
        double[] mf = new double[PHENOTYPE_DIM];

        double[] difference = scaling_factors[candidates.get(0)].clone();
        for(int n = 0; n < PHENOTYPE_DIM; n++){
            for(int i=1; i<candidates.size(); i++){
                if(i < (candidates.size()/2)){
                    difference[n] += scaling_factors[candidates.get(i)][n];
                }else{
                    difference[n] -= scaling_factors[candidates.get(i)][n];
                }
            }
            double d = scaling_factors[base_indiv][n] + SCALING_FACTOR * difference[n];
            if ( d < 0){
                d = -1.0*d;
            }
            if(d > 1){
                d = d % 1.0;
            }

            mf[n] = Math.max(Math.min(d,1.0),0.1);
        }

        return mf;





    }


    private List<? extends Object> get_mutant_vector(double[][] pop, double[] fitness_values, double[][] scaling_factors, int evals){ // Efi
        // function to create a new mutant population based on pop
        // result dims should be [POP_SIZE][PHENOTYPE_DIM]
        // create placeholder for population
        double[][] mutants = new double[POP_SIZE][PHENOTYPE_DIM];
        double[][] mutants_f = new double[POP_SIZE][PHENOTYPE_DIM];

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
            double tau = rnd_.nextFloat();
            double [] mf = scaling_factors[j].clone();

            // if(tau<0.1) { // evolve scaling factors
            mf = get_mutant_f(pop, indexIndividual1, randomCandidateList, scaling_factors, fitness_values, evals);
            // }

            mutants_f[j] = mf.clone();


            // compute mutant input values
            double[] difference = pop[randomCandidateList.get(0)].clone();
            for(int n = 0; n < PHENOTYPE_DIM; n++){
                for(int i=1; i<randomCandidateList.size(); i++){
                    if(i < (randomCandidateList.size()/2)){
                        difference[n] += pop[randomCandidateList.get(i)][n];
                    }else{
                        difference[n] -= pop[randomCandidateList.get(i)][n];
                    }
                }
                // rescale to [-5,5] if necessary
                double d = individual1[n] + mf[n] * difference[n];
                if (d < -5){
                    d = -1.0*d;
                    d = d % 5;
                    d = -1.0*d;
                }
                if (d > 5){
                    d = d % 5;
                }
                mutants[j][n] = d;
            }

        }
        return Arrays.asList(mutants, mutants_f);
    }

    private List<? extends Object> get_trial_vector(double[][] pop, double[][] scaling_factors, double[][] mutants, double[][] mutants_f){ // Hans
        // function to create the trial vectors T
        // apply crossover based on CROSSOVER_SCHEME
        // return T
        double[][] T = new double[POP_SIZE][PHENOTYPE_DIM];
        double[][] Tf = new double[POP_SIZE][PHENOTYPE_DIM];
        for(int i = 0; i < POP_SIZE; i++) {
            List res = crossover(pop[i], scaling_factors[i], mutants[i], mutants_f[i]);
            T[i] = (double[])res.get(0);
            Tf[i] = (double[])res.get(1);
        }

        return Arrays.asList(T,Tf);
    }

    private List<? extends Object> crossover(double[] p1, double[] f1, double[] p2, double[] f2){
        // function to perform crossover between two parents, return new child
        double[] child = new double[PHENOTYPE_DIM];
        double[] child_f = new double[PHENOTYPE_DIM];
        int R = rnd_.nextInt(PHENOTYPE_DIM); //index on which at least crossover is applied

        if(CROSSOVER_SCHEME=="bin") {
            for (int i = 0; i < PHENOTYPE_DIM; i++) {
                double r_i = rnd_.nextDouble();
                if(r_i < CROSSOVER_PROB || i == R){
                    child[i] = p2[i];
                    child_f[i] = f2[i];
                }else{
                    child[i] = p1[i];
                    child_f[i] = f1[i];
                }
            }
        } else if(CROSSOVER_SCHEME=="exp") {
            int i = 0;
            while(rnd_.nextDouble()<CROSSOVER_PROB && i < PHENOTYPE_DIM) {
                child[i] = p2[i];
                child_f[i] = f2[i];
                i++;
            }
            for (int j = i; i < PHENOTYPE_DIM; i++) {
                child[j] = p1[j];
                child_f[j] =  f1[j];
            }
            child[R] = p2[R];
            child_f[R] = f2[R];

        }

        return Arrays.asList(child, child_f);
    }

    private List<? extends Object> survival_selection(double[][] parents, double[][] parents_f, double[] parents_fitness, double[][] children, double [][] children_f){ // Niels
        // function to select new surivors based on fitness
        double[][] survivors = new double[POP_SIZE][PHENOTYPE_DIM];
        double[][] survivors_f = new double[POP_SIZE][PHENOTYPE_DIM];
        double[] survivor_fitness = new double[POP_SIZE];
        // evaluate children
        double[] child_fitness = eval_pop(children);
        // compare each parent to child and save fittest in new population
        for(int i = 0; i < parents_fitness.length; i++){
            if(parents_fitness[i] > child_fitness[i]){
                survivors[i] = parents[i].clone();
                survivors_f[i] = parents_f[i].clone();
                survivor_fitness[i] = parents_fitness[i];
            }else{
                survivors[i] = children[i].clone();
                survivors_f[i] = children_f[i].clone();
                survivor_fitness[i] = child_fitness[i];
            }
        }
        // return list of objects, to be converted back outside function
        return Arrays.asList(survivors, survivor_fitness, survivors_f);
    }

    // MAIN FUNCTION
    public void run()
    {

        int evals = POP_SIZE;
        // init population
        double[][] pop = init_population(POP_SIZE, PHENOTYPE_DIM, DIM_LOWER_BOUND, DIM_UPPER_BOUND);
        double[][] scaling_factors = init_population(POP_SIZE, PHENOTYPE_DIM, 0.5,1.0);

        double diversity;
        // calculate fitness
        double[] fitness_scores = eval_pop(pop);
        while(evals<evaluations_limit_){


            diversity = getDiversity(pop, fitness_scores);
            //System.out.println(fitness_scores[getFittest(fitness_scores)]);
            System.out.println(diversity);

            // Select parents
            // Apply crossover / mutation operators
            List m = get_mutant_vector(pop, fitness_scores, scaling_factors, evals);
            double[][] mutant_pop = (double[][])m.get(0);
            double[][] mutants_f = (double[][])m.get(1);


            List T = get_trial_vector(pop, scaling_factors, mutant_pop, mutants_f);
            double[][] trial_pop = (double[][])T.get(0);
            double[][] trial_f = (double[][])T.get(1);
            // Selection
            List res = survival_selection(pop, scaling_factors, fitness_scores, trial_pop, trial_f);
            pop = (double[][])res.get(0);
            fitness_scores = (double[])res.get(1);
            scaling_factors = (double[][])res.get(2);
            evals += POP_SIZE;

        }
    }

    public void main(String[] argv){
        run();
    }
}