import java.util.Random;
import java.util.Arrays;
import java.util.*;
class test{

    static Random rnd_;
    // Static for assignment
    public static int PHENOTYPE_DIM = 2; // each phenotype has 10 dimensions in our assignment
    public int DIM_LOWER_BOUND = -5; // each dimension ranges from [-5, 5]
    public int DIM_UPPER_BOUND = 5;

    // Changeable params
    public static int POP_SIZE = 2; // mu
    public double SCALING_FACTOR = 0.5; //F
    public double CROSSOVER_PROB = 0.5; // Cr
    public double CROSSOVER_RATE = 0.1; // Pc

    // Params for DE operators (different versions of algorithm)
    public int NR_PERTURBATION_VECTORS = 1;
    public String BASE_VECTOR = "random";
    public String CROSSOVER_SCHEME = "bin";

    private static double[][] init_population(int pop_size, int phenotype_dim){
            // function to intialize a random population of size pop_size

            double[][] pop = new double[pop_size][phenotype_dim];

            for(int i = 0; i < pop_size;i++){
                // init phenotype
                double[] pheno = new double[phenotype_dim];

                for(int j = 0; j < phenotype_dim; j++){
                    // initialize values randomly
                    pheno[j] = rnd_.nextDouble();
                    }

                pop[i] = pheno;

            }

            return pop;
    }

    private static double[] eval_pop(double[][] pop){ //Ronald
        // function to evaluate each phenotype in a population

        double[] dummy = new double[POP_SIZE];
        return  dummy;
    }


    private static List<Object> survival_selection(double[][] parents, double[] parents_fitness, double[][] children){ // Niels
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

    public static void run(String[] argv){
        rnd_ = new Random(); // for testing purposes
        rnd_.setSeed(1);
        double[][] pop = init_population(2, 2);
        // System.out.println(Arrays.deepToString(pop));

        // System.out.println("TEST");

        // double[][] new_pop = init_population(3, 10);
        // System.out.println(Arrays.deepToString(new_pop));
        double[] fitness = {0.1, 0.1}
        System.out.println(Arrays.toString(fitness));

        List res = survival_selection(pop, fitness, pop);
        double[][] survivors = (double[][])res.get(0);
        double[] surv_fitness = (double[])res.get(1);

        System.out.println(Arrays.deepToString(survivors));

        System.out.println(Arrays.toString(surv_fitness));

    }


    public static void main(String[] args){
        run(args);
    }
}