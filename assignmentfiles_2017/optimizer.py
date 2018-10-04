from bayes_opt import BayesianOptimization
import os
import sys
import pandas as pd
from io import StringIO
import numpy as np
from subprocess import run
from subprocess import PIPE

class Optimizer:

    def __init__(self, optimization_folder, nr_iterations, iteration_chunck_size, nr_init_points, nr_folds, algo_name, obj_fun):

        # Set static variables
        self.INTERMEDIATE_RESULTS_FOLDER = optimization_folder
        self.FINAL_RESULTS_FOLDER = optimization_folder
        self.NR_ITERATIONS = nr_iterations
        self.ITERATION_CHUNCK_SIZE = iteration_chunck_size
        self.NR_INIT_POINTS = nr_init_points
        self.NR_OF_FOLDS = nr_folds

        self.OBJ_FUN = obj_fun
        self.ALGO_NAME = algo_name



        # Boundaries between which to explore the input space
        self.param_boundaries = {
            'pop_size' : (20,500), 
            'scaling_factor' : (0,1), 
            'crossover_prob' : (0,1),  
            'nr_perturbation_vectors' : (1,3), 
            'base_vector' : (0,1), # rand / best
            'crossover_scheme':(0,1) # binomial / exponential
        }
        # Set points on which to evaluate the model for exploration of the solution space
        self.explore_points = {
            'pop_size' : [25, 100, 500], 
            'scaling_factor' : [0.5, 0.9, 0.1], 
            'crossover_prob' : [0.3, 0.05, 0.8],  
            'nr_perturbation_vectors' : [1,2,3], 
            'base_vector' : [0,0,1], # rand / best
            'crossover_scheme':[1,0,0] # binomial / exponential
        }



        self.bo = BayesianOptimization(self._eval_fun, self.param_boundaries)

    ######################################################
    # Helper functions
    ######################################################

    def _change_params(self, filename, pop_size = 50, scaling_factor = 0.1, crossover_prob = 0.2, nr_perturbation_vectors = 1, base_vector = '"random1"', crossover_scheme='"bin1"'):


        vals = [pop_size, scaling_factor, crossover_prob, nr_perturbation_vectors, base_vector, crossover_scheme]
        result = ''
        with open(filename, 'r') as f:

            param_decl = ['public int POP_SIZE','public double SCALING_FACTOR','public double CROSSOVER_PROB',  
                            'public int NR_PERTURBATION_VECTORS',
                            'public String BASE_VECTOR', 'public String CROSSOVER_SCHEME']

            param_val = list(zip(param_decl, vals))
            for line in f.readlines():
                l = ''
                for decl, val in param_val:
                    if decl in line: # replace param values
                        l = '\t{} = {};\n'.format(decl, val)



                if l == '': #no param in current line
                    l = line

                result += l


        with open(filename, 'w') as f:
            f.write(result)

    def _cast_params(self, params):

        # convert integers
        for param in ['pop_size', 'nr_perturbation_vectors']:
            params[param] = int(round(params[param],0))

        # cast categorical variablest
        params['base_vector'] = '"best"' if params['base_vector'] > 0.5 else '"rand"'
        params['crossover_scheme'] = '"bin"' if params['crossover_scheme'] > 0.5 else '"exp"'


        return params

    def _eval_fun(self, pop_size, scaling_factor, crossover_prob, nr_perturbation_vectors, base_vector, crossover_scheme):

        
        #TODO: Include default settings
        # params to evaluate:
        params = {
            'pop_size' : pop_size, 
            'scaling_factor' : scaling_factor, 
            'crossover_prob' : crossover_prob,  
            'nr_perturbation_vectors' : nr_perturbation_vectors, 
            'base_vector' : base_vector, # rand / best
            'crossover_scheme':crossover_scheme # binomial / exponential
        }

        params = self._cast_params(params)

        score = 0
        for i in range(self.NR_OF_FOLDS):
            # change params
            self._change_params(self.ALGO_NAME, **params)
            # compile algo
            run('javac -cp contest.jar {}'.format(self.ALGO_NAME), shell=True)

            # run algorithm and capture output
            submission = self.ALGO_NAME.replace('.java','')
            seed = 1234567*i
            out = run('java -jar testrun.jar -submission={} -evaluation={} -seed={}'.format(submission, self.OBJ_FUN, seed), shell=True, stdout=PIPE)
            # parse output
            out = str(out)
            out = out[out.find('Score:'):]

            score += float(out.split('\\n')[0][7:])

        #with Capturing() as output:
        return score


    #############################################################################
    # Public functions
    #############################################################################

    def init_csv(self, path):
        """
            Initializes KNOWN target values stored in *.csv format to optimizer object
        """

        df = pd.read_csv(path, header = 0)
        df = df.rename(columns=lambda x: x.strip('#')) # strip hashtag since it is automatically added when saving results
        df = df.rename(columns=lambda x: x.strip())
        self.bo.initialize_df(df)

    def get_best_params(self):
        """
            Return the best found parameter settings after optimizing in dictionary form
        """
        assert self.bo is not None, "Parameters can only be retrieved AFTER optimizing"

        print(self.bo.res['max'])

        params = self.bo.res['max']['max_params']
        params = self._cast_params(params)

        return params

    

    def optimize(self):
        """
        Main function for optimization
        """

        # Explore the input and target space on predefined points
        self.bo.explore(self.explore_points)
        
        
        # Set parameters for Gaussian Process
        #TODO: read in into GPs
        gp_params = {}
        #{'kernel': None, 'alpha': 1e-5}

        # Iteratively maximize, store values, reinitialize
        for i in range(0, self.NR_ITERATIONS, self.ITERATION_CHUNCK_SIZE):
            # Only initialize extra random points in first epoch
            nr_init_points = self.NR_INIT_POINTS if i ==0 else 1

            # Perform maximization
            self.bo.maximize(init_points = nr_init_points, n_iter=self.ITERATION_CHUNCK_SIZE, acq='ei', **gp_params)

            # Store intermediate results
            file_path = os.path.join(self.INTERMEDIATE_RESULTS_FOLDER, '{}-intermediate.csv'.format(i))
            self.bo.points_to_csv(file_path)

        # Store and upload results
        filepath = os.path.join(self.FINAL_RESULTS_FOLDER, 'final_optimization_results.csv')
        self.bo.points_to_csv(filepath)


Opt = Optimizer('optimizer_results', 250, 25, 25, 5, 'DENF.java', 'KatsuuraEvaluation')

# Reinitialize from *.csv file in case of interruption. IMPORTANT: also disable line line "self.bo.explore(self.explore_points)"
#Opt.init_csv('optimizer_results/0-intermediate.csv')

Opt.optimize()