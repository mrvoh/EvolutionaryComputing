
from subprocess import run
from subprocess import PIPE


# out = run('java -jar testrun.jar -submission=player26 -evaluation=SphereEvaluation -seed=1', shell=True, stdout=PIPE)

# print('hier: {}'.format(out.stdout))


def change_params(filename, pop_size = 50, scaling_factor = 0.1, crossover_prob = 0.2, crossover_rate = 0.3, nr_perturbation_vectors = 1, base_vector = '"random1"', crossover_scheme='"bin1"'):


    vals = [pop_size, scaling_factor, crossover_prob, crossover_rate, nr_perturbation_vectors, base_vector, crossover_scheme]
    result = ''
    with open(filename, 'r') as f:

        param_decl = ['public int POP_SIZE','public double SCALING_FACTOR','public double CROSSOVER_PROB',  
                        'public double CROSSOVER_RATE', 'public int NR_PERTURBATION_VECTORS',
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


    with open('test.java', 'w') as f:
        f.write(result)

params = {
    'pop_size' : 50000, 
    'scaling_factor' : 0.111, 
    'crossover_prob' : 0.2, 
    'crossover_rate' : 0.3, 
    'nr_perturbation_vectors' : 1, 
    'base_vector' : '"random1rthtyrh"', 
    'crossover_scheme':'"bin1"'
}

change_params('DifferentialEvolution.java', **params)