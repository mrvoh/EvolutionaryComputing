
from subprocess import run
from subprocess import PIPE


out = run('java -jar testrun.jar -submission=player26 -evaluation=SphereEvaluation -seed=1', shell=True, stdout=PIPE)

print('hier: {}'.format(out.stdout))
