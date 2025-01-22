import sys

from importlib.resources import files
def path(package, file):
    return files(package).joinpath(file)

acps = []
non_acps = []
cpps = []
non_cpps = []

with open(path('accpps.data', 'acps_positive.txt')) as fff:
    acps = fff.read().splitlines()

with open(path('accpps.data', 'acps_negative.txt')) as fff:
    non_acps = fff.read().splitlines()

with open(path('accpps.data', 'cpps_positive.txt')) as fff:
    cpps = fff.read().splitlines()

with open(path('accpps.data', 'cpps_negative.txt')) as fff:
    non_cpps = fff.read().splitlines()

