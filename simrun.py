# Python script for cleaning up data folders, generating input files, and
# running a simulation

import csv
import numpy as np
import os

"""Function for writing input file
Inputs:
outfile - file where parameter file is to be exported
templatefile - file that will be used as a template
               any parameters unspecified will be kept
               from the template
params - dictionary of param name strings and associated values
         (use consistent caps to be consistent with template)
"""


def write_input(outfile, templatefile, params):

    # Define space delimited dialect
    csv.register_dialect('mydialect', delimiter=' ')
    # Define space delimited dialect
    param_dict = {}
    # Open the template file and read in as a dictionary of
    # paramater name, value pairs
    with open(templatefile, 'r') as template:
        reader = csv.reader(filter(lambda row: row[0] != '#', template),
                            dialect='mydialect')

        # Loop over paramname,value pairs in template
        for param, value in reader:
            param_dict[param] = value
    # Update param_dict with the new params
    for param in params:
        param_dict[param] = params[param]
    # Open the target output file for writing
    with open(outfile, 'w') as inputfile:
        writer = csv.writer(inputfile, dialect='mydialect')
        # Loop over the param dictionary and save as key,value pairs
        for param in param_dict:
            writer.writerow([param, param_dict[param]])


# ################################################
# Set up input files and data folders
# ################################################


# Set the LKs for the simulation
LKs = [0, 1, 2, 3]
np.savetxt('input/LKs', LKs, fmt='%i')


# Clean up data folder
os.system('mv data/* trash/')

# Generate data folders for each LK
for lk in LKs:
    os.system('mkdir data/LK_%s' % lk)

NINIT = 100
L = 1000.
NB = 100
NSTEP = 100
NREPLICAEXCHANGE = 100
INDMAX = 100
PTON = 'T'
LK = 0
params = {'NINIT': NINIT, 'L': L, 'INDMAX': INDMAX, 'NSTEP': NSTEP,
          'PTON': PTON, 'LK': LK, 'NREPLICAEXCHANGE': NREPLICAEXCHANGE,
          'NB': NB}

write_input('input/input', 'input/template', params)

# #################################################
# Compile and run the simulation
# #################################################

# compile the simulation
os.system('./compile.sh')

# run the simulation with the appropriate number of processors
os.system('mpirun -np %i ./wlcsim' % (len(LKs)+1))
