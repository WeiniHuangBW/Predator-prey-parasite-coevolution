##### THIS SCRIPT DOES STOCHASTIC SIMULATIONS USING THE GILLESPIE ALGORITHM #####
##### THIS SCRIPT RECORDS THE COEXISTENCE OF PREDATOR, PREY AND PARASITES IN A PREDATOR-PREY-PARASITE SYSTEM WHERE THE PARASITE IS TRANSMITTED TROPHICALLY FROM INFECTED PREY TO PREDATORS #####

## Import relevant packages and modules ##
import pandas as pd, math, statistics, random, numpy.random as nura, numpy as np, array as arr, matplotlib.pyplot as plt, matplotlib.patches as mpatches, sys, getopt, time

# Create .csv output file for recording the average abundance, genotypic diversity, genotypic lifespan, emergence of super-types and coexistence
# str(float(sys.argv[4])) position 4 in command line to assign a number of the corresponding realisation to avoid overwriting the correspondent files
output_file = open("output_file_" + str(float(sys.argv[4])) + ".csv", "w") # This file will provide information about average abundance, average genotypes, emergence of super-types and coexistence from a single realisation
output_file.write("rp,rx,prey individuals,parasite individuals,predator individuals,prey genotypes,parasite genotypes,predator genotypes,prey genotypic lifespan,parasite genotypic lifespan,predator genotypic lifespan,super resistant prey,super infective parasite,super resistant predator,all coexist,predator and prey coexist,prey coexist,all extinct\n")

## Function for appending rows of new genotypes(for simulation purposes)
def loop_to_compare_array(the_array, new_genotype, time):
    for jj in range(0, the_array.shape[0]): # Go through all the rows in array_prey
        if(np.array_equal(new_genotype, the_array[jj, 2])):
            the_array[jj, 1] += 1 # if the new genotype is the same as an existing one (dead or alive) sum one individual
            return the_array
            break
        if jj == the_array.shape[0] - 1:
            the_array = np.append(the_array,[[jj+1, 1, new_genotype]], axis=0) # if the new genotype is different than any of the existing ones (dead or alive) add another row with one individual of that genotype
            return the_array

## Function for appending new genotypes (for storage purposes)
def loop_to_store_array(the_array, new_genotype, time):
    for jj in range(0, the_array.shape[0]): # Go through all the rows in array_prey
        if(np.array_equal(new_genotype, the_array[jj, 2])):
            the_array[jj, 1] += 1 # if the new genotype is the same as an existing one (dead or alive) sum one individual
            return the_array
            break
        if jj == the_array.shape[0] - 1:
            the_array = np.append(the_array,[[jj+1, 1, new_genotype, 0, time]], axis=0) # if the new genotype is different than any of the existing ones (dead or alive) add another row with one individual of that genotype
            return the_array

## Function to remove individuals (storage)
def loop_to_remove_array(the_array, new_genotype):
    for jj in range(0, the_array.shape[0]): # Go through all the rows in array_prey
        if(np.array_equal(new_genotype, the_array[jj, 2]) and (the_array[jj, 1] != 0)):
            the_array[jj, 1] -= 1 # if the new genotype is the same as an existing one (dead or alive) sum one individual
            return the_array
            break

## Function for checking mutations in all loci
def my_mutation(x,mu): # give random number (between 0 and 1) and mutation rate
    return x < mu # condition for using in the function below (the mutation will only happen if the random number is lower that mutation rate) 

def my_mutation_loci(n_loci, mu, initial_array): # give number of loci in genotype, mutation rate, and genotype that is reproducing and may mutate
    mutation_loci = np.zeros(n_loci) # all the posible positions (loci) in which the genotype can mutate
    sum_mutation = 0
    for ii in range(0, n_loci): # check all loci (positions in genotype)
        mutation_loci[ii] = random.uniform(0,1) # throw a random number (between 0 and 1) for all the positions in genotype
        sum_mutation = sum(1 for x in mutation_loci if my_mutation(x,mu)) # sum all the positions of the genotype were the random number was lower than the mutation rate
    temporal_array = np.array(initial_array) # asign a temporal name to the genotype that may (or may not) mutate
    if(sum_mutation != 0): # if any of the positions had a random number lower that the mutation rate, then mutation will happen
        mut = np.random.choice(n_loci, sum_mutation, replace = False) # pick randomly those loci that will mutate, where sum_mutation is the number of loci mutating, replace = False to avoid using the same locus twice.
        for ii in range(0, sum_mutation): # for all the loci mutating check whether it is a 0 or a 1
            if(temporal_array[mut[ii]] == 0): # if it is a 0, change it to 1
                temporal_array[mut[ii]] = 1
            else: # if it is a 1, change it to 0
                temporal_array[mut[ii]] = 0
    return temporal_array

## Function to calculate the probability of successful infection of all the prey/parasite combinations (prey-parasite interactions)
def infection_of_prey(genotype_prey, genotype_parasite):
    part1 = genotype_prey[genotype_parasite == 0] # interaction of current prey and parasite genotypes
    part2 = part1[part1 == 1] # which loci of the prey genotype are resistant to the parasite genotype (1 in the prey matching the same position as a 0 in the parasite)
    r = part2.sum() # some all those 1 in prey matching 0 in parasite (calculate resistance)
    Qx = float(sigma_value_prey)**float(r) # calculate infection rate according to resistance
    return Qx

## Function to calculate the probability of successful infection of all the predator/parasite combinations (predator-parasite interactions)
def infection_of_predator(genotype_predator, genotype_parasite):
    part1 = genotype_predator[genotype_parasite == 0] # interaction of current prey and parasite genotypes
    part2 = part1[part1 == 1] # which loci of the predator genotype are resistant to the parasite genotype (1 in the predator matching the same position as a 0 in the parasite)
    r = part2.sum()  # some all those 1 in predator matching 0 in parasite (calculate resistance)
    Qy = float(sigma_value_predator)**float(r) # calculate infection rate according to resistance
    return Qy

## Function to append new infected prey
def loop_infection(array_infected, host_genotype, parasite_genotype):
    for ii in range(0, array_infected.shape[0]): # Go through all the rows in array of infected prey
        if(np.array_equal(array_infected[ii, 2], host_genotype) and np.array_equal(array_infected[ii, 3], parasite_genotype)):
            array_infected[ii, 1] += 1 # if we already have that prey genotype infected by that parasite genotype then just add one individual to that row
            return array_infected
            break
        if ii == array_infected.shape[0] - 1: # else add the new combination of interaction by adding a row with that prey and parasite genotype (one individual)
            array_infected = np.append(array_infected,[[ii+1, 1, host_genotype, parasite_genotype]], axis=0)
            return array_infected

### PARAMETERS ###
# initial population sizes
pop_uninfected_prey = 800 # uninfected prey
pop_parasite = 1000 # free-living parasites
pop_uninfected_predator = 100 # uninfected predator
pop_infected_prey = 0 # at the beginning of the simulation, all prey are uninfected
pop_infected_predator = 0 # at the beginning of the simulation, all predator are uninfected

# prey parameters 
gx = 2.0 # growth rate of the prey
rx = float(sys.argv[1]) # first position in command line to vary the parameter values of the reproduction cost of the prey due to parasite infection
dx = 0.1 # intrinsic death rate of the prey
pop_limit = 2000 # population limit (competition between prey for carrying capacity)

# parasite parameters
n_z = 6 # number of parasite offspring per reproduction
dz = 0.09 # intrinsic death rate of the free-living parasite
S = 0.0005 # scaling factor for prey-parasite population sizes

#predator parameters
fy = 0.01 # predation rate
ky = 0.2 # reproduction rate of the predator
re = 1.0 # reproduction cost of the predator due to parasite exposure
rp = float(sys.argv[2]) # second position in command line to vary the parameter values of the reproduction cost of the predator due to parasite infection
dy = 1.0 # intrinsic death rate of the predator

# probability of successful infection
sigma_value_predator = 0.85 # to calculate probability of infection according to predator resistance
sigma_value_prey = 0.85 # to calculate probability of infection according to prey resistance
n_loci = 10 # total number of loci in genotypes
mx = 0.00002 # mutation rate prey (probability of mutation per locus per reproduction event)
mz = 0.000006 # mutation rate parasite (probability of mutation per locus per reproduction event)
my = 0.00005 # mutation rate predator (probability of mutation per locus per reproduction event)

# Genotypes length
length_genotype_prey = np.zeros(n_loci) # initial prey genotype (all loci are initially susceptible, i.e. zeros)
length_genotype_parasite = np.zeros(n_loci) # initial parasite genotype (all loci are initially non-infective, i.e. zeros)
length_genotype_predator = np.zeros(n_loci) # initial predator genotype (all loci are initially susceptible, i.e. zeros)

# POPULATION ARRAYS
np.warnings.filterwarnings('ignore',category=np.VisibleDeprecationWarning) # to avoid warning messages that result due to python version conflict
# Store information of uninfected predator, infected prey and free-living parasites
# Position 0 corresponds to the index. Position 1 corresponds to the number of individuals. Position 2 corresponds to the host genotype (string of binary numbers).
prey = np.array([[0,pop_uninfected_prey,length_genotype_prey]], dtype = object)
parasite = np.array([[0,pop_parasite,length_genotype_parasite]], dtype = object)
predator = np.array([[0,pop_uninfected_predator,length_genotype_predator]], dtype = object)

# Store information of infected predator and prey
# Position 0 corresponds to the index. Position 1 corresponds to the number of individuals. Position 2 corresponds to the host genotype (string of binary numbers). Position 3 corresponds to the parasite genotype (string of binary numbers).
infected_prey = np.array([[0,pop_infected_prey,length_genotype_prey,length_genotype_parasite]], dtype = object)
infected_predator = np.array([[0,pop_infected_predator,length_genotype_predator,length_genotype_parasite]], dtype = object)

# store all information of predators and prey (infected and uninfected), and free-living parasites
# Position 0 corresponds to the index. Position 1 corresponds to the number of individuals. Position 2 corresponds to the genotype (string of binary numbers). 
# Position 3 corresponds to whether the genotype reached at least 2% of the total prey population size (0 means "no" and 1 means "yes"). 
# Position 4 corresponds to the time in which a new genotype emerged (for calculating genotypic lifespan)
store_prey = np.array([[0,pop_uninfected_prey,length_genotype_prey,1,0]], dtype = object)
store_parasite = np.array([[0,pop_parasite,length_genotype_parasite,1,0]], dtype = object)
store_predator = np.array([[0,pop_uninfected_predator,length_genotype_predator,1,0]], dtype = object)

# Variables for recording time of emergence of all-resistant hosts and all-infective parasites (simulation purposes)
emergence_st_prey = False
emergence_st_predator = False
emergence_st_parasite = False

# Variables to record extinctions (simulation puporses)
extinction_prey = False
extinction_predator = False
extinction_parasite = False

# Assign extinctions (output purposes; this is done at the end of the script after the algorithm finishes)
coexistence = "0"
coexistence_predator_and_prey = "0"
coexistence_prey = "0"
extinction = "0"

## if there was no emergence of super-types, then we assign "0"
emergence_prey = "0"
emergence_parasite = "0"
emergence_predator = "0"

# Empty arrays to store the abundance of each entity in each time step
store_sum_prey = []
store_sum_parasite = []
store_sum_predator = []

# Empty arrays to store the number of effective genotypes of each entity in each time step
store_prey_genotype = []
store_parasite_genotype = []
store_predator_genotype = []

# Empty arrays to store all genotypic lifespans
store_prey_genotypic_lifespan = []
store_parasite_genotypic_lifespan = []
store_predator_genotypic_lifespan = []

# Variables to record average abundances
av_prey_individuals = 0
av_parasite_individuals = 0
av_predator_individuals = 0

# Variables to record average genotypes
av_prey_genotypes = 0
av_parasite_genotypes = 0
av_predator_genotypes = 0

# Variables to record average genotypic lifespans
av_prey_genotypic_lifespan = 0
av_parasite_genotypic_lifespan = 0
av_predator_genotypic_lifespan = 0


# Time variables
max_time = float(sys.argv[3]) # third position in command line to vary the total time to run the simulation (Gillespie algorithm)
Time = 0 # initial Gillespie time (we add dt_next_event and run until reaching max_time)
dt_next_event = 0 # random time step after event occurs (following the Gillespie algorithm). This quantity is summed to the total time (continuos time simulation)
n = 0 # number of time points across simulations in which we record abundances of subpopulations (in units of one, i.e., n+1)

### GILLESPIE ALGORITHM ###
### SIMULATION STARTS ###
while Time < max_time: # SIMULATION STARTS: repeat simulation until reaching max time

    # Optimize simulation arrays (remove extinct genotypes for speeding simulations):
    if(prey.shape[0] != 1):
        prey = prey[prey[:,1] != 0]
    if(infected_prey.shape[0] != 1):
        infected_prey = infected_prey[infected_prey[:,1] != 0]
    if(parasite.shape[0] != 1):
        parasite = parasite[parasite[:,1] != 0]
    if(predator.shape[0] != 1):
        predator = predator[predator[:,1] != 0]
    if(infected_predator.shape[0] != 1):
        infected_predator = infected_predator[infected_predator[:,1] != 0]
    
    # Optimize storing arrays (keep only those that are still alive or were effective)
    if(store_prey.shape[0] != 1):
        store_prey = store_prey[np.logical_not(np.logical_and(store_prey[:,1] == 0, store_prey[:,3] == 0))]
    if(store_parasite.shape[0] != 1):
        store_parasite = store_parasite[np.logical_not(np.logical_and(store_parasite[:,1] == 0, store_parasite[:,3] == 0))]
    if(store_predator.shape[0] != 1):
        store_predator = store_predator[np.logical_not(np.logical_and(store_predator[:,1] == 0, store_predator[:,3] == 0))]

## STEP 1 ##
## CALCULATE VALUES FOR EVERY POSSIBLE EVENT ACCORDING TO CURRENT CONDITIONS ##      
###### events uninfected prey ######
    prey_growth = prey[:,1] * gx # uninfected prey reproduces
    prey_death = prey[:,1] * dx # uninfected prey dies
    prey_competition = prey[:,1] * (sum(prey[:,1]) + sum(infected_prey[:,1])) * (1 /pop_limit) # uninfected prey dies due to competition

###### events infected prey ######
    infected_prey_growth = infected_prey[:,1] * gx * rx # infected prey reproduces
    infected_prey_death = infected_prey[:,1] * dx # infected prey dies
    infected_prey_competition = infected_prey[:,1] * (sum(infected_prey[:,1]) + sum(prey[:,1])) * (1 /pop_limit) # infected prey dies due to competition

###### events free-living parasite ######
    infection_prey = np.zeros([parasite.shape[0],prey.shape[0]], dtype = float) # storage array for event
    non_infection_prey = np.zeros([parasite.shape[0],prey.shape[0]], dtype = float) # storage array for event
    for i in range(0,parasite.shape[0]):
        for j in range(0,prey.shape[0]):
            infection_prey[i,j] = parasite[i,1] * prey[j,1] * infection_of_prey(prey[j,2],parasite[i,2]) * S # free-living parasite infects prey
            non_infection_prey[i,j] = parasite[i,1] * prey[j,1] * (1-infection_of_prey(prey[j,2],parasite[i,2])) * S # free-living parasite fails infecting prey
    parasite_death = parasite[:,1] * dz # free-living parasite dies

###### events uninfected predator ######
    predator_growth = np.zeros([predator.shape[0],prey.shape[0]], dtype = float) # storage array for event
    predator_non_growth = np.zeros([predator.shape[0],prey.shape[0]], dtype = float) # storage array for event
    for i in range(0,predator.shape[0]):
        for j in range(0,prey.shape[0]):
            predator_growth[i,j] = predator[i,1] * prey[j,1] * fy * ky # uninfected predator reproduces after feeding
            predator_non_growth[i,j] = predator[i,1] * prey[j,1] * fy * (1-ky) # uninfected predator does not reproduce after feeding

    predator_exposure_growth = np.zeros([predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    predator_exposure_non_growth = np.zeros([predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    predator_infection_growth = np.zeros([predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    predator_infection_non_growth = np.zeros([predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    for i in range(0,predator.shape[0]):
        for j in range(0,infected_prey.shape[0]):
            predator_exposure_growth[i,j] = predator[i,1] * infected_prey[j,1] * fy * (1-infection_of_predator(predator[i,2],infected_prey[j,3])) * re * ky # uninfected predator exposed to parasite reproduces
            predator_exposure_non_growth[i,j] = predator[i,1] * infected_prey[j,1] * fy * (1-infection_of_predator(predator[i,2],infected_prey[j,3])) * (1 - (re * ky)) # uninfected predator exposed to parasite reproduces
            predator_infection_growth[i,j] = predator[i,1] * infected_prey[j,1] * fy * infection_of_predator(predator[i,2],infected_prey[j,3]) * rp * ky # uninfected predator infected by parasite reproduces
            predator_infection_non_growth[i,j] = predator[i,1] * infected_prey[j,1] * fy * infection_of_predator(predator[i,2],infected_prey[j,3]) * (1 - (rp * ky)) # uninfected predator infected by parasite does not reproduce
    predator_death = predator[:,1] * dy # predator intrinsic death
    
    # events infected predator
    infected_predator_growth = np.zeros([infected_predator.shape[0],prey.shape[0]], dtype = float) # storage array for event
    infected_predator_non_growth = np.zeros([infected_predator.shape[0],prey.shape[0]], dtype = float) # storage array for event
    for i in range(0,infected_predator.shape[0]):
        for j in range(0,prey.shape[0]):
            infected_predator_growth[i,j] = infected_predator[i,1] * prey[j,1] * fy * rp * ky # infected predator reproduces after feeding
            infected_predator_non_growth[i,j] = infected_predator[i,1] * prey[j,1] * fy * (1 - (rp * ky)) # infected predator does not reproduce after feeding
    
    infected_predator_exposure_growth = np.zeros([infected_predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    infected_predator_exposure_non_growth = np.zeros([infected_predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    infected_predator_infection_growth = np.zeros([infected_predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    infected_predator_infection_non_growth = np.zeros([infected_predator.shape[0],infected_prey.shape[0]], dtype = float) # storage array for event
    for i in range(0,infected_predator.shape[0]):
        for j in range(0,infected_prey.shape[0]):
            infected_predator_exposure_growth[i,j] = infected_predator[i,1] * infected_prey[j,1] * fy * (1-infection_of_predator(infected_predator[i,2],infected_prey[j,3])) * re * rp * ky # infected predator exposed to the parasite reproduces
            infected_predator_exposure_non_growth[i,j] = infected_predator[i,1] * infected_prey[j,1] * fy * (1-infection_of_predator(infected_predator[i,2],infected_prey[j,3])) * (1 - (re * rp * ky))  # infected predator exposed to the parasite does not reproduce
            infected_predator_infection_growth[i,j] = infected_predator[i,1] * infected_prey[j,1] * fy * infection_of_predator(infected_predator[i,2],infected_prey[j,3]) * rp * rp * ky # infected predator infected by parasite reproduces
            infected_predator_infection_non_growth[i,j] = infected_predator[i,1] * infected_prey[j,1] * fy * infection_of_predator(infected_predator[i,2],infected_prey[j,3]) * (1 - (rp * rp * ky)) # infected predator infected by parasite does not reproduce
    infected_predator_death = infected_predator[:,1] * dy # infected predator intrinsic death
    
    # Store the sum of all the events
    sum_events = (prey_growth.sum() + prey_death.sum() + prey_competition.sum() + 
    infected_prey_growth.sum() + infected_prey_death.sum() + infected_prey_competition.sum() + 
    infection_prey.sum() + non_infection_prey.sum() + parasite_death.sum() + 
    predator_growth.sum() + predator_non_growth.sum() + 
    predator_exposure_growth.sum() + predator_exposure_non_growth.sum() + 
    predator_infection_growth.sum() + predator_infection_non_growth.sum() +
    predator_death.sum() + 
    infected_predator_growth.sum() + infected_predator_non_growth.sum() + 
    infected_predator_exposure_growth.sum() + infected_predator_exposure_non_growth.sum()  + 
    infected_predator_infection_growth.sum() + infected_predator_infection_non_growth.sum() + 
    infected_predator_death.sum())

    ## STEP 2 ##
    ## CALCULATE NEXT TIME STEP AND NEXT EVENT ##
    ## GENERATE RANDOM NUMBERS ##
    
    # Next time step
    dt_next_event = np.random.exponential(scale=1/sum_events)

    # Next event
    URN = random.uniform(0,1) # unit-interval uniform random number generator for next event
    P = 0 # for doing cummulative sum in picking the next event

    ## STEP 3 ##
    ## EVENT HAPPENS, UPDATE POPULATION SIZES AND ADD TIME STEP TO TOTAL TIME ##
    
    ############### Uninfected prey ####################
    occurrence = False
    while not occurrence:
        for i in range(0,prey.shape[0]):
            if URN > P and URN <= P + prey_growth[i]/sum_events:
                bx = 1 # uninfected prey increases by one
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbx = i # row number in uninfected prey array
                occurrence = True
                break
            P += prey_growth[i]/sum_events #! use += to modify in place
            
            if URN > P and URN <= P + prey_death[i]/sum_events:
                bx = 0 # uninfected prey decreases by one
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbx = i # row number in uninfected prey array
                occurrence = True
                break
            P += prey_death[i]/sum_events
        
            if URN > P and URN <= P + prey_competition[i]/sum_events:
                bx = 0 # uninfected prey decreases by one
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbx = i # row number in uninfected prey array
                occurrence = True
                break
            P += prey_competition[i]/sum_events
        
        if occurrence:
            break

        ############### Infected prey #################
        for i in range(0,infected_prey.shape[0]):
            if URN > P and URN <= P + infected_prey_growth[i]/sum_events:
                bx = 6 # infected prey reproduces
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbi = i # row number in infected prey array
                occurrence = True
                break
            P += infected_prey_growth[i]/sum_events

            if URN > P and URN <= P + infected_prey_death[i]/sum_events:
                bx = 5 # infected prey decreases by one
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbi = i # row number in infected prey array
                occurrence = True
                break
            P += infected_prey_death[i]/sum_events

            if URN > P and URN <= P + infected_prey_competition[i]/sum_events:
                bx = 5 # infected prey decreases by one
                bz = 2 # nothing happens to free-living parasite
                by = 2 # nothing happens to predator
                mbi = i # row number in infected prey array
                occurrence = True
                break
            P += infected_prey_competition[i]/sum_events
        
        if occurrence:
            break
                                            
        ################ Free-living parasite ####################
        for i in range(0,parasite.shape[0]):
            for j in range(0,prey.shape[0]):
                if URN > P and URN <= P + infection_prey[i,j]/sum_events:
                    bx = 3 # prey is carrying a parasite (now it is infected)
                    bz = 0 # parasite gets in the prey
                    by = 2 # nothing happens to predator
                    mbx = j # row number in uninfected prey array
                    mbz = i # row number in free-living parasite array
                    occurrence = True
                    break
                P += infection_prey[i,j]/sum_events
            
                if URN > P and URN <= P + non_infection_prey[i,j]/sum_events:
                    bx = 2 # nothing happens to prey (it is not infected)
                    bz = 2 # nothing happens to free-living parasite
                    by = 2 # nothing happens to predator
                    mbx = j # row number in uninfected prey array
                    mbz = i # row number in free-living parasite array
                    occurrence = True
                    break
                P += non_infection_prey[i,j]/sum_events

        if occurrence:
            break

        for i in range(0,parasite.shape[0]):
            if URN > P and URN <= P + parasite_death[i]/sum_events:
                bx = 2 # nothing happens to prey
                bz = 0 # free-living parasite decreases by one
                by = 2 # nothing happens to predator
                mbz = i # row number in free-living parasite array
                occurrence = True
                break
            P += parasite_death[i]/sum_events
        
        if occurrence:
            break

        ################ Uninfected predator #####################
        for i in range(0,predator.shape[0]):
            for j in range(0,prey.shape[0]):
                if URN > P and URN <= P + predator_growth[i,j]/sum_events:
                    bx = 0 # uninfected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 1 # uninfected predator increases by one
                    mby = i # row number in uninfected predator array
                    mbx = j # row number in uninfected prey array
                    occurrence = True
                    break
                P += predator_growth[i,j]/sum_events

                if URN > P and URN <= P + predator_non_growth[i,j]/sum_events:
                    bx = 0 # uninfected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 2 # nothing happens to predator
                    mby = i # row number in uninfected predator array
                    mbx = j # row number in uninfected prey array
                    occurrence = True
                    break
                P += predator_non_growth[i,j]/sum_events

        if occurrence:
            break
        
        for i in range(0,predator.shape[0]):          
            if URN > P and URN <= P + predator_death[i]/sum_events:
                bx = 2 # nothing happens to prey
                bz = 2 # nothing happens to free-living parasite
                by = 0 # uninfected predator decreases by one  
                mby = i # row number in uninfected predator array
                occurrence = True
                break
            P += predator_death[i]/sum_events
        
        if occurrence:
            break

        for i in range(0,predator.shape[0]):
            for j in range(0,infected_prey.shape[0]):
                if URN > P and URN <= P + predator_exposure_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 1 # uninfected predator increases by one
                    mby = i # row number in uninfected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += predator_exposure_growth[i,j]/sum_events

                if URN > P and URN <= P + predator_exposure_non_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 2 # nothing happens to predator
                    mby = i # row number in uninfected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += predator_exposure_non_growth[i,j]/sum_events

                if URN > P and URN <= P + predator_infection_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 1 # free-living parasite increases by nz (parasite inside the prey reproduces in predator)
                    by = 3 # predator is infected and increases by one
                    mby = i # row number in uninfected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += predator_infection_growth[i,j]/sum_events
                            
                if URN > P and URN <= P + predator_infection_non_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 1 # free-living parasite increases by nz (parasite inside the prey reproduces in predator)
                    by = 4 # predator is infected and does not reproduce
                    mby = i # row number in uninfected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += predator_infection_non_growth[i,j]/sum_events

        if occurrence:
            break
                            
        ################ Infected predator #####################
        for i in range(0,infected_predator.shape[0]):
            for j in range(0,prey.shape[0]):
                if URN > P and URN <= P + infected_predator_growth[i,j]/sum_events:
                    bx = 0 # ancestral prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 6 # predator increases by one
                    mby = i # row number in infected predator array
                    mbx = j # row number in uninfected prey array
                    occurrence = True
                    break
                P += infected_predator_growth[i,j]/sum_events
                                
                if URN > P and URN <= P + infected_predator_non_growth[i,j]/sum_events:
                    bx = 0 # ancestral prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 2 # nothing happens to predator
                    mby = i # row number in infected predator array
                    mbx = j # row number in uninfected prey array
                    occurrence = True
                    break
                P += infected_predator_non_growth[i,j]/sum_events
        
        if occurrence:
            break

        for i in range(0,infected_predator.shape[0]):
            for j in range(0,infected_prey.shape[0]):
                if URN > P and URN <= P + infected_predator_exposure_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 6 # uninfected predator increases by one
                    mby = i # row number in infected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += infected_predator_exposure_growth[i,j]/sum_events
                                    
                if URN > P and URN <= P + infected_predator_exposure_non_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 2 # nothing happens to free-living parasite
                    by = 2 # nothing happens to predator
                    mby = i # row number in infected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += infected_predator_exposure_non_growth[i,j]/sum_events
                                        
                if URN > P and URN <= P + infected_predator_infection_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 1 # free-living parasite increases by nz (parasite inside the prey reproduces in predator)
                    by = 6 # uninfected predator increases by one
                    mby = i # row number in infected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += infected_predator_infection_growth[i,j]/sum_events
                                        
                if URN > P and URN <= P + infected_predator_infection_non_growth[i,j]/sum_events:
                    bx = 5 # infected prey decreases by one
                    bz = 1 # free-living parasite increases by nz (parasite inside the prey reproduces in predator)
                    by = 2 # nothing happens to predator
                    mby = i # row number in infected predator array
                    mbi = j # row number in infected prey array
                    occurrence = True
                    break
                P += infected_predator_infection_non_growth[i,j]/sum_events

        if occurrence:
            break

        for i in range(0,infected_predator.shape[0]):
            if URN > P and URN <= P + infected_predator_death[i]/sum_events:
                bx = 2 # nothing happens to prey
                bz = 2 # nothing happens to free-living parasite
                by = 5 # infected predator decreases by one  
                mby = i # row number in infected predator array
                occurrence = True
                break
            P += infected_predator_death[i]/sum_events
    
        if occurrence:
            break

##################### PREY EVENTS #########################
    if(bx == 1): # prey reproduces
        new_prey = np.array(my_mutation_loci(n_loci, mx, prey[mbx,2])) # new genotype that results after reproduction (may have a mutation or not)
        prey = loop_to_compare_array(prey, new_prey, Time) # either append new genotype or add individual if already exists
        store_prey = loop_to_store_array(store_prey, new_prey, Time) # update the store arrays of all the prey genotypes
        temporal_prey = str(new_prey) # check if offspring was super-resistant
        if(temporal_prey == str(np.ones(n_loci)) and not emergence_st_prey):
            emergence_st_prey = True

    elif(bx == 0): # prey dies
        prey[mbx,1] -= 1 # decrease one prey individual
        store_prey = loop_to_remove_array(store_prey, prey[mbx,2])
    
    elif(bx == 3):# parasite infects prey
        prey[mbx,1] -= 1  # prey gets infected
        infected_prey = loop_infection(infected_prey, prey[mbx,2], parasite[mbz,2]) # either append new genotype infected or add individual if already exists

    elif(bx == 5): # infected prey dies
        infected_prey[mbi,1] -= 1 # decrease one infected prey from that particular row (where event happened)
        store_prey = loop_to_remove_array(store_prey, infected_prey[mbi,2])

    elif(bx == 6): # infected prey reproduces
        new_infected_prey = np.array(my_mutation_loci(n_loci, mx, infected_prey[mbi,2])) # new genotype that results after reproduction (may have a mutation or not)
        prey = loop_to_compare_array(prey, new_infected_prey, Time) # either append new genotype or add individual if already exists (it is added to uninfected prey array beacuse parasite is not transmitted from parent to offspring
        store_prey = loop_to_store_array(store_prey, new_infected_prey, Time) # update the store arrays of all the prey genotypes
        temporal_prey = str(new_infected_prey) # check if offspring was super-resistant
        if(temporal_prey == str(np.ones(n_loci)) and not emergence_st_prey):
            emergence_st_prey = True

################### PARASITE EVENTS #######################
    if(bz == 0): # free-living parasite dies
        parasite[mbz,1] -= 1
        store_parasite = loop_to_remove_array(store_parasite, parasite[mbz,2])
    
    elif(bz == 1): # parasite infects a predator
        new_parasite = np.zeros(n_z)
        for i in range(0,n_z): # repeat this for all parasite offspring
            new_parasite = np.array(my_mutation_loci(n_loci, mz, infected_prey[mbi,3])) # check for mutations
            parasite = loop_to_compare_array(parasite, new_parasite, Time) # add new genotype (may contain mutation or may not)
            store_parasite = loop_to_store_array(store_parasite, new_parasite, Time) # update the store arrays of all the parasite genotypes
        # If all-resistant or all-infective genotypes emerged, record the id and emergence time
        temporal_parasite = str(new_parasite) # check if offspring was super-infective
        if(temporal_parasite == str(np.ones(n_loci)) and not emergence_st_parasite):
            emergence_st_parasite = True
    
################### PREDATOR EVENTS ########################
    if(by == 1): # predator reproduces
        new_predator = np.array(my_mutation_loci(n_loci, my, predator[mby,2])) # check for mutation in predator offspring
        predator = loop_to_compare_array(predator, new_predator, Time) # add new genotype (may contain mutation or may not)
        store_predator = loop_to_store_array(store_predator, new_predator, Time) # update the store arrays of all the predator genotypes
        temporal_predator = str(new_predator) # check if offspring was super-resistant
        if(temporal_predator == str(np.ones(n_loci)) and not emergence_st_predator):
            emergence_st_predator = True

    elif(by == 0): # predator dies
        predator[mby,1] -= 1 # uninfected predator of that genotype (row in uninfected predator array) decreases by one
        store_predator = loop_to_remove_array(store_predator, predator[mby,2])

    elif(by == 3): # predator gets infected and reproduces
        # Infection part
        predator[mby,1] -= 1
        infected_predator = loop_infection(infected_predator, predator[mby,2], infected_prey[mbi,3]) 
        # Reproduction part
        new_predator = np.array(my_mutation_loci(n_loci, my, predator[mby,2])) # the predator is also reproducing, create new genotype (may or may not have mutation)
        predator = loop_to_compare_array(predator, new_predator, Time) # add new genotype to the uninfected predator array (parasite is not transmitted from parent to offspring)
        store_predator = loop_to_store_array(store_predator, new_predator, Time) # update the store arrays of all the predator genotypes
        temporal_predator = str(new_predator) # check if offspring was super-resistant
        if(temporal_predator == str(np.ones(n_loci)) and not emergence_st_predator):
            emergence_st_predator = True

    elif(by == 4): # predator gets infected and does not reproduce
        # Infection part
        predator[mby,1] -= 1
        infected_predator = loop_infection(infected_predator, predator[mby,2], infected_prey[mbi,3])

    elif(by == 5): # infected predator dies
        infected_predator[mby,1] -= 1
        store_predator = loop_to_remove_array(store_predator, infected_predator[mby,2])
        
    elif(by == 6): # infected predator reproduces
        # Reproduction part
        new_infected_predator = np.array(my_mutation_loci(n_loci, my, infected_predator[mby,2])) # check for mutations in genotype of offspring
        predator = loop_to_compare_array(predator, new_infected_predator, Time) # add genotype to uninfected predator array
        store_predator = loop_to_store_array(store_predator, new_infected_predator, Time) # update the store arrays of all the predator genotypes
        temporal_predator = str(new_infected_predator) # check if offspring was super-resistant
        if(temporal_predator == str(np.ones(n_loci)) and not emergence_st_predator):
            emergence_st_predator = True

    # Record effective genotypes (genotypes that reached at least 2% of the total population size)
    # If they have a "1" in column 4 of array, they are effective genotypes
    for i in range(0, store_prey.shape[0]):
        if sum(store_prey[:,1]) >= 1:
            if(store_prey[i,1] >= (sum(store_prey[:,1]) * 0.02) and store_prey[i,3] == 0): # if the genotype reached at least 2% of the total population size and it has not been recorded
                store_prey[i,3] = 1 # becomes effective genotype

    for i in range(0, store_parasite.shape[0]):
        if sum(store_parasite[:,1]) >= 1:
            if(store_parasite[i,1] >= (sum(store_parasite[:,1]) * 0.02) and store_parasite[i,3] == 0): # if the genotype reached at least 2% of the total population size and it has not been recorded
                store_parasite[i,3] = 1 # becomes effective genotype

    for i in range(0, store_predator.shape[0]):
        if sum(store_predator[:,1]) >= 1:
            if(store_predator[i,1] >= (sum(store_predator[:,1]) * 0.02) and store_predator[i,3] == 0): # if the genotype reached at least 2% of the total population size and it has not been recorded
                store_predator[i,3] = 1 # becomes effective genotype

    # Record lifespan of the extinct effective genotypes
    for i in range(0, store_prey.shape[0]):
        if store_prey[i,1] == 0 and store_prey[i,3] == 1: # if genotype went extinct and was an effective genotype (more than 2% of total population size)
            store_prey_genotypic_lifespan.append(Time - store_prey[i,4]) # Time when genotype went extinct minus time when the genotype emerged (store_prey[i,4] is the time it first emerged)

    for i in range(0, store_parasite.shape[0]):
        if store_parasite[i,1] == 0 and store_parasite[i,3] == 1: # if genotype went extinct and was an effective genotype (more than 2% of total population size)
           store_parasite_genotypic_lifespan.append(Time - store_parasite[i,4]) # Time when genotype went extinct minus time when the genotype emerged (store_parasite[i,4] is the time it first emerged)

    for i in range(0, store_predator.shape[0]):
        if store_predator[i,1] == 0 and store_predator[i,3] == 1: # if genotype went extinct and was an effective genotype (more than 2% of total population size)
            store_predator_genotypic_lifespan.append(Time - store_predator[i,4]) # Time when genotype went extinct minus time when the genotype emerged (store_predator[i,4] is the time it first emerged)

    # Update population sizes
    # free-living individuals
    sum_uninfected_prey = np.sum(prey[:,1]) # uninfected prey
    sum_parasite = np.sum(store_parasite[:,1]) # free-living parasites
    sum_uninfected_predator = np.sum(predator[:,1]) # uninfected predator
    # infected hosts
    sum_infected_prey = np.sum(infected_prey[:,1]) # infected prey
    sum_infected_predator = np.sum(infected_predator[:,1]) # infected predator
    # total
    sum_prey = sum_uninfected_prey + sum_infected_prey # sum prey (uninfected and infected)
    sum_all_parasite = sum_parasite + sum_infected_prey # sum parasites for extinction purposes (free-living and those inside prey)
    sum_predator = sum_uninfected_predator + sum_infected_predator # sum predator (uninfected and infected)
    
    # Update number of genotypes
    prey_genotype = 0
    for i in range(0, store_prey.shape[0]):
        if store_prey[i,1] != 0:
            prey_genotype += store_prey[i,3]
    parasite_genotype = 0
    for i in range(0, store_parasite.shape[0]):
        if store_parasite[i,1] != 0:
            parasite_genotype += store_parasite[i,3]
    predator_genotype = 0
    for i in range(0, store_predator.shape[0]):
        if store_predator[i,1] != 0:
            predator_genotype += store_predator[i,3]

# if prey, parasite, and predator go extinct record which one went extinct
    if(sum_prey <= 0 and not extinction_prey):
        extinction_prey = True

    if(sum_all_parasite <= 0 and not extinction_parasite):
        extinction_parasite = True

    if(sum_predator <= 0 and not extinction_predator):
        extinction_predator = True
    
    ## RECORD ABUNDANCE OF ALL ENTITIES IN EACH TIME STEP ##
    if Time > n:

        # STORE NUMBER OF INDIVIDUALS
        store_sum_prey.append(sum_prey)
        store_sum_parasite.append(sum_parasite)
        store_sum_predator.append(sum_predator)

        # STORE NUMBER OF GENOTYPES
        store_prey_genotype.append(prey_genotype)
        store_parasite_genotype.append(parasite_genotype)
        store_predator_genotype.append(predator_genotype)

        n += 0.1

    if(sum_parasite <= 0 and sum_predator <= 0): # break while loop if only the prey remains
        break

    # Advance in Gillespie time step
    Time += dt_next_event # continuous time simulation
        
### SIMULATION FINISHES WHEN REACHING MAX TIME ###
# Record coexistence/extinctions
if not extinction_prey and not extinction_parasite and not extinction_predator:
    coexistence = "1"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "0"
    extinction = "0"
elif not extinction_prey and extinction_parasite and not extinction_predator:
    coexistence = "0"
    coexistence_predator_and_prey = "1"
    coexistence_prey = "0"
    extinction = "0"
elif not extinction_prey and extinction_parasite and extinction_predator:
    coexistence = "0"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "1"
    extinction = "0"     
else:
    coexistence = "0"
    coexistence_predator_and_prey = "0"
    coexistence_prey = "0"
    extinction = "1"

# Record emergence super-types
if emergence_st_prey:
    emergence_prey = "1"
if emergence_st_parasite:
    emergence_parasite = "1"
if emergence_st_predator:
    emergence_predator = "1"

# Record genotypic lifespan of effective genotypes (those still alive after the simulation finished and have more than 2% of total population size)
for i in range(0, store_prey.shape[0]):
    if store_prey[i,1] != 0 and store_prey[i,3] == 1: # if genotype is still alive and it is an effective genotype (more than 2% of total population size)
        store_prey_genotypic_lifespan.append(max_time - store_prey[i,4]) # Time when genotype went extinct minus time when the genotype emerged (store_prey[i,4] is the time it first emerged)

for i in range(0, store_parasite.shape[0]):
    if store_parasite[i,1] != 0 and store_parasite[i,3] == 1: # if genotype is still alive and it is an effective genotype (more than 2% of total population size)
        store_parasite_genotypic_lifespan.append(max_time - store_parasite[i,4]) # Time when genotype emerged minus time when genotype went extinct (store_parasite[i,4] is the time it first emerged)

for i in range(0, store_predator.shape[0]):
    if store_predator[i,1] != 0 and store_predator[i,3] == 1: # if genotype is still alive and it is an effective genotype (more than 2% of total population size)
        store_predator_genotypic_lifespan.append(max_time - store_predator[i,4]) # Time when genotype emerged minus time when genotype went extinct (store_predator[i,4] is the time it first emerged)

# Record average abundances
if store_sum_prey != []: #avoid division by zero
    av_prey_individuals = sum(store_sum_prey[:]) / len(store_sum_prey)
if store_sum_parasite != []: #avoid division by zero
    av_parasite_individuals = sum(store_sum_parasite[:]) / len(store_sum_parasite)
if store_sum_predator != []: #avoid division by zero
    av_predator_individuals = sum(store_sum_predator[:]) / len(store_sum_predator)
# Record average genotypes
if store_prey_genotype != []: #avoid division by zero
    av_prey_genotypes = sum(store_prey_genotype[:]) / len(store_prey_genotype)
if store_parasite_genotype != []: #avoid division by zero
    av_parasite_genotypes = sum(store_parasite_genotype[:]) / len(store_parasite_genotype)
if store_predator_genotype != []: #avoid division by zero
    av_predator_genotypes = sum(store_predator_genotype[:]) / len(store_predator_genotype)
# Record average genotypic lifespans
if store_prey_genotypic_lifespan != []: #avoid division by zero
    av_prey_genotypic_lifespan = sum(store_prey_genotypic_lifespan[:]) / len(store_prey_genotypic_lifespan)
if store_parasite_genotypic_lifespan != []: #avoid division by zero
    av_parasite_genotypic_lifespan = sum(store_parasite_genotypic_lifespan[:]) / len(store_parasite_genotypic_lifespan)
if store_predator_genotypic_lifespan != []: #avoid division by zero
    av_predator_genotypic_lifespan = sum(store_predator_genotypic_lifespan[:]) / len(store_predator_genotypic_lifespan)

# Save output extinctions and emergence of super-types (all-resistant host and all-infective parasite genotypes)
output_file.write(str(rp) + "," + str(rx) + "," + str(av_prey_individuals) + "," + str(av_parasite_individuals) + "," + str(av_predator_individuals) + "," + str(av_prey_genotypes) + "," + str(av_parasite_genotypes) + "," + str(av_predator_genotypes) + "," + str(av_prey_genotypic_lifespan) + "," + str(av_parasite_genotypic_lifespan) + "," + str(av_predator_genotypic_lifespan) + "," + str(emergence_prey) + "," + str(emergence_parasite) + "," + str(emergence_predator) + "," + str(coexistence) + "," + str(coexistence_predator_and_prey) + "," + str(coexistence_prey) + "," + str(extinction) + "\n")

# Close files
output_file.close() # This file is for average abundance, average genotypes, emergence of super-types and coexistence from a single realisation