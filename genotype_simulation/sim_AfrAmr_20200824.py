# Simulate genotype data of African American
# NOTE: this model will cause admixture from Asian
# African American: European with African admixture 16 generations ago
import msprime, sys
from math import log
from math import exp
#seed = int(sys.argv[1])
seed = 20

mu=1.25e-8 # mutation rate per bp
rho=1e-8 # recombination rate per bp
nbp=1e8 # generate 100Mb
N0=7310 # initial population size
Thum=5920 # time (gens) of advent of modern humans
Naf=14474 # size of african population
Tooa=2040 # number of generations back to Out of Africa
Nb=1861 # size of out of africa population
mafb=1.5e-4 # migration rate Africa and Out-of-Afica
Teu=920 # number generations back to Asia-Euroupe split
Neu=1032; Nas=554 # bottleneck population sizes
mafeu=2.5e-5;mafas=7.8e-6;meuas=3.11e-5 #mig.rates
reu=0.0038 #growth rate per generation in Europe
ras=0.0048 #growth rate per generation in Asia
Tadmix=16 # time of admixture
Nadmix=30000 # initial size of admixed population
radmix=.05 # growth rate of admixed population

# pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is admixed
refsamplesize=200
admsamplesize=1000 #sample size of admixted pop
pop_config = [
    msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Naf, growth_rate=0.0), #100
    msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Neu*exp(reu*Teu), growth_rate=reu), #100
    msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Nas*exp(ras*Teu), growth_rate=ras), #100
    msprime.PopulationConfiguration(sample_size=admsamplesize, initial_size=Nadmix*exp(radmix*Tadmix), growth_rate=radmix)] #500
mig_mat=[
    [0,mafeu,mafas,0],
    [mafeu,0,meuas,0],
    [mafas,meuas,0,0],
    [0,0,0,0]]

# Admixture event, 0.8 Africa, 0.2 Europe
admixture_event = [msprime.MassMigration(time=Tadmix, source=3, destination=0, proportion=0.8), # African came to America
                  msprime.MassMigration(time=Tadmix+0.0001, source=3, destination=1, proportion=1.0)] # European came to America

# Asia and Europe split
eu_event = [
    msprime.MigrationRateChange(time=Teu, rate=0.0),
    msprime.MassMigration(time=Teu+0.0001, source=2, destination=1, proportion=1.0),
    msprime.PopulationParametersChange(time=Teu+0.0002, initial_size=Nb, growth_rate=0.0, population_id=1),
    msprime.MigrationRateChange(time=Teu+0.0003, rate=mafb, matrix_index=(0,1)),
    msprime.MigrationRateChange(time=Teu+0.0003, rate=mafb, matrix_index=(1,0))
]

# Out of Africa event
ooa_event = [
    msprime.MigrationRateChange(time=Tooa, rate=0.0),
    msprime.MassMigration(time=Tooa+0.0001, source=1, destination=0, proportion=1.0)
]

# initial population size
init_event=[msprime.PopulationParametersChange(time=Thum, initial_size=N0, population_id=0)]

events = admixture_event + eu_event + ooa_event + init_event

treeseq = msprime.simulate(
    population_configurations=pop_config,
    migration_matrix=mig_mat,
    demographic_events=events,
    length=nbp,
    recombination_rate=rho, 
    mutation_rate=mu,
    random_seed=seed)

#print((treeseq.num_trees, treeseq.num_sites))

# output genotype data (diploid)
# deal with Plink error "Sample ID ends with "_0", which induces an invalid IID of '0'." 
n_dip_indv = int(treeseq.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
with open("AfricanAmerican.vcf", "w") as vcf_file:
    treeseq.write_vcf(vcf_file, ploidy=2, individual_names=indv_names)
