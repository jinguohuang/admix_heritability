# simulate a rather simple admixed population with equal parental pop size

import msprime 
import tskit
import numpy as np
import argparse
# variables
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")
# required input: segment number(seed number) 
req_grp.add_argument("--seg",dest="seg",help="segment number",type=str,required=True)
args=parser.parse_args()
print(args)

seed=args.seg #random seed
length=1_000_000  # segment length, set as 1M
ancestry_seed=seed # random seed for ancestry simulation
mutation_seed=seed # random seed for mutation simulation
admsamplesize=10_000 #sample size for admix population, set as 10k
Asamplesize=1_000 # sample size for parental pop A, set as 1k
Bsamplesize=1_000 # sample size for parental pop B, set as 1k

# the demography
demography = msprime.Demography()
demography.add_population(name="A", initial_size=10_000)
demography.add_population(name="B", initial_size=10_000)
demography.add_population(name="ADMIX", initial_size=10_000)
demography.add_population(
    name="ANC",
    description="Ancestral equilibrium",
    initial_size=10_000,
)
demography.set_symmetric_migration_rate(["A", "B"], 2.5e-5)
demography.add_admixture(
    time=16, 
    derived="ADMIX", 
    ancestral=["A", "B"], 
    proportions=[0.5, 0.5],
)

demography.add_census(time=17)
demography.add_population_split(
    time=2040, 
    derived=["A", "B"], 
    ancestral="ANC")
#demography.debug()
demography.sort_events() # sort events

# simulate
ts = msprime.sim_ancestry(
    {"ADMIX": admsamplesize, "A": Asamplesize, "B": Bsamplesize}, 
    demography=demography, 
    recombination_rate=1e-08,
    sequence_length=length,
    random_seed=ancestry_seed
)
# add mutations
ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=mutation_seed)
#SVG(ts.draw_svg(y_axis=True, time_scale="log_time"))

print("writing genotype to vcf file")
#export to diploid to validate simulated global anc
filename=f'sim_seg{seed}.vcf'
n_dip_indv = int(ts.num_samples/2)
#indv_names = [f"sim_{str(i)}indv" for i in range(n_dip_indv)]
indv_names = [f"adm_{str(i)}indv" for i in range(admsamplesize)] +\
    [f"A_{str(i)}indv" for i in range(Asamplesize)] +\
    [f"B_{str(i)}indv" for i in range(Bsamplesize)]
with open(filename, "w") as vcf_file:
    ts.write_vcf(vcf_file, 
                 contig_id=seed, 
                 individual_names=indv_names)

##### output the local ancestry chunk of each person
import tskit
import numpy as np
import msprime 

# replace the parent to population label
def get_population_id(node, ts):
    return ts.tables.nodes.population[node]

def local_ancestry(Tcensus, ts, filename):
    # select nodes that are Tcensus generations ago
    ancestors_Tcensus_gens = np.where(ts.tables.nodes.time == Tcensus)[0]
    # which of your samples descend from which ancestors
    #ancestrytable_Tcensus_gens = ts.tables.link_ancestors(
    #samples=ts.samples(), ancestors=ancestors_Tcensus_gens)
    ### loop on a small set 
    a_list = ts.samples() # samples list
    length = len(a_list) #total number of sample size
    lu=100 #loop unit:loop 100 at a time
    ln=length//lu +1 #loop number
    x=0  #initial number for the start index
    for i in range(1,ln):
        mysample = a_list[x:(x+lu)]
        x=lu*i #new start index
        ancestrytable_Tcensus_gens = ts.tables.link_ancestors(
        samples=mysample, ancestors=ancestors_Tcensus_gens)
        import pandas as pd
        ancestry_table = pd.DataFrame(
            data = {
                'left': ancestrytable_Tcensus_gens.left,
                'right': ancestrytable_Tcensus_gens.right,
                'populations': [get_population_id(u, ts) for u in ancestrytable_Tcensus_gens.parent],
                'child': ancestrytable_Tcensus_gens.child
            }
        )
        #return ancestry_table
        # output the ancestry table to a file, append mode so it will accumulate for each loop
        ancestry_table.to_csv(f'lanc_{filename}.txt', sep='\t', encoding='utf-8', index=False, header=False, mode='a')

# run the function
print("writing local ancestry file")
local_ancestry(17, ts, filename)

