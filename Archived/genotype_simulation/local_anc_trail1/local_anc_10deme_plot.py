# functions to get local ancestry with msprime
#%matplotlib inline
#%config InlineBackend.figure_format = 'svg'
import random
import collections
import msprime
import numpy as np
import seaborn as sns
import multiprocessing
import matplotlib.pyplot as plt

#from IPython.display import SVG
#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:90% !important; }</style>"))
# new imports
import sys
import itertools
from matplotlib import collections  as mc
from intervaltree import Interval, IntervalTree # thanks for the suggestion

# merge interval
def merge_intervals(intervals):
    """takes a list of (population-specific) ancestry intervals and combines them into contiguous intervals"""
    merged_intervals = []
    start = None
    stop = None
    for newstart, newstop in intervals:
        if start is None: # start first
            start = newstart
            stop = newstop
        elif np.allclose(newstart, stop): # extend
            stop = newstop
        else:  # end previous and start new 
            merged_intervals.append((start, stop))
            start = newstart
            stop = newstop
    # get terminal interval
    merged_intervals.append((start, newstop))
    return(np.array(merged_intervals))

# plot the merged segments
def plot_merged_segments(merged_segments_from_pop):
    """rough plot of local ancestry, plots the output of find_local_ancestry()"""
    fig, ax = plt.subplots(figsize=(10,2))
    for anc_pop, ms in merged_segments_from_pop.items():
        lines = zip(zip(ms[:,0], itertools.repeat(anc_pop)), zip(ms[:,1], itertools.repeat(anc_pop)))
        lc = mc.LineCollection(lines, linewidths=8)
        ax.add_collection(lc)
    ax.margins(0.1)
    maxpop = max(merged_segments_from_pop.keys())
    plt.ylim(-.2, maxpop+.2)
    plt.yticks(range(maxpop+1), range(maxpop+1))
    plt.ylabel('source population')
    plt.xlabel('bp position')

# find the local ancestry
def find_local_ancestry(sample, time, ts):
    """returns a dict.
    keys are the ancestral populations.
    values are the contiguous tracks of ancestry for the sample inherited from that pop.  
    Pops are defined at time """
    # put relevant migrations on an interval tree (otherwise too slow)
    mig_int_tree = IntervalTree()
    for migration in ts.migrations():
        if migration.time < time:
            mig_int_tree[migration.left:migration.right] = migration
    
    # for each tree, find the oldest node prior to [time]
    ancestor_before_timex_of_tree = dict()
    for tree in ts.trees():
        target = sample
        node_time = tree.time(target)
        parent_node = tree.parent(target)
        if parent_node != msprime.NULL_NODE:
            parent_time = tree.time(tree.parent(target))
        else:
            parent_time = time+1 
        while parent_time < time:
            node_time = parent_time
            target = tree.parent(target)
            parent_node = tree.parent(target)
            if parent_node != msprime.NULL_NODE:
                parent_time = tree.time(tree.parent(target))
            else:
                parent_time = time+1 
        ancestor_before_timex_of_tree[tree.index] = target
    
    # loop over trees and their relevant ancestors, find the pop that contributes to this sample 
    pop_at_time_of_parent = dict()
    pop_at_time_of_tree = dict()

    for tree in ts.trees():
        parent_node = ancestor_before_timex_of_tree[tree.index]
        # the original population of the node
        pop_at_time_of_parent[parent_node] = tree.population(parent_node)
        # trace the node back until [time], accounting for relevant migrations
        overlapping_intervals = mig_int_tree[tree.interval[0]]
        if len(overlapping_intervals) > 0:
            overlapping_migrations = [interval.data for interval in overlapping_intervals if interval.data.node == parent_node]
            if len(overlapping_migrations) > 0:
                #print len(overlapping_migrations)
                overlapping_migrations = sorted(overlapping_migrations, key = lambda x : x.time)
                #last_mig = overlapping_migrations[-1]#.pop()
                last_mig = overlapping_migrations.pop()
                #assert (pop_at_time_of_parent[parent_node] == mig.source), (parent_node, mig.time, pop_at_time_of_parent[parent_node], mig.source, mig.dest)
                pop_at_time_of_parent[parent_node] = last_mig.dest
        pop_at_time_of_tree[tree.index] = pop_at_time_of_parent[parent_node]
    
    # for each tree, record the interval and pop
    intervals_of_tree = dict()
    for tree in ts.trees():
        intervals_of_tree[tree.index] = tree.interval
    segments_from_pop = collections.defaultdict(list)
    for ti, anc_pop in pop_at_time_of_tree.items():
        segments_from_pop[anc_pop].append(intervals_of_tree[ti])

    # merge adjacent intervals from the same population
    merged_segments_from_pop = dict()
    for anc_pop in segments_from_pop:
        merged_segments_from_pop[anc_pop] = merge_intervals(segments_from_pop[anc_pop])  
    return(merged_segments_from_pop)


# my function to simulate 
# African American: European with African admixture 16 generations ago
# Note, this version deleted EUR-ASN split and all ASN-related migration event
import msprime, sys
from math import log
from math import exp
#seed = int(sys.argv[1])
#seed = 20
def run_simulation(random_seed=None):
    # 10 demes model
    mu=1.25e-8 # mutation rate per bp
    rho=1e-8 # recombination rate per bp
    nbp=1e8 # generate 100Mb
    #nbp=1e7 # generate 10Mb
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
    mdeme=0.05/2 # migration rate between neighboring demes

    # pop0 is Africa, pop1 is Europe,  pop2-10 is admixed
    refsamplesize=200
    #admsamplesize=1000 #sample size of admixted pop
    # sample size of each admixed pop
    adm1=30
    adm2=0
    adm3=15
    adm4=30
    adm5=15
    adm6=91
    adm7=258
    adm8=592
    adm9=742
    adm10=227
    pop_config = [
        msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Naf), #AFR
        msprime.PopulationConfiguration(sample_size=refsamplesize, initial_size=Neu*exp(reu*Teu)), #EUR
        msprime.PopulationConfiguration(sample_size=adm1, initial_size=Nadmix*exp(radmix*Tadmix)), #adm1
        msprime.PopulationConfiguration(sample_size=adm2, initial_size=Nadmix*exp(radmix*Tadmix)), #adm2
        msprime.PopulationConfiguration(sample_size=adm3, initial_size=Nadmix*exp(radmix*Tadmix)), #adm3
        msprime.PopulationConfiguration(sample_size=adm4, initial_size=Nadmix*exp(radmix*Tadmix)), #adm4
        msprime.PopulationConfiguration(sample_size=adm5, initial_size=Nadmix*exp(radmix*Tadmix)), #adm5
        msprime.PopulationConfiguration(sample_size=adm6, initial_size=Nadmix*exp(radmix*Tadmix)), #adm6
        msprime.PopulationConfiguration(sample_size=adm7, initial_size=Nadmix*exp(radmix*Tadmix)), #adm7
        msprime.PopulationConfiguration(sample_size=adm8, initial_size=Nadmix*exp(radmix*Tadmix)), #adm8
        msprime.PopulationConfiguration(sample_size=adm9, initial_size=Nadmix*exp(radmix*Tadmix)), #adm9
        msprime.PopulationConfiguration(sample_size=adm10, initial_size=Nadmix*exp(radmix*Tadmix)), #adm10 
    ]
    mig_mat=[
        [    0,mafeu,0,0,0,0,0,0,0,0,0,0], #AFR
        [mafeu,0,    0,0,0,0,0,0,0,0,0,0], #EUR
        [    0,0,0,mdeme,0,0,0,0,0,0,0,0], #ADM1
        [0,0,mdeme,0,mdeme,0,0,0,0,0,0,0], #ADM2
        [0,0,0,mdeme,0,mdeme,0,0,0,0,0,0], #ADM3
        [0,0,0,0,mdeme,0,mdeme,0,0,0,0,0], #ADM4
        [0,0,0,0,0,mdeme,0,mdeme,0,0,0,0], #ADM5
        [0,0,0,0,0,0,mdeme,0,mdeme,0,0,0], #ADM6
        [0,0,0,0,0,0,0,mdeme,0,mdeme,0,0], #ADM7
        [0,0,0,0,0,0,0,0,mdeme,0,mdeme,0], #ADM8
        [0,0,0,0,0,0,0,0,0,mdeme,0,mdeme], #ADM9
        [0,0,0,0,0,0,0,0,0,0,mdeme,0], #ADM10
    ]

    # Admixture event, 0.8 Africa, 0.2 Europe in total
    # 0 Afr, 1 Eur, source 2-12 admixed 
    admixture_event = [
        msprime.MigrationRateChange(time=Tadmix, rate=0.0), #newly added, change the migration rate to 0
        msprime.MassMigration(time=Tadmix, source=2, destination=0, proportion=0.05), # Afr 0.05
        msprime.MassMigration(time=Tadmix+0.0001, source=2, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=3, destination=0, proportion=0.15), # Afr 0.15
        msprime.MassMigration(time=Tadmix+0.0001, source=3, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=4, destination=0, proportion=0.25), # Afr 0.25
        msprime.MassMigration(time=Tadmix+0.0001, source=4, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=5, destination=0, proportion=0.35), # Afr 0.35
        msprime.MassMigration(time=Tadmix+0.0001, source=5, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=6, destination=0, proportion=0.45), # Afr 0.45
        msprime.MassMigration(time=Tadmix+0.0001, source=6, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=7, destination=0, proportion=0.55), # Afr 0.55
        msprime.MassMigration(time=Tadmix+0.0001, source=7, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=8, destination=0, proportion=0.65), # Afr 0.65
        msprime.MassMigration(time=Tadmix+0.0001, source=8, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=9, destination=0, proportion=0.75), # Afr 0.75
        msprime.MassMigration(time=Tadmix+0.0001, source=9, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=10, destination=0, proportion=0.85), # Afr 0.85
        msprime.MassMigration(time=Tadmix+0.0001, source=10, destination=1, proportion=1.0), # Eur
        msprime.MassMigration(time=Tadmix, source=11, destination=0, proportion=0.95), # Afr 0.95
        msprime.MassMigration(time=Tadmix+0.0001, source=11, destination=1, proportion=1.0), # Eur                   
    ] 

    # Out of Africa event
    ooa_event = [
        msprime.MigrationRateChange(time=Tooa, rate=0.0),
        msprime.MassMigration(time=Tooa+0.0001, source=1, destination=0, proportion=1.0)
    ]

    # initial population size
    init_event=[msprime.PopulationParametersChange(time=Thum, initial_size=N0, population_id=0)]

    #events = admixture_event + eu_event + ooa_event + init_event
    events = admixture_event + ooa_event + init_event 

    ts = msprime.simulate(
        population_configurations=pop_config,
        migration_matrix=mig_mat,
        demographic_events=events,
        length=nbp,
        recombination_rate=rho, 
        mutation_rate=mu,
        record_migrations=True, # newly added, Needed for tracking segments.
        random_seed=random_seed)
    return ts


#my function to calculate global anc from local anc results
def get_global_anc(x):
    j=0
    if (x==None) is True:
        j=j
    else:
        for i in range(len(x)):
            j=j+x[i,1]-x[i,0] 
    return j/1e8

# simulate!
ts = run_simulation(1)
# plot
for sample in [100, 300, 1000, 2000]:
    for time in [20, 2030]:
        la = find_local_ancestry(sample = sample, time = time, ts=ts)
        plot_merged_segments(la)
        plt.title('Sample: {}, Time: {}'.format(sample, time))
        plt.savefig('Sample{}Time{}'.format(sample, time)+ '.jpg', bbox_inches='tight')


