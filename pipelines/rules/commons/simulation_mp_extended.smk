import numpy as np
import pandas as pd
import pybedtools
import tskit
import math
import pyranges as pr
from multiprocessing import Process, Manager



### CONFIG ###
#configfile: "config/test_gorilla.config.yaml"

params_set = "test"

feature_id = config["feature_id"]
feature_config = config["feature_config"]
nfeature = config["nfeature"]
nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
geno_state_list = config["geno_states"]
win_step = config["win_step"]
cutoff_num = config["cutoff_num"]
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])

output_prefix = config["output_prefix"]
nrep = config["nrep"][params_set]
seq_len = config["seq_len"]
demog_id = config["demog_id"][params_set]
demes_file = config["demes"][params_set]
mut_rate = config["mut_rate"][params_set]
rec_rate = config["rec_rate"][params_set]
ploidy = config["ploidy"]
ref_id = config["ref_id"][params_set]
tgt_id = config["tgt_id"][params_set]
src_id = config["src_id"][params_set]


np.random.seed(config["seed"])
seed_list = np.random.random_integers(1, 2**31, nrep)

#output_dir = output_dir + f'/{demog_id}/nref_{nref}/ntgt_{ntgt}'
output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'



def create_final_tracts(tract_file):
    import pyranges as pr
    tract_file.columns=["Chromosome", "Start", "End", "Sample"]
    pr_tracts = pr.from_dict(tract_file)
    pr_tracts = pr_tracts.merge(strand=False, by='Sample')
    return pr_tracts


def partition_into_partitions(lst, num_partitions):
    partition_size = math.ceil(len(lst) / num_partitions)
    partitions = [lst[i:i+partition_size] for i in range(0, len(lst), partition_size)]
    return partitions


def tree_batch_processing(args):

    tree_batches, tree_object, migration_tree_object, bed_list, tgt_id = args

    ts = tree_object.value

    migration_ts = migration_tree_object.value

    for m in migration_ts.migrations():
        # Tree-sequences are sorted by the left ends of the intervals
        # Can skip those tree-sequences are not overlapped with the interval of i.

        for tree_index in tree_batches:
            t = ts.at_index(tree_index)

            if m.left >= t.interval.right: continue
            if m.right <= t.interval.left: break # [l, r)
            for n in ts.samples(tgt_id):
                if t.is_descendant(n, m.node):
                    left = m.left if m.left > t.interval.left else t.interval.left
                    right = m.right if m.right < t.interval.right else t.interval.right

                    bed_list.append([1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}'])

    return True


def _get_truth_tracts(ts, tgt_id, src_id, ploidy):
    """
    Description:
        Helper function to obtain ground truth introgressed tracts at the haploid level from tree-sequence.

    Arguments:
        ts tskit.TreeSqueuece: Tree-sequence containing ground truth introgressed tracts.
        tgt_id str: Name of the target population.
        src_id str: Name of the source population.
        ploidy int: Ploidy of the genomes.

    Returns:
        tracts pandas.DataFrame: Ground truth introgressed fragments from a given source population to a given target populations.
    """
    tracts = pd.DataFrame(columns=['chrom', 'start', 'end', 'sample'])

    src_id = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
    tgt_id = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

    for m in ts.migrations():
        if (m.dest==src_id) and (m.source==tgt_id):
            # For simulations with a long sequence, large sample size, and/or deep generation time
            # This function may become slow
            # May parallelize this function when necessary
            for t in ts.trees():
                # Tree-sequences are sorted by the left ends of the intervals
                # Can skip those tree-sequences are not overlapped with the interval of i.
                if m.left >= t.interval.right: continue
                if m.right <= t.interval.left: break # [l, r)
                for n in ts.samples(tgt_id):
                    if t.is_descendant(n, m.node):
                        left = m.left if m.left > t.interval.left else t.interval.left
                        right = m.right if m.right < t.interval.right else t.interval.right
                        tracts.loc[len(tracts.index)] = [1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}']

    tracts = tracts.sort_values(by=['sample', 'chrom', 'start', 'end'])

    return tracts

### RULES ###


rule all:
    input:
        expand(output_dir + "/{seed}/{output_prefix}.truth.tracts.bed",
               output_prefix=output_prefix, seed=seed_list),
        expand(output_dir + "/{seed}/{output_prefix}.vcf.gz",
               output_prefix=output_prefix, seed=seed_list),


rule simulate_data:
    input:
    output:
        ts = output_dir + "/{seed}/{output_prefix}.ts",
        vcf = output_dir + "/{seed}/{output_prefix}.vcf",
        ref = output_dir + "/{seed}/{output_prefix}.ref.ind.list",
        tgt = output_dir + "/{seed}/{output_prefix}.tgt.ind.list",
    run:
        import demes, msprime
        demo_graph = demes.load(demes_file)
        demography = msprime.Demography.from_demes(demo_graph)
        samples = [
            msprime.SampleSet(nref, ploidy=ploidy, population=ref_id),
            msprime.SampleSet(ntgt, ploidy=ploidy, population=tgt_id),
        ]

        ts = msprime.sim_ancestry(
            recombination_rate=rec_rate,
            sequence_length=seq_len,
            samples=samples,
            demography=demography,
            record_migrations=True,
            random_seed=wildcards.seed,
        )
        ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=wildcards.seed, model=msprime.BinaryMutationModel())
        ts.dump(output.ts)
        with open(output.vcf, 'w') as o: ts.write_vcf(o)

        with open(output.ref, 'w') as f:
            for i in range(nref):
                f.write(f'tsk_{i}\n')

        with open(output.tgt, 'w') as f:
            for i in range(i+1, nref + ntgt):
                f.write(f'tsk_{i}\n')


rule compress_vcf:
    input:
        vcf = rules.simulate_data.output.vcf,
    output:
        vcf = output_dir + "/{seed}/{output_prefix}.vcf.gz",
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {input.vcf}
        """


rule get_truth_tracts:
    input:
        ts = rules.simulate_data.output.ts,
    output:
        bed = output_dir + "/{seed}/{output_prefix}.truth.tracts.bed",
    resources:
        partition="himem", time=120, nodes=1, mem_gb=2000, cpus=128,
    run:


        ts = tskit.load(input.ts)

        src_name = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
        tgt_name = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

        tracts = pd.DataFrame(columns=['chrom', 'start', 'end', 'sample'])
        ploidy=2


        #now we create reduced tree sequence objects
        ts_dump_mig = ts.dump_tables()
        migtable = ts_dump_mig.migrations
        from copy import deepcopy
        migtable2 = deepcopy(migtable)
        migtable2.clear()

        #we search for all rows involving source and target
        for mrow in migtable:
            if (mrow.dest==src_name) and (mrow.source==tgt_name):
                migtable2.append(mrow)


        #the new tree sequence stores only the relevant migrations (involving source and target)
        ts_dump_mig.migrations.replace_with(migtable2)
        ts_dump_sequence_mig = ts_dump_mig.tree_sequence()

        #in the other tree sequence, we delete all migration events
        ts_dump = ts.dump_tables()
        ts_dump.migrations.clear()
        ts_dump_sequence = ts_dump.tree_sequence()

        #we search for all nodes involving the relevant populations
        populations_not_to_remove = [src_name, tgt_name]
        individuals_not_to_remove = []
        for ind in ts.nodes():
            if ind.population in populations_not_to_remove:
                individuals_not_to_remove.append(ind.id)

        #the tree sequence object without migrations can be simplified
        #the simplification contains only the relevant (source-target involving) information
        ts_dump_sequence_simplified = ts_dump_sequence.simplify(individuals_not_to_remove, filter_populations=False, filter_individuals=False, filter_sites=False, filter_nodes=False)

        tracts.columns=["Chromosome", "Start", "End", "Sample"]

        #we can either use concurrent.futures or multiprocessing library
        from concurrent.futures import ProcessPoolExecutor as Pool
        import multiprocessing as mp

        #create manager
        manager = mp.Manager()

        #batches
        tree_inds = list(range(ts_dump_sequence_simplified.num_trees))

        #in case we also want to batch over migrations
        migration_inds = list(range(ts_dump_sequence_mig.num_migrations))


        #managers to share both tree sequences
        #manager_trees = manager.Value(tskit.trees.TreeSequence, ts)
        manager_trees = manager.Value(tskit.trees.TreeSequence, ts_dump_sequence_simplified)

        manager_migrations = manager.Value(tskit.trees.TreeSequence, ts_dump_sequence_mig)

        #manager to write shared file
        bed_list_tree_batches = mp.Manager().list()

        #version for tree batch iteration
        #args = [(curr_index, manager_trees, manager_migrations, bed_list_tree_batches, tgt_name) for curr_index in tree_inds_batches]


        #initialization of multiprocessing pool
        #pool = Pool()
        pool = mp.Pool()

        if len(migration_inds) > 0:
            tree_inds_batches =  partition_into_partitions(tree_inds, mp.cpu_count())
            migration_inds_batches =  partition_into_partitions(migration_inds, mp.cpu_count())


            args = [(manager_trees, curr_index, manager_migrations, bed_list_tree_batches, tgt_name) for curr_index in migration_inds_batches]

            #all_done = pool.map(tree_batch_processing, args)
            #migrations
            all_done =pool.map(migration_batch_processing, args)

            all_done = list(all_done)
            #for concurrent.futures close and join do not exist, so in this case one would have to comment the next two lines
            pool.close()
            pool.join()


            #computation of tracts accomplished

            #now we have to convert the shared list into a pd dataframe
            list_bed_list_tree_batches =  list(bed_list_tree_batches)
            bed_list_tree_df = pd.DataFrame(list_bed_list_tree_batches)
            bed_list_tree_df = create_final_tracts(bed_list_tree_df)
            bed_list_tree_df = bed_list_tree_df.df

            #writing as usual
            open(output.bed, 'w').close()
            for s in bed_list_tree_df['Sample'].unique():
                #in case the simplified tree is used, there could be individuals set to -1
                if "-" not in s:
                    sample_tracts = bed_list_tree_df[bed_list_tree_df['Sample'] == s]
                    sample_tracts = pybedtools.BedTool.from_dataframe(sample_tracts).sort().merge().to_dataframe()
                    sample_tracts['Sample'] = s
                    sample_tracts.to_csv(output.bed, sep="\t", mode='a', header=False, index=False)

        else:
            print("no migration events found")
            df = pd.DataFrame()
            df.to_csv(output.bed, sep="\t", mode='a', header=False, index=False)



#rule obtain_ind_truth_tracts:
#    input:
#        bed = rules.simulate_data.output.bed,
#    output:
#        bed = output_dir + "/{seed}/0/{output_prefix}.0.ind.truth.tracts.bed",
#    run:
#        import pandas as pd
#        import pybedtools
#
#        try:
#            truth_tracts = pd.read_csv(input.bed, sep="\t", header=None)
#        except pd.errors.EmptyDataError:
#            open(output.bed, 'w').close()
#        else:
#            truth_tracts.columns = ['chrom', 'start', 'end', 'sample']
#            truth_tracts['sample'] = truth_tracts.apply(lambda row: "_".join(row['sample'].split("_")[0:2]), axis=1)
#
#            open(output.bed, 'w').close()
#            for s in truth_tracts['sample'].unique():
#                sample_tracts = truth_tracts[truth_tracts['sample'] == s]
#                sample_tracts = pybedtools.BedTool.from_dataframe(sample_tracts).sort().merge().to_dataframe()
#                sample_tracts['sample'] = s
#                sample_tracts.to_csv(output.bed, sep="\t", mode='a', header=False, index=False)
