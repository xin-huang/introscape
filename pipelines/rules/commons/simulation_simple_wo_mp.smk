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
        partition="basic", time=2000, mem_gb=100, cpus=16,
    run:

        ts = tskit.load(input.ts)


        src_name = [p.id for p in ts.populations() if p.metadata['name']==src_id][0]
        tgt_name = [p.id for p in ts.populations() if p.metadata['name']==tgt_id][0]

        tracts = pd.DataFrame(columns=['chrom', 'start', 'end', 'sample'])
        ploidy=2


        truth_tracts = _get_truth_tracts(ts, tgt_id, src_id, ploidy)

        print(truth_tracts)

        #computation of tracts accomplished
        open(output.bed, 'w').close()
        for s in truth_tracts['sample'].unique():
            sample_tracts = truth_tracts[truth_tracts['sample'] == s]
            sample_tracts = pybedtools.BedTool.from_dataframe(sample_tracts).sort().merge().to_dataframe()
            sample_tracts['sample'] = s
            sample_tracts.to_csv(output.bed, sep="\t", mode='a', header=False, index=False)




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
