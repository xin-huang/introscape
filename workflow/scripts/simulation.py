import sys
sys.stderr = open(snakemake.log[0], "w")

import demes, msprime


nref = int(snakemake.wildcards.nref)
ntgt = int(snakemake.wildcards.ntgt)
seed = int(snakemake.wildcards.test_seed)

demo_graph = demes.load(snakemake.input.demes_file)
demography = msprime.Demography.from_demes(demo_graph)

samples = [
    msprime.SampleSet(nref, ploidy=snakemake.params.ploidy, population=snakemake.params.ref_id),
    msprime.SampleSet(ntgt, ploidy=snakemake.params.ploidy, population=snakemake.params.tgt_id),
]

ts = msprime.sim_ancestry(
    recombination_rate=snakemake.params.rec_rate,
    sequence_length=snakemake.params.seq_len,
    samples=samples,
    demography=demography,
    record_migrations=True,
    random_seed=seed,
)

ts = msprime.sim_mutations(ts, rate=snakemake.params.mut_rate, 
                           random_seed=seed, model=msprime.BinaryMutationModel())
ts.dump(snakemake.output.ts)

with open(snakemake.output.vcf, 'w') as o: ts.write_vcf(o)

with open(snakemake.output.ref, 'w') as f:
    for i in range(nref):
        f.write(f'tsk_{i}\n')

with open(snakemake.output.tgt, 'w') as f:
    for i in range(i+1, nref + ntgt):
        f.write(f'tsk_{i}\n')
