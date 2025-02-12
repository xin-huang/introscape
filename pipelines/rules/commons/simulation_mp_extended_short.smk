import os

#-----------------------------------------------------------------------------------------------------------------------
# functions needed for rule get_truth_tracts - moved to pipelines/rules/commons/simulation_mp.py
# run section of rule get_truth_tracts - moved to pipelines/rules/commons/simulation_mp.py
#-----------------------------------------------------------------------------------------------------------------------

### RULES ###

'''
rule all:
    input:
        expand(output_dir + "/{seed}/{output_prefix}.truth.tracts.bed",
               output_prefix=output_prefix, seed=seed_list),
        expand(output_dir + "/{seed}/{output_prefix}.vcf.gz",
               output_prefix=output_prefix, seed=seed_list),
'''

rule simulate_data:
    output:
        ts = output_dir + "/{seed}/{output_prefix}.ts",
        vcf = output_dir + "/{seed}/{output_prefix}.vcf",
        ref = output_dir + "/{seed}/{output_prefix}.ref.ind.list",
        tgt = output_dir + "/{seed}/{output_prefix}.tgt.ind.list",
    log:
        "logs/sim.{seed}.{output_prefix}.log",
    benchmark:
        "benchmarks/sim.{seed}.{output_prefix}.benchmark.txt",
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

        if binary:
            ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=wildcards.seed, model=msprime.BinaryMutationModel())
        else:
            ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=wildcards.seed)

        ts.dump(output.ts)
        with open(output.vcf, 'w') as o: ts.write_vcf(o)

        if not binary:
            new_vcf = output_dir + f"/{wildcards.seed}/{output_prefix}.vcf.recode.vcf"
            shell("vcftools --vcf {output.vcf} --min-alleles 2 --max-alleles 2 --recode --out {output.vcf}")
            shell("mv {new_vcf} {output.vcf}")

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

rule mp_get_truth_tracts:
    input:
        ts = rules.simulate_data.output.ts,
    output:
        bed = output_dir + "/{seed}/{output_prefix}.truth.tracts.bed",
    log:
        "logs/tract.{seed}.{output_prefix}.log",
    benchmark:
        "benchmarks/tract.{seed}.{output_prefix}.benchmark.txt",
    params:
        tgt_id = config["tgt_id"],
        src_id = config["src_id"],
        ploidy = config["ploidy"]
    resources:
        partition="himem", time=120, nodes=1, mem_gb=2000, cpus=128,
    shell:
        """    
        python pipelines/rules/commons/simulation_mp.py \
         --inputts '{input.ts}' --tgtid '{params.tgt_id}' \
         --srcid '{params.src_id}'  --ploidy {params.ploidy} --output '{output.bed}'
        """
