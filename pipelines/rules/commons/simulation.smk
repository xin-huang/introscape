import numpy as np

feature_id = config["feature_id"]
feature_config = config["feature_config"]
nfeature = config["nfeature"]
nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
geno_state_list = config["geno_states"]
output_prefix = config["output_prefix"]
win_step = config["win_step"]
cutoff_num = config["cutoff_num"]
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])

nrep = {}
seq_len = {}
demog_id = {}
demes = {}
mut_rate = {}
rec_rate = {}
ref_id = {}
tgt_id = {}
src_id = {}
seed_list = {}
output_dir = {}

for k in ["train", "test"]:
    np.random.seed(config["seed"])
    nrep[k] = config["nrep"][k]
    seq_len[k] = config["seq_len"][k]
    demog_id[k] = config["demog_id"][k]
    demes[k] = config["demes"][k]
    mut_rate[k] = config["mut_rate"][k]
    rec_rate[k] = config["rec_rate"][k]
    ref_id[k] = config["ref_id"][k]
    tgt_id[k] = config["tgt_id"][k]
    src_id[k] = config["src_id"][k]
    seed_list[k] = np.random.randint(1, 2**31, nrep[k])
    output_dir[k] = f"results/data/{k}/{demog_id[k]}/nref_{nref}/ntgt_{ntgt}"

test_seed = seed_list["test"]
output_dir = output_dir["test"]

rule all:
    input:
        expand(output_dir + "/{seed}/{output_prefix}.phased.true.tracts.bed",
               output_prefix=output_prefix, seed=test_seed),
        expand(output_dir + "/{seed}/{output_prefix}.vcf.gz",
               output_prefix=output_prefix, seed=test_seed),


rule simulate_test_data:
    input:
        demes_file = demes["test"],
    output:
        ts = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.ts",
        vcf = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.vcf",
        ref = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.ref.ind.list",
        tgt = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.tgt.ind.list",
    log:
        "logs/simulate_test_data/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    params:
        ploidy = ploidy,
        seq_len = seq_len["test"],
        mut_rate = mut_rate["test"],
        rec_rate = rec_rate["test"],
        ref_id = ref_id["test"],
        tgt_id = tgt_id["test"],
        src_id = src_id["test"],
    script:
        "../../scripts/simulation.py"


rule compress_vcf:
    input:
        vcf = rules.simulate_test_data.output.vcf,
    output:
        vcf = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.vcf.gz",
    log:
        "logs/compress_vcf/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        rm {input.vcf} 2>> {log}
        """


rule get_phased_true_tracts:
    input:
        ts = rules.simulate_test_data.output.ts,
    output:
        bed = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.phased.true.tracts.bed",
    log:
        "logs/get_phased_true_tracts/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.log",
    params:
        ploidy = ploidy,
        tgt_id = tgt_id["test"],
        src_id = src_id["test"],
    resources:
        partition = "himem", 
        time = 60,
        mem_gb = 2000,
        cpus = 128,
    run:
        import tskit
        import pyranges as pr
        from multiprocessing import Process, Manager

        def worker_func(in_queue, out_queue, **kwargs):
            while True:
                trees, migration = in_queue.get()
                tgt_name = kwargs['tgt_name']
                ts = kwargs['ts']
                ploidy = kwargs['ploidy']

                res = ''

                try:
                    for t in trees:
                        for n in ts.samples(tgt_name):
                            if t.is_descendant(n, migration.node):
                                 left = migration.left if migration.left > t.interval.left else t.interval.left
                                 right = migration.right if migration.right < t.interval.right else t.interval.right
                                 res += f'1\t{int(left)}\t{int(right)}\ttsk_{ts.node(n).individual}_{int(n%ploidy+1)}\n'
                    out_queue.put(res)
                except Exception as e:
                # Handle or log the exception as needed
                    print(f"Error in worker: {e}")

        ts = tskit.load(input.ts)

        src_name = [p.id for p in ts.populations() if p.metadata['name']==params.src_id][0]
        tgt_name = [p.id for p in ts.populations() if p.metadata['name']==params.tgt_id][0]

        res = "Chromosome\tStart\tEnd\tSample\n"
        with Manager() as manager:
            in_queue  = manager.Queue()
            out_queue = manager.Queue()
            keywords = {'tgt_name': tgt_name, 'ts': ts, 'ploidy': params.ploidy}
            workers = [
                Process(target=worker_func, args=(in_queue, out_queue), kwargs=keywords) for i in range(resources.cpus)
            ]

            num_introgression = 0

            for m in ts.migrations():
                if (m.dest==src_name) and (m.source==tgt_name):
                    in_queue.put((ts.trees(left=m.left, right=m.right), m))
                    num_introgression += 1

            for w in workers: w.start()

            try:
                for i in range(num_introgression):
                    item = out_queue.get()
                    res += item
                for w in workers: w.terminate()
            except Exception as e:
                print(f"Error in manager: {e}")

        res = pr.from_string(res)
        res = res.merge(strand=False, by='Sample')
        if not res.empty: res.to_csv(output.bed, sep="\t", header=False)
        else: open(output.bed, 'w').close()


rule get_unphased_true_tracts:
    input:
        bed = rules.get_phased_true_tracts.output.bed,
    output:
        bed = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.unphased.true.tracts.bed",
    log:
        "logs/get_unphased_true_tracts/{test_demog}/nref_{nref}/ntgt_{ntgt}/{test_seed}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    script:
        "../../scripts/get_unphased_true_tracts.py"
