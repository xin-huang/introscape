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
        time = 60, 
        mem_gb = 128, 
        cpus = 16,
    run:
        import tskit
        import numpy as np
        import pyranges as pr
        from multiprocessing import Process, Manager

        def worker_func(in_queue, out_queue, shared_trees, shared_samples, **kwargs):
            while True:
                m_left, m_right, m_node = in_queue.get()
                ploidy = kwargs['ploidy']

                res = ''

                try:
                    for t in shared_trees.value.trees():
                        if m_left >= t.interval.right: continue
                        if m_right <= t.interval.left: break # [l, r)
                        for n in shared_samples.value:
                            if t.is_descendant(n, m_node):
                                 left = max(m_left, t.interval.left)
                                 right = min(m_right, t.interval.right)
                                 res += f'1\t{int(left)}\t{int(right)}\ttsk_{ts.node(n).individual}_{int(n%ploidy+1)}\n'
                    out_queue.put(res)
                except Exception as e:
                # Handle or log the exception as needed
                    print(f"Error in worker: {e}")

        def simplify_ts(ts, tgt_id, src_id):
            from copy import deepcopy

            #now we create reduced tree sequence objects
            ts_dump_mig = ts.dump_tables()
            migtable = ts_dump_mig.migrations
            migtable2 = deepcopy(migtable)
            migtable2.clear()

            for mrow in migtable:
                if (mrow.dest==src_id) and (mrow.source==tgt_id):
                    migtable2.append(mrow)

            ts_dump_mig.migrations.replace_with(migtable2)
            ts_dump_sequence_mig = ts_dump_mig.tree_sequence()

            ts_dump = ts.dump_tables()
            ts_dump.migrations.clear()
            ts_dump_sequence = ts_dump.tree_sequence()

            populations_not_to_remove = [src_id, tgt_id]
            individuals_not_to_remove = []
            for ind in ts.nodes():
                if ind.population in populations_not_to_remove:
                    individuals_not_to_remove.append(ind.id)

            ts_dump_sequence_simplified = ts_dump_sequence.simplify(
                individuals_not_to_remove, 
                filter_populations=False, 
                filter_individuals=False, 
                filter_sites=False, 
                filter_nodes=False
            )

            return ts_dump_sequence_simplified, migtable2

        ts = tskit.load(input.ts)

        src_name = [p.id for p in ts.populations() if p.metadata['name']==params.src_id][0]
        tgt_name = [p.id for p in ts.populations() if p.metadata['name']==params.tgt_id][0]

        tgt_samples = ts.samples(tgt_name)

        ts, migtable = simplify_ts(ts=ts, tgt_id=tgt_name, src_id=src_name)

        res = "Chromosome\tStart\tEnd\tSample\n"
        with Manager() as manager:
            in_queue  = manager.Queue()
            out_queue = manager.Queue()

            num_introgression = 0
            shared_trees = manager.Value(tskit.trees.TreeSequence, ts)
            shared_samples = manager.Value(np.ndarray, tgt_samples)

            keywords = {'tgt_name': tgt_name, 'ts': ts, 'ploidy': params.ploidy}
            workers = [
                Process(target=worker_func, args=(in_queue, out_queue, shared_trees, shared_samples), kwargs=keywords) for i in range(resources.cpus)
            ]


            for m in migtable:
                in_queue.put((m.left, m.right, m.node))
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
