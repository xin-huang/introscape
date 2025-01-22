rule plot:
    input:
        performance = rules.summary.output.summary,
    output:
        pdf = "results/plots/train_{train_demog}_test_{test_demog}/{output_prefix}.performance.pdf",
    log:
        "logs/plot/train_{train_demog}_test_{test_demog}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    script:
        "../../scripts/plot.py"
