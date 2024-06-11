rule plot:
    input:
        performance = rules.summary.output.summary,
    output:
        pdf = "results/plots/{output_prefix}.performance.pdf",
    log:
        "logs/plot/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    script:
        "../../scripts/plot.py"
