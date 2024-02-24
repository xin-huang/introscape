rule plot:
    input:
        performance = rules.summary.output.summary,
    output:
        png = "results/plots/{output_prefix}.performance.png",
    log:
        "logs/plot/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
    script:
        "../scripts/plot.py"
