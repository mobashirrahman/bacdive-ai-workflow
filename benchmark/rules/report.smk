"""Report: figures plus publication of final tables and figures."""


rule report:
    input:
        tables=RESULT_TABLES,
    output:
        figures=RESULT_FIGURES,
        docs_tables=DOC_TABLES,
        docs_figures=DOC_FIGURES,
    log:
        "logs/report.log",
    conda:
        "../envs/stats.yaml"
    params:
        repo=str(REPO),
    script:
        "../scripts/report.py"
