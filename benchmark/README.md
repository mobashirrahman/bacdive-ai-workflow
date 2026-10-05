# Benchmark: how far can BacDive-AI predictions be trusted?

Analyses on selected isolate genomes: leakage-controlled accuracy (A),
annotation drift vs published features (B), completeness degradation (C),
annotation robustness (D). Design and results in `docs/BENCHMARK.md`.

```bash
snakemake -s benchmark/Snakefile --directory <data_root>/run \
  --workflow-profile benchmark/profiles/local --sdm conda --cores 16 \
  --resources ips=4
```

InterProScan runs offline and version-pure (`-dp`); `--resources ips=4`
caps concurrent InterProScan jobs at 4 (each uses threads 4).

Config: `benchmark/config/config.yaml`. Selection `genomes.tsv` and
`label_mapping.yaml` are frozen and committed; logic lives in
`bacdive_workflow/bench/` with tests in `tests/test_bench.py`.

Offline check (real models, fixture annotation, ~1 min):

```bash
R=$(mktemp -d) && benchmark/.test/stage.sh "$R"
snakemake -s benchmark/Snakefile --configfile benchmark/.test/config.yaml \
  --directory "$R" --cores 2 --sdm conda
```

InterProScan runs offline and version-pure (`-dp`); at most 4 run at
once (`resources: ips=1`, threads 4). Seeds derive from config via
SHA-256; same config in, same bytes out.
