# Contributing

Use Python 3.12 and install the project with `python -m pip install -e '.[dev,workflow]'`.

Before submitting changes:

```bash
ruff check bacdive_workflow tests scripts predict.py
ruff format --check bacdive_workflow tests scripts predict.py
pytest -q
python scripts/reproduce_case_study.py --check
```

For workflow changes, also prepare the upstream models, run the real annotation
example with Snakemake, and run `python scripts/check_example.py
results/predictions/1120941.3.json`. Explain any changed predictions and record
the model/database/environment versions. Do not silently change the published
feature encoding, E-value threshold, or reference-label interpretation.

Bug reports should include the command, package versions, relevant stderr, and
a minimal input that can be shared. Do not include credentials, private paths,
or unrelated genome data. Pull requests should explain the user-visible problem
and include regression coverage for changed behavior.
