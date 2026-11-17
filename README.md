# ToMMo NMR Quantification

A common pipeline for converting plasma NMR BIN spectra into a curated concentration table (µM).

**BIN columns → features → frozen models → reference 10SD exclusion → low-value exclusion.**

[日本語](README.ja.md) · [Run instructions](docs/setup.md) · [All 50 column/model definitions](config/model_columns.tsv) · [Compound-specific adjustments](docs/curation.md)

## Calculation

Every input row contains an ID followed by 37,000 spectral values on the supplied [ppm axis](reference/BinTable_Header.tsv). All compounds use inclusive start/end BIN indices from one [model catalogue](config/models.json):

```python
region = spectrum[bin_start:bin_end + 1]
feature = max(region)                 # peak height
# or feature = serial_sum(region)    # BIN integral, without a ppm multiplier
concentration = slope * feature + intercept
```

The [column table](config/model_columns.tsv) gives the actual first/last TSV columns, ppm endpoints, feature and coefficients for every compound. **TSV columns are numbered from 1, including the ID column.** In the code, the spectral vector excludes the ID and starts at index 0: `TSV column = BIN index + 2`. Both endpoints are included.

There are 45 ordinary scalar linear models, four frozen MLPs (TCN000007/24/32/34), and the TCN000042 v6 model. Succinate subtracts a reference-region contribution before regression. TCN000042 uses a batch-specific reference peak to adjust its feature automatically, followed by one fixed regression. MLP weights are included; inference does not retrain models.

## Filtering

After any designated batch correction, apply the following rules **in order**:

1. Set a concentration to `NA` when it is strictly outside the fixed reference `mean ± 10 × SD` for that metabolite. Values at the boundary are retained.
2. Set a remaining concentration to `NA` when it is at or below the metabolite's [operational reporting cutoff](config/reporting_cutoffs.json).

The default [reference table](config/reference_stats.tsv) supplies all 50 means and SDs in µM at the reported two-decimal precision. These stored values define the exclusion bounds; they are not recalculated from the incoming samples. Existing missing values stay missing. Outputs retain the table after 10SD filtering, the final table and per-stage exclusion flags.

## Run

```bash
conda env create -f environment.yml
conda activate tommo-nmr-quantification
bash scripts/run_synthetic_demo.sh
```

For cohort spectra:

```bash
python -m nmr_cohort run \
  --bin data/spectra.tsv --format header \
  --batch-map data/batch_map.tsv --groups data/groups.json \
  --out runs/cohort
```

See [setup](docs/setup.md) for the input tables and optional frozen batch offsets. Synthetic examples contain no cohort observations. Coefficients retain their supplied precision; calculations use float64, with rounding only in the separate display table. Model-development and performance details belong to the accompanying manuscript.

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

Source code and bundled model parameters and weights are provided under the [MIT License](LICENSE). Copyright (c) 2026 ToMMo NMR Quantification Project Contributors. Cohort participant data remain subject to their separate data-access terms.
