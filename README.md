# ToMMo NMR Quantification

A pipeline for quantifying **50 plasma metabolites** from NMR BIN spectra and producing a curated concentration table in µM.

**BIN spectra → metabolite-specific features → concentration estimates → designated batch adjustments → 10SD exclusion → low-value exclusion.**

[日本語](README.ja.md) · [Run instructions](docs/setup.md) · [All 50 model definitions](config/model_columns.tsv) · [Feature and batch adjustments](docs/curation.md)

## Quantification

Each input row contains a sample ID and 37,000 spectral intensities on the supplied [ppm axis](reference/BinTable_Header.tsv). Input spectra must use the same axis and intensity scale as the model inputs.

For each metabolite, the pipeline selects the specified spectral region and applies its corresponding model:

- **Linear regression for 46 metabolites**, using a BIN intensity sum, peak height, or an adjusted peak-height feature.
- **Multilayer perceptrons (MLPs) for four metabolites**: methionine, glycerol, valine and alanine, using the selected spectral intensities as input.

For a linear model, `concentration = slope × feature + intercept`. The BIN intensity sum is calculated without an additional ppm-width factor. Where required, feature extraction accounts for reference signals, peak displacement or overlapping signals. Models are applied with the supplied coefficients or weights; prediction does not retrain them.

The [model table](config/model_columns.tsv) lists the method, input columns, ppm endpoints, feature and coefficients or weight file for every metabolite. TSV columns are numbered from **1, including the sample ID column**, and both region endpoints are included. Executable definitions are stored in [models.json](config/models.json).

## Data curation

Selected batches can receive configured spectral-region adjustments or additive median corrections. Batch assignments and settings are supplied explicitly; see [adjustments](docs/curation.md).

After any designated concentration correction, the pipeline applies these rules in order:

1. Set values strictly outside the metabolite's fixed reference **mean ± 10 × SD** to `NA`. Values at the boundary are retained.
2. Set remaining values at or below the metabolite's [operational reporting cutoff](config/reporting_cutoffs.json) to `NA`.

The [reference table](config/reference_stats.tsv) provides the means and SDs; they are not estimated from the incoming samples. Intermediate tables and exclusion flags record each stage.

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

See [run instructions](docs/setup.md) for required inputs and optional batch corrections. The example uses synthetic spectra. Calculations use float64; rounding is confined to a separate display table.

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

Source code and bundled model parameters and weights are provided under the [MIT License](LICENSE). Copyright (c) 2026 ToMMo NMR Quantification Project Contributors. Cohort participant data remain subject to their separate data-access terms.
