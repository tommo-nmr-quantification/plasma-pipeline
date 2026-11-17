# Feature extraction and batch adjustments

The [model table](../config/model_columns.tsv) describes all 50 metabolites using the same columns. `linear_regression` applies the listed slope and intercept to a scalar feature; `multilayer_perceptron` applies the supplied network weights to a spectral vector.

| Feature | Calculation |
|---|---|
| `sum` | Sum of BIN intensities within the specified columns, without an additional ppm-width factor. |
| `max` | Maximum BIN intensity within the specified columns. |
| `max_minus_reference` | Peak height minus a scaled reference-region peak height. |
| `adjusted_peak_height` | Reference-guided peak detection and adjustment for overlapping signals. |
| `vector` | Selected spectral intensities supplied to an MLP. |

For scalar features, `concentration = slope × feature + intercept`. Input regions, coefficients and weight files are recorded per metabolite.

## Signal adjustments

Succinate uses a reference-region contribution in its feature definition. N,N-dimethylglycine uses a reference peak position selected by inspecting spectra for each processing batch. Peak-pattern detection and overlap adjustment then produce a peak-height feature for one shared linear regression. The reference setting is supplied as an input; subsequent feature extraction and prediction are automated.

## Batch-specific settings

For designated processing groups, additional settings can address shifted peaks or differences in the concentration distribution:

| Setting | Application |
|---|---|
| Spectral-region override | Change the selected columns for the specified metabolite and batch. |
| Additive median correction | Apply a fixed concentration offset estimated from designated reference and target groups; threonine is an example in this cohort workflow. |

The additive correction is:

`corrected = predicted + median(reference group) − median(target group)`

Offsets are estimated from the designated calibration groups and stored before application. Batch membership and correction settings are supplied explicitly. This adjustment assumes that the designated median difference represents a technical offset; it is not estimated automatically from each new input subset.

Batch corrections precede reference-distribution 10SD exclusion and operational low-value exclusion. Sample-to-batch mappings and applicable processing settings accompany cohort data provided through the data-request process. See [run instructions](setup.md) for the configuration interface.
