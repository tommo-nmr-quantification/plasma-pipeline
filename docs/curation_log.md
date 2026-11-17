# Metabolite-specific curation record

Automated quantification results were manually reviewed for metabolite-specific and batch-associated patterns. This review included spectral regions and concentration distributions, with particular attention to differences in sample buffers and pH-sensitive peak positions. Where needed, expert-guided adjustments were applied to designated processing groups. Buffer and pH differences were considered during review; the entries below describe the actions taken without assigning every observed difference to a single cause.

The record summarizes standard feature adjustments and selected-group curation during dataset preparation. It complements the [current model catalogue](../config/model_columns.tsv) and [executable adjustment settings](curation.md); it is not a list of additional corrections to apply to every sample.

| Metabolite | Review focus | Recorded action |
|---|---|---|
| TCN000001 — Succinate | Contribution from a reference signal | Subtract one fifth of the reference-region peak height before regression; included in the standard feature definition. |
| TCN000004 — Formate | Individual quantification results | Set review-selected concentration values to missing. |
| TCN000014 — Cysteine | Placement and coverage of the integration region | Revise the region position and subsequently expand the integration interval. |
| TCN000015 — Glutamate | Signal capture in a designated processing group | Use an alternative integration region for that group. |
| TCN000017 — Asparagine | Concentration calibration and selected low values | Recalibrate estimates in an earlier processing workflow and set review-selected values to missing. The earlier recalibration is not added to the current regression. |
| TCN000021 — Ornithine | Placement of the peak-height region | Revise the spectral interval used to obtain the maximum intensity. |
| TCN000025 — Histidine | Peak selection and individual quantification results | Use an alternative peak-height region for selected groups and apply review-based exclusions, including a rule linked to missing Asparagine values. |
| TCN000028 — Serine | Signal capture and unusually high concentrations | Shift the integration region for a designated group; identify high values for review, with missing-value exclusion in a release-specific curation step. |
| TCN000030 — Threonine | Differences in concentration distributions between processing groups | Apply an additive offset equal to the reference-group median minus the target-group median. |
| TCN000035 — Glutamine | Peak selection within the reference region | Narrow the peak-height interval for a designated group. |
| TCN000042 — N,N-Dimethylglycine | Peak displacement and overlapping signals | Select a reference peak position manually for each processing batch, then automatically adjust the peak-height feature and apply one shared regression. |
| TCN000043 — Uridine | Placement of the peak-height region in a designated group | Use an alternative peak-height interval for that group. |

Group-specific region changes are made during feature extraction. Designated concentration corrections precede the common post-quantification rules: exclusion strictly outside the fixed reference mean ± 10 SD, followed by exclusion at or below the metabolite-specific [operational reporting cutoff](../config/reporting_cutoffs.json). These common rules are implemented in [postprocess.py](../nmr_cohort/postprocess.py).

The current catalogue defines the default models. Group-specific actions are associated with the relevant processing groups and dataset releases. Sample-to-group mappings, applicable settings and review decisions are maintained with the cohort data and supplied through its data-request process. The standard pipeline does not infer individual review decisions from sample identifiers.
