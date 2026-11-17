# Compound-specific processing notes

The standard models share one column/feature catalogue. The following records summarize additional processing choices found in the cohort workflow. Historical or group-specific adjustments are associated with their designated processing versions and sample groups.

| Compound | Recorded adjustment |
|---|---|
| TCN000001 Succinate | Subtract one fifth of a reference-region peak height before regression; included in the standard feature definition. |
| TCN000004 Formate | Exclusion of review-selected observations. |
| TCN000005 2-Aminobutyrate | Updated peak-height calibration; included in the standard catalogue. |
| TCN000014 Cysteine | Revision and expansion of the integration region. |
| TCN000015 Glutamate | Alternative integration region for a selected processing group. |
| TCN000017 Asparagine | Historical concentration recalibration and review-based exclusions in designated groups. |
| TCN000021 Ornithine | Revision of the peak-height region in the main processing route. |
| TCN000025 Histidine | Alternative peak-height region and review-based exclusions for selected groups. |
| TCN000028 Serine | Shifted integration region for a selected group; review of high concentrations. |
| TCN000030 Threonine | Additive median correction between processing groups. |
| TCN000032 Valine / TCN000034 Alanine | MLP predictions selected for the final integrated table. |
| TCN000033 Glycine | Alternative integration region for a selected processing group. |
| TCN000035 Glutamine | Narrowed peak-height region for a selected processing group. |
| TCN000042 N,N-Dimethylglycine | Manually assigned group reference peak, automated interference adjustment and one shared v6 regression. |
| TCN000043 Uridine | Alternative peak-height region for a selected processing group. |

For additive median correction: `corrected = predicted + median(reference group) − median(target group)`. Offsets are frozen from the designated population before application. Only designated groups receive a region change or correction; the standard pipeline does not infer group membership from sample IDs.

The executable workflow applies the current model catalogue, explicitly supplied batch settings, reference 10SD exclusion and operational low-value cutoffs. Review-based decisions and historical recalibrations are recorded here as processing history. Sample-to-batch mappings and the relevant processing settings accompany approved data requests.
