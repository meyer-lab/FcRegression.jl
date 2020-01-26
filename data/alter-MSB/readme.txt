Input data for G. Alter et al., "High-resolution definition of humoral correlates of effective immunity against HIV", Molecular Systems Biology, 2018

data: the experimental data
* for each, the rows are subjects (id in first column) and the columns are features
* some subjects are entirely missing and are not included at all 
* "NA" indicates a single missing measurement
* files (Table S1 provides details)
  - functions: effector function assays
  - glycan-gp120: gp120-specific glycan abundances
  - luminex-igg: Fc Array, IgG detection reagents
  -  luminex: Fc Array, non-IgG detection reagents

meta: metadata for analysis
- antigens: groups of Fc Array antigen names
- colors: names of colors to use when plotting various properties
- detections: groups of Fc Array detection reagent names
- glycans: summary characteristics of each glycan feature: _b_isected GlcNAc, _s_ialylated, level of _g_alactosylation, and _f_ucosylated
- subjects: groups by clinical status, with original 4-class labels (EC/VC/TP/UP) also broken into two 2-class labels (controller/progressor, viremic/nonviremic)
