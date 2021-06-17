# DNHiPSE
This R-package is aimed to illustrate the methodology presented in "Modeling the natural history of human diseases".

## Installation
```r
devtools::install_github("eric40065/DNHiPSE")
```

## Usage
The following code can reproduce everything mentioned in our paper. It will take 2 hours to reproduce the all the simulation and analyses of real data application.
```r
library(DNHiPSE)
# To reproduce Table 1
Table11 = DNH_simulation("coverage", "null", sample_size = 300)
Table12 = DNH_simulation("coverage", "null", sample_size = 1000)
Table13 = DNH_simulation("coverage", "alternative", sample_size = 1000)
Table14 = DNH_simulation("coverage", "alternative", sample_size = 1000)

# To reproduce Figure 3
DNH_simulation("unbiasedness", "null", sample_size = 300)
DNH_simulation("unbiasedness", "null", sample_size = 1000)

# To reproduce Figure 4
DNH_simulation("unbiasedness", "alternative", sample_size = 300)
DNH_simulation("unbiasedness", "alternative", sample_size = 1000)

# To reproduce Figure 5
result_HBV = DNHiPSE(DNH_REVEAL_HBV, plot_result = TRUE, bootstrap = TRUE, plot_unit = 365)

# To reproduce Figure 6 and proportion of mediation from Section 7
result_HCV = DNHiPSE(DNH_REVEAL_HCV, plot_result = TRUE, bootstrap = TRUE, plot_unit = 365, PM = TRUE)

# To reproduce Figure S4 of supplement material
result_HBV_supp = DNHiPSE(DNH_REVEAL_HBV, plot_result = TRUE, plot_unit = 365, sensitivity_analysis = TRUE)

# To reproduce Figure S5 of supplement material
result_HCV_supp = DNHiPSE(DNH_REVEAL_HCV, plot_result = TRUE, plot_unit = 365, sensitivity_analysis = TRUE)
```

### REVEAL-HBV and REVEAL-HCV
#### Abstract
- REVEAL is a community-based prospective cohort study conducted in Taiwan.
- we focused on male participants with age at cohort entry less than 55 years and without alcohol consumption (m = 7,782 and 6,305, respectively, for hepatitis B and C analyses).
#### Dictionary
REVEAL-HBV/HCV contain 7 columns.
- T1: time to liver cirrhosis or censored time (days).
- d1: whether or not the liver cirrhosis is observed (0: no vs. 1:yes).
- T2: time to liver cancer incidence or censored time (days).
- d2: whether or not the liver cancer incidence is observed (0: no vs. 1:yes).
- T3: time to death or censored time (days).
- d3: whether or not death is observed (no: 0 vs. yes: 1).
- Z: hepatitis B surface antigen/anti-hepatitis C antibody
