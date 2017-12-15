# microneedle-montecarlo
MATLAB code for Monte Carlo simulation of microneedle insertion / capillary strikes

Blicharz et al. "Microneedle-based device for the one-step painless collection of capillary blood samples"
Nature Biomedical Engineering in press


| Parameter | Definition | Values |
| --------- | ---------- | ------ |
| N | Number of trials | 100000 |
| needleN | Number of needles in array | 1,4,8,30 |
| needleXspacing | Spacing between needles in X dimension (mm) | 0.2 |
| needleYspacing | Spacing between needles in Y dimension (mm) | 0.85 |
| needleXlength | Length of needle in X dimension (“width”) (mm) | 0.35 |
| needleYlength | Length of needle in Y dimension (“thickness”) (mm) | 0.05 |
| needleRows | Number of rows in needle array | 1,2,2,6 |
| capillaryDensity | Capillary density (per mm2) | 14 |
| capillaryDiameter | Capillary diameter (mm) | 0.015 |
| maxiter | Maximum number of attempts to generate a valid capillary array | 500 |
| capillaryMinimumDistance | Minimum intercapillary distance (mm) | 0.08 |

