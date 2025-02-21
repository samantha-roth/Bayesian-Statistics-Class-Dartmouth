## Overview
This codes includes a toy example for implementing a workflow to
- fit a model to data
- perform a positive control test
- perform a method comparison (Frequentist Bootstrap vs Bayesian calibration via MCMC) 
- assess convergence and 
- perform two methods for sensitivity analysis

This code also includes an example for implementing a workflow to compare Bayesian calibration via MCMC to precalibration.

## How to run:
1. Save the file in a directory
2. Open R-Studio and navigate to the folder with the file
3. Make this directory the work directory
4. Open the file
5. Source the file

## Main document
Workflow_example_SR.R

## Additional examples
- precalibration_vs_mcmc_example.R is not included in the main document. This code compares Bayesian calibration via MCMC to precalibration.
- labX_SR.R is also not included in the main document. This code illustrates how to identify lack of convergence by considering multiple additional convergence metrics.
- 
## Contributors
Atieh Alipour (atieh.alipour@dartmouth.edu), Klaus Keller (Klaus@dartmouth.edu), Haochen Ye (hxy46@psu.edu), Prabhat Hegde (Prabhat.Hegde.TH@dartmouth.edu), Sitara Baboolal (Sitara.D.Baboolal@dartmouth.edu), and Samantha Roth (samantha.m.roth@dartmouth.edu)

## Versions
Format: Version / last changes / by whom / what
1. Jan 20 2023 / Contributors from above / write initial version
2. Dec 27 2023 / Klaus Keller / update header, clean up code, rename file
3. Jan 25 2024 / Klaus Keller / clean up comments and format
4. Feb 05 2024 / Klaus Keller / add comments
5. Feb 21 2025 / Samantha Roth / add comments, labX_SR.R, precalibration_vs_mcmc_example.R
