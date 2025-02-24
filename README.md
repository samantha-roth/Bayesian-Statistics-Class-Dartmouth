## Overview
This codes includes a toy example for implementing a workflow to
- fit a model to data
- perform a positive control test
- perform a method comparison (Frequentist Bootstrap vs Bayesian calibration via MCMC) 
- assess convergence and 
- perform two methods for sensitivity analysis

This code also includes an example for implementing a workflow to compare Bayesian calibration via MCMC to precalibration.

## How to run:
1. Download this repository in a zip file or clone it.
2. Unzip the folder corresponding to the respository if needed.
3. Open R-Studio and navigate to the folder corresponding to this repository.
4. Make this directory the work directory.
5. Open each R script.
6. Source each R script.
7. Compare the figures produced by each R script to the figures with the same name and folder in the `compare` folder.

## Main document
Workflow_example_SR.R

## Additional examples
- precalibration_vs_mcmc_example.R is not included in the main document. This code compares Bayesian calibration via MCMC to precalibration.
- Reading.Rmd is an R markdown file that produces ReadingCorrection.pdf.
- Reading.pdf describes an example implementing MCMC and checking for lack of convergence and provides code snippets. This pdf was produced by knitting ReadingCorrection.Rmd to a pdf.

## Contributors
Atieh Alipour (atieh.alipour@dartmouth.edu), Klaus Keller (Klaus@dartmouth.edu), Haochen Ye (hxy46@psu.edu), Prabhat Hegde (Prabhat.Hegde.TH@dartmouth.edu), Sitara Baboolal (Sitara.D.Baboolal@dartmouth.edu), and Samantha Roth (samantha.m.roth@dartmouth.edu).
Kelsey Ruckert, Tony Wong, Yawen Guan, Patrick Applegate, and Murali Haran provided the original code and example in Reading.Rmd. Samantha Roth later edited this example and provided additional code.

## Versions
Format: Version / last changes / by whom / what
1. Jan 20 2023 / Contributors from above / write initial version
2. Dec 27 2023 / Klaus Keller / update header, clean up code, rename file
3. Jan 25 2024 / Klaus Keller / clean up comments and format
4. Feb 05 2024 / Klaus Keller / add comments
5. Feb 21 2025 / Samantha Roth / add comments, labX_SR.R, precalibration_vs_mcmc_example.R
