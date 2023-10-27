# Europe_4p1000
Implementation of multiple SOC models statistically calibrated to predict C input requirements to reach a 4‰ SOC stock increase at the European scale.
##Project
Project founded by CLAND (https://cland.lsce.ipsl.fr/) and Holisoils (https://holisoils.eu/)
##Summary
Several SOC models are run over the LUCAS cropland sites (https://esdac.jrc.ec.europa.eu/projects/lucas) to estimate the C input requirements to reach a 4‰ SOC stock increase target in Europe (https://4p1000.org/).
Model parameters are estimated using a statistical approach that links calibrated parameters from Bruni et al. (2022) to environmental variables.

The repository includes:
- Model scripts [Roth-C (Coleman and Jenkinson, 1996), ICBM (Andren and Katterer, 1997), AMG (Andriulo et al., 1999)]
- Scripts to perform the multiple linear regressions for the estimation of model parameters
- Input data (C input forcings and LUCAS soil data for the cropland sites). Climate forcing was derived from the ISIMIP open repository (https://www.isimip.org/outputdata/isimip-repository/)

For questions, comments, or inquieries, please contact Elisa Bruni: ebruni93@gmail.com
