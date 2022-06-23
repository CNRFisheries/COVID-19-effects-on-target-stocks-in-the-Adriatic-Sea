# The potential effects of COVID-19 lockdown and the following restrictions on the status of eight target stocks in the Adriatic Sea

This repo accompanies the manuscript _"The potential effects of COVID-19 lockdown and the following restrictions on the status of eight target stocks in the Adriatic Sea"_.

All the required data for catches and biomass are stored in the _data_ folder.
To run the analysis proceed as follow:
* clone the repository on your machine
* clean out the results folder: **_it must be empty_**
* open the R project first
* open the script ```step1_CMSY++12c.R``` and run it twice: once having _final_year_ parameter (line 57) set as 2019 and once having _final_year_ parameter set as 2020.

**NOTE:** the script applies the surplus production model CMSY++, from [Froese et al., 2017](https://enriconarmelloni.github.io/SOLEA/) and also available [at this repo](https://github.com/SISTA16/cmsyPlusPlus), on the selected target stocks in the Adriatic Sea. The algorithm is the original version, slight modifications were applied to the script to allow BSM results extraction. For the purpose of the analysis, two stock assessments per stock are done: one on the data up to 2019 and one on the data up to 2020. To do so the script have to be run twice, changing the _final_year_ parameter (line 57).

* open the script ```step2_handle_estimations.R``` and run it

**NOTE:** the script merge toghether the BSM model estimations for all the stocks to be used by the next script.

* open the script ```step3_bootstrap_analysis.R``` and run it

**NOTE:** the script perform a bootstrap analysis to test difference between BSM model estimations for 2019 and 2020.

* open the script ```step4_manuscript_figures.R``` and run it

**NOTE:** the script extract various results and generates the figures.
