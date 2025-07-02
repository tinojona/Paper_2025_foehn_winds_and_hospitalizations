# The risk of hospitalization associated with foehn winds and heat in the mountainous region of Switzerland

Authors: Tino Schneidewind [1,2], Sujung Lee [1,2], Ana Maria Vicedo-Cabrera [1,2], Apolline Saucy [1,2]

1 Institute for Social and Preventive Medicine, University of Bern, Switzerland

2 Oeschger Centre for Climate Change Research, University of Bern, Switzerland

<br>
<!--
## Abstract
Foehn winds are intense warm winds, common in mountain regions, but their health impacts and potential to exacerbate existing heat-related risks remain poorly understood. We investigated the independent and combined association of foehn winds and temperature with cause-specific emergency hospitalizations in Switzerland. We found that foehn winds daily intensity showed small and no consistent association with hospitalizations in temperature-adjusted and non-adjusted models. However, foehn winds amplified heat-related hospitalization risk with a 14% increase in risk at the 99th temperature percentile on foehn days, compared to -2% on non-foehn days. The association was larger for females, older adults, and for hospitalizations due to respiratory and mental health causes. While foehn winds did not directly impact hospitalizations, they may contribute to an amplification of heat-related health risks, especially for females and older adults.
<p align="center">
  <img src="output/figures/Figure4.png" alt="Figure 1" width="500"/>
</p>
*(left) Cumulative relative risk (Model 4) for subpopulations at -8.9°C (1st percentile) with 95% confidence intervals and (right) cumulative relative risk (Model 4) for subpopulations at 24.7°C (99th percentile) with 95% confidence intervals.*
-->

<br>

## Aim of this repository

This `repository` presents the workflow and the analysis of the paper, which developed out of my master thesis. 

<br>

## Data sensitivity issues 

Due to the sensitivity of the hospitalization data, I cannot publish the full data of the analysis. The data that is not part of this repository but can be applied for using the following procedure. The hospitalization data was provided by the [Swiss Federal Office for Statistics](https://www.bfs.admin.ch/bfs/de/home.html) and the data can be accessed as described on their [webpage](https://www.bfs.admin.ch/bfs/de/home/statistiken/gesundheit/erhebungen/ms.html). The non-sensitive meteorological data was provided by the [Swiss Federal Office for Meteorology and Climatology](https://www.meteoschweiz.admin.ch/#tab=forecast-map) and is part of this repository.

<br>

## Folders

### The analysis folder

The `analysis` folder contains the crossbasis definition methods for both foehn winds and temperature exposure (A01, A02), the procedure by which the descriptive statistics were calculated (B01, B02) and both research questions were answered (C01, C02).


### The data folder

The `data` folder includes all the environmental data in both raw and processed form, as well as the necessary shapefiles for the buffer size calculation and subsequent Medstat region selection.


### The functions folder

The `functions` folder contains R functions, which are stored there so that they can be accessed in other R scripts.


### The output folder

The `output` folder has both all figures and tables of the paper that were created in the analysis. 


### The vignettes folder

The `vignettes` folder presents how the buffersize influences the Medstat selection (01) and an explorative investigation into the foehn wind data (02), the hospitalization data (03) and for both Altdorf (04) and Lugano (05) an overview over the processed hospitalization-foehnwind-temperature data set.
