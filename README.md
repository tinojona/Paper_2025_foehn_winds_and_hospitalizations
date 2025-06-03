# Master Thesis: *The risk of hospitalization associated with foehn winds and heat in the mountainous region of Switzerland*

author: Tino Schneidewind

Disclamer:

This repository presents my workflow and my results of my master thesis. Due to the sensitivity of the hospitalization data, I cannot publish the full data of the analysis. The data that is not part of this repository but can be applied for using the following procedure. The meteorological data was provided by the [Swiss Federal Office for Meteorology and Climatology](https://www.meteoschweiz.admin.ch/#tab=forecast-map) and is accessable through their own distribution platform [IDAweb](https://www.meteoschweiz.admin.ch/service-und-publikationen/service/wetter-und-klimaprodukte/datenportal-fuer-lehre-und-forschung.html). The hospitalization data was provided by the [Swiss Federal Office for Statistics](https://www.bfs.admin.ch/bfs/de/home.html) and the data can be accessed as described on their [webpage](https://www.bfs.admin.ch/bfs/de/home/statistiken/gesundheit/erhebungen/ms.html).



## Folders

### The analysis folder

The `analysis` folder contains the crossbasis definition methods for both foehn winds and temperature exposure (A01, A02). The procedure by which the descriptive statistics were calculated (B01, B02) and both research questions were answered (C01, C02).


### The data folder

The `data` folder includes all the environmental data in both raw and processed form, as well as the necessary shapefiles for the buffer size calculation and subsequent Medstat region selection.


### The functions folder

The `functions` folder contains R functions, which are stored there so that they can be accessed
in other R scripts.


### The output folder

The `output` folder has both all figures and tables of the paper that were created in the analysis. 


### The vignettes folder

The `vignettes` folder presents how the buffersize influences the Medstat selection (01) and an explorative investigation into the foehn wind data (02), the hospitalization data (03) and for both Altdorf and Lugano an overview over the processed hospitalization-foehnwind-temperature data set.
