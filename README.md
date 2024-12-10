# CDS_heat_stress_indicators

This repository provides the equations that were used to calculate heat stress indicators for CMIP6 models. The resulting data is going to be published on the [C3S Climate Data Store](https://cds.climate.copernicus.eu). The following heat stress indicators are included:
<ol>
<li>Heat Index</li>
<li>Humidex</li>
<li>Universal Thermal Climate Index</li>
<li>Wet-Bulb Temperature</li>
<li>Wet-Bulb Globe Temperature</li>
</ol>

More details about the heat stress indicators can be found in [Schwingshackl et al. (2021)](https://doi.org/10.1029/2020EF001885).

### Update (10.12.2024):
The calculation of Wet-Bulb Temperature (WBT) in Schwingshackl et al. (2021) has been updated after the publication of the paper. The update calculates WBT based on the iterative solution from [Buzan et al. (2015)](https://doi.org/10.5194/gmd-8-151-2015). The new WBT formula can be found on [https://github.com/jrbuzan/WetBulb.m](https://github.com/jrbuzan/WetBulb.m).
The text in Schwingshackl et al. (2021) and all figures have been updated accordingly on 18 October 2022.
