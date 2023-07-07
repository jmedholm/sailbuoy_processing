# Sailbuoy SO-CHIC
Sailbuoy processing for the SO-CHIC mission
[![DOI](https://zenodo.org/badge/534208314.svg)](https://zenodo.org/badge/latestdoi/534208314)

These sets of notebooks deals with the ADCP sensor on the Sailbuoy, and process the data, and combines it with the output from the other sensors.

Different outputs:
1. Datalogger - Compiled output from all sensors
2. Autopilot - Navigational variables
3. DCPS - Output from the DCPS-sensor
4. Airmar - Output from the Airmar, incl. relative humidity

# Sailbuoy Kringla DCPS

Kringla was deployed in the Southern Ocean for the SO-CHIC project, analyzing surface characteristics in a couple of subregions. Settings for the DCPS were:
- **Basic output + beams**
- **Cell size** = 2 m
- **Cell spacing** = 2 m
- **Number of cells** = 40
- **Number of pings** = 150

We chose the above settings to have a good compromise between range, accuracy, and power consumption. The extra output allows us to do a rigorous QC  

And for the datalogger:
- **$LRMIN 15**

Minutes between samples. This was later (8 June) changed to 10 minutes.
- **$DCPS.MIN 3**  

Minutes that the sensor was on.
- **$DCPS.AVG 6**  

Number of cells to average and send over Iridium.
- **$DCPS.DEV -21**  

Local magnetic declination, this was added after 16 Feb, and then continuously updated based on latitude. Does not apply to the DCPS.txt.
- **$STRDL 300**  

Startup delay, 5 minutes. Was added 22 Feb to let the SB complete each tack, which could take 2 minutes. This was later (25 March) changed to 3 minutes. This does affect DCPS.txt as without it the SB -could- have had a tack while measuring, which affects the calculated speed and heading.
- **$DCPS.RC 4**  

Run counter. I.e., let the sensor open every 4th time the datalogger turns on. This was done 24th May to conserve battery. Does affect DCPS.txt since it controls how often the sensor samples. This was changed later (8 June) to 6, to comply with an increase in frequency to improve T/S sampling. It was also changed 3 July to 3, to increase sampling frequency to 30 minutes, and changed to 1 (10 minutes) on 8 July.


