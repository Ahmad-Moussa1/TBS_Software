This software conists mainly of 5 scripts developed in python to optimize sailing time by choosing optimimum waypoints based on weather
conditions. it is a dynamic optimization solution proposed through TBS algorithm and TBS software implementation.

programming langauge: Mainly python 

accompanying libraries: NumPy, Pandas, Matplotlib, xarray, and Cartopy 
all these libraries must be installed on your system with python.

All the commands will be executed in command prompet with the name of the script.py using python keyword before it

***************************************************************************
The software structure consists of 5 Python scripts as follows: 
•	Parameters.py
•	Get_data_new.py
•	Methods.py
•	Main1.py
•	Main2.py 
The input variables are consolidated in a single file, parameters.py, which contains all the necessary input parameters
for the simulation. These parameters are logged by the user. The code execution is structured into three sequential steps: 
downloading data, Phase 1, and Phase 2.
The script get_data_new.py downloads sea surface wave variables, specifically significant wave height (meters), 
wave direction (degrees measured clockwise from north), and wave period (seconds). 
This data is organized in netCDF format (.nc) files based on time, latitudes, and longitudes and is used to calculate 
involuntary speed reduction using various mathematical models and optimization techniques. 
The methods.py script includes all methods and functions used to complete the simulation.
Phase 1, executed by main1.py, validates the input parameters and outputs a map and report on these constraints. 
Phase 2, executed by main2.py, runs the TBS algorithm using optimization techniques to select the optimum waypoints for the voyage,
resulting in a voyage map and two reports: simulation metadata and a simulation output report that includes the chosen waypoints, 
total sailing time, and total distance covered.

Main_initial_trials directory contains outputs of initial trials of TBS software for phase 1 and phase 2
Validation directory contains comparisons between TBS software and SIMROUTE software to validate results.

*************************************************************
execution sequecnce:

1- input the parameters used for simulation in parameters.py file
2- run get_data_new.py script to download data
3- run main1.py to check and validate your inputs
4- check output report and if all checks are validated move to phase 2
5- run main2.py to get the main output maps and text files demonstrating chosen waypoints and total sailing time and distance

**************************************************************

sample input (scenario metadata): (Alex, Egypt) to (Almeria, Spain)

Simulation name	TBS_Alex_Almeria
Paths number 	5
Start point coord.	(31.25, 29.8800659)
End point coord.	(36.8299675, -2.47345209)
Paths coordinates	[(31.25, 29.8800659), (35, 15), (37, 12), (38, 10), (36.4, 0), (36.8299675, -2.47345209)]
ETD	2024-03-16T14:00:00
Wave data zone	Mediterranean Sea
Time resolution	Hourly 
Distance per head.	9
Reduction factor	0.3
Gain factor	1.1
Degree steps	7
OBC diameter	63
Math model used	1
Service speed	16
Max. speed	20
TBS radius	20
LBP	250
DWT	56000



*****************************************************************

TBS check report (output of phase 1) for TBS simulation of (Alex, Almeria) case

             PHASE I report:     TBS_Alex_Almeria
=========================================================================
1st Check: OBC Diameter is greater than TBS radius && smaller than distance between paths coordinates
** number of checks equal number of paths 
** All values must be "OBC Diameter is acceptable"

path 1: OBC Diameter is acceptable
path 2: OBC Diameter is acceptable
path 3: OBC Diameter is acceptable
path 4: OBC Diameter is acceptable
path 5: OBC Diameter is acceptable

 2nd Check: Check OBC Diameter (specially Last OBC in each path) is smaller than TBS radius
** number of checks equal number of OBCs
** All Values must be "Last OBC Diameter is acceptable"

OBC number 1:  OBC diameter 2nd Check is acceptable
OBC number 2:  OBC diameter 2nd Check is acceptable
OBC number 3:  OBC diameter 2nd Check is acceptable
OBC number 4:  OBC diameter 2nd Check is acceptable
OBC number 5:  OBC diameter 2nd Check is acceptable
OBC number 6:  OBC diameter 2nd Check is acceptable
OBC number 7:  OBC diameter 2nd Check is acceptable
OBC number 8:  OBC diameter 2nd Check is acceptable
OBC number 9:  OBC diameter 2nd Check is acceptable
OBC number 10:  OBC diameter 2nd Check is acceptable
OBC number 11:  OBC diameter 2nd Check is acceptable
OBC number 12:  OBC diameter 2nd Check is acceptable
OBC number 13:  OBC diameter 2nd Check is acceptable
OBC number 14:  OBC diameter 2nd Check is acceptable
OBC number 15:  OBC diameter 2nd Check is acceptable
OBC number 16:  OBC diameter 2nd Check is acceptable
OBC number 17:  OBC diameter 2nd Check is acceptable
OBC number 18:  OBC diameter 2nd Check is acceptable
OBC number 19:  OBC diameter 2nd Check is acceptable
OBC number 20:  OBC diameter 2nd Check is acceptable
OBC number 21:  OBC diameter 2nd Check is acceptable
OBC number 22:  OBC diameter 2nd Check is acceptable
OBC number 23:  OBC diameter 2nd Check is acceptable
OBC number 24:  OBC diameter 2nd Check is acceptable
OBC number 25:  OBC diameter 2nd Check is acceptable
OBC number 26:  OBC diameter 2nd Check is acceptable
OBC number 27:  OBC diameter 2nd Check is acceptable
OBC number 28:  OBC diameter 2nd Check is acceptable

 3rd Check: OBC coordinates status (is_land OR is_water) 
** number of checks equal number of OBCs +1 
** All values of status must be is_water

Latitude      Longitude       Status    
  31.250       29.880         ['is_water'] 
  31.616       28.728         ['is_water'] 
  31.983       27.571         ['is_water'] 
  32.349       26.409         ['is_water'] 
  32.715       25.243         ['is_water'] 
  33.081       24.072         ['is_water'] 
  33.448       22.896         ['is_water'] 
  33.814       21.715         ['is_water'] 
  34.179       20.528         ['is_water'] 
  34.545       19.337         ['is_water'] 
  34.911       18.140         ['is_water'] 
  35.277       16.938         ['is_water'] 
  35.642       15.731         ['is_water'] 
  35.000       15.000         ['is_water'] 
  35.676       14.016         ['is_water'] 
  36.351       13.023         ['is_water'] 
  37.000       12.000         ['is_water'] 
  37.564       10.888         ['is_water'] 
  38.000       10.000         ['is_water'] 
  37.841        8.685         ['is_water'] 
  37.682        7.373         ['is_water'] 
  37.524        6.064         ['is_water'] 
  37.365        4.758         ['is_water'] 
  37.206        3.454         ['is_water'] 
  37.048        2.153         ['is_water'] 
  36.889        0.855         ['is_water'] 
  36.400        0.000         ['is_water'] 
  36.628       -1.274         ['is_water'] 
  36.830       -2.473         ['is_water'] 
--------------------End of Report------------





TBS simulation report (waypoints) for TBS simulation of (Alex, Almeria) case

PHASE II Simulation Report:     TBS_Alex_Almeria
=========================================================================

   (Lat1,Long1)         (Lat2,Long2)             heading            distance           time    

  31.250,  29.880       31.396,  29.875          358.302             8.760             0.547 
  31.396,  29.875       31.521,  29.792          330.396             8.634             0.540 
  31.521,  29.792       31.646,  29.667          319.606             9.859             0.616 
  31.646,  29.667       31.646,  29.542          270.033             6.389             0.399 
  31.646,  29.542       31.729,  29.417          308.110             8.113             0.507 
  31.729,  29.417       31.616,  28.728          259.299            35.848             2.241 
  31.616,  28.728       31.646,  28.542          280.559             9.674             0.605 
  31.646,  28.542       31.688,  28.375          286.413             8.877             0.555 
  31.688,  28.375       31.688,  28.208          270.044             8.515             0.532 
  31.688,  28.208       31.688,  28.042          270.044             8.515             0.532 
  31.688,  28.042       31.729,  27.917          291.429             6.857             0.429 
  31.729,  27.917       31.983,  27.571          310.899            23.303             1.456 
  31.983,  27.571       32.021,  27.417          286.259             8.173             0.511 
  32.021,  27.417       32.062,  27.375          319.726             3.280             0.205 
  32.062,  27.375       32.146,  27.333          337.055             5.434             0.340 
  32.146,  27.333       32.146,  27.250          270.022             4.236             0.265 
  32.146,  27.250       32.188,  27.167          300.592             4.919             0.307 
  32.188,  27.167       32.349,  26.409          284.364            39.659             2.479 
  32.349,  26.409       32.479,  26.375          347.487             7.999             0.500 
  32.479,  26.375       32.604,  26.292          330.686             8.609             0.538 
  32.604,  26.292       32.646,  26.250          319.905             3.271             0.204 
  32.646,  26.250       32.688,  26.208          319.920             3.270             0.204 
  32.688,  26.208       32.729,  26.125          300.743             4.897             0.306 
  32.729,  26.125       32.715,  25.243          269.170            44.565             2.785 
  32.715,  25.243       32.854,  25.208          348.191             8.517             0.532 
  32.854,  25.208       32.979,  25.167          344.379             7.793             0.487 
  32.979,  25.167       33.104,  25.083          330.824             8.598             0.537 
  33.104,  25.083       33.229,  24.958          320.102             9.787             0.612 
  33.229,  24.958       33.271,  24.875          300.896             4.875             0.305 
  33.271,  24.875       33.081,  24.072          254.490            41.938             2.621 
  33.081,  24.072       33.146,  23.917          296.405             8.703             0.544 
  33.146,  23.917       33.271,  23.792          320.116             9.786             0.612 
  33.271,  23.792       33.396,  23.667          320.156             9.780             0.611 
  33.396,  23.667       33.521,  23.542          320.197             9.774             0.611 
  33.521,  23.542       33.604,  23.375          301.012             9.724             0.608 
  33.604,  23.375       33.448,  22.896          248.728            25.770             1.611 
  33.448,  22.896       33.604,  22.875          353.734             9.461             0.591 
  33.604,  22.875       33.729,  22.833          344.506             7.789             0.487 
  33.729,  22.833       33.812,  22.792          337.442             5.418             0.339 
  33.812,  22.792       33.812,  22.708          270.023             4.157             0.260 
  33.812,  22.708       33.854,  22.625          301.069             4.851             0.303 
  33.854,  22.625       33.814,  21.715          267.177            45.472             2.842 
  33.814,  21.715       33.938,  21.625          329.072             8.680             0.543 
  33.938,  21.625       34.062,  21.500          320.375             9.749             0.609 
  34.062,  21.500       34.146,  21.458          337.521             5.415             0.338 
  34.146,  21.458       34.188,  21.417          320.408             3.247             0.203 
  34.188,  21.417       34.271,  21.250          301.210             9.669             0.604 
  34.271,  21.250       34.179,  20.528          261.493            36.242             2.265 
  34.179,  20.528       34.312,  20.417          325.274             9.725             0.608 
  34.312,  20.417       34.438,  20.375          344.629             7.784             0.486 
  34.438,  20.375       34.562,  20.292          331.238             8.563             0.535 
  34.562,  20.292       34.562,  20.125          270.047             8.241             0.515 
  34.562,  20.125       34.646,  20.000          309.041             7.949             0.497 
  34.646,  20.000       34.545,  19.337          259.747            33.320             2.083 
  34.545,  19.337       34.562,  19.250          283.555             4.425             0.277 
  34.562,  19.250       34.688,  19.208          344.673             7.782             0.486 
  34.688,  19.208       34.771,  19.167          337.673             5.409             0.338 
  34.771,  19.167       34.771,  19.083          270.024             4.110             0.257 
  34.771,  19.083       34.854,  18.958          309.113             7.937             0.496 
  34.854,  18.958       34.911,  18.140          275.076            40.433             2.527 
  34.911,  18.140       35.062,  18.125          355.263             9.127             0.570 
  35.062,  18.125       35.146,  17.958          301.479             9.594             0.600 
  35.146,  17.958       35.229,  17.917          337.786             5.405             0.338 
  35.229,  17.917       35.354,  17.792          320.814             9.688             0.605 
  35.354,  17.792       35.396,  17.625          287.093             8.534             0.533 
  35.396,  17.625       35.277,  16.938          258.191            34.387             2.149 
  35.277,  16.938       35.438,  16.917          353.747             9.714             0.607 
  35.438,  16.917       35.521,  16.750          301.598             9.562             0.598 
  35.521,  16.750       35.562,  16.708          320.877             3.225             0.202 
  35.562,  16.708       35.562,  16.625          270.024             4.070             0.254 
  35.562,  16.625       35.646,  16.458          301.638             9.551             0.597 
  35.646,  16.458       35.642,  15.731          269.867            35.499             2.219 
  35.642,  15.731       35.688,  15.542          286.459             9.616             0.601 
  35.688,  15.542       35.729,  15.375          287.162             8.502             0.531 
  35.729,  15.375       35.729,  15.208          270.049             8.123             0.508 
  35.729,  15.208       35.604,  15.083          219.128             9.670             0.604 
  35.604,  15.083       35.562,  15.042          219.131             3.225             0.202 
  35.562,  15.042       35.000,  15.000          183.473            33.834             2.115 
  35.000,  15.000       35.146,  14.917          334.960             9.666             0.604 
  35.146,  14.917       35.229,  14.750          301.507             9.587             0.599 
  35.229,  14.750       35.354,  14.750            0.000             7.505             0.469 
  35.354,  14.750       35.479,  14.625          320.858             9.682             0.605 
  35.479,  14.625       35.604,  14.500          320.901             9.676             0.605 
  35.604,  14.500       35.676,  14.016          280.431            24.012             1.501 
  35.676,  14.016       35.729,  13.833          289.927             9.460             0.591 
  35.729,  13.833       35.854,  13.792          344.882             7.775             0.486 
  35.854,  13.792       35.938,  13.625          301.732             9.526             0.595 
  35.938,  13.625       35.979,  13.583          321.025             3.219             0.201 
  35.979,  13.583       36.021,  13.542          321.038             3.218             0.201 
  36.021,  13.542       36.351,  13.023          308.432            32.010             2.001 
  36.351,  13.023       36.479,  13.000          351.737             7.769             0.486 
  36.479,  13.000       36.521,  12.917          301.905             4.736             0.296 
  36.521,  12.917       36.604,  12.750          301.953             9.468             0.592 
  36.604,  12.750       36.646,  12.583          287.351             8.412             0.526 
  36.646,  12.583       36.688,  12.542          321.280             3.207             0.200 
  36.688,  12.542       37.000,  12.000          305.950            32.084             2.005 
  37.000,  12.000       37.021,  11.833          278.946             8.088             0.505 
  37.021,  11.833       37.021,  11.667          270.050             7.990             0.499 
  37.021,  11.667       37.021,  11.583          270.025             3.995             0.250 
  37.021,  11.583       37.062,  11.417          287.442             8.370             0.523 
  37.062,  11.417       37.104,  11.250          287.451             8.366             0.523 
  37.104,  11.250       37.564,  10.888          328.019            32.562             2.035 
  37.564,  10.888       37.521,  10.750          248.622             7.039             0.440 
  37.521,  10.750       37.604,  10.750            0.000             5.004             0.313 
  37.604,  10.750       37.646,  10.583          287.568             8.311             0.519 
  37.646,  10.583       37.771,  10.500          332.217             8.485             0.530 
  37.771,  10.500       37.812,  10.417          302.349             4.679             0.292 
  37.812,  10.417       38.000,  10.000          299.825            22.723             1.420 
  38.000,  10.000       37.979,   9.833          261.039             7.985             0.499 
  37.979,   9.833       38.021,   9.792          321.773             3.185             0.199 
  38.021,   9.792       37.979,   9.625          252.450             8.273             0.517 
  37.979,   9.625       37.938,   9.542          237.644             4.671             0.292 
  37.938,   9.542       37.896,   9.417          247.132             6.428             0.402 
  37.896,   9.417       37.841,   8.685          264.818            34.825             2.177 
  37.841,   8.685       37.979,   8.583          329.818             9.587             0.599 
  37.979,   8.583       38.062,   8.500          321.795             6.369             0.398 
  38.062,   8.500       38.062,   8.292          270.064             9.848             0.616 
  38.062,   8.292       37.979,   8.167          229.799             7.745             0.484 
  37.979,   8.167       37.938,   8.042          247.121             6.425             0.402 
  37.938,   8.042       37.682,   7.373          244.424            35.215             2.201 
  37.682,   7.373       37.688,   7.250          273.033             5.861             0.366 
  37.688,   7.250       37.729,   7.083          287.588             8.303             0.519 
  37.729,   7.083       37.688,   6.958          247.190             6.443             0.403 
  37.688,   6.958       37.688,   6.875          270.025             3.959             0.247 
  37.688,   6.875       37.688,   6.792          270.025             3.959             0.247 
  37.688,   6.792       37.524,   6.064          254.357            35.984             2.249 
  37.524,   6.064       37.562,   5.875          284.589             9.294             0.581 
  37.562,   5.875       37.562,   5.792          270.025             3.966             0.248 
  37.562,   5.792       37.604,   5.708          302.276             4.688             0.293 
  37.604,   5.708       37.562,   5.625          237.775             4.688             0.293 
  37.562,   5.625       37.562,   5.542          270.025             3.966             0.248 
  37.562,   5.542       37.365,   4.758          252.630            39.205             2.450 
  37.365,   4.758       37.312,   4.750          186.530             3.170             0.198 
  37.312,   4.750       37.188,   4.667          207.979             8.496             0.531 
  37.188,   4.667       37.271,   4.500          302.177             9.408             0.588 
  37.271,   4.500       37.312,   4.333          287.496             8.345             0.522 
  37.312,   4.333       37.312,   4.250          270.025             3.979             0.249 
  37.312,   4.250       37.206,   3.454          260.727            38.574             2.411 
  37.206,   3.454       37.062,   3.375          203.642             9.424             0.589 
  37.062,   3.375       36.938,   3.250          218.650             9.605             0.600 
  36.938,   3.250       36.979,   3.167          302.061             4.716             0.295 
  36.979,   3.167       36.896,   3.042          230.206             7.812             0.488 
  36.896,   3.042       36.854,   2.917          247.418             6.504             0.406 
  36.854,   2.917       37.048,   2.153          287.822            38.443             2.403 
  37.048,   2.153       37.146,   2.000          308.865             9.400             0.587 
  37.146,   2.000       37.188,   1.833          287.469             8.357             0.522 
  37.188,   1.833       37.062,   1.708          218.604             9.599             0.600 
  37.062,   1.708       37.104,   1.542          287.451             8.366             0.523 
  37.104,   1.542       37.062,   1.375          252.650             8.366             0.523 
  37.062,   1.375       36.889,   0.855          247.515            27.041             1.690 
  36.889,   0.855       36.812,   0.708          236.851             8.403             0.525 
  36.812,   0.708       36.688,   0.667          194.966             7.768             0.486 
  36.688,   0.667       36.646,   0.542          247.471             6.519             0.407 
  36.646,   0.542       36.521,   0.417          218.801             9.625             0.602 
  36.521,   0.417       36.396,   0.292          218.847             9.631             0.602 
  36.396,   0.292       36.400,   0.000          271.104            14.098             0.881 
  36.400,   0.000       36.438,  -0.042          318.212             3.020             0.189 
  36.438,  -0.042       36.479,  -0.208          287.317             8.428             0.527 
  36.479,  -0.208       36.479,  -0.292          270.025             4.023             0.251 
  36.479,  -0.292       36.562,  -0.458          301.938             9.471             0.592 
  36.562,  -0.458       36.604,  -0.542          301.935             4.733             0.296 
  36.604,  -0.542       36.628,  -1.274          272.587            35.334             2.208 
  36.628,  -1.274       36.646,  -1.333          290.109             3.031             0.189 
  36.646,  -1.333       36.771,  -1.417          331.902             8.510             0.532 
  36.771,  -1.417       36.854,  -1.583          302.037             9.446             0.590 
  36.854,  -1.583       36.938,  -1.750          302.064             9.438             0.590 
  36.938,  -1.750       36.979,  -1.875          292.682             6.498             0.406 
  36.979,  -1.875       36.830,  -2.473          252.863            30.096             1.881 

Number of Way points: 169
Total Voyage Distance: 1959.027680742692
Total Sailing Time: 122.43923004641825
--------------------End of Report----------------------


 

