"""
Parameters file in which planner will input main data for Algorithm execution,
All values are default and the planner will modify these values regarding the intended voyage

"""
import numpy as np


#Simulation name
name_Simu='TBS_Muscat_JebelALI'

# number of Paths: (int)
paths_num = 2

# coordinates (latits,longs) of start, intermediate and end points over the paths. paths_coord elements = paths_num + 1.
paths_coord = [(23.75,58.5),(26.55,56.5),(25.0236,54.9)]
#(26.55,56.5),(26.25,56)
#Data to donwnload wave files from CMEMS server:
ETD = '2024-05-08T06:00:00'        #[year-month-dayT00:00:00]   # Estimated Time Deperature ETD, Always make seconds 00
ETA = '2024-05-10T23:59:59'        # Estimated Time Arrival ETA, it will be maximum value to ensure all data downloaded contains the time zones of required data




# folder where the data downloaded will be saved and you can make subfolders for more specification, so make it TBS_Trial_2 and _3 and etc
dir_waveDATA='D:\Study\Masters of Marine Engineering\Thesis\weather routing\Journal paper\TBS_software\Validation\case_5\TBS\storeWaves'

dir_out_phaseI ='D:\Study\Masters of Marine Engineering\Thesis\weather routing\Journal paper\TBS_software\Validation\case_5\TBS\output_phase_1'

dir_out_phaseII='D:\Study\Masters of Marine Engineering\Thesis\weather routing\Journal paper\TBS_software\Validation\case_5\TBS\output_phase_2'


# number of Zones: (int) they are the weather zones where the voyage will pass through, based on paths coordinates start, end and intermediate points
Zones_num= 1

# Zones Wave models number for each given zone (it will be the same as zones number and the values will be of the CMEMS product you need)
zones = [1]

#Time resolution of CMEMS product for each Zone, also it is the TCV limit for TBS based on the zones:
time_res = [3]             # In hours: 1 for all regions except the GLOBAL it will be 3 hours

# for each zone, these data must be obtained: (lonMin, lonMax, LatitMin, LatitMax) to obtain data in this range from the server,
# you have to ensure that your paths coordinates bounded by these values not outside it.
# these values are obtained from the coordinates of the paths and the zones of waves
# in the following arrays, the elements of each array will be = no of zones.
# for maximum latitudes, add OBC diameter to maximum and minimum latitudes in path coordinates
# each degree in latitudes approximatelly equals 60 nautical miles and for longitudes maximum will be 60 at equator and 0 at poles
# Adjust the following parameters till you got status for each point of OBC (is_land/is_water) to go to next phase

lonMin= [50]
lonMax= [60]
LatitMin= [20]
LatitMax= [30]
# ensure you have extra boundaries when you get the data to ensure all coordinates are bounded inside
dx= 0.5


# distance/heading (nautical miles / waypoint), it is an approximate value because real chosen way point may be any value less than or equal it: int
# distance/heading will be appropraite to be used in all the zones taking into consideration their data resolution, it is better to choose dist/head based
# on the lower resolution.
# this is a maximum value where way points distance may be this or less
# distance/ heading cannot be more than the distance covered by the ship with its service speed if the resolution is one hour,
# the way points distance will be between min: 70% of dist/head and max dist/head, and this will be a loss of 30% of total distance when
# divided over all TBCS it will be smaller than 1 mile
dist_per_head = 9
value = np.power(dist_per_head,2) + np.power(dist_per_head,2)
dist_per_head_diagonal = np.sqrt(value)
reduction_factor = 0.3 # this value indicates the min distance of way point chosen
min_dist_per_head = reduction_factor * dist_per_head
gain_factor = 1.4  # this value indicates max distance of way points chosen
max_dist_per_head = gain_factor*dist_per_head

# Degree steps, for outer boundary circles, it will be multiples or halves or thirds or quarters of degree of 90.(30, 45, 90, 180, ....) (step)
deg_step = 9

# diameter of outer boundary circles(OBC) = ditance / heading * degree steps in nautical miles:
OBC_diameter = dist_per_head * deg_step

# Mathematical Model used to calculate involuntary speed reduction (1= without_khokhlov, 2 = with_Khokhlov)
# Khokhlov method has 2 main constraints to be used:
# * service speed >= 9 knots & service speed <= 20 knots
# * DWT >= 4000 and DWT <= 20000
# so, be cautious to fulfill these 2 constraints before using it
math_model = 1

# service speed of the ship in knots (it will be constant all over the voyage with variable power consumption to optimize sailing time)
serv_speed = 18

# Maximum speed of the ship in knots (for TBS construction)
Max_speed = 22

time_res_max = np.max(time_res)
TBS_radius = time_res_max * Max_speed

# length between perpendiculars of the ship in meters
LBP = 125

# Dead weight of the ship in tonnes
DWT = 38000

