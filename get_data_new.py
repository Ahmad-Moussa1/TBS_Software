import copernicusmarine
import parameters as p
import os
import threading
import sys
import time

# when you run the script you will be prompted to input your user and password of copernicusmarine services


"""
Full name : Global Ocean Waves Analysis and Forecast
Product ID : GLOBAL_ANALYSISFORECAST_WAV_001_027
dataset_id : cmems_mod_glo_wav_anfc_0.083deg_PT3H-i
Spatial extent : Global OceanLat -80° to 90°Lon -180° to 179.92°
Spatial resolution : 0.083° × 0.083°
Temporal resolution : 3 Hours

this is the joker zone, please make sure to choose this zone at boundaries between zones
"""

"""
Full name : Mediterranean Sea Waves Analysis and Forecast
Product ID : MEDSEA_ANALYSISFORECAST_WAV_006_017
dataset_id : cmems_mod_med_wav_anfc_4.2km_PT1H-i
Spatial extent : Mediterranean Sea Lat 30.19° to 45.98° Lon -18.12° to 36.29°
Spatial resolution : 0.042° × 0.042°
Temporal resolution : Hourly
"""

"""
Full name : Atlantic-Iberian Biscay Irish- Ocean Wave Analysis and Forecast
Product ID : IBI_ANALYSISFORECAST_WAV_005_005
dataset_id : cmems_mod_ibi_wav_anfc_0.05deg_PT1H-i
Spatial extent : Atlantic-Iberia-Biscay-IrelandAtlantic: NorthLat 26° to 56°Lon -19° to 5°
Spatial resolution : 0.05° × 0.05°
Temporal resolution : Hourly
"""

"""
Full name : Black Sea Significant Wave Height extreme from Reanalysis
Product ID: BLKSEA_OMI_SEASTATE_extreme_var_swh_mean_and_anomaly
dataset_id : not_ready_yet_from_copernicusmarine, alternatively use global or download .nc file manually
dataset_id : cmems_mod_glo_wav_anfc_0.083deg_PT3H-i
Spatial extent : Black SeaLat 40.86° to 46.8°Lon 27.37° to 41.96°
Spatial resolution : 0.028° × 0.037°
Temporal resolution: Hourly
"""

"""
Full name : Baltic Sea Wave Analysis and Forecast
Product ID : BALTICSEA_ANALYSISFORECAST_WAV_003_010
dataset_id : cmems_mod_bal_wav_anfc_PT1H-i
Spatial extent : Baltic SeaLat 53.01° to 65.91°Lon 9.01° to 30.21°
Spatial resolution : 2 × 2 km
Temporal resolution : Hourly
"""

"""
Full name : Arctic Ocean Wave Analysis and Forecast
Product ID : ARCTIC_ANALYSIS_FORECAST_WAV_002_014
dataset_id : not_ready_yet_from_copernicusmarine, alternatively use global or download .nc file manually
dataset_id : cmems_mod_glo_wav_anfc_0.083deg_PT3H-i
Spatial extent : Arctic OceanLat 63° to 90°Lon -180° to 180°
Spatial resolution : 3 × 3 km
Temporal resolution : Hourly
"""

class Spinner:
    def __init__(self, delay=0.1):
        self.spinner_cycle = ['|', '/', '-', '\\']
        self.delay = delay
        self.running = False
        self.thread = threading.Thread(target=self._spin)

    def start(self):
        self.running = True
        self.thread.start()

    def _spin(self):
        while self.running:
            for char in self.spinner_cycle:
                if not self.running:
                    break
                sys.stdout.write('\r' + char)
                sys.stdout.flush()
                time.sleep(self.delay)

    def stop(self):
        self.running = False
        self.thread.join()
        sys.stdout.write('\r')
        sys.stdout.flush()


directory_path = p.dir_waveDATA
# Initialize an empty list to store file paths
file_paths = []


"""

for i in range(0,p.Zones_num,1):
    wave_prod = p.zones[i]
    if wave_prod == 1:
        ### GLOBAL OCEAN
        dataset = "cmems_mod_glo_wav_anfc_0.083deg_PT3H-i"
        prod = "GLOBAL"

    elif wave_prod == 2:
        ### Service MEDSEA
        dataset = "cmems_mod_med_wav_anfc_4.2km_PT1H-i"
        prod = "MEDSEA"

    elif wave_prod == 3:
        ### Service IBI
        dataset = "cmems_mod_ibi_wav_anfc_0.05deg_PT1H-i"
        prod = "IBI"

    elif wave_prod == 4:
        ### Service Blacksea
        dataset = "cmems_mod_glo_wav_anfc_0.083deg_PT3H-i"
        prod = "BLCKSEA"

    elif wave_prod == 5:
        ### Service Balticsea
        dataset = "cmems_mod_bal_wav_anfc_PT1H-i"
        prod = "BALTIC"

    elif wave_prod == 6:
        ### Arctic Sea
        dataset = "cmems_mod_glo_wav_anfc_0.083deg_PT3H-i"
        prod = "ARCTIC"

    copernicusmarine.subset(
        dataset_id= dataset,
        variables=["VHM0", "VMDR","VTPK"],
        minimum_longitude=p.lonMin[i]-p.dx,
        maximum_longitude=p.lonMax[i] + p.dx,
        minimum_latitude=p.LatitMin[i] - p.dx,
        maximum_latitude=p.LatitMax[i] + p.dx,
        start_datetime=p.ETD,
        end_datetime=p.ETA,
        minimum_depth=0,
        maximum_depth=1,
        output_filename = 'Waves_'+str(i)+'.nc',
        output_directory = directory_path)
"""


# Loop through the directory and append file paths to the list
for filename in os.listdir(directory_path):
    if os.path.isfile(os.path.join(directory_path, filename)):
        file_paths.append(os.path.join(directory_path, filename))

file_paths_comp = file_paths
"""
# check the following code snippet
# Loop over the list of file paths
for file_path in file_paths:
    print(f"Processing file: {file_path}")
"""
"""
if __name__ == "__main__":
    spinner = Spinner()

    print("Start downloading.......")
    spinner.start()

    main()

    spinner.stop()
    print("Downloading completed!")

"""