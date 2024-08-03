import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import get_data_new as data
import parameters as param
import numpy as np
import xarray as xr
from collections import Counter



def calc_distance(lat1,lon1,lat2,lon2):
    """
    calculate distance between 2 points in nautical miles using spherical law of cosines
    """
# Convert latitude and longitude from degrees to radians
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)

# Calculate the angular distance between the two points using the spherical law of cosines
    angular_distance = np.arccos(np.sin(lat1) * np.sin(lat2) +
                              np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1))

# Radius of the Earth in nautical miles
    earth_radius = 3440.065

# Calculate the distance using the formula: distance = angular_distance * earth_radius
    distance = angular_distance * earth_radius

   # r1 = 0.5*distance

    return distance

def calculate_initial_bearing(lat1, lon1, lat2, lon2):
    """
    Calculate the initial bearing between two points.
    Inputs:
    lat1, lon1: Latitude and longitude of the first point (in degrees).
    lat2, lon2: Latitude and longitude of the second point (in degrees).
    Returns:
    Initial bearing in degrees (clockwise from north).
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    delta_lon = lon2 - lon1

    y = np.sin(delta_lon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon)

    initial_bearing = np.arctan2(y, x)

    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def destination_point(lat1, lon1, distance, initial_bearing):
    """
    Calculate the destination point given a start point, distance, and bearing.

    Parameters:
    lat1 (float): Latitude of the starting point in degrees.
    lon1 (float): Longitude of the starting point in degrees.
    distance (float): Distance to the destination point in Nautical miles.
    bearing (float): Bearing (azimuth) in degrees from the starting point.

    Returns:
    (float, float): Latitude and longitude of the destination point in degrees.
    """
    # Convert latitude and longitude from degrees to radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    bearing_1 = np.radians(initial_bearing)

    # Earth's radius in Nautical Miles
    R = 3440.065

    # Calculate the destination point's latitude
    lat2 = np.arcsin(np.sin(lat1) * np.cos(distance / R) +
                     np.cos(lat1) * np.sin(distance / R) * np.cos(bearing_1))

    # Calculate the destination point's longitude
    lon2 = lon1 + np.arctan2(np.sin(bearing_1) * np.sin(distance / R) * np.cos(lat1),
                             np.cos(distance / R) - np.sin(lat1) * np.sin(lat2))

    # Convert latitude and longitude back from radians to degrees
    destiny_lat = np.degrees(lat2)
    destiny_lon = np.degrees(lon2)

    return destiny_lat, destiny_lon

def plot_path(outer_bound,boundaries_coordinates):
    """
    Plot OBCs on a map using Cartopy.

    Parameters:
    waypoints (list): List of tuples containing latitude and longitude of waypoints.
    """
    # Extract latitude and longitude from waypoints
    lats, lons = zip(*outer_bound)
    lats_1,lons_1 = zip(*boundaries_coordinates)
    fig = plt.figure(figsize=(10, 10))
    m6 = plt.axes(projection=ccrs.PlateCarree())
    #m6 = plt.contourf(lons,lats,transform=ccrs.PlateCarree())
    # (x0, x1, y0, y1)
    m6.set_extent([45, 65, 15, 35], ccrs.PlateCarree())
    #m6.set_global()

    #m6.coastlines()
    m6.add_feature(cfeature.LAND)
    m6.add_feature(cfeature.OCEAN)
    m6.add_feature(cfeature.COASTLINE)
    m6.gridlines(draw_labels=True)
    m6.add_feature(cfeature.BORDERS, linestyle=':')
    m6.scatter(lons, lats, linestyle='-', color='black', transform=ccrs.PlateCarree())
    m6.scatter(lons_1, lats_1, color='red', label='start & End', marker='*', s=30,
               linewidths=1)
    m6.legend()
    plt.title('outer boundaries circles')
    plt.show()

def plot_path2(coordinates):
    """
    plot paths on a map using cartopy
    :return:
    """
    # Separate latitudes and longitudes


    lats_1,lons_1 = zip(*coordinates)

    # Create a new plot
    fig = plt.figure(figsize=(10, 10))

    # Create a map projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([45, 65, 15, 35], ccrs.PlateCarree())

    # Add features to the map
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines(draw_labels=True)

    # Plot the points
    ax.scatter(lons_1[0], lats_1[0], color='green', label= 'start', marker='*', s=40, linewidths=2)
    ax.scatter(lons_1[-1], lats_1[-1], color='yellow', label='End', marker='*', s=40, linewidths=2)
    # Plot the lines connecting the points
    ax.plot(lons_1, lats_1, color='red', linestyle='-', transform=ccrs.PlateCarree())


    # Show the plot
    ax.legend()
    plt.title('MAP')
    plt.show()


def retrieve_data(target_lat,target_lon,target_time):
    for i in range(0, param.Zones_num, 1):
        if target_lat >= param.LatitMin[i] - param.dx and target_lat <= param.LatitMax[i] + param.dx and target_lon >= param.lonMin[i] - param.dx and target_lon <= param.lonMax[i] + param.dx:
            dm = xr.open_dataset(data.file_paths[0])
            target_lat_calc = target_lat
            target_lon_calc = target_lon
            target_time_calc = target_time
            selected_data = dm.sel(latitude=target_lat_calc, longitude=target_lon_calc, time=target_time_calc,method='pad')
            VHM0_at_point = selected_data['VHM0'].values
            VMDR_at_point = selected_data['VMDR'].values

            print(
                f"Hs at latitude {target_lat_calc} and longitude {target_lon_calc} and time {target_time_calc}: {VHM0_at_point}")
            print(
                f"theta at latitude {target_lat_calc} and longitude {target_lon_calc} and time {target_time_calc}: {VMDR_at_point}")
            break
        else:
            print(f"not found in : {i} file")
            if i == param.Zones_num - 1:
                print("Out of Range, Given Positions are out of Zones")

def check_land(target_lat,target_lon,target_time):
    for i in range(0, param.Zones_num, 1):
        status = []
        if target_lat >= param.LatitMin[i] - param.dx and target_lat <= param.LatitMax[i] + param.dx and target_lon >= param.lonMin[i] - param.dx and target_lon <= param.lonMax[i] + param.dx:
            dm = xr.open_dataset(data.file_paths_comp[0])
            target_lat_calc = target_lat
            target_lon_calc = target_lon
            target_time_calc = target_time
            selected_data = dm.sel(latitude=target_lat_calc, longitude=target_lon_calc, time=target_time_calc,method='nearest')
            VHM0_at_point = selected_data['VHM0'].values
            VMDR_at_point = selected_data['VMDR'].values
            if np.isnan(VHM0_at_point) or np.isnan(VMDR_at_point):
                print("is_land")
                status.append("is_land")
                #print(f"check land at latitude {target_lat_calc} and longitude {target_lon_calc}: is Land")
                break
            else:
                print("is_water")
                status.append("is_water")
                #print(f"check land at latitude {target_lat_calc} and longitude {target_lon_calc}: is Water")
                break

        else:
            print(f"not found in : {i} file")
            if i == param.Zones_num - 1:
                print("Out of Range, Given Position out of Zones")


    return status

"""
transform matrix method is used to change the heading between start and end of each OBC into:
1- Divergence array
2- Convergence array
3- TBS array 
these arrays are used in calculating and excluding possible candidate points to be chosen as waypoints.
"""
def transform_matrix(outer_boundary):

    heading_array = []
    TBS_array = []
    div_array = []
    conv_array= []


    for m in range(0, len(outer_boundary) - 1, 1):
        lat3 = outer_boundary[m][0]
        lon3 = outer_boundary[m][1]
        lat4 = outer_boundary[m + 1][0]
        lon4 = outer_boundary[m + 1][1]
        heading = calculate_initial_bearing(lat3, lon3, lat4, lon4)
        heading_array.append(heading)
        #print(heading)

        if heading == 0:
            div_val_min = 270  # the same used for TBS
            div_val_max = 90
            Div_deg_arr = []
            conv_val_max = 225
            conv_val_min = 135
            Conv_deg_arr = []
            TBS = [270, 0, 90]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading == 90:
            div_val_min = 0
            div_val_max = 180
            Div_deg_arr = []
            conv_val_max = 315
            conv_val_min = 225
            Conv_deg_arr = []
            TBS = [0, 90, 180]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading == 180:
            div_val_min = 90
            div_val_max = 270
            Div_deg_arr = []
            conv_val_max = 45
            conv_val_min = 315
            Conv_deg_arr = []
            TBS = [90, 180, 270]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading == 270:
            div_val_min = 180
            div_val_max = 360
            Div_deg_arr = []
            conv_val_max = 135
            conv_val_min = 45
            Conv_deg_arr = []
            TBS = [180, 270, 360]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading > 0 and heading < 90:
            div_val_min = heading -90
            div_val_max = heading +90
            Div_deg_arr = []
            conv_val_max = 225 + heading
            conv_val_min = 135 + heading
            Conv_deg_arr = []
            TBS = [heading -90, heading-45, heading+45, heading +90]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading > 90 and heading < 180:
            div_val_min = heading - 90
            div_val_max = heading + 90
            Div_deg_arr = []
            conv_val_max = heading + 225
            conv_val_min = heading + 135
            Conv_deg_arr = []
            TBS = [heading - 90, heading-45, heading+45, heading + 90]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading > 180 and heading < 270:
            div_val_min = heading - 90
            div_val_max = heading + 90
            Div_deg_arr = []
            conv_val_max = heading - 135
            conv_val_min = heading - 225
            Conv_deg_arr = []
            TBS = [heading - 90, heading-45, heading+45, heading + 90]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1


        elif heading > 270 and heading < 360:
            div_val_min = heading - 90
            div_val_max = heading + 90
            Div_deg_arr = []
            conv_val_max = heading - 135
            conv_val_min = heading - 225
            Conv_deg_arr = []
            TBS = [heading - 90, heading-45, heading+45, heading + 90]
            delta = 90 / param.deg_step
            for i in range(0, int(param.deg_step / 2), 1):
                div_val_min = div_val_min + delta
                Div_deg_arr.append(div_val_min)
                div_val_max = div_val_max - delta
                Div_deg_arr.append(div_val_max)
                conv_val_max = conv_val_max - delta
                Conv_deg_arr.append(conv_val_max)
                conv_val_min = conv_val_min + delta
                Conv_deg_arr.append(conv_val_min)
                i += 1

        TBS_array.append(TBS)
        div_array.append(Div_deg_arr)
        conv_array.append(Conv_deg_arr)
        m += 1

    return heading_array,TBS_array,div_array,conv_array






"""
the following methods are intended to employ different math models and their combinations to calculate different involuntary speed reductions, at
different coordinates
"""

def bowditch(cleaned_df, v0):
# argument dataframe that contain cleaned data ready for implementation and V0
    cleaned_df.loc[:, 'bowditch_red'] = 0

    for index, row in cleaned_df.iterrows():
        # Access value of specific column
        column_value = row['adj_theta']
        # Perform operations based on conditions
        if ((column_value >= 0 and column_value <= 45) or (column_value >= 270 and column_value <= 360)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0083
        elif ((column_value > 45) and (column_value < 135)) or ((column_value > 225) and (column_value < 270)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0165
        elif (column_value >= 135) and (column_value <= 225):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0248

    cleaned_df['Bowditch_model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Bowditch_model'] = v0 - row['bowditch_red'] * row['VHM0'] * row['VHM0'] * 3.2808 * 3.2808



def artessen(cleaned_df,v0,LBP):
    # argument dataframe that contain cleaned data ready for implementation, v0 and LBP
    cleaned_df['Artessen_model'] = 0

    for index, row in cleaned_df.iterrows():
        column_value_1 = row['adj_theta']
        column_value_2 = row['VHM0']

        if (0 <= column_value_2 and column_value_2 < 2.5):
            m = 0
            n = 0
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)
        elif (2.5 <= column_value_2 and column_value_2 < 4.0):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 100
                n = 0
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 350
                n = 1
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 700
                n = 2
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 900
                n = 2
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (4.0 <= column_value_2 and column_value_2 < 5.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 200
                n = 1
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 500
                n = 3
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1000
                n = 5
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 1300
                n = 6
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (5.5 <= column_value_2 and column_value_2 < 7.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 400
                n = 2
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 700
                n = 5
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1400
                n = 8
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 2100
                n = 11
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (7.5 <= column_value_2):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 700
                n = 3
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 1000
                n = 7
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 2300
                n = 12
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 3600
                n = 18

            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)



def khoklov(cleaned_df,v0,DWT):
    # argument dataframe that contain cleaned data ready for implementation, v0 and DWT

    cleaned_df['Khokhlov_Model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Khokhlov_Model'] = v0 - (
                    (0.745 * row['VHM0'] - 0.245 * row['VHM0'] * row['adj_theta_rad']) * (1.0 - 1.35e-6 * DWT * v0))


def bow_khokh(cleaned_df, v0, DWT):
    cleaned_df.loc[:, 'bowditch_red'] = 0

    for index, row in cleaned_df.iterrows():
        # Access value of specific column
        column_value = row['adj_theta']
        # Perform operations based on conditions
        if ((column_value >= 0 and column_value <= 45) or (column_value >= 270 and column_value <= 360)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0083
        elif ((column_value > 45) and (column_value < 135)) or ((column_value > 225) and (column_value < 270)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0165
        elif (column_value >= 135) and (column_value <= 225):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0248

    cleaned_df['Bowditch_model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Bowditch_model'] = v0 - row['bowditch_red'] * row['VHM0'] * row['VHM0'] * 3.2808 * 3.2808

    cleaned_df['Khokhlov_Model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Khokhlov_Model'] = v0 - (
                (0.745 * row['VHM0'] - 0.245 * row['VHM0'] * row['adj_theta_rad']) * (1.0 - 1.35e-6 * DWT * v0))


    cleaned_df['Bowditch_Khokhlov'] = (cleaned_df['Bowditch_model'] + cleaned_df['Khokhlov_Model'])/2


def bow_khokh_artess(cleaned_df,v0,DWT,LBP):
    cleaned_df.loc[:, 'bowditch_red'] = 0

    for index, row in cleaned_df.iterrows():
        # Access value of specific column
        column_value = row['adj_theta']
        # Perform operations based on conditions
        if ((column_value >= 0 and column_value <= 45) or (column_value >= 270 and column_value <= 360)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0083
        elif ((column_value > 45) and (column_value < 135)) or ((column_value > 225) and (column_value < 270)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0165
        elif (column_value >= 135) and (column_value <= 225):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0248

    cleaned_df['Bowditch_model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Bowditch_model'] = v0 - row['bowditch_red'] * row['VHM0'] * row['VHM0'] * 3.2808 * 3.2808

    cleaned_df['Khokhlov_Model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Khokhlov_Model'] = v0 - (
                (0.745 * row['VHM0'] - 0.245 * row['VHM0'] * row['adj_theta_rad']) * (1.0 - 1.35e-6 * DWT * v0))

    cleaned_df['Artessen_model'] = 0

    for index, row in cleaned_df.iterrows():
        column_value_1 = row['adj_theta']
        column_value_2 = row['VHM0']

        if (0 <= column_value_2 and column_value_2 < 2.5):
            m = 0
            n = 0
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)
        elif (2.5 <= column_value_2 and column_value_2 < 4.0):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 100
                n = 0
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 350
                n = 1
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 700
                n = 2
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 900
                n = 2
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (4.0 <= column_value_2 and column_value_2 < 5.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 200
                n = 1
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 500
                n = 3
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1000
                n = 5
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 1300
                n = 6
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (5.5 <= column_value_2 and column_value_2 < 7.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 400
                n = 2
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 700
                n = 5
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1400
                n = 8
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 2100
                n = 11
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (7.5 <= column_value_2):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 700
                n = 3
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 1000
                n = 7
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 2300
                n = 12
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 3600
                n = 18

            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

    cleaned_df['Bowditch_Khokhlov_Artessen'] = (cleaned_df['Bowditch_model'] + cleaned_df['Khokhlov_Model']+ cleaned_df['Artessen_model']) / 3

    cleaned_df.drop('Artessen_model', axis=1, inplace=True)


def with_khokhlov(cleaned_df,v0,DWT,LBP):
    cleaned_df.loc[:, 'bowditch_red'] = 0

    for index, row in cleaned_df.iterrows():
        # Access value of specific column
        column_value = row['adj_theta']
        # Perform operations based on conditions
        if ((column_value >= 0 and column_value <= 45) or (column_value >= 270 and column_value <= 360)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0083
        elif ((column_value > 45) and (column_value < 135)) or ((column_value > 225) and (column_value < 270)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0165
        elif (column_value >= 135) and (column_value <= 225):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0248

    cleaned_df['Bowditch_model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Bowditch_model'] = v0 - row['bowditch_red'] * row['VHM0'] * row['VHM0'] * 3.2808 * 3.2808

    cleaned_df['Khokhlov_Model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Khokhlov_Model'] = v0 - (
                (0.745 * row['VHM0'] - 0.245 * row['VHM0'] * row['adj_theta_rad']) * (1.0 - 1.35e-6 * DWT * v0))

    cleaned_df['Artessen_model'] = 0

    for index, row in cleaned_df.iterrows():
        column_value_1 = row['adj_theta']
        column_value_2 = row['VHM0']

        if (0 <= column_value_2 and column_value_2 < 2.5):
            m = 0
            n = 0
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)
        elif (2.5 <= column_value_2 and column_value_2 < 4.0):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 100
                n = 0
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 350
                n = 1
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 700
                n = 2
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 900
                n = 2
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (4.0 <= column_value_2 and column_value_2 < 5.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 200
                n = 1
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 500
                n = 3
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1000
                n = 5
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 1300
                n = 6
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (5.5 <= column_value_2 and column_value_2 < 7.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 400
                n = 2
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 700
                n = 5
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1400
                n = 8
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 2100
                n = 11
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (7.5 <= column_value_2):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 700
                n = 3
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 1000
                n = 7
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 2300
                n = 12
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 3600
                n = 18

            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

    cleaned_df['Bowditch_Artessen'] = (cleaned_df['Bowditch_model'] + cleaned_df['Artessen_model']) / 2
    cleaned_df['Bowditch_Khokhlov'] = (cleaned_df['Bowditch_model'] + cleaned_df['Khokhlov_Model']) / 2
    cleaned_df['Artessen_Khokhlov'] = (cleaned_df['Artessen_model'] + cleaned_df['Khokhlov_Model']) / 2
    cleaned_df['Bowditch_Khokhlov_Artessen'] = (cleaned_df['Bowditch_model'] + cleaned_df['Khokhlov_Model'] + cleaned_df['Artessen_model']) / 3




def without_khokhlov(cleaned_df,v0,LBP):

    cleaned_df.loc[:, 'bowditch_red'] = 0

    for index, row in cleaned_df.iterrows():
        # Access value of specific column
        column_value = row['adj_theta']
        # Perform operations based on conditions
        if ((column_value >= 0 and column_value <= 45) or (column_value >= 270 and column_value <= 360)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0083
        elif ((column_value > 45) and (column_value < 135)) or ((column_value > 225) and (column_value < 270)):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0165
        elif (column_value >= 135) and (column_value <= 225):
            cleaned_df.loc[index, 'bowditch_red'] = 0.0248

    cleaned_df['Bowditch_model'] = 0

    for index, row in cleaned_df.iterrows():
        cleaned_df.loc[index, 'Bowditch_model'] = v0 - row['bowditch_red'] * row['VHM0'] * row['VHM0'] * 3.2808 * 3.2808

    cleaned_df['Artessen_model'] = 0

    for index, row in cleaned_df.iterrows():
        column_value_1 = row['adj_theta']
        column_value_2 = row['VHM0']

        if (0 <= column_value_2 and column_value_2 < 2.5):
            m = 0
            n = 0
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)
        elif (2.5 <= column_value_2 and column_value_2 < 4.0):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 100
                n = 0
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 350
                n = 1
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 700
                n = 2
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 900
                n = 2
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (4.0 <= column_value_2 and column_value_2 < 5.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 200
                n = 1
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 500
                n = 3
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1000
                n = 5
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 1300
                n = 6
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (5.5 <= column_value_2 and column_value_2 < 7.5):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 400
                n = 2
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 700
                n = 5
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 1400
                n = 8
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 2100
                n = 11
            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

        elif (7.5 <= column_value_2):
            if (column_value_1 >= 0 and column_value_1 <= 30) or (column_value_1 >= 330 and column_value_1 <= 360):
                m = 700
                n = 3
            if (column_value_1 > 30 and column_value_1 <= 60) or (column_value_1 >= 300 and column_value_1 < 330):
                m = 1000
                n = 7
            if (column_value_1 > 60 and column_value_1 <= 150) or (column_value_1 >= 210 and column_value_1 < 300):
                m = 2300
                n = 12
            if (column_value_1 > 150 and column_value_1 < 210):
                m = 3600
                n = 18

            cleaned_df.loc[index, 'Artessen_model'] = v0 - ((m / LBP) + n) * (v0 / 100)

    cleaned_df['Bowditch_Artessen'] = (cleaned_df['Bowditch_model'] + cleaned_df['Artessen_model']) / 2


"""
Optimization Method: this method is employed to get the optimum point between candidate waypoints based on repetitions of maximum values
upon different columns
"""

def optimize_without_khokhlov(cleaned_df):

    columns = ['Bowditch_model', 'Artessen_model', 'Bowditch_Artessen']

    max_indices = []

    # Iterate over the specified columns and get the indices of the maximum values for each column
    for col in columns:
        max_value = cleaned_df[col].max()
        indices = cleaned_df.index[cleaned_df[col] == max_value].tolist()
        max_indices.extend(indices)

    print("Indices of maximum values:", max_indices)

    # Create a Counter object to count the occurrences of each element
    counter = Counter(max_indices)

    # Find the highest frequency
    max_count = max(counter.values())

    # Find all elements with the highest frequency
    most_common_elements = [elem for elem, count in counter.items() if count == max_count]

    print("Most common elements:", most_common_elements)
    print("Count:", max_count)

    if len(most_common_elements) == 1:
        row_index = most_common_elements[0]
        # Column name
        column_name_1 = 'latitude'
        column_name_2 = 'longitude'
        column_name_3 = 'ang_ship'
        column_name_4 = 'distance'

        # Get the value using .iloc
        value_1 = cleaned_df.iloc[row_index][column_name_1]
        value_2 = cleaned_df.iloc[row_index][column_name_2]
        value_3 = cleaned_df.iloc[row_index][column_name_3]
        value_4 = cleaned_df.iloc[row_index][column_name_4]


    else:
    # here priorites will be applied priority in without khokhlov is "Bowditch" method
        row_index = cleaned_df['Bowditch_model'].idxmax()
        # Column name
        column_name_1 = 'latitude'
        column_name_2 = 'longitude'
        column_name_3 = 'ang_ship'
        column_name_4 = 'distance'

    # Get the value using .iloc
        value_1 = cleaned_df.iloc[row_index][column_name_1]
        value_2 = cleaned_df.iloc[row_index][column_name_2]
        value_3 = cleaned_df.iloc[row_index][column_name_3]
        value_4 = cleaned_df.iloc[row_index][column_name_4]

    return value_1,value_2,value_3,value_4


def optimize_with_khokhlov(cleaned_df):

    columns = ['Bowditch_model', 'Khokhlov_Model','Artessen_model', 'Bowditch_Artessen', 'Bowditch_Khokhlov', 'Artessen_Khokhlov', 'Bowditch_Khokhlov_Artessen']

    max_indices = []

    # Iterate over the specified columns and get the indices of the maximum values for each column
    for col in columns:
        max_value = cleaned_df[col].max()
        indices = cleaned_df.index[cleaned_df[col] == max_value].tolist()
        max_indices.extend(indices)

    print("Indices of maximum values:", max_indices)

    # Create a Counter object to count the occurrences of each element
    counter = Counter(max_indices)

    # Find the highest frequency
    max_count = max(counter.values())

    # Find all elements with the highest frequency
    most_common_elements = [elem for elem, count in counter.items() if count == max_count]

    print("Most common elements:", most_common_elements)
    print("Count:", max_count)

    if len(most_common_elements) == 1:
        row_index = most_common_elements[0]
        # Column name
        column_name_1 = 'latitude'
        column_name_2 = 'longitude'
        column_name_3 = 'ang_ship'
        column_name_4 = 'distance'

        # Get the value using .iloc
        value_1 = cleaned_df.iloc[row_index][column_name_1]
        value_2 = cleaned_df.iloc[row_index][column_name_2]
        value_3 = cleaned_df.iloc[row_index][column_name_3]
        value_4 = cleaned_df.iloc[row_index][column_name_4]


    else:
    # here priorites will be applied priority in with khokhlov is to choose the maximum value between three columns: "Bowditch", "Khokhlov", "Bowditch-Khokhlov"
        maxima_val=[]
        maxima_index=[]
        max_value_bow = cleaned_df['Bowditch_model'].max()
        maxima_val.append(max_value_bow)
        max_index_bow = cleaned_df['Bowditch_model'].idxmax()
        maxima_index.append(max_index_bow)
        max_value_khokh = cleaned_df['Khokhlov_Model'].max()
        maxima_val.append(max_value_khokh)
        max_index_khokh = cleaned_df['Khokhlov_Model'].idxmax()
        maxima_index.append(max_index_khokh)
        max_value_bow_khokh = cleaned_df['Bowditch_Khokhlov'].max()
        maxima_val.append(max_value_bow_khokh)
        max_index_bow_khokh = cleaned_df['Bowditch_Khokhlov'].idxmax()
        maxima_index.append(max_index_bow_khokh)

        def max_value_and_index(arr):
            # Using max() with enumerate() to get both value and index
            max_index, max_value = max(enumerate(arr), key=lambda x: x[1])
            return max_value, max_index

        max_value, max_index = max_value_and_index(maxima_val)

        row_index = maxima_index[max_index]
        # Column name
        column_name_1 = 'latitude'
        column_name_2 = 'longitude'
        column_name_3 = 'ang_ship'
        column_name_4 = 'distance'

    # Get the value using .iloc
        value_1 = cleaned_df.iloc[row_index][column_name_1]
        value_2 = cleaned_df.iloc[row_index][column_name_2]
        value_3 = cleaned_df.iloc[row_index][column_name_3]
        value_4 = cleaned_df.iloc[row_index][column_name_4]

    return value_1,value_2,value_3,value_4