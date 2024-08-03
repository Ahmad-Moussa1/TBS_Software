"""
Phase 2 will start by the main input which is "outer_boundary" array which contains Start_End coordinates of OBCs for simulation
Divide and Conquer paradigm.

main constraint is : all TBSs are perpendicular on the diameter between Start and End of OBC

Transformation phase is a must to ensure: the angles (bearings) of TBSs & Divergence - Convergence constraints to be settled
with any inclinations of heading between Start and End of OBC

Transformation matrix for divergence and convergence:
Divergence:  min --> heading - 90
             max --> heading + 90
    in the convergence: if heading > 180, subtract 180 and if heading < 180 , add 180.
Convergence: min --> heading -/+ 180 -/+ 45
All bearings are calculated from north clockwise
"""
import methods
import xarray as xr
import parameters as param
import numpy as np
import get_data_new as data
import sys
import time
import threading
from datetime import datetime, timedelta


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


def main():
    global start_time_TBS, end_time_TBS, modified_date, time_res_used, lat_bound_2, long_bound_2, long_bound_1, lat_bound_1, ds, start_time_TBS_1, end_time_TBS_1, filtered_data
    paths_dist = []  # distance of each given path, number of elements = number of paths
    paths_bearing = []
    paths_number = []
    for i in range(1, param.paths_num + 1, 1):
        paths_number.append(i)

    for i in range(0, param.paths_num, 1):
        lat1 = param.paths_coord[i][0]
        lon1 = param.paths_coord[i][1]
        lat2 = param.paths_coord[i + 1][0]
        lon2 = param.paths_coord[i + 1][1]
        distance1 = methods.calc_distance(lat1, lon1, lat2, lon2)
        bearing1 = methods.calculate_initial_bearing(lat1, lon1, lat2, lon2)
        paths_dist.append(distance1)
        paths_bearing.append(bearing1)
        i += 1

# print(paths_dist)
# print(paths_bearing)

    outer_boundary = []
    obcd_check_1 = []

    for i in range(0, param.paths_num, 1):
        n = int(paths_dist[i] // param.OBC_diameter)
        last_OBC_diameter = paths_dist[i] - (n * param.OBC_diameter)
        destiny_latitude = param.paths_coord[i][0]
        destiny_longitude = param.paths_coord[i][1]
        distance2 = param.OBC_diameter
        bearing2 = paths_bearing[i]
        destiny_tuple_initial = (destiny_latitude, destiny_longitude)
        outer_boundary.append(destiny_tuple_initial)
        # print(param.OBC_diameter)
        # print(last_OBC_diameter)
        # print(bearing2)
        # print(i)
        # print(n)
        if param.OBC_diameter > param.TBS_radius and param.OBC_diameter < paths_dist[i]:
            # print("OBC Diameter is acceptable")
            obcd_check_1.append("OBC Diameter is acceptable")

            for j in range(0, n, 1):
                destiny_latit, destiny_long = methods.destination_point(outer_boundary[-1][0], outer_boundary[-1][1],
                                                                        distance2, bearing2)
                destiny_tuple = (destiny_latit, destiny_long)
                outer_boundary.append(destiny_tuple)
                j += 1

        else:
            """print("values for degree steps and distance/heading need to be corrected, \n these values leads to OBC diameters that exceeding paths distance or smaller than TBS radius")"""
            obcd_check_1.append(
                "values for degree steps and distance/heading need to be corrected, \n these values leads to OBC diameters that exceeding paths distance or smaller than TBS radius")
    #print('1st check is done')
    #print(f'number of paths is: {len(paths_number)}')
    #print(len(obcd_check_1))
    destiny_latitude_final = param.paths_coord[-1][0]
    destiny_longitude_final = param.paths_coord[-1][1]
    destiny_tuple_final = (destiny_latitude_final, destiny_longitude_final)
    outer_boundary.append(destiny_tuple_final)

    print(f'Outer_Boundary: {outer_boundary}')

# calculate bearing between start and end of each OBC, this bearing will be transformed to fit angles of TBSs and Divergence - convergence constraints

    heading_array,TBS_array, div_array, conv_array = methods.transform_matrix(outer_boundary)
    print(f'heading array: {heading_array}')
    print(f'TBS array: {TBS_array}')
    TBS_array_1 = np.array(TBS_array)
    print("Shape of TBS_array:", TBS_array_1.shape)
    print(f'divergance_array: {div_array}')
    div_array_1 = np.array(div_array)
    print("Shape of div_array:", div_array_1.shape)
    print(f'Convergence array: {conv_array}')
    conv_array_1 = np.array(conv_array)
    print("Shape of conv_array:", conv_array_1.shape)


    way_points_out =[] # tuples of latitudes and longitudes of waypoints chosen (Main output)
    distances_out  =[] # this used to output all distances covered between waypoints (one of the outputs)
    headings_out   =[] # this used to output all headings values between each point (one of the outputs)
    TCV_array_out  =[] # this used to output all values time between each point (one of the outputs)
    TCV_array_tmp  =[] # this used for every time resolution and get zero when exceeds time resolution values
    TCV_array_tmp_1=[] # this is used for initial steps if minutes are not equal to zero
    v0 = param.serv_speed
    LBP = param.LBP
    DWT = param.DWT
    H = 0   # counter used for hours in TBS identification for start time calculations
    result = 0
    X_Z = 0 # indicator if the algorithm executes specific IF condition
    X_Y = 0 # indicator if the algorithm executes specific IF condition
    zh = 0  # counter for way_points_out






    for i in range(0, len(outer_boundary)-1, 1):
        print(f'i: {i}')
        lat10  = outer_boundary[i][0]
        long10 = outer_boundary[i][1]
        lat11  = outer_boundary[i+1][0]
        long11 = outer_boundary[i+1][1]

        # Define Bounding Box boundaries to get into the waves files and return data
        # Bounding box values will be the start point + distance/heading in 3 direction of TBS array
        # you will loop over each step

        way_points_tuples = (lat10, long10)
        way_points_out.append(way_points_tuples)

        ml = 0  # counter for div_array_1
        nl = 0  # counter for conv_array_1

        for step in range(0, int(div_array_1.shape[1]/2), 1):

            print(f'counter for way points: {zh}')
            if heading_array[i] == 0:
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] == 90:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] == 180:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] == 270:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] > 0 and heading_array[i] < 90:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                lat_bound_1 = destiny_TBS_lat

            elif heading_array[i] > 90 and heading_array[i] < 180:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] > 180 and heading_array[i] < 270:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] > 270 and heading_array[i] < 360:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                long_bound_2 = destiny_TBS_long

            nod_lat = way_points_out[-1][0]
            nod_lon = way_points_out[-1][1]

            print(f'loop OBC is: {i}')
            print(f'counter for way points: {zh}')

            print(f'Step is: {step}')

            print(f'lat_bound_1: {lat_bound_1}')
            print(f'lat_bound_2: {lat_bound_2}')
            print(f'long_bound_1: {long_bound_1}')
            print(f'long_bound_2: {long_bound_2}')

            for h in range(len(param.zones)):
                if lat_bound_1 >= param.LatitMin[h] - param.dx and lat_bound_2 >= param.LatitMin[h] - param.dx \
                        and lat_bound_1 <= param.LatitMax[h] + param.dx and lat_bound_2 <= param.LatitMax[h] + param.dx \
                        and long_bound_1 >= param.lonMin[h] - param.dx and long_bound_2 >= param.lonMin[h] - param.dx \
                        and long_bound_1 <= param.lonMax[h] + param.dx and long_bound_2 <= param.lonMax[h] + param.dx:
                    print(f'h is : {h}')
                    print(f'file path: {data.file_paths_comp[h]}')
                    ds = xr.open_dataset(data.file_paths_comp[h])
                    time_res_used = param.time_res[h]
                    print(f'time_Res_used: {time_res_used}')
                    break
                else:
                    print(f"not found in : {h} file")
                    if h == param.Zones_num - 1:
                        print("Out of Range, Given Positions are out of Zones")

            # Print the dataset to understand its structure (optional)
            #print(ds)

            sum_TCV_tmp = np.nansum(TCV_array_tmp_1)
            sum_TCV = np.nansum(TCV_array_tmp)

            print(f'divergence array 1 : {div_array_1[i][ml]}')
            print(f'divergence array 1+1: {div_array_1[i][ml+1]}')


            if zh == 0:
                date_str = param.ETD
            # Convert string to datetime object
                date_obj = datetime.fromisoformat(date_str)
            # Extract minutes from datetime object
                minutes = date_obj.minute
                hours = date_obj.hour
                Days = date_obj.day
                months = date_obj.month
                if minutes != 0:
                    result = (60 - minutes) / 60.0 # time in hours
                    modified_date = date_obj - timedelta(minutes= minutes)
                    dist_check = result * v0
                    if dist_check < param.min_dist_per_head:
                        modified_date = (date_obj - timedelta(minutes= minutes)) + timedelta(hours=1)
                        X_Y = 2



            zh += 1




            if sum_TCV_tmp < result and X_Y == 0:
                """
                
                
                X_Z = 2
                start_time_TBS = modified_date.isoformat()
                if time_res_used == 1:
                    original_datetime = datetime.fromisoformat(start_time_TBS)
                    hours_to_add = 0
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS = end_time.isoformat()

                elif time_res_used == 3:
                    original_datetime = datetime.fromisoformat(start_time_TBS)
                    hours_to_add = 2
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS = end_time.isoformat()



                subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2), longitude=slice(long_bound_1, long_bound_2),
                                time=slice(start_time_TBS, end_time_TBS))

                print("minutes are not equal zero")
                print(f'subset: {subset}')


                subset['nod_lat'] = nod_lat
                subset['nod_lon'] = nod_lon
                var1 = subset['latitude']
                var2 = subset['longitude']
                var3 = subset['nod_lat']
                var4 = subset['nod_lon']
                var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                subset['ang_ship'] = var5
                var12 = methods.calc_distance(var3, var4, var1, var2)
                subset['distance'] = var12

                print(f'div array element: {div_array_1[i][ml]}')
                print(f'div array element: {div_array_1[i][ml + 1]}')
                print(subset['ang_ship'])
                filtered_data = subset.where(
                    (subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),
                    drop=True)


                ml +=2
                filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (filtered_data['distance'] >= param.min_dist_per_head) , drop=True)

                var6 = filtered_data_2['ang_ship']
                var7 = filtered_data_2['VMDR']
                var8 = filtered_data_2['VHM0']

                # ship to wave relative direction (VMDR - ang_ship)
                theta = var7 - var6
                filtered_data_2['theta_rel'] = theta

                adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                    filtered_data_2['theta_rel'])
                filtered_data_2['adj_theta'] = adjusted_theta

                df = filtered_data_2.to_dataframe().reset_index()

                cleaned_df = df.dropna(axis=0, how='any')

                cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                cleaned_df_reset = cleaned_df.reset_index(drop=True)

                if cleaned_df_reset.empty:
                    way_points_tuples_var = (np.nan, np.nan)
                    way_points_out.append(way_points_tuples_var)
                    distances_out.append(np.nan)
                    headings_out.append(np.nan)
                    TCV_array_out.append(np.nan)
                    TCV_array_tmp_1.append(np.nan)


                else:

                    if param.math_model == 1:
                        methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out= methods.optimize_without_khokhlov(cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out/v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp_1.append(time_var)

                    elif param.math_model == 2:
                        methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp_1.append(time_var)
                        
                """




            elif sum_TCV < time_res_used:
                if X_Z == 0:
                    start_time = param.ETD

                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()


                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()


                    print(f'start time: {start_time_TBS}')
                    print(f'End Time: {end_time_TBS}')
                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes = 00")
                    print(f'subset: {subset}')

                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    print(subset['nod_lat'])
                    print(subset['nod_lon'])
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var3, var4, var1, var2)
                    subset['distance'] = var12
                    print(subset['distance'])

                    print(f'div array element: {div_array_1[i][ml]}')
                    print(f'div array element: {div_array_1[i][ml + 1]}')
                    print(subset['ang_ship'])



                    
                    filtered_data = subset.where(
                        (subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),
                        drop=True)

                    print(filtered_data['ang_ship'])
                    print(filtered_data['distance'])


                    ml +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head) & (filtered_data['distance'] >= param.min_dist_per_head), drop=True)


                    print(filtered_data_2['distance'])
                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                          filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    print(f'df : {df}')

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    print(f'cleaned_df : {cleaned_df}')

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)

                    print(f'cleaned_df_reset : {cleaned_df_reset}')

                    if cleaned_df_reset.empty:
                        print('chosen waypoints maybe land or empty dataframe due to empty distance or angle ship')

                        if df.empty:
                            print('empty dataframe due to empty distance')


                            filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head), drop=True)

                            print(filtered_data_2['distance'])
                            var6 = filtered_data_2['ang_ship']
                            var7 = filtered_data_2['VMDR']
                            var8 = filtered_data_2['VHM0']

                            # ship to wave relative direction (VMDR - ang_ship)
                            theta = var7 - var6
                            filtered_data_2['theta_rel'] = theta

                            adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0,
                                                      filtered_data_2['theta_rel'] + 360,
                                                      filtered_data_2['theta_rel'])
                            filtered_data_2['adj_theta'] = adjusted_theta

                            df = filtered_data_2.to_dataframe().reset_index()

                            print(f'df : {df}')

                            cleaned_df = df.dropna(axis=0, how='any')

                            cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                            print(f'cleaned_df : {cleaned_df}')

                            cleaned_df_reset = cleaned_df.reset_index(drop=True)

                            print(f'cleaned_df_reset : {cleaned_df_reset}')

                            if param.math_model == 1:
                                methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                                latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                    cleaned_df_reset)
                                way_points_tuples_var = (latitude_out, longitude_out)
                                way_points_out.append(way_points_tuples_var)
                                distances_out.append(displacement_out)
                                headings_out.append(bearing_out)
                                time_var = displacement_out / v0
                                TCV_array_out.append(time_var)
                                TCV_array_tmp.append(time_var)

                            elif param.math_model == 2:
                                methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                                latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                    cleaned_df_reset)
                                way_points_tuples_var = (latitude_out, longitude_out)
                                way_points_out.append(way_points_tuples_var)
                                distances_out.append(displacement_out)
                                headings_out.append(bearing_out)
                                time_var = displacement_out / v0
                                TCV_array_out.append(time_var)
                                TCV_array_tmp.append(time_var)


                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)



                elif X_Z == 2 and X_Y == 0:

                    print("minutes  != 0 and first round of small TCV is over")

                    """
                    
                    
                    start_time = modified_date + timedelta(hours= 1)
                    start_time_tmp = start_time.isoformat()
                    original_datetime = datetime.fromisoformat(start_time_tmp)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()

                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),
                                    longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes  != 0 and first round of small TCV is over")
                    print(f'subset: {subset}')


                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var3, var4, var1, var2)
                    subset['distance'] = var12

                    print(f'div array element: {div_array_1[i][ml]}')
                    print(f'div array element: {div_array_1[i][ml + 1]}')
                    print(subset['ang_ship'])


                    filtered_data = subset.where(
                        (subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),
                        drop=True)
                        
                    

                    ml +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (
                            filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                      filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)


                    if cleaned_df_reset.empty:
                        way_points_tuples_var = (np.nan, np.nan)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(np.nan)
                        headings_out.append(np.nan)
                        TCV_array_out.append(np.nan)
                        TCV_array_tmp.append(np.nan)


                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)
                            
                    """



                elif X_Z == 2 and X_Y == 2:

                    print("minutes  != 0 and first round of small TCV is over")


                    """
                    
                    
                    start_time = modified_date
                    start_time_tmp = start_time.isoformat()
                    original_datetime = datetime.fromisoformat(start_time_tmp)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()

                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),
                                    longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes  != 0 and first round of small TCV is over")
                    print(f'subset: {subset}')



                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var3, var4, var1, var2)
                    subset['distance'] = var12

                    print(f'div array element: {div_array_1[i][ml]}')
                    print(f'div array element: {div_array_1[i][ml + 1]}')
                    print(subset['ang_ship'])

                    filtered_data = subset.where(
                        (subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),
                        drop=True)




                    ml +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (
                            filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                            filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)

                    if cleaned_df_reset.empty:
                        way_points_tuples_var = (np.nan, np.nan)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(np.nan)
                        headings_out.append(np.nan)
                        TCV_array_out.append(np.nan)
                        TCV_array_tmp.append(np.nan)


                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)
                            
                    """


            elif sum_TCV >= time_res_used:

                if time_res_used == 1:
                    H += 1
                    TCV_array_tmp.clear()
                    start_time = start_time_TBS
                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS_1 = start_time_used.isoformat()

                    original_datetime = datetime.fromisoformat(start_time_TBS_1)
                    hours_to_add = 0
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS_1 = end_time.isoformat()

                elif time_res_used == 3:
                    H += 3
                    TCV_array_tmp.clear()
                    start_time = start_time_TBS
                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS_1 = start_time_used.isoformat()

                    original_datetime = datetime.fromisoformat(start_time_TBS_1)
                    hours_to_add = 2
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS_1 = end_time.isoformat()

                # Select the data within the bounding box and time range
                subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2), longitude=slice(long_bound_1, long_bound_2), time=slice(start_time_TBS_1, end_time_TBS_1))

                # Print the subset to check (optional)
                print(f'subset: {subset}')
                print("TCV limit Exceeded")


                subset['nod_lat'] = nod_lat
                subset['nod_lon'] = nod_lon
                var1 = subset['latitude']
                var2 = subset['longitude']
                var3 = subset['nod_lat']
                var4 = subset['nod_lon']
                var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                subset['ang_ship'] = var5
                var12 = methods.calc_distance(var3, var4, var1, var2)
                subset['distance'] = var12

                print(f'div array element: {div_array_1[i][ml]}')
                print(f'div array element: {div_array_1[i][ml + 1]}')
                print(subset['ang_ship'])

                
                
                filtered_data = subset.where(
                    (subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),
                    drop=True)

                print(subset['ang_ship'])

                print(filtered_data['distance'])


                ml +=2
                filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head) & (
                        filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                print(filtered_data_2['distance'])

                var6 = filtered_data_2['ang_ship']
                var7 = filtered_data_2['VMDR']
                var8 = filtered_data_2['VHM0']

                # ship to wave relative direction (VMDR - ang_ship)
                theta = var7 - var6
                filtered_data_2['theta_rel'] = theta

                adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                          filtered_data_2['theta_rel'])
                filtered_data_2['adj_theta'] = adjusted_theta

                df = filtered_data_2.to_dataframe().reset_index()

                print(f'df : {df}')

                cleaned_df = df.dropna(axis=0, how='any')

                cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                print(f'cleaned_df : {cleaned_df}')
                cleaned_df_reset = cleaned_df.reset_index(drop=True)

                print(f'cleaned_df_reset : {cleaned_df_reset}')

                if cleaned_df_reset.empty:
                    print('chosen waypoints maybe land or empty dataframe due to empty distance or angle ship')

                    if df.empty:
                        print('empty dataframe due to empty distance')

                        filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head),
                                                              drop=True)

                        print(filtered_data_2['distance'])
                        var6 = filtered_data_2['ang_ship']
                        var7 = filtered_data_2['VMDR']
                        var8 = filtered_data_2['VHM0']

                        # ship to wave relative direction (VMDR - ang_ship)
                        theta = var7 - var6
                        filtered_data_2['theta_rel'] = theta

                        adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0,
                                                  filtered_data_2['theta_rel'] + 360,
                                                  filtered_data_2['theta_rel'])
                        filtered_data_2['adj_theta'] = adjusted_theta

                        df = filtered_data_2.to_dataframe().reset_index()

                        print(f'df : {df}')

                        cleaned_df = df.dropna(axis=0, how='any')

                        cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                        print(f'cleaned_df : {cleaned_df}')

                        cleaned_df_reset = cleaned_df.reset_index(drop=True)

                        print(f'cleaned_df_reset : {cleaned_df_reset}')

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)


                else:

                    if param.math_model == 1:
                        methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                            cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp.append(time_var)

                    elif param.math_model == 2:
                        methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                            cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp.append(time_var)

            step +=1

            print(f': way points out : {way_points_out}')
            print(f': distances : {distances_out}')
            print(f': headings: {headings_out}')
            print(f': TCV : {TCV_array_out}')

            print(f' counter OBC: {i}')
            print(f' counter Conv: {nl}')

            result = 0



        for step in range(0, int(conv_array_1.shape[1]/2)-1, 1):

            print(f'counter for way points: {zh}')

            if heading_array[i] == 0:
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] == 90:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] == 180:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] == 270:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][1])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] > 0 and heading_array[i] < 90:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                lat_bound_1 = destiny_TBS_lat

            elif heading_array[i] > 90 and heading_array[i] < 180:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                long_bound_2 = destiny_TBS_long

            elif heading_array[i] > 180 and heading_array[i] < 270:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                long_bound_2 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                lat_bound_2 = destiny_TBS_lat

            elif heading_array[i] > 270 and heading_array[i] < 360:

                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][0])
                lat_bound_1 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][1])
                long_bound_1 = destiny_TBS_long
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head_diagonal,TBS_array[i][2])
                lat_bound_2 = destiny_TBS_lat
                destiny_TBS_lat, destiny_TBS_long = methods.destination_point(way_points_out[-1][0], way_points_out[-1][1], param.dist_per_head,TBS_array[i][3])
                long_bound_2 = destiny_TBS_long

            nod_lat = way_points_out[-1][0]
            nod_lon = way_points_out[-1][1]

            print(f'loop OBC is: {i}')
            print(f'counter for way points: {zh}')

            print(f'Step is: {step}')

            print(f'lat_bound_1: {lat_bound_1}')
            print(f'lat_bound_2: {lat_bound_2}')
            print(f'long_bound_1: {long_bound_1}')
            print(f'long_bound_2: {long_bound_2}')

            for h in range(0, param.Zones_num, 1):
                if lat_bound_1 >= param.LatitMin[h] - param.dx and lat_bound_2 >= param.LatitMin[h] - param.dx \
                        and lat_bound_1 <= param.LatitMax[h] + param.dx and lat_bound_2 <= param.LatitMax[h] + param.dx \
                        and long_bound_1 >= param.lonMin[h] - param.dx and long_bound_2 >= param.lonMin[h] - param.dx \
                        and long_bound_1 <= param.lonMax[h] + param.dx and long_bound_2 <= param.lonMax[h] + param.dx:
                    print(f'h is : {h}')
                    print(f'file path: {data.file_paths_comp[h]}')
                    ds = xr.open_dataset(data.file_paths_comp[h])
                    time_res_used = param.time_res[h]
                    print(f'time_Res_used: {time_res_used}')
                    break
                else:
                    print(f"not found in : {h} file")
                    if h == param.Zones_num - 1:
                        print("Out of Range, Given Positions are out of Zones")

            # Print the dataset to understand its structure (optional)
            #print(ds)

            sum_TCV_tmp = np.nansum(TCV_array_tmp_1)
            sum_TCV = np.nansum(TCV_array_tmp)

            print(f'convergence array 1 : {conv_array_1[i][nl]}')
            print(f'convergence array 1+1 : {conv_array_1[i][nl+1]}')


            zh += 1
            if sum_TCV_tmp < result and X_Y == 0:

                """
                
                
                X_Z = 2
                start_time_TBS = modified_date.isoformat()
                if time_res_used == 1:
                    original_datetime = datetime.fromisoformat(start_time_TBS)
                    hours_to_add = 0
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS = end_time.isoformat()

                elif time_res_used == 3:
                    original_datetime = datetime.fromisoformat(start_time_TBS)
                    hours_to_add = 2
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS = end_time.isoformat()

                subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2), longitude=slice(long_bound_1, long_bound_2),
                                time=slice(start_time_TBS, end_time_TBS))

                print("minutes are not equal zero")
                print(subset)


                subset['nod_lat'] = nod_lat
                subset['nod_lon'] = nod_lon

                var1 = subset['latitude']
                var2 = subset['longitude']
                var3 = subset['nod_lat']
                var4 = subset['nod_lon']
                var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                subset['ang_ship'] = var5
                var12 = methods.calc_distance(var3, var4, var1, var2)
                subset['distance'] = var12

                filtered_data = subset.where((subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),drop=True)
                ml +=2
                filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (filtered_data['distance'] >= param.min_dist_per_head) , drop=True)

                var6 = filtered_data_2['ang_ship']
                var7 = filtered_data_2['VMDR']
                var8 = filtered_data_2['VHM0']

                # ship to wave relative direction (VMDR - ang_ship)
                theta = var7 - var6
                filtered_data_2['theta_rel'] = theta

                adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                    filtered_data_2['theta_rel'])
                filtered_data_2['adj_theta'] = adjusted_theta

                df = filtered_data_2.to_dataframe().reset_index()

                cleaned_df = df.dropna(axis=0, how='any')

                cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                cleaned_df_reset = cleaned_df.reset_index(drop=True)

                if cleaned_df_reset.empty:
                    way_points_tuples_var = (np.nan, np.nan)
                    way_points_out.append(way_points_tuples_var)
                    distances_out.append(np.nan)
                    headings_out.append(np.nan)
                    TCV_array_out.append(np.nan)
                    TCV_array_tmp_1.append(np.nan)

                else:

                    if param.math_model == 1:
                        methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out= methods.optimize_without_khokhlov(cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out/v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp_1.append(time_var)

                    elif param.math_model == 2:
                        methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp_1.append(time_var)
                
                
                """



            elif sum_TCV < time_res_used:
                if X_Z == 0:
                    start_time = param.ETD

                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()


                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    print(f'start time: {start_time_TBS}')
                    print(f'End Time: {end_time_TBS}')

                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes = 00")
                    #print(subset)

                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    print(subset['nod_lat'])
                    print(subset['nod_lon'])
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var1, var2, var3, var4)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var1, var2, var3, var4)
                    subset['distance'] = var12


                    print(f'conv array element: {conv_array_1[i][nl]}')
                    print(f'conv array element: {conv_array_1[i][nl + 1]}')
                    print(subset['ang_ship'])
                    filtered_data = subset.where(
                        (subset['ang_ship'] <= conv_array_1[i][nl]) & (subset['ang_ship'] >= conv_array_1[i][nl + 1]),
                        drop=True)

                    print(filtered_data['ang_ship'])
                    print(filtered_data['distance'])

                    


                    nl +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head) & (filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                    print(filtered_data_2['distance'])

                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                          filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    print(f'df : {df}')

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    print(f'cleaned_df : {cleaned_df}')

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)

                    print(f'cleaned_df_reset : {cleaned_df_reset}')

                    if cleaned_df_reset.empty:
                        print('chosen waypoints maybe land or empty dataframe due to empty distance or angle ship')

                        if df.empty:
                            print('empty dataframe due to empty distance')

                            filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head),
                                                                  drop=True)

                            print(filtered_data_2['distance'])
                            var6 = filtered_data_2['ang_ship']
                            var7 = filtered_data_2['VMDR']
                            var8 = filtered_data_2['VHM0']

                            # ship to wave relative direction (VMDR - ang_ship)
                            theta = var7 - var6
                            filtered_data_2['theta_rel'] = theta

                            adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0,
                                                      filtered_data_2['theta_rel'] + 360,
                                                      filtered_data_2['theta_rel'])
                            filtered_data_2['adj_theta'] = adjusted_theta

                            df = filtered_data_2.to_dataframe().reset_index()

                            print(f'df : {df}')

                            cleaned_df = df.dropna(axis=0, how='any')

                            cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                            print(f'cleaned_df : {cleaned_df}')

                            cleaned_df_reset = cleaned_df.reset_index(drop=True)

                            print(f'cleaned_df_reset : {cleaned_df_reset}')

                            if param.math_model == 1:
                                methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                                latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                    cleaned_df_reset)
                                way_points_tuples_var = (latitude_out, longitude_out)
                                way_points_out.append(way_points_tuples_var)
                                distances_out.append(displacement_out)
                                headings_out.append(bearing_out)
                                time_var = displacement_out / v0
                                TCV_array_out.append(time_var)
                                TCV_array_tmp.append(time_var)

                            elif param.math_model == 2:
                                methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                                latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                    cleaned_df_reset)
                                way_points_tuples_var = (latitude_out, longitude_out)
                                way_points_out.append(way_points_tuples_var)
                                distances_out.append(displacement_out)
                                headings_out.append(bearing_out)
                                time_var = displacement_out / v0
                                TCV_array_out.append(time_var)
                                TCV_array_tmp.append(time_var)


                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)






                elif X_Z == 2 and X_Y == 0:

                    """
                    start_time = modified_date + timedelta(hours= 1)
                    start_time_tmp = start_time.isoformat()
                    original_datetime = datetime.fromisoformat(start_time_tmp)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()

                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),
                                    longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes  != 0 and first round of small TCV is over")
                    print(subset)


                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var3, var4, var1, var2)
                    subset['distance'] = var12


                    filtered_data = subset.where((subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),drop=True)

                    ml +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (
                            filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                      filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)


                    if cleaned_df_reset.empty:
                        way_points_tuples_var = (np.nan, np.nan)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(np.nan)
                        headings_out.append(np.nan)
                        TCV_array_out.append(np.nan)
                        TCV_array_tmp.append(np.nan)

                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)
                            
                    """




                elif X_Z == 2 and X_Y == 2:

                    """
                    start_time = modified_date
                    start_time_tmp = start_time.isoformat()
                    original_datetime = datetime.fromisoformat(start_time_tmp)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS = start_time_used.isoformat()

                    if time_res_used == 1:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 0
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    elif time_res_used == 3:
                        original_datetime = datetime.fromisoformat(start_time_TBS)
                        hours_to_add = 2
                        minutes_to_add = 59
                        seconds_to_add = 59
                        delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                        end_time = original_datetime + delta
                        end_time_TBS = end_time.isoformat()

                    subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2),
                                    longitude=slice(long_bound_1, long_bound_2),
                                    time=slice(start_time_TBS, end_time_TBS))

                    print("minutes  != 0 and first round of small TCV is over")
                    print(subset)



                    subset['nod_lat'] = nod_lat
                    subset['nod_lon'] = nod_lon
                    var1 = subset['latitude']
                    var2 = subset['longitude']
                    var3 = subset['nod_lat']
                    var4 = subset['nod_lon']
                    var5 = methods.calculate_initial_bearing(var3, var4, var1, var2)
                    subset['ang_ship'] = var5
                    var12 = methods.calc_distance(var3, var4, var1, var2)
                    subset['distance'] = var12



                    filtered_data = subset.where((subset['ang_ship'] >= div_array_1[i][ml]) & (subset['ang_ship'] <= div_array_1[i][ml + 1]),drop=True)

                    ml +=2
                    filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.dist_per_head) & (
                            filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                    var6 = filtered_data_2['ang_ship']
                    var7 = filtered_data_2['VMDR']
                    var8 = filtered_data_2['VHM0']

                    # ship to wave relative direction (VMDR - ang_ship)
                    theta = var7 - var6
                    filtered_data_2['theta_rel'] = theta

                    adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                            filtered_data_2['theta_rel'])
                    filtered_data_2['adj_theta'] = adjusted_theta

                    df = filtered_data_2.to_dataframe().reset_index()

                    cleaned_df = df.dropna(axis=0, how='any')

                    cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                    cleaned_df_reset = cleaned_df.reset_index(drop=True)

                    if cleaned_df_reset.empty:
                        way_points_tuples_var = (np.nan, np.nan)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(np.nan)
                        headings_out.append(np.nan)
                        TCV_array_out.append(np.nan)
                        TCV_array_tmp.append(np.nan)

                    else:

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)
                            
                    """


            elif sum_TCV >= time_res_used:

                if time_res_used == 1:
                    H += 1
                    TCV_array_tmp.clear()
                    start_time = start_time_TBS
                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS_1 = start_time_used.isoformat()

                    original_datetime = datetime.fromisoformat(start_time_TBS_1)
                    hours_to_add = 0
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS_1 = end_time.isoformat()

                elif time_res_used == 3:
                    H += 3
                    TCV_array_tmp.clear()
                    start_time = start_time_TBS
                    original_datetime = datetime.fromisoformat(start_time)
                    hours_to_add = H
                    # Create a timedelta object for the hours and minutes to add
                    delta = timedelta(hours=hours_to_add)

                    # Add the timedelta to the original datetime
                    start_time_used = original_datetime + delta

                    # Format the modified datetime back into ISO format
                    start_time_TBS_1 = start_time_used.isoformat()

                    original_datetime = datetime.fromisoformat(start_time_TBS_1)
                    hours_to_add = 2
                    minutes_to_add = 59
                    seconds_to_add = 59
                    delta = timedelta(hours=hours_to_add, minutes=minutes_to_add, seconds=seconds_to_add)
                    end_time = original_datetime + delta
                    end_time_TBS_1 = end_time.isoformat()

                # Select the data within the bounding box and time range
                subset = ds.sel(latitude=slice(lat_bound_1, lat_bound_2), longitude=slice(long_bound_1, long_bound_2), time=slice(start_time_TBS_1, end_time_TBS_1))

                # Print the subset to check (optional)
                print(f'subset: {subset}')
                print("TCV limit Exceeded")

                subset['nod_lat'] = nod_lat
                subset['nod_lon'] = nod_lon
                var1 = subset['latitude']
                var2 = subset['longitude']
                var3 = subset['nod_lat']
                var4 = subset['nod_lon']
                var5 = methods.calculate_initial_bearing(var1, var2, var3, var4)
                subset['ang_ship'] = var5
                var12 = methods.calc_distance(var1, var2, var3, var4)
                subset['distance'] = var12

                print(subset['distance'])

                print(f'conv array element: {conv_array_1[i][nl]}')
                print(f'conv array element: {conv_array_1[i][nl + 1]}')
                print(subset['ang_ship'])

                filtered_data = subset.where(
                    (subset['ang_ship'] <= conv_array_1[i][nl]) & (subset['ang_ship'] >= conv_array_1[i][nl + 1]),
                    drop=True)

                print(subset['ang_ship'])

                print(filtered_data['distance'])

                nl +=2
                filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head) & (filtered_data['distance'] >= param.min_dist_per_head), drop=True)

                print(filtered_data_2['distance'])
                var6 = filtered_data_2['ang_ship']
                var7 = filtered_data_2['VMDR']
                var8 = filtered_data_2['VHM0']

                # ship to wave relative direction (VMDR - ang_ship)
                theta = var7 - var6
                filtered_data_2['theta_rel'] = theta

                adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0, filtered_data_2['theta_rel'] + 360,
                                          filtered_data_2['theta_rel'])
                filtered_data_2['adj_theta'] = adjusted_theta

                df = filtered_data_2.to_dataframe().reset_index()

                print(f'df : {df}')

                cleaned_df = df.dropna(axis=0, how='any')

                cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                print(f'cleaned_df : {cleaned_df}')
                cleaned_df_reset = cleaned_df.reset_index(drop=True)

                print(f'cleaned_df_reset : {cleaned_df_reset}')

                if cleaned_df_reset.empty:
                    print('chosen waypoints maybe land or empty dataframe due to empty distance or angle ship')

                    if df.empty:
                        print('empty dataframe due to empty distance')

                        filtered_data_2 = filtered_data.where((filtered_data['distance'] <= param.max_dist_per_head),
                                                              drop=True)

                        print(filtered_data_2['distance'])
                        var6 = filtered_data_2['ang_ship']
                        var7 = filtered_data_2['VMDR']
                        var8 = filtered_data_2['VHM0']

                        # ship to wave relative direction (VMDR - ang_ship)
                        theta = var7 - var6
                        filtered_data_2['theta_rel'] = theta

                        adjusted_theta = xr.where(filtered_data_2['theta_rel'] < 0,
                                                  filtered_data_2['theta_rel'] + 360,
                                                  filtered_data_2['theta_rel'])
                        filtered_data_2['adj_theta'] = adjusted_theta

                        df = filtered_data_2.to_dataframe().reset_index()

                        print(f'df : {df}')

                        cleaned_df = df.dropna(axis=0, how='any')

                        cleaned_df['adj_theta_rad'] = cleaned_df['adj_theta'].apply(np.radians)

                        print(f'cleaned_df : {cleaned_df}')

                        cleaned_df_reset = cleaned_df.reset_index(drop=True)

                        print(f'cleaned_df_reset : {cleaned_df_reset}')

                        if param.math_model == 1:
                            methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)

                        elif param.math_model == 2:
                            methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                            latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                                cleaned_df_reset)
                            way_points_tuples_var = (latitude_out, longitude_out)
                            way_points_out.append(way_points_tuples_var)
                            distances_out.append(displacement_out)
                            headings_out.append(bearing_out)
                            time_var = displacement_out / v0
                            TCV_array_out.append(time_var)
                            TCV_array_tmp.append(time_var)


                else:

                    if param.math_model == 1:
                        methods.without_khokhlov(cleaned_df_reset, v0, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_without_khokhlov(
                            cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp.append(time_var)

                    elif param.math_model == 2:
                        methods.with_khokhlov(cleaned_df_reset, v0, DWT, LBP)
                        latitude_out, longitude_out, bearing_out, displacement_out = methods.optimize_with_khokhlov(
                            cleaned_df_reset)
                        way_points_tuples_var = (latitude_out, longitude_out)
                        way_points_out.append(way_points_tuples_var)
                        distances_out.append(displacement_out)
                        headings_out.append(bearing_out)
                        time_var = displacement_out / v0
                        TCV_array_out.append(time_var)
                        TCV_array_tmp.append(time_var)

            step += 1

            print(f': way points out : {way_points_out}')
            print(f': distances : {distances_out}')
            print(f': headings: {headings_out}')
            print(f': TCV : {TCV_array_out}')

            print(f' counter OBC: {i}')
            print(f' counter Conv: {nl}')

            result = 0

    way_points_tuples = (lat11, long11)
    way_points_out.append(way_points_tuples)
    # after convergence for loop, you will append latitude 11 and longitude 11 to coordinates and calculate its distance, heading and TCV and append
    # them.


# out of OBCS loop you will print all 4 output arrays



    distances_to_be_printed = []
    headings_to_be_printed = []
    TCV_to_be_printed = []

    for i in range(0,len(way_points_out)-1,1):
        lat_x = way_points_out[i][0]
        lon_y = way_points_out[i][1]
        lat_z = way_points_out[i + 1][0]
        lon_l = way_points_out[i + 1][1]
        distance_xh = methods.calc_distance(lat_x, lon_y, lat_z, lon_l)
        distances_to_be_printed.append(distance_xh)
        bearing_xh = methods.calculate_initial_bearing(lat_x,lon_y, lat_z, lon_l)
        headings_to_be_printed.append(bearing_xh)
        time_var_xh= distance_xh/v0
        TCV_to_be_printed.append(time_var_xh)





    print(f' way points out : {way_points_out}')
    print(f'number of elements in way points array : {len(way_points_out)}')
    print(f' distances : {distances_to_be_printed}')
    print(f'number of elements in distances array : {len(distances_to_be_printed)}')
    print(f' headings: {headings_to_be_printed}')
    print(f'number of elements in headings array : {len(headings_to_be_printed)}')
    print(f' TCV : {TCV_to_be_printed}')
    print(f'number of elements in TCV array : {len(TCV_to_be_printed)}')

    total_distance = sum(distances_to_be_printed)
    print(f"Total Voyage Distance: {total_distance}")
    total_time = sum(TCV_to_be_printed)
    print(f"Total Sailing Time: {total_time}")

    with open(
            f'{param.dir_out_phaseII}'+ '\Simulation_Report_' + param.name_Simu + '.txt',
            'w') as file:

        st = '              PHASE II Simulation Report:     ' + param.name_Simu + '\n'
        file.write(st)
        st = '=========================================================================\n\n'
        file.write(st)

        st = '   (Lat1,Long1)         (Lat2,Long2)             heading            distance           time    \n\n'
        file.write(st)

        for z in range(0, len(way_points_out)-1, 1):
            file.write('{:8.3f},{:8.3f}     {:8.3f},{:8.3f}         {:8.3f}          {:8.3f}          {:8.3f} \n'.format(way_points_out[z][0],way_points_out[z][1], way_points_out[z+1][0],way_points_out[z+1][1], headings_to_be_printed[z], distances_to_be_printed[z], TCV_to_be_printed[z]))

        file.write('\n')
        file.write(f'Number of Way points: {len(way_points_out)}\n\n')
        file.write(f'Total Voyage Distance: {total_distance}\n\n')
        file.write(f'Total Sailing Time: {total_time}\n\n')


        file.write("--------------------End of Report----------------------")



    with open(
            f'{param.dir_out_phaseII}'+ '\Simulation_Metadata_' + param.name_Simu + '.txt',
            'w') as file:

        st = '              Simulation MetaData:     ' + param.name_Simu + '\n'
        file.write(st)
        st = '=========================================================================\n\n'
        file.write(st)

        file.write('\n')
        file.write(f'Simulation name : {param.name_Simu}\n\n')
        file.write(f'Paths number : {param.paths_num}\n\n')
        file.write(f'Paths coordinates : {param.paths_coord}\n\n')
        file.write(f'Estimated Time Departure : {param.ETD}\n\n')
        file.write(f'Longitude Minimum : {param.lonMin}\n\n')
        file.write(f'Longitude Maximum : {param.lonMax}\n\n')
        file.write(f'Latitude Minimum : {param.LatitMin}\n\n')
        file.write(f'Latitude Maximum : {param.LatitMax}\n\n')
        file.write(f'Distance per Heading : {param.dist_per_head}\n\n')
        file.write(f'Reduction factor : {param.reduction_factor}\n\n')
        file.write(f'Gain factor : {param.gain_factor}\n\n')
        file.write(f'Degree Steps : {param.deg_step}\n\n')
        file.write(f'OBC diameter : {param.OBC_diameter}\n\n')
        file.write(f'Math Model Used : {param.math_model}\n\n')
        file.write(f'Service Speed : {param.serv_speed}\n\n')
        file.write(f'Max. Speed : {param.Max_speed}\n\n')
        file.write(f'TBS radius : {param.TBS_radius}\n\n')
        file.write(f'LBP : {param.LBP}\n\n')
        file.write(f'DWT : {param.DWT}\n\n')

        file.write("--------------------End of Report----------------------")









    data_array = np.array(way_points_out)

    # Create a mask to filter out tuples with NaN values
    mask = ~np.isnan(data_array).any(axis=1)

    # Apply the mask to get the cleaned array of tuples
    cleaned_data = data_array[mask]

    # Convert the cleaned array back to a list of tuples if needed
    cleaned_data_tuples = [tuple(x) for x in cleaned_data]

    print(cleaned_data_tuples)

    methods.plot_path2(cleaned_data_tuples)






if __name__ == "__main__":
    spinner = Spinner()

    print("Start calculating.......")
    spinner.start()

    main()

    spinner.stop()
    print("calculations completed!")







