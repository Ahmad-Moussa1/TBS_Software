import parameters as param
import methods
import sys
import time
import threading

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
    obcd_check_2 = []

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
    print('1st check is done')
    print(f'number of paths is: {len(paths_number)}')
    # print(len(obcd_check_1))
    destiny_latitude_final = param.paths_coord[-1][0]
    destiny_longitude_final = param.paths_coord[-1][1]
    destiny_tuple_final = (destiny_latitude_final, destiny_longitude_final)
    outer_boundary.append(destiny_tuple_final)

    print(outer_boundary)
    # print(len(outer_boundary))

    distance4 = []
    radius4 = []
    bearing4 = []
    centers_OBC = []
    for i in range(0, len(outer_boundary) - 1, 1):
        lat3 = outer_boundary[i][0]
        lon3 = outer_boundary[i][1]
        lat4 = outer_boundary[i + 1][0]
        lon4 = outer_boundary[i + 1][1]
        distance3 = methods.calc_distance(lat3, lon3, lat4, lon4)
        distance4.append(distance3)
        radius3 = distance3 / 2
        radius4.append(radius3)
        bearing3 = methods.calculate_initial_bearing(lat3, lon3, lat4, lon4)
        bearing4.append(bearing3)
        center_latit, center_long = methods.destination_point(outer_boundary[i][0], outer_boundary[i][1], radius4[i],
                                                              bearing4[i])
        center_tuple = (center_latit, center_long)
        centers_OBC.append(center_tuple)
        i += 1

    # print(radius4)
    # print(bearing4)
    # print(centers_OBC)
    # print(len(distance4))
    distance4_num = []
    for i in range(1, len(distance4) + 1, 1):
        distance4_num.append(i)

    for i in range(0, len(distance4), 1):
        if distance4[i] > param.TBS_radius:
            """print("last OBC diameter is acceptable")"""
            obcd_check_2.append(" OBC diameter 2nd Check is acceptable")
        else:
            """print("last OBC in each path need to be checked because its diameter is smaller than TBS. \n Please optimize distance/heading or deg. steps ")"""
            obcd_check_2.append(
                "last OBC in each path need to be checked because its diameter is smaller than TBS. \n Please optimize distance/heading or deg. steps")
    print('2nd Check is done')
    print(f'number of OBCs is: {len(outer_boundary) - 1}')
    # you have to check if all outer boundary coordinates (Start , End) not on land areas before plotting.
    status_print = []
    for i in range(0, len(outer_boundary), 1):
        status_case = methods.check_land(outer_boundary[i][0], outer_boundary[i][1], param.ETD)
        status_print.append(status_case)
    print('3rd check is done')
    print(f'number of OBCs coordinates is: {len(outer_boundary)}')
    print(status_print)
    # if you got out of zones points, adjust lons and latitudes in parameters and download another data_set till you got status for each point

    with open(
            f'{param.dir_out_phaseI}'+ '\Report_check_' + param.name_Simu + '.txt',
            'w') as file:

        st = '              PHASE I report:     ' + param.name_Simu + '\n'
        file.write(st)
        st = '=========================================================================\n\n'
        file.write(st)

        st = '1st Check: OBC Diameter is greater than TBS radius && smaller than distance between paths coordinates\n' \
             '** number of checks equal number of paths \n' \
             '** All values must be "OBC Diameter is acceptable"\n\n'
        file.write(st)
        for i, value in zip(paths_number, obcd_check_1):
            file.write(f'path {i}: {value}\n')

        st = '\n 2nd Check: Check OBC Diameter (specially Last OBC in each path) is smaller than TBS radius\n' \
             '** number of checks equal number of OBCs\n' \
             '** All Values must be "Last OBC Diameter is acceptable"\n\n'
        file.write(st)
        for i, value in zip(distance4_num, obcd_check_2):
            file.write(f'OBC number {i}: {value}\n')

        st = '\n 3rd Check: OBC coordinates status (is_land OR is_water) \n' \
             '** number of checks equal number of OBCs +1 \n' \
             '** All values of status must be is_water\n\n'
        file.write(st)

        st = 'Latitude      Longitude       Status    \n'
        file.write(st)
        for i in range(0, len(outer_boundary), 1):
            file.write('{:8.3f}     {:8.3f}         {:9} \n'.format(outer_boundary[i][0], outer_boundary[i][1],
                                                                    str(status_print[i])))

        file.write("--------------------End of Report----------------------")

    outer_bound_plot = []
    for i in range(0, len(outer_boundary) - 1, 1):
        for bear in range(0, 361, 1):
            outer_bound_lat, outer_bound_lon = methods.destination_point(centers_OBC[i][0], centers_OBC[i][1],
                                                                         radius4[i], bear)
            outer_bound_plot.append((outer_bound_lat, outer_bound_lon))

    # print(outer_bound_plot)
    methods.plot_path(outer_bound_plot, param.paths_coord)

    #print(f'outer_bound: {outer_bound_plot}')
    #print(f'paths_coord: {param.paths_coord}')

    print('Done Plotting')
    print('End of Phase I')
    #return outer_boundary

if __name__ == "__main__":
    spinner = Spinner()

    print("Start calculating.......")
    spinner.start()

    main()

    spinner.stop()
    print("calculations completed!")



