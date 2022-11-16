import WaterShed

import pandas as pd
import numpy as np
from numpy import unravel_index
import networkx as nx
from tqdm import tqdm
from os import path



class Graph():

    def __init__(self, dem=None, fdir=None, acc=None, compression=1):
        self.dem = dem
        self.fdir = fdir
        self.acc = acc
        self.compression = compression


    def compute_height(self, file_path, DEMs_path, min_acc, save_path):
        df = pd.read_csv(file_path, sep=';', decimal=',')
        df['x_lon_int'] = df['x_lon'].values.astype(int)
        df['y_lat_int'] = df['y_lat'].values.astype(int)
        
        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis = 0, ascending = True, inplace = True, na_position = "first")
        
        self.df_new = pd.DataFrame(columns=['hstation_id', 'x_lon', 'y_lat', 'height', 'error'])
        x_lon_past, y_lat_past = None, None

        for i, row in df.iterrows():
            print(f'{i+1}/{df.shape[0]} hydropost...')
            hstation_id  = row[0]
            x_lon, y_lat = row[1], row[2]
            coordinate = (x_lon, y_lat)

            # Define coordinate of map to download 
            lng_num, lat_num = int(x_lon), int(y_lat)

            # Check if this coordinates weren't calculated
            if (x_lon_past != lng_num) or (y_lat_past != lat_num):
                x_lon_past, y_lat_past = lng_num, lat_num
                # Set acc values as None to calulate them later
                self.acc_slice = None
                self.acc_Graph = None

                self.tif_pathes = []
                for i in range(lat_num-1, lat_num+2):
                    for j in range(lng_num-1, lng_num+2):
                        lat = str(i)
                        lng = ''.join((['0'] + list(str(int(j))))[-3:])
                        self.tif_pathes.append(f'{DEMs_path}/n{lat}_e{lng}_1arc_v3.tif')
                        
                # check if files 'exisits'
                success_list = []
                for tif_path in self.tif_pathes:
                    if path.exists(tif_path) == False:
                        print(f'{tif_path} is not exist in path {DEMs_path}')
                        success_list.append(False)

                # Download DEM and preprocess it
                if len(success_list) == 0:
                    self.compute_DEM(self.tif_pathes, lng_num, lat_num)
                else:
                    # Temporary while I'm thinking what to with others frames DEMs
                    self.compression = 1
                    self.acc = None
                    self.dem = None
                    self.fdir = None
            

            # Calculate Heights
            top_left = (lng_num-1, lat_num+2) if len(self.tif_pathes) == 9 else (lng_num, lat_num+1)
            bottom_right = (lng_num+2, lat_num-1) if len(self.tif_pathes) == 9 else (lng_num+1, lat_num)

            if self.dem.all() != None:
                height, success = self.compute_height_differance(coordinate, top_left, bottom_right, 10000, min_acc)
                error = 1 if success == False else 0

                dct = {
                    'hstation_id': hstation_id, 
                    'x_lon':x_lon, 
                    'y_lat':y_lat, 
                    'height':height, 
                    'error':error
                }
                self.df_new = self.df_new.append(dct, ignore_index=True)
                self.df_new.to_csv(f'{save_path}/hydroposts_height_calculated.csv', sep=';')

    

    def compute_DEM(self, pathes, lng_num, lat_num):
        if len(pathes) == 9:
            shed = WaterShed.WaterSheds(files_pathes=pathes, compute_acc=True, compression=3)
            self.compression = 3
            self.acc = shed.acc
            self.dem = shed.dem
            self.fdir = shed.fdir
        else:
            #  Use only central tile
            lat = str(lat_num)
            lng = ''.join((['0'] + list(str(int(lng_num))))[-3:])
            shed = WaterShed.WaterSheds(file_path=f'n{lat}_e{lng}_1arc_v3.tif', compute_acc=True)
            self.compression = 1
            self.acc = shed.acc
            self.dem = shed.dem
            self.fdir = shed.fdir


    def compute_height_differance(self, coordinate, top_left, bottom_right, lenth, min_acc):
        # In case to not create acc_graph every time for same lon & lat
        if (self.acc_slice == None) or (self.acc_Graph == None):
            acc_slice = self.acc.copy()
            # Filter river cells
            self.create_acc_graph(acc_slice, self.fdir, min_acc)

        point = self.coordinate2point(coordinate, top_left, bottom_right)
        river_pathes_nodes_DOWN, river_pathes_nodes_UP, success = self.compute_river_path(point, lenth)

        start, end = river_pathes_nodes_DOWN[-1], river_pathes_nodes_UP[-1]
        start_height, end_height = self.dem[start], self.dem[end]
        return abs(start_height - end_height), success


    def compute_river_path(self, point, lenth):
        """
        Function returns all river nodes down and up from the given point.
        * lenth - in meteres from start to up and down.
        """
        point_lenth = self.compression * 35       # around 35 meteres shape for each point

        # DOWN
        river_pathes_lenght_DOWN = int(lenth/point_lenth)        # количество затопленных клеток реки вниз по течению
        river_pathes_nodes_DOWN = [point]        # Стартовая точка тоже входит.
        node = point
        for _ in tqdm(range(river_pathes_lenght_DOWN)):
            try:
                new_node = self.out_node(node)
                river_pathes_nodes_DOWN.append(new_node)
                node = new_node
            except:
                print(f"Out nodes not definded.\nLast node: {river_pathes_nodes_DOWN[-1]}\nCollected {len(river_pathes_nodes_DOWN)} nodes")
                break
        
        # UP
        river_pathes_lenght_UP = int(lenth/point_lenth)        # количество затопленных клеток реки вверх по течению
        river_pathes_nodes_UP = [point]       # Стартовая точка тоже входит.
        node = point
        for _ in tqdm(range(river_pathes_lenght_UP)):
            try:
                new_node = self.in_node(node)
                river_pathes_nodes_UP.append(new_node)
                node = new_node
            except:
                print(f"Out nodes not definded.\nLast node: {river_pathes_nodes_UP[-1]}\nCollected {len(river_pathes_nodes_UP)} nodes")
                break
                
        success = True if (len(river_pathes_nodes_DOWN) == river_pathes_lenght_DOWN + 1) and (len(river_pathes_nodes_UP) == river_pathes_lenght_UP + 1) else False
        return river_pathes_nodes_DOWN, river_pathes_nodes_UP, success


    def coordinate2point(self, coordinate, top_left, bottom_right):
        # lng - horizontal, lat - vertical
        # shape[0] - vertical 
        # shape[1] - horizontal
        lng, lat = coordinate
        lng_left, lng_right = top_left[0], bottom_right[0]
        lat_top, lat_bottom = top_left[1], bottom_right[1]
        lng = abs(lng_left - lng) / abs(lng_left - lng_right)
        lat = 1 - (abs(lat_bottom - lat) / abs(lat_top - lat_bottom))

        x_path, y_path = 1/self.dem.shape[1], 1/self.dem.shape[0]
        x_coordinate = int(round(lng/x_path, 0))
        y_coordinate = int(round(lat/y_path, 0))

        # Check that coordinate contains right cell
        # There can be error with rounding
        if self.accumulation_for_point((y_coordinate, x_coordinate)) < 1000:
            acc_near = self.acc[
                y_coordinate-4 : y_coordinate+5,
                x_coordinate-4 : x_coordinate+5
            ]
            max_index = unravel_index(acc_near.argmax(), acc_near.shape)
            y_coordinate = y_coordinate + max_index[0]-4
            x_coordinate = x_coordinate + max_index[1]-4
            return (y_coordinate, x_coordinate)
        else:
            return (y_coordinate, x_coordinate)


    def create_acc_graph(self, acc_slice, dir_slice, min_acc):
        acc_slice[acc_slice < min_acc] = 0
        self.acc_slice = acc_slice # For future nodes matching
        self.acc_Graph = nx.DiGraph()
        for i in tqdm(range(1, acc_slice.shape[0]-1)):
            for j in range(1, acc_slice.shape[1]-1):
                if acc_slice[i, j] != 0:
                    dir = dir_slice[i, j]
                    dir = dir if dir >= 1 else 1
                    start = (i, j)
                    target = self.fdir_coordinate(start, dir)
                    #print(start, target)
                    self.acc_Graph.add_edge(start, target)


    def fdir_coordinate(self, point, dir):
        (row, column) = point
        if dir == 64: #Up
            return (row - 1, column)
        elif dir == 128: #Up-Right
            return (row - 1, column + 1)
        elif dir == 1: #Right
            return (row, column + 1)
        elif dir == 2: #Down-Right
            return (row + 1, column + 1)
        elif dir == 4: #Down
            return (row + 1, column)
        elif dir == 8: #Down-Left
            return (row + 1, column - 1)
        elif dir == 16: #Left
            return (row, column - 1)
        elif dir == 32: #Up-Left
            return (row - 1, column - 1)

    
    def out_node(self, node):
        return [node for node in self.acc_Graph.edges(node)][0][1]

    def out_node_G(self, node):
        return [node for node in self.G.edges(node)][0][1]

    def in_node(self, node):
        in_nodes = [edge[0] for edge in self.acc_Graph.in_edges(node)]
        nodes_accumulation = [self.acc_slice[node[0], node[1]] for node in in_nodes]
        return in_nodes[nodes_accumulation.index(max(nodes_accumulation))]

    def in_nodes(self, node):
        return [node[0] for node in self.G.in_edges(node)]


    def accumulation_for_point(self, point):
        return self.acc[point[0], point[1]]



    
    # Flood Part

    def compute_flood(self, coordinate, top_left, bottom_right, lenth, target_h):
        point = self.coordinate2point(coordinate, top_left, bottom_right)
        y, x = point[0], point[1]
        print(point)
        lenth_with_offset = int(lenth/35 + 40)  # For offset
        flood_area_fdir = self.fdir[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        flood_area_dem = self.dem[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        flood_area_acc = self.acc.copy()
        flood_area_acc = flood_area_acc[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        new_x, new_y = lenth_with_offset, lenth_with_offset
        new_point = (new_y, new_x)

        # Create temporary acc graph
        self.create_acc_graph(flood_area_acc, flood_area_fdir, 200)

        # Get river pathes
        river_pathes_nodes_DOWN, river_pathes_nodes_UP = self.compute_river_path(new_point, lenth)

        # Graph for specified area
        self.G = nx.DiGraph()
        shape = flood_area_fdir.shape

        for row in tqdm(range(1, shape[0]-1)):
            for column in range(1, shape[1]-1):
                dir = flood_area_fdir[row, column]
                start = (row, column)
                target = self.fdir_coordinate(start, dir)
                self.G.add_edge(start, target)

        # Make flood
        h = flood_area_dem[new_y, new_x] + target_h

        flooded_nodes_down = []
        all_out_nodes = [] # For both up and down

        for i, node in tqdm(enumerate(river_pathes_nodes_DOWN[::-1])): # Начинаем с последней затопленной клетки
            all_nodes = [node]
            nodes = [node]
            out_nodes_log = []

            while len(nodes) > 0:
                node_ = nodes[0]
                if node_ == river_pathes_nodes_DOWN[0]:
                    nodes.pop(0)
                    break

                nodes.pop(0)
                in_nodes = self.in_nodes(node_)
                #print(node, in_nodes)
                if len(in_nodes) == 0:
                    break

                intersection = set(river_pathes_nodes_DOWN).intersection(set(in_nodes))
                if len(intersection) > 0:
                    in_nodes.remove(list(intersection)[0]) # Удаление участков реки ниже. Чтобы обрабатывать только области у рек.
                    

                in_nodes_ = [node for node in in_nodes if flood_area_dem[node[0], node[1]] <= h]
                    
                all_nodes += in_nodes_
                nodes.append(in_nodes_)

                # adding in-edge parts
                out_nodes = [node for node in in_nodes if flood_area_dem[node[0], node[1]] > h]
                
                for out_node in out_nodes:
                    start, end = out_node, self.out_node_G(out_node)
                    start_height, end_height = flood_area_dem[start[0], start[1]], flood_area_dem[end[0], end[1]]
                    height_rise = abs(end_height - h)
                    height_difference = abs(start_height - end_height)
                    meter_path = height_difference / 35
                    point_path = round(height_rise / meter_path, 0) # * out of 35 (in meteres). From end (lower) to upper.
                    out_nodes_log.append((out_node, end, point_path))
            
            if len(all_nodes) == 0:
                continue
            flooded_nodes_down += all_nodes
            all_out_nodes += out_nodes_log

        flooded_nodes_up = []
        for i, node in tqdm(enumerate(river_pathes_nodes_UP)): # Начинаем с первой клетки
            all_nodes = [node]
            nodes = [node]
            out_nodes_log = []

            while len(nodes) > 0:
                node_ = nodes[0]
                if node_ == river_pathes_nodes_UP[-1]:
                    nodes.pop(0)
                    break

                nodes.pop(0)
                in_nodes = self.in_nodes(node_)
                if len(in_nodes) == 0:
                    break

                intersection = set(river_pathes_nodes_UP).intersection(set(in_nodes))
                if len(intersection) > 0:
                    in_nodes.remove(list(intersection)[0]) # Удаление участков реки ниже. Чтобы обрабатывать только области у рек.
                    
                in_nodes_ = [node for node in in_nodes if flood_area_dem[node[0], node[1]] <= h]
                    
                all_nodes += in_nodes_
                nodes.append(in_nodes_)

                # adding in-edge parts
                out_nodes = [node for node in in_nodes if flood_area_dem[node[0], node[1]] > h]
                
                for out_node in out_nodes:
                    start, end = out_node, self.out_node_G(out_node)
                    start_height, end_height = flood_area_dem[start[0], start[1]], flood_area_dem[end[0], end[1]]
                    height_rise = abs(end_height - h)
                    height_difference = abs(start_height - end_height)
                    meter_path = height_difference / 35
                    point_path = round(height_rise / meter_path, 0) # * out of 35 (in meteres). From end (lower) to upper.
                    out_nodes_log.append((out_node, end, point_path))

            
            if len(all_nodes) == 0:
                continue
            flooded_nodes_up += all_nodes
            all_out_nodes += out_nodes_log

        return flood_area_acc.shape, flooded_nodes_down, flooded_nodes_up, all_out_nodes, flood_area_dem, h, point, lenth_with_offset
