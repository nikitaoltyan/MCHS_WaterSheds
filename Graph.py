import WaterShed

import pandas as pd
import numpy as np
from numpy import unravel_index
import pytz
from datetime import datetime
import networkx as nx
from tqdm import tqdm
from os import path
import os



class Graph():

    def __init__(self, dem=None, fdir=None, acc=None, compression=1):
        self.dem = dem
        self.fdir = fdir
        self.acc = acc
        self.compression = compression
        self.point_size_meteres = 35


    def compute_height(self, file_path, DEMs_path, min_acc, save_path):
        df = pd.read_csv(file_path, sep=';', decimal=',')
        df['x_lon_int'] = df['x_lon'].values.astype(int)
        df['y_lat_int'] = df['y_lat'].values.astype(int)
        
        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis = 0, ascending = True, inplace = True, na_position = "first")
        
        dt_string = datetime.now(pytz.timezone('Europe/Moscow')).strftime("%d_%m_%Y__%H:%M")
        self.df_new = pd.DataFrame(columns=['hstation_id', 'x_lon', 'y_lat', 'height', 'distance_m', 'error'])
        self.df_new = self.df_new.astype(dtype= {'hstation_id':'int64', 'height':'int64', 'distance_m':'int64', 'error':'int64'})
        x_lon_past, y_lat_past = None, None

        for i, row in df.iterrows():
            print(f'{i+1}/{df.shape[0]} hydropost...')
            hstation_id  = int(row[0])
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
                
                if lat_num+1 < 60:
                    self.point_size_meteres = 35
                else:
                    self.point_size_meteres = 65

                self.tif_pathes = []
                for i in range(lat_num-1, lat_num+2):
                    for j in range(lng_num-1, lng_num+2):
                        lat = str(i)
                        lng = ''.join((['0'] + list(str(int(j))))[-3:])
                        file_name = f'n{lat}_e{lng}_1arc_v3.tif' if lat_num+1 < 60 else f'n{lat}_e{lng}_1arc_v3_1201x1201.tif'
                        self.tif_pathes.append(f'{DEMs_path}/{file_name}')
                        
                # check if files 'exisits'
                success_list = []
                for tif_path in self.tif_pathes:
                    if path.exists(tif_path) == False:
                        print(f'{tif_path} is not exist in path {DEMs_path}')
                        success_list.append(False)

                 # Download DEM and preprocess it
                if len(success_list) == 0:
                    print('All required DEMs exist')
                    self.compute_DEM(self.tif_pathes, lng_num, lat_num)
                else:
                    # Temporary while I'm thinking what to with others frames DEMs
                    print('Not all required DEMs exist')
                    self.compression = 1
                    self.acc = None
                    self.dem = None
                    self.fdir = None
            

            # Calculate Heights
            top_left = (lng_num-1, lat_num+2) if len(self.tif_pathes) == 9 else (lng_num, lat_num+1)
            bottom_right = (lng_num+2, lat_num-1) if len(self.tif_pathes) == 9 else (lng_num+1, lat_num)

            if self.dem is not None:
                height, distance, success = self.compute_height_differance(coordinate, top_left, bottom_right, 10000, min_acc)
                error = 1 if success == False else 0

                dct = {
                    'hstation_id': hstation_id, 
                    'x_lon': x_lon, 
                    'y_lat': y_lat, 
                    'height': int(height),
                    'distance_m': int(distance),
                    'error': int(error)
                }
                self.df_new = self.df_new.append(dct, ignore_index=True)
                self.df_new.to_csv(f'{save_path}/{dt_string}_hydroposts_height_calculated.csv', sep=';', decimal=',', index=False)
            else:
                # here is save for error hydropost (corner hydropost without all DEMs)
                dct = {
                    'hstation_id': hstation_id, 
                    'x_lon': x_lon, 
                    'y_lat': y_lat, 
                    'height': 0,
                    'distance_m': 0,
                    'error': 1
                }
                self.df_new = self.df_new.append(dct, ignore_index=True)
                self.df_new.to_csv(f'{save_path}/{dt_string}_hydroposts_height_calculated.csv', sep=';', decimal=',', index=False)

    

    def compute_DEM(self, pathes, lng_num, lat_num, compression=3):
        if len(pathes) == 9:
            shed = WaterShed.WaterSheds(files_pathes=pathes, compute_acc=True, compression=compression)
            self.compression = compression
            self.acc = shed.acc
            self.dem = shed.dem
            self.fdir = shed.fdir
        else:
            #  Use only central tile
            lat = str(lat_num)
            lng = ''.join((['0'] + list(str(int(lng_num))))[-3:])
            shed = WaterShed.WaterSheds(file_path=f'n{lat}_e{lng}_1arc_v3.tif', compute_acc=True)
            self.compression = compression
            self.acc = shed.acc
            self.dem = shed.dem
            self.fdir = shed.fdir


    def compute_height_differance(self, coordinate, top_left, bottom_right, lenth, min_acc):
        # In case to not create acc_graph every time for same lon & lat
        if (self.acc_slice is None) or (self.acc_Graph is None):
            acc_slice = self.acc.copy()
            # Filter river cells
            self.create_acc_graph(acc_slice, self.fdir, min_acc)

        point = self.coordinate2point(coordinate, top_left, bottom_right)
        river_pathes_nodes_DOWN, river_pathes_nodes_UP, distance, success = self.compute_river_path(point, lenth)

        start, end = river_pathes_nodes_DOWN[-1], river_pathes_nodes_UP[-1]
        start_height, end_height = self.dem[start], self.dem[end]
        return abs(start_height - end_height), distance, success


    def compute_river_path(self, point, lenth):
        """
        Function returns all river nodes down and up from the given point.
        * lenth - in meteres from start to up and down.
        """
        point_lenth = self.compression * self.point_size_meteres       # around 35 meteres shape for each point

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
                
        distance = (len(river_pathes_nodes_DOWN) + len(river_pathes_nodes_UP)) * point_lenth
        success = True if (len(river_pathes_nodes_DOWN) == river_pathes_lenght_DOWN + 1) and (len(river_pathes_nodes_UP) > (0.5*river_pathes_lenght_UP)) else False
        return river_pathes_nodes_DOWN, river_pathes_nodes_UP, distance, success


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
        for i in range(1, acc_slice.shape[0]-1):
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
        else: #Up-Left
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

    def compute_flood(self, coordinate, top_left, bottom_right, lenth, target_h, uniform_flooding=False):
        point = self.coordinate2point(coordinate, top_left, bottom_right)
        y, x = point[0], point[1]
        lenth_with_offset = int(lenth/(self.point_size_meteres * self.compression) + 40)  # For offset
        flood_area_fdir = self.fdir[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        flood_area_dem = self.dem[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        flood_area_acc = self.acc.copy()
        flood_area_acc = flood_area_acc[y-lenth_with_offset:y+lenth_with_offset, x-lenth_with_offset:x+lenth_with_offset]
        new_x, new_y = lenth_with_offset, lenth_with_offset
        new_point = (new_y, new_x)

        # Create temporary acc graph
        self.create_acc_graph(flood_area_acc, flood_area_fdir, 200)

        # Get river pathes
        river_pathes_nodes_DOWN, river_pathes_nodes_UP, _, success = self.compute_river_path(new_point, lenth)

        # Graph for specified area
        self.G = nx.DiGraph()
        shape = flood_area_fdir.shape

        for row in range(1, shape[0]-1):
            for column in range(1, shape[1]-1):
                dir = flood_area_fdir[row, column]
                start = (row, column)
                target = self.fdir_coordinate(start, dir)
                self.G.add_edge(start, target)

        # Make flood
        if uniform_flooding == False:
            self.h = flood_area_dem[new_y, new_x] + target_h

        flooded_nodes_down = []
        all_out_nodes = [] # For both up and down
        for i, node in enumerate(river_pathes_nodes_DOWN[::-1]): # Начинаем с последней затопленной клетки
            all_nodes = [node]
            nodes = [node]
            out_nodes_log = []
            if uniform_flooding:
                self.h = flood_area_dem[node[0], node[1]] + target_h

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
                    

                in_nodes_ = [node for node in in_nodes if flood_area_dem[node[0], node[1]] <= self.h]
                    
                all_nodes += in_nodes_
                nodes.append(in_nodes_)

                # adding in-edge parts
                out_nodes = [node for node in in_nodes if flood_area_dem[node[0], node[1]] > self.h]
                
                for out_node in out_nodes:
                    start, end = out_node, self.out_node_G(out_node)
                    start_height, end_height = flood_area_dem[start[0], start[1]], flood_area_dem[end[0], end[1]]
                    height_rise = abs(end_height - self.h)
                    height_difference = abs(start_height - end_height)
                    meter_path = height_difference / self.point_size_meteres
                    point_path = round(height_rise / meter_path, 0) # * out of self.point_size_meteres (in meteres). From end (lower) to upper.
                    out_nodes_log.append((out_node, end, point_path))
            
            if len(all_nodes) == 0:
                continue
            flooded_nodes_down += all_nodes
            all_out_nodes += out_nodes_log

        flooded_nodes_up = []
        for i, node in enumerate(river_pathes_nodes_UP): # Начинаем с первой клетки
            all_nodes = [node]
            nodes = [node]
            out_nodes_log = []
            if uniform_flooding:
                self.h = flood_area_dem[node[0], node[1]] + target_h

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
                    
                in_nodes_ = [node for node in in_nodes if flood_area_dem[node[0], node[1]] <= self.h]
                    
                all_nodes += in_nodes_
                nodes.append(in_nodes_)

                # adding in-edge parts
                out_nodes = [node for node in in_nodes if flood_area_dem[node[0], node[1]] > self.h]
                
                for out_node in out_nodes:
                    start, end = out_node, self.out_node_G(out_node)
                    start_height, end_height = flood_area_dem[start[0], start[1]], flood_area_dem[end[0], end[1]]
                    height_rise = abs(end_height - self.h)
                    height_difference = abs(start_height - end_height)
                    meter_path = height_difference / self.point_size_meteres
                    point_path = round(height_rise / meter_path, 0) # * out of self.point_size_meteres (in meteres). From end (lower) to upper.
                    out_nodes_log.append((out_node, end, point_path))

            
            if len(all_nodes) == 0:
                continue
            flooded_nodes_up += all_nodes
            all_out_nodes += out_nodes_log

        self.h = flood_area_dem[new_y, new_x] + target_h
        return flood_area_acc.shape, flooded_nodes_down, flooded_nodes_up, all_out_nodes, flood_area_dem, self.h, point, lenth_with_offset

    
    def get_step(self, y_delta, x_delta):
        """
        Calculate step for future river slice
        """
        if ((y_delta > 3) and (x_delta > 3)) \
            or ((y_delta < -3) and (x_delta < -3)):
            return (1, 1)
        elif ((y_delta > 3) and (x_delta < -3)) \
            or ((y_delta < -3) and (x_delta >3)):
            return (1, -1)
        elif x_delta <= 3:
            return (0, 1)
        else:
            return (1, 0)

        
    def get_river_slice(self, file_path, DEMs_path, save_path):
        df = pd.read_csv(file_path, sep=';', decimal=',')
        df['x_lon_int'] = df['x_lon'].values.astype(int)
        df['y_lat_int'] = df['y_lat'].values.astype(int)
        
        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis = 0, ascending = True, inplace = True, na_position = "first")
        
        # Creat new df to save successes of river slices
        self.df_new = pd.DataFrame(columns=['hstation_id', 'success'])
        
        x_lon_past, y_lat_past = None, None
        
        for i, row in df.iterrows():
            print(f'{i+1}/{df.shape[0]} hydropost...')
            hstation_id  = int(row[0])
            x_lon, y_lat = row[1], row[2]
            coordinate = (x_lon, y_lat)
        
            # Define coordinate of map to download 
            lng_num, lat_num = int(x_lon), int(y_lat)

            # Check if this coordinates weren't calculated
            if (x_lon_past != lng_num) or (y_lat_past != lat_num):
                x_lon_past, y_lat_past = lng_num, lat_num

                self.tif_pathes = []
                for i in range(lat_num-1, lat_num+2):
                    for j in range(lng_num-1, lng_num+2):
                        lat = str(i)
                        lng = ''.join((['0'] + list(str(int(j))))[-3:])
                        file_name = f'n{lat}_e{lng}_1arc_v3.tif' if lat_num+1 < 60 else f'n{lat}_e{lng}_1arc_v3_1201x1201.tif'
                        self.tif_pathes.append(f'{DEMs_path}/{file_name}')
                        
                # check if files 'exisits'
                success_list = []
                for tif_path in self.tif_pathes:
                    if path.exists(tif_path) == False:
                        print(f'{tif_path} is not exist in path {DEMs_path}')
                        success_list.append(False)

                # Download DEM and preprocess it
                if len(success_list) == 0:
                    print('All required DEMs exist')
                    self.compute_DEM(self.tif_pathes, lng_num, lat_num, compression=1)
                else:
                    # Temporary while I'm thinking what to with others frames DEMs
                    print('Not all required DEMs exist')
                    self.compression = 1
                    self.acc = None
                    self.dem = None
                    self.fdir = None
        
        
            # Calculate Heights
            top_left = (lng_num-1, lat_num+2) if len(self.tif_pathes) == 9 else (lng_num, lat_num+1)
            bottom_right = (lng_num+2, lat_num-1) if len(self.tif_pathes) == 9 else (lng_num+1, lat_num)

            if self.dem is not None:
                point = self.coordinate2point(coordinate, top_left, bottom_right)
            
                new_top_left, new_bottom_right = (point[0]-70, point[1]-70), (point[0]+70, point[1]+70)
                cut_fdir = self.fdir[new_top_left[0]:new_bottom_right[0], new_top_left[1]:new_bottom_right[1]]
                cut_acc = self.acc[new_top_left[0]:new_bottom_right[0], new_top_left[1]:new_bottom_right[1]]
                cut_dem = self.dem[new_top_left[0]:new_bottom_right[0], new_top_left[1]:new_bottom_right[1]]

                # Graph for specified area
                G = nx.DiGraph()
                shape = cut_fdir.shape

                try:
                    for row in tqdm(range(1, shape[0]-1)):
                        for column in range(1, shape[1]-1):
                            dir = cut_fdir[row, column]
                            start = (row, column)
                            target = self.fdir_coordinate(start, dir)
                            G.add_edge(start, target)
                except:
                    dct = {
                        'hstation_id': hstation_id, 
                        'success': 0
                    }
                    self.df_new = self.df_new.append(dct, ignore_index=True)
                    self.df_new.to_csv(f'{save_path}/river_slice_success_table.csv', sep=';', decimal=',', index=False)
                    continue

                real_coord = (point[0] - new_top_left[0], point[1] - new_top_left[1])
                
                def out_node_G(node):
                    return [node for node in G.edges(node)][0][1]
    
                def in_node_G(node):
                    in_nodes = [edge[0] for edge in G.in_edges(node)]
                    nodes_accumulation = [cut_acc[node[0], node[1]] for node in in_nodes]
                    return in_nodes[nodes_accumulation.index(max(nodes_accumulation))]

                outs = [real_coord]
                while len(outs) < 5:
                    out = out_node_G(outs[-1])
                    outs.append(out)

                ins = [real_coord]
                while len(ins) < 5:
                    in_ = in_node_G(ins[-1])
                    ins.append(in_)

                start, end = ins[-1], outs[-1]
                y_delta, x_delta = end[0] - start[0], end[1] - start[1]

                step = self.get_step(y_delta, x_delta)
                target_height = cut_dem[real_coord] + 15
                start = real_coord

                right_heights = [cut_dem[real_coord]]
                right_coords = [real_coord]
                while (max(right_heights) < target_height) and (len(right_heights) < 65):
                    start = [sum(x) for x in zip(start, step)]
                    height = cut_dem[start[0], start[1]]
                    right_heights.append(height)
                    right_coords.append(start)

                start = real_coord
                step = [-i for i in step]
                left_heights = [cut_dem[real_coord]]
                left_coords = [real_coord]
                while (max(left_heights) < target_height) and (len(left_heights) < 65):
                    start = [sum(x) for x in zip(start, step)]
                    height = cut_dem[start[0], start[1]]
                    left_heights.append(height)
                    left_coords.append(start)

                river_slice = left_heights[::-1] + right_heights[1:]
                coords_slice = left_coords[::-1] + right_coords[1:]
                coords_bin_slice = [0 if type(coord) != tuple else 1 for coord in coords_slice]

                if os.path.exists(f'{save_path}/{hstation_id}') == False:
                    os.makedirs(f'{save_path}/{hstation_id}')

                path_to_csv = f'{save_path}/{hstation_id}/{hstation_id}_river_slice.csv'
                final_df = pd.DataFrame({'HEIGHTS': river_slice, 'WaterpostFlag': coords_bin_slice})
                final_df.to_csv(path_to_csv, index=False, sep=';')
                
                dct = {
                    'hstation_id': hstation_id, 
                    'success': 1
                }
                self.df_new = self.df_new.append(dct, ignore_index=True)
                self.df_new.to_csv(f'{save_path}/river_slice_success_table.csv', sep=';', decimal=',', index=False)
                
            else:
                # here is save for error hydropost (corner hydropost without all DEMs)
                dct = {
                    'hstation_id': hstation_id, 
                    'success': 0
                }
                self.df_new = self.df_new.append(dct, ignore_index=True)
                self.df_new.to_csv(f'{save_path}/river_slice_success_table.csv', sep=';', decimal=',', index=False)
