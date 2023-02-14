import sys
import subprocess
import os
from os import path

import numpy as np
from numpy import unravel_index
import pandas as pd
from PIL import Image
import struct
from tqdm import tqdm
import random

import skimage.morphology
import skimage.measure

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Circle
import seaborn as sns
from scipy.stats import mode
from scipy.ndimage import convolve

from osgeo import gdal, osr, ogr # Python bindings for GDAL

from glob import glob

import networkx as nx

# Try to import pysheds
try:
    from pysheds.grid import Grid
    print('pysheds imported')
except:
    print('intalling pysheds')
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'pysheds'])
    from pysheds.grid import Grid
    print('pysheds imported')
    
    
# Upgrade tbb
print('Upgrading tbb')
subprocess.check_call([sys.executable, "-m", "pip", "install", '--upgrade', 'tbb'])

# Upgrade numba
print('Upgrading numba')
subprocess.check_call([sys.executable, "-m", "pip", "install", '--upgrade', 'numba'])
    

class WaterSheds():
    
    def __init__(self, 
                 file_path=None, 
                 files_pathes=None, 
                 compute_acc=False,
                 compression=1):
        assert file_path != None or files_pathes != None, "Insert file path!"
        self.dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
        self.box = ((104, 105), (51, 52)) # Should be difined automatically.
        self.compression_coef = compression

        if files_pathes == None:
            self.file_path = file_path
        else:
            self.file_path = "merged.tif"
            self.__merge_tifs(files_pathes)
            

        self.__grid()
        self.__dem()
        if compute_acc:
            self.compute_flow_directions()
            self.compute_accumulation()
        print('All data prepared')

    def __merge_tifs(self, pathes):
        """ Merging multiple geo-files together.
        """
        # Delete merged tif file if exists
        try:
            os.remove(f'/content/{self.file_path}')
        except:
            pass

        print('Merging...')
        cmd = f"gdal_merge.py -o {self.file_path}"
        subprocess.call(cmd.split()+pathes)
        print('Merged')

        # Compress if needed:
        if self.compression_coef != 1:
            self.__compress()


    def __grid(self):
        self.grid = Grid.from_raster(self.file_path)

    def __dem(self):
        dem = self.grid.read_raster(self.file_path)
        self.dem = dem

    def __compress(self):
        # Compressing .tif file for file_name
        print("Compressing...")

        # Open tif
        dataset = gdal.Open(self.file_path, gdal.GA_ReadOnly)
        projection = dataset.GetProjection()
        rast_band = dataset.GetRasterBand(1)
        img_array = rast_band.ReadAsArray()

        shape = img_array.shape
        coeff = self.compression_coef
        print(f'Shape before compressing: {shape}')

        # Collecting rows & columns to delete
        rows_to_delete = [i for i in range(shape[0]) if i % coeff != 0]
        columns_to_delete = [i for i in range(shape[1]) if i % coeff != 0]

        # Deleting rows & columns
        img_array = np.delete(img_array, columns_to_delete, 1) # delete column
        img_array = np.delete(img_array, rows_to_delete, 0) # delete rows
        
        print(f'Shape after compressing: {img_array.shape}')

        # Save merged tif
        driver = gdal.GetDriverByName("GTiff")
        driver.Register()
        out_ds = driver.Create(self.file_path, img_array.shape[1], img_array.shape[0], 2, gdal.GDT_Float32)
        out_ds.SetProjection(projection)
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(img_array)
        out_band.SetNoDataValue(np.nan)
        out_band.FlushCache()
        out_ds=None
        print('Compressed')


    def __resolve_ground(self):
        pit_filled_dem = self.grid.fill_pits(self.dem)
        print("1/4 computed!")
        flooded_dem = self.grid.fill_depressions(pit_filled_dem)
        print("2/4 computed!")
        inflated_dem = self.grid.resolve_flats(flooded_dem)
        print("3/4 computed!")
        return inflated_dem

    def compute_flow_directions(self):
        print('Computing flow directions...')
        inflated_dem = self.__resolve_ground()
        self.fdir = self.grid.flowdir(inflated_dem, dirmap=self.dirmap)
        print("4/4 computed!")
        print('Done.')
        return self.fdir

    def compute_accumulation(self, do_previous_steps=False):
        print('Computing water accumulation...')
        if do_previous_steps:
            self.compute_flow_directions()
        self.acc = self.grid.accumulation(self.fdir, dirmap=self.dirmap)
        print('Done.')
        return self.acc

    def accumulation_for_point(self, point):
        return self.acc[point[0], point[1]]

    def extract_river_network(self, min_acc=5000, do_previous_steps=False):
        if do_previous_steps:
            self.compute_accumulation(do_previous_steps=True)
        self.branches = self.grid.extract_river_network(
            self.fdir, 
            self.acc > min_acc, 
            dirmap=self.dirmap
        )
        return self.branches

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

    def create_DiGraph(self, do_previous_steps=False):
        # TODO: Add simplified version (with less resolution).
        if do_previous_steps:
            self.compute_flow_directions()
        G = nx.DiGraph()
        shape = self.fdir.shape
        for row in tqdm(range(1, shape[0]-1)):
            for column in range(1, shape[1]-1):
                dir = self.fdir[row, column]
                start = (row, column)
                target = self.fdir_coordinate(start, dir)
                G.add_edge(start, target)
        self.Graph = G
        print("Graph created!")
        return G

    def __in_nodes(self, edges):
        return [edge[0] for edge in edges]

    def __get_pixel_position(self, M):
        #TODO: Compute x_path, y_path dinamicaly.
        _x_path, _y_path = 1/1800, 1/3601
        x, y = M[:,0] % 1, M[:,1] % 1
        return np.array((3601 - np.fix(y/_y_path).astype(int), np.fix(x/_x_path).astype(int))).T[0]

    def __define_maximum_near(self, M):
        max_coord = M
        row, column = M
        pathes = [
              (row-1, column),
              (row-1, column+1),
              (row, column+1),
              (row+1, column+1),
              (row+1, column),
              (row+1, column-1),
              (row, column-1),
              (row-1, column-1)
        ]

        # TODO: Set 4 anchor points and get coordinate of the maximum point.
        # That's it!
        for path in pathes:
            try:
                if self.acc[path[0], path[1]] > self.acc[max_coord[0], max_coord[1]]:
                    max_coord = path
            except:
                pass
        return max_coord

    def extract_watersheds_positions(self, do_previous_steps=False):
        if do_previous_steps:
            self.extract_river_network(do_previous_steps=True)
        position = []
        box = self.box
        for branch in self.branches['features']:
            line = np.asarray(branch['geometry']['coordinates'])
            start, end = line[0], line[-1]
            if (end[0] > box[0][0]) and (end[0] < box[0][1]) and (end[1] > box[1][0]) and (end[1] < box[1][1]):
                start_pos, end_pos = tuple(self.__get_pixel_position(np.array([start]))), tuple(self.__get_pixel_position(np.array([end])))
                start_pos = self.__define_maximum_near(start_pos)
                end_pos = self.__define_maximum_near(end_pos)
                position.append([start_pos, end_pos])
        return position
    
    def __watershed(self, start, break_point, do_previous_steps=False):
        all_nodes = []
        nodes = [break_point]
        break_point = start

        while len(nodes) > 0:
            for node in nodes:
                nodes.pop(0)
                in_nodes = self.__in_nodes(self.Graph.in_edges(node))
                if len(in_nodes) == 0:
                    break

                try:
                    in_nodes.remove(break_point)
                except:
                    pass
                
                all_nodes += in_nodes
                nodes.append(in_nodes)

        return all_nodes

    def extract_watersheeds(self, do_previous_steps=False):
        if do_previous_steps:
            self.create_DiGraph(do_previous_steps=True)
        sheds = []
        positions = self.extract_watersheds_positions(do_previous_steps=do_previous_steps)
        for pos in tqdm(positions):
            print(f"Position: {pos}")
            print(f"Accumulation: {self.acc[pos[0], pos[1]]}")
            sheds.append(self.__watershed(pos[0], pos[1]))
        self.sheds = sheds
        return sheds

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


    def compute_watershed_area(self, coordinate, top_left, bottom_right):
        point_index = self.coordinate2point(coordinate, top_left, bottom_right)
        return self.watershed_area(point_index)


    def watershed_area(self, point):
        """
        Calculating watershed area for specific point.
        #Area of 1 point is defined in merging (based on compression).

        Returns: 
            * area of watershed in point (sq.km)
        """
        x_size_point = 30.8 * self.compression_coef
        y_size_point = 37.3 * self.compression_coef
        area_per_point_km = (x_size_point * y_size_point) / 1000000
        accumulated = self.acc[point[0], point[1]]
        return area_per_point_km * accumulated



    # Drawing part
    def draw_dem(self, figsize=(11,7)):
        fig, ax = plt.subplots(figsize=figsize)
        fig.patch.set_alpha(0)

        plt.imshow(self.dem, cmap='terrain', zorder=1)
        plt.colorbar(label='Elevation (m)')

        plt.title('Digital elevation map', size=10)
        plt.tight_layout()

    def draw_fdir(self, figsize=(11,7)):
        fig = plt.figure(figsize=figsize)
        fig.patch.set_alpha(0)

        plt.imshow(self.fdir, cmap='viridis', zorder=2)
        boundaries = ([0] + sorted(list(self.dirmap)))
        plt.colorbar(boundaries= boundaries, values=sorted(self.dirmap))
        plt.title('Flow direction grid', size=10)
        plt.grid(zorder=-1)
        plt.tight_layout()

    def draw_acc(self, figsize=(11,7)):
        fig, ax = plt.subplots(figsize=figsize)
        fig.patch.set_alpha(0)
        plt.grid('on', zorder=0)

        im = ax.imshow(self.acc, zorder=2, cmap='cubehelix',
                      norm=colors.LogNorm(1, self.acc.max()), interpolation='bilinear')

            
        plt.colorbar(im, ax=ax, label='Upstream Cells')
        plt.title('Flow Accumulation', size=10)
        plt.tight_layout()

    def draw_branches(self, figsize=(11,11)):
        sns.set_palette('husl')
        fig, ax = plt.subplots(figsize=figsize)

        plt.xlim(self.grid.bbox[0], self.grid.bbox[2])
        plt.ylim(self.grid.bbox[1], self.grid.bbox[3])

        for branch in self.branches['features']:
            line = np.asarray(branch['geometry']['coordinates'])
            plt.plot(line[:, 0], line[:, 1])
            
        _ = plt.title('Rivers network', size=10)
