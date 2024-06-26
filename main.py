import WaterShed
import Graph
import guard

from osgeo import gdal, osr, ogr # Python bindings for GDAL
import pytz
from datetime import datetime
from tqdm import tqdm
import numpy as np
import pandas as pd
from os import path
import os
import sys



class Main():
  
    def __init__(self):
        self.shed = None
        self.DEMs_path = None
        self.lng_num = None
        self.lat_num = None
        
        
    def _compute_shed(self, x_lon, y_lat, compression):
        """Compute shed if it isn't exists
        """
        lng_num, lat_num = int(x_lon), int(y_lat)
        
        if (self.lng_num != lng_num) or (self.lat_num != lat_num):
            self.lng_num = lng_num
            self.lat_num = lat_num
          
            self.tif_pathes = []
            for i in range(lat_num-1, lat_num+2):
                for j in range(lng_num-1, lng_num+2):
                    lat = str(i)
                    lng = ''.join((['0'] + list(str(int(j))))[-3:])
                    file_name = f'n{lat}_e{lng}_1arc_v3.tif' if lat_num+1 < 60 else f'n{lat}_e{lng}_1arc_v3_1201x1201.tif'
                    self.tif_pathes.append(f'{self.DEMs_path}/{file_name}')
                    
            # check if files 'exists'
            success_list = []
            for tif_path in self.tif_pathes:
                if path.exists(tif_path) == False:
                    print(f'{tif_path} is not exist in path {self.DEMs_path}')
                    success_list.append(False)

            # Download DEM and preprocess it
            if len(success_list) == 0:
                print('All required DEMs exist')
                self.shed = WaterShed.WaterSheds(files_pathes=self.tif_pathes, compute_acc=True, compression=compression)
                return self.shed
            else:
                return None
        else:
            return self.shed
          
          
    def __frequency_to_name(self, frequency):
        if frequency == 0.5:
            return '005'
        elif frequency == 1.0:
            return '01'
        elif frequency == 5.0:
            return '05'
        elif frequency == 10.0:
            return '10'
        else:
            return '20'
  
  
    def compute_shape(self, save_path, data_dict, uniform_flooding=False):
        coordinate = data_dict['coordinate']
        freq_name = data_dict['frequency_name']
        frequency = data_dict['frequency']
        waterpost_name = int(data_dict['hsts_id'])
        wtrdepth = data_dict['wtrdepth']
        wtrlvltime = data_dict['wtrlvltime']
        compression = data_dict['compression']
        lenth = data_dict['lenth']
        target_h = data_dict['wtrdepth']
        
        # Define coordinate of map to download 
        (x_lon, y_lat) = coordinate
      
        # Get WaterShed (precomputed or computed on demand)
        shed = self._compute_shed(x_lon, y_lat, compression)
        
        lng_num, lat_num = int(x_lon), int(y_lat)
        if shed is not None:
            top_left = (lng_num-1, lat_num+2)
            bottom_right = (lng_num+2, lat_num-1)

            GraphClass = Graph.Graph(dem=shed.dem, fdir=shed.fdir, acc=shed.acc, compression=shed.compression_coef)
            (shape, 
                flooded_nodes_down, 
                flooded_nodes_up,
                all_out_nodes,
                dem, 
                height, 
                point_coord, 
                offset
            ) = GraphClass.compute_flood(coordinate, top_left, bottom_right, lenth, target_h, uniform_flooding=uniform_flooding)

            new_space_no_interpol = np.zeros((shape[0], shape[1]), dtype=np.uint8)
            for down_node in flooded_nodes_down:
                y, x = down_node[0], down_node[1]
                new_space_no_interpol[y, x] = 1

            for up_node in flooded_nodes_up:
                y, x = up_node[0], up_node[1]
                new_space_no_interpol[y, x] = 1

            flooded_area = (len(flooded_nodes_down) + len(flooded_nodes_up)) * ((0.03*compression)**2)
            flood_outer = new_space_no_interpol

            main_point = point_coord
            top_left, bottom_right = top_left, bottom_right
            main_shape = shed.dem.shape

            y_path, x_path = (top_left[1]-bottom_right[1])/main_shape[0], (bottom_right[0] - top_left[0])/main_shape[1]

            cut_top_left = (top_left[1]-(main_point[0]-offset)*y_path, top_left[0]+(main_point[1]-offset)*x_path)
            cut_bottom_right = (top_left[1]-(main_point[0]+offset)*y_path, top_left[0]+(main_point[1]+offset)*x_path)

            y_path, x_path = (cut_top_left[0]-cut_bottom_right[0])/flood_outer.shape[0], (cut_bottom_right[1] - cut_top_left[1])/flood_outer.shape[1]

            # ----------------------------
            print('Flood computed!')
            print('Saving it as .tif...')
            # ----------------------------



            nlines = flood_outer.shape[0]
            ncolumns = flood_outer.shape[1]
            data = flood_outer

            def getGeoTransform(extent, nlines, ncols):
                resx = (extent[2] - extent[0]) / ncols
                resy = (extent[3] - extent[1]) / nlines
                return [extent[0], resx, 0, extent[3] , 0, -resy]

            # Define the data extent (min. lon, min. lat, max. lon, max. lat)
            extent = [cut_top_left[1], cut_bottom_right[0], cut_bottom_right[1], cut_top_left[0]]

            # Export the test array to GeoTIFF ================================================

            # Get GDAL driver GeoTiff
            driver = gdal.GetDriverByName('GTiff')

            # Get dimensions
            nlines = data.shape[0]
            ncols = data.shape[1]
            nbands = len(data.shape)
            data_type = gdal.GDT_Int16 # gdal.GDT_Float32

            # Create a temp grid
            grid_data = driver.Create('grid_data', ncols, nlines, 1, data_type)

            # Write data for each bands
            grid_data.GetRasterBand(1).WriteArray(data)

            # Lat/Lon WSG84 Spatial Reference System
            srs = osr.SpatialReference()
            srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

            # Setup projection and geo-transform
            grid_data.SetProjection(srs.ExportToWkt())
            grid_data.SetGeoTransform(getGeoTransform(extent, nlines, ncols))

            # Save the file
            file_name = f'{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}.tif'
            try:
                os.makedirs(f"{save_path}/{waterpost_name}")
                driver.CreateCopy(file_name, grid_data, 0)
                print(f'Generated GeoTIFF: {file_name}')
            except:
                driver.CreateCopy(file_name, grid_data, 0)


            # Close the file
            driver = None
            grid_data = None

            # Delete the temp grid               
            os.remove('grid_data')
            
            dct = {
                'hstst_id': int(waterpost_name), 
                'lat_y': y_lat,
                'lon_x': x_lon,
                'frequency': frequency, 
                'wtrdepth': wtrdepth, 
                'wtrlvltime': wtrlvltime,
                'lenth': lenth,
                'flooded_area': flooded_area,
                'success': 1
            }
            self.df_new = self.df_new.append(dct, ignore_index=True)
            self.df_new.to_csv(f'{save_path}/{self.dt_string}_tifs_for_shape_calculated.csv', sep=';', decimal=',', index=False)

            # ----------------------------
            print(f'Tif for hstation_id: {waterpost_name}, with freq_name:{freq_name} saved')
            # ----------------------------
        else:
            dct = {
                'hstst_id': int(waterpost_name), 
                'lat_y': y_lat,
                'lon_x': x_lon,
                'frequency': frequency, 
                'wtrdepth': wtrdepth, 
                'wtrlvltime': wtrlvltime,
                'lenth': lenth,
                'flooded_area': 0,
                'success': 0
            }
            self.df_new = self.df_new.append(dct, ignore_index=True)
            self.df_new.to_csv(f'{save_path}/{self.dt_string}_tifs_for_shape_calculated.csv', sep=';', decimal=',', index=False)


    # New function for computing flood shape
    def compute_shape_with_flood_height(self, save_path, data_dict, uniform_flooding=False):
        coordinate = data_dict['coordinate']
        freq_name = data_dict['frequency_name']
        frequency = data_dict['frequency']
        waterpost_name = int(data_dict['hsts_id'])
        wtrdepth = data_dict['wtrdepth']
        wtrlvltime = data_dict['wtrlvltime']
        compression = data_dict['compression']
        lenth = data_dict['lenth']
        target_h = data_dict['wtrdepth']

        # Define coordinate of map to download
        (x_lon, y_lat) = coordinate

        # Get WaterShed (precomputed or computed on demand)
        shed = self._compute_shed(x_lon, y_lat, compression)

        lng_num, lat_num = int(x_lon), int(y_lat)
        if shed is not None:
            top_left = (lng_num - 1, lat_num + 2)
            bottom_right = (lng_num + 2, lat_num - 1)

            GraphClass = Graph.Graph(dem=shed.dem, fdir=shed.fdir, acc=shed.acc, compression=shed.compression_coef)
            (shape,
             flooded_nodes_down,
             flooded_nodes_down_height,
             flooded_nodes_up,
             flooded_nodes_up_height,
             all_out_nodes,
             dem,
             height,
             point_coord,
             offset
             ) = GraphClass.compute_flood_with_height(coordinate, top_left, bottom_right, lenth, target_h,
                                          uniform_flooding=uniform_flooding)

            new_space_no_interpol = np.zeros((shape[0], shape[1]), dtype=np.uint8)
            for (down_node, _height) in zip(flooded_nodes_down, flooded_nodes_down_height):
                y, x = down_node[0], down_node[1]
                new_space_no_interpol[y, x] = _height

            for (up_node, _height) in zip(flooded_nodes_up, flooded_nodes_up_height):
                y, x = up_node[0], up_node[1]
                new_space_no_interpol[y, x] = _height

            flooded_area = (len(flooded_nodes_down) + len(flooded_nodes_up)) * ((0.03 * compression) ** 2)
            flood_outer = new_space_no_interpol

            main_point = point_coord
            top_left, bottom_right = top_left, bottom_right
            main_shape = shed.dem.shape

            y_path, x_path = (top_left[1] - bottom_right[1]) / main_shape[0], (bottom_right[0] - top_left[0]) / \
                             main_shape[1]

            cut_top_left = (
            top_left[1] - (main_point[0] - offset) * y_path, top_left[0] + (main_point[1] - offset) * x_path)
            cut_bottom_right = (
            top_left[1] - (main_point[0] + offset) * y_path, top_left[0] + (main_point[1] + offset) * x_path)

            y_path, x_path = (cut_top_left[0] - cut_bottom_right[0]) / flood_outer.shape[0], (
                        cut_bottom_right[1] - cut_top_left[1]) / flood_outer.shape[1]

            # ----------------------------
            print('Flood computed!')
            print('Saving it as .tif...')
            # ----------------------------

            nlines = flood_outer.shape[0]
            ncolumns = flood_outer.shape[1]
            data = flood_outer

            def getGeoTransform(extent, nlines, ncols):
                resx = (extent[2] - extent[0]) / ncols
                resy = (extent[3] - extent[1]) / nlines
                return [extent[0], resx, 0, extent[3], 0, -resy]

            # Define the data extent (min. lon, min. lat, max. lon, max. lat)
            extent = [cut_top_left[1], cut_bottom_right[0], cut_bottom_right[1], cut_top_left[0]]

            # Export the test array to GeoTIFF ================================================

            # Get GDAL driver GeoTiff
            driver = gdal.GetDriverByName('GTiff')

            # Get dimensions
            nlines = data.shape[0]
            ncols = data.shape[1]
            nbands = len(data.shape)
            data_type = gdal.GDT_Int16  # gdal.GDT_Float32

            # Create a temp grid
            grid_data = driver.Create('grid_data', ncols, nlines, 1, data_type)

            # Write data for each bands
            grid_data.GetRasterBand(1).WriteArray(data)

            # Lat/Lon WSG84 Spatial Reference System
            srs = osr.SpatialReference()
            srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

            # Setup projection and geo-transform
            grid_data.SetProjection(srs.ExportToWkt())
            grid_data.SetGeoTransform(getGeoTransform(extent, nlines, ncols))

            # Save the file
            file_name = f'{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}.tif'
            try:
                os.makedirs(f"{save_path}/{waterpost_name}")
                driver.CreateCopy(file_name, grid_data, 0)
                print(f'Generated GeoTIFF: {file_name}')
            except:
                driver.CreateCopy(file_name, grid_data, 0)

            # Close the file
            driver = None
            grid_data = None

            # Delete the temp grid
            os.remove('grid_data')

            dct = {
                'hstst_id': int(waterpost_name),
                'lat_y': y_lat,
                'lon_x': x_lon,
                'frequency': frequency,
                'wtrdepth': wtrdepth,
                'wtrlvltime': wtrlvltime,
                'lenth': lenth,
                'flooded_area': flooded_area,
                'success': 1
            }
            self.df_new = self.df_new.append(dct, ignore_index=True)
            self.df_new.to_csv(f'{save_path}/{self.dt_string}_tifs_for_shape_calculated.csv', sep=';', decimal=',',
                               index=False)

            # ----------------------------
            print(f'Tif for hstation_id: {waterpost_name}, with freq_name:{freq_name} saved')
            # ----------------------------
        else:
            dct = {
                'hstst_id': int(waterpost_name),
                'lat_y': y_lat,
                'lon_x': x_lon,
                'frequency': frequency,
                'wtrdepth': wtrdepth,
                'wtrlvltime': wtrlvltime,
                'lenth': lenth,
                'flooded_area': 0,
                'success': 0
            }
            self.df_new = self.df_new.append(dct, ignore_index=True)
            self.df_new.to_csv(f'{save_path}/{self.dt_string}_tifs_for_shape_calculated.csv', sep=';', decimal=',',
                               index=False)
        
    
    
    def compute_tifs_for_shapes(self, csv_data_path, DEMs_path, save_path, uniform_flooding=False):
        # ---- Guard ---- 
        # TODO: Perform this shape functions
        # guard.data_is_not_none(data)
        # guard.data_contains_values(data

        # ---- Compute ----
        self.DEMs_path = DEMs_path
        
        df = pd.read_csv(csv_data_path, sep=';', decimal=',')
        df['x_lon_int'] = df['lon_x'].values.astype(int)
        df['y_lat_int'] = df['lat_y'].values.astype(int)
        
        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis = 0, ascending = True, inplace = True, na_position = "first")
        
        # Create success df
        self.df_new = pd.DataFrame(columns=['hstst_id', 'lat_y', 'lon_x', 'frequency', 'wtrdepth', 'wtrlvltime', 'lenth', 'flooded_area', 'success'])
        self.df_new = self.df_new.astype(dtype= {'hstst_id': 'int64', 
                                                 'lat_y': 'float64', 
                                                 'lon_x': 'float64',
                                                 'frequency': 'float64', 
                                                 'wtrdepth': 'float64', 
                                                 'wtrlvltime': 'float64',
                                                 'lenth': 'int64',
                                                 'flooded_area': 'float64',
                                                 'success': 'int64'})
    
        # For future save
        self.dt_string = datetime.now(pytz.timezone('Europe/Moscow')).strftime("%d_%m_%Y__%H:%M")
    
        unique_id = df['hstst_id'].unique()
        for id in tqdm(unique_id):
            temp_df = df[df['hstst_id'] == id]

            for i, row in temp_df.iterrows():
                data_dict = {
                    'hsts_id': id,
                    'coordinate': (row[2], row[1]),
                    'frequency_name': self.__frequency_to_name(round(row[3], 1)),
                    'frequency': round(row[3], 1),
                    "wtrdepth": round(row[5], 2),
                    'wtrlvltime': round(row[6], 2),
                    'compression': int(row[8]),
                    'lenth': int(row[7])
                }
                self.compute_shape(save_path, data_dict, uniform_flooding=uniform_flooding)
                
        # ----------------------------
        print('DONE')
        # ----------------------------

    # New function to store flood as flooding height, not flooding bool
    def compute_tifs_for_shapes_with_height(self, csv_data_path, DEMs_path, save_path, uniform_flooding=False):
        # ---- Guard ----
        # TODO: Perform this shape functions
        # guard.data_is_not_none(data)
        # guard.data_contains_values(data

        # ---- Compute ----
        self.DEMs_path = DEMs_path

        df = pd.read_csv(csv_data_path, sep=';', decimal=',')
        df['x_lon_int'] = df['lon_x'].values.astype(int)
        df['y_lat_int'] = df['lat_y'].values.astype(int)

        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis=0, ascending=True, inplace=True, na_position="first")

        # Create success df
        self.df_new = pd.DataFrame(
            columns=['hstst_id', 'lat_y', 'lon_x', 'frequency', 'wtrdepth', 'wtrlvltime', 'lenth', 'flooded_area',
                     'success'])
        self.df_new = self.df_new.astype(dtype={'hstst_id': 'int64',
                                                'lat_y': 'float64',
                                                'lon_x': 'float64',
                                                'frequency': 'float64',
                                                'wtrdepth': 'float64',
                                                'wtrlvltime': 'float64',
                                                'lenth': 'int64',
                                                'flooded_area': 'float64',
                                                'success': 'int64'})

        # For future save
        self.dt_string = datetime.now(pytz.timezone('Europe/Moscow')).strftime("%d_%m_%Y__%H:%M")

        unique_id = df['hstst_id'].unique()
        for id in tqdm(unique_id):
            temp_df = df[df['hstst_id'] == id]

            for i, row in temp_df.iterrows():
                data_dict = {
                    'hsts_id': id,
                    'coordinate': (row[2], row[1]),
                    'frequency_name': self.__frequency_to_name(round(row[3], 1)),
                    'frequency': round(row[3], 1),
                    "wtrdepth": round(row[5], 2),
                    'wtrlvltime': round(row[6], 2),
                    'compression': int(row[8]),
                    'lenth': int(row[7])
                }
                self.compute_shape_with_flood_height(save_path, data_dict, uniform_flooding=uniform_flooding)

        # ----------------------------
        print('DONE')
        # ----------------------------
                

        
    def add_fields(self, dst_layer):
        # Создаю записи филдов, которые должны быть в шейп файле и задаю их тип
        field_name = ogr.FieldDefn("hsts_id", ogr.OFTInteger)
        field_name2 = ogr.FieldDefn("lat_y", ogr.OFTReal)
        field_name2.SetPrecision(6)
        field_name3 = ogr.FieldDefn("lon_x", ogr.OFTReal)
        field_name3.SetPrecision(6)
        field_name4 = ogr.FieldDefn("frequency", ogr.OFTReal)
        field_name4.SetPrecision(2)
        field_name5 = ogr.FieldDefn("wtrdepth", ogr.OFTReal)
        field_name5.SetPrecision(2)
        field_name6 = ogr.FieldDefn("wtrlvltime", ogr.OFTReal)
        field_name6.SetPrecision(2)

        # Создаю эти филды
        dst_layer.CreateField(field_name)
        dst_layer.CreateField(field_name2)
        dst_layer.CreateField(field_name3)
        dst_layer.CreateField(field_name4)
        dst_layer.CreateField(field_name5)
        dst_layer.CreateField(field_name6)

        
        
    def _concat_shapes(self, df, tifs_path, save_path):
        print('computing shape...')
        
        # this allows GDAL to throw Python Exceptions
        gdal.UseExceptions()
        dst_layername = 'result_union'
        dst_drv = ogr.GetDriverByName("ESRI Shapefile")
        dst_file = f'{save_path}/{self.dt_string}_{dst_layername}.shp'

        if os.path.exists(dst_file):
            dst_drv.DeleteDataSource(dst_file)

        dst_ds = dst_drv.CreateDataSource(dst_file)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dst_layer = dst_ds.CreateLayer(dst_layername, geom_type=ogr.wkbMultiPolygon, srs = srs )
        self.add_fields(dst_layer)

        dst_layer_defn = dst_layer.GetLayerDefn()

        tmp_layername = 'temp'
        mem_drv = ogr.GetDriverByName("MEMORY")
        
        for _, row in tqdm(df.iterrows()):
            hsts_id = int(row[0])
            frequency_name = self.__frequency_to_name(round(row[3], 1))
            tmp_ds = mem_drv.CreateDataSource('mem_temp_data')
            tmp_layer = tmp_ds.CreateLayer(tmp_layername, geom_type=ogr.wkbPolygon, srs = srs )
            
            raster_file = f'{tifs_path}/{hsts_id}/{hsts_id}_{frequency_name}.tif'
            
            if os.path.exists(raster_file) == False:
                continue
            
            src_ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
            srcband = src_ds.GetRasterBand(1)
            gdal.Polygonize(srcband, srcband, tmp_layer, -1, [], callback=None )
            
            multi_poly = ogr.Geometry(ogr.wkbMultiPolygon)

            for poly in tmp_layer:
                geom = poly.GetGeometryRef()
                if geom.IsValid():
                    multi_poly = multi_poly.Union(geom)
                else:
                    clean = geom.Buffer(0)
                    multi_poly = multi_poly.Union(clean)

            out_feature = ogr.Feature(dst_layer_defn)
            out_feature.SetGeometry(multi_poly)

            out_feature.SetField('hsts_id', hsts_id)
            out_feature.SetField('lat_y', row[1])
            out_feature.SetField('lon_x', row[2])
            out_feature.SetField('frequency', round(row[3], 1))
            out_feature.SetField('wtrdepth', round(row[5], 2))
            out_feature.SetField('wtrlvltime', round(row[6], 2))

            dst_layer.CreateFeature(out_feature)

            del out_feature
            del multi_poly
            del tmp_layer
            del tmp_ds

        del dst_layer
        del dst_ds

        # ----------------------------
        print('DONE')
        # ----------------------------
        
        
        
    def compute_union_shape(self, csv_data_path, tifs_path, save_path):
        # ---- Guard ---- 
        # TODO: Perform this shape functions
        # guard.data_is_not_none(data)
        # guard.data_contains_values(data

        # ---- Compute ----
        self.tifs_path = tifs_path
        self.dt_string = datetime.now(pytz.timezone('Europe/Moscow')).strftime("%d_%m_%Y__%H:%M")
        
        df = pd.read_csv(csv_data_path, sep=';', decimal=',')
        df['x_lon_int'] = df['lon_x'].values.astype(int)
        df['y_lat_int'] = df['lat_y'].values.astype(int)
        
        # Sort df by x_lon and y_lat for future reduction of DEM computing
        df.sort_values(['x_lon_int', 'y_lat_int'], axis = 0, ascending = True, inplace = True, na_position = "first")
        
        self._concat_shapes(df, tifs_path, save_path)
        
        # ----------------------------
        print('DONE')
        # ----------------------------
        
        
        
    def compute_heights(self, excel_path, DEMs_path, save_path, min_acc=1000):
        # ---- Guard ---- 
        df = pd.read_csv(excel_path, sep=';', decimal=',')
        guard.height_data_contains_columns(df)
          
        # ---- Compute ----
        GraphClass = Graph.Graph()
        GraphClass.compute_height(excel_path, DEMs_path, min_acc, save_path)
        
        # ----------------------------
        print('DONE')
        # ----------------------------
        
        
        
    def compute_river_slices(self, excel_path, DEMs_path, save_path):
        # ---- Guard ---- 
#         df = pd.read_csv(excel_path, sep=';', decimal=',')
#         guard.height_data_contains_columns(df)
          
        # ---- Compute ----
        GraphClass = Graph.Graph()
        GraphClass.get_river_slice(excel_path, DEMs_path, save_path)
        
        # ----------------------------
        print('DONE')
        # ----------------------------
        
        
        
    def compute_height_for_coordinate(self, DEMs_path, coordinate):
        # ---- Compute ----
        (x_lon, y_lat) = coordinate
        lng_num, lat_num = int(x_lon), int(y_lat)
        lat = str(lat_num)
        lng = ''.join((['0'] + list(str(int(lng_num))))[-3:])
        file_name = f'n{lat}_e{lng}_1arc_v3.tif' if lat_num+1 < 60 else f'n{lat}_e{lng}_1arc_v3_1201x1201.tif'
        path = f'{DEMs_path}/{file_name}'
        
        Shed = WaterShed.WaterSheds(file_path=path, compute_acc=True)
        top_left = (lng_num, lat_num+1)
        bottom_right = (lng_num+1, lat_num)
        
        return Shed.dem[Shed.coordinate2point(coordinate, top_left, bottom_right)]
      

      
    def compute_watershed_area(self, DEMs_path, csv_data_path=None, save_path=None, coordinate=None, top_left=None, bottom_right=None, compress_coef=1):

        if (coordinate is not None) and (top_left is not None) and (bottom_right is not None):
            # ---- Compute path if only one point passed ---- 
            x_path, y_path = bottom_right[0] - top_left[0], top_left[1] - bottom_right[1]

            tif_pathes = []
            for i in range(y_path):
                for j in range(x_path):
                    lat = str(bottom_right[1] + i)
                    lng = ''.join((['0'] + list(str(int(top_left[0] + j))))[-3:])
                    file_name = f'n{lat}_e{lng}_1arc_v3.tif' if top_left[1]+1 < 60 else f'n{lat}_e{lng}_1arc_v3_1201x1201.tif'
                    tif_pathes.append(f'{DEMs_path}/{file_name}')

            success_list = []
            for tif_path in tif_pathes:
                if path.exists(tif_path) == False:
                    print(f'{tif_path} is not exist in path {DEMs_path}')
                    success_list.append(False)

            if len(success_list) != 0:
                return None

            # ---- Compute ---- 
            # Download DEM and preprocess it
            shed = WaterShed.WaterSheds(files_pathes=tif_pathes, compute_acc=True, compression=compress_coef)

            watershed_area = shed.compute_watershed_area(coordinate, top_left, bottom_right)

            return f'{watershed_area} км^2'
          
          
        elif (csv_data_path is not None) and (save_path is not None):
            df = pd.read_csv(csv_data_path, sep=';', decimal=',')
            
            # Create success df
            df_new = pd.DataFrame(columns=['hstst_id', 'lat', 'lon', 'watershed_area', 'success'])
            df_new = df_new.astype(dtype= {'hstst_id': 'int64', 
                                                'lat': 'float64', 
                                                'lon': 'float64',
                                                'watershed_area': 'float64',
                                                'success': 'int64'})

            # For future save
            dt_string = datetime.now(pytz.timezone('Europe/Moscow')).strftime("%d_%m_%Y__%H:%M")
            
            # iterate over points to calculate water area
            df['top_left_bottom_right'] = df['top_left'] + df['bottom_right']
            unique_areas = df['top_left_bottom_right'].unique()
            for id in tqdm(unique_areas):
                temp_df = df[df['top_left_bottom_right'] == id]
                dem_computed = False
                
                for i, row in temp_df.iterrows():
                    hstst_id = row[0]
                    lat = float(row[1])
                    lon = float(row[2])
                    top_left = tuple(map(int, row[3].replace('(', '').replace(')', '').split(', ')))
                    bottom_right = tuple(map(int, row[4].replace('(', '').replace(')', '').split(', ')))
                    compress_coef = row[5]
                    coordinate = (lat, lon)
                    
                    print(f'Processing hstst_id {hstst_id}...')
                    
                    # ---- Compute path for point ----
                    if dem_computed == False:
                        x_path, y_path = bottom_right[0] - top_left[0], top_left[1] - bottom_right[1]

                        tif_pathes = []
                        for i in range(y_path):
                            for j in range(x_path):
                                lat = str(bottom_right[1] + i)
                                lon = ''.join((['0'] + list(str(int(top_left[0] + j))))[-3:])
                                file_name = f'n{lat}_e{lon}_1arc_v3.tif' if top_left[1]+1 <= 60 else f'n{lat}_e{lon}_1arc_v3_1201x1201.tif'
                                tif_pathes.append(f'{DEMs_path}/{file_name}')

                        success_list = []
                        for tif_path in tif_pathes:
                            if path.exists(tif_path) == False:
                                print(f'{tif_path} is not exist in path {DEMs_path}')
                                success_list.append(False)

                        if len(success_list) != 0:
                            return None

                        # ---- Compute ---- 
                        # Download DEM and preprocess it
                        self.shed = WaterShed.WaterSheds(files_pathes=tif_pathes, compute_acc=True, compression=compress_coef)
                        dem_computed = True

                    
                    watershed_area = self.shed.compute_watershed_area(coordinate, top_left, bottom_right)
                    print(f'Watershed area for hstst_id {hstst_id} was calculated')

                    # Save result
                    dct = {
                        'hstst_id': int(hstst_id), 
                        'lat': lat,
                        'lon': lon,
                        'watershed_area': watershed_area,
                        'success': 1
                    }
                    df_new = df_new.append(dct, ignore_index=True)
                    df_new.to_csv(f'{save_path}/{dt_string}_watershed_area.csv', sep=';', decimal=',', index=False)
                    
                del self.shed, dem_computed
                
#             for i, row in df.iterrows():
#                 hstst_id = row[0]
#                 lat = float(row[1])
#                 lon = float(row[2])
#                 top_left = tuple(map(int, row[3].replace('(', '').replace(')', '').split(', ')))
#                 bottom_right = tuple(map(int, row[4].replace('(', '').replace(')', '').split(', ')))
#                 compress_coef = row[5]
#                 coordinate = (lat, lon)
                
#                 print(f'Processing hstst_id {hstst_id}...')
                
#                 # ---- Compute path for point ---- 
#                 x_path, y_path = bottom_right[0] - top_left[0], top_left[1] - bottom_right[1]

#                 tif_pathes = []
#                 for i in range(y_path):
#                     for j in range(x_path):
#                         lat = str(bottom_right[1] + i)
#                         lon = ''.join((['0'] + list(str(int(top_left[0] + j))))[-3:])
#                         file_name = f'n{lat}_e{lon}_1arc_v3.tif' if top_left[1]+1 <= 60 else f'n{lat}_e{lon}_1arc_v3_1201x1201.tif'
#                         tif_pathes.append(f'{DEMs_path}/{file_name}')

#                 success_list = []
#                 for tif_path in tif_pathes:
#                     if path.exists(tif_path) == False:
#                         print(f'{tif_path} is not exist in path {DEMs_path}')
#                         success_list.append(False)

#                 if len(success_list) != 0:
#                     return None

#                 # ---- Compute ---- 
#                 # Download DEM and preprocess it
#                 shed = WaterShed.WaterSheds(files_pathes=tif_pathes, compute_acc=True, compression=compress_coef)

#                 watershed_area = shed.compute_watershed_area(coordinate, top_left, bottom_right)
#                 print(f'Watershed area for hstst_id {hstst_id} was calculated')
                
#                 # Save result
#                 dct = {
#                     'hstst_id': int(hstst_id), 
#                     'lat': lat,
#                     'lon': lon,
#                     'watershed_area': watershed_area,
#                     'success': 1
#                 }
#                 df_new = df_new.append(dct, ignore_index=True)
#                 df_new.to_csv(f'{save_path}/{dt_string}_watershed_area.csv', sep=';', decimal=',', index=False)
                
            
        else:
            print('Not all veriables passed')
            print('Try using csv_data_path={csv_path}, save_path={save_path} for file with miltiple inputs')
            print('Try using coordinate={point_coordinate}, top_left={top_left_coordinate}, bottom_right={bottom_right_coordinate} for one point input')
