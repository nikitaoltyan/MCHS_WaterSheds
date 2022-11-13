import WaterShed
import Graph

from osgeo import gdal, osr, ogr # Python bindings for GDAL
import numpy as np
import os

class Main():
  
    def compute_shape(self, tif_path, save_path, data_dict):
        shed = WaterShed.WaterSheds(tif_path, compute_acc=True)
        
        coordinate = data_dict['coordinate']
        freq_name = data_dict['frequency_name']
        waterpost_name = data_dict['hsts_id']
        top_left = (int(coordinate[0]), int(coordinate[1])+1)
        bottom_right = (int(coordinate[0])+1, int(coordinate[1]))
        lenth = data_dict['lenth']
        target_h = data_dict['target_h']

        print('Computing flood for...')
        GraphClass = Graph.Graph(dem=shed.dem, fdir=shed.fdir, acc=shed.acc)
        (shape, 
            flooded_nodes_down, 
            flooded_nodes_up,
            all_out_nodes,
            dem, 
            height, 
            point_coord, 
            offset
        ) = GraphClass.compute_flood(coordinate, top_left, bottom_right, lenth, target_h)
        
        new_space_no_interpol = np.zeros((shape[0], shape[1]), dtype=np.uint8)
        for down_node in flooded_nodes_down:
            y, x = down_node[0], down_node[1]
            new_space_no_interpol[y, x] = 1

        for up_node in flooded_nodes_up:
            y, x = up_node[0], up_node[1]
            new_space_no_interpol[y, x] = 1
            
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
            driver.CreateCopy(file_name, grid_data, 0)
            print(f'Generated GeoTIFF: {file_name}')
        except:
            try:
                os.mkdir(f"{save_path}/{waterpost_name}")
            except:
                pass
            os.mkdir(f"{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}")
            driver.CreateCopy(file_name, grid_data, 0)
            print(f'Generated GeoTIFF in new folder: {file_name}')


        # Close the file
        driver = None
        grid_data = None

        # Delete the temp grid               
        os.remove('grid_data')
        
        # ----------------------------
        print('Tif saved')
        # ----------------------------

        
    def add_fields(self, dst_layer):
        # Создаю записи филдов, которые должны быть в шейп файле и задаю их тип
        field_name = ogr.FieldDefn("region_id", ogr.OFTInteger)
        field_name2 = ogr.FieldDefn("LAT_Y", ogr.OFTReal)
        field_name2.SetPrecision(6)
        field_name3 = ogr.FieldDefn("LON_X", ogr.OFTReal)
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
