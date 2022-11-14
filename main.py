import WaterShed
import Graph

from osgeo import gdal, osr, ogr # Python bindings for GDAL
import numpy as np
import os
import sys

class Main():
  
    def __init__(self):
        self.shed = None
        
        
    def _compute_shed(self, file_path):
        """Compute shed if it isn't exists
        """
        if self.shed == None:
            self.shed = WaterShed.WaterSheds(file_path, compute_acc=True)
            return self.shed
        else:
            return self.shed
  
  
    def compute_shape(self, tif_path, save_path, data_dict):
        shed = self._compute_shed(tif_path)
        
        coordinate = data_dict['coordinate']
        freq_name = data_dict['frequency_name']
        waterpost_name = data_dict['hsts_id']
        top_left = (int(coordinate[0]), int(coordinate[1])+1)
        bottom_right = (int(coordinate[0])+1, int(coordinate[1]))
        lenth = data_dict['lenth']
        target_h = data_dict['wtrdepth']

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
        file_name = f'{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}/{waterpost_name}_{freq_name}.tif'
        try:
            os.makedirs(f"{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}")
            driver.CreateCopy(file_name, grid_data, 0)
            print(f'Generated GeoTIFF: {file_name}')
        except:
            driver.CreateCopy(file_name, grid_data, 0)


        # Close the file
        driver = None
        grid_data = None

        # Delete the temp grid               
        os.remove('grid_data')
        
        # ----------------------------
        print('Tif saved')
        # ----------------------------
        
    
    def compute_shapes(self, tif_path, save_path, data):
        for d in data:
            print(d)
            self.compute_shape(tif_path, save_path, data[d])
        self._concat_shapes(data, save_path)

        
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

        
    def _concat_shapes(self, data, save_path):
        print('computing shape...')
        waterpost_name = data[list(data.keys())[0]]['hsts_id']
        
        # this allows GDAL to throw Python Exceptions
        gdal.UseExceptions()
        dst_layername = 'result_union'
        dst_drv = ogr.GetDriverByName("ESRI Shapefile")
        dst_file = f'{save_path}/{waterpost_name}/{waterpost_name}_{dst_layername}.shp'

        if os.path.exists(dst_file):
            dst_drv.DeleteDataSource(dst_file)

        dst_ds = dst_drv.CreateDataSource(dst_file)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dst_layer = dst_ds.CreateLayer(dst_layername, geom_type=ogr.wkbMultiPolygon, srs = srs )
        add_fields(dst_layer)

        dst_layer_defn = dst_layer.GetLayerDefn()

        tmp_layername = 'temp'
        mem_drv = ogr.GetDriverByName("MEMORY")

        for file_name, file_data in data.items():
            print(file_name)
            waterpost_name = file_data['hsts_id']
            freq_name = file_data['frequency_name']
            tmp_ds = mem_drv.CreateDataSource('mem_temp_data')
            tmp_layer = tmp_ds.CreateLayer(tmp_layername, geom_type=ogr.wkbPolygon, srs = srs )

            raster_file = f'{save_path}/{waterpost_name}/{waterpost_name}_{freq_name}/{file_name}'
            src_ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
            srcband = src_ds.GetRasterBand(1)
            gdal.Polygonize(srcband, srcband, tmp_layer, -1, [], callback=None )

            multi_poly = ogr.Geometry(ogr.wkbMultiPolygon)

            for poly in tmp_layer:
                multi_poly = multi_poly.Union(poly.geometry())

            out_feature = ogr.Feature(dst_layer_defn)
            out_feature.SetGeometry(multi_poly)

            out_feature.SetField('region_id', waterpost_name)
            out_feature.SetField('LAT_Y', file_data['coordinate'][1])
            out_feature.SetField('LON_X', file_data['coordinate'][0])
            out_feature.SetField('frequency', file_data['frequency'])
            out_feature.SetField('wtrdepth', file_data['wtrdepth'])
            out_feature.SetField('wtrlvltime', file_data['wtrlvltime'])

            dst_layer.CreateFeature(out_feature)

            del out_feature
            del multi_poly
            del tmp_layer
            del tmp_ds

        del dst_layer
        del dst_ds
