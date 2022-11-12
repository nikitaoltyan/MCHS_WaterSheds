import WaterShed
import Graph

class Main():
  
    def compute_shape(self, file_path, data_dict):
        shed = WaterShed(file_path, compute_acc=True)
        
        coordinate = data_dict['coordinate']
        top_left = data_dict['top_left']
        bottom_right = data_dict['bottom_right']
        lenth = data_dict['lenth']
        target_h = data_dict['target_h']

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
        
        print('Done')
