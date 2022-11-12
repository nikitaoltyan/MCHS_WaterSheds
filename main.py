import WaterShed
import Graph

class Main():
  
  def compute_for_all_freqencies(file_path):
      shed = WaterShed(file_path, compute_acc=True)
      print(shed.shape)
