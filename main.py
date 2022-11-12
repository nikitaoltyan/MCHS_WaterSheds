import WaterShed
import Graph

class Main
  
  def compute_for_all_freqencies(file_path):
      shed = WaterSheds(file_path, compute_acc=True)
      print(shed.shape)
