import pandas as pd


def data_is_not_none(data):
    assert data != None, 'Data is not specified'

    
def data_contains_values(data_dict):
    """Check if given dict of dicts contains all required keys
    """
    for data in data_dict:
        assert 'coordinate' in data, f'Missing coordinate for {data}'
        assert 'frequency_name' in data, f'Missing frequency_name for {data}'
        assert 'frequency' in data, f'Missing frequency for {data}'
        assert 'hsts_id' in data, f'Missing hsts_id for {data}'
        assert 'lenth' in data, f'Missing lenth for {data}'
        assert 'wtrdepth' in data, f'Missing wtrdepth for {data}'
        assert 'wtrlvltime' in data, f'Missing wtrlvltime for {data}'
    

def height_data_contains_columns(file_path):
    df = pd.read_csv(file_path, sep=';', decimal=',')
    columns = df.columns
    dtypes = df.dtypes
    assert len(columns) == 3, "Must be 3 columns: ['hstation_id', 'x_lon', 'y_lat']"
    assert columns[0] == 'hstation_id', "Column 1 should have name 'hstation_id'"
    assert columns[1] == 'x_lon', "Column 3 should have name 'x_lon'"
    assert columns[2] == 'y_lat', "Column 4 should have name 'y_lat'"
    # Check long and lat dtypes for future sorting
    assert dtypes[1] == 'float64', "x_lon should have a float64 type"
    assert dtypes[2] == 'float64', "y_lat should have a float64 type"
    
