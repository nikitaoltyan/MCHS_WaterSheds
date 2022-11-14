def data_is_not_none(data):
    assert data != None, 'Data not specified'

    
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
    
