from tracer_tools import tools
from tracer_tools import app_tools
from pprint import pprint
from page_design import map_face
import os

def main():
    polygon = [[[151.301454, -33.700754],
            [151.311453, -33.702075],
            [151.307237, -33.739761],
            [151.294220, -33.736329],
            [151.301454, -33.700754]]]
    
    polygon = tools.smallest_rectangle(polygon)
    dates = ['2024-12-31', '2025-06-09']
    sitename = 'NARRA'
    prj_id = 'gee-playground-jas13'
    filepath_data = os.path.join(os.getcwd(), 'data')
    
    inputs = {
    'polygon': polygon,
    'dates': dates,
    'sat_list': ['S2'],
    'sitename': sitename,
    'ee_project_id': prj_id,
    'filepath': filepath_data,
    # 'LandsatWRS': '089083',
    # 'S2tile': '56HLH',
        }
    
    img = tools.check_images_available(inputs) #returns t1 and t2 images in a tuple
    s2Img = img[0]['S2']

    #the structure goes like this:
    # LIST of IMGS
    #   - DICT of IMG Metadata
    #       - bands (list of bands)
    #       - id (earth engine dataset id)
    #       - properties (image params)
    #       - type (datatype e.g. Image)
    #       - version ID (version of the image)

    for img in s2Img:
        print(app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d'))

def main2():
    M = map_face(13.198839454428844, 122.30218706556167)
    
if __name__ == "__main__":
    main()