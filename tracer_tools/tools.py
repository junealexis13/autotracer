from datetime import datetime
from shapely import geometry
from scipy import stats, interpolate
import numpy as np
import pandas as pd
import pyproj
import matplotlib.pyplot as plt
import ee
import os
import pickle
import pytz
import pyproj



def authenticate_and_initialize(prject_id=None):
    """
    Authenticates and initializes the Earth Engine API.
    This function handles the authentication and initialization process:
        1. Try to use existing token to initialize
        2. If 1 fails, try to refresh the token using Application Default Credentials
        3. If 2 fails, authenticate manually via the web browser
    """
    # first try to initialize connection with GEE server with existing token
    try: 
        ee.Initialize(project=prject_id)
        print('GEE initialized (existing token).')
    except:
        # if token is expired, try to refresh it
        # based on https://stackoverflow.com/questions/53472429/how-to-get-a-gcp-bearer-token-programmatically-with-python
        try:
            import google.auth
            import google.auth.transport.requests
            creds, project = google.auth.default()
            # creds.valid is False, and creds.token is None
            # refresh credentials to populate those
            auth_req = google.auth.transport.requests.Request()
            creds.refresh(auth_req)
            # initialise GEE session with refreshed credentials
            ee.Initialize(creds)
            print('GEE initialized (refreshed token).')
        except:
            # get the user to authenticate manually and initialize the sesion
            ee.Authenticate()
            ee.Initialize(project=prject_id)
            print('GEE initialized (manual authentication).')

def check_images_available(inputs):
    """
    Scan the GEE collections to see how many images are available for each
    satellite mission (L5,L7,L8,L9,S2), collection (C02) and tier (T1,T2).
    
    Note: Landsat Collection 1 (C01) is deprecated. Users should migrate to Collection 2 (C02).
    For more information, visit: https://developers.google.com/earth-engine/landsat_c1_to_c2

    KV WRL 2018

    Arguments:
    -----------
    inputs: dict
        inputs dictionnary

    Returns:
    -----------
    im_dict_T1: list of dict
        list of images in Tier 1 and Level-1C
    im_dict_T2: list of dict
        list of images in Tier 2 (Landsat only)
    """

    dates = [datetime.strptime(_,'%Y-%m-%d') for _ in inputs['dates']]
    dates_str = inputs['dates']
    polygon = inputs['polygon']
    
    # check if dates are in chronological order
    if  dates[1] <= dates[0]:
        raise Exception('Verify that your dates are in the correct chronological order')

    # check if EE was initialised or not
    authenticate_and_initialize(prject_id=inputs['ee_project_id'])
        
    print('Number of images available between %s and %s:'%(dates_str[0],dates_str[1]), end='\n')
    
    # get images in Landsat Tier 1 as well as Sentinel Level-1C
    # print('- In Landsat Tier 1 & Sentinel-2 Level-1C:')
    im_dict_T1 = dict([])
    sum_img = 0

    # filter S2 collection only - removed Landsat imges from COAST SAT Algo
    im_list = get_image_info('COPERNICUS/S2_HARMONIZED',polygon,dates_str)

    # for S2, filter collection to only keep images with same UTM Zone projection 
    # (there duplicated images in different UTM projections)
    im_list = filter_S2_collection(im_list)

    sum_img = sum_img + len(im_list)
    print('     %s: %d images'%('S2',len(im_list)))
    im_dict_T1['S2'] = im_list          
        
    print('  Total to download: %d images'%sum_img)

    # check if images already exist  
    # print('\nLooking for existing imagery...')
    filepath = os.path.join(inputs['filepath'],inputs['sitename'])
    if os.path.exists(filepath):
        metadata_existing = get_metadata(inputs)
        for satname in inputs['sat_list']:
            # remove from download list the images that are already existing
            if satname in metadata_existing:
                if len(metadata_existing[satname]['dates']) > 0:
                    # get all the possible availabe dates for the imagery requested
                    avail_date_list = [datetime.fromtimestamp(image['properties']['system:time_start'] / 1000, tz=pytz.utc).replace( microsecond=0) for image in im_dict_T1[satname]]
                    # if no images are available, skip this loop
                    if len(avail_date_list) == 0:
                        print(f'{satname}:There are {len(avail_date_list)} images available, {len(metadata_existing[satname]["dates"])} images already exist, {len(avail_date_list)} to download')
                        continue
                    # get the dates of the images that are already downloaded
                    downloaded_dates = metadata_existing[satname]['dates']
                    # if no images are already downloaded, skip this loop and use whats already in im_dict_T1[satname]
                    if len(downloaded_dates) == 0:
                        print(f'{satname}:There are {len(avail_date_list)} images available, {len(downloaded_dates)} images already exist, {len(avail_date_list)} to download')
                        continue
                    # get the indices of the images that are not already downloaded 
                    idx_new = np.where([ not avail_date in downloaded_dates for avail_date in avail_date_list])[0]
                    im_dict_T1[satname] = [im_dict_T1[satname][index] for index in idx_new]
                    print('%s: %d/%d images already exist, %s to be downloaded'%(satname,
                                                                            len(downloaded_dates),
                                                                            len(avail_date_list),
                                                                            len(idx_new)))

    # if only S2 is in sat_list, stop here as no Tier 2 for Sentinel
    if len(inputs['sat_list']) == 1 and inputs['sat_list'][0] == 'S2':
        return im_dict_T1, []

    # if user also requires Tier 2 images, check the T2 collections as well
    col_names_T2 = {'L5':'LANDSAT/LT05/C02/T2_TOA',
                    'L7':'LANDSAT/LE07/C02/T2_TOA',
                    'L8':'LANDSAT/LC08/C02/T2_TOA'}
    print('- In Landsat Tier 2 (not suitable for time-series analysis):', end='\n')
    im_dict_T2 = dict([])
    sum_img = 0
    for satname in inputs['sat_list']:
        if satname in ['L9','S2']: continue # no Tier 2 for Sentinel-2 and Landsat 9
        im_list = get_image_info(col_names_T2[satname],polygon,dates_str)
        sum_img = sum_img + len(im_list)
        print('     %s: %d images'%(satname,len(im_list)))
        im_dict_T2[satname] = im_list

    print('  Total Tier 2: %d images'%sum_img)
    
    return im_dict_T1, im_dict_T2

def get_image_info(collection,polygon,dates):
    """
    Reads info about EE images for the specified collection, sensor and dates

    KV WRL 2022

    Arguments:
    -----------
    satname: str 
        name of the satellite mission [deprecated from CoastSat Library]
    polygon: list
        coordinates of the polygon in lat/lon
    dates: list of str
        start and end dates (e.g. '2022-01-01')

    Returns:
    -----------
    im_list: list of ee.Image objects
        list with the info for the images
    """
    while True:
        try:
            # get info about images
            ee_col = ee.ImageCollection(collection)
            # if S2tile is specified    
            col = ee_col.filterBounds(ee.Geometry.Polygon(polygon))\
                        .filterDate(dates[0],dates[1])
            # convert to dict
            im_list = col.getInfo().get('features')
            break
        except:
            continue
    # remove very cloudy images (>95% cloud cover)
    im_list = remove_cloudy_images(im_list, prc_cloud_cover=95) #this function now only works for s2 and only accepts imlist
    return im_list

def smallest_rectangle(polygon):
    """
    Converts a polygon to the smallest rectangle polygon with sides parallel
    to coordinate axes.
     
    KV WRL 2020

    Arguments:
    -----------
    polygon: list of coordinates 
        pair of coordinates for 5 vertices, in clockwise order,
        first and last points must match     
                
    Returns:    
    -----------
    polygon: list of coordinates
        smallest rectangle polygon
        
    """
    
    multipoints = geometry.Polygon(polygon[0])
    polygon_geom = multipoints.envelope
    coords_polygon = np.array(polygon_geom.exterior.coords)
    polygon_rect = [[[_[0], _[1]] for _ in coords_polygon]]
    return polygon_rect

def filter_S2_collection(im_list):
    """
    Removes duplicates from the EE collection of Sentinel-2 images (many duplicates)
    Finds the images that were acquired at the same time but have different utm zones.

    KV WRL 2018

    Arguments:
    -----------
    im_list: list
        list of images in the collection

    Returns:
    -----------
    im_list_flt: list
        filtered list of images
    """

    # get datetimes
    timestamps = [datetime.fromtimestamp(_['properties']['system:time_start']/1000,
                                         tz=pytz.utc) for _ in im_list]
    # get utm zone projections
    utm_zones = np.array([int(_['bands'][0]['crs'][5:]) for _ in im_list])
    if len(np.unique(utm_zones)) == 1:
        return im_list
    else:
        idx_max = np.argmax([np.sum(utm_zones == _) for _ in np.unique(utm_zones)])
        utm_zone_selected =  np.unique(utm_zones)[idx_max]
        # find the images that were acquired at the same time but have different utm zones
        idx_all = np.arange(0,len(im_list),1)
        idx_covered = np.ones(len(im_list)).astype(bool)
        idx_delete = []
        i = 0
        while 1:
            same_time = np.abs([(timestamps[i]-_).total_seconds() for _ in timestamps]) < 60*60*24
            idx_same_time = np.where(same_time)[0]
            same_utm = utm_zones == utm_zone_selected
            # get indices that have the same time (less than 24h apart) but not the same utm zone
            idx_temp = np.where([same_time[j] == True and same_utm[j] == False for j in idx_all])[0]
            idx_keep = idx_same_time[[_ not in idx_temp for _ in idx_same_time]]
            # if more than 2 images with same date and same utm, drop the last ones
            if len(idx_keep) > 2:
               idx_temp = np.append(idx_temp,idx_keep[-(len(idx_keep)-2):])
            for j in idx_temp:
                idx_delete.append(j)
            idx_covered[idx_same_time] = False
            if np.any(idx_covered):
                i = np.where(idx_covered)[0][0]
            else:
                break
        # update the collection by deleting all those images that have same timestamp
        # and different utm projection
        im_list_flt = [x for k,x in enumerate(im_list) if k not in idx_delete]
        # print('%d S2 duplicates removed'%(len(idx_delete)))
    return im_list_flt


def get_metadata(inputs):
    """
    Gets the metadata from the downloaded images by parsing .txt files located
    in the meta subfolder.

    KV WRL 2018

    Arguments:
    -----------
    inputs: dict with the following fields
        'sitename': str
            name of the site
        'filepath_data': str
            filepath to the directory where the images are downloaded

    Returns:
    -----------
    metadata: dict
        contains the information about the satellite images that were downloaded:
        date, filename, georeferencing accuracy and image coordinate reference system

    """
    # directory containing the images
    filepath = os.path.join(inputs['filepath'],inputs['sitename'])
    # initialize metadata dict
    metadata = dict([])
    # loop through the satellite missions
    for satname in ['L5','L7','L8','L9','S2']:
        # if a folder has been created for the given satellite mission
        if satname in os.listdir(filepath):
            # update the metadata dict
            metadata[satname] = {'filenames':[],'dates':[],'tilename':[],'epsg':[],'acc_georef':[],
                                 'im_quality':[],'im_dimensions':[]}
            # directory where the metadata .txt files are stored
            filepath_meta = os.path.join(filepath, satname, 'meta')
            # get the list of filenames and sort it chronologically
            filenames_meta = os.listdir(filepath_meta)
            filenames_meta.sort()
            # loop through the .txt files
            for im_meta in filenames_meta:
                # read the file and extract the metadata info
                df = pd.read_csv(os.path.join(filepath_meta, im_meta),sep='\t', 
                                 names=['property','value'])
                df.set_index('property', inplace=True)
                # if statement to be compatible with older version without tilename
                filename = df.at['filename','value']
                if 'tile' in df.index:
                    tilename = df.at['tile','value']
                else:
                    tilename = 'NA'
                epsg = df.at['epsg','value']    
                acc_georef = df.at['acc_georef','value'] 
                im_quality = df.at['image_quality','value'] 
                im_width = df.at['im_width','value'] 
                im_height = df.at['im_height','value'] 
                date_str = filename[0:19]
                date = pytz.utc.localize(datetime(int(date_str[:4]),int(date_str[5:7]),
                                                  int(date_str[8:10]),int(date_str[11:13]),
                                                  int(date_str[14:16]),int(date_str[17:19])))
                # check if they are quantitative values (Landsat) or Pass/Fail flags (Sentinel-2)
                try: acc_georef = float(acc_georef)
                except: acc_georef = str(acc_georef)
                try: im_quality = float(im_quality)
                except: im_quality = str(im_quality)
                # store the information in the metadata dict
                metadata[satname]['filenames'].append(filename)
                metadata[satname]['dates'].append(date)
                metadata[satname]['tilename'].append(tilename)
                metadata[satname]['epsg'].append(epsg)
                metadata[satname]['acc_georef'].append(acc_georef)
                metadata[satname]['im_quality'].append(im_quality)
                metadata[satname]['im_dimensions'].append([im_height,im_width])

    # save a .pkl file containing the metadata dict
    with open(os.path.join(filepath, inputs['sitename'] + '_metadata' + '.pkl'), 'wb') as f:
        pickle.dump(metadata, f)

    return metadata


def remove_cloudy_images(im_list, prc_cloud_cover=95):
    """
    Removes from the EE collection very cloudy images (>95% cloud cover)

    KV WRL 2018

    Arguments:
    -----------
    im_list: list
        list of images in the collection
    satname:
        name of the satellite mission
    prc_cloud_cover: int
        percentage of cloud cover acceptable on the images

    Returns:
    -----------
    im_list_upt: list
        updated list of images
    """

    # refactoring to use only for sentinel 2
    cloud_property = 'CLOUDY_PIXEL_PERCENTAGE'
    cloud_cover = [_['properties'][cloud_property] for _ in im_list]
    if np.any([_ > prc_cloud_cover for _ in cloud_cover]):
        idx_delete = np.where([_ > prc_cloud_cover for _ in cloud_cover])[0]
        im_list_upt = [x for k,x in enumerate(im_list) if k not in idx_delete]
    else:
        im_list_upt = im_list

    return im_list_upt