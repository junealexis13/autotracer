from datetime import datetime
from shapely import geometry
from scipy import stats, interpolate
from osgeo import gdal
import numpy as np
import pandas as pd
import pyproj
import matplotlib.pyplot as plt
import ee
import os
import pickle
import pytz
import pyproj
import streamlit as st
import requests
import time
import zipfile

def retrieve_images(inputs, filtered_images: list):
    """
    Downloads all images from Sentinel-2 covering the area of interest and acquired between the specified dates.
    The downloaded images are in .TIF format and organised in subfolders, divided
    by satellite mission. The bands are also subdivided by pixel resolution.

    Credits to KV WRL 2018 for this Awesome Algorithm!

    Arguments:
    -----------
    inputs: dict with the following keys
        'sitename': str
            name of the site
        'polygon': list
            polygon containing the lon/lat coordinates to be extracted,
            longitudes in the first column and latitudes in the second column,
            there are 5 pairs of lat/lon with the fifth point equal to the first point:
            ```
            polygon = [[[151.3, -33.7],[151.4, -33.7],[151.4, -33.8],[151.3, -33.8],
            [151.3, -33.7]]]
            ```
        'dates': list of str
            list that contains 2 strings with the initial and final dates in
            format 'yyyy-mm-dd':
            ```
            dates = ['1987-01-01', '2018-01-01']
            ```
        'sat_list': list of str
            list that contains the names of the satellite missions to include:
            ```
            sat_list = ['L5', 'L7', 'L8', 'S2']
            ```
        'filepath_data': str
            filepath to the directory where the images are downloaded
    
        'ee_project_id': str
            Earth Engine project ID, if not specified, the default project will be used

    filtered_images: list of dict
        list of images that were filtered based on the dates and polygon.

    Returns:
    -----------
    metadata: dict
        contains the information about the satellite images that were downloaded:
        date, filename, georeferencing accuracy and image coordinate reference system

    """
    # initialise connection with GEE server
    authenticate_and_initialize(prject_id=inputs['ee_project_id'])


    if 'S2' in inputs['sat_list'] and len(filtered_images)>0:
        im_dict_s2cloudless = get_s2cloudless(filtered_images, inputs)
    
    
    # create a new directory for this site with the name of the site
    im_folder = os.path.join(inputs['filepath'],inputs['sitename'])
    if not os.path.exists(im_folder): os.makedirs(im_folder)
    
    #define the variables / only needed S2 vars
    qa_band_S2 = 'QA60'
    suffix = '.tif'
    band_ids = ['B2','B3','B4','B8','s2cloudless','B11',qa_band_S2]
    
    # Sentinel-2 products don't provide a georeferencing accuracy (RMSE as in Landsat)
    # but they have a flag indicating if the geometric quality control was PASSED or FAILED
    # if passed a value of 1 is stored if failed a value of -1 is stored in the metadata
    # check which flag name is used for the image as it changes for some reason in the archive
    #FLAG NAMES were set of properties in the metadata dict
    # 'GEOMETRIC_QUALITY_FLAG', 'GEOMETRIC_QUALITY', 'quality_check', 'GENERAL_QUALITY','GENERAL_QUALITY_FLAG'

    # download the images and set the folder structure
    st.toast('Downloading images from Earth Engine...')
    filepaths = create_folder_structure(im_folder, 'S2')

    all_names = [] # list for detecting duplicates
    for i, img in enumerate(filtered_images):

        img_date = convert_unix_to_dt(
            img['properties']['system:time_start']).strftime('%Y-%m-%d-%H-%M-%S')
        
        # check the accuracy
        flag_names = ['GEOMETRIC_QUALITY_FLAG', 'GEOMETRIC_QUALITY', 'quality_check', 'GENERAL_QUALITY','GENERAL_QUALITY_FLAG']
        for key in flag_names: 
            if key in img['properties'].keys(): 
                break # use the first flag that is found
        if len(key) > 0:
            acc_georef = img['properties'][key]
        else:
            st.warning('WARNING: could not find Sentinel-2 geometric quality flag,'+ 
                        ' raise an issue at https://github.com/kvos/CoastSat/issues'+
                        ' and add you inputs in text (not a screenshot pls).')
            acc_georef = 'PASSED'
        
        # add the radiometric image quality ['PASSED' or 'FAILED']
        flag_names = ['RADIOMETRIC_QUALITY', 'RADIOMETRIC_QUALITY_FLAG']
        key = []
        for key in flag_names: 
            if key in img['properties'].keys(): 
                break # use the first flag that is found
        if len(key) > 0:
            rad_quality = img['properties'][key]
        else:
            print('WARNING: could not find Sentinel-2 geometric quality flag,'+ 
                    ' raise an issue at https://github.com/kvos/CoastSat/issues'+
                    ' and add your inputs in text (not a screenshot pls).')
            rad_quality = 'PASSED'

        tilename = img['properties']['MGRS_TILE']
        im_epsg = int(img['bands'][0]['crs'][5:])


        #select the image by ID
        image_EE = ee.Image(img['id'])

        #for cloudless query
        if len(im_dict_s2cloudless[i]) == 0:
            print('Warning: S2cloudless mask for image %s is not available yet, try again tomorrow.'%img_date)
            continue
        im_cloud = ee.Image(im_dict_s2cloudless[i]['id'])
        cloud_prob = im_cloud.select('probability').rename('s2cloudless')
        image_ee = image_ee.addBands(cloud_prob)

        # download the images as .tif files
        bands = dict([])
        im_fn = dict([])
        # first delete dimensions key from dictionnary
        # otherwise the entire image is extracted (don't know why)
        im_bands = image_ee.getInfo()['bands']
        for j in range(len(im_bands)):
            if 'dimensions' in im_bands[j].keys():
                del im_bands[j]['dimensions']

        fp_ms = filepaths[1]
        fp_swir = filepaths[2]
        fp_mask = filepaths[3]    
        # select bands (10m ms RGB+NIR+s2cloudless, 20m SWIR1, 60m QA band)
        bands['ms'] = [im_bands[_] for _ in range(len(im_bands)) if im_bands[_]['id'] in band_ids[:5]]
        bands['swir'] = [im_bands[_] for _ in range(len(im_bands)) if im_bands[_]['id'] in band_ids[5:6]]
        bands['mask'] = [im_bands[_] for _ in range(len(im_bands)) if im_bands[_]['id'] in band_ids[-1:]]
        # adjust polygon for both ms and pan bands
        proj_ms = image_ee.select('B1').projection()
        proj_swir = image_ee.select('B11').projection()
        proj_mask = image_ee.select('QA60').projection()
        ee_region_ms = adjust_polygon(inputs['polygon'],proj_ms)
        ee_region_swir = adjust_polygon(inputs['polygon'],proj_swir)
        ee_region_mask = adjust_polygon(inputs['polygon'],proj_mask)
        # download the ms, swir and QA bands from EE
        count = 0
        while True:
            try:    
                fn_ms = download_tif(image_ee,ee_region_ms,bands['ms'],fp_ms,'S2')
                fn_swir = download_tif(image_ee,ee_region_swir,bands['swir'],fp_swir,'S2')
                fn_QA = download_tif(image_ee,ee_region_mask,bands['mask'],fp_mask,'S2')
                break
            except Exception as e:
                print('\nDownload failed, trying again...')
                print(f"The error that occurred: {e}")
                time.sleep(20)
                count += 1
                if count > 100:
                    raise Exception('Too many attempts, crashed while downloading image %s'%img['id'])
                else:
                    continue       

         # create filename for the three images (ms, swir and mask)
        for key in bands.keys():
            im_fn[key] = img_date + '_' + "S2" + '_' + tilename + '_' + inputs['sitename'] + '_' + key + suffix
        # if multiple images taken at the same date add 'dupX' to the name (duplicate)
        duplicate_counter = 0
        while im_fn['ms'] in all_names:
            duplicate_counter += 1
            for key in bands.keys():
                im_fn[key] = img_date + '_' + "S2" + '_' + tilename + '_' \
                    + inputs['sitename'] + '_' + key \
                    + '_dup%d'%duplicate_counter + suffix
        filename_ms = im_fn['ms']
        all_names.append(im_fn['ms']) 
        
        # resample the 20m swir band to the 10m ms band with bilinear interpolation
        fn_in = fn_swir
        fn_target = fn_ms
        fn_out = os.path.join(fp_swir, im_fn['swir'])
        warp_image_to_target(fn_in,fn_out,fn_target,double_res=False,resampling_method='bilinear')             
        
        # resample 60m QA band to the 10m ms band with nearest-neighbour interpolation
        fn_in = fn_QA
        fn_target = fn_ms
        fn_out = os.path.join(fp_mask, im_fn['mask'])
        warp_image_to_target(fn_in,fn_out,fn_target,double_res=False,resampling_method='near')
        
        # delete original downloads
        for _ in [fn_swir,fn_QA]: os.remove(_)  
        # rename the multispectral band file
        os.rename(fn_ms,os.path.join(fp_ms, im_fn['ms']))
            
        # get image dimensions (width and height)
        image_path = os.path.join(fp_ms,im_fn['ms'])
        width, height = get_image_dimensions(image_path)
        # write metadata in a text file for easy access
        filename_txt = im_fn['ms'].replace('_ms','').replace('.tif','')
        metadict = {'filename':filename_ms,'tile':tilename,'epsg':im_epsg,
                    'acc_georef':acc_georef,'image_quality':rad_quality,
                    'im_width':width,'im_height':height}
        with open(os.path.join(filepaths[0],filename_txt + '.txt'), 'w') as f:
            for key in metadict.keys():
                f.write('%s\t%s\n'%(key,metadict[key]))
        # print percentage completion for user
        print('\r%d%%' %int((i+1)/len(img)*100), end='')


    # once all images have been downloaded, load metadata from .txt files
    metadata = get_metadata(inputs)
    
    # merge overlapping images (necessary only if the polygon is at the boundary of an image)
    # if 'S2' in metadata.keys():
    #     print("\n Called merge_overlapping_images\n")
    #     try:
    #         metadata = merge_overlapping_images(metadata,inputs)
    #     except:
    #         print('WARNING: there was an error while merging overlapping S2 images,'+
    #               ' please open an issue on Github at https://github.com/kvos/CoastSat/issues'+
    #               ' and include your script so we can find out what happened.')

    # save metadata dict
    with open(os.path.join(im_folder, inputs['sitename'] + '_metadata' + '.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
    print('Satellite images downloaded from GEE and save in %s'%im_folder)
    return metadata

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

def get_s2cloudless(im_list, inputs):
    "Match the list of S2 images with the corresponding s2cloudless images"
    # get s2cloudless collection
    dates = [datetime.strptime(_,'%Y-%m-%d') for _ in inputs['dates']]
    polygon = inputs['polygon']
    collection = 'COPERNICUS/S2_CLOUD_PROBABILITY'
    s2cloudless_col = ee.ImageCollection(collection).filterBounds(ee.Geometry.Polygon(polygon))\
                                                    .filterDate(dates[0],dates[1])
    im_list_cloud = s2cloudless_col.getInfo().get('features')
    # get image ids
    indices_cloud = [_['properties']['system:index'] for _ in im_list_cloud]
    # match with S2 images
    im_list_cloud_matched = []
    for i in range(len(im_list)):
        index = im_list[i]['properties']['system:index'] 
        if index in indices_cloud:
            k = np.where([_ == index for _ in indices_cloud])[0][0]
            im_list_cloud_matched.append(im_list_cloud[k])
        else: # put an empty list if no match
            im_list_cloud_matched.append([])
    return im_list_cloud_matched

def adjust_polygon(polygon,proj):
    """
    Adjust polygon of ROI to fit exactly with the pixels of the underlying tile

    KV WRL 2022

    Arguments:
    -----------
    polygon: list
        polygon containing the lon/lat coordinates to be extracted,
        longitudes in the first column and latitudes in the second column,
        there are 5 pairs of lat/lon with the fifth point equal to the first point:
        ```
        polygon = [[[151.3, -33.7],[151.4, -33.7],[151.4, -33.8],[151.3, -33.8],
        [151.3, -33.7]]]
        ```
    proj: ee.Proj
        projection of the underlying tile

    Returns:
    -----------
    ee_region: ee
        updated list of images
    """    
    # adjust polygon to match image coordinates so that there is no resampling
    polygon_ee = ee.Geometry.Polygon(polygon)    
    # convert polygon to image coordinates
    polygon_coords = np.array(ee.List(polygon_ee.transform(proj, 1).coordinates().get(0)).getInfo())
    # make it a rectangle
    xmin = np.min(polygon_coords[:,0])
    ymin = np.min(polygon_coords[:,1])
    xmax = np.max(polygon_coords[:,0])
    ymax = np.max(polygon_coords[:,1])
    # round to the closest pixels
    rect = [np.floor(xmin), np.floor(ymin), 
            np.ceil(xmax),  np.ceil(ymax)]
    # convert back to epsg 4326
    ee_region = ee.Geometry.Rectangle(rect, proj, True, False).transform("EPSG:4326")
    
    return ee_region

def download_tif(image, polygon, bands, filepath, satname):
    """
    Downloads a .TIF image from the ee server. The image is downloaded as a
    zip file then moved to the working directory, unzipped and stacked into a
    single .TIF file. Any QA band is saved separately.

    KV WRL 2018

    Arguments:
    -----------
    image: ee.Image
        Image object to be downloaded
    polygon: list
        polygon containing the lon/lat coordinates to be extracted
        longitudes in the first column and latitudes in the second column
    bands: list of dict
        list of bands to be downloaded
    filepath: str
        location where the temporary file should be saved
    satname: str
        name of the satellite missions ['L5','L7','L8','S2']
    Returns:
    -----------
    Downloads an image in a file named data.tif

    """
     
    # crop and download
    download_id = ee.data.getDownloadId({'image': image,
                                         'region': polygon,
                                         'bands': bands,
                                         'filePerBand': True,
                                         'name': 'image'})
    response = requests.get(ee.data.makeDownloadUrl(download_id))  
    fp_zip = os.path.join(filepath,'temp.zip')
    with open(fp_zip, 'wb') as fd:
      fd.write(response.content) 
    # unzip the individual bands
    with zipfile.ZipFile(fp_zip) as local_zipfile:
        for fn in local_zipfile.namelist():
            local_zipfile.extract(fn, filepath)
        fn_all = [os.path.join(filepath,_) for _ in local_zipfile.namelist()]
    os.remove(fp_zip)
    # now process the individual bands:
    # - for Landsat
    if satname in ['L5','L7','L8','L9']:
        # if there is only one band, it's the panchromatic
        if len(fn_all) == 1:
            # return the filename of the .tif
            return fn_all[0]
        # otherwise there are multiple multispectral bands so we have to merge them into one .tif
        else:
            # select all ms bands except the QA band (which is processed separately)
            fn_tifs = [_ for _ in fn_all if not 'QA' in _]
            filename = 'ms_bands.tif'
            # build a VRT and merge the bands (works the same with pan band)
            outds = gdal.BuildVRT(os.path.join(filepath,'temp.vrt'),
                                  fn_tifs, separate=True)
            outds = gdal.Translate(os.path.join(filepath,filename), outds) 
            # remove temporary files
            os.remove(os.path.join(filepath,'temp.vrt'))
            for _ in fn_tifs: os.remove(_)
            if os.path.exists(os.path.join(filepath,filename+'.aux.xml')):
                os.remove(os.path.join(filepath,filename+'.aux.xml'))
            # return file names (ms and QA bands separately)
            fn_image = os.path.join(filepath,filename)
            fn_QA = [_ for _ in fn_all if 'QA' in _][0]
            return fn_image, fn_QA
    # - for Sentinel-2
    if satname in ['S2']:
        # if there is only one band, it's either the SWIR1 or QA60
        if len(fn_all) == 1:
            # return the filename of the .tif
            return fn_all[0]
        # otherwise there are multiple multispectral bands so we have to merge them into one .tif
        else:
            # select all ms bands except the QA band (which is processed separately)
            fn_tifs = fn_all
            filename = 'ms_bands.tif'
            # build a VRT and merge the bands (works the same with pan band)
            outds = gdal.BuildVRT(os.path.join(filepath,'temp.vrt'),
                                  fn_tifs, separate=True)
            outds = gdal.Translate(os.path.join(filepath,filename), outds) 
            # remove temporary files
            os.remove(os.path.join(filepath,'temp.vrt'))
            for _ in fn_tifs: os.remove(_)
            if os.path.exists(os.path.join(filepath,filename+'.aux.xml')):
                os.remove(os.path.join(filepath,filename+'.aux.xml'))
            # return filename of the merge .tif file
            fn_image = os.path.join(filepath,filename)
            return fn_image           
        
def warp_image_to_target(fn_in,fn_out,fn_target,double_res=True,resampling_method='bilinear'):
    """
    Resample an image on a new pixel grid based on a target image using gdal_warp.
    This is used to align the multispectral and panchromatic bands, as well as just downsample certain bands.

    KV WRL 2022

    Arguments:
    -----------
    fn_in: str
        filepath of the input image (points to .tif file)
    fn_out: str
        filepath of the output image (will be created)
    fn_target: str
        filepath of the target image
    double_res: boolean
        this function can be used to downsample images by settings the input and target 
        filepaths to the same imageif the input and target images are the same and settings
        double_res = True to downsample by a factor of 2
    resampling_method: str
        method using to resample the image on the new pixel grid. See gdal_warp documentation
        for options (https://gdal.org/programs/gdalwarp.html)

    Returns:
    -----------
    Creates a new .tif file (fn_out)

    """    
    # get output extent from target image
    im_target = gdal.Open(fn_target, gdal.GA_ReadOnly)
    georef_target = np.array(im_target.GetGeoTransform())
    xres =georef_target[1]
    yres = georef_target[5]
    if double_res:
        xres = int(georef_target[1]/2)
        yres = int(georef_target[5]/2)      
    extent_pan = get_image_bounds(fn_target)
    extent_coords = np.array(extent_pan.exterior.coords)
    xmin = np.min(extent_coords[:,0])
    ymin = np.min(extent_coords[:,1])
    xmax = np.max(extent_coords[:,0])
    ymax = np.max(extent_coords[:,1])
    
    # use gdal_warp to resample the inputon the target image pixel grid
    options = gdal.WarpOptions(xRes=xres, yRes=yres,
                               outputBounds=[xmin, ymin, xmax, ymax],
                               resampleAlg=resampling_method,
                               targetAlignedPixels=False)
    gdal.Warp(fn_out, fn_in, options=options)
    
    # check that both files have the same georef and size (important!)
    im_target = gdal.Open(fn_target, gdal.GA_ReadOnly)
    im_out = gdal.Open(fn_out, gdal.GA_ReadOnly)
    georef_target = np.array(im_target.GetGeoTransform())
    georef_out = np.array(im_out.GetGeoTransform())
    size_target = np.array([im_target.RasterXSize,im_target.RasterYSize])
    size_out = np.array([im_out.RasterXSize,im_out.RasterYSize])
    if double_res: size_target = size_target*2
    if np.any(np.nonzero(georef_target[[0,3]]-georef_out[[0,3]])): 
        raise Exception('Georef of pan and ms bands do not match for image %s'%fn_out)
    if np.any(np.nonzero(size_target-size_out)): 
        raise Exception('Size of pan and ms bands do not match for image %s'%fn_out)

def get_image_bounds(fn):
    """
    Returns a polygon with the bounds of the image in the .tif file
     
    KV WRL 2020

    Arguments:
    -----------
    fn: str
        path to the image (.tif file)         
                
    Returns:    
    -----------
    bounds_polygon: shapely.geometry.Polygon
        polygon with the image bounds
        
    """
    
    # nested functions to get the extent 
    # copied from https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
    def GetExtent(gt,cols,rows):
        'Return list of corner coordinates from a geotransform'
        ext=[]
        xarr=[0,cols]
        yarr=[0,rows]
        for px in xarr:
            for py in yarr:
                x=gt[0]+(px*gt[1])+(py*gt[2])
                y=gt[3]+(px*gt[4])+(py*gt[5])
                ext.append([x,y])
            yarr.reverse()
        return ext
    
    # load .tif file and get bounds
    if not os.path.exists(fn):
        raise FileNotFoundError(f"{fn}")
    data = gdal.Open(fn, gdal.GA_ReadOnly)
    # Check if data is null meaning the open failed
    if data is None:
        print("TIF file: ",fn, "cannot be opened" )
        os.remove(fn)
        raise AttributeError
    else:
        gt = data.GetGeoTransform()
        cols = data.RasterXSize
        rows = data.RasterYSize
        ext = GetExtent(gt,cols,rows)
    
    return geometry.Polygon(ext)

def get_image_bounds(fn):
    """
    Returns a polygon with the bounds of the image in the .tif file
     
    KV WRL 2020

    Arguments:
    -----------
    fn: str
        path to the image (.tif file)         
                
    Returns:    
    -----------
    bounds_polygon: shapely.geometry.Polygon
        polygon with the image bounds
        
    """
    
    # nested functions to get the extent 
    # copied from https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
    def GetExtent(gt,cols,rows):
        'Return list of corner coordinates from a geotransform'
        ext=[]
        xarr=[0,cols]
        yarr=[0,rows]
        for px in xarr:
            for py in yarr:
                x=gt[0]+(px*gt[1])+(py*gt[2])
                y=gt[3]+(px*gt[4])+(py*gt[5])
                ext.append([x,y])
            yarr.reverse()
        return ext
    
    # load .tif file and get bounds
    if not os.path.exists(fn):
        raise FileNotFoundError(f"{fn}")
    data = gdal.Open(fn, gdal.GA_ReadOnly)
    # Check if data is null meaning the open failed
    if data is None:
        print("TIF file: ",fn, "cannot be opened" )
        os.remove(fn)
        raise AttributeError
    else:
        gt = data.GetGeoTransform()
        cols = data.RasterXSize
        rows = data.RasterYSize
        ext = GetExtent(gt,cols,rows)
    
    return geometry.Polygon(ext)

def get_image_dimensions(image_path):
    "function to get image dimensions with GDAL"
    dataset = gdal.Open(image_path, gdal.GA_ReadOnly)
    if dataset is None:
        raise Exception("Failed to open the image file %s"%image_path)
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    dataset = None

    return width, height


def convert_unix_to_dt(unix_ts):
    """
    Convert a Unix timestamp in milliseconds to a UTC datetime object.
    
    :param unix_ts: Unix timestamp in milliseconds
    :return: UTC datetime object
    """
    # convert ms to s
    timestamp_s = unix_ts / 1000

    # conv to datetime (UTC)
    dt = datetime.fromtimestamp(timestamp_s, tz=pytz.utc)
    return dt

def create_folder_structure(im_folder, satname):
    """
    Create the structure of subfolders for each satellite mission

    KV WRL 2018

    Arguments:
    -----------
    im_folder: str
        folder where the images are to be downloaded
    satname:
        name of the satellite mission

    Returns:
    -----------
    filepaths: list of str
        filepaths of the folders that were created
    """

    # one folder for the metadata (common to all satellites)
    filepaths = [os.path.join(im_folder, satname, 'meta')]
    # subfolders depending on satellite mission

    filepaths.append(os.path.join(im_folder, satname, 'ms'))
    filepaths.append(os.path.join(im_folder, satname, 'swir'))
    filepaths.append(os.path.join(im_folder, satname, 'mask'))

    for fp in filepaths:
        if not os.path.exists(fp): os.makedirs(fp)

    return filepaths