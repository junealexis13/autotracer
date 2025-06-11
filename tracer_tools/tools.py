from datetime import datetime, timedelta
from shapely import geometry
from scipy import stats, interpolate
from osgeo import gdal
from osgeo import gdal_array
import skimage.morphology as morphology
import skimage.exposure as exposure
import skimage.transform as transform
from skimage import img_as_ubyte
from skimage.io import imsave
from astropy.convolution import convolve
from dateutil.relativedelta import relativedelta
import geopandas as gpd
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
import re

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
        image_ee= ee.Image(img['id'])

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
                time.sleep(10)
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
        st.toast('Status: \r%d%%' %int((i+1)/len(img)*100))


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
    st.toast('Satellite images downloaded from GEE and save in %s'%im_folder)
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
        st.toast('GEE initialized (existing token).')
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
            st.toast('GEE initialized (refreshed token).')
        except:
            # get the user to authenticate manually and initialize the sesion
            ee.Authenticate()
            ee.Initialize(project=prject_id)
            st.toast('GEE initialized (manual authentication).')

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
        
    st.toast('Number of images available between %s and %s:'%(dates_str[0],dates_str[1]))
    
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
    st.toast('     %s: %d images'%('S2',len(im_list)))
    im_dict_T1['S2'] = im_list          
        
    st.toast('  Total to download: %d images'%sum_img)

    # check if images already exist  
    # st.toast('\nLooking for existing imagery...')
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
                        st.toast(f'{satname}:There are {len(avail_date_list)} images available, {len(metadata_existing[satname]["dates"])} images already exist, {len(avail_date_list)} to download')
                        continue
                    # get the dates of the images that are already downloaded
                    downloaded_dates = metadata_existing[satname]['dates']
                    # if no images are already downloaded, skip this loop and use whats already in im_dict_T1[satname]
                    if len(downloaded_dates) == 0:
                        st.toast(f'{satname}:There are {len(avail_date_list)} images available, {len(downloaded_dates)} images already exist, {len(avail_date_list)} to download')
                        continue
                    # get the indices of the images that are not already downloaded 
                    idx_new = np.where([ not avail_date in downloaded_dates for avail_date in avail_date_list])[0]
                    im_dict_T1[satname] = [im_dict_T1[satname][index] for index in idx_new]
                    st.toast('%s: %d/%d images already exist, %s to be downloaded'%(satname,
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
    st.toast('- In Landsat Tier 2 (not suitable for time-series analysis):')
    im_dict_T2 = dict([])
    sum_img = 0
    for satname in inputs['sat_list']:
        if satname in ['L9','S2']: continue # no Tier 2 for Sentinel-2 and Landsat 9
        im_list = get_image_info(col_names_T2[satname],polygon,dates_str)
        sum_img = sum_img + len(im_list)
        st.toast('     %s: %d images'%(satname,len(im_list)))
        im_dict_T2[satname] = im_list

    st.toast('  Total Tier 2: %d images'%sum_img)
    
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
                                         'name': 'image',
                                         'scale': 20 })
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
        st.toast("TIF file: ",fn, "cannot be opened" )
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
        st.toast("TIF file: ",fn, "cannot be opened" )
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

def get_filepath(inputs,satname):
    """
    Create filepath to the different folders containing the satellite images.
    
    KV WRL 2018

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
        'sat_list': ["S2"]
            Only set to S2 - not changeable
            ```
            sat_list = ['L5', 'L7', 'L8', 'L9', 'S2']
            ```
        'filepath': str
            filepath to the directory where the images are downloaded
    satname: str
        short name of the satellite mission ('L5','L7','L8','S2')
        NOTE: This was removed and set only to S2
                
    Returns:    
    -----------
    filepath: str or list of str
        contains the filepath(s) to the folder(s) containing the satellite images
    
    """     
    
    sitename = inputs['sitename']
    filepath_data = inputs['filepath']
    # access the images

    # access downloaded Sentinel 2 images
    fp_ms = os.path.join(filepath_data, sitename, satname, 'ms')
    fp_swir = os.path.join(filepath_data, sitename, satname, 'swir')
    fp_mask = os.path.join(filepath_data, sitename, satname, 'mask')
    filepath = [fp_ms, fp_swir, fp_mask]
        
    return filepath

def get_filenames(filename, filepath, satname='S2'):
    """
    Creates filepath + filename for all the bands belonging to the same image.
    
    KV WRL 2018

    Arguments:
    -----------
    filename: str
        name of the downloaded satellite image as found in the metadata
    filepath: str or list of str
        contains the filepath(s) to the folder(s) containing the satellite images
    satname: str
        SET TO S2 only - CoastSat enables the user to choose between it. But are disabled.       
        
    Returns:    
    -----------
    fn: str or list of str
        contains the filepath + filenames to access the satellite image
        
    """ 
    fn_swir = filename.replace('_ms','_swir')
    fn_mask = filename.replace('_ms','_mask')
    fn = [os.path.join(filepath[0], filename),
            os.path.join(filepath[1], fn_swir),
            os.path.join(filepath[2], fn_mask)]
        
    return fn

def create_cloud_mask(im_QA, satname, cloud_mask_issue):
    """
    Creates a cloud mask using the information contained in the QA band.

    KV WRL 2018

    Arguments:
    -----------
    im_QA: np.array
        Image containing the QA band
    satname: string
        short name for the satellite: ```'L5', 'L7', 'L8' or 'S2'```
    cloud_mask_issue: boolean
        True if there is an issue with the cloud mask and sand pixels are being
        erroneously masked on the images
        
    Returns:
    -----------
    cloud_mask : np.array
        boolean array with True if a pixel is cloudy and False otherwise

    """

    # 1024 = dense cloud, 2048 = cirrus clouds
    cloud_values = [1024, 2048] 
    cloud_mask = np.isin(im_QA, cloud_values)
    if sum(sum(cloud_mask)) > 0 and sum(sum(~cloud_mask)) > 0:
        cloud_mask = morphology.remove_small_objects(cloud_mask, min_size=40, connectivity=1)

    if cloud_mask_issue:
        cloud_mask = np.zeros_like(im_QA, dtype=bool)
        for value in cloud_values:
            cloud_mask_temp = np.isin(im_QA, value)         
            elem = morphology.footprint_rectangle(6) # use a square of width 6 pixels
            cloud_mask_temp = morphology.binary_opening(cloud_mask_temp, elem) # perform image opening            
            cloud_mask_temp = morphology.remove_small_objects(cloud_mask_temp, min_size=100, connectivity=1)
            cloud_mask = np.logical_or(cloud_mask, cloud_mask_temp)

    return cloud_mask

def create_s2cloudless_mask(cloud_prob, s2cloudless_prob):
    """
    Creates a cloud mask using the s2cloudless band.

    KV WRL 2023

    Arguments:
    -----------
    cloud_prob: np.array
        Image containing the s2cloudless cloud probability
        
    Returns:
    -----------
    cloud_mask : np.array
        boolean array with True if a pixel is cloudy and False otherwise

    """
    # find which pixels have bits corresponding to cloud values
    cloud_mask = cloud_prob > s2cloudless_prob
    # dilate cloud mask
    elem = morphology.footprint_rectangle(6) # use a square of width 6 pixels
    cloud_mask = morphology.binary_opening(cloud_mask,elem) # perform image opening

    return cloud_mask

def preprocess_single(fn, satname, cloud_mask_issue, pan_off, s2cloudless_prob=40):
    """
    Reads the image and outputs the pansharpened/down-sampled multispectral bands,
    the georeferencing vector of the image (coordinates of the upper left pixel),
    the cloud mask, the QA band and a no_data image.
    For Landsat 7-8 it also outputs the panchromatic band and for Sentinel-2 it
    also outputs the 20m SWIR band.

    KV WRL 2018

    Arguments:
    -----------
    fn: str or list of str
        filename of the .TIF file containing the image. For L7, L8 and S2 this
        is a list of filenames, one filename for each band at different
        resolution (30m and 15m for Landsat 7-8, 10m, 20m, 60m for Sentinel-2)
    satname: str
        name of the satellite mission (e.g., 'L5')
    cloud_mask_issue: boolean
        True if there is an issue with the cloud mask and sand pixels are being masked on the images
    pan_off : boolean
        if True, disable panchromatic sharpening and ignore pan band
    s2cloudless_prob: float [0,100)
        threshold to identify cloud pixels in the s2cloudless probability mask
        
    Returns:
    -----------
    im_ms: np.array
        3D array containing the pansharpened/down-sampled bands (B,G,R,NIR,SWIR1)
    georef: np.array
        vector of 6 elements [Xtr, Xscale, Xshear, Ytr, Yshear, Yscale] defining the
        coordinates of the top-left pixel of the image
    cloud_mask: np.array
        2D cloud mask with True where cloud pixels are
    im_extra : np.array
        2D array containing the 20m resolution SWIR band for Sentinel-2 and the 15m resolution
        panchromatic band for Landsat 7 and Landsat 8. This field is empty for Landsat 5.
    im_QA: np.array
        2D array containing the QA band, from which the cloud_mask can be computed.
    im_nodata: np.array
        2D array with True where no data values (-inf) are located

    """
    
    if isinstance(fn, list):
        fn_to_split = fn[0]
    elif isinstance(fn, str):
        fn_to_split = fn
    # split by os.sep and only get the filename at the end then split again to remove file extension
    fn_to_split=fn_to_split.split(os.sep)[-1].split('.')[0]
    # search for the year the tif was taken with regex and convert to int
    year = int(re.search('[0-9]+',fn_to_split).group(0))

    # read 10m bands (R,G,B,NIR)
    fn_ms = fn[0]
    data = gdal.Open(fn_ms, gdal.GA_ReadOnly)
    georef = np.array(data.GetGeoTransform())
    bands = [data.GetRasterBand(k + 1).ReadAsArray() for k in range(data.RasterCount-1)]
    im_ms = np.stack(bands, 2)
    im_ms = im_ms/10000 # TOA scaled to 10000
    # read s2cloudless cloud probability (last band in ms image)
    cloud_prob = data.GetRasterBand(data.RasterCount).ReadAsArray()

    # image size
    nrows = im_ms.shape[0]
    ncols = im_ms.shape[1]
    # if image contains only zeros (can happen with S2), skip the image
    if sum(sum(sum(im_ms))) < 1:
        im_ms = []
        georef = []
        # skip the image by giving it a full cloud_mask
        cloud_mask = np.ones((nrows,ncols)).astype('bool')
        return im_ms, georef, cloud_mask, [], [], []

    # read 20m band (SWIR1)
    fn_swir = fn[1]
    data = gdal.Open(fn_swir, gdal.GA_ReadOnly)
    bands = [data.GetRasterBand(k + 1).ReadAsArray() for k in range(data.RasterCount)]
    im_swir = bands[0]
    im_swir = im_swir/10000 # TOA scaled to 10000
    im_swir = np.expand_dims(im_swir, axis=2)

    # append down-sampled SWIR1 band to the other 10m bands
    im_ms = np.append(im_ms, im_swir, axis=2)

    # create cloud mask using 60m QA band (not as good as Landsat cloud cover)
    fn_mask = fn[2]
    data = gdal.Open(fn_mask, gdal.GA_ReadOnly)
    bands = [data.GetRasterBand(k + 1).ReadAsArray() for k in range(data.RasterCount)]
    im_QA = bands[0]
    # compute cloud mask using QA60 band
    cloud_mask_QA60 = create_cloud_mask(im_QA, satname, cloud_mask_issue)
    # compute cloud mask using s2cloudless probability band
    cloud_mask_s2cloudless = create_s2cloudless_mask(cloud_prob, s2cloudless_prob)
    # combine both cloud masks
    cloud_mask = np.logical_or(cloud_mask_QA60,cloud_mask_s2cloudless)
    
    # check if -inf or nan values on any band and create nodata image
    im_nodata = np.zeros(cloud_mask.shape).astype(bool)
    for k in range(im_ms.shape[2]):
        im_inf = np.isin(im_ms[:,:,k], -np.inf)
        im_nan = np.isnan(im_ms[:,:,k])
        im_nodata = np.logical_or(np.logical_or(im_nodata, im_inf), im_nan)
    # add the edges of the SWIR1 band that contains only 0's to the nodata image
    # these are created when reprojecting the SWIR1 20 m band onto the 10m pixel grid
    im_nodata = pad_edges(im_swir, im_nodata)        
    # check if there are pixels with 0 intensity in the Green, NIR and SWIR bands and add those
    # to the cloud mask as otherwise they will cause errors when calculating the NDWI and MNDWI
    im_zeros = np.ones(im_nodata.shape).astype(bool)
    im_zeros = np.logical_and(np.isin(im_ms[:,:,1],0), im_zeros) # Green
    im_zeros = np.logical_and(np.isin(im_ms[:,:,3],0), im_zeros) # NIR
    im_zeros = np.logical_and(np.isin(im_ms[:,:,4],0), im_zeros) # SWIR
    # add to im_nodata
    im_nodata = np.logical_or(im_zeros, im_nodata)
    # dilate if image was merged as there could be issues at the edges
    if 'merged' in fn_ms:
        im_nodata = morphology.dilation(im_nodata,morphology.footprint_rectangle(5))

    # update cloud mask with all the nodata pixels
    cloud_mask = np.logical_or(cloud_mask, im_nodata)

    # no extra image
    im_extra = []

    return im_ms, georef, cloud_mask, im_extra, im_QA, im_nodata

def save_jpg(metadata, settings: dict, use_matplotlib=False):
    """
    Saves a .jpg image for all the images contained in metadata.

    KV WRL 2018

    Arguments:
    -----------
    metadata: dict
        contains all the information about the satellite images that were downloaded
    settings: dict with the following keys
        'inputs': dict
            input parameters (sitename, filepath, polygon, dates, sat_list)
        'cloud_thresh': float
            value between 0 and 1 indicating the maximum cloud fraction in
            the cropped image that is accepted
        'cloud_mask_issue': boolean
            True if there is an issue with the cloud mask and sand pixels
            are erroneously being masked on the images
        's2cloudless_prob': float [0,100)
            threshold to identify cloud pixels in the s2cloudless probability mask
        'use_matplotlib': boolean
            False to save a .jpg and True to save as matplotlib plots

    Returns:
    -----------
    Stores the images as .jpg in a folder named /preprocessed

    """

    sitename = settings['inputs']['sitename']
    cloud_thresh = settings['cloud_thresh']
    s2cloudless_prob = settings['s2cloudless_prob']
    filepath_data = settings['inputs']['filepath']
    
    # create subfolder to store the jpg files
    filepath_jpg = os.path.join(filepath_data, sitename, 'jpg_files', 'preprocessed')
    if not os.path.exists(filepath_jpg):
            os.makedirs(filepath_jpg)

    # loop through satellite list
    st.toast('Saving images as jpg:')
    for satname in metadata.keys():
        
        filepath = get_filepath(settings['inputs'],satname)
        filenames = metadata[satname]['filenames']
        st.toast('%s: %d images'%(satname,len(filenames)))
        # loop through images
        for i in range(len(filenames)):
            st.toast('Status:  \r%d%%' %int((i+1)/len(filenames)*100))
            # image filename
            fn = get_filenames(filenames[i],filepath, satname)
            # read and preprocess image
            im_ms, georef, cloud_mask, im_extra, im_QA, im_nodata = preprocess_single(fn, satname, settings['cloud_mask_issue'],
                                                                                      settings['pan_off'], s2cloudless_prob)

            # compute cloud_cover percentage (with no data pixels)
            cloud_cover_combined = np.divide(sum(sum(cloud_mask.astype(int))),
                                    (cloud_mask.shape[0]*cloud_mask.shape[1]))
            if cloud_cover_combined > 0.99: # if 99% of cloudy pixels in image skip
                continue

            # remove no data pixels from the cloud mask (for example L7 bands of no data should not be accounted for)
            cloud_mask_adv = np.logical_xor(cloud_mask, im_nodata)
            # compute updated cloud cover percentage (without no data pixels)
            cloud_cover = np.divide(sum(sum(cloud_mask_adv.astype(int))),
                                    (sum(sum((~im_nodata).astype(int)))))
            # skip image if cloud cover is above threshold
            if cloud_cover > cloud_thresh or cloud_cover == 1:
                continue
            # save .jpg with date and satellite in the title
            date = filenames[i][:19]
            plt.ioff()  # turning interactive plotting off
            create_jpg(im_ms, cloud_mask, date, satname, filepath_jpg, use_matplotlib)
    # print the location where the images have been saved
    st.toast('Satellite images saved as .jpg in ' + os.path.join(filepath_data, sitename,
                                                    'jpg_files', 'preprocessed'))

def rescale_image_intensity(im, cloud_mask, prob_high):
    """
    Rescales the intensity of an image (multispectral or single band) by applying
    a cloud mask and clipping the prob_high upper percentile. This functions allows
    to stretch the contrast of an image, only for visualisation purposes.

    KV WRL 2018

    Arguments:
    -----------
    im: np.array
        Image to rescale, can be 3D (multispectral) or 2D (single band)
    cloud_mask: np.array
        2D cloud mask with True where cloud pixels are
    prob_high: float
        probability of exceedence used to calculate the upper percentile

    Returns:
    -----------
    im_adj: np.array
        rescaled image
    """

    # lower percentile is set to 0
    prc_low = 0

    # reshape the 2D cloud mask into a 1D vector
    vec_mask = cloud_mask.reshape(im.shape[0] * im.shape[1])

    # if image contains several bands, stretch the contrast for each band
    if len(im.shape) > 2:
        # reshape into a vector
        vec =  im.reshape(im.shape[0] * im.shape[1], im.shape[2])
        # initiliase with NaN values
        vec_adj = np.ones((len(vec_mask), im.shape[2])) * np.nan
        # loop through the bands
        for i in range(im.shape[2]):
            # find the higher percentile (based on prob)
            prc_high = np.percentile(vec[~vec_mask, i], prob_high)
            # clip the image around the 2 percentiles and rescale the contrast
            vec_rescaled = exposure.rescale_intensity(vec[~vec_mask, i],
                                                      in_range=(prc_low, prc_high))
            vec_adj[~vec_mask,i] = vec_rescaled
        # reshape into image
        im_adj = vec_adj.reshape(im.shape[0], im.shape[1], im.shape[2])

    # if image only has 1 bands (grayscale image)
    else:
        vec =  im.reshape(im.shape[0] * im.shape[1])
        vec_adj = np.ones(len(vec_mask)) * np.nan
        prc_high = np.percentile(vec[~vec_mask], prob_high)
        vec_rescaled = exposure.rescale_intensity(vec[~vec_mask], in_range=(prc_low, prc_high))
        vec_adj[~vec_mask] = vec_rescaled
        im_adj = vec_adj.reshape(im.shape[0], im.shape[1])

    return im_adj

def pad_edges(im_swir: np.ndarray, im_nodata: np.ndarray) -> np.ndarray:
    """
    Adds 0's located along the edges of im_swir to the nodata array.

    Fixes the issue where 0s are created along the edges of the SWIR1 band caused by reprojecting the 20 m band onto the 10m pixel grid (with bilinear interpolation in GDAL)

    Args:
        im_swir (np.ndarray): The SWIR image.
        im_nodata (np.ndarray): The nodata array.

    Returns:
        np.ndarray: The nodata array with padded edges.
    """
    top_pad, bottom_pad, left_pad, right_pad = find_edge_padding(im_swir)
    # Apply this padding to your masks or other arrays as needed

    # if bottom pad is 0 the entire image gets set to True
    if bottom_pad > 0:
        im_nodata[-bottom_pad:, :] = True
    # if right pad is 0 the entire image gets set to True
    if right_pad > 0:
        im_nodata[:, -right_pad:] = True

    im_nodata[:, :left_pad] = True
    im_nodata[:top_pad, :] = True
    return im_nodata

def create_jpg(im_ms, cloud_mask, date, satname, filepath, use_matplotlib=True):
    """
    Saves a .jpg file with the RGB image as well as the NIR and SWIR1 grayscale images.
    This functions can be modified to obtain different visualisations of the
    multispectral images.

    KV WRL 2018

    Arguments:
    -----------
    im_ms: np.array
        3D array containing the pansharpened/down-sampled bands (B,G,R,NIR,SWIR1)
    cloud_mask: np.array
        2D cloud mask with True where cloud pixels are
    date: str
        string containing the date at which the image was acquired
    satname: str
        name of the satellite mission (e.g., 'L5')
    filepath: str
        directory in which to save the images
    use_matplotlib: boolean
        False to save a .jpg and True to save as matplotlib plots

    Returns:
    -----------
        Saves a .jpg image corresponding to the preprocessed satellite image

    """
    # rescale image intensity for display purposes
    im_RGB = rescale_image_intensity(im_ms[:,:,[2,1,0]], cloud_mask, 99.9)
    im_NIR = rescale_image_intensity(im_ms[:,:,3], cloud_mask, 99.9)
    im_SWIR = rescale_image_intensity(im_ms[:,:,4], cloud_mask, 99.9)
    
    # creates raw jpg files that can be used for ML applications
    if not use_matplotlib:
        # convert images to bytes so they can be saved
        im_RGB = img_as_ubyte(im_RGB)
        im_NIR = img_as_ubyte(im_NIR)
        im_SWIR = img_as_ubyte(im_SWIR)
        # Save each kind of image with skimage.io
        file_types = ["RGB","SWIR","NIR"]
        # create folders RGB, SWIR, and NIR to hold each type of image
        for ext in file_types:
            ext_filepath = filepath + os.sep + ext
            if not os.path.exists(ext_filepath):
                os.mkdir(ext_filepath)
            # location to save image rgb image would be in sitename/RGB/sitename.jpg
            fname=os.path.join(ext_filepath, date + '_'+ ext +'_' + satname + '.jpg')
            if ext == "RGB":
                imsave(fname, im_RGB, quality=100)
            if ext == "SWIR":
                imsave(fname, im_SWIR, quality=100)
            if ext == "NIR":
                imsave(fname, im_NIR, quality=100)
                
    # if use_matplotlib=True, creates a nicer plot
    else:
        fig = plt.figure()
        fig.set_size_inches([18,9])
        fig.set_tight_layout(True)
        ax1 = fig.add_subplot(111)
        ax1.axis('off')
        ax1.imshow(im_RGB)
        ax1.set_title(date + '   ' + satname, fontsize=16)
        
        # choose vertical or horizontal based on image size
        # if im_RGB.shape[1] > 2.5*im_RGB.shape[0]:
        #     ax1 = fig.add_subplot(311)
        #     ax2 = fig.add_subplot(312)
        #     ax3 = fig.add_subplot(313)
        # else:
        #     ax1 = fig.add_subplot(131)
        #     ax2 = fig.add_subplot(132)
        #     ax3 = fig.add_subplot(133)
        # # RGB
        # ax1.axis('off')
        # ax1.imshow(im_RGB)
        # ax1.set_title(date + '   ' + satname, fontsize=16)
        # # NIR
        # ax2.axis('off')
        # ax2.imshow(im_NIR, cmap='seismic')
        # ax2.set_title('Near Infrared', fontsize=16)
        # # SWIR
        # ax3.axis('off')
        # ax3.imshow(im_SWIR, cmap='seismic')
        # ax3.set_title('Short-wave Infrared', fontsize=16)
    
        # save figure
        fig.savefig(os.path.join(filepath, date + '_' + satname + '.jpg'), dpi=150)

### OTHERS
def find_edge_padding(im_band: np.ndarray) -> np.ndarray:
    """
    Finds the padding required for each edge of an image band based on the presence of data.

    Parameters:
    im_band (numpy.ndarray): The image band.

    Returns:
    tuple: A tuple containing the top, bottom, left, and right padding values.
    """
    # Assuming non-data values are zeros. Adjust the condition if needed.
    is_data = im_band != 0

    # Function to find padding for one edge
    def find_edge_data(is_data_along_edge):
        for idx, has_data in enumerate(is_data_along_edge):
            if has_data:
                return idx
        return len(is_data_along_edge)  # Return full length if no data found

    # Calculate padding for each side
    top_padding = find_edge_data(np.any(is_data, axis=1))
    bottom_padding = find_edge_data(np.any(is_data, axis=1)[::-1])
    left_padding = find_edge_data(np.any(is_data, axis=0))
    right_padding = find_edge_data(np.any(is_data, axis=0)[::-1])

    return top_padding, bottom_padding, left_padding, right_padding

def nd_index(im1, im2, cloud_mask):
    """
    Computes normalised difference index on 2 images (2D), given a cloud mask (2D).

    KV WRL 2018

    Arguments:
    -----------
    im1: np.array
        first image (2D) with which to calculate the ND index
    im2: np.array
        second image (2D) with which to calculate the ND index
    cloud_mask: np.array
        2D cloud mask with True where cloud pixels are

    Returns:    
    -----------
    im_nd: np.array
        Image (2D) containing the ND index
        
    """

    # reshape the cloud mask
    vec_mask = cloud_mask.reshape(im1.shape[0] * im1.shape[1])
    # initialise with NaNs
    vec_nd = np.ones(len(vec_mask)) * np.nan
    # reshape the two images
    vec1 = im1.reshape(im1.shape[0] * im1.shape[1])
    vec2 = im2.reshape(im2.shape[0] * im2.shape[1])
    # compute the normalised difference index
    temp = np.divide(vec1[~vec_mask] - vec2[~vec_mask],
                     vec1[~vec_mask] + vec2[~vec_mask])
    vec_nd[~vec_mask] = temp
    # reshape into image
    im_nd = vec_nd.reshape(im1.shape[0], im1.shape[1])

    return im_nd

def image_std(image, radius):
    """
    Calculates the standard deviation of an image, using a moving window of 
    specified radius. Uses astropy's convolution library'
    
    Arguments:
    -----------
    image: np.array
        2D array containing the pixel intensities of a single-band image
    radius: int
        radius defining the moving window used to calculate the standard deviation. 
        For example, radius = 1 will produce a 3x3 moving window.
        
    Returns:    
    -----------
    win_std: np.array
        2D array containing the standard deviation of the image
        
    """  
    
    # convert to float
    image = image.astype(float)
    # first pad the image
    image_padded = np.pad(image, radius, 'reflect')
    # window size
    win_rows, win_cols = radius*2 + 1, radius*2 + 1
    # calculate std with uniform filters
    win_mean = convolve(image_padded, np.ones((win_rows,win_cols)), boundary='extend',
                        normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)
    win_sqr_mean = convolve(image_padded**2, np.ones((win_rows,win_cols)), boundary='extend',
                        normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)
    win_var = win_sqr_mean - win_mean**2
    win_std = np.sqrt(win_var)
    # remove padding
    win_std = win_std[radius:-radius, radius:-radius]

    return win_std

def convert_epsg(points, epsg_in, epsg_out):
    """
    Converts from one spatial reference to another using the epsg codes
    
    KV WRL 2018

    Arguments:
    -----------
    points: np.array or list of np.ndarray
        array with 2 columns (rows first and columns second)
    epsg_in: int
        epsg code of the spatial reference in which the input is
    epsg_out: int
        epsg code of the spatial reference in which the output will be            
                
    Returns:    
    -----------
    points_converted: np.array or list of np.array 
        converted coordinates from epsg_in to epsg_out
        
    """
    
    # define transformer
    proj = pyproj.Transformer.from_crs(epsg_in, epsg_out, always_xy=True)
    
    # transform points
    if type(points) is list:
        points_converted = []
        # iterate over the list
        for i, arr in enumerate(points): 
            x,y = proj.transform(arr[:,0], arr[:,1])
            arr_converted = np.transpose(np.array([x,y]))
            points_converted.append(arr_converted)
    elif type(points) is np.ndarray:
        x,y = proj.transform(points[:,0], points[:,1])
        points_converted = np.transpose(np.array([x,y]))
    else:
        raise Exception('invalid input type')

    return points_converted

def convert_world2pix(points, georef):
    """
    Converts world projected coordinates (X,Y) to image coordinates 
    (pixel row and column) performing an affine transformation.
    
    KV WRL 2018

    Arguments:
    -----------
    points: np.array or list of np.array
        array with 2 columns (X,Y)
    georef: np.array
        vector of 6 elements [Xtr, Xscale, Xshear, Ytr, Yshear, Yscale]
                
    Returns:    
    -----------
    points_converted: np.array or list of np.array 
        converted coordinates (pixel row and column)
    
    """
    
    # make affine transformation matrix
    aff_mat = np.array([[georef[1], georef[2], georef[0]],
                       [georef[4], georef[5], georef[3]],
                       [0, 0, 1]])
    # create affine transformation
    tform = transform.AffineTransform(aff_mat)
    
    # if list of arrays
    if type(points) is list:
        points_converted = []
        # iterate over the list
        for i, arr in enumerate(points): 
            points_converted.append(tform.inverse(points))
            
    # if single array    
    elif type(points) is np.ndarray:
        points_converted = tform.inverse(points)
        
    else:
        print('invalid input type')
        raise
        
    return points_converted

def convert_pix2world(points, georef):
    """
    Converts pixel coordinates (pixel row and column) to world projected 
    coordinates performing an affine transformation.
    
    KV WRL 2018

    Arguments:
    -----------
    points: np.array or list of np.array
        array with 2 columns (row first and column second)
    georef: np.array
        vector of 6 elements [Xtr, Xscale, Xshear, Ytr, Yshear, Yscale]
                
    Returns:    
    -----------
    points_converted: np.array or list of np.array 
        converted coordinates, first columns with X and second column with Y
        
    """
    
    # make affine transformation matrix
    aff_mat = np.array([[georef[1], georef[2], georef[0]],
                       [georef[4], georef[5], georef[3]],
                       [0, 0, 1]])
    # create affine transformation
    tform = transform.AffineTransform(aff_mat)

    # if list of arrays
    if type(points) is list:
        points_converted = []
        # iterate over the list
        for i, arr in enumerate(points): 
            tmp = arr[:,[1,0]]
            points_converted.append(tform(tmp))
          
    # if single array
    elif type(points) is np.ndarray:
        tmp = points[:,[1,0]]
        points_converted = tform(tmp)
        
    else:
        raise Exception('invalid input type')
        
    return points_converted

def merge_output(output):
    """
    Function to merge the output dictionnary, which has one key per satellite mission
    into a dictionnary containing all the shorelines and dates ordered chronologically.
    
    Arguments:
    -----------
    output: dict
        contains the extracted shorelines and corresponding dates, organised by 
        satellite mission
    
    Returns:    
    -----------
    output_all: dict
        contains the extracted shorelines in a single list sorted by date
    
    """     
    
    # initialize output dict
    output_all = dict([])
    satnames = list(output.keys())
    for key in output[satnames[0]].keys():
        output_all[key] = []
    # create extra key for the satellite name
    output_all['satname'] = []
    # fill the output dict
    for satname in list(output.keys()):
        for key in output[satnames[0]].keys():
            output_all[key] = output_all[key] + output[satname][key]
        output_all['satname'] = output_all['satname'] + [_ for _ in np.tile(satname,
                  len(output[satname]['dates']))]
    # sort chronologically
    idx_sorted = sorted(range(len(output_all['dates'])), key=output_all['dates'].__getitem__)
    for key in output_all.keys():
        output_all[key] = [output_all[key][i] for i in idx_sorted]

    return output_all

def remove_duplicates(output):
    """
    Function to remove from the output dictionnary entries containing shorelines for 
    the same date and satellite mission. This happens when there is an overlap 
    between adjacent satellite images.
    
    KV WRL 2020
    
    Arguments:
    -----------
        output: dict
            contains output dict with shoreline and metadata
        
    Returns:    
    -----------
        output_no_duplicates: dict
            contains the updated dict where duplicates have been removed
        
    """
    # remove duplicates
    dates = output['dates'].copy()
    # find the pairs of images that are within 5 minutes of each other
    time_delta = 5*60 # 5 minutes in seconds
    pairs = []
    for i,date in enumerate(dates):
        # dummy value so it does not match it again
        dates[i] = pytz.utc.localize(datetime(1,1,1) + timedelta(days=i+1))
        # calculate time difference
        time_diff = np.array([np.abs((date - _).total_seconds()) for _ in dates])
        # find the matching times and add to pairs list
        boolvec = time_diff <= time_delta
        if np.sum(boolvec) == 0:
            continue
        else:
            idx_dup = np.where(boolvec)[0][0]
            pairs.append([i,idx_dup])
            
    # if there are duplicates, only keep the longest shoreline
    if len(pairs) > 0:
        # initialise variables
        output_no_duplicates = dict([])
        idx_remove = []
        # for each pair
        for pair in pairs:
            # check if any of the shorelines are empty
            empty_bool = [(len(output['shorelines'][_]) < 2) for _ in pair]
            if np.all(empty_bool): # if both empty remove both
                idx_remove.append(pair[0])
                idx_remove.append(pair[1])
            elif np.any(empty_bool): # if one empty remove that one
                idx_remove.append(pair[np.where(empty_bool)[0][0]])
            else: # remove the shorter shoreline and keep the longer one
                satnames = [output['satname'][_] for _ in pair]
                # keep Landsat 9 if it duplicates Landsat 7
                if 'L9' in satnames and 'L7' in satnames: 
                    idx_remove.append(pair[np.where([_ == 'L7' for _ in satnames])[0][0]])
                else: # keep the longest shorelines
                    sl0 = geometry.LineString(output['shorelines'][pair[0]]) 
                    sl1 = geometry.LineString(output['shorelines'][pair[1]])
                    if sl0.length >= sl1.length: idx_remove.append(pair[1])
                    else: idx_remove.append(pair[0])
        # create a new output structure with all the duplicates removed
        idx_remove = sorted(idx_remove)
        idx_all = np.linspace(0, len(dates)-1, len(dates)).astype(int)
        idx_keep = list(np.where(~np.isin(idx_all,idx_remove))[0])        
        for key in output.keys():
            output_no_duplicates[key] = [output[key][i] for i in idx_keep]
        print('%d duplicates' % len(idx_remove))
        return output_no_duplicates 
    else: 
        print('0 duplicates')
        return output

def remove_inaccurate_georef(output, accuracy):
    """
    Function to remove from the output dictionnary entries containing shorelines 
    that were mapped on images with inaccurate georeferencing:
        - RMSE > accuracy for Landsat images
        - failed geometric test for Sentinel images (flagged with -1)

    Arguments:
    -----------
        output: dict
            contains the extracted shorelines and corresponding metadata
        accuracy: int
            minimum horizontal georeferencing accuracy (metres) for a shoreline to be accepted

    Returns:
    -----------
        output_filtered: dict
            contains the updated dictionnary

    """

    # find indices of shorelines to be removed
    idx = []
    for i in range(len(output['geoaccuracy'])):
        geoacc = output['geoaccuracy'][i]
        if geoacc in ['PASSED','FAILED']:
            if geoacc == 'PASSED':
                idx.append(i)
        else:
            if geoacc <= accuracy:
                idx.append(i)
    # idx = np.where(~(np.array(output['geoaccuracy']) >= accuracy))[0]
    output_filtered = dict([])
    for key in output.keys():
        output_filtered[key] = [output[key][i] for i in idx]
    print('%d bad georef' % (len(output['geoaccuracy']) - len(idx)))
    return output_filtered

def output_to_gdf(output, geomtype):
    """
    Saves the mapped shorelines as a gpd.GeoDataFrame    
    
    KV WRL 2018

    Arguments:
    -----------
    output: dict
        contains the coordinates of the mapped shorelines + attributes
    geomtype: str
        'lines' for LineString and 'points' for Multipoint geometry      
                
    Returns:    
    -----------
    gdf_all: gpd.GeoDataFrame
        contains the shorelines + attirbutes
  
    """    
     
    # loop through the mapped shorelines
    counter = 0
    gdf_all = None
    for i in range(len(output['shorelines'])):
        # skip if there shoreline is empty 
        if len(output['shorelines'][i]) < 2:
            continue
        else:
            # save the geometry depending on the linestyle [removed for points only needed lines]
            if geomtype == 'lines':
                geom = geometry.LineString(output['shorelines'][i]) 
            else:
                raise Exception('geomtype %s is not an option, choose between lines or points'%geomtype)
            # save into geodataframe with attributes
            gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(geom))
            gdf.index = [i]
            gdf.loc[i,'date'] = output['dates'][i].strftime('%Y-%m-%d %H:%M:%S')
            gdf.loc[i,'satname'] = output['satname'][i]
            gdf.loc[i,'geoaccuracy'] = output['geoaccuracy'][i]
            gdf.loc[i,'cloud_cover'] = output['cloud_cover'][i]
            # store into geodataframe
            if counter == 0:
                gdf_all = gdf
            else:
                gdf_all = pd.concat([gdf_all, gdf])
            counter = counter + 1
            
    return gdf_all


## specialized functions
def is_safe_directory_name(name: str) -> bool:
    return re.fullmatch(r'[A-Za-z0-9_]+', name) is not None

def is_more_than_6_months(start: datetime, end: datetime) -> bool:
    return end >= start + relativedelta(months=+6)