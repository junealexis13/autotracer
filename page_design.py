import streamlit as st
from streamlit_folium import st_folium
from folium import Map, Marker
from folium.plugins import Draw
from tracer_tools import tools
import shoreline_tracer
import os, time, traceback
import json
from datetime import datetime, timedelta, timezone

class UI:

    @staticmethod
    def _set_session_states():
        if 'session_objects' not in st.session_state:
            st.session_state['session_objects'] = None
        if 'last_inputs' not in st.session_state:
            st.session_state['last_inputs'] = None

    @staticmethod
    def map_face(lat, lon):
        m = Map(location=[lat, lon], zoom_start=5, tiles="esri-worldimagery")
        
        Draw(
            export=True,
            filename='drawn_data.geojson',
            position='topleft',
            draw_options={
                'polyline': True,
                'polygon': True,
                'rectangle': True,
                'circle': True,
                'marker': True,
                'circlemarker': False,
            },
            edit_options={'edit': True, 'remove': True},
        ).add_to(m)

        st_data = st_folium(m, use_container_width=True, height=800)
        return st_data


    @staticmethod
    @st.fragment
    def process_request(folium_data):
        objects_store = {}
        
        with st.container(border=True):
            st.title('Map Objects')
            if folium_data["all_drawings"] is not None:
                for i, objects in enumerate(folium_data["all_drawings"]):
                    with st.container(border=True):
                        st.subheader(f"Object {i+1}")
                        colb1, colb2 = st.columns([1,3])
                        with colb1:
                            st.write("__:red[Object Type:]__")
                            st.write("__:red[Object Vertices Coordinates:]__")
                        with colb2:
                            st.text(objects['geometry']['type'])
                            st.text(objects['geometry']['coordinates'])

                            objects_store[f"Object {i+1}"] = objects
            else:
                st.write("Start drawing some vectors on the map")



        with st.form("my_form"):
            st.subheader("Generate a Coastal Trace")
            st.text('Required theres a valid object on the map')
            object_selection = st.selectbox('Select Objects',options=objects_store.keys(), index=0)


            ### checklist ###
            # a set of Checklist to allow users proceed
            check = []
            #################
            from_, to_ = st.columns(2)
            with from_:
                from_date = st.date_input("Select start date", value=datetime.now(timezone.utc).date() - timedelta(days=30))
            with to_:
                to_date = st.date_input("Select end date", value=datetime.now(timezone.utc).date())

            # VALIDATE DATES
            if from_date > to_date:
                st.error("Start date must be before end date.")
                check.append(False)
            else:
                check.append(True)

            if tools.is_more_than_6_months(from_date, to_date):
                st.error('This program was only limited to 6 months max. Please change the date accordingly. If you need more than 6 months trace, start querying again separately for that date range.')
                check.append(False)
            else:
                check.append(True)


            ### parameters
            satellite_missions = ['S2']

            sitename = st.text_input("Site Name (cannot be empty). Take note: Space and '\' or '/' is not allowed")
            
            #check if sitename is valid
            if all([not tools.is_safe_directory_name(sitename), sitename is not None, sitename != ""]):
                st.error('Sitename Invalid! Input characters are not permitted.')
                check.append(False)
            else:
                check.append(True)

            with st.container(border=True):
                st.subheader("Configure Earth engine variables")
                st.caption("If you don't have a Google Earth Engine account, login with your Google Account to console.google.com and create a new project. Or login directly to earthengine.google.com to establish your GEE account and project id.")
                project_id = st.text_input("Earth Engine Project ID", placeholder="your-project-id-in-console")

            ####
            try:
                polygon_shape = objects_store[object_selection]['geometry']['coordinates']
                check.append(True)
            except KeyError:
                st.error('Cannot proceed! You should specify an object first by drawing it.')
                check.append(False)

            submit_button = st.form_submit_button(label="Submit")
            if submit_button:
                if all(check):

                    
                    date_transformed = [from_date.strftime('%Y-%m-%d'), to_date.strftime('%Y-%m-%d')]
        

                    inputs = {
                        'polygon': polygon_shape,
                        'dates': date_transformed,
                        'sat_list': satellite_missions,
                        'sitename': sitename,
                        'ee_project_id': project_id,
                        'filepath': os.getcwd(),
                            }
                    

                    with st.spinner("Checking for available images..."):
                        try:
                            img = tools.check_images_available(inputs) #returns t1 and t2 images in a tuple
                        except Exception as e:
                            st.error(f"Error checking images: {e}")
                            return None # cut here
                        
                    s2Img = img[0]['S2']

                    st.toast("Images fetched successfully!")
                    st.info(f"There are {len(img[0]['S2'])} images available for the selected date range and polygon.")
                    
                    #####################################
                    #view list
                    st.subheader('Available Images:')
                    with st.expander("View Image Dates", expanded=False):
                        for img in s2Img:
                            st.write(tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S'))
                    #####################################

                    # Filter the images based on the selected dates
                    #create a filtered image list 

                    # filtered_images = st.multiselect(
                    #     'Select Images to Process',
                    #     options=[tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') for img in s2Img],
                    #     default=[tools.convert_unix_to_dt(s2Img[0]['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S')]
                    # )

                    # # Filter the images based on the selected dates
                    # filtered_s2Img = [img for img in s2Img if app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') in filtered_images]

                    st.session_state['session_objects'] = s2Img
                    st.session_state['last_inputs'] = inputs
                    st.rerun()

                else:
                    st.error('You have issues with your input. Kindly address them first.')
    @staticmethod
    def display_image_states():
        with st.container(border=True):
            st.subheader("Sentinel-2 Images")
            st.write(f'There are {len(st.session_state["session_objects"])} images available for processing.')


            if st.session_state['session_objects'] is not None:
                metadata_json = json.dumps(st.session_state['session_objects'], indent=1, default=str)
                with st.expander("View Image Dates", expanded=False):
                    for img in st.session_state['session_objects']:
                        st.write(tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S'))


            #create a multiselect to select images to process
            filtered_images = st.multiselect(
                'Select Images to Process',
                options=[tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') for img in st.session_state['session_objects']],
                default=[tools.convert_unix_to_dt(st.session_state['session_objects'][0]['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S')]
            )

            # Filter the images based on the selected dates
            filtered_s2Img = [img for img in st.session_state['session_objects'] if tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') in filtered_images]

            return filtered_s2Img, metadata_json
        

    @staticmethod
    @st.fragment
    def settings():
        with st.container(border=True):
            st.title('Settings')
            st.subheader('General Settings')
            CLOUD_TRESHOLD = st.slider(
                        'Cloud Threshold (%)',
                        min_value=0.0,
                        max_value=1.0,
                        value=0.3,
                        step=0.01,
                        help='threshold on maximum cloud cover'
                    )
            MAX_CLOUD_DIST = st.number_input(
                        'Maximum Cloud Distance (m)',
                        min_value=0,
                        max_value=1000,
                        value=300,
                        help='Maximum distance around clouds where shoreline cant be mapped'
                    )
            OUTPUT_EPSG = st.number_input(
                        'Output EPSG Code',
                        min_value=0,
                        max_value=999999,
                        value=4326,
                        help='EPSG code for the output coordinate system (default is Mercator WGS84 - 4326)'
                )
            
            st.subheader('QC Controls')
            CHECK_DETECTION = st.checkbox(
                        'Check Detection',
                        value=False,
                        help='If checked, shows each shoreline detection to the user for validation. [Disabled]',
                        disabled=True
                    )
            ADJUST_DETECTION = st.checkbox(
                        'Adjust Detection',
                        value=False,
                        help='If checked, allows user to adjust the postion of each shoreline by changing the threhold'
                    )
            SAVE_FIGURE = st.checkbox(
                        'Save Figure',
                        value=True,
                        help='If checked, saves the figure showing the mapped shoreline for each image'
                    )
            
            st.subheader('Advanced Settings - Optional')
            MIN_BEACH_AREA = st.number_input(
                        'Minimum Beach Area (m²)',
                        min_value=0,
                        max_value=100000000,
                        value=1000,
                        help='Minimum area of the beach to be considered for shoreline mapping'
                    )
            MIN_SHORELINE_LENGTH = st.number_input(
                        'Minimum Shoreline Length (m)',
                        min_value=0,
                        max_value=100000,
                        value=500,
                        help='Minimum length of the shoreline to be considered for mapping'
                    )
            S2CLOUD_PROBABILITY = st.slider(
                        'S2 Cloud Probability Threshold',
                        min_value=0,
                        max_value=100,
                        value=40,
                        step=1,
                        help='Threshold for Sentinel-2 cloud probability'
                    )
            SAND_COLOR = st.selectbox(
                        'Sand Color',
                        options=['default', 'latest', 'dark', 'bright'],
                        index=0,
                        help='Color of the sand to be used for shoreline mapping. Pertains to the color of the sand in the images.'
                    )
            CLOUD_MASKING = st.checkbox(
                        'Cloud Mask Issue',
                        value=False,
                        help='if Checked, set this if sand pixels are masked (in black) on many images'
                    )
            LANDSAT_PANSHARPENING = st.checkbox(
                        'Landsat Pansharpening',
                        value=False,
                        disabled=True,  # Disabled until Landsat support is implemented
                        help='if Checked, uses pansharpened Landsat images for shoreline mapping'
                    )
            
            return {
                # general parameters:
                'cloud_thresh': CLOUD_TRESHOLD,       
                'dist_clouds': MAX_CLOUD_DIST,        
                'output_epsg': OUTPUT_EPSG,      
                # quality control:
                'check_detection': CHECK_DETECTION,   
                'adjust_detection': ADJUST_DETECTION,  
                'save_figure': SAVE_FIGURE,        
                # [ONLY FOR ADVANCED USERS] 
                'min_beach_area': MIN_BEACH_AREA,     
                'min_length_sl': MIN_SHORELINE_LENGTH,      
                'cloud_mask_issue': CLOUD_MASKING,  
                'sand_color': SAND_COLOR,   
                'pan_off': LANDSAT_PANSHARPENING,      #disabled for this setup     
                's2cloudless_prob': S2CLOUD_PROBABILITY,     
                # add the inputs defined previously
                'inputs': st.session_state['last_inputs'],
                    }
        



    @staticmethod
    def retrieve_data(inputs, imlist, settings: dict,):
        #download the data
        with st.spinner('Fetching Data in progress', show_time=True):
            data = tools.retrieve_images(inputs,imlist)
            data = tools.get_metadata(inputs=inputs)

        st.toast(f'[PROCESS] Data collected: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

        tools.save_jpg(data,settings,use_matplotlib=True)
        st.toast("[PROCESS] Images saved")
        time.sleep(5)
        st.toast('[PROCESS] Attempting to do Shoreline Trace')

        try:
            OUTPUT = shoreline_tracer.extract_shorelines(data, settings)
            OUTPUTS = tools.remove_duplicates(OUTPUT)
            OUTPUTS = tools.remove_inaccurate_georef(OUTPUTS, 10)
            shoreline_tracer.generate_GIS_output()
        except Exception as e:
            st.error(e)
            st.write(traceback.print_exc())

        st.toast('[UDPATE] Trace Complete!')
