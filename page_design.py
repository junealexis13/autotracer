import streamlit as st
from streamlit_folium import st_folium
from folium import Map, Marker
from folium.plugins import Draw
from tracer_tools import tools
from tracer_tools import app_tools
import os

from datetime import datetime, timedelta, timezone

class UI:

    @staticmethod
    def _set_session_states():
        if 'session_objects' not in st.session_state:
            st.session_state['session_objects'] = None

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

            from_, to_ = st.columns(2)
            with from_:
                from_date = st.date_input("Select start date", value=datetime.now(timezone.utc).date() - timedelta(days=30))
            with to_:
                to_date = st.date_input("Select end date", value=datetime.now(timezone.utc).date())

            if from_date > to_date:
                st.error("Start date must be before end date.")

            ### parameters
            satellite_missions = ['S2']
            sitename = st.text_input("Site Name (cannot be empty)")
            
            with st.container(border=True):
                st.subheader("Configure Earth engine variables")
                project_id = st.text_input("Earth Engine Project ID", placeholder="your-project-id-in-console")

            ####


            submit_button = st.form_submit_button(label="Submit")
            if submit_button:

                polygon_shape = objects_store[object_selection]['geometry']['coordinates']
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
                        st.write(app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S'))
                #####################################

                # Filter the images based on the selected dates
                #create a filtered image list 

                # filtered_images = st.multiselect(
                #     'Select Images to Process',
                #     options=[app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') for img in s2Img],
                #     default=[app_tools.convert_unix_to_dt(s2Img[0]['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S')]
                # )

                # # Filter the images based on the selected dates
                # filtered_s2Img = [img for img in s2Img if app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') in filtered_images]

                st.session_state['session_objects'] = s2Img
                st.rerun()


    @staticmethod
    def display_image_states():
        with st.container(border=True):
            st.subheader("Sentinel-2 Images")
            st.write(f'There are {len(st.session_state["session_objects"])} images available for processing.')


            if st.session_state['session_objects'] is not None:
                with st.expander("View Image Metadata", expanded=False):
                    for i, img in enumerate(st.session_state['session_objects']):
                        st.write(f"Image {i+1}:")
                        st.json(img)
                with st.expander("View Image Dates", expanded=False):
                    for img in st.session_state['session_objects']:
                        st.write(app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S'))


            #create a multiselect to select images to process
            filtered_images = st.multiselect(
                'Select Images to Process',
                options=[app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') for img in st.session_state['session_objects']],
                default=[app_tools.convert_unix_to_dt(st.session_state['session_objects'][0]['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S')]
            )

            # Filter the images based on the selected dates
            filtered_s2Img = [img for img in st.session_state['session_objects'] if app_tools.convert_unix_to_dt(img['properties']['system:time_start']).strftime('%Y-%m-%d %H:%M:%S') in filtered_images]

            return filtered_s2Img