from tracer_tools import tools
from pprint import pprint
from page_design import UI
import os

#
import streamlit as st

#the structure goes like this:
# LIST of IMGS
#   - DICT of IMG Metadata
#       - bands (list of bands)
#       - id (earth engine dataset id)
#       - properties (image params)
#       - type (datatype e.g. Image)
#       - version ID (version of the image)


def main2():
    UI._set_session_states()
    MAP = UI.map_face(13.198839454428844, 122.30218706556167)
    INPUTS = UI.process_request(MAP)

    if st.session_state['session_objects'] is None:
        st.info('Consider drawing a polygon on the map to proceed.')
    else:
        with st.form('create-trace'):
            FILTER_IMG = UI.display_image_states()
            SETTINGS = UI.settings()

            submit = st.form_submit_button('Run')
            if submit:
                IMAGES_DOWNLOAD = UI.retrieve_data(st.session_state.last_inputs,
                                                   FILTER_IMG[0], SETTINGS)

    try:
        st.download_button(
            label="Download Filtered Image Metadata (JSON)",
            data=FILTER_IMG[1],
            file_name="sentinel2_image_metadata.json",
            mime="application/json"
        )
    except UnboundLocalError: # data unavailable
        st.toast('Download JSON Metadata disabled.')
        
if __name__ == "__main__":
    main2()