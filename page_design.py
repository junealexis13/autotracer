import streamlit as st
from streamlit_folium import st_folium
from folium import Map, Marker
from folium.plugins import Draw

def map_face(lat, lon):
    m = Map(location=[lat, lon], zoom_start=12)
    st_data = st_folium(m, width=700, height=500)
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
    return st_data