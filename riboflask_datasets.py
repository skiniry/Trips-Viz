import sys
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
import os
import re
import subprocess
import shelve
import mpld4
from mpld4 import plugins,utils
import collections
from mpld4.utils import get_id
import pandas as pd
from fetch_shelve_reads2 import get_reads
import numpy as np

from bokeh.plotting import figure, show, output_file, reset_output
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette
from bokeh.models import CustomJS
from bokeh.layouts import column


from bokeh.layouts import widgetbox
from bokeh.models.widgets import Select

from bokeh.io import show
from bokeh.models import (
	TapTool,
	OpenURL,
	LogTicker,
	ColumnDataSource,
	HoverTool,
	LinearColorMapper,
	LogColorMapper,
	BasicTicker,
	PrintfTickFormatter,
	ColorBar
)




# Define some CSS to control our custom labels
point_tooltip_css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #000000;
  background-color: #d2d4d8;
}
td
{
  background-color: #ffffff;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 0px solid black;
  text-align: left;
}
"""

def merge_dict(dict1,dict2):
	master_dict = dict1
	for key in dict2:
		if key in master_dict:
			master_dict[key] += dict2[key]
		else:
			master_dict[key] = dict2[key]
	return master_dict


def set_axis_color(axis, color, alpha=None):
	"""Sets the spine color of all sides of an axis (top, right, bottom, left)."""
	for side in ('top', 'right',  'bottom', 'left'):
		spine = axis.spines[side]
		spine.set_color(color)
		if alpha is not None:
			spine.set_alpha(alpha)


def get_color_palette(scheme):
	"""Return colors for a given scheme. Default colors are returned for an item
	if undefined in scheme.

	"""
	color_schemes = {
		'default': {
			'frames': ['#FF4A45', '#64FC44', '#5687F9'], 'background': '#ffffff',
			'color': '#616161', 'ticks': '#757575', 'start': '#ffffff', 'stop': '#909090',
			'rna': '#BFBFBF', 'axis': '#e0e0e0', 'grey': '#bdbdbd'
		},
		'colorbrewer': {
			'frames': ['#fc8d62', '#66c2a5', '#8da0cb']
		},
		'rgb': {
			'frames': ['red', 'green', 'blue']
		},
		'greyorfs': {}
	}

	colors = {}
	for k, v in color_schemes['default'].items():
		try:
			vals = color_schemes[scheme][k]
		except KeyError:
			vals = v
		colors[k] = vals
	return colors


def get_near_cog_starts(seq):
	near_cog_starts = {0:[],1:[],2:[]}
	for i in range(0,len(seq)):
		codon = seq[i:i+3]
		if codon == "CTG":
			near_cog_starts[(i+1)%3].append(i+1)
	return near_cog_starts

def plot_profile(xlist,ylist,filenames,file_descs,studies,raw_reads,controls, cell_lines,control_colors, study_colors, cell_line_colors, transcript, start, stop):
	reset_output()
	TOOLS="hover,pan,wheel_zoom,reset"
	full_title = "Dataset breakdown for {} from {} to {}".format(transcript, start, stop)
	
	p = figure(plot_width=1200,plot_height=550,tools=TOOLS, y_axis_label = "Normalised read count", title=full_title)
	p.title.align="center"
	p.title.text_font_size = "18pt"
	p.xaxis.axis_label_text_font_size = "26pt"
	p.xaxis.major_label_text_font_size = "26pt"
	p.yaxis.axis_label_text_font_size = "26pt"
	p.yaxis.major_label_text_font_size = "26pt"
	p.background_fill_color = "#e0e1e2"
	p.xaxis.major_tick_line_color = None
	p.xaxis.minor_tick_line_color = None
	p.xaxis.major_label_text_font_size = '0pt' 
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"

	csource = ColumnDataSource({'x': xlist,'y':ylist,'filenames':filenames, 'file_descs':file_descs, "studies":studies, 
							   "raw_reads":raw_reads,"controls":controls, "cell_lines":cell_lines,"control_colors":control_colors,
							   "study_colors":study_colors, "cell_line_colors":cell_line_colors})

	sct = p.scatter('x', 'y',source=csource,size=12,fill_color="study_colors",line_color="study_colors")
	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Filename", "@filenames"),("Desc","@file_descs"),("Study","@studies"),("Normalised read count","@y"), ("Raw read count","@raw_reads"),("Control","@controls"),("Cell line","@cell_lines")]
	hover.mode = 'mouse'

	cb_cselect = CustomJS(args=dict(sct=sct,csource=csource), code ="""
		var selected_color = cb_obj.value;
		sct.glyph.line_color.field = selected_color;
		sct.glyph.fill_color.field = selected_color;
		sct.glyph.change.emit();
	""")
	color_select = Select(title="Choose color pallette:", value="study_colors", 
						options = ["study_colors","control_colors","cell_line_colors"], callback = cb_cselect)

	

	output_file("color_scatter.html", title="color_scatter.py example")

	#show(p)  # open a browser
	layout = column(color_select, p)
	#graph = "<link href='https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.12.14.min.css' rel='stylesheet' type='text/css'><script src='https://cdn.pydata.org/bokeh/release/bokeh-0.12.14.min.js'></script><script src='https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.12.14.min.js'></script>"
	graph = file_html(layout,CDN)
	reset_output()
	return graph
