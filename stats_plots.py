import sys
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import StringIO
import os 
import re 
import subprocess 
import shelve 
import mpld3 
from mpld3 import plugins,utils 
import collections 
from mpld3.utils import get_id
import pandas as pd
import numpy as np
from flask import make_response



from bokeh.plotting import figure, show, output_file
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette
from bokeh.palettes import all_palettes
from bokeh.models.glyphs import Text

from bokeh.io import show
from bokeh.models import (
    Range1d,
    Label,
    FuncTickFormatter,
    LogTicker, 
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    LogColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar
)




redhex="#FF5F5B"
greenhex="#90E090"
bluehex="#9ACAFF"
yellowhex="#FFFF91"
background_col = "#EEEEF4"


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
  border: 2px solid black;
  text-align: right;
}
"""








def year_dist(final_year_dict): 
    df = pd.DataFrame(final_year_dict, columns = ['no_studies', 'year'])
    pos = list(range(len(df['no_studies'])))
    width = 0.25
    fig, ax = plt.subplots( figsize=(23,12))

    
    # Create a bar with riboseq files data in position pos,
    plt.bar(pos,df['no_studies'], width, alpha=1,color=redhex,linewidth=0,label="Number of studies")#df['readlengths'][0])

    # Set the y axis label
    ax.set_ylabel('Number of studies', labelpad=100,fontsize="21")
    ax.set_xlabel('Publication year',fontsize="21")
    # Set the chart's title
    ax.set_title("No. of studies per year",y=1.05,fontsize="25")

    # Set the position of the x ticks
    ax.set_xticks([p + 1 * width for p in pos])
    
    # Set the labels for the x ticks
    ax.set_xticklabels(df['year'])

    # Setting the x-axis and y-axis limits
    plt.xlim(min(pos)-width, max(pos)+width*4)
    plt.ylim([0, max(df['no_studies'])*1.1])


    ax.set_axis_bgcolor(background_col)
    ax.tick_params('both', labelsize=16)
    plt.grid(color="white", linewidth=2,linestyle="solid")
    plugins.connect(fig, plugins.TopToolbar(xoffset=-13, yoffset=115),plugins.DownloadPNG(returnstr="Study distribution by year"))
    return mpld3.fig_to_html(fig)

def org_breakdown_plot(read_dict): 
    df = pd.DataFrame(read_dict, columns = ['riboseq_files', 'rnaseq_files', 'organisms'])
    pos = list(range(len(df['riboseq_files'])))
    width = 0.25
    fig, ax = plt.subplots( figsize=(23,12))

    
    # Create a bar with riboseq files data in position pos,
    plt.bar(pos,df['riboseq_files'], width, alpha=1,color=redhex,linewidth=0,label="Riboseq files")#df['readlengths'][0])
    # Create a bar with rnaseq files data in position pos + some width buffer,
    plt.bar([p + width for p in pos],df['rnaseq_files'],width,alpha=1,color=greenhex,linewidth=0,label="Rnaseq files")

    # Set the y axis label
    ax.set_ylabel('Count', labelpad=100,fontsize="21")
    ax.set_xlabel('Organisms',fontsize="21")
    # Set the chart's title
    ax.set_title("No. of files per organism",y=1.05,fontsize="25")

    # Set the position of the x ticks
    ax.set_xticks([p + 1 * width for p in pos])
    
    # Set the labels for the x ticks
    ax.set_xticklabels(df['organisms'])

    # Setting the x-axis and y-axis limits
    plt.xlim(min(pos)-width, max(pos)+width*4)
    plt.ylim([0, max(max(df['riboseq_files']),max(df['rnaseq_files']))*1.1] )

    # Adding the legend and showing the plot
    leg = plt.legend(['Riboseq files', 'Rnaseq files'], loc='upper right')
    leg.get_frame().set_edgecolor('#D2D2EB')
    ax.set_axis_bgcolor(background_col)
    ax.tick_params('both', labelsize=16)
    plt.grid(color="white", linewidth=2,linestyle="solid")
    plugins.connect(fig, plugins.TopToolbar(xoffset=-13, yoffset=115),plugins.DownloadPNG(returnstr="Organism breakdown"))
    return mpld3.fig_to_html(fig)
    
