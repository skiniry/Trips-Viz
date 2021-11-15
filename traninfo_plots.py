import sys
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import os
import re
import subprocess
import shelve
import mpld3
from mpld3 import plugins,utils
from new_plugins import InteractiveLegendPlugin,TopToolbar,DownloadProfile,DownloadPNG
import collections
from mpld3.utils import get_id
import pandas as pd
import numpy as np
from flask import make_response
from scipy.stats.stats import spearmanr,pearsonr
import matplotlib.cm as cm
from bokeh.plotting import figure, show, output_file
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette
from bokeh.palettes import inferno
from bokeh.palettes import all_palettes
from bokeh.models.glyphs import Text
import bokeh.models as bmo
from bokeh.io import show
from bokeh.models import (
	TapTool,
	OpenURL,
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
#ViennaRNA can be installed from here https://github.com/ViennaRNA/ViennaRNA
try:
	import RNA
	vienna_rna = True
except:
	print ("Could not import the RNA module, ViennaRNA needs to be installed (https://github.com/ViennaRNA/ViennaRNA), MFE will not be plotted on traninfo plot")
	vienna_rna = False


redhex="#FF5F5B"
greenhex="#90E090"
bluehex="#9ACAFF"
yellowhex="#FFFF91"



# Define some CSS to control our custom labels
line_tooltip_css = """
.tooltip
{
  color: #000000;
  background-color: #d2d4d8;
  font-family:Arial, Helvetica, sans-serif;
  text-align: left;
}

"""


def nuc_freq_plot(master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size,filename):
	labels = ["A","T","G","C"]
	returnstr = "Position,A,T,G,C\n"
	minpos = min(master_dict.keys())
	maxpos = max(master_dict.keys())
	x_pos = []
	a_counts = []
	t_counts = []
	g_counts = []
	c_counts = []
	for i in range(minpos,maxpos):
		returnstr += "{},{:.2f},{:.2f},{:.2f},{:.2f}\n".format(i,master_dict[i]["A"],master_dict[i]["T"],master_dict[i]["G"],master_dict[i]["C"])
		x_pos.append(i)
		a_counts.append(master_dict[i]["A"])
		t_counts.append(master_dict[i]["T"])
		g_counts.append(master_dict[i]["G"])
		c_counts.append(master_dict[i]["C"])
		
	fig, ax = plt.subplots( figsize=(13,12))
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Position (nucleotides)',fontsize=axis_label_size)

	#if nuc_comp_type == "nuc_comp_per":
	#	ax.set_ylim(0,1)
	#	ax.set_ylabel('Percent',fontsize=axis_label_size,labelpad=50)
	#elif nuc_comp_type == "nuc_comp_count":
	maxheight = max(max(a_counts), max(t_counts),max(g_counts),max(c_counts))
	ax.set_ylim(0,maxheight)
	ax.set_ylabel('Count',fontsize=axis_label_size,labelpad=100)
	#if nuc_comp_direction == "nuc_comp_five":
	#	ax.set_xlim(0,maxreadlen)
	#elif nuc_comp_direction == "nuc_comp_three":
	#	ax.set_xlim((maxreadlen*-1),-1)
	width = 0.95
	#plot it
	ax = plt.subplot(111)
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str, y=1.05,fontsize=title_size)
	a_line = ax.plot(x_pos, a_counts, label=labels, color="blue", linewidth=4)
	t_line = ax.plot(x_pos, t_counts, label=labels, color="red", linewidth=4)
	g_line = ax.plot(x_pos, g_counts, label=labels, color="green", linewidth=4)
	c_line = ax.plot(x_pos, c_counts, label=labels, color="orange", linewidth=4)
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")
	ilp = InteractiveLegendPlugin([a_line, t_line, g_line, c_line], ["A","T","G","C"], alpha_unsel=0,alpha_sel=1,start_visible=True)
	plugins.connect(fig, ilp,TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)	
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as fasta file</b></button></a> </div>".format(short_code,filename)
	graph += mpld3.fig_to_html(fig)
	return graph




def nuc_comp_single(tran, master_dict,title,short_code,background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size,traninfo):
	labels = ["Exon Junctions","CDS markers"]
	start_visible=[True,True]
	stop_codons = ["TAG","TAA","TGA"]
	frame_orfs = {1:[],2:[],3:[]}
	color_dict = {'frames': ['#FF4A45', '#64FC44', '#5687F9']}
	try:
		traninfo["stop_list"] = [int(x) for x in traninfo["stop_list"]]
	except:
		traninfo["stop_list"] = []

	try:
		traninfo["start_list"] = [int(x) for x in traninfo["start_list"]]
	except:
		traninfo["start_list"] = []

	if str(traninfo["exon_junctions"][0]) != "":
		traninfo["exon_junctions"] = [int(x) for x in traninfo["exon_junctions"]]
	else:
		traninfo["exon_junctions"] = []
	gene = traninfo["gene"]
	tranlen = traninfo["length"]
	cds_start = traninfo["cds_start"]
	cds_stop = traninfo["cds_stop"]
	if cds_start == "NULL" or cds_start == None:
		cds_start = 0
	if cds_stop == "NULL" or cds_stop == None:
		cds_stop = 0
	all_starts = traninfo["start_list"]
	all_stops = {"TAG":[],"TAA":[],"TGA":[]}
	exon_junctions = traninfo["exon_junctions"]
	seq = traninfo["seq"].upper()
	for i in range(0,len(seq)):
		if seq[i:i+3] in stop_codons:
			all_stops[seq[i:i+3]].append(i+1)
	# Error occurs if one of the frames is empty for any given start/stop, so we initialise with -5 as this won't be seen by user and will prevent the error
	start_stop_dict = {1:{"starts":[-5], "stops":{"TGA":[-5],"TAG":[-5],"TAA":[-5]}},
				   2:{"starts":[-5], "stops":{"TGA":[-5],"TAG":[-5],"TAA":[-5]}},
				   3:{"starts":[-5], "stops":{"TGA":[-5],"TAG":[-5],"TAA":[-5]}}}
	for start in all_starts:
		rem = ((start-1)%3)+1
		start_stop_dict[rem]["starts"].append(start)
	for stop in all_stops:
		for stop_pos in all_stops[stop]:
			rem = ((stop_pos-1)%3)+1
			start_stop_dict[rem]["stops"][stop].append(stop_pos)
	#find all open reading frames
	for frame in [1,2,3]:
		for start in start_stop_dict[frame]["starts"]:
			best_stop_pos = 10000000
			for stop in start_stop_dict[frame]["stops"]:
				for stop_pos in start_stop_dict[frame]["stops"][stop]:
					if stop_pos > start and stop_pos < best_stop_pos:
						best_stop_pos = stop_pos
			if best_stop_pos != 10000000:
				frame_orfs[frame].append((start, best_stop_pos))
	y_max = 100
	fig = plt.figure(figsize=(13,8))

	ax_main = plt.subplot2grid((30,1), (0,0),rowspan=22)
	ax_main.spines['bottom'].set_visible(False)
	#ax_main.set_ylabel(label,  fontsize=axis_label_size, labelpad=30)
	label = 'Position (nucleotides)'
	ax_main.set_xlabel(label, fontsize=axis_label_size)
	ax_main.set_ylim(0, y_max)
	cds_markers = ax_main.plot((cds_start+1,cds_start+1), (0, y_max), color="black",linestyle = 'solid', linewidth=2)
	cds_markers += ax_main.plot((cds_stop+1,cds_stop+1), (0, y_max), color="black",linestyle = 'solid', linewidth=2)
	ax_f1 = plt.subplot2grid((30,1), (26,0),rowspan=1,sharex=ax_main)
	ax_f1.set_facecolor(color_dict['frames'][0])
	ax_f2 = plt.subplot2grid((30,1), (27,0),rowspan=1,sharex=ax_main)
	ax_f2.set_facecolor(color_dict['frames'][1])
	ax_f3 = plt.subplot2grid((30,1), (28,0),rowspan=1,sharex=ax_main)
	ax_f3.set_facecolor(color_dict['frames'][2])
	ax_nucseq = plt.subplot2grid((30,1), (29,0),rowspan=1,sharex=ax_main)
	ax_nucseq.set_xlabel('Transcript: {} Length: {} nt'.format(tran, tranlen), fontsize=subheading_size, labelpad=10)

	#plot a dummy exon junction at postion -1, needed in cases there are no exon junctions, this wont be seen
	allexons = ax_main.plot((-1,-1), (0, 1), alpha=0.01,color='black',linestyle = '-.', linewidth=2)
	for exon in exon_junctions:
		allexons += ax_main.plot((exon,exon), (0, y_max), alpha=0.95,color='black',linestyle = '-.', linewidth=3)
	xy = 0
	ax_nucseq.set_facecolor(background_col)
	mrnaseq = seq.replace("T","U")
	color_list = ["#FF4A45","#64FC44","#5687F9"]
	char_frame = 0
	for char in mrnaseq:
		ax_nucseq.text((xy+1)-0.1,0.2,mrnaseq[xy],fontsize=20,color=color_list[char_frame%3])
		xy += 1
		char_frame += 1
	for axisname in (ax_f1, ax_f2, ax_f3,ax_nucseq):
		axisname.tick_params(top=False, bottom=False, labelleft=False, labelright=False, labelbottom=False)
	for label in ax_main.xaxis.get_majorticklabels():
		label.set_fontsize(36)
	for axis, frame in ((ax_f1, 1), (ax_f2, 2), (ax_f3, 3)):
		axis.set_xlim(1, tranlen)
		starts = [(item, 1) for item in start_stop_dict[frame]['starts']]
		uag_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TAG']]
		uaa_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TAA']]
		uga_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TGA']]
		axis.broken_barh(starts, (0.30, 1),color="white", zorder=2,linewidth=7)
		axis.broken_barh(uag_stops, (0, 1), color="black", zorder=2, linewidth=4)
		axis.broken_barh(uaa_stops, (0, 1), color="black", zorder=2, linewidth=4)
		axis.broken_barh(uga_stops, (0, 1), color="black", zorder=2, linewidth=4)
		axis.set_ylim(0, 1)
		axis.set_ylabel('{}'.format(frame), labelpad=10, verticalalignment='center',rotation="horizontal",color="black")
	title_str = '{} ({})'.format(gene,short_code)
	plt.title(title_str, fontsize=title_size,y=36)
	line_collections = [allexons,cds_markers]

	plot_gc = True
	plot_mfe = True
	if plot_mfe == True and vienna_rna == True:
		step_size = 2
		window_size = 60
		mfe_dict = collections.OrderedDict()
		for i in range(0,len(seq)-(window_size),step_size):
			seq_window = str(seq[i:i+window_size])
			(ss, mfe) = RNA.fold(seq_window)
			mfe_dict[i+(window_size/2)] = abs(mfe)
	else:
		mfe_dict = {}
	if plot_gc == True:
		step_size = 2
		window_size = 60
		a_dict = collections.OrderedDict()
		t_dict = collections.OrderedDict()
		g_dict = collections.OrderedDict()
		c_dict = collections.OrderedDict()
		gc_dict = collections.OrderedDict()
		for i in range(0,len(seq)-(window_size),step_size):
			a_count = 0.0
			t_count = 0.0
			g_count = 0.0
			c_count = 0.0
			for x in range(i,i+window_size):
				if seq[x] == "A":
					a_count += 1
				elif seq[x] == "T":
					t_count += 1
				elif seq[x] == "G":
					g_count += 1
				elif seq[x] == "C":
					c_count += 1

			gc_count = g_count + c_count
			norm_a = a_count/window_size
			norm_t = t_count/window_size
			norm_g = g_count/window_size
			norm_c = c_count/window_size
			norm_gc = gc_count/window_size

			final_a = norm_a*y_max
			final_t = norm_t*y_max
			final_g = norm_g*y_max
			final_c = norm_c*y_max
			final_gc = norm_gc*y_max
			a_dict[i+(window_size/2)] = final_a
			t_dict[i+(window_size/2)] = final_t
			g_dict[i+(window_size/2)] = final_g
			c_dict[i+(window_size/2)] = final_c
			gc_dict[i+(window_size/2)] = final_gc
		a_plot = ax_main.plot(a_dict.keys(), a_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][0], linewidth=4)
		t_plot = ax_main.plot(t_dict.keys(), t_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][1], linewidth=4)
		g_plot = ax_main.plot(g_dict.keys(), g_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][2], linewidth=4)
		c_plot = ax_main.plot(c_dict.keys(), c_dict.values(), alpha=0.01, label = labels, zorder=1, color='#ffff99', linewidth=4)
		gc_plot = ax_main.plot(gc_dict.keys(), gc_dict.values(), alpha=1, label = labels, zorder=1, color='grey', linewidth=4)
		mfe_plot = ax_main.plot(mfe_dict.keys(), mfe_dict.values(), alpha=0.01, label = labels, zorder=1,color='#df8500', linewidth=4)
		for item,lbl,viz in [(a_plot,"A%",False),(t_plot,"T%",False),(g_plot,"G%",False),(c_plot,"C%",False),(gc_plot,"GC%",True),(mfe_plot,"MFE",False)]:
			line_collections.append(item)
			labels.append(lbl)
			start_visible.append(viz)

	leg_offset = (30-17)*5
	if leg_offset <0:
		leg_offset = 0

	ilp = InteractiveLegendPlugin(line_collections, labels, alpha_unsel=0,alpha_sel=0.85,start_visible=start_visible,fontsize=30,xoffset=leg_offset)
	htmllabels = {1:[],2:[],3:[]}
	all_start_points = {1:[],2:[],3:[]}

	points1 =ax_f1.plot(all_start_points[1], [0.75]*len(all_start_points[1]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)
	points2 =ax_f2.plot(all_start_points[2], [0.75]*len(all_start_points[2]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)
	points3 =ax_f3.plot(all_start_points[3], [0.75]*len(all_start_points[3]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)

	ax_f3.axes.get_yaxis().set_ticks([])
	ax_f2.axes.get_yaxis().set_ticks([])
	ax_f1.axes.get_yaxis().set_ticks([])

	plugins.connect(fig, ilp, TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=""),DownloadPNG(returnstr=tran))

	ax_main.set_facecolor(background_col)
	# This changes the size of the tick markers, works on both firefox and chrome.
	ax_main.tick_params('both', labelsize=marker_size)
	ax_main.xaxis.set_major_locator(plt.MaxNLocator(3))
	ax_main.yaxis.set_major_locator(plt.MaxNLocator(3))
	ax_main.grid(True, color="white", linewidth=30,linestyle="solid")
	#Without this style tag the markers sizes will appear correct on browser but be original size when downloaded via png
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph +=  mpld3.fig_to_html(fig)
	return graph

def gc_metagene(title, short_code, background_col, readlength_col, title_size, axis_label_size, subheading_size, marker_size, traninfo):
	labels = ["CDS markers"]
	start_visible=[True]
	stop_codons = ["TAG","TAA","TGA"]
	frame_orfs = {1:[],2:[],3:[]}
	color_dict = {'frames': ['#FF4A45', '#64FC44', '#5687F9']}
	gene = ""
	y_max = 100
	fig = plt.figure(figsize=(13,8))

	ax_main = plt.subplot2grid((30,1), (0,0),rowspan=22)
	ax_main.spines['bottom'].set_visible(False)
	ax_main.set_ylabel("%",  fontsize=axis_label_size, labelpad=30)
	ax_main.set_ylim(0, y_max)
	ax_main.set_xlim(0, 1500)
	cds_markers = ax_main.plot((500,500), (0, y_max-3), color="black",linestyle = 'solid', linewidth=2)
	cds_markers += ax_main.plot((1000,1000), (0, y_max-3), color="black",linestyle = 'solid', linewidth=2)
	#ax_nucseq.set_xlabel('Transcript: {} Length: {} nt'.format(tran, tranlen), fontsize=subheading_size, labelpad=10)
	xy = 0
	color_list = ["#FF4A45","#64FC44","#5687F9"]
	for label in ax_main.xaxis.get_majorticklabels():
		label.set_fontsize(36)

	title_str = '{} ({})'.format(gene,short_code)
	plt.title(title_str, fontsize=title_size,y=36)
	line_collections = [cds_markers]

	plot_gc = True
	plot_mfe = False
	if plot_mfe == True and vienna_rna == True:
		step_size = 10
		window_size = 60
		mfe_dict = collections.OrderedDict()
		for item in traninfo:
			transcript = item[0]
			cds_start = float(item[1])
			cds_stop = float(item[2])
			seq = item[3]
			seqlen = len(seq)
			for i in range(0,len(seq)-(window_size),step_size):
				seq_window = str(seq[i:i+window_size])
				(ss, mfe) = RNA.fold(seq_window)
				if i < cds_start:
					per = (i+(window_size/2)/cds_start)*5
				if i >= cds_start and i <= cds_stop:
					per = 500+((i+(window_size/2)/(cds_stop-cds_start))*5)
				if i > cds_stop:
					per = 1000+((i+(window_size/2)/(seqlen-cds_stop))*5)
				if per not in mfe_dict:
					mfe_dict[per] = [abs(mfe)]
				else:
					mfe_dict[per].append(abs(mfe))
		for per in mfe_dict:
			mfe_dict[per] = sum(mfe_dict[per])/len(mfe_dict[per])
			
			
			
	if plot_gc == True:
		step_size = 10
		window_size = 60
		a_dict = {}
		t_dict = {}
		g_dict = {}
		c_dict = {}
		gc_dict = {}
		sorted_a_dict = collections.OrderedDict()
		sorted_t_dict = collections.OrderedDict()
		sorted_g_dict = collections.OrderedDict()
		sorted_c_dict = collections.OrderedDict()
		sorted_gc_dict = collections.OrderedDict()
		#for i in range(0,1503):
		#	a_dict[i] = [0]
		#	t_dict[i] = [0]
		#	g_dict[i] = [0]
		#	c_dict[i] = [0]
		#	gc_dict[i] = [0]
		for item in traninfo:
			transcript = item[0]
			cds_start = float(item[1])
			cds_stop = float(item[2])
			seq = item[3]
			seqlen = len(seq)
			for i in range(0,len(seq)-(window_size),step_size):
				mid_window = i+(window_size/2)
				a_count = 0.0
				t_count = 0.0
				g_count = 0.0
				c_count = 0.0
				for x in range(i,i+window_size):
					if seq[x] == "A":
						a_count += 1
					elif seq[x] == "T":
						t_count += 1
					elif seq[x] == "G":
						g_count += 1
					elif seq[x] == "C":
						c_count += 1

				gc_count = g_count + c_count
				norm_a = a_count/window_size
				norm_t = t_count/window_size
				norm_g = g_count/window_size
				norm_c = c_count/window_size
				norm_gc = gc_count/window_size

				final_a = norm_a*y_max
				final_t = norm_t*y_max
				final_g = norm_g*y_max
				final_c = norm_c*y_max
				final_gc = norm_gc*y_max
				if mid_window < cds_start:
					per = mid_window/cds_start
					per = per*500
					per = int(per)
				if mid_window >= cds_start and mid_window <= cds_stop:
					per = (mid_window-cds_start)/(cds_stop-cds_start)
					per = per*500
					per += 500 
					per = int(per)
				if mid_window > cds_stop:
					per = (mid_window-cds_stop)/(seqlen-cds_stop)
					per = per*500
					per += 1000
					per = int(per)
				if per not in a_dict:
					a_dict[per] = [final_a]
				else:
					a_dict[per].append(final_a)
					
				if per not in t_dict:
					t_dict[per] = [final_t]
				else:
					t_dict[per].append(final_t)
				
				if per not in g_dict:
					g_dict[per] = [final_g]
				else:
					g_dict[per].append(final_g)
				
				if per not in c_dict:
					c_dict[per] = [final_c]
				else:
					c_dict[per].append(final_c)
					
				if per not in gc_dict:
					gc_dict[per] = [final_gc]
				else:
					gc_dict[per].append(final_gc)
					
		for per in a_dict:
			a_dict[per] = sum(a_dict[per])/len(a_dict[per])
		for per in t_dict:
			t_dict[per] = sum(t_dict[per])/len(t_dict[per])
		for per in g_dict:
			g_dict[per] = sum(g_dict[per])/len(g_dict[per])
		for per in c_dict:
			c_dict[per] = sum(c_dict[per])/len(c_dict[per])
		for per in sorted(gc_dict.keys()):
			sorted_gc_dict[per] = sum(gc_dict[per])/len(gc_dict[per])
			
			
		a_plot = ax_main.plot(a_dict.keys(), a_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][0], linewidth=4)
		t_plot = ax_main.plot(t_dict.keys(), t_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][1], linewidth=4)
		g_plot = ax_main.plot(g_dict.keys(), g_dict.values(), alpha=0.01, label = labels, zorder=1, color=color_dict['frames'][2], linewidth=4)
		c_plot = ax_main.plot(c_dict.keys(), c_dict.values(), alpha=0.01, label = labels, zorder=1, color='#ffff99', linewidth=4)
		gc_plot = ax_main.plot(sorted_gc_dict.keys(), sorted_gc_dict.values(), alpha=1, label = labels, zorder=1, color='grey', linewidth=4)
		if plot_mfe == True:
			mfe_plot = ax_main.plot(mfe_dict.keys(), mfe_dict.values(), alpha=0.01, label = labels, zorder=1,color='#df8500', linewidth=4)
		#for item,lbl,viz in [(a_plot,"A%",False),(t_plot,"T%",False),(g_plot,"G%",False),(c_plot,"C%",False),(gc_plot,"GC%",True),(mfe_plot,"MFE",False)]:
		for item,lbl,viz in [(a_plot,"A%",False),(t_plot,"T%",False),(g_plot,"G%",False),(c_plot,"C%",False),(gc_plot,"GC%",True)]:
			line_collections.append(item)
			labels.append(lbl)
			start_visible.append(viz)

	leg_offset = (30-17)*5
	if leg_offset <0:
		leg_offset = 0

	ilp = InteractiveLegendPlugin(line_collections, labels, alpha_unsel=0,alpha_sel=0.85,start_visible=start_visible,fontsize=30,xoffset=leg_offset)
	htmllabels = {1:[],2:[],3:[]}


	plugins.connect(fig, ilp, TopToolbar(yoffset=420,xoffset=130))

	ax_main.set_facecolor(background_col)
	# This changes the size of the tick markers, works on both firefox and chrome.
	ax_main.tick_params('both', labelsize=marker_size)
	ax_main.xaxis.set_major_locator(plt.MaxNLocator(3))
	ax_main.yaxis.set_major_locator(plt.MaxNLocator(3))
	ax_main.grid(True, color="white", linewidth=30,linestyle="solid")
	ax_main.text(500,y_max*0.97,"CDS start",fontsize=18,color="black",ha="center")
	ax_main.text(1000,y_max*0.97,"CDS stop",fontsize=18,color="black",ha="center")
	#hide x axis set_ticks
	ax_main.axes.get_xaxis().set_ticks([])
	#Without this style tag the markers sizes will appear correct on browser but be original size when downloaded via png
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph +=  mpld3.fig_to_html(fig)
	return graph


def nuc_comp_scatter(master_dict,filename,title_size,axis_label_size,marker_size,nucleotide,short_code):
	x_values = []
	gc_list = master_dict[1]["lengths"]
	tran_list = master_dict[1]["trans"]
	for i in range(1,len(gc_list)+1):
		x_values.append(i)
	source = ColumnDataSource({'x': x_values,'y':gc_list,'trans':tran_list})	
	x_num = len(gc_list)	
		
	x_values2 = []
	gc_list2 = master_dict[2]["lengths"]
	tran_list2 = master_dict[2]["trans"]
	for i in range(1,len(gc_list2)+1):
		x_values2.append(x_num+i)
	source2 = ColumnDataSource({'x': x_values2,'y':gc_list2,'trans':tran_list2})
	x_num += len(gc_list2)
	
	x_values3 = []
	gc_list3 = master_dict[3]["lengths"]
	tran_list3 = master_dict[3]["trans"]
	for i in range(1,len(gc_list3)+1):
		x_values3.append(x_num+i)
	source3 = ColumnDataSource({'x': x_values3,'y':gc_list3,'trans':tran_list3})
	x_num += len(gc_list3)
	
	x_values4 = []
	gc_list4 = master_dict[4]["lengths"]
	tran_list4 = master_dict[4]["trans"]
	for i in range(1,len(gc_list4)+1):
		x_values4.append(x_num+i)
	source4 = ColumnDataSource({'x': x_values4,'y':gc_list4,'trans':tran_list4})
	x_num += len(gc_list4)	
	full_title = "{}% ({})".format(nucleotide,short_code)
	y_lab = "{} %".format(nucleotide)

	p = figure(plot_width=1300, plot_height=1300,x_axis_label="", title=full_title,y_axis_label=y_lab,toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	p.title.align="center"
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	#p.title.text_font_size = title_size
	#p.xaxis.axis_label_text_font_size = axis_label_size
	#p.xaxis.major_label_text_font_size = marker_size
	#p.yaxis.axis_label_text_font_size = axis_label_size
	#p.yaxis.major_label_text_font_size = marker_size
	#p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='green')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source2,fill_color='red')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source3,fill_color='blue')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source4,fill_color='yellow')
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	hover.tooltips = [("GC%", "@y"),("Count","@x"),("Transcript","@trans")]
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph


def lengths_scatter(master_dict,filename,title_size,axis_label_size,marker_size,short_code):
	x_values = []
	gc_list = master_dict[1]["lengths"]
	tran_list = master_dict[1]["trans"]
	for i in range(1,len(gc_list)+1):
		x_values.append(i)
	source = ColumnDataSource({'x': x_values,'y':gc_list,'trans':tran_list})	
	x_num = len(gc_list)	
		
	x_values2 = []
	gc_list2 = master_dict[2]["lengths"]
	tran_list2 = master_dict[2]["trans"]
	for i in range(1,len(gc_list2)+1):
		x_values2.append(x_num+i)
	source2 = ColumnDataSource({'x': x_values2,'y':gc_list2,'trans':tran_list2})
	x_num += len(gc_list2)
	
	x_values3 = []
	gc_list3 = master_dict[3]["lengths"]
	tran_list3 = master_dict[3]["trans"]
	for i in range(1,len(gc_list3)+1):
		x_values3.append(x_num+i)
	source3 = ColumnDataSource({'x': x_values3,'y':gc_list3,'trans':tran_list3})
	x_num += len(gc_list3)
	
	x_values4 = []
	gc_list4 = master_dict[4]["lengths"]
	tran_list4 = master_dict[4]["trans"]
	for i in range(1,len(gc_list4)+1):
		x_values4.append(x_num+i)
	source4 = ColumnDataSource({'x': x_values4,'y':gc_list4,'trans':tran_list4})
	x_num += len(gc_list4)	
	full_title = "Lengths ({})".format(short_code)
	y_lab = "Length"

	p = figure(plot_width=1300, plot_height=1300,x_axis_label="", title=full_title,y_axis_label=y_lab,toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	p.title.align="center"
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	#p.title.text_font_size = title_size
	#p.xaxis.axis_label_text_font_size = axis_label_size
	#p.xaxis.major_label_text_font_size = marker_size
	#p.yaxis.axis_label_text_font_size = axis_label_size
	#p.yaxis.major_label_text_font_size = marker_size
	#p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='green')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source2,fill_color='red')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source3,fill_color='blue')
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source4,fill_color='yellow')
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	hover.tooltips = [("GC%", "@y"),("Count","@x"),("Transcript","@trans")]
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph



def nuc_comp_box(master_dict,filename,nucleotide,title_size,box_colour,axis_label_size,marker_size,short_code):
	gc_list = master_dict[1]["gc"]
	gc_list2 = master_dict[2]["gc"]
	gc_list3 = master_dict[3]["gc"]
	gc_list4 = master_dict[4]["gc"]
	
	a_list = master_dict[1]["a"]
	a_list2 = master_dict[2]["a"]
	a_list3 = master_dict[3]["a"]
	a_list4 = master_dict[4]["a"]
	
	t_list = master_dict[1]["t"]
	t_list2 = master_dict[2]["t"]
	t_list3 = master_dict[3]["t"]
	t_list4 = master_dict[4]["t"]
	
	g_list = master_dict[1]["g"]
	g_list2 = master_dict[2]["g"]
	g_list3 = master_dict[3]["g"]
	g_list4 = master_dict[4]["g"]
	
	c_list = master_dict[1]["c"]
	c_list2 = master_dict[2]["c"]
	c_list3 = master_dict[3]["c"]
	c_list4 = master_dict[4]["c"]
	
	cats = []
	grouplist = []
	gclist = []
	
	
	if nucleotide == "A":
		cats.append("Group 1")
		for item in a_list:
			grouplist.append("Group 1")
			gclist.append(item)
		if a_list2 != []:
			cats.append("Group 2")
			for item in a_list2:
				grouplist.append("Group 2")
				gclist.append(item)
		if a_list3 != []:
			cats.append("Group 3")
			for item in a_list3:
				grouplist.append("Group 3")
				gclist.append(item)
		if a_list4 != []:
			cats.append("Group 4")
			for item in a_list4:
				grouplist.append("Group 4")
				gclist.append(item)
	if nucleotide == "C":
		cats.append("Group 1")
		for item in c_list:
			grouplist.append("Group 1")
			gclist.append(item)
		if c_list2 != []:
			cats.append("Group 2")
			for item in c_list2:
				grouplist.append("Group 2")
				gclist.append(item)
		if c_list3 != []:
			cats.append("Group 3")
			for item in c_list3:
				grouplist.append("Group 3")
				gclist.append(item)
		if c_list4 != []:
			cats.append("Group 4")
			for item in c_list4:
				grouplist.append("Group 4")
				gclist.append(item)
	if nucleotide == "GC":
		cats.append("Group 1")
		for item in gc_list:
			grouplist.append("Group 1")
			gclist.append(item)
		if gc_list2 != []:
			cats.append("Group 2")
			for item in gc_list2:
				grouplist.append("Group 2")
				gclist.append(item)
		if gc_list3 != []:
			cats.append("Group 3")
			for item in gc_list3:
				grouplist.append("Group 3")
				gclist.append(item)
		if gc_list4 != []:
			cats.append("Group 4")
			for item in gc_list4:
				grouplist.append("Group 4")
				gclist.append(item)
	if nucleotide == "G":
		cats.append("Group 1")
		for item in g_list:
			grouplist.append("Group 1")
			gclist.append(item)
		if g_list2 != []:
			cats.append("Group 2")
			for item in g_list2:
				grouplist.append("Group 2")
				gclist.append(item)
		if g_list3 != []:
			cats.append("Group 3")
			for item in g_list3:
				grouplist.append("Group 3")
				gclist.append(item)
		if g_list4 != []:
			cats.append("Group 4")
			for item in g_list4:
				grouplist.append("Group 4")
				gclist.append(item)	
	if nucleotide == "T":
		cats.append("Group 1")
		for item in t_list:
			grouplist.append("Group 1")
			gclist.append(item)
		if t_list2 != []:
			cats.append("Group 2")
			for item in t_list2:
				grouplist.append("Group 2")
				gclist.append(item)
		if t_list3 != []:
			cats.append("Group 3")
			for item in t_list3:
				grouplist.append("Group 3")
				gclist.append(item)
		if t_list4 != []:
			cats.append("Group 4")
			for item in t_list4:
				grouplist.append("Group 4")
				gclist.append(item)
	df = pd.DataFrame({"group":grouplist,"gc":gclist})
	groups = df.groupby('group')
	q1 = groups.quantile(q=0.25)
	q2 = groups.quantile(q=0.5)
	q3 = groups.quantile(q=0.75)
	iqr = q3 - q1
	upper = q3 + 1.5*iqr
	lower = q1 - 1.5*iqr

	# find the outliers for each category
	def outliers(group):
		cat = group.name
		return group[(group.gc > upper.loc[cat]['gc']) | (group.gc < lower.loc[cat]['gc'])]['gc']
	out = groups.apply(outliers).dropna()
		
	# find the outliers for each category
	def outliers(group):
		cat = group.name
		return group[(group.gc > upper.loc[cat]['gc']) | (group.gc < lower.loc[cat]['gc'])]['gc']
			
	# prepare outlier data for plotting, we need coordinates for every outlier.
	if not out.empty:
		outx = []
		outy = []
		for keys in out.index:
			try:
				outx.append(keys[0])
				outy.append(out.loc[keys[0]].loc[keys[1]])
			except:
				pass
	full_title = "{}% ({})".format(nucleotide,short_code)
	y_lab = '{} %'.format(nucleotide)
	p = figure(plot_width=1300, plot_height=1300,tools="reset,pan,box_zoom,save,hover,tap", background_fill_color="#efefef", x_range=cats, toolbar_location="below",title=full_title,y_axis_label=y_lab)
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	#p.background_fill_color = background_color
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"
	# if no outliers, shrink gcs of stems to be no longer than the minimums or maximums
	qmin = groups.quantile(q=0.00)
	qmax = groups.quantile(q=1.00)
	upper.gc = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'gc']),upper.gc)]
	lower.gc = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'gc']),lower.gc)]

	# stems
	p.segment(cats, upper.gc, cats, q3.gc, line_color="black")
	p.segment(cats, lower.gc, cats, q1.gc, line_color="black")

	# boxes
	p.vbar(cats, 0.7, q2.gc, q3.gc, fill_color=box_colour, line_color="black")
	p.vbar(cats, 0.7, q1.gc, q2.gc, fill_color=box_colour, line_color="black")

	# whiskers (almost-0 height rects simpler than segments)
	p.rect(cats, lower.gc, 0.2, 0.01, line_color="black")
	p.rect(cats, upper.gc, 0.2, 0.01, line_color="black")
	if not out.empty:
		p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)

	p.xgrid.grid_line_color = None
	p.ygrid.grid_line_color = "white"
	p.grid.grid_line_width = 2
	p.xaxis.major_label_text_font_size="12pt"
	
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph


def lengths_box(master_dict,filename,box_colour,short_code,title_size,marker_size, axis_label_size):
	gc_list = master_dict[1]["lengths"]
	gc_list2 = master_dict[2]["lengths"]
	gc_list3 = master_dict[3]["lengths"]
	gc_list4 = master_dict[4]["lengths"]
	cats = []
	grouplist = []
	gclist = []
		
	cats.append("Lengths_G1")
	for item in gc_list:
		grouplist.append("Lengths_G1")
		gclist.append(item)
	if gc_list2 != []:
		cats.append("Lengths_G2")
		for item in gc_list2:
			grouplist.append("Lengths_G2")
			gclist.append(item)
	if gc_list3 != []:
		cats.append("Lengths_G3")
		for item in gc_list3:
			grouplist.append("Lengths_G3")
			gclist.append(item)
	if gc_list4 != []:
		cats.append("Lengths_G4")
		for item in gc_list4:
			grouplist.append("Lengths_G4")
			gclist.append(item)

	df = pd.DataFrame({"group":grouplist,"lengths":gclist})
	groups = df.groupby('group')
	q1 = groups.quantile(q=0.25)
	q2 = groups.quantile(q=0.5)
	q3 = groups.quantile(q=0.75)
	iqr = q3 - q1
	upper = q3 + 1.5*iqr
	lower = q1 - 1.5*iqr

	# find the outliers for each category
	def outliers(group):
		cat = group.name
		return group[(group.lengths > upper.loc[cat]['lengths']) | (group.lengths < lower.loc[cat]['lengths'])]['lengths']
	out = groups.apply(outliers).dropna()
		
	# find the outliers for each category
	def outliers(group):
		cat = group.name
		return group[(group.lengths > upper.loc[cat]['lengths']) | (group.lengths < lower.loc[cat]['lengths'])]['lengths']
			
	# prepare outlier data for plotting, we need coordinates for every outlier.
	if not out.empty:
		outx = []
		outy = []
		for keys in out.index:
			outx.append(keys[0])
			outy.append(out.loc[keys[0]].loc[keys[1]])
		
	full_title = "Lengths ({})".format(short_code)
	p = figure(plot_width=1300, plot_height=1300,tools="reset,pan,box_zoom,save,hover,tap",title=full_title, background_fill_color="#efefef", x_range=cats, toolbar_location="below")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	# if no outliers, shrink gcs of stems to be no longer than the minimums or maximums
	qmin = groups.quantile(q=0.00)
	qmax = groups.quantile(q=1.00)
	upper.lengths = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'lengths']),upper.lengths)]
	lower.lengths = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'lengths']),lower.lengths)]

	# stems
	p.segment(cats, upper.lengths, cats, q3.lengths, line_color="black")
	p.segment(cats, lower.lengths, cats, q1.lengths, line_color="black")

	# boxes
	p.vbar(cats, 0.7, q2.lengths, q3.lengths, fill_color=box_colour, line_color="black")
	p.vbar(cats, 0.7, q1.lengths, q2.lengths, fill_color=box_colour, line_color="black")

	# whiskers (almost-0 height rects simpler than segments)
	p.rect(cats, lower.lengths, 0.2, 0.01, line_color="black")
	p.rect(cats, upper.lengths, 0.2, 0.01, line_color="black")
	if not out.empty:
		p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)

	p.xgrid.grid_line_color = None
	p.ygrid.grid_line_color = "white"
	p.grid.grid_line_width = 2
	p.xaxis.major_label_text_font_size="12pt"
	
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph




def codon_usage(codon_dict,short_code,title_size, axis_label_size, marker_size,filename):
	allxvals = []
	allyvals = []
	alllabels = []
	amino_acids = []
	aa_dict = {"TTT":"Phenylalanine", "TTC":"Phenylalanine", "TTA":"Leucine", "TTG":"Leucine",
    "TCT":"Serine", "TCC":"Serine", "TCA":"Serine", "TCG":"Serine",
    "TAT":"Tyrosine", "TAC":"Tyrosine", "TAA":"*", "TAG":"*",
    "TGT":"Cysteine", "TGC":"Cysteine", "TGA":"*", "TGG":"Tryptophan",
    "CTT":"Leucine", "CTC":"Leucine", "CTA":"Leucine", "CTG":"Leucine",

    "CCT":"Proline", "CCC":"Proline", "CCA":"Proline", "CCG":"Proline",
    "CAT":"Histidine", "CAC":"Histidine", "CAA":"Glutamine", "CAG":"Glutamine",
    "CGT":"Arginine", "CGC":"Arginine", "CGA":"Arginine", "CGG":"Arginine",
    "ATT":"Isoleucine", "ATC":"Isoleucine", "ATA":"Isoleucine", "ATG":"Methionine",

    "ACT":"Threonine", "ACC":"Threonine", "ACA":"Threonine", "ACG":"Threonine",
    "AAT":"Asparagine", "AAC":"Asparagine", "AAA":"Lysine", "AAG":"Lysine",
    "AGT":"Serine", "AGC":"Serine", "AGA":"Arginine", "AGG":"Arginine",
    "GTT":"Valine", "GTC":"Valine", "GTA":"Valine", "GTG":"Valine",

    "GCT":"Alanine", "GCC":"Alanine", "GCA":"Alanine", "GCG":"Alanine",
    "GAT":"Aspartic Acid", "GAC":"Aspartic Acid", "GAA":"Glutamic Acid", "GAG":"Glutamic Acid",
    "GGT":"Glycine", "GGC":"Glycine", "GGA":"Glycine", "GGG":"Glycine"}

	codon_list = ["ATG","TTT","TTC","CTT", "CTC", "CTA", "CTG","TTA","TTG", "AGT", "AGC","TCT", "TCC", "TCA", "TCG","TAT", "TAC","TGT", "TGC","TGG","CCT", "CCC", "CCA",
			   "CCG","CAT", "CAC", "CAA", "CAG","AGA", "AGG","CGT", "CGC", "CGA", "CGG","ATT", "ATC", "ATA","ACT", "ACC", "ACA", "ACG","AAT", "AAC", "AAA", "AAG","GTT",
			   "GTC", "GTA", "GTG","GCT", "GCC", "GCA", "GCG","GAT", "GAC", "GAA", "GAG","GGT", "GGC", "GGA", "GGG","TAG","TAA","TGA"]

	curr_count = 0
	for codon in codon_list:
		curr_count+=1
		allxvals.append(curr_count)
		if codon in codon_dict:
			allyvals.append(codon_dict[codon])
		else:
			allyvals.append(0)
		alllabels.append(codon)
		amino_acids.append(aa_dict[codon])
	full_title = "Codon usage ({})".format(short_code)
	x_lab = ''
	y_lab = 'Count'
	min_y = min(0,min(allyvals))-.02
	max_y = max(allyvals)*1.05
	p = figure(plot_width=1300, plot_height=1300,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap,save",y_range=(min_y,max_y))
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = "14pt"
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	#p.background_fill_color = background_color
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"
	color_palette_list = []
	colormap =cm.get_cmap("gist_rainbow") #choose any matplotlib colormap here
	start_val = 0.75
	for i in range(0,21):
		start_val -= 0.0293
		rgb = colormap(start_val)[:3]
		hexval = matplotlib.colors.rgb2hex(rgb)
		color_palette_list.append(hexval)

	color_map = bmo.CategoricalColorMapper(factors=["Phenylalanine","Leucine","Serine","Tyrosine","*","Cysteine","Tryptophan","Proline","Histidine","Glutamine","Arginine","Isoleucine","Methionine","Threonine","Asparagine","Lysine","Valine","Alanine","Aspartic Acid","Glutamic Acid","Glycine"],
									palette=color_palette_list)
	p.quad(top=[max_y,max_y,max_y,max_y,max_y,max_y,max_y,max_y,max_y,max_y,max_y],
				bottom=[min_y,min_y,min_y,min_y,min_y,min_y,min_y,min_y,min_y,min_y,min_y],
				left=[0.5,3.5,15.5,19.5,24.5,28.5,37.5,43.5,49.5,55.5,61.5],
				right=[1.5,9.5,17.5,20.5,26.5,34.5,41.5,45.5,53.5,57.5,64.5],
				color="#e0e0e0")
	source = ColumnDataSource({'x': allxvals,'y':allyvals,'labels':alllabels,'amino_acids':amino_acids})
	p.scatter('x','y',source=source, alpha=1,color={'field': 'amino_acids', 'transform': color_map},size=16,line_color="black")
	p.xaxis.ticker = [1,2.5,6.5,12.5,16.5,18.5,20,22.5,25.5,27.5,31.5,36,39.5,42.5,44.5,47.5,51.5,54.5,56.5,59.5,63]
	p.xaxis.major_label_overrides = {1:"Met",2.5:"Phe",6.5:"Leu",12.5:"Ser",16.5:"Tyr",18.5:"Cys",20:"Trp",22.5:"Pro",25.5:"His",27.5:"Gln",31.5:"Arg",
								 36:"Ile",39.5:"Thr",42.5:"Asn",44.5:"Lys",47.5:"Val",51.5:"Ala",54.5:"Asp",56.5:"Glu",59.5:"Gly",63:"Stop"}
	#p.vbar(x=[1], width=1,bottom=0,color="gray",top=[1])

	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Count", "@y"),("Codon","@labels"),("Amino acid","@amino_acids")]

	output_file("scatter10k.html", title="Codon usage")
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(filename)
	graph += file_html(p,CDN)
	return graph




def make_autopct(values):
	def my_autopct(pct):
		total = sum(values)
		val = int(round(pct*total/100.0))
		return '{p:.2f}% \nCount: ({v:d})'.format(p=pct,v=val)
	return my_autopct

def gene_count(short_code, background_col, title_size, axis_label_size, subheading_size, marker_size,coding,noncoding):
	labels = ["","Genes","Transcripts",""]
	fig, ax = plt.subplots(figsize=(13,8))
	N = 4
	bar_width = 0.35
	bar_l = [i for i in range(N)]
	tick_pos = [i+(bar_width/2) for i in bar_l]
	ind = np.arange(N)    # the x locations for the groups
	all_reads_count = 0
	title_str = "Reads breakdown ({})".format(short_code)
	plt.title(title_str,fontsize=title_size)
	plt.xticks(tick_pos, labels)
	#plt.yticks(np.arange(0, 81, 10))
	#plt.legend((p1[0], p2[0], p3[0],p4[0],p5[0],p6[0]), ('Cutadapt removed', 'rRNA removed', 'Unmapped', 'Ambiguous','Mapped noncoding', 'Mapped coding'))

	ax.set_facecolor(background_col)
	ax.tick_params('y', labelsize=marker_size)
	if len(labels) > 12:
		marker_size = int(marker_size/(len(labels)/8))
	ax.tick_params('x', labelsize=marker_size)


	plt.grid(color="white", linewidth=2,linestyle="solid")
	totals = []
	for i in range(0,N):
		curr_total = 0
		curr_total += coding[i]
		curr_total += noncoding[i]
		if curr_total > 0:
			all_reads_count += curr_total
			totals.append(float(curr_total))
		else:
			totals.append(1)

	for i in range(0,len(coding)):
		per = (coding[i]/totals[i])*100

	for i in range(0,len(noncoding)):
		per = (noncoding[i]/totals[i])*100

	plt.ylabel('Count',fontsize=axis_label_size,labelpad=100)
	p1 = plt.bar(ind, coding, bar_width, color='#80ff80',linewidth=0)
	p2 = plt.bar(ind, noncoding, bar_width,color='#ff7275', bottom=coding,linewidth=0)
	p8 = plt.bar(ind, totals, bar_width, color='#5e0003',linewidth=0,alpha=0)

	#Dummy plot point so we can add total reads to the legend
	total_reads_plot = plt.plot(0,0,alpha=0)
	plt.legend(
			(p2[0],p1[0]),
			(	'Non Coding: {:,}'.format(sum(noncoding)),
				'Coding {:,}'.format(sum(coding))
				))
	for i, bar in enumerate(p8.get_children()):
		lab1 = ""
		per = int(round((noncoding[i]/totals[i])*100,0))
		lab1 += '<br><b>Non-coding:</b> {:,}  ({}%)'.format(noncoding[i], per)
		per = int(round((coding[i]/totals[i])*100,0))
		lab1 += '<br><b>Coding:</b> {:,}  ({}%)'.format(coding[i], per)
		tooltip1 = plugins.LineHTMLTooltip(bar, lab1,voffset=10, hoffset=30,css=line_tooltip_css)
		plugins.connect(fig, tooltip1)
	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph






















def letterAt(letter, x, y, yscale=1, ax=None):
	fp = FontProperties(family="Arial", weight="bold")
	LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
				"G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
				"A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
				"C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
	text = LETTERS[letter]

	t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
		mpl.transforms.Affine2D().translate(x,y) + ax.transData
	p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
	if ax != None:
		ax.add_artist(p)




def ticker(tick):
	return "{:.0f} + {:.2f}".format(tick, tick % 1)



