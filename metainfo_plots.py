
import sys
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import os
import re
import subprocess
import shelve
import mpld3
import operator
import logging
from sqlitedict import SqliteDict
from math import log
from mpld3 import plugins,utils
import collections
from mpld3.utils import get_id
import pandas as pd
import numpy as np
from flask import make_response
from new_plugins import InteractiveLegendPlugin,PointHTMLTooltip,TopToolbar,DownloadProfile,DownloadPNG
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





def mismatches(master_dict, title, short_code, background_col,title_size, axis_label_size, subheading_size,marker_size):
	fig, ax = plt.subplots( figsize=(13,8))
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Position', fontsize="26")
	ax.set_ylabel('Count',fontsize="26",labelpad=50)
	if master_dict.values():
		ax.set_ylim(0,max(master_dict.values())*1.25)
	else:
		ax.set_ylim(0,1)
	width = 0.90
	#plot it
	ax = plt.subplot(111)
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str,y=1.05,fontsize=title_size)
	ax.bar(master_dict.keys(), master_dict.values(), width, color="green", linewidth=0, align="center")
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=18)
	ax.set_ylabel('Count', labelpad=100,fontsize="26")
	ax.set_xlabel('Readlength',fontsize="26")
	ax.tick_params('both', labelsize=marker_size)
	#ax.xaxis.set_major_locator(plt.MaxNLocator(3))
	#ax.yaxis.set_major_locator(plt.MaxNLocator(3))
	#plt.rc('axes', linewidth=40,edgecolor="green")
	plt.grid(color="white", linewidth=2,linestyle="solid")
	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph

def calc_factor(master_dict):
	maxval = max(master_dict.values())
	string_maxval = str(maxval)
	zeroes = len(string_maxval)-1
	zeroes_string = "0"*zeroes
	factor = int("1"+zeroes_string)
	for key in master_dict:
		master_dict[key] = float(master_dict[key])/factor
	return master_dict, zeroes
	



def readlen_dist(master_dict,title,short_code,background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size):
	master_dict, factor = calc_factor(master_dict)
	returnstr = "Readlen,Count\n"
	for key in master_dict:
		returnstr += "{},{}\n".format(key, master_dict[key])
	fig, ax = plt.subplots( figsize=(13,8))
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Read Length', fontsize="26")
	ax.set_ylabel('Count',fontsize="26",labelpad=50)
	if master_dict.values():
		ax.set_ylim(0,max(master_dict.values())*1.25)
	else:
		ax.set_ylim(0,1)
	width = 0.90
	#plot it
	ax = plt.subplot(111)
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str,y=1.05,fontsize=title_size)
	#logging.warn("Width is ", width)
	read_length_list = test_list = [int(i) for i in master_dict.keys()]
	ax.bar(read_length_list, master_dict.values(), width, color=readlength_col, linewidth=0, align="center")
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	#ax.set_ylabel('Count', labelpad=100,fontsize=axis_label_size)
	ax.set_ylabel('Count (x10 {})'.format(factor), labelpad=100,fontsize=axis_label_size)
	ax.set_xlabel('Readlength',labelpad=-15,fontsize=axis_label_size)
	#ax.xaxis.set_major_locator(plt.MaxNLocator(3))
	#ax.yaxis.set_major_locator(plt.MaxNLocator(3))
	#plt.rc('axes', linewidth=40,edgecolor="green")
	plt.grid(color="white", linewidth=2,linestyle="solid")
	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph



def mismatch_pos(master_dict,title,short_code,background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size):
	fig, ax = plt.subplots( figsize=(13,8))
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Read Length', fontsize="26")
	ax.set_ylabel('Count',fontsize="26",labelpad=50)
	
	if master_dict.values():
		ax.set_ylim(0,max(master_dict.values())*1.25)
	else:
		ax.set_ylim(0,1)
	width = 0.90
	#plot it
	ax = plt.subplot(111)
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str,y=1.05,fontsize=title_size)
	ax.bar(master_dict.keys(), master_dict.values(), width, color=readlength_col, linewidth=0, align="center")
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	ax.set_ylabel('Count', labelpad=100,fontsize=axis_label_size)
	ax.set_xlabel('Position',fontsize=axis_label_size)
	ax.set_xlim(-1,max(master_dict.keys())+1)
	#ax.xaxis.set_major_locator(plt.MaxNLocator(3))
	#ax.yaxis.set_major_locator(plt.MaxNLocator(3))
	#plt.rc('axes', linewidth=40,edgecolor="green")
	plt.grid(color="white", linewidth=2,linestyle="solid")
	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph



def nuc_comp(master_dict,maxreadlen,title, nuc_comp_type,nuc_comp_direction,short_code,background_col,a_col,t_col,g_col,c_col,title_size, axis_label_size, subheading_size,marker_size,legend_size):

	labels = ["A","T","G","C"]
	returnstr = "Position,A,T,G,C\n"
	
	for i in range(0,maxreadlen):
		returnstr += "{},{:.2f},{:.2f},{:.2f},{:.2f}\n".format(i,master_dict["A"][i],master_dict["T"][i],master_dict["G"][i],master_dict["C"][i])
	fig, ax = plt.subplots( figsize=(13,8))
	ax = plt.subplot(111)
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Position (nucleotides)',labelpad=-10,fontsize=axis_label_size)
	if nuc_comp_type == "nuc_comp_per":
		ax.set_ylim(0,100)
		ax.set_ylabel('Percent %',fontsize=axis_label_size,labelpad=50)
	elif nuc_comp_type == "nuc_comp_count":
		maxheight = max(max(master_dict["A"].values()),max(master_dict["T"].values()),max(master_dict["G"].values()),max(master_dict["C"].values()))
		ax.set_ylim(0,maxheight)
		ax.set_ylabel('Count',fontsize=axis_label_size,labelpad=0)
	if nuc_comp_direction == "nuc_comp_five":
		ax.set_xlim(0,maxreadlen)
	elif nuc_comp_direction == "nuc_comp_three":
		ax.set_xlim((maxreadlen*-1),-1)
	width = 0.95
	#plot it
	
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str, y=1.05,fontsize=title_size)
	a_line = ax.plot(master_dict["A"].keys(), master_dict["A"].values(), label=labels, color=a_col, linewidth=6)
	t_line = ax.plot(master_dict["T"].keys(), master_dict["T"].values(), label=labels, color=t_col, linewidth=6)
	g_line = ax.plot(master_dict["G"].keys(), master_dict["G"].values(), label=labels, color=g_col, linewidth=6)
	c_line = ax.plot(master_dict["C"].keys(), master_dict["C"].values(), label=labels, color=c_col, linewidth=6)
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")
	ilp = InteractiveLegendPlugin([a_line, t_line, g_line, c_line], ["A","T","G","C"], alpha_unsel=0,alpha_sel=1,start_visible=True,fontsize=legend_size)
	plugins.connect(fig, ilp,TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	#return graph
	return graph




def mrna_dist_readlen(mrna_dist_dict,mrna_readlen_per,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,legend_size):

	
	if mrna_readlen_per == False:
		mrna_dist_dict, factor = calc_mrnadist_factor(mrna_dist_dict)
	returnstr = "Readlength,5_leader,start_codon,cds,stop_codon,3_trailer\n"
	
	
	labels = ["5_leader","start_codon","cds","stop_codon","3_trailer"]
	
	for i in range(0,100):
		returnstr += "{},".format(i)
		for label in labels:
			if i in mrna_dist_dict[label]:
				returnstr += "{},".format(mrna_dist_dict[label][i])
			else:
				returnstr += "0,"
		returnstr += "\n"
	fig, ax = plt.subplots( figsize=(13,8))
	ax = plt.subplot(111)
	#rects1 = ax.bar([20,21,22,23,24,25,26,27,28], [100,200,100,200,100,200,100,200,100], 0.1, color='r',align='center')
	ax.set_xlabel('Readlength',fontsize=axis_label_size,labelpad=-10)

	maxheight = 0
	for pos in mrna_dist_dict:
		for readlen in mrna_dist_dict[pos]:
			if mrna_dist_dict[pos][readlen] > maxheight:
				maxheight = mrna_dist_dict[pos][readlen]

	if mrna_readlen_per == False:
		ax.set_ylim(0,maxheight*1.1)
		ax.set_ylabel('Count (x 10 {})'.format(factor),fontsize=axis_label_size,labelpad=100)
	elif mrna_readlen_per == True:
		ax.set_ylim(0,100)
		ax.set_ylabel('Percent %',fontsize=axis_label_size,labelpad=100)
	ax.set_xlim(15,100)
	width = 0.95
	#plot it
	
	title_str = "mRNA distribution vs Readlengths ({})".format(short_code)
	ax.set_title(title_str, y=1.05,fontsize=title_size)
	five_leader_line = ax.plot(mrna_dist_dict["5_leader"].keys(), mrna_dist_dict["5_leader"].values(),label=labels, color='#ff6d6d', linewidth=4)
	start_codon_line = ax.plot(mrna_dist_dict["start_codon"].keys(), mrna_dist_dict["start_codon"].values(),label=labels, color='#a38b22', linewidth=4)
	cds_line = ax.plot(mrna_dist_dict["cds"].keys(), mrna_dist_dict["cds"].values(),label=labels, color='#00d8cd', linewidth=4)
	stop_codon_line = ax.plot(mrna_dist_dict["stop_codon"].keys(), mrna_dist_dict["stop_codon"].values(),label=labels, color='#05bc27', linewidth=4)
	three_trailer_line = ax.plot(mrna_dist_dict["3_trailer"].keys(), mrna_dist_dict["3_trailer"].values(),label=labels, color='#ff77d4', linewidth=4)
	#total_line = ax.plot(mrna_dist_dict["total"].keys(), mrna_dist_dict["total"].values(),label=labels, color='grey', linewidth=4)
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")
	ilp = InteractiveLegendPlugin([five_leader_line, start_codon_line, cds_line, stop_codon_line, three_trailer_line], labels, alpha_unsel=0,alpha_sel=1,start_visible=True)
	plugins.connect(fig, ilp,TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph



def dinuc_bias(master_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size):
	master_count_dict,factor = calc_factor(master_count_dict)
	fig, ax = plt.subplots(figsize=(13,8))
	N = 16
	bar_width = 0.35
	bar_l = [i for i in range(16)]
	tick_pos = [i+(bar_width/2) for i in bar_l]
	ind = np.arange(N)    # the x locations for the groups
	p1 = plt.bar(ind, master_count_dict.values(), bar_width, color='#9ACAFF',linewidth=4,edgecolor='#9ACAFF')
	plt.ylabel('Count (x 10 {})'.format(factor),fontsize=axis_label_size,labelpad=100)
	plt.xlabel('Dinculeotide',fontsize=axis_label_size,labelpad=-10)
	title_str = "Dinucleotide composition ({})".format(short_code)
	plt.title(title_str,fontsize=title_size)
	plt.xticks(tick_pos, master_count_dict.keys())
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plugins.connect(fig,TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph


def calc_meta_factor(inlist):
	maxval = max(inlist)
	string_maxval = str(maxval)
	zeroes = len(string_maxval)-1
	zeroes_string = "0"*zeroes
	factor = int("1"+zeroes_string)
	for i in range(0,len(inlist)):
		inlist[i] = float(inlist[i])/factor
	return inlist, zeroes
	



def metagene_plot(readlen_list, fiveprime_list, threeprime_list,metagene_type,title,minreadlen, maxreadlen,short_code,background_col,metagene_fiveprime_col,metagene_threeprime_col,title_size, axis_label_size, subheading_size,marker_size,metagene_end, metagene_aggregate):
	fig, ax = plt.subplots( figsize=(13,8))
	ind = np.array(readlen_list)
	file_colors = ["#FF4A45","#4286f4","#42f450","#f4f142","#ff9e16","#a800aa"]

	if metagene_aggregate == True:	
		width = 0.35       # the width of the bars
		returnstr = "Position,Count\n"
		if metagene_end == "metagene_five":
			fiveprime_list,factor = calc_meta_factor(fiveprime_list)
			rects1 = ax.bar(ind, fiveprime_list, width, color=metagene_fiveprime_col, linewidth=2,edgecolor=metagene_fiveprime_col)
			start_pos = -300
			for count in fiveprime_list:
				returnstr += "{},{}\n".format(start_pos, count)
				start_pos += 1
		if metagene_end == "metagene_three":
			threeprime_list,factor = calc_meta_factor(threeprime_list)
			rects2 = ax.bar(ind + width, threeprime_list, width, color=metagene_threeprime_col, linewidth=2,edgecolor=metagene_threeprime_col)
			start_pos = -300
			for count in threeprime_list:
				returnstr += "{},{}\n".format(start_pos, count)
				start_pos += 1
		#if metagene_five == True and metagene_three == True:
		#	leg = ax.legend((rects1[0], rects2[0]), ('5 prime', '3 prime'),fontsize="26")
		if metagene_end == "metagene_five":
			leg = ax.legend((rects1), ['5 prime'],fontsize="26")
		elif metagene_end == "metagene_three":
			leg = ax.legend((rects2), ['3 prime'],fontsize="26")
		leg.get_frame().set_edgecolor('#D2D2EB')
		
	else:
		line_collections = []
		labels = []
		file_ind = 0
		returnstr = "Position,"
		return_dict = {}
		if metagene_end == "metagene_five":
			
			count_dict = fiveprime_list
			labelend = " 5' ends"
		elif metagene_end == "metagene_three":
			count_dict = threeprime_list
			labelend = " 3' ends"
		#Create the header for the output file
		for file_id in count_dict:
			returnstr += "{} counts,".format(file_id)
		returnstr += "\n"
		
		for file_id in count_dict:
			value_list = []
			pos_list = []
			i = (len(count_dict[file_id])/2)*(-1)
			for val in count_dict[file_id]:
				value_list.append(val)
				pos_list.append(i)
				if i not in return_dict:
					return_dict[i] = ""
				return_dict[i] += str(val)+","
				i+=1
			plotcounts = ax.plot(pos_list, value_list, alpha=0.75, label = file_id, zorder=2, color=file_colors[file_ind], linewidth=2)
			file_ind += 1
			line_collections.append(plotcounts)
			labels.append(file_id+labelend)
		ilp = InteractiveLegendPlugin(line_collections, labels, alpha_unsel=0,alpha_sel=0.85,start_visible=True,fontsize=19,xoffset=200)
		for i in range(-600,600):
			if i in return_dict:
				returnstr += "{},{}\n".format(i,return_dict[i])
			else:
				returnstr += "{},0\n".format(i)
	#ilp = InteractiveLegendPlugin([rects1,rects2], ["5'","3'"], alpha_unsel=0,alpha_sel=0.75)
	# add some text for labels, title and axes ticks
	if metagene_aggregate == True:	
		ax.set_ylabel('Count x (10 x{})'.format(factor), labelpad=25,fontsize=axis_label_size)
	else:
		ax.set_ylabel('Count', labelpad=25,fontsize=axis_label_size)
	if metagene_type == "metagene_start":
		ax.set_xlabel("Position relative to cds start (nucleotides)",labelpad=-10,fontsize=axis_label_size)
	elif metagene_type == "metagene_stop":
		ax.set_xlabel("Position relative to cds stop (nucleotides)",labelpad=-10,fontsize=axis_label_size)
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str, y=1,fontsize=title_size)
	ax.tick_params('both', labelsize=marker_size)
	plt.xlim(minreadlen, maxreadlen)
	
	ax.set_facecolor(background_col)
	if metagene_aggregate == True:
		plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str),DownloadProfile(returnstr=returnstr))
	else:
		plugins.connect(fig, ilp, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str),DownloadProfile(returnstr=returnstr))
	plt.grid(color="white", linewidth=2,linestyle="solid")
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph

def calc_trip_factor(read_dict):
	maxval = 0
	for i in range(0,len(read_dict["frame1"])):
		count_list = [read_dict["frame1"][i],read_dict["frame2"][i],read_dict["frame3"][i]]
		for count in count_list:
			if count > maxval:
				maxval = int(count)
	string_maxval = str(maxval)
	zeroes = len(string_maxval)-1
	zeroes_string = "0"*zeroes
	factor = int("1"+zeroes_string)
	for i in range(0,len(read_dict["frame1"])):
		read_dict["frame1"][i] = float(read_dict["frame1"][i])/factor
		read_dict["frame2"][i] = float(read_dict["frame2"][i])/factor
		read_dict["frame3"][i] = float(read_dict["frame3"][i])/factor
	return read_dict, zeroes




def trip_periodicity_plot(read_dict,title,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,legend_size):
	read_dict, factor = calc_trip_factor(read_dict)
	tot_high_count = 0.01
	tot_low_count = 0.0
	returnstr = "Readlength, Frame 1 Count, Frame 2 count, Frame 3 count\n"
	for i in range(0,len(read_dict["frame1"])):
		count_list = [read_dict["frame1"][i],read_dict["frame2"][i],read_dict["frame3"][i]]
		returnstr += "{},{},{},{}\n".format(read_dict["readlengths"][i],read_dict["frame1"][i],read_dict["frame2"][i],read_dict["frame3"][i])
		sorted_count_list = sorted(count_list)
		tot_high_count += sorted_count_list[2]
		tot_low_count += sorted_count_list[1]

	trip_periodicity_score = round(1-(tot_low_count/tot_high_count),2)
	df = pd.DataFrame(read_dict, columns = ['frame1', 'frame2', 'frame3', 'readlengths'])
	pos = list(range(len(df['frame1'])))
	width = 0.25
	fig, ax = plt.subplots( figsize=(13,8))

	plt.plot(0,0,alpha=0,label="score")
	# Create a bar with frame1 data in position pos,
	plt.bar(pos,df['frame1'], width, alpha=1,color=redhex,linewidth=0,label="frame 0")#df['readlengths'][0])
	# Create a bar with frame2 data in position pos + some width buffer,
	plt.bar([p + width for p in pos],df['frame2'],width,alpha=1,color=greenhex,linewidth=0,label="frame 1")
	# Create a bar with frame3 data in position pos + some width buffer,
	plt.bar([p + width*2 for p in pos], df['frame3'],width,alpha=1,color=bluehex,linewidth=0, label="frame 2")

	# Set the y axis label
	#ax.set_ylabel('Count', labelpad=38,fontsize=axis_label_size)
	ax.set_ylabel('Count (x10 {})'.format(factor), labelpad=100,fontsize=axis_label_size)
	ax.set_xlabel('Readlength',labelpad=-15,fontsize=axis_label_size)
	# Set the chart's title
	title_str = "{} ({})".format(title,short_code)
	ax.set_title(title_str,y=1.05,fontsize=title_size)

	# Set the position of the x ticks
	ax.set_xticks([p + 1.5 * width for p in pos])

	# Set the labels for the x ticks
	ax.set_xticklabels(df['readlengths'])

	# Setting the x-axis and y-axis limits
	plt.xlim(min(pos)-width, max(pos)+width*4)
	plt.ylim([0, max(max(df['frame1']),max(df['frame2']),max(df['frame3']))*1.1] )

	# Adding the legend and showing the plot
	leg = plt.legend(["Score: {}".format(trip_periodicity_score),'Frame 1', 'Frame 2', 'Frame 3'], loc='upper right', fontsize=legend_size)
	leg.get_frame().set_edgecolor('#D2D2EB')
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")
	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str),DownloadProfile(returnstr=returnstr))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph

def make_autopct(values):
	def my_autopct(pct):
		total = sum(values)
		val = int(round(pct*total/100.0))
		return '{p:.2f}% \nCount: ({v:d})'.format(p=pct,v=val)
	return my_autopct




def mapped_reads_plot(unmapped, mapped_coding, mapped_noncoding, labels, ambiguous,cutadapt_removed, rrna_removed,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size, breakdown_per, pcr_duplicates,legend_size):
	fig, ax = plt.subplots(figsize=(13,8))
	N = len(unmapped)
	bar_width = 0.35
	bar_l = [i for i in range(len(unmapped))]
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
	for i in range(0,len(unmapped)):
		curr_total = 0
		curr_total += unmapped[i]
		curr_total += mapped_coding[i]
		curr_total += mapped_noncoding[i]
		curr_total += ambiguous[i]
		curr_total += cutadapt_removed[i]
		curr_total += rrna_removed[i]
		curr_total += pcr_duplicates[i]
		if curr_total > 0:
			all_reads_count += curr_total
			totals.append(float(curr_total))
		else:
			totals.append(1)
	
	cutadapt_removed_per = []
	rrna_removed_per = []
	unmapped_per = []
	ambiguous_per= []
	mapped_noncoding_per = []
	mapped_coding_per = []
	pcr_duplicates_per = []
	totals_per = []
	
	
	for i in range(0,len(cutadapt_removed)):
		per = (cutadapt_removed[i]/totals[i])*100
		cutadapt_removed_per.append(per)
		totals_per.append(100)

	for i in range(0,len(rrna_removed)):
		per = (rrna_removed[i]/totals[i])*100
		rrna_removed_per.append(per)

	for i in range(0,len(unmapped)):
		per = (unmapped[i]/totals[i])*100
		unmapped_per.append(per)

	for i in range(0,len(ambiguous)):
		per = (ambiguous[i]/totals[i])*100
		ambiguous_per.append(per)

	for i in range(0,len(mapped_noncoding)):
		per = (mapped_noncoding[i]/totals[i])*100
		mapped_noncoding_per.append(per)

	for i in range(0,len(mapped_coding)):
		per = (mapped_coding[i]/totals[i])*100
		mapped_coding_per.append(per)

	for i in range(0,len(pcr_duplicates)):
		per = (pcr_duplicates[i]/totals[i])*100
		pcr_duplicates_per.append(per)


	if breakdown_per == False:
		plt.ylabel('Count',fontsize=axis_label_size,labelpad=100)
		p1 = plt.bar(ind, cutadapt_removed, bar_width, color='#5e0003',linewidth=0)
		p2 = plt.bar(ind, rrna_removed, bar_width,color='#b21c1c', bottom=cutadapt_removed,linewidth=0)
		p3 = plt.bar(ind, unmapped, bar_width,color='#e5584b', bottom=[i+j for i,j in zip(cutadapt_removed,rrna_removed)],linewidth=0)
		p4 = plt.bar(ind, ambiguous, bar_width,color='#F1B592', bottom=[i+j+z for i,j,z in zip(cutadapt_removed,rrna_removed, unmapped)],linewidth=0)
		p5 = plt.bar(ind, mapped_noncoding, bar_width,color='#C3E188', bottom=[i+j+z+q for i,j,z,q in zip(cutadapt_removed,rrna_removed, unmapped, ambiguous)],linewidth=0)
		p6 = plt.bar(ind, mapped_coding, bar_width,color='#42951F', bottom=[i+j+z+q+b for i,j,z,q,b in zip(cutadapt_removed,rrna_removed, unmapped, ambiguous, mapped_noncoding)],linewidth=0)
		p7 = plt.bar(ind, pcr_duplicates, bar_width,color='#01541e', bottom=[i+j+z+q+b+y for i,j,z,q,b,y in zip(cutadapt_removed,rrna_removed, unmapped, ambiguous, mapped_noncoding,mapped_coding)],linewidth=0)
		p8 = plt.bar(ind, totals, bar_width, color='#5e0003',linewidth=0,alpha=0)
	else:
		plt.ylabel('Percent %',fontsize=axis_label_size,labelpad=100)
		p1 = plt.bar(ind, cutadapt_removed_per, bar_width, color='#5e0003',linewidth=0)
		p2 = plt.bar(ind, rrna_removed_per, bar_width,color='#b21c1c', bottom=cutadapt_removed_per,linewidth=0)
		p3 = plt.bar(ind, unmapped_per, bar_width,color='#e5584b', bottom=[i+j for i,j in zip(cutadapt_removed_per,rrna_removed_per)],linewidth=0)
		p4 = plt.bar(ind, ambiguous_per, bar_width,color='#F1B592', bottom=[i+j+z for i,j,z in zip(cutadapt_removed_per,rrna_removed_per, unmapped_per)],linewidth=0)
		p5 = plt.bar(ind, mapped_noncoding_per, bar_width,color='#C3E188', bottom=[i+j+z+q for i,j,z,q in zip(cutadapt_removed_per,rrna_removed_per, unmapped_per, ambiguous_per)],linewidth=0)
		p6 = plt.bar(ind, mapped_coding_per, bar_width,color='#42951F', bottom=[i+j+z+q+b for i,j,z,q,b in zip(cutadapt_removed_per,rrna_removed_per, unmapped_per, ambiguous_per, mapped_noncoding_per)],linewidth=0)
		p7 = plt.bar(ind, pcr_duplicates_per, bar_width,color='#01541e', bottom=[i+j+z+q+b+y for i,j,z,q,b,y in zip(cutadapt_removed_per,rrna_removed_per, unmapped_per, ambiguous_per, mapped_noncoding_per,mapped_coding_per)],linewidth=0)
		p8 = plt.bar(ind, totals_per, bar_width, color='#5e0003',linewidth=0,alpha=0)
	#Dummy plot point so we can add total reads to the legend
	total_reads_plot = plt.plot(0,0,alpha=0)
	plt.legend(
			(p7[0],p6[0], p5[0], p4[0],p3[0],p2[0],p1[0],total_reads_plot[0]),
			('PCR Duplicates: {:,}'.format(sum(pcr_duplicates)),
				'Mapped coding: {:,}'.format(sum(mapped_coding)),
				'Mapped noncoding: {:,}'.format(sum(mapped_noncoding)),
				'Ambiguous: {:,}'.format(sum(ambiguous)),
				'Unmapped: {:,}'.format(sum(unmapped)),
				'rRNA removed: {:,}'.format(sum(rrna_removed)),
				'Cutadapt removed {:,}'.format(sum(cutadapt_removed)),
				'Total reads: {:,}'.format(all_reads_count)
				),fontsize=legend_size/1.5)


	for i, bar in enumerate(p8.get_children()):
		per = int(round((pcr_duplicates[i]/totals[i])*100,0))
		lab1 = '<div class="tooltip"><span class="tooltiptext">'
		lab1 += '<b>Filename:</b> {}<br><br><b>Pcr duplicates:</b> {:,}  ({}%)'.format(labels[i],pcr_duplicates[i], per)
		per = int(round((mapped_coding[i]/totals[i])*100,0))
		lab1 += '<br><b>Mapped coding:</b> {:,}  ({}%)'.format(mapped_coding[i], per)
		per = int(round((mapped_noncoding[i]/totals[i])*100,0))
		lab1 += '<br><b>Mapped noncoding:</b> {:,}  ({}%)'.format(mapped_noncoding[i], per)
		per = int(round((ambiguous[i]/totals[i])*100,0))
		lab1 += '<br><b>Ambiguous:</b> {:,}  ({}%)'.format(ambiguous[i], per)
		per = int(round((unmapped[i]/totals[i])*100,0))
		lab1 += '<br><b>Unmapped:</b> {:,}  ({}%)'.format(unmapped[i], per)
		per = int(round((rrna_removed[i]/totals[i])*100,0))
		lab1 += '<br><b>rRNA removed:</b> {:,}  ({}%)'.format(rrna_removed[i], per)
		per = int(round((cutadapt_removed[i]/totals[i])*100,0))
		lab1 += '<br><b>Cutadapt_removed:</b> {:,}  ({}%)'.format(cutadapt_removed[i], per)
		#lab1 = (str(lab1).to_html())
		tooltip1 = plugins.LineHTMLTooltip(bar, lab1,voffset=10, hoffset=30,css=line_tooltip_css)
		plugins.connect(fig, tooltip1)

	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph




def single_tran_de(single_tran_de_transcript, sorted_master_list, study_master_list, organism, transcriptome,single_tran_de_study):
	xvals = []
	yvals = []
	labels = []
	range1counts = []
	file_ids = []
	mapped_reads = []
	file_descs = []
	study_names = []
	
	if single_tran_de_study == False:
		for tup in sorted_master_list:
			file_id = tup[0]
			filename = tup[1]
			r1 = float(tup[2])
			r2 = float(tup[3])
			file_desc = tup[4]
			study = tup[5]
			
			ratio = r1/r2
			avg = (r1+r2)/2
			
			xvals.append(avg)
			yvals.append(ratio)
			labels.append(filename)
			range1counts.append(r1)
			file_ids.append(file_id)
			mapped_reads.append(r2)
			file_descs.append(file_desc)
			study_names.append(study)
	else:
		for tup in study_master_list:
			file_id = tup[0]
			study = tup[1]
			r1 = float(tup[2])
			r2 = float(tup[3])

			
			ratio = log(r1/r2,2)
			avg = (r1+r2)/2
			
			xvals.append(avg)
			yvals.append(ratio)
			labels.append(study)
			range1counts.append(r1)
			file_ids.append(file_id)
			mapped_reads.append(r2)
			file_descs.append(study)
			study_names.append(study)
		
		
		
		
	

	full_title = "ORF TPMs ({})"
	x_lab = ''
	y_lab = 'TPM'
	p = figure(plot_width=1300, plot_height=750,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap",logo=None)
	p.title.align="center"
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"

	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Ratio", "@y"),("Max count","@x"),("File name","@labels"),("Range 1 count","@range1counts"),("Mapped reads","@mapped_reads"),("Description","@file_descs"),("Study","@study")]
	source = ColumnDataSource({'x':xvals,'y':yvals,'labels':labels,"range1counts":range1counts,"file_id":file_ids,"mapped_reads":mapped_reads,"file_descs":file_descs,"study":study_names})
	p.scatter('x','y',source=source, alpha=1,color="grey",size=9)
	output_file("scatter10k.html", title="Single transcript differential translation")
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'

	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,%23ff1f00_228,%233BFF00_231,%23ffffff_232,%23000000_&transcript=YLR162W&normalize=F&cov=T&ambig=F&minread=25&maxread=100
	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,228,229,%23ff1f00_230,%233BFF00&transcript=YDR003W&normalize=T&cov=T&ambig=T&minread=18&maxread=45
	#url = "http://143.239.109.139/tripsviz/{}/{}/interactive_plot/?transcript=@labels".format(organism, transcriptome)
	#file_string = ""
	#label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"


	# remove the trailing _ in file_string if it's been populated
	#if file_string:
	#	file_string = file_string[:len(file_string)-1]

	url = "http://trips.ucc.ie/{}/{}/interactive_plot/?files=@file_id&tran={}&ambig=F&minread=25&maxread=150".format(organism, transcriptome,single_tran_de_transcript)

	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)


	#TODO FIX HARDCODED TMP FILE LINK
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)

	graph = file_html(p,CDN)
	return graph








def tran_corr(tran_corr_transcript1, tran_corr_transcript2,sorted_master_list,organism, transcriptome):
	xvals = []
	yvals = []
	labels = []
	range1counts = []
	range2counts = []
	file_ids = []

	for tup in sorted_master_list:
		xvals.append(tup[2])
		yvals.append(tup[3])
		labels.append(tup[1])
		range1counts.append(tup[4])
		range2counts.append(tup[3])
		file_ids.append(tup[0])

	full_title = "Transcript correllation ({})"


	x_lab = '{} count (log2)'.format(tran_corr_transcript1)
	y_lab = '{} count (log2)'.format(tran_corr_transcript2)

	p = figure(plot_width=2000, plot_height=1000,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap",logo=None)
	p.title.align="center"

	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"


	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Ratio", "@y"),("Max count","@x"),("File name","@labels"),("Range 1 count","@range1counts"),("Range 2 count","@range2counts")]
	source = ColumnDataSource({'x':xvals,'y':yvals,'labels':labels,"range1counts":range1counts,"range2counts":range2counts,"file_id":file_ids})
	p.scatter('x','y',source=source, alpha=1,color="grey",size=12)
	output_file("scatter10k.html", title="Single transcript differential translation")
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'

	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,%23ff1f00_228,%233BFF00_231,%23ffffff_232,%23000000_&transcript=YLR162W&normalize=F&cov=T&ambig=F&minread=25&maxread=100
	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,228,229,%23ff1f00_230,%233BFF00&transcript=YDR003W&normalize=T&cov=T&ambig=T&minread=18&maxread=45
	#url = "http://143.239.109.139/tripsviz/{}/{}/interactive_plot/?transcript=@labels".format(organism, transcriptome)
	#file_string = ""
	#label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"


	# remove the trailing _ in file_string if it's been populated
	#if file_string:
	#	file_string = file_string[:len(file_string)-1]

	url = "http://trips.ucc.ie/{}/{}/interactive_plot/?files=@file_id&tran={}&ambig=F&minread=25&maxread=150".format(organism, transcriptome,tran_corr_transcript1)

	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)


	#TODO FIX HARDCODED TMP FILE LINK
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)

	graph = file_html(p,CDN)
	return graph










def explore_offsets(f0_counts, f1_counts, f2_counts, labels,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size):
	width = 0.25
	pos = list(range(len(f0_counts)))
	fig, ax = plt.subplots(figsize=(13,8))
	N = len(f0_counts)
	bar_width = 0.35
	bar_l = [i for i in range(len(f0_counts))]
	tick_pos = [i+(bar_width/2) for i in bar_l]
	ind = np.arange(N)
	for i in range(0,len(f0_counts)):
		if i%3 == 0:
			width = 0.4
		elif i%3 == 1:
			width = 0.0
		else:
			width = -0.4
		p1 = plt.bar(ind[i]+width, f0_counts[i], bar_width, color=greenhex,linewidth=0)
		p2 = plt.bar([pos[i] + width], f1_counts[i], bar_width,color=redhex, bottom=f0_counts[i],linewidth=0)
	plt.ylabel('Count',fontsize=axis_label_size,labelpad=100)
	title_str = "Frame breakdown ({})".format(short_code)
	plt.title(title_str,fontsize=title_size)
	plt.xticks(tick_pos, labels)
	plt.legend((p2[0],p1[0]), ('CDS: out of frame', 'CDS: in frame'))
	ax.set_facecolor(background_col)
	ax.tick_params('both', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")
	totals = []
	for i in range(0,len(f0_counts)):
		curr_total = 0
		curr_total += f0_counts[i]
		curr_total += f1_counts[i]
		curr_total += f2_counts[i]
		if curr_total > 0:
			totals.append(float(curr_total))
		else:
			totals.append(1)
	for i, bar in enumerate(p1.get_children()):
		per = int(round((f0_counts[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Frame 0: {:,}  ({}%)".format(f0_counts[i],per))
		plugins.connect(fig, tooltip1)
	for i, bar in enumerate(p2.get_children()):
		per = int(round((f1_counts[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Frame 1: {:,}  ({}%)".format(f1_counts[i], per))
		plugins.connect(fig, tooltip1)
	for i, bar in enumerate(p3.get_children()):
		per = int(round((f2_counts[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Frame 2: {:,}  ({}%)".format(f2_counts[i], per))
		plugins.connect(fig, tooltip1)

	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadPNG(returnstr=title_str))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph



def replicate_comp(labels, transcript_dict,min_log_val,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,corr_type):
	list_dict = {}
	corr_dict = {}
	for lbl in labels:
		list_dict[lbl] = []
		corr_dict[lbl] = {}
	for transcript in transcript_dict:
		for label in transcript_dict[transcript]:
			list_dict[label].append(transcript_dict[transcript][label])
	for label in list_dict:
		for label2 in list_dict:
			list1 = list_dict[label]
			list2 = list_dict[label2]
			if corr_type == "pearson":
				corr = pearsonr(list1, list2)
			elif corr_type == "spearman":
				corr = spearmanr(list1,list2)
			final_corr = round(corr[0],2)
			corr_dict[label][label2] = final_corr
	tablehtml = "<table class='corr_table'><thead><tr><th></th>"
	sorted_list = sorted(corr_dict.keys())
	for label in sorted_list:
		tablehtml += "<th>{}</th>".format(label)
	tablehtml += "</tr></thead><tbody>"
	for label in sorted_list:
		tablehtml += "<tr><td><b>{}</b></td>".format(label)
		for label2 in sorted_list:
			val = corr_dict[label][label2]
			if val >= 0.98:
				clss = "col1"
			elif val >= 0.96:
				clss = "col2"
			elif val >= 0.94:
				clss = "col3"
			elif val >= 0.92:
				clss = "col4"
			elif val >= 0.90:
				clss = "col5"
			elif val >= 0.88:
				clss = "col6"
			elif val >= 0.86:
				clss = "col7"
			else:
				clss = "col8"
			tablehtml += "<td class='{}'>{}</td>".format(clss, val)
		tablehtml += "</tr>"
	tablehtml += "</tbody></table>"
	return  tablehtml


def most_freq_unmapped(file_paths_dict,short_code):
	master_dict = {}
	for filetype in file_paths_dict:
		for file_id in file_paths_dict[filetype]:
			filepath = file_paths_dict[filetype][file_id]
			if os.path.isfile(filepath):
				sqlite_db = SqliteDict(filepath, autocommit=False)
			else:
				return_str =  "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				return return_str
			if "frequent_unmapped_reads" not in sqlite_db:
				return_str =  "No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				return return_str
			#unmapped reads list is a list of tuples of length 100, first item in tuple is a sequence second is a count
			unmapped_reads_list = sqlite_db["frequent_unmapped_reads"]
			sqlite_db.close()

			for tup in unmapped_reads_list:
				if tup[0] in master_dict:
					master_dict[tup[0]] += tup[1]
				else:
					master_dict[tup[0]] = tup[1]
	title = "Most frequent unmapped reads ({})".format(short_code)
	studyname = ""
	top_reads = (sorted(master_dict.items(), key=operator.itemgetter(1)))[-50:]
	html_table = "<h1><center>{}</center></h1>".format(title)
	html_table += """<table class="unmapped_table">
	<thead><tr><th>Sequence</th><th>Frequency</th><th>Blast Link</th></tr></thead>"""
	for tup in top_reads[::-1]:
		html_table += ("<tr><td>{0}</td><td>{1}</td>    <td><a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY=%3E{2}_unmapped_sequence%0A{0}' target='_blank'>Blast</a></td></tr>".format(tup[0], tup[1],studyname))
	html_table += ("</table>")
	return html_table



def te_table(table_str):
	 return table_str


def heatplot(min_readlen, max_readlen, min_pos, max_pos, positions,readlengths,count_list,heatmap_metagene_type,title,reverse_scale,color_palette,short_code,background_col,maxscaleval,title_size, axis_label_size, subheading_size,marker_size):
	xlabs = []
	ylabs = []
	for i in range(min_readlen,max_readlen+1):
		ylabs.append(i)
	for i in range(min_pos, max_pos):
		xlabs.append(i)
	df = pd.DataFrame(np.column_stack([positions, readlengths, count_list]), columns=['positions',"readlengths","counts"])
	source = ColumnDataSource(df)

	if color_palette == "gist_rainbow":
		color_pallete_list = []
		colormap =cm.get_cmap("gist_rainbow") #choose any matplotlib colormap here
		start_val = 0.75
		for i in range(0,256):
			start_val -= 0.00293
			rgb = colormap(start_val)[:3]
			hexval = matplotlib.colors.rgb2hex(rgb)
			color_pallete_list.append(hexval)
	else:
		color_pallete_list = all_palettes[color_palette][max(all_palettes[color_palette].keys())]
	
	colormap =cm.get_cmap("gist_rainbow") #choose any matplotlib colormap here
	rgb = colormap(0.5)[:3]
	hexval = matplotlib.colors.rgb2hex(rgb)

	if reverse_scale == True:
		color_pallete_list = color_pallete_list[::-1]
	
	print (count_list)
	print (max(filter(lambda v: v is not None, count_list)))
	
	if maxscaleval == "None":
		color_mapper = LogColorMapper(palette=color_pallete_list, low=1, high=max(filter(lambda v: v is not None, count_list)),nan_color=background_col)
	else:
		color_mapper = LogColorMapper(palette=color_pallete_list, low=1, high=maxscaleval,nan_color=background_col)
	p = figure(plot_width=1300, plot_height=750,x_range=(min_pos,max_pos), y_range=(min_readlen,max_readlen),title="Heatmap ({})".format(short_code),
			   y_axis_label="Readlengths", toolbar_location="below",tools = "reset,pan,box_zoom")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_col
	if heatmap_metagene_type == "metagene_start":
		p.xaxis.axis_label = ('Position relative to start')
	elif heatmap_metagene_type == "metagene_stop":
		p.xaxis.axis_label = ('Position relative to stop')
	p.rect(x="positions", y="readlengths", width=1, height=1,source=source,fill_color={'field': 'counts',"transform":color_mapper},line_color=None)
	color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),title="Counts", border_line_color=None, location=(0,0))
	p.add_layout(color_bar, 'right')
	output_file("scatter10k.html", title="Heatmap ({})".format(short_code))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += file_html(p,CDN)
	return graph



def ticker(tick):
	return "{:.0f} + {:.2f}".format(tick, tick % 1)




def rust_dwell(codon_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size):
	aa_color_dict = {"gly":"white",
			"arg":"blue",
			"ser":"green",
			"trp":"lightgreen",
			"cys":"white",
			"glu":"red",
			"asp":"red",
			"lys":"blue",
			"asn":"green",
			"gln":"green",
			"his":"blue",
			"tyr":"lightgreen",
			"ala":"grey",
			"thr":"green",
			"pro":"white",
			"val":"grey",
			"met":"grey",
			"ile":"grey",
			"leu":"grey",
			"phe":"lightgreen"}
	
	min_count = float(min(codon_count_dict.values()))+0.000001
	max_count = float(max(codon_count_dict.values()))
	max_rust_ratio = (max_count/(min_count))
	aa_dict = {"gly":["GGT","GGC","GGA","GGG"],
			   "arg":["AGA","AGG","CGT","CGC","CGA","CGG"],
			   "ser":["AGT","AGC","TCT","TCC","TCA","TCG"],
			   "trp":["TGG"],
			   "cys":["TGT","TGC"],
			   "glu":["GAA","GAG"],
			   "asp":["GAT","GAC"],
			   "lys":["AAA","AAG"],
			   "asn":["AAT","AAC"],
			   "gln":["CAA","CAG"],
			   "his":["CAT","CAC"],
			   "tyr":["TAT","TAC"],
			   "ala":["GCT","GCC","GCA","GCG"],
			   "thr":["ACT","ACC","ACA","ACG"],
			   "pro":["CCT","CCC","CCA","CCG"],
			   "val":["GTT","GTC","GTA","GTG"],
			   "met":["ATG"],
			   "ile":["ATA","ATT","ATC"],
			   "leu":["CTT","CTC","CTA","CTG","TTA","TTG"],
			   "phe":["TTT","TTC"]}
	avg_dict = {}
	for aa in aa_dict:
		tot = 0.0
		for codon in aa_dict[aa]:
			tot += codon_count_dict[codon]
		avg = tot/len(aa_dict[aa])
		avg_dict[aa] = avg
	sorted_list = sorted(avg_dict.items(), key=lambda x: x[1])
	p = figure(plot_width=1300, plot_height=750, y_axis_label='RUST ratio (relative to min)',title="RUST: Relative codon dwell times ({})".format(short_code),toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover",logo=None)#,y_range=Range1d(bounds=(0, max_rust_ratio)))
	p.y_range = Range1d(0,max_rust_ratio)
	p.x_range = Range1d(0,23)
	p.xaxis.visible = False
	p.title.align="center"
	p.title.text_font_size = "24pt"
	p.xaxis.axis_label_text_font_size = "22pt"
	p.xaxis.major_label_text_font_size = "17pt"
	p.yaxis.axis_label_text_font_size = "22pt"
	p.yaxis.major_label_text_font_size = "17pt"
	p.background_fill_color = background_col
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"
	x_vals = []
	y_vals = []
	codons = []
	x_val = 1
	for tup in sorted_list:
		aa = tup[0]
		for codon in aa_dict[aa]:
			codons.append(codon)
			x_vals.append(x_val)
			rust_ratio = (float(codon_count_dict[codon]+0.000001)/min_count)
			y_vals.append(rust_ratio)
			mytext = Label(x=x_val-0.2, y=-10, text=aa,background_fill_color=aa_color_dict[aa],text_font_size="18pt")
			p.add_layout(mytext)
		x_val += 1
	p.legend.label_text_size = "30pt"
	p.scatter(-1,-1,color="red",legend="Acidic")
	p.scatter(-1,-1,color="lightgreen",legend="Aromatic")
	p.scatter(-1,-1,color="green",legend="Polar uncharged")
	p.scatter(-1,-1,color="grey",legend="Aliphatic")
	p.scatter(-1,-1,color="blue",legend="Basic")
	source = ColumnDataSource({'x': x_vals,'y':y_vals,'codons':codons})
	p.scatter('x','y',source=source, color="grey",size=16)
	p.y_range=Range1d(-10, max(y_vals)*1.1)
	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Rust ratio", "@y"),
					  ("Codon","@codons")]
	output_file("scatter10k.html", title="RUST: Relative codon dwell times ({})".format(short_code))
	graph = "<style>.mpld3-xaxis {{font-size: {0}px;}} .mpld3-yaxis {{font-size: {0}px;}}</style>".format(marker_size)
	graph += "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += file_html(p,CDN)
	return graph




def codon_usage(codon_dict,short_code,title_size, axis_label_size, marker_size):
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
	
	codon_list = ["ATG","TTT","TTC","CTT", "CTC", "CTA", "CTG","TTA","TTG", "AGT", "AGC","TCT", "TCC", "TCA", "TCG","TAT", "TAC","TGT", "TGC","TGG","CCT", "CCC", "CCA", "CCG",
			   "CAT", "CAC", "CAA", "CAG","AGA", "AGG","CGT", "CGC", "CGA", "CGG","ATT", "ATC", "ATA","ACT", "ACC", "ACA", "ACG","AAT", "AAC", "AAA", "AAG", 
				"GTT", "GTC", "GTA", "GTG","GCT", "GCC", "GCA", "GCG","GAT", "GAC", "GAA", "GAG","GGT", "GGC", "GGA", "GGG","TAG","TAA","TGA"]

	curr_count = 0
	for codon in codon_list:
		curr_count+=1
		allxvals.append(curr_count)
		allyvals.append(codon_dict[codon]["ribo_count"]/codon_dict[codon]["codon_count"])
		alllabels.append(codon)
		amino_acids.append(aa_dict[codon])
	full_title = "Codon usage ({})".format(short_code)
	x_lab = ''
	y_lab = 'Normalised Count'
	min_y = min(0,min(allyvals))-.02
	max_y = max(allyvals)+.02
	p = figure(plot_width=1300, plot_height=750,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap",logo=None,y_range=(min_y,max_y))
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

	graph = file_html(p,CDN)
	return graph



def calc_mrnadist_factor(mrna_dist_dict):
	maxval = 0
	for key in mrna_dist_dict:
		for region in mrna_dist_dict[key]:
			count = mrna_dist_dict[key][region]
			if count > maxval:
				maxval = int(count)
	string_maxval = str(maxval)
	zeroes = len(string_maxval)-1
	zeroes_string = "0"*zeroes
	factor = int("1"+zeroes_string)
	for key in mrna_dist_dict:
		for region in mrna_dist_dict[key]:
			mrna_dist_dict[key][region] = float(mrna_dist_dict[key][region])/factor
	return mrna_dist_dict, zeroes





def mrna_dist(mrna_dist_dict, short_code, background_col, title_size, axis_label_size, subheading_size, marker_size, mrna_dist_per,md_start,md_stop,legend_size):
	if mrna_dist_per == False:
		mrna_dist_dict,factor = calc_mrnadist_factor(mrna_dist_dict)
	fig, ax = plt.subplots(figsize=(13,8))
	returnstr = "Sample, 5' leader, Start codon, CDS, Stop codon, 3' trailer\n"
	
	#Add two because we plot two empty bars at the beginning and end for aesthethics
	N = len(mrna_dist_dict.keys())+2
	bar_width = 0.35
	bar_l = [i for i in range(N)]
	tick_pos = [i+(bar_width/2) for i in bar_l]
	ind = np.arange(N)    # the x locations for the groups

	labels = [""]
	five_leaders = [0]
	start_codons = [0]
	cds = [0]
	stop_codons = [0]
	three_trailers = [0]

	five_leaders_per = [0]
	start_codons_per = [0]
	cds_per = [0]
	stop_codons_per = [0]
	three_trailers_per = [0]



	for key in mrna_dist_dict:
		#try:
		#	total = float(mrna_dist_dict[key]["total"])
		#except:
		total = float(mrna_dist_dict[key]["5_leader"]+mrna_dist_dict[key]["start_codon"]+mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["stop_codon"]+mrna_dist_dict[key]["3_trailer"])
		labels.append(key)
		five_leaders.append(mrna_dist_dict[key]["5_leader"])
		five_leaders_per.append((mrna_dist_dict[key]["5_leader"]/total)*100)
		three_trailers.append(mrna_dist_dict[key]["3_trailer"])
		three_trailers_per.append((mrna_dist_dict[key]["3_trailer"]/total)*100)
		if md_start == True and md_stop == True:
			start_codons.append(mrna_dist_dict[key]["start_codon"])
			start_codons_per.append((mrna_dist_dict[key]["start_codon"]/total)*100)
			cds.append(mrna_dist_dict[key]["cds"])
			cds_per.append((mrna_dist_dict[key]["cds"]/total)*100)
			stop_codons.append(mrna_dist_dict[key]["stop_codon"])
			stop_codons_per.append((mrna_dist_dict[key]["stop_codon"]/total)*100)
		elif md_start == True and md_stop == False:
			start_codons.append(mrna_dist_dict[key]["start_codon"])
			start_codons_per.append((mrna_dist_dict[key]["start_codon"]/total)*100)
			cds.append(mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["stop_codon"])
			cds_per.append(((mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["stop_codon"])/total)*100)
			stop_codons.append(0)
			stop_codons_per.append(0)
		elif md_start == False and md_stop == True:
			stop_codons.append(mrna_dist_dict[key]["stop_codon"])
			stop_codons_per.append((mrna_dist_dict[key]["stop_codon"]/total)*100)
			cds.append(mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["start_codon"])
			cds_per.append(((mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["start_codon"])/total)*100)
			start_codons.append(0)
			start_codons_per.append(0)
		elif md_start == False and md_stop == False:
			cds.append(mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["start_codon"]+mrna_dist_dict[key]["stop_codon"])
			cds_per.append(((mrna_dist_dict[key]["cds"]+mrna_dist_dict[key]["start_codon"]+mrna_dist_dict[key]["stop_codon"])/total)*100)	
			stop_codons.append(0)
			stop_codons_per.append(0)
			start_codons.append(0)
			start_codons_per.append(0)			
		returnstr += "{},{},{},{},{},{}\n".format(key, mrna_dist_dict[key]["5_leader"],mrna_dist_dict[key]["start_codon"],mrna_dist_dict[key]["cds"],mrna_dist_dict[key]["stop_codon"],mrna_dist_dict[key]["3_trailer"])
	labels.append("")
	five_leaders.append(0)
	start_codons.append(0)
	cds.append(0)
	stop_codons.append(0)
	three_trailers.append(0)
	five_leaders_per.append(0)
	start_codons_per.append(0)
	cds_per.append(0)
	stop_codons_per.append(0)
	three_trailers_per.append(0)
	if mrna_dist_per == False:
		p1 = plt.bar(ind, five_leaders, bar_width, color='#ff6d6d',linewidth=0)
		p2 = plt.bar(ind, start_codons, bar_width,color='#a38b22', bottom=five_leaders,linewidth=0)
		p3 = plt.bar(ind, cds, bar_width,color='#00d8cd', bottom=[i+j for i,j in zip(five_leaders,start_codons)],linewidth=0)
		p4 = plt.bar(ind, stop_codons, bar_width,color='#05bc27', bottom=[i+j+z for i,j,z in zip(five_leaders,start_codons, cds)],linewidth=0)
		p5 = plt.bar(ind, three_trailers, bar_width,color='#ff77d4', bottom=[i+j+z+q for i,j,z,q in zip(five_leaders,start_codons, cds, stop_codons)],linewidth=0)
	else:
		p1 = plt.bar(ind, five_leaders_per, bar_width, color='#ff6d6d',linewidth=0)
		p2 = plt.bar(ind, start_codons_per, bar_width,color='#a38b22', bottom=five_leaders_per,linewidth=0)
		p3 = plt.bar(ind, cds_per, bar_width,color='#00d8cd', bottom=[i+j for i,j in zip(five_leaders_per,start_codons_per)],linewidth=0)
		p4 = plt.bar(ind, stop_codons_per, bar_width,color='#05bc27', bottom=[i+j+z for i,j,z in zip(five_leaders_per,start_codons_per, cds_per)],linewidth=0)
		p5 = plt.bar(ind, three_trailers_per, bar_width,color='#ff77d4', bottom=[i+j+z+q for i,j,z,q in zip(five_leaders_per,start_codons_per, cds_per, stop_codons_per)],linewidth=0)	
	if mrna_dist_per == False:
		plt.ylabel('Count (x 10 {})'.format(factor),fontsize=axis_label_size,labelpad=100)
	else:
		plt.ylabel('Percent %',fontsize=axis_label_size,labelpad=100)
	title_str = "Reads breakdown ({})".format(short_code)
	plt.title(title_str,fontsize=title_size)
	if len(labels) <= 5:
		xlabel_size = marker_size
	else:
		xlabel_size = marker_size/(len(labels)/5.0)

	plt.xticks(tick_pos, labels,fontsize=xlabel_size)
	if md_start == True and md_stop == True:
		plt.legend((p5[0], p4[0],p3[0],p2[0],p1[0]), ('Three trailers','Stop codons','Cds','Start codons','Five leaders'),fontsize=legend_size)
	elif md_start == True and md_stop == False:
		plt.legend((p5[0],p3[0],p2[0],p1[0]), ('Three trailers','Cds','Start codons','Five leaders'),fontsize=legend_size)
	elif md_start == False and md_stop == True:
		plt.legend((p5[0], p4[0],p3[0],p1[0]), ('Three trailers','Stop codons','Cds','Five leaders'),fontsize=legend_size)
	elif md_start == False and md_stop == False:
		plt.legend((p5[0],p3[0],p1[0]), ('Three trailers','Cds','Five leaders'),fontsize=legend_size)

	ax.set_facecolor(background_col)
	ax.tick_params('y', labelsize=marker_size)
	plt.grid(color="white", linewidth=2,linestyle="solid")

	totals = []
	for i in range(0,len(cds)):
		curr_total = 0
		curr_total += cds[i]
		curr_total += three_trailers[i]
		curr_total += stop_codons[i]
		curr_total += five_leaders[i]
		curr_total += start_codons[i]
		if curr_total > 0:
			totals.append(float(curr_total))
		else:
			totals.append(1)


	for i, bar in enumerate(p1.get_children()):
		per = int(round((five_leaders[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="5' leaders: {:,}  ({}%)".format(five_leaders[i],per))
		#tooltip1 = PointHTMLTooltip(p1[i], p1labels[i],voffset=10, hoffset=10, css=point_tooltip_css)
		plugins.connect(fig, tooltip1)

	for i, bar in enumerate(p2.get_children()):
		per = int(round((start_codons[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Start codons: {:,}  ({}%)".format(start_codons[i], per))
		#tooltip1 = PointHTMLTooltip(p1[i], p1labels[i],voffset=10, hoffset=10, css=point_tooltip_css)
		plugins.connect(fig, tooltip1)

	for i, bar in enumerate(p3.get_children()):
		per = int(round((cds[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Cds: {:,}  ({}%)".format(cds[i], per))
		#tooltip1 = PointHTMLTooltip(p1[i], p1labels[i],voffset=10, hoffset=10, css=point_tooltip_css)
		plugins.connect(fig, tooltip1)

	for i, bar in enumerate(p4.get_children()):
		per = int(round((stop_codons[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="Stop codons: {:,}  ({}%)".format(stop_codons[i], per))
		#tooltip1 = PointHTMLTooltip(p1[i], p1labels[i],voffset=10, hoffset=10, css=point_tooltip_css)
		plugins.connect(fig, tooltip1)

	for i, bar in enumerate(p5.get_children()):
		per = int(round((three_trailers[i]/totals[i])*100,0))
		tooltip1 = plugins.LineLabelTooltip(bar, label="3' trailers: {:,}  ({}%)".format(three_trailers[i], per))
		#tooltip1 = PointHTMLTooltip(p1[i], p1labels[i],voffset=10, hoffset=10, css=point_tooltip_css)
		plugins.connect(fig, tooltip1)

	plugins.connect(fig, TopToolbar(yoffset=750,xoffset=600),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	graph = "<div style='padding-left: 55px;padding-top: 22px;'> <a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a> </div>".format(short_code)
	graph += mpld3.fig_to_html(fig)
	return graph
