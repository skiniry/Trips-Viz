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
import collections
import pandas as pd
from math import log
from bokeh.plotting import figure, show, output_file
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette
from bokeh.layouts import column
from bokeh.models.widgets import TextInput
from scipy.stats.stats import spearmanr,pearsonr
from bokeh.io import show
from bokeh.models import (
	TapTool,
	CustomJS,
	OpenURL,
	Label,
	LogTicker,
	ColumnDataSource,
	HoverTool,
	LinearColorMapper,
	LogColorMapper,
	BasicTicker,
	PrintfTickFormatter,
	ColorBar
)

import logging

#import seaborn


#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
#seaborn.set_style("whitegrid", {'axes.grid' : False})




#Defaults
#lite = "y"
#min_read = 25
#max_read = 35
#ribo = "a"
#rna = "a"
#subcodon = "d"
#tran = ""
#ambig = "u"
#user_ribo_files  = []
#user_rna_files  = []
#offset_dict = {}



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

def generate_plot(sorted_min_exp_list,bin_list,organism, label,transcriptome,riboseq1,riboseq2,rnaseq1,rnaseq2,background_color, short_code,normalized,filename,no_groups,title_size, axis_label_size, subheading_size,marker_size,ambiguous,gene_list):
	#Convert gene_list from string to list
	logging.debug("generate plot called")
	#logging.debug("sorted_min_exp_list: {}".format(sorted_min_exp_list))
	logging.debug("bin_list: {}".format(bin_list))
	if gene_list != "":
		gene_list = gene_list.replace(","," ").replace("\t"," ")
		split_list = gene_list.split(" ")
		gene_list = []
		for item in split_list:
			gene_list.append(item.strip(" ").upper())
	posallxvals = []
	posallyvals = []
	posalllabels = []
	posallgenes = []
	poszscores = []

	negallxvals = []
	negallyvals = []
	negalllabels = []
	negallgenes = []
	negzscores = []

	nonde_xvals = []
	nonde_yvals = []
	nonde_labels = []
	nonde_allgenes = []
	nonde_zscores = []
	


	upper_thresholds_y = [bin_list[0][2]]
	upper_thresholds_x = [(sorted_min_exp_list[0][1])*0.99]
	lower_thresholds_y = [bin_list[0][3]]
	lower_thresholds_x = [(sorted_min_exp_list[0][1])*0.99]

	# Easiest way to plot threshold lines is to just plot the threshold and min_exp for each bin, but this will result in slanted lines, i.ie if my first bin threshold is 8 and my second is 4
	# there would be a slanted line between 8 and 4, but this makes no sense as the threshold is the same across an entire bin, so instead we plot two y values at every step (except first and last)

	#d((transcript, log(min_exp,2), fold_change,gene))
	cur_count = 0
	bin_count = 0
	for i in range(0,len(sorted_min_exp_list)):
		cur_count+=1
		if cur_count == 300:
			#To x we add the log2(min exp) of the 300th (or multiple of) item in min exp list
			bin_count += 1
			try:
				tst = bin_list[bin_count]
			except:
				bin_count -= 1
			cur_count = 0
			upper_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			upper_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			upper_thresholds_y.append(bin_list[bin_count-1][2])
			upper_thresholds_y.append(bin_list[bin_count][2])
			lower_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			lower_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			lower_thresholds_y.append(bin_list[bin_count-1][3])
			lower_thresholds_y.append(bin_list[bin_count][3])
		if bin_list[bin_count][1] == 0.0:
			continue
		z_score = (sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1])
		
		if gene_list == "":
			if sorted_min_exp_list[i][2] <= bin_list[bin_count][2] and sorted_min_exp_list[i][2] >= bin_list[bin_count][3]:
				nonde_xvals.append(sorted_min_exp_list[i][1])
				nonde_yvals.append(sorted_min_exp_list[i][2])
				nonde_labels.append(sorted_min_exp_list[i][0])
				nonde_allgenes.append(sorted_min_exp_list[i][3])
				nonde_zscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
			elif sorted_min_exp_list[i][2] >  bin_list[bin_count][2]:
				posallxvals.append(sorted_min_exp_list[i][1])
				posallyvals.append(sorted_min_exp_list[i][2])
				posalllabels.append(sorted_min_exp_list[i][0])
				posallgenes.append(sorted_min_exp_list[i][3])
				poszscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
			elif sorted_min_exp_list[i][2] <  (bin_list[bin_count][3]):
				negallxvals.append(sorted_min_exp_list[i][1])
				negallyvals.append(sorted_min_exp_list[i][2])
				negalllabels.append(sorted_min_exp_list[i][0])
				negallgenes.append(sorted_min_exp_list[i][3])
				negzscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
		else:
			#If the user passes a gene list only highlight those genes red/green depending on if fold changes is >/< 0
			if sorted_min_exp_list[i][3] not in gene_list:
				nonde_xvals.append(sorted_min_exp_list[i][1])
				nonde_yvals.append(sorted_min_exp_list[i][2])
				nonde_labels.append(sorted_min_exp_list[i][0])
				nonde_allgenes.append(sorted_min_exp_list[i][3])
				nonde_zscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
			elif sorted_min_exp_list[i][3] in gene_list and sorted_min_exp_list[i][2] > 0:
				posallxvals.append(sorted_min_exp_list[i][1])
				posallyvals.append(sorted_min_exp_list[i][2])
				posalllabels.append(sorted_min_exp_list[i][0])
				posallgenes.append(sorted_min_exp_list[i][3])
				poszscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
			elif sorted_min_exp_list[i][3] in gene_list and sorted_min_exp_list[i][2] < 0:
				negallxvals.append(sorted_min_exp_list[i][1])
				negallyvals.append(sorted_min_exp_list[i][2])
				negalllabels.append(sorted_min_exp_list[i][0])
				negallgenes.append(sorted_min_exp_list[i][3])
				negzscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))


	# Add z score thresholds for the last bin, which wouldn't have been covered in the above loop
	upper_thresholds_x.append((sorted_min_exp_list[-1][1])*0.99)
	upper_thresholds_y.append(bin_list[bin_count][2])
	lower_thresholds_x.append((sorted_min_exp_list[-1][1])*0.99)
	lower_thresholds_y.append(bin_list[bin_count][3])

	full_title = "Differential translation ({})".format(short_code)

	if no_groups == 1:
		x_lab = 'Geometric mean (log2)'
		y_lab = 'Fold change (log2)'
	else:
		x_lab = 'Average Geometric mean (log2)'
		y_lab = 'Average Fold change (log2)'
	p = figure(plot_width=1300, plot_height=750,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_color
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"
	if gene_list == "":
		p.line(upper_thresholds_x, upper_thresholds_y, color="yellow",line_width=3)
		p.line(lower_thresholds_x, lower_thresholds_y, color="yellow",line_width=3)
	source = ColumnDataSource({'x': nonde_xvals,'y':nonde_yvals,'labels':nonde_labels, 'genes':nonde_allgenes})
	sct = p.scatter('x','y',source=source, alpha=1,color="grey",size=5)
	#source = ColumnDataSource({'x': posallxvals,'y':posallyvals,'labels':posalllabels, 'genes':posallgenes,"zscores":poszscores})
	source = ColumnDataSource({'x': posallxvals,'y':posallyvals,'labels':posalllabels, 'genes':posallgenes})
	p.scatter('x','y',source=source, alpha=1,color="green",size=10)
	hover = p.select(dict(type=HoverTool))
	#hover.tooltips = [("Fold change", "@y"),("Geometric mean","@x"),("Transcript","@labels"),("Genes","@genes"), ("Z score","@zscores")]
	hover.tooltips = [("Fold change", "@y"),("Geometric mean","@x"),("Transcript","@labels"),("Genes","@genes")]
	#source = ColumnDataSource({'x': negallxvals,'y':negallyvals,'labels':negalllabels, 'genes':negallgenes,"zscores":negzscores})
	source = ColumnDataSource({'x': negallxvals,'y':negallyvals,'labels':negalllabels, 'genes':negallgenes})
	p.scatter('x','y',source=source, alpha=1,color="red",size=10)
	#source = ColumnDataSource({'x': nonde_xvals,'y':nonde_yvals,'labels':nonde_labels, 'genes':nonde_allgenes,"zscores":nonde_zscores})

	output_file("scatter10k.html", title="Differential translation")
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'

	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,%23ff1f00_228,%233BFF00_231,%23ffffff_232,%23000000_&transcript=YLR162W&normalize=F&cov=T&ambig=F&minread=25&maxread=100
	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,228,229,%23ff1f00_230,%233BFF00&transcript=YDR003W&normalize=T&cov=T&ambig=T&minread=18&maxread=45
	#url = "http://143.239.109.139/tripsviz/{}/{}/interactive_plot/?transcript=@labels".format(organism, transcriptome)
	file_string = ""
	label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"
	if riboseq1:
		if riboseq1[0] != "":
			for file_id in riboseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%23007a02_")

	if riboseq2:
		if riboseq2[0] != "":
			for file_id in riboseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23960000_")
	if rnaseq1:
		if rnaseq1[0] != "":
			for file_id in rnaseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%2374ed76_")
	if rnaseq2:
		if rnaseq2[0] != "":
			for file_id in rnaseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23ff6d6d_")

	# remove the trailing _ in file_string if it's been populated
	if file_string:
		file_string = file_string[:len(file_string)-1]

	if ambiguous == True:
		ambig = "T"
	else:
		ambig = "F"

	url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@labels&normalize={}&cov=T&ambig={}&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0],ambig)


	'''
	hili_gene = CustomJS(args=dict(sct=sct,source=source), code="""
			console.log("Called hili gene");
			var f = cb_obj['value']
			var x = source["genes"][f]
			console.log("x is" + x);
			console.log("F is "+f);
			sct.glyph.fill_color = '#4286f4';
		""")
	'''

	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)
	#text_input = TextInput(value="B2M",title="Gene",callback=hili_gene)


	#TODO FIX HARDCODED TMP FILE LINK
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	logging.debug("Returning plot")
	return graph




def ribo_vs_rna(ribo_rna_dict,organism,transcriptome,riboseq1,riboseq2,rnaseq1,rnaseq2,background_col,short_code,normalized,filename,no_groups,title_size, axis_label_size, subheading_size,marker_size,ambiguous,gene_list,label):
	#Convert gene_list from string to list
	if gene_list != "":
		gene_list = gene_list.replace(","," ").replace("\t"," ")
		split_list = gene_list.split(" ")
		gene_list = []
		for item in split_list:
			gene_list.append(item.strip(" ").upper())
	x_values = []
	y_values = []
	genes = []
	trans = []
	
	hili_x_values = []
	hili_y_values = []
	hili_genes = []
	hili_trans = []
	for gene in ribo_rna_dict:
		if label == "TE":
			if gene not in gene_list:
				y_values.append(log(ribo_rna_dict[gene]["ribo2"]/ribo_rna_dict[gene]["ribo1"],2))
				x_values.append(log(ribo_rna_dict[gene]["rna2"]/ribo_rna_dict[gene]["rna1"],2))
				genes.append(gene)
				trans.append(ribo_rna_dict[gene]["tran"])
			else:
				hili_y_values.append(log(ribo_rna_dict[gene]["ribo2"]/ribo_rna_dict[gene]["ribo1"],2))
				hili_x_values.append(log(ribo_rna_dict[gene]["rna2"]/ribo_rna_dict[gene]["rna1"],2))
				hili_genes.append(gene)
				hili_trans.append(ribo_rna_dict[gene]["tran"])
		elif label == "Riboseq":
			if gene not in gene_list:
				if ribo_rna_dict[gene]["ribo1"] >= 1:
					y_values.append(log(ribo_rna_dict[gene]["ribo1"],2))
				else:
					y_values.append(0)
				if ribo_rna_dict[gene]["ribo2"] >=  1:
					x_values.append(log(ribo_rna_dict[gene]["ribo2"],2))
				else:
					x_values.append(0)
				genes.append(gene)
				trans.append(ribo_rna_dict[gene]["tran"])
			else:
				if ribo_rna_dict[gene]["ribo1"] >= 1:
					hili_y_values.append(log(ribo_rna_dict[gene]["ribo1"],2))
				else:
					hili_y_values.append(0)
				if ribo_rna_dict[gene]["ribo2"] >= 1:
					hili_x_values.append(log(ribo_rna_dict[gene]["ribo2"],2))
				else:
					hili_x_values.append(0)
				hili_genes.append(gene)
				hili_trans.append(ribo_rna_dict[gene]["tran"])
		elif label == "Rnaseq":
			if gene not in gene_list:
				if ribo_rna_dict[gene]["rna1"] != 0:
					y_values.append(log(ribo_rna_dict[gene]["rna1"],2))
				else:
					y_values.append(0)
				if ribo_rna_dict[gene]["rna2"] != 0:
					x_values.append(log(ribo_rna_dict[gene]["rna2"],2))
				else:
					x_values.append(0)
				genes.append(gene)
				trans.append(ribo_rna_dict[gene]["tran"])
			else:
				if ribo_rna_dict[gene]["rna1"] != 0:
					hili_y_values.append(log(ribo_rna_dict[gene]["rna1"],2))
				else:
					hili_y_values.append(0)
				if ribo_rna_dict[gene]["rna2"] != 0:
					hili_x_values.append(log(ribo_rna_dict[gene]["rna2"],2))
				else:
					hili_x_values.append(0)
				hili_genes.append(gene)
				hili_trans.append(ribo_rna_dict[gene]["tran"])
	source = ColumnDataSource({'x': x_values,'y':y_values,'trans':trans, 'genes':genes})
	if label == "TE":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq FC (log2)",  y_axis_label='Ribo-Seq FC (log2)',title="Ribo-Seq FC vs RNA-Seq FC ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Riboseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="Ribo-Seq Cond2 count (log2)",  y_axis_label='Ribo-Seq Cond1 count (log2)',title="Ribo-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Rnaseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq Cond2 count (log2)",  y_axis_label='RNA-Seq Cond1 count (log2)',title="RNA-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"


	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='grey')
	source = ColumnDataSource({'x':hili_x_values, 'y':hili_y_values, 'trans':hili_trans,'genes':hili_genes})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#4286f4')
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	
	if label == "TE":
		p.line([-8,8], [-8,8], color="#cccccc",line_width=1)
		hover.tooltips = [("Ribo fc", "@y"),("RNA fc","@x"),("Genes","@genes"),("Transcript","@trans")]
	elif label == "Riboseq":
			p.line([0,16], [0,16], color="#cccccc",line_width=1)
			hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans")]
			#corr = spearmanr(x_values, y_values)
			pearson_corr = pearsonr(x_values, y_values)
			mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
			p.add_layout(mytext)
	else:
			p.line([0,16], [0,16], color="#cccccc",line_width=1)
			hover.tooltips = [("Rna-seq Cond 1 count (log2)", "@y"),("Rna-seq Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans")]
			hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans")]
			#corr = spearmanr(x_values, y_values)
			pearson_corr = pearsonr(x_values, y_values)
			mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
			p.add_layout(mytext)
	file_string = ""
	label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"
	
	if riboseq1:
		if riboseq1[0] != "":
			for file_id in riboseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%23007a02_")

	if riboseq2:
		if riboseq2[0] != "":
			for file_id in riboseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23960000_")
	if rnaseq1:
		if rnaseq1[0] != "":
			for file_id in rnaseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%2374ed76_")
	if rnaseq2:
		if rnaseq2[0] != "":
			for file_id in rnaseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23ff6d6d_")
	# remove the trailing _ in file_string if it's been populated
	if file_string:
		file_string = file_string[:len(file_string)-1]
	if ambiguous == True:
		ambig = "T"
	else:
		ambig = "F"
	url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@trans&normalize={}&cov=T&ambig={}&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0],ambig)
	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph






def deseq2_plot(ribo_rna_dict,organism,transcriptome,riboseq1,riboseq2,rnaseq1,rnaseq2,background_col,short_code,normalized,filename,no_groups,title_size, axis_label_size, subheading_size,marker_size,ambiguous,gene_list,label,minzscore):
	#Convert gene_list from string to list
	if gene_list != "":
		gene_list = gene_list.replace(","," ").replace("\t"," ")
		split_list = gene_list.split(" ")
		gene_list = []
		for item in split_list:
			gene_list.append(item.strip(" ").upper())
	x_values = []
	y_values = []
	basemeans = []
	lfcses = []
	genes = []
	trans = []
	
	highlight_x_values = []
	highlight_y_values = []
	highlight_basemeans = []
	highlight_lfcses = []
	highlight_genes = []
	highlight_trans = []
	
	teup_x_values = []
	teup_y_values = []
	teup_basemeans = []
	teup_lfcses = []
	teup_genes = []
	teup_trans = []
	tedown_x_values = []
	tedown_y_values = []
	tedown_basemeans = []
	tedown_lfcses = []
	tedown_genes = []
	tedown_trans = []	
	
	rnaup_x_values = []
	rnaup_y_values = []
	rnaup_basemeans = []
	rnaup_lfcses = []
	rnaup_genes = []
	rnaup_trans = []
	rnadown_x_values = []
	rnadown_y_values = []
	rnadown_basemeans = []
	rnadown_lfcses = []
	rnadown_genes = []
	rnadown_trans = []	

	largest_fc = 0

	for gene in ribo_rna_dict:
		if label == "TE":
			x = ribo_rna_dict[gene]["rna_fc"]
			y = ribo_rna_dict[gene]["ribo_fc"]
			basemean = ribo_rna_dict[gene]["te_basemean"]
			lfcSE = ribo_rna_dict[gene]["te_lfcSE"]
			te_padj = ribo_rna_dict[gene]["te_padj"]
			ribo_padj = ribo_rna_dict[gene]["ribo_padj"]
			rna_padj = ribo_rna_dict[gene]["rna_padj"]
			
			if x != "NA" and y != "NA":
				x = float(x)
				y = float(y)
				if abs(x) > largest_fc:
					largest_fc = abs(x)
				if abs(y) > largest_fc:
					largest_fc = abs(y)
				if "NA" in te_padj:
					te_padj = 1
				if "NA" in ribo_padj:
					ribo_padj = 1
				if "NA" in rna_padj:
					rna_padj = 1
				te_padj = float(te_padj)
				ribo_padj = float(ribo_padj)
				rna_padj = float(rna_padj)
				if gene not in gene_list:
					if te_padj > (minzscore/100) and ribo_padj > (minzscore/100) and rna_padj > (minzscore/100):
						y_values.append(y)
						x_values.append(x)
						basemeans.append(basemean)
						lfcses.append(lfcSE)
						genes.append(gene)
						trans.append(ribo_rna_dict[gene]["tran"])
					elif te_padj < (minzscore/100) and x < y :
						teup_y_values.append(y)
						teup_x_values.append(x)
						teup_basemeans.append(basemean)
						teup_lfcses.append(lfcSE)
						teup_genes.append(gene)
						teup_trans.append(ribo_rna_dict[gene]["tran"])
					elif te_padj < (minzscore/100) and y < x :
						tedown_y_values.append(y)
						tedown_x_values.append(x)
						tedown_basemeans.append(basemean)
						tedown_lfcses.append(lfcSE)
						tedown_genes.append(gene)
						tedown_trans.append(ribo_rna_dict[gene]["tran"])
					elif rna_padj < (minzscore/100):
						if x > 0:
							rnaup_y_values.append(y)
							rnaup_x_values.append(x)
							rnaup_basemeans.append(basemean)
							rnaup_lfcses.append(lfcSE)
							rnaup_genes.append(gene)
							rnaup_trans.append(ribo_rna_dict[gene]["tran"])
						elif x < 0 :
							rnadown_y_values.append(y)
							rnadown_x_values.append(x)
							rnadown_basemeans.append(basemean)
							rnadown_lfcses.append(lfcSE)
							rnadown_genes.append(gene)
							rnadown_trans.append(ribo_rna_dict[gene]["tran"])
				else:
					highlight_y_values.append(y)
					highlight_x_values.append(x)
					highlight_basemeans.append(basemean)
					highlight_lfcses.append(lfcSE)
					highlight_genes.append(gene)
					highlight_trans.append(ribo_rna_dict[gene]["tran"])
		elif label == "Riboseq":
			ribo_padj = ribo_rna_dict[gene]["ribo_padj"]
			basemean = ribo_rna_dict[gene]["ribo_basemean"]
			lfcSE = ribo_rna_dict[gene]["ribo_lfcSE"]
			
			if "NA" in ribo_padj:
				ribo_padj = 1
			ribo_padj = float(ribo_padj)
			if gene not in gene_list:
				if ribo_padj > (minzscore/100):
					if ribo_rna_dict[gene]["ribo1"] != 0:
						y_values.append(log(ribo_rna_dict[gene]["ribo1"],2))
					else:
						y_values.append(0)
					if ribo_rna_dict[gene]["ribo2"] != 0:
						x_values.append(log(ribo_rna_dict[gene]["ribo2"],2))
					else:
						x_values.append(0)
					genes.append(gene)
					trans.append(ribo_rna_dict[gene]["tran"])
					basemeans.append(basemean)
					lfcses.append(lfcSE)
				else:
					if ribo_rna_dict[gene]["ribo1"] != 0:
						teup_y_values.append(log(ribo_rna_dict[gene]["ribo1"],2))
					else:
						teup_y_values.append(0)
						
					if ribo_rna_dict[gene]["ribo2"] != 0:
						teup_x_values.append(log(ribo_rna_dict[gene]["ribo2"],2))
					else:
						teup_x_values.append(0)
					teup_genes.append(gene)
					teup_trans.append(ribo_rna_dict[gene]["tran"])
					teup_basemeans.append(basemean)
					teup_lfcses.append(lfcSE)
			else:
				if ribo_rna_dict[gene]["ribo1"] != 0:
					highlight_y_values.append(log(ribo_rna_dict[gene]["ribo1"],2))
				else:
					highlight_y_values.append(0)
				if ribo_rna_dict[gene]["ribo2"] != 0:
					highlight_x_values.append(log(ribo_rna_dict[gene]["ribo2"],2))
				else:
					highlight_x_values.append(0)
				highlight_basemeans.append(basemean)
				highlight_lfcses.append(lfcSE)
				highlight_genes.append(gene)
				highlight_trans.append(ribo_rna_dict[gene]["tran"])
		elif label == "Rnaseq":
			rna_padj = ribo_rna_dict[gene]["rna_padj"]
			basemean = ribo_rna_dict[gene]["rna_basemean"]
			lfcSE = ribo_rna_dict[gene]["rna_lfcSE"]
			
			if "NA" in rna_padj:
				rna_padj = 1
			rna_padj = float(rna_padj)
			if gene not in gene_list:
				if rna_padj > (minzscore/100):
					basemeans.append(basemean)
					lfcses.append(lfcSE)
					if ribo_rna_dict[gene]["rna1"] != 0:
						y_values.append(log(ribo_rna_dict[gene]["rna1"],2))
					else:
						y_values.append(0)
					if ribo_rna_dict[gene]["rna2"] != 0:
						x_values.append(log(ribo_rna_dict[gene]["rna2"],2))
					else:
						x_values.append(0)

					genes.append(gene)
					trans.append(ribo_rna_dict[gene]["tran"])
				else:
					if ribo_rna_dict[gene]["rna1"] != 0:
						teup_y_values.append(log(ribo_rna_dict[gene]["rna1"],2))
					else:
						teup_y_values.append(0)
					if ribo_rna_dict[gene]["rna2"] != 0:
						teup_x_values.append(log(ribo_rna_dict[gene]["rna2"],2))
					else:
						teup_x_values.append(ribo_rna_dict[gene]["rna2"])
					teup_genes.append(gene)
					teup_trans.append(ribo_rna_dict[gene]["tran"])
					teup_basemeans.append(basemean)
					teup_lfcses.append(lfcSE)
			else:
				if ribo_rna_dict[gene]["rna1"] != 0:
					highlight_y_values.append(log(ribo_rna_dict[gene]["rna1"],2))
				else:
					highlight_y_values.append(0)
				if ribo_rna_dict[gene]["rna2"] != 0:
					highlight_x_values.append(log(ribo_rna_dict[gene]["rna2"],2))
				else:
					highlight_x_values.append(0)
				highlight_basemeans.append(basemean)
				highlight_lfcses.append(lfcSE)
				highlight_genes.append(gene)
				highlight_trans.append(ribo_rna_dict[gene]["tran"])
				
	source = ColumnDataSource({'x': x_values,'y':y_values,'trans':trans, 'genes':genes,'basemeans':basemeans,'lfcses':lfcses})
	if label == "TE":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq FC (log2)",  y_axis_label='Ribo-Seq FC (log2)',title="Ribo-Seq FC vs RNA-Seq FC ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Riboseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="Ribo-Seq Cond2 count (log2)",  y_axis_label='Ribo-Seq Cond1 count (log2)',title="Ribo-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Rnaseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq Cond2 count (log2)",  y_axis_label='RNA-Seq Cond1 count (log2)',title="RNA-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='gray')


	if highlight_genes != []:
		source = ColumnDataSource({'x':highlight_x_values, 'y':highlight_y_values, 'trans':highlight_trans,'genes':highlight_genes,'basemeans':highlight_basemeans,'lfcses':highlight_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#000000',legend="Highlighted Genes")
	
	if label == "TE":
		source = ColumnDataSource({'x':teup_x_values, 'y':teup_y_values, 'trans':teup_trans,'genes':teup_genes,'basemeans':teup_basemeans,'lfcses':teup_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Translation up")
		source = ColumnDataSource({'x':tedown_x_values, 'y':tedown_y_values, 'trans':tedown_trans,'genes':tedown_genes,'basemeans':tedown_basemeans,'lfcses':tedown_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#ff5050',legend="Translation down")

		source = ColumnDataSource({'x':rnaup_x_values, 'y':rnaup_y_values, 'trans':rnaup_trans,'genes':rnaup_genes,'basemeans':rnaup_basemeans,'lfcses':rnaup_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#cc0099',legend="mRNA up")
		source = ColumnDataSource({'x':rnadown_x_values, 'y':rnadown_y_values, 'trans':rnadown_trans,'genes':rnadown_genes,'basemeans':rnadown_basemeans,'lfcses':rnadown_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#ffff00',legend="mRNA down")
	elif label == "Riboseq":
		source = ColumnDataSource({'x':teup_x_values, 'y':teup_y_values, 'trans':teup_trans,'genes':teup_genes,'basemeans':teup_basemeans,'lfcses':teup_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Significant genes")
		source = ColumnDataSource({'x':tedown_x_values, 'y':tedown_y_values, 'trans':tedown_trans,'genes':tedown_genes,'basemeans':tedown_basemeans,'lfcses':tedown_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Significant genes")
	elif label == "Rnaseq":
		source = ColumnDataSource({'x':teup_x_values, 'y':teup_y_values, 'trans':teup_trans,'genes':teup_genes,'basemeans':teup_basemeans,'lfcses':teup_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Significant genes")
		source = ColumnDataSource({'x':rnadown_x_values, 'y':rnadown_y_values, 'trans':rnadown_trans,'genes':rnadown_genes,'basemeans':rnadown_basemeans,'lfcses':rnadown_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Significant genes")
	

	p.legend.location = "top_left"
	p.legend.click_policy="hide"
	p.legend.label_text_font_size = "28px"
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	plot_limit = largest_fc*1.1
	if label == "TE":
		p.line([plot_limit*-1,plot_limit], [plot_limit*-1,plot_limit], color="#cccccc",line_width=1)
		hover.tooltips = [("Ribo fc", "@y"),("RNA fc","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
	elif label == "Riboseq":
		p.line([0,16], [0,16], color="#cccccc",line_width=1)
		hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		#corr = spearmanr(x_values, y_values)
		#pearson_corr = pearsonr(x_values, y_values)
		#mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
		#p.add_layout(mytext)
	else:
		p.line([0,16], [0,16], color="#cccccc",line_width=1)
		hover.tooltips = [("Rna-seq Cond 1 count (log2)", "@y"),("Rna-seq Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		#corr = spearmanr(x_values, y_values)
		#pearson_corr = pearsonr(x_values, y_values)
		#mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
		#p.add_layout(mytext)
	
	file_string = ""
	label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"
	
	if riboseq1:
		if riboseq1[0] != "":
			for file_id in riboseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%23007a02_")
	if riboseq2:
		if riboseq2[0] != "":
			for file_id in riboseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23960000_")
	if rnaseq1:
		if rnaseq1[0] != "":
			for file_id in rnaseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%2374ed76_")
	if rnaseq2:
		if rnaseq2[0] != "":
			for file_id in rnaseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23ff6d6d_")

	# remove the trailing _ in file_string if it's been populated
	if file_string:
		file_string = file_string[:len(file_string)-1]

	if ambiguous == True:
		ambig = "T"
	else:
		ambig = "F"
	
	url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@trans&normalize={}&cov=T&ambig={}&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0],ambig)
	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)

	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download DESeq2 output files</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph









def anota2seq_plot(ribo_rna_dict,organism,transcriptome,riboseq1,riboseq2,rnaseq1,rnaseq2,background_col,short_code,normalized,filename,no_groups,title_size, axis_label_size, subheading_size,marker_size,ambiguous,gene_list,label,minzscore,sig_translated,sig_rna,sig_buffering):
	#Convert gene_list from string to list
	# ("gene list passed to deseq2_plot", gene_list)
	if gene_list != "":
		gene_list = gene_list.replace(","," ").replace("\t"," ")
		split_list = gene_list.split(" ")
		gene_list = []
		for item in split_list:
			gene_list.append(item.strip(" ").upper())
	x_values = []
	y_values = []
	basemeans = []
	lfcses = []
	genes = []
	trans = []
	
	highlight_x_values = []
	highlight_y_values = []
	highlight_basemeans = []
	highlight_lfcses = []
	highlight_genes = []
	highlight_trans = []
	
	teup_x_values = []
	teup_y_values = []
	teup_basemeans = []
	teup_lfcses = []
	teup_genes = []
	teup_trans = []
	tedown_x_values = []
	tedown_y_values = []
	tedown_basemeans = []
	tedown_lfcses = []
	tedown_genes = []
	tedown_trans = []	
	
	bufferedup_x_values = []
	bufferedup_y_values = []
	bufferedup_basemeans = []
	bufferedup_lfcses = []
	bufferedup_genes = []
	bufferedup_trans = []
	buffereddown_x_values = []
	buffereddown_y_values = []
	buffereddown_basemeans = []
	buffereddown_lfcses = []
	buffereddown_genes = []
	buffereddown_trans = []		
	
	rnaup_x_values = []
	rnaup_y_values = []
	rnaup_basemeans = []
	rnaup_lfcses = []
	rnaup_genes = []
	rnaup_trans = []
	rnadown_x_values = []
	rnadown_y_values = []
	rnadown_basemeans = []
	rnadown_lfcses = []
	rnadown_genes = []
	rnadown_trans = []	
	
	largest_fc = 0
	for gene in ribo_rna_dict:
		if label == "TE":
			x = ribo_rna_dict[gene]["rna_fc"]
			y = ribo_rna_dict[gene]["ribo_fc"]
			basemean = ribo_rna_dict[gene]["te_basemean"]
			lfcSE = ribo_rna_dict[gene]["te_lfcSE"]
			te_padj = ribo_rna_dict[gene]["te_padj"]
			ribo_padj = ribo_rna_dict[gene]["ribo_padj"]
			rna_padj = ribo_rna_dict[gene]["rna_padj"]
			
			if x != "NA" and y != "NA":
				x = float(x)
				y = float(y)
				if abs(x) > largest_fc:
					largest_fc = abs(x)
				if abs(y) > largest_fc:
					largest_fc = abs(y)
				if "NA" in te_padj:
					te_padj = 1
				if "NA" in ribo_padj:
					ribo_padj = 1
				if "NA" in rna_padj:
					rna_padj = 1
				te_padj = float(te_padj)
				ribo_padj = float(ribo_padj)
				rna_padj = float(rna_padj)
				if gene not in gene_list:
					if gene not in sig_translated and gene not in sig_rna and gene not in sig_buffering:
						y_values.append(y)
						x_values.append(x)
						basemeans.append(basemean)
						lfcses.append(lfcSE)
						genes.append(gene)
						trans.append(ribo_rna_dict[gene]["tran"])
					elif gene in sig_translated:
						if x < y :
							teup_y_values.append(y)
							teup_x_values.append(x)
							teup_basemeans.append(basemean)
							teup_lfcses.append(lfcSE)
							teup_genes.append(gene)
							teup_trans.append(ribo_rna_dict[gene]["tran"])
						elif y < x :
							tedown_y_values.append(y)
							tedown_x_values.append(x)
							tedown_basemeans.append(basemean)
							tedown_lfcses.append(lfcSE)
							tedown_genes.append(gene)
							tedown_trans.append(ribo_rna_dict[gene]["tran"])
					elif gene in sig_rna:
							if x > 0:
								rnaup_y_values.append(y)
								rnaup_x_values.append(x)
								rnaup_basemeans.append(basemean)
								rnaup_lfcses.append(lfcSE)
								rnaup_genes.append(gene)
								rnaup_trans.append(ribo_rna_dict[gene]["tran"])
							elif x < 0 :
								rnadown_y_values.append(y)
								rnadown_x_values.append(x)
								rnadown_basemeans.append(basemean)
								rnadown_lfcses.append(lfcSE)
								rnadown_genes.append(gene)
								rnadown_trans.append(ribo_rna_dict[gene]["tran"])
					elif gene in sig_buffering:
						#If no change in ribo, then this is buffered
						if x > 0:
							bufferedup_y_values.append(y)
							bufferedup_x_values.append(x)
							bufferedup_basemeans.append(basemean)
							bufferedup_lfcses.append(lfcSE)
							bufferedup_genes.append(gene)
							bufferedup_trans.append(ribo_rna_dict[gene]["tran"])
						else:
							buffereddown_y_values.append(y)
							buffereddown_x_values.append(x)
							buffereddown_basemeans.append(basemean)
							buffereddown_lfcses.append(lfcSE)
							buffereddown_genes.append(gene)
							buffereddown_trans.append(ribo_rna_dict[gene]["tran"])
				else:
					highlight_y_values.append(y)
					highlight_x_values.append(x)
					highlight_basemeans.append(basemean)
					highlight_lfcses.append(lfcSE)
					highlight_genes.append(gene)
					highlight_trans.append(ribo_rna_dict[gene]["tran"])
	source = ColumnDataSource({'x': x_values,'y':y_values,'trans':trans, 'genes':genes,'basemeans':basemeans,'lfcses':lfcses})
	if label == "TE":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq FC (log2)",  y_axis_label='Ribo-Seq FC (log2)',title="Ribo-Seq FC vs RNA-Seq FC ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Riboseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="Ribo-Seq Cond2 count (log2)",  y_axis_label='Ribo-Seq Cond1 count (log2)',title="Ribo-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	elif label == "Rnaseq":
		p = figure(plot_width=1300, plot_height=1300,x_axis_label="RNA-Seq Cond2 count (log2)",  y_axis_label='RNA-Seq Cond1 count (log2)',title="RNA-Seq correlation ({})".format(short_code),toolbar_location="below",
			tools = "reset,pan,box_zoom,save,hover,tap")
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_col
	p.xgrid.grid_line_color = "#cccccc"
	p.ygrid.grid_line_color = "#cccccc"
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='gray')

	
	if highlight_genes != []:
		source = ColumnDataSource({'x':highlight_x_values, 'y':highlight_y_values, 'trans':highlight_trans,'genes':highlight_genes,'basemeans':highlight_basemeans,'lfcses':highlight_lfcses})
		p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#000000',legend="Highlighted Genes")
	
	source = ColumnDataSource({'x':teup_x_values, 'y':teup_y_values, 'trans':teup_trans,'genes':teup_genes,'basemeans':teup_basemeans,'lfcses':teup_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#00ff99',legend="Translation up")
	source = ColumnDataSource({'x':tedown_x_values, 'y':tedown_y_values, 'trans':tedown_trans,'genes':tedown_genes,'basemeans':tedown_basemeans,'lfcses':tedown_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#ff5050',legend="Translation down")
	
	source = ColumnDataSource({'x':bufferedup_x_values, 'y':bufferedup_y_values, 'trans':bufferedup_trans,'genes':bufferedup_genes,'basemeans':bufferedup_basemeans,'lfcses':bufferedup_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#6699ff',legend="Buffered up")
	source = ColumnDataSource({'x':buffereddown_x_values, 'y':buffereddown_y_values, 'trans':buffereddown_trans,'genes':buffereddown_genes,'basemeans':buffereddown_basemeans,'lfcses':buffereddown_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#ffcc66',legend="Buffered down")
	
	source = ColumnDataSource({'x':rnaup_x_values, 'y':rnaup_y_values, 'trans':rnaup_trans,'genes':rnaup_genes,'basemeans':rnaup_basemeans,'lfcses':rnaup_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#cc0099',legend="mRNA abundance up")
	source = ColumnDataSource({'x':rnadown_x_values, 'y':rnadown_y_values, 'trans':rnadown_trans,'genes':rnadown_genes,'basemeans':rnadown_basemeans,'lfcses':rnadown_lfcses})
	p.scatter('x','y', alpha=0.2,color="black",fill_alpha=1,size=12,source=source,fill_color='#ffff00',legend="mRNA abundance down")
	
	#legend_label=name)

	p.legend.location = "top_left"
	p.legend.click_policy="hide"
	p.legend.label_text_font_size = "28px"
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'
	plot_limit = largest_fc*1.1
	if label == "TE":
		p.line([plot_limit*-1,plot_limit], [plot_limit*-1,plot_limit], color="#cccccc",line_width=1)
		hover.tooltips = [("Ribo fc", "@y"),("RNA fc","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
	elif label == "Riboseq":
		p.line([0,16], [0,16], color="#cccccc",line_width=1)
		hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		#corr = spearmanr(x_values, y_values)
		#pearson_corr = pearsonr(x_values, y_values)
		#mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
		#p.add_layout(mytext)
	else:
		p.line([0,16], [0,16], color="#cccccc",line_width=1)
		hover.tooltips = [("Rna-seq Cond 1 count (log2)", "@y"),("Rna-seq Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		hover.tooltips = [("Ribo Cond 1 count (log2)", "@y"),("Ribo Cond 2 count (log2)","@x"),("Genes","@genes"),("Transcript","@trans"),("Basemean","@basemeans"),("lfcSE","@lfcses")]
		#corr = spearmanr(x_values, y_values)
		#pearson_corr = pearsonr(x_values, y_values)
		#mytext = Label(x=0.1,y=max(y_values),text="Pearson correlation: {}".format(round(pearson_corr[0],2)),background_fill_color="white",text_font_size="13pt")
		#p.add_layout(mytext)
	
	file_string = ""
	label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"
	
	if riboseq1:
		if riboseq1[0] != "":
			for file_id in riboseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%23007a02_")
	if riboseq2:
		if riboseq2[0] != "":
			for file_id in riboseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23960000_")
	if rnaseq1:
		if rnaseq1[0] != "":
			for file_id in rnaseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%2374ed76_")
	if rnaseq2:
		if rnaseq2[0] != "":
			for file_id in rnaseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23ff6d6d_")

	# remove the trailing _ in file_string if it's been populated
	if file_string:
		file_string = file_string[:len(file_string)-1]

	if ambiguous == True:
		ambig = "T"
	else:
		ambig = "F"
	
	url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@trans&normalize={}&cov=T&ambig={}&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0],ambig)
	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)

	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download DESeq2 output files</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)
	#layout = column(text_input, p)
	graph += file_html(p,CDN)
	return graph































