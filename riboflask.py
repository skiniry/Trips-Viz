import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
import mpld3
from mpld3 import plugins,utils
import collections
from sqlitedict import SqliteDict
import pandas as pd
from fetch_shelve_reads2 import get_reads,get_seq_var,get_readlength_breakdown
import sqlite3
import os
import config
from new_plugins import InteractiveLegendPlugin,PointHTMLTooltip,TopToolbar,DownloadProfile,DownloadPNG


# CSS for popup tables that appear when hovering over aug codons
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

color_dict = {'frames': ['#FF4A45', '#64FC44', '#5687F9']}

def get_user_defined_seqs(seq,seqhili):
	iupac_dict = {"A":["A"],"U":["U"],"G":["G"],"C":["C"],"R":["A","G"],"Y":["C","U"],"S":["G","C"],"W":["A","U"],"K":["G","U"],
			   "M":["A","C"],"B":["C","G","U"],"D":["A","G","U"],"D":["A","G","U"],"H":["A","C","U"],"V":["A","C","G"],"N":["A","U","G","C"]}
	signalhtml = {0:[],1:[],2:[]}
	seq = seq.replace("T","U")
	near_cog_starts = {0:[],1:[],2:[]}
	for i in range(0,len(seq)):
		for subseq in seqhili:
			subseq = subseq.upper()
			subseq = subseq.replace("T","U").replace(" ","")
			partial_seq = list(seq[i:i+len(subseq)])
			if len(partial_seq) != len(subseq):
				continue
			x= 0
			for x in range(0,len(subseq)):
				char = subseq[x]
				if partial_seq[x] in iupac_dict[char]:
					partial_seq[x] = char
			partial_seq = "".join(partial_seq)
			if partial_seq == subseq:
				near_cog_starts[(i)%3].append(i+1)
				datadict = {'sequence': [subseq]}
				df = pd.DataFrame(datadict, columns=(["sequence"]))
				label = df.ix[[0], :].T
				label.columns = ["Position: {}".format(i)]
				signalhtml[(i)%3].append(str(label.to_html()))
	return near_cog_starts,signalhtml


def merge_dicts(dict1,dict2):
	for nuc in dict2:
		#print "nuc", nuc
		#print dict1
		#print dict2[nuc]
		if nuc not in dict1:
			dict1[nuc] = dict2[nuc]
		else:
			for pos in dict2[nuc]:
				if pos not in dict1[nuc]:
					dict1[nuc][pos] = dict2[nuc][pos]
				else:
					dict1[nuc][pos] += dict2[nuc][pos]
	return dict1

def generate_plot(tran, ambig, min_read, max_read,lite,ribocoverage,organism,readscore, noisered, primetype, minfiles,nucseq, user_hili_starts, user_hili_stops,uga_diff,file_paths_dict, short_code, color_readlen_dist, background_col,uga_col, uag_col, uaa_col,advanced,trips_annotation_location,seqhili,seq_rules,title_size,
subheading_size,axis_label_size,marker_size, transcriptome, trips_uploads_location,cds_marker_size,cds_marker_colour,legend_size,ribo_linewidth, secondary_readscore,pcr,mismatches, hili_start, hili_stop):
	
	if lite == "n" and ribocoverage == True:
		return "Error: Cannot display Ribo-Seq Coverage when 'Line Graph' is turned off"
	labels = ["Frame 1 profiles","Frame 2 profiles","Frame 3 profiles","RNA", "Exon Junctions"]
	start_visible=[True, True, True, True, False]
	if mismatches == True:
		labels.append("Mismatches A")
		labels.append("Mismatches T")
		labels.append("Mismatches G")
		labels.append("Mismatches C")
		start_visible.append(False)
		start_visible.append(False)
		start_visible.append(False)
		start_visible.append(False)
	start_visible.append(True)
	labels.append("CDS markers")
	#This is a list of booleans that decide if the interactive legends boxes are filled in or not.Needs to be same length as labels
	stop_codons = ["TAG","TAA","TGA"]
	frame_orfs = {1:[],2:[],3:[]}
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
	if owner == 1:
		if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
			transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
		else:
			return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	else:
		transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(trips_uploads_location,owner,organism,transcriptome))
	cursor = transhelve.cursor()
	cursor.execute("SELECT * from transcripts WHERE transcript = '{}'".format(tran))
	result = cursor.fetchone()
	traninfo = {"transcript":result[0] , "gene":result[1], "length":result[2] , "cds_start":result[3] , "cds_stop":result[4] , "seq":result[5] ,
			 "strand":result[6], "stop_list":result[7].split(","),"start_list":result[8].split(","), "exon_junctions":result[9].split(","),
			 "tran_type":result[10], "principal":result[11]}
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

	all_cds_regions = []
	# Check if the 'coding_regions' table exists
	cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='coding_regions';")
	result = cursor.fetchone()
	#print ("CODING REGION RESULT",result)
	if result != None:
		print ("result is not empty")
		cursor.execute("SELECT * from coding_regions WHERE transcript = '{}'".format(tran))
		result = cursor.fetchall()
		for row in result:
			all_cds_regions.append((row[1],row[2]))
	transhelve.close()
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
	all_rna_reads, rna_seqvar_dict = get_reads(ambig, min_read, max_read, tran, file_paths_dict,tranlen,True, organism, False,noisered, primetype,"rnaseq",readscore,pcr,get_mismatches=mismatches)

	all_subcodon_reads,ribo_seqvar_dict = get_reads(ambig, min_read, max_read, tran, file_paths_dict,tranlen,ribocoverage, organism, True,noisered, primetype,"riboseq",readscore,secondary_readscore,pcr,get_mismatches=mismatches)
	#print "riboseq var dict",ribo_seqvar_dict
	seq_var_dict = merge_dicts(ribo_seqvar_dict, rna_seqvar_dict)
	#print "all subcodon reads", all_subcodon_reads
	try:
		rnamax = max(all_rna_reads.values())
	except:
		rnamax = 0
	try:
		subcodonmax = max(all_subcodon_reads.values())
	except:
		subcodonmax = 0

	y_max = max(1,rnamax, subcodonmax)*1.1
	fig = plt.figure(figsize=(23,12))
	ax_main = plt.subplot2grid((30,1), (0,0),rowspan=22)
	ax_main.spines['bottom'].set_visible(False)
	for s in ['bottom', 'left','top','right']:
		ax_main.spines[s].set_linewidth(15)
		ax_main.spines[s].set_color("red")
	alt_seq_type_vars = []
	# Plot any alternative sequence types if there are any
	for seq_type in file_paths_dict:
		if seq_type != "riboseq" and seq_type != "rnaseq":
			#print "seq_type", seq_type
			if seq_rules[seq_type]["frame_breakdown"] == 1:
				frame_breakdown = True
			else:
				frame_breakdown = False
			alt_sequence_reads,empty_seqvar_dict = get_reads(ambig, min_read, max_read, tran, file_paths_dict,tranlen,True, organism, frame_breakdown,noisered, primetype,seq_type,readscore)

			if frame_breakdown == False:
				alt_seq_plot = ax_main.plot(alt_sequence_reads.keys(), alt_sequence_reads.values(), alpha=1, label = seq_type, zorder=2, color='#5c5c5c', linewidth=2)
				labels.append(seq_type)
				start_visible.append(True)
				alt_seq_type_vars.append(alt_seq_plot)
			else:
				alt_frame_counts = {0: collections.OrderedDict(), 1: collections.OrderedDict(), 2: collections.OrderedDict()}
				for key in alt_sequence_reads:
					start = key
					rem = start % 3
					if rem == 1:    # frame 1
						frame = 2
					elif rem == 2:  # frame 2
						frame = 0
					elif rem == 0:  # frame 3
						frame = 1
					alt_frame_counts[frame][key] = alt_sequence_reads[key]
				frame0_altseqplot = ax_main.plot(alt_frame_counts[0].keys(), alt_frame_counts[0].values(), alpha=0.75, label = seq_type+"frame0", zorder=2, color= "#FF4A45", linewidth=2)
				frame1_altseqplot = ax_main.plot(alt_frame_counts[1].keys(), alt_frame_counts[1].values(), alpha=0.75, label = seq_type+"frame1", zorder=2, color= "#64FC44", linewidth=2)
				frame2_altseqplot = ax_main.plot(alt_frame_counts[2].keys(), alt_frame_counts[2].values(), alpha=0.75, label = seq_type+"frame2*", zorder=2, color= "#5687F9", linewidth=2)
				labels.append(seq_type+"frame 1")
				labels.append(seq_type+"frame 2")
				labels.append(seq_type+"frame 3")
				start_visible.append(True)
				start_visible.append(True)
				start_visible.append(True)
				alt_seq_type_vars.append(frame0_altseqplot)
				alt_seq_type_vars.append(frame1_altseqplot)
				alt_seq_type_vars.append(frame2_altseqplot)
			if max(alt_sequence_reads.values()) > y_max:
				y_max = max(alt_sequence_reads.values())

	label = 'Read count'
	ax_main.set_ylabel(label,  fontsize=axis_label_size, labelpad=30)
	label = 'Position (nucleotides)'
	ax_main.set_xlabel(label, fontsize=axis_label_size,labelpad=-10)
	ax_main.set_ylim(0, y_max)

	if lite == "n":
		rna_bars = ax_main.bar(all_rna_reads.keys(), all_rna_reads.values(), alpha=1, label = labels, zorder=1,color='lightgray', linewidth=0, width=1)
	else:
		rna_bars = ax_main.plot(all_rna_reads.keys(), all_rna_reads.values(), alpha=1, label = labels, zorder=1,color='#a7adb7', linewidth=4)

	cds_markers = ax_main.plot((cds_start,cds_start), (0, y_max*0.97), color=cds_marker_colour,linestyle = 'solid', linewidth=cds_marker_size)
	ax_main.text(cds_start,y_max*0.97,"CDS start",fontsize=18,color="black",ha="center")
	#ax_main.annotate('axes fraction',xy=(3, 1), xycoords='data',xytext=(0.8, 0.95), textcoords='axes fraction',arrowprops=dict(facecolor='black', shrink=0.05),horizontalalignment='right', verticalalignment='top')
	#trans = blended_transform_factory(ax_main.transData, ax_main.transAxes)
	#ax_main.annotate('CDS RELATIVE START',(100,100),transform=trans)
	#tform = blended_transform_factory(ax_main.transData, ax_main.transAxes)
	#r=10
	#ax_main.text(cds_start, 0.9, "CDS START OR WHATEVER", fontsize='xx-large', color='r', transform=tform)
	cds_markers += ax_main.plot((cds_stop+1,cds_stop+1), (0, y_max*0.97), color=cds_marker_colour,linestyle = 'solid', linewidth=cds_marker_size)
	ax_main.text(cds_stop,y_max*0.97,"CDS stop",fontsize=18,color="black",ha="center")
	ax_cds = plt.subplot2grid((31,1), (26,0),rowspan=1,sharex=ax_main)
	ax_cds.set_axis_bgcolor("white")
	ax_cds.set_ylabel('Merged CDS', labelpad=4, verticalalignment='center',horizontalalignment="right",rotation="horizontal",color="black",fontsize=axis_label_size-5)
	ax_f1 = plt.subplot2grid((31,1), (27,0),rowspan=1,sharex=ax_main)
	ax_f1.set_axis_bgcolor(color_dict['frames'][0])
	ax_f2 = plt.subplot2grid((31,1), (28,0),rowspan=1,sharex=ax_main)
	ax_f2.set_axis_bgcolor(color_dict['frames'][1])
	ax_f3 = plt.subplot2grid((31,1), (29,0),rowspan=1,sharex=ax_main)
	ax_f3.set_axis_bgcolor(color_dict['frames'][2])
	ax_nucseq = plt.subplot2grid((31,1), (30,0),rowspan=1,sharex=ax_main)
	ax_nucseq.set_xlabel('Transcript: {} Length: {} nt'.format(tran, tranlen), fontsize=subheading_size)


	for tup in all_cds_regions:
		ax_cds.fill_between([tup[0],tup[1]], [1, 1],zorder=0, alpha=1, color="#001285")


	#plot a dummy exon junction at postion -1, needed in cases there are no exon junctions, this wont be seen
	allexons = ax_main.plot((-1,-1), (0, 1), alpha=0.01,color='black',linestyle = '-.', linewidth=2)
	for exon in exon_junctions:
		allexons += ax_main.plot((exon,exon), (0, y_max), alpha=0.01,color='black',linestyle = '-.', linewidth=3)

	#dictionary for each frame in which the keys are the posistions and the values are the counts
	frame_counts = {0: collections.OrderedDict(), 1: collections.OrderedDict(), 2: collections.OrderedDict()}
	for key in all_subcodon_reads:
		rem = key % 3
		if rem == 1:    # frame 1
			frame = 2
		elif rem == 2:  # frame 2
			frame = 0
		elif rem == 0:  # frame 3
			frame = 1
		frame_counts[frame][key] = all_subcodon_reads[key]
		if lite == "n":
			frame_counts[frame][key+1] = 0
			frame_counts[frame][key+2] = 0

	if lite == "n":
		frame0subpro = ax_main.bar(frame_counts[0].keys(), frame_counts[0].values(), alpha=0.75, label = labels, zorder=2, color= "#FF4A45", width=1, linewidth=0)
		frame1subpro = ax_main.bar(frame_counts[1].keys(), frame_counts[1].values(), alpha=0.75, label = labels, zorder=2, color= "#64FC44", width=1, linewidth=0)
		frame2subpro = ax_main.bar(frame_counts[2].keys(), frame_counts[2].values(), alpha=0.75, label = labels, zorder=2, color= "#5687F9", width=1, linewidth=0)
	else:
		frame0subpro = ax_main.plot(frame_counts[0].keys(), frame_counts[0].values(), alpha=0.75, label = labels, zorder=2, color= "#FF4A45", linewidth=ribo_linewidth)
		frame1subpro = ax_main.plot(frame_counts[1].keys(), frame_counts[1].values(), alpha=0.75, label = labels, zorder=2, color= "#64FC44", linewidth=ribo_linewidth)
		frame2subpro = ax_main.plot(frame_counts[2].keys(), frame_counts[2].values(), alpha=0.75, label = labels, zorder=2, color= "#5687F9", linewidth=ribo_linewidth)
	if mismatches == True:
		a_mismatches = ax_main.plot(seq_var_dict["A"].keys(), seq_var_dict["A"].values(),alpha=0.01, label = labels, zorder=2, color= "purple", linewidth=2)
		t_mismatches = ax_main.plot(seq_var_dict["T"].keys(), seq_var_dict["T"].values(),alpha=0.01, label = labels, zorder=2, color= "yellow", linewidth=2)
		g_mismatches = ax_main.plot(seq_var_dict["G"].keys(), seq_var_dict["G"].values(),alpha=0.01, label = labels, zorder=2, color= "orange", linewidth=2)
		c_mismatches = ax_main.plot(seq_var_dict["C"].keys(), seq_var_dict["C"].values(),alpha=0.01, label = labels, zorder=2, color= "pink", linewidth=2)

	xy = 0
	if nucseq == True:
		ax_nucseq.set_axis_bgcolor(background_col)
		mrnaseq = seq.replace("T","U")
		color_list = ["#FF4A45","#64FC44","#5687F9"]
		char_frame = 0
		for char in mrnaseq:
			ax_nucseq.text((xy+1)-0.1,0.2,mrnaseq[xy],fontsize=20,color=color_list[char_frame%3])
			xy += 1
			char_frame += 1

	# If the user passed a list of sequences to highlight, find and plot them here.
	if seqhili != ['']:
		near_cog_starts,signalhtml = get_user_defined_seqs(seq, seqhili)
		for slip in near_cog_starts[0]:
			try:
				hili_sequences += ax_f1.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4,color='black', linewidth=5)
			except Exception as e:
				hili_sequences = ax_f1.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4, color='black', linewidth=5)
		for slip in near_cog_starts[1]:
			try:
				hili_sequences += ax_f2.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4,color='black', linewidth=5)
			except:
				hili_sequences = ax_f2.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4,color='black',linewidth=5)
		for slip in near_cog_starts[2]:
			try:
				hili_sequences += ax_f3.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4,color='black', linewidth=5)
			except:
				hili_sequences = ax_f3.plot((slip, slip),(0,0.5), alpha=1, label = labels, zorder=4,color='black', linewidth=5)

		#Plot sequence identifiers which will create a popup telling user what the subsequence is (useful if they have passed multiple subsequences)
		frame1_subsequences = ax_f1.plot(near_cog_starts[0], [0.25]*len(near_cog_starts[0]), 'o', color='b',mec='k', ms=12, mew=1, alpha=0, zorder=4)
		frame2_subsequences = ax_f2.plot(near_cog_starts[1], [0.25]*len(near_cog_starts[1]), 'o', color='b',mec='k', ms=12, mew=1, alpha=0, zorder=4)
		frame3_subsequences = ax_f3.plot(near_cog_starts[2], [0.25]*len(near_cog_starts[2]), 'o', color='b',mec='k', ms=12, mew=1, alpha=0, zorder=4)

		#Attach the labels to the subsequences plotted above
		signaltooltip1 = PointHTMLTooltip(frame1_subsequences[0], signalhtml[0], voffset=10, hoffset=10, css=point_tooltip_css)
		signaltooltip2 = PointHTMLTooltip(frame2_subsequences[0], signalhtml[1], voffset=10, hoffset=10, css=point_tooltip_css)
		signaltooltip3 = PointHTMLTooltip(frame3_subsequences[0], signalhtml[2], voffset=10, hoffset=10, css=point_tooltip_css)
	for axisname in (ax_f1, ax_f2, ax_f3,ax_nucseq,ax_cds):
		axisname.tick_params(top=False, bottom=False, labelleft=False, labelright=False, labelbottom=False)
	for label in ax_main.xaxis.get_majorticklabels():
		label.set_fontsize(36)
	for axis, frame in ((ax_f1, 1), (ax_f2, 2), (ax_f3, 3)):
		axis.set_xlim(1, tranlen)
		starts = [(item, 1) for item in start_stop_dict[frame]['starts']]
		uag_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TAG']]
		uaa_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TAA']]
		uga_stops = [(item, 1) for item in start_stop_dict[frame]['stops']['TGA']]
		axis.broken_barh(starts, (0.5, 1),color="white", zorder=2)
		axis.broken_barh(uag_stops, (0, 1), color=uag_col, zorder=2, linewidth=2)
		axis.broken_barh(uaa_stops, (0, 1), color=uaa_col, zorder=2, linewidth=2)
		axis.broken_barh(uga_stops, (0, 1), color=uga_col, zorder=2, linewidth=2)
		axis.set_ylim(0, 1)
		axis.set_ylabel('Frame {}'.format(frame), labelpad=4, verticalalignment='center',horizontalalignment="right",rotation="horizontal",color="black",fontsize=axis_label_size-5)
	title_str = '{} ({})'.format(gene,short_code)
	plt.title(title_str, fontsize=title_size,y=38)
	line_collections = [frame0subpro, frame1subpro, frame2subpro, rna_bars, allexons]

	if mismatches == True:
		line_collections.append(a_mismatches)
		line_collections.append(t_mismatches)
		line_collections.append(g_mismatches)
		line_collections.append(c_mismatches)
	line_collections.append(cds_markers)

	if not (hili_start == 0 and hili_stop == 0):
		hili_start = int(hili_start)
		hili_stop = int(hili_stop)
		hili = ax_main.fill_between([hili_start,hili_stop], [y_max, y_max],zorder=0, alpha=0.75, color="#fffbaf")
		labels.append("Highligted region")
		start_visible.append(True)
		line_collections.append(hili)

	for alt_plot in alt_seq_type_vars:
		line_collections.append(alt_plot)
	if 'hili_sequences' in locals():
		labels.append("Highligted sequences")
		start_visible.append(True)
		line_collections.append(hili_sequences)
	if user_hili_starts != [] and user_hili_stops != []:
		for i in range(0,len(user_hili_starts)):
			user_hili_start = int(user_hili_starts[i])
			user_hili_stop = int(user_hili_stops[i])
			try:
				hili += ax_main.fill_between([user_hili_start,user_hili_stop],[y_max, y_max], alpha=0.75,color="#fffbaf")
			except:
				hili = ax_main.fill_between([user_hili_start,user_hili_stop],[y_max, y_max], alpha=0.75,color="#fffbaf")
		labels.append("Highligter")
		start_visible.append(True)
		line_collections.append(hili)

	leg_offset = (legend_size-17)*5
	if leg_offset <0:
		leg_offset = 0

	ilp = InteractiveLegendPlugin(line_collections, labels, alpha_unsel=0,alpha_sel=0.85,start_visible=start_visible,fontsize=legend_size,xoffset=leg_offset)
	htmllabels = {1:[],2:[],3:[]}
	all_start_points = {1:[],2:[],3:[]}
	try:
		con_scores = SqliteDict("{0}homo_sapiens/score_dict.sqlite".format(trips_annotation_location))
	except Exception as e:
		print "Couldn't open conservation scores "+e
		con_scores = []
	for frame in [1,2,3]:
		orf_list = frame_orfs[frame]
		for tup in orf_list:
			orf_ribo = 0.0
			outframe_ribo = 0.0
			orf_rna = 0.0001
			start = tup[0]
			try:
				context = (seq[start-7:start+4].upper()).replace("T","U")
			except Exception as e:
				con_score = "?"
			if len(context) != 11 or context[6:9] != "AUG":
				con_score = "?"
			else:
				try:
					con_score = con_scores[context.upper()]
				except Exception as e:
					con_score = "?"
			all_start_points[frame].append(start-1)
			stop = tup[1]
			other_ribo = 0.0
			otherother_ribo = 0.0
			for i in range(start+2, stop,3):
				for subframe in [0,1,2]:
					if i in frame_counts[subframe]:
						orf_ribo += frame_counts[subframe][i]

			for i in range(start, stop,3):
				for subframe in [0,1,2]:
					if i in frame_counts[subframe]:
						outframe_ribo += frame_counts[subframe][i]

			for i in range(start+1, stop,3):
				for subframe in [0,1,2]:
					if i in frame_counts[subframe]:
						outframe_ribo += frame_counts[subframe][i]

			for i in range(start, stop+1):
				if i in all_rna_reads:
					orf_rna += all_rna_reads[i]

			orf_te = float(orf_ribo)/float(orf_rna)
			orf_len = int(stop-start)

			try:
				in_out_ratio = orf_ribo/outframe_ribo
			except:
				in_out_ratio = "Null"

			datadict = {'inframe ribo': [orf_ribo],
						'outframe ribo':[outframe_ribo],
						'in/out ratio':[in_out_ratio],
						'rna':  [orf_rna],
						'te':   [orf_te],
						'len':  [orf_len],
						'context_score':[str(con_score)+"/150"]}
			df = pd.DataFrame(datadict, columns=(["inframe ribo", "outframe ribo", "in/out ratio","rna","te","len","context_score"]))
			label = df.ix[[0], :].T
			label.columns = ["Start pos: {}".format(start-1)]
			htmllabels[frame].append(str(label.to_html()))

	points1 =ax_f1.plot(all_start_points[1], [0.75]*len(all_start_points[1]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)
	points2 =ax_f2.plot(all_start_points[2], [0.75]*len(all_start_points[2]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)
	points3 =ax_f3.plot(all_start_points[3], [0.75]*len(all_start_points[3]), 'o', color='b',mec='k', ms=13, mew=1, alpha=0, zorder=3)

	tooltip1 = PointHTMLTooltip(points1[0], htmllabels[1],voffset=10, hoffset=10, css=point_tooltip_css)
	tooltip2 = PointHTMLTooltip(points2[0], htmllabels[2],voffset=10, hoffset=10, css=point_tooltip_css)
	tooltip3 = PointHTMLTooltip(points3[0], htmllabels[3],voffset=10, hoffset=10, css=point_tooltip_css)

	ax_f3.axes.get_yaxis().set_ticks([])
	ax_f2.axes.get_yaxis().set_ticks([])
	ax_f1.axes.get_yaxis().set_ticks([])

	returnstr = "Position,Sequence,Frame 1,Frame 2,Frame 3,RNA-Seq\n"
	for i in range(0,len(seq)):
		f1_count = 0
		f2_count = 0
		f3_count = 0
		rna_count = 0
		if i+1 in frame_counts[0]:
			f1_count = frame_counts[0][i+1]
		elif i+1 in frame_counts[1]:
			f2_count = frame_counts[1][i+1]
		elif i+1 in frame_counts[2]:
			f3_count = frame_counts[2][i+1]
		if i+1 in all_rna_reads:
			rna_count = all_rna_reads[i+1]
		returnstr += "{},{},{},{},{},{}\n".format(i+1,seq[i],f1_count, f2_count, f3_count,rna_count)

	if seqhili == ['']:
		plugins.connect(fig, ilp, tooltip1, tooltip2, tooltip3, TopToolbar(yoffset=100),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))
	else:
		plugins.connect(fig, ilp, tooltip1, tooltip2, tooltip3, signaltooltip1,signaltooltip2,signaltooltip3, TopToolbar(yoffset=100),DownloadProfile(returnstr=returnstr),DownloadPNG(returnstr=title_str))

	ax_main.set_axis_bgcolor(background_col)
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
