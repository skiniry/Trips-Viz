from flask import Blueprint, render_template, abort, request
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import time
import re
import operator
from math import log
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,build_profile
import metainfo_plots
import collections
from flask_login import current_user
import subprocess




metainfo_plotpage_blueprint = Blueprint("metainfo_plotpage", __name__, template_folder="templates")
@metainfo_plotpage_blueprint.route('/<organism>/<transcriptome>/metainfo_plot/')
def metainfo_plotpage(organism, transcriptome):
	global user_short_passed
	user_short_passed = True
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except Exception as e:
		print "Error ", e
		user = None
	print "user", user
	organism = str(organism)
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict, accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)

	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())
	studyinfo_dict = fetch_study_info(organism)
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	
	
	
	
	
	
	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				"user_plot_type":str(request.args.get('plot')),
				"nuc_comp_direction":str(request.args.get('nc_dir')),
				"nuc_comp_type":str(request.args.get('nc_type')),
				"nuc_comp_min_readlen":str(request.args.get('nc_minreadlen')),
				"nuc_comp_max_readlen":str(request.args.get('nc_maxreadlen')),
				"trip_periodicity_min_readlen":str(request.args.get('tp_minreadlen')),
				"trip_periodicity_max_readlen":str(request.args.get('tp_maxreadlen')),
				"heatmap_dir":str(request.args.get('hm_dir')),
				"heatmap_log_scale":str(request.args.get('hm_log')),
				"heatmap_reverse":str(request.args.get('hm_rev')),
				"heatmap_position":str(request.args.get('hm_pos')),
				"heatmap_minreadlen":str(request.args.get('hm_minreadlen')),
				"heatmap_maxreadlen":str(request.args.get('hm_maxreadlen')),
				"heatmap_start":str(request.args.get('hm_start')),
				"heatmap_stop":str(request.args.get('hm_stop')),
				"heatmap_colour":str(request.args.get('hm_col')),
				"metagene_pos":str(request.args.get('mg_pos')),
				"metagene_minreadlen":str(request.args.get('mg_minreadlen')),
				"metagene_maxreadlen":str(request.args.get('mg_maxreadlen')),
				"replicate_minreads":str(request.args.get('rp_minreads')),
				"mrna_dist_readlen_per":str(request.args.get('mdr_per')),
				"transcriptome":str(transcriptome),
				"maxscaleval":str(request.args.get('maxscaleval')),
				"mrna_dist_readlen_smooth":str(request.args.get('mdr_smooth')),
				"te_minreads":str(request.args.get('te_minreads')),
				"include_first":str(request.args.get('include_first')),
				"include_last":str(request.args.get('include_last')),
				"exclude_first":str(request.args.get('exclude_first')),
				"exclude_last":str(request.args.get('exclude_last')),
				"custom_seq_list":str(request.args.get('custom_seq_list')),
				"exclude_first_val":str(request.args.get('exclude_first_val')),
				"exclude_last_val":str(request.args.get('exclude_last_val')),
				"include_first_val":str(request.args.get('include_first_val')),
				"include_last_val":str(request.args.get('include_last_val')),
				"metagene_tranlist":str(request.args.get('metagene_tranlist')),
				"te_tranlist":str(request.args.get('te_tranlist'))}
	user_files = request.args.get('files')
	if user_files != None:
		user_files = user_files.split(",")
		html_args["user_files"] = [str(x) for x in user_files]
	else:
		html_args["user_files"] = []
	#print"html args", html_args
	user_ribo_studies = request.args.get('ribo_studies')
	if user_ribo_studies != None:
		user_ribo_studies = user_ribo_studies.split(",")
		html_args["user_ribo_studies"] = [str(x) for x in user_ribo_studies]
	else:
		html_args["user_ribo_studies"] = []
	user_rna_studies = request.args.get('rna_studies')
	if user_rna_studies != None:
		user_rna_studies = user_rna_studies.split(",")
		html_args["user_rna_studies"] = [str(x) for x in user_rna_studies]
	else:
		html_args["user_rna_studies"] = []
	connection.close()
	return render_template('metainfo_index.html', gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,transcriptome=transcriptome,organism=organism,default_tran=default_tran,current_username=user,local=local,
						   studies_dict=accepted_studies,accepted_files=accepted_files,html_args=html_args,user_files=user_files,user_ribo_studies=user_ribo_studies, user_rna_studies=user_rna_studies,
						   studyinfo_dict=studyinfo_dict,seq_types=seq_types)



# Used to create custom metagene plots on the metainformation plot page
def create_custom_metagene(custom_seq_list, exclude_first_val, exclude_last_val, include_first_val, include_last_val, custom_search_region, exclude_first, exclude_last, include_first, include_last,sqlite_db, organism,metagene_tranlist,metagene_frame,transcriptome):
	print "create custom metagene called"
	custom_metagene_id = "cmgc_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}".format((custom_seq_list.upper()).replace(" ","").replace("T","U").replace(",","_"), custom_search_region, exclude_first, exclude_last, include_first, include_last,exclude_first_val,exclude_last_val,include_first_val,include_last_val,metagene_tranlist)
	#try:
	#	mgc = sqlite_db[custom_metagene_id]
	#except:
	#	pass
	offsets = sqlite_db["offsets"]
	#print"offsets",offsets
	if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
		transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
	else:
		return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	cursor = transhelve.cursor()
	if metagene_tranlist == "":
		#print"selecting principal transcripts"
		cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1")
	else:
		metagene_tranlist = metagene_tranlist.split(",")
		cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE transcript IN ({})".format(str(metagene_tranlist).strip("[]")))

	result = cursor.fetchall()
	mgc = {"fiveprime":{},"threeprime":{}}
	iupac_dict =   {"R":"[AG]","Y":"[CU]","S":"[GC]","W":"[AU]","K":"[GU]",
					"M":"[AC]","B":"[CGU]","D":"[AGU]","D":"[AGU]",
					"H":"[ACU]","V":"[ACG]","N":"[AUGC]"}

	subseq_list = []
	for subseq in custom_seq_list.split(","):
		subseq = subseq.upper()
		subseq = subseq.replace("T","U").replace(" ","")
		for code in iupac_dict:
			subseq = subseq.replace(code,iupac_dict[code])
		subseq_list.append(subseq)
	#print "result", result
	for row in result:
		#print "row", row
		tran = row[0]
		#print "tran", tran
		if custom_search_region != "whole_gene":
			try:
				cds_start = int(row[1])
				cds_stop = int(row[2])
			except:
				continue
		seq = row[3].replace("T","U")
		if custom_search_region == "whole_gene":
			min_pos = 0
			max_pos = len(seq)
		if custom_search_region == "five_leader":
			min_pos = 0
			max_pos = cds_start
		if custom_search_region == "cds":
			min_pos = cds_start
			max_pos = cds_stop+3
		if custom_search_region == "three_trailer":
			min_pos = cds_stop
			max_pos = len(seq)
		if include_first == True:
			max_pos = min_pos+include_first_val
		if include_last == True:
			min_pos = max_pos-include_last_val
		if exclude_first == True:
			min_pos = min_pos + exclude_first_val
		if exclude_last == True:
			max_pos = max_pos-exclude_last_val

		if cds_start != "NULL" and cds_start != None:
			seq_positions = []
			for subseq in subseq_list:
				#print "subseq",subseq
				pattern = re.compile(r"{}".format(subseq))
				m = ""
				search_pos = min_pos-1
				
				while m != None:
					#print "search pos, max_pos", search_pos, max_pos, seq[search_pos:max_pos]
					m = pattern.search(seq,search_pos,max_pos)
					if m != None:
						position = m.span()[0]
						seq_positions.append(position)
						search_pos = position+1
						#print"cds_start, search_pos",cds_start, search_pos
		if seq_positions != []:
			#print "seq positions", seq_positions
			profile = {"fiveprime":{},"threeprime":{}}
			try:
				tran_reads = sqlite_db[tran]["unambig"]
			except:
				tran_reads = {"unambig":{}}
			#print "tran reads", tran_reads[25]
			
			for readlen in tran_reads:
				if readlen in offsets["fiveprime"]["offsets"]:
					offset = offsets["fiveprime"]["offsets"][readlen]
				else:
					offset = 15
				#print "offset", offset
				profile["fiveprime"][readlen] = {}
				profile["threeprime"][readlen] = {}
				for pos in tran_reads[readlen]:
					relative_frame = ((pos+offset)-cds_start)%3
					#print "pos, relative frame", pos, relative_frame
					#print "pos, cds_start, relative_frame", pos, cds_start, relative_frame
					#if pos != 0:
					#	relative_frame = seq_position%pos
					#else:
					#	relative_frame = 0
					#print "relative frame", relative_frame
					#print "metagene frame", metagene_frame
					if metagene_frame == "minus_frame":
						if relative_frame != 1:
							continue
					if metagene_frame == "in_frame":
						if relative_frame != 2:
							#print "its not 0", pos%3
							continue
					if metagene_frame == "plus_frame":
						if relative_frame != 0:
							continue
					#print "including"
					three_pos = pos+readlen
					count = tran_reads[readlen][pos]
					try:
						profile["fiveprime"][readlen][pos] += count
					except:
						profile["fiveprime"][readlen][pos] = 0
						profile["fiveprime"][readlen][pos] += count
					try:
						profile["threeprime"][readlen][three_pos] += count
					except:
						profile["threeprime"][readlen][three_pos] = 0
						profile["threeprime"][readlen][three_pos] += count
			for readlen in profile["fiveprime"]:
				if readlen not in mgc["fiveprime"]:
					mgc["fiveprime"][readlen] = {}
					for i in range(-600,601):
						mgc["fiveprime"][readlen][i] = 0
				for pos in profile["fiveprime"][readlen]:
					#print "pos", pos
					count = profile["fiveprime"][readlen][pos]
					for seq_position in seq_positions:
						#print "seq position ", seq_position
						#if pos != 0:
						#	relative_frame = seq_position%pos
						#else:
						#	relative_frame = 0
						#print "relative frame", relative_frame
						#print "metagene frame", metagene_frame
						#if metagene_frame == "minus_frame":
						#	if relative_frame != 1:
						#		continue
						#if metagene_frame == "in_frame":
						#	if relative_frame != 0:
						#		#print "its not 0", pos%3
						#		continue
						#if metagene_frame == "plus_frame":
						#	if relative_frame != 2:
						#		continue
						#print pos%3
						if seq_position >= pos-600 and seq_position <= pos+600:
							relative_seq_position = pos-seq_position
							#if metagene_frame == "minus_frame":
							#	if relative_seq_position%3 != 1:
							#		count = 0
							#if metagene_frame == "in_frame":
							#	if relative_seq_position%3 != 0:
							#		count = 0
							#if metagene_frame == "plus_frame":
							#	if relative_seq_position%3 != 2:
							#		count = 0
							#print "rsp, count", relative_seq_position%3, count
							mgc["fiveprime"][readlen][relative_seq_position] += count

			for readlen in profile["threeprime"]:
				if readlen not in mgc["threeprime"]:
					mgc["threeprime"][readlen] = {}
					for i in range(-600,601):
						mgc["threeprime"][readlen][i] = 0
				for pos in profile["threeprime"][readlen]:
					count = profile["threeprime"][readlen][pos]
					for seq_position in seq_positions:
						if seq_position >= pos-600 and seq_position <= pos+600:
							relative_seq_position = pos-seq_position
							mgc["threeprime"][readlen][relative_seq_position] += count
		else:
			print "No valid sequence positions"
	sqlite_db[custom_metagene_id] = mgc
	sqlite_db.commit()
	#print"returning mgc"
	return mgc







metainfoquery_blueprint = Blueprint("metainfoquery", __name__, template_folder="templates")
@metainfoquery_blueprint.route('/metainfoquery', methods=['POST'])
def metainfoquery():
	global user_short_passed
	tran_dict = {}
	gene_dict = {}
	data = ast.literal_eval(request.data)
	plottype = data["plottype"]
	minreadlen = int(data['minreadlen'])
	maxreadlen = int(data['maxreadlen'])
	
	sw_diff_group1 = (data['sw_diff_group1'].strip(",")).split(",")
	sw_diff_group2 = (data['sw_diff_group2'].strip(",")).split(",")
	sw_diff_window_size = int(data['sw_diff_window_size'])
	sw_diff_step_size = int(data['sw_diff_step_size'])
	sw_diff_min_diff = int(data['sw_diff_min_diff'])
	
	
	single_tran_de_transcript = data['single_tran_de_transcript']

	tran_corr_transcript1 = data['tran_corr_transcript1']
	tran_corr_transcript2 = data['tran_corr_transcript2']

	single_tran_de_range1 = data['single_tran_de_range1']
	custom_seq_list = data["custom_seq_list"]
	exclude_first_val = int(data["exclude_first_val"])
	exclude_last_val = int(data["exclude_last_val"])
	include_first_val = int(data["include_first_val"])
	include_last_val = int(data["include_last_val"])
	metagene_tranlist = (data["metagene_tranlist"].strip(" ")).replace(" ",",")
	trip_minreadlen = int(data['trip_minreadlen'])
	trip_maxreadlen = int(data['trip_maxreadlen'])
	mismatch_minreadcount = int(data['mismatch_minreadcount'])
	mismatch_minper = int(data['mismatch_minper'])
	mismatch_maxper = int(data['mismatch_maxper'])
	mismatch_maxhit = int(data['mismatch_maxhit'])
	mismatch_region = data['mismatch_region']
	nuc_minreadlen = int(data['nuc_minreadlen'])
	nuc_maxreadlen = int(data['nuc_maxreadlen'])
	mismatch_minreadlen = int(data['mismatch_minreadlen'])
	mismatch_maxreadlen = int(data['mismatch_maxreadlen'])
	heatmap_minreadlen = int(data['heatmap_minreadlen'])
	heatmap_maxreadlen = int(data['heatmap_maxreadlen'])
	heatmap_startpos = int(data['heatmap_startpos'])
	heatmap_endpos = int(data['heatmap_endpos'])
	maxscaleval = data['maxscaleval']
	
	if maxscaleval != "None" and maxscaleval != "":
		try:
			maxscaleval = int(maxscaleval)
		except:
			maxscalval = "None"
	color_palette = (data["color_palette"])
	minimum_reads = int(data['minimum_reads'])
	te_minimum_reads = data['te_minimum_reads']
	raw_te_tranlist = data['te_tranlist']
	raw_te_tranlist = raw_te_tranlist.replace(" ",",")
	te_tranlist = []
	for item in raw_te_tranlist.split(","):
		te_tranlist.append(item)
	organism = data['organism']
	transcriptome = data['transcriptome']
	metagene_type = data["metagene_type"]
	metagene_end = data["metagene_end"]
	custom_search_region = data["custom_search_region"]
	metagene_frame = data["metagene_frame"]
	heatmap_metagene_type = data["heatmap_metagene_type"]
	contaminant_organism = data["contaminant_organism"]
	nuc_comp_type = data["nuc_comp_type"]
	nuc_comp_frame = data["nuc_comp_frame"]
	nuc_comp_direction = data["nuc_comp_direction"]
	heatmap_direction = data["heatmap_direction"]
	html_args = data["html_args"]
	all_seq_types = data["all_seq_types"]
	nuccomp_reads = data["nuccomp_reads"]
	corr_type = data["corr_type"]
	count_type = data["count_type"]

	try:
		user = current_user.name

	except:
		user = None
	connection = sqlite3.connect('/home/DATA/www/tripsviz/tripsviz/trips.sqlite')
	connection.text_factory = str
	cursor = connection.cursor()

	background_col = config.BACKGROUND_COL
	readlength_col = config.READLENGTH_COL
	metagene_fiveprime_col = config.METAGENE_FIVEPRIME_COL
	metagene_threeprime_col = config.METAGENE_THREEPRIME_COL
	a_col = config.A_COL
	t_col = config.T_COL
	g_col = config.G_COL
	c_col = config.C_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE

	#get a list of organism id's this user can access
	if user != None:
		#get user_id
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		cursor.execute("SELECT background_col,readlength_col,metagene_fiveprime_col,metagene_threeprime_col,nuc_comp_a_col,nuc_comp_t_col,nuc_comp_g_col,nuc_comp_c_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		readlength_col = result[1]
		metagene_fiveprime_col = result[2]
		metagene_threeprime_col = result[3]
		a_col = result[4]
		t_col = result[5]
		g_col = result[6]
		c_col = result[7]
		title_size = result[8]
		subheading_size = result[9]
		axis_label_size = result[10]
		marker_size = result[11]

	if "readlen_ambig" in data:
		readlen_ambig = True
	else:
		readlen_ambig = False

	if "breakdown_per" in data:
		breakdown_per = True
	else:
		breakdown_per = False

	if "metagene_offsets" in data:
		metagene_offsets = True
	else:
		metagene_offsets = False

	if "metagene_aggregate" in data:
		metagene_aggregate = True
	else:
		metagene_aggregate = False

	if "metagene_normalise" in data:
		metagene_normalise = True
	else:
		metagene_normalise = False

	if "exclude_first" in data:
		exclude_first = True
	else:
		exclude_first = False

	if "exclude_last" in data:
		exclude_last = True
	else:
		exclude_last = False

	if "include_first" in data:
		include_first = True
	else:
		include_first = False

	if "include_last" in data:
		include_last = True
	else:
		include_last = False
	if "normalise" in data:
		normalise = True
	else:
		normalise = False

	total_files = len(data["file_list"])
	# Send file_list (a list of integers intentionally encoded as strings due to javascript), to be converted to a dictionary with riboseq/rnaseq lists of file paths.
	file_paths_dict = fetch_file_paths(data["file_list"],organism)

	if "log_scale" in data:
		log_scale = True
	else:
		log_scale = False

	if "mrna_dist_per" in data:
		mrna_dist_per = True
	else:
		mrna_dist_per = False

	if "md_start" in data:
		md_start = True
	else:
		md_start = False

	if "md_stop" in data:
		md_stop = True
	else:
		md_stop = False


	if "mrna_readlen_per" in data:
		mrna_readlen_per = True
	else:
		mrna_readlen_per = False
	if "count_agg" in data:
		count_agg = True
	else:
		count_agg = False
	if "apply_offset" in data:
		apply_offset = True
	else:
		apply_offset = False
	if "mismatch_agg" in data:
		mismatch_agg = True
	else:
		mismatch_agg = False

	if "reverse_scale" in data:
		reverse_scale = True
	else:
		reverse_scale = False

	if html_args["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"metainfo_plot")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True
	if plottype == "sw_diff":
		curr_time = time.time()
		filename = "SW_Diff_{}".format(curr_time)
		outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		outfile.write("Gene,Tran,Window_start,Grp1_count,Grp2_count,Ratio,Difference\n")
		trandict = {}
		traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{2}.sqlite".format(organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		codon_dict = {}
		principal_transcripts = {}
		traninfo_cursor.execute("SELECT transcript,cds_start,cds_stop,length,gene FROM transcripts WHERE principal = 1;")
		result = traninfo_cursor.fetchall()
		gene_dict = {}
		for row in result:
			trandict[row[0]] = {"cds_start":row[1],"cds_stop":row[2],"length":row[3]}
			gene_dict[row[0]] = row[4]
					
					
		grp1_profile_dict = {}
		grp2_profile_dict = {}
		grp1_file_paths_dict = fetch_file_paths(sw_diff_group1,organism)
		#printgrp1_file_paths_dict
		grp2_file_paths_dict = fetch_file_paths(sw_diff_group2,organism)
		for file_id in sw_diff_group1:
			#printgrp1_file_paths_dict["riboseq"]
			sqlite_db = SqliteDict(grp1_file_paths_dict["riboseq"][int(file_id)])
			offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
			#print"offsets", offsets
			for transcript in trandict.keys():
				if transcript not in grp1_profile_dict:
					grp1_profile_dict[transcript] = {}
				try:
					counts = sqlite_db[transcript]
				except:
					continue
				subprofile = build_profile(counts, offsets,False)
				for pos in subprofile:
					try:
						grp1_profile_dict[transcript][pos] += subprofile[pos]
					except:
						grp1_profile_dict[transcript][pos] = subprofile[pos]
			sqlite_db.close()
			
			
			
		for file_id in sw_diff_group2:
			
			sqlite_db = SqliteDict(grp2_file_paths_dict["riboseq"][int(file_id)])
			offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
			#print"offsets", offsets
			for transcript in trandict.keys():
				if transcript not in grp2_profile_dict:
					grp2_profile_dict[transcript] = {}
				try:
					counts = sqlite_db[transcript]
				except:
					continue
				subprofile = build_profile(counts, offsets,False)
				for pos in subprofile:
					try:
						grp2_profile_dict[transcript][pos] += subprofile[pos]
					except:
						grp2_profile_dict[transcript][pos] = subprofile[pos]
			sqlite_db.close()
		for tran in grp1_profile_dict:
			#if tran != "ENSMUST00000018470":
			#	continue
			#print grp1_profile_dict[tran]
			#print "\n\n\n"
			#print grp2_profile_dict[tran]
			if tran in grp2_profile_dict:
				tranlen = trandict[tran]["length"]
				startpoint = 0
				endpoint = tranlen
				if custom_search_region != "whole_gene":
					try:
						cds_start = int(trandict[tran]["cds_start"])
						cds_stop = int(trandict[tran]["cds_stop"])
					except:
						continue
					if custom_search_region == "five_leader":
						startpoint = 0
						endpoint = cds_start
					elif custom_search_region == "cds":
						startpoint = cds_start
						endpoint = cds_stop
					elif custom_search_region == "three_trailer":
						startpoint = cds_stop
						endpoint = tranlen
				for i in range(startpoint,endpoint, sw_diff_step_size):
					#if i != 690:
					#	continue
					grp1_count = 1.0
					grp2_count = 1.0
					for x in range(i,i+sw_diff_window_size):
						#print "x",x
						if x in grp1_profile_dict[tran]:
							#print "adding for grp1",  grp1_profile_dict[tran][x]
							grp1_count += grp1_profile_dict[tran][x]
						if x in grp2_profile_dict[tran]:
							#print "adding fro grp2", grp2_profile_dict[tran][x]
							grp2_count += grp2_profile_dict[tran][x]
					diff = abs(grp1_count-grp2_count)
					if diff > sw_diff_min_diff:
						#print "grp1 count grp2 coutn, i, "
						ratio = grp1_count/grp2_count
						gene = gene_dict[tran]
						outfile.write("{},{},{},{},{},{},{}\n".format(gene,tran,i,grp1_count,grp2_count,ratio,(grp1_count-grp2_count)))

		
		return filename

	if plottype == "bulk_dl":
		curr_time = time.time()
		transcripts = data["region"]
		traninfo_dict = {}
		traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{2}.sqlite".format(organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		if transcripts == "custom":
			traninfo_cursor.execute("SELECT transcript,gene,cds_start,cds_stop,length FROM transcripts WHERE transcript IN ({})".format(str(te_tranlist).strip("[]").replace('"','')))
		elif transcripts == "principal":
			traninfo_cursor.execute("SELECT transcript,gene,cds_start,cds_stop,length FROM transcripts WHERE principal = 1")
		elif transcripts == "all":
			traninfo_cursor.execute("SELECT transcript,gene,cds_start,cds_stop,length FROM transcripts")
		result = traninfo_cursor.fetchall()
		for row in result:
			#print "length", row[4]
			traninfo_dict[row[0]] = {"gene":row[1], "cds_start":row[2], "cds_stop":row[3], "length":int(row[4])}
		identifiers = []
		filepaths = []
		map_types = ["unambig"]
		if readlen_ambig == True:
			map_types.append("ambig")
		for transcript in traninfo_dict:
			outfile = open("{}/static/tmp/{}_{}".format(config.SCRIPT_LOC,transcript,curr_time),"w")
			filepaths.append("{}_{}".format(transcript,curr_time))
			count_dict = {}
			for filetype in file_paths_dict:
				for file_id in file_paths_dict[filetype]:
					filepath = file_paths_dict[filetype][file_id]
					cursor.execute("SELECT file_name, file_description FROM files WHERE file_id = {}".format(file_id))
					result = cursor.fetchone()
					file_name = result[0]
					file_desc = result[1]
					identifier = "{}_({})".format(file_name, file_desc)
					identifiers.append(identifier)
					count_dict[identifier] = {}
					for i in range(0,traninfo_dict[transcript]["length"]):
						count_dict[identifier][i] = 0
					if os.path.isfile(filepath):
						sqlite_db = SqliteDict(filepath, autocommit=False)
						offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					else:
						return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
					if transcript in sqlite_db:
						for map_type in map_types:
							counts = sqlite_db[transcript][map_type]
							for readlen in counts:
								if readlen > minreadlen and readlen < maxreadlen:
									if readlen in offsets:
										offset = offsets[readlen]
									else:
										offset = 15
									for pos in counts[readlen]:
										if apply_offset == True:
											try:
												count_dict[identifier][pos+offset] += counts[readlen][pos]
											except:
												pass
										else:
											count_dict[identifier][pos] += counts[readlen][pos]
			if traninfo_dict[transcript]["cds_start"] != None:
				outfile.write("#{}_{}_{}_{}\n".format(transcript, traninfo_dict[transcript]["gene"],traninfo_dict[transcript]["cds_start"],traninfo_dict[transcript]["cds_stop"]))
			else:
				outfile.write("#{}_{}_noncoding\n".format(transcript, traninfo_dict[transcript]["gene"]))
			outfile.write("Pos,")
			for identifier in identifiers:
				outfile.write("{},".format(identifier))
			outfile.write("Total\n")
			for i in range(0,traninfo_dict[transcript]["length"]):
				outfile.write(str(i)+",")
				pos_total = 0 
				for identifier in identifiers:
					outfile.write("{},".format(count_dict[identifier][i]))
					pos_total += count_dict[identifier][i]
				outfile.write("{}\n".format(pos_total))
			
		#print "tar -czvf {1}/static/tmp/bulk_dl_{0}.tar.gz -C  {2}".format(curr_time,config.SCRIPT_LOC, str(filepaths).strip("[]").replace(",","").replace("'",""))
		subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/bulk_dl_{0}.tar.gz {2}".format(curr_time,config.SCRIPT_LOC, str(filepaths).strip("[]").replace(",","").replace("'","")),shell=True)
				
		return  "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/static/tmp/bulk_dl_{1}.tar.gz' target='_blank' ><button class='button centerbutton' type='submit'><b>Download result</b></button></a> </div>".format(config.SCRIPT_LOC, curr_time)
		
	if plottype == "readlen_dist":
		#print"readlength dist called, custom_search_region is ", custom_search_region
		#print"metagene_tranlist is", metagene_tranlist
		master_dict = {}
		if metagene_tranlist != "":
			metagene_tranlist = metagene_tranlist.split(",")
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
				# If no transcripts given and no region specified, get the precomputed read lengths (all transcripts, entire gene)
				if metagene_tranlist == "" and custom_search_region == "whole_gene":
					if readlen_ambig == True:
						if "read_lengths" not in sqlite_db:
							return "No readlength distribution data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
						else:
							read_lengths = sqlite_db["read_lengths"]
						sqlite_db.close()
						for i in read_lengths:
							if i in master_dict:
								master_dict[i] += read_lengths[i]
							else:
								master_dict[i] = read_lengths[i]
					elif readlen_ambig == False:
						if "unambig_read_lengths" not in sqlite_db:
							return "No unambiguous readlength distribution data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
						else:
							read_lengths = sqlite_db["unambig_read_lengths"]
						sqlite_db.close()
						for i in read_lengths:
							if i in master_dict:
								master_dict[i] += read_lengths[i]
							else:
								master_dict[i] = read_lengths[i]
				else:
					traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{2}.sqlite".format(organism,transcriptome))
					traninfo_cursor = traninfo_connection.cursor()
					if metagene_tranlist == "":
						#print"metagene tranlist is", metagene_tranlist
						traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE principal = 1;")
						result = traninfo_cursor.fetchall()
					else:
						#print"metagene tranlist is", metagene_tranlist

						#print "SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE transcript IN ({})".format(str(metagene_tranlist).strip("[]"))
						traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE transcript IN ({})".format(str(metagene_tranlist).strip("[]").replace('"','')))
						result = traninfo_cursor.fetchall()
						#print "result", result
						for row in result:
							tran = row[0]
							seq = row[1]
							cds_start = row[2]
							cds_stop = row[3]
							if tran in sqlite_db:
								counts = sqlite_db[tran]["unambig"]
								for readlen in counts:
									if readlen not in master_dict:
										master_dict[readlen] = 0
									for pos in counts[readlen]:
										cnt = counts[readlen][pos]
										if custom_search_region == "whole_gene":
											master_dict[readlen] += cnt
										elif custom_search_region == "five_leader":
											if pos < cds_start:
												master_dict[readlen] += cnt
										elif custom_search_region == "cds":
											if pos > cds_start and pos < cds_stop:
												master_dict[readlen] += cnt
										elif custom_search_region == "three_trailer":
											if pos > cds_stop:
												master_dict[readlen] += cnt

		title = "Readlength distribution"
		connection.close()
		return metainfo_plots.readlen_dist(master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size)
	if plottype == "mismatch_pos":
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)

				if "global_mismatches" not in sqlite_db:
					return "No mismatch data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				else:
					mismatches = sqlite_db["global_mismatches"]
				sqlite_db.close()
				for readlen in mismatches:
					if readlen < mismatch_minreadlen or readlen > mismatch_maxreadlen:
						continue
					for pos in mismatches[readlen]:
						if pos in master_dict:
							master_dict[pos] += mismatches[readlen][pos]
						else:
							master_dict[pos] = mismatches[readlen][pos]
		title = "Mismatch positions"
		return metainfo_plots.mismatch_pos(master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "single_tran_de":
		range1 = single_tran_de_range1.split("_")
		range1_kbp = (float(int(range1[1]) - int(range1[0])))/1000
		master_list = []
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				range1_count = 1.001
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				cursor.execute("SELECT file_description from files WHERE file_id = {};".format(file_id))
				result = (cursor.fetchone())
				file_desc = result[0]
				study = filepath.split("/")[-2]
				if study not in master_dict:
					master_dict[study] = {"range1_count":0}
				if os.path.isfile(filepath):
					#Add the counts to the profile
					sqlite_db = SqliteDict(filepath, autocommit=False)
					mapped_reads = 0
					try:
						mapped_reads += float(sqlite_db["noncoding_counts"])
					except:
						pass
					try:				
						mapped_reads += float(sqlite_db["coding_counts"])
					except:
						pass
					if mapped_reads < 100000:
						return "ERROR: File {} has less than 100,000 mapped reads".format(filename)
					if "offsets" in sqlite_db:
						offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					else:
						offsets = {}
					profile = {}
					if single_tran_de_transcript in sqlite_db:
						sqlite_db_tran = sqlite_db[single_tran_de_transcript]
						for readlen in sqlite_db_tran["unambig"]:
							if readlen in offsets:
								offset = offsets[readlen]
							else:
								offset = 15
							for pos in sqlite_db_tran["unambig"][readlen]:
								count = sqlite_db_tran["unambig"][readlen][pos]
								offset_pos = offset+pos
								if offset_pos not in profile:
									profile[offset_pos] = 0
								profile[offset_pos] += count
					for x in range(int(range1[0]), int(range1[1])):
						if x in profile:
							range1_count += profile[x]
					#print"range1_count, range1_len, mapped reads",range1_count, range1_kbp, mapped_reads
					range1_tpm = (range1_count/range1_kbp)
					range1_tpm = range1_tpm/(mapped_reads/1000000)
					#print"range1_tpm", range1_tpm
				master_dict[study]["range1_count"] += range1_count
				master_list.append((file_id, filename, range1_tpm,mapped_reads,file_desc))
		study_master_list = []
		for study in master_dict:
			range1_count = master_dict[study]["range1_count"]
			study_master_list.append((0, study, range1_count))
		sorted_master_list = sorted(master_list, key=lambda x:x[1])
		return metainfo_plots.single_tran_de(single_tran_de_transcript, sorted_master_list,study_master_list,organism, transcriptome)
	elif plottype == "codon_usage":
		traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{2}.sqlite".format(organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		codon_dict = {}
		principal_transcripts = {}
		traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE principal = 1;")
		result = traninfo_cursor.fetchall()
		for row in result:
			if row[2] != "None" and row[2] != "" and row[2] != None:
				principal_transcripts[str(row[0])] = {"seq":str(row[1]),"cds_start":int(row[2]),"cds_stop":int(row[3])}

		if file_paths_dict["riboseq"] == {} and file_paths_dict["rnaseq"] == {}:
			flash("Error no files selected")
			return "Error no files selected"
		all_values = []
		offset_dict = {}
		for file_id in file_paths_dict["riboseq"]:
			sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
			try:
				offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
				offset_dict[file_id] = offsets
			except:
				offset_dict[file_id] = {}
			sqlite_db.close()
		tran_count = 0


		for file_id in file_paths_dict["riboseq"]:
			sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
			if "codon_usage_dict2" in sqlite_db:
				codon_usage_dict = sqlite_db["codon_usage_dict"]
				for codon in codon_usage_dict:
					if codon not in codon_dict:
						codon_dict[codon] = {"ribo_count":0,"codon_count":0.0}
					codon_dict[codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
					codon_dict[codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
			else:
				#codon_dict is the main dict that holds counts from all files, codon_usage_dict is file specific and will be saved for quick access later.
				codon_usage_dict = {}
				offsets = offset_dict[file_id]
				#print"old offsets", offsets
				poffsets = {}
				for offset in offsets:
					value = offsets[offset]
					new_value = value - 3
					poffsets[offset] = new_value
				#print"new offsets", poffsets
				for transcript in principal_transcripts:
					tran_count += 1
					if tran_count%1000 == 0:
						print"tran_count", tran_count
					profile = {}
					if transcript not in sqlite_db:
						continue
					
					subprofile = build_profile(sqlite_db[transcript], poffsets,"unambig")
					for pos in subprofile:
						if pos not in profile:
							profile[pos] = 0
						profile[pos] += subprofile[pos]
					seq = principal_transcripts[transcript]["seq"]
					for i in range(principal_transcripts[transcript]["cds_start"]-1,principal_transcripts[transcript]["cds_start"]+30):
						codon = seq[i:i+3]
						if len(codon) != 3:
							continue
						count = 0
						if i in profile:
							count += profile[i]
						#if i+1 in profile:
						#	count += profile[i+1]
						#if i+2 in profile:
						#	count += profile[i+2]
						if codon not in codon_dict:
							codon_dict[codon] = {"ribo_count":0,"codon_count":0.0}
						if codon not in codon_usage_dict:
							codon_usage_dict[codon] = {"ribo_count":0,"codon_count":0.0}
						codon_usage_dict[codon]["ribo_count"] += count
						codon_usage_dict[codon]["codon_count"] += 1
				for codon in codon_usage_dict:
					codon_dict[codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
					codon_dict[codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
				sqlite_db["codon_usage_dict"] = codon_usage_dict
				sqlite_db.commit()
			sqlite_db.close()
		return metainfo_plots.codon_usage(codon_dict,short_code,str(title_size)+"pt", str(axis_label_size)+"pt", str(marker_size)+"pt")
	elif plottype == "diff_codon_usage":
		traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.{2}.sqlite".format(organism,transcriptome))
		traninfo_cursor = traninfo_connection.cursor()
		codon_dict_cond = {"condition1":{},"condition2":{}}
		diff_codon_dict = {}
		principal_transcripts = {}
		condition_totals = {"condition1":0,"condition2":0}
		traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop FROM transcripts WHERE principal = 1;")
		result = traninfo_cursor.fetchall()
		for row in result:
			if row[2] != "None" and row[2] != "" and row[2] != None:
				principal_transcripts[str(row[0])] = {"seq":str(row[1]),"cds_start":int(row[2]),"cds_stop":int(row[3])}

		if file_paths_dict["riboseq"] == {} and file_paths_dict["rnaseq"] == {}:
			flash("Error no files selected")
			return "Error no files selected"

		condition_dict = {"condition1":[file_paths_dict["riboseq"].keys()[0]], "condition2":[file_paths_dict["riboseq"].keys()[1]]}

		all_values = []
		offset_dict = {}
		for condition in condition_dict:
			for file_id in condition_dict[condition]:
				sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
				try:
					offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					offset_dict[file_id] = offsets
				except:
					offset_dict[file_id] = {}
				sqlite_db.close()
			tran_count = 0

			for file_id in condition_dict[condition]:
				sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
				if "codon_usage_dict" in sqlite_db:
					codon_usage_dict = sqlite_db["codon_usage_dict"]
					for codon in codon_usage_dict:
						if codon not in codon_dict_cond[condition]:
							codon_dict_cond[condition][codon] = {"ribo_count":0,"codon_count":0.0}
						codon_dict_cond[condition][codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
						condition_totals[condition] += codon_usage_dict[codon]["ribo_count"]
						codon_dict_cond[condition][codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
				else:
					#codon_dict_cond is the main dict that holds counts from all files, codon_usage_dict is file specific and will be saved for quick access later.
					codon_usage_dict = {}
					for transcript in principal_transcripts:
						tran_count += 1
						if tran_count%1000 == 0:
							print "tran_count", tran_count
						profile = {}
						if transcript not in sqlite_db:
							continue
						offsets = offset_dict[file_id]
						subprofile = build_profile(sqlite_db[transcript], offsets,"unambig")
						for pos in subprofile:
							if pos not in profile:
								profile[pos] = 0
							profile[pos] += subprofile[pos]
						seq = principal_transcripts[transcript]["seq"]
						for i in range(principal_transcripts[transcript]["cds_start"], principal_transcripts[transcript]["cds_stop"],3):
							codon = seq[i:i+3]
							if len(codon) != 3:
								continue
							count = 0
							if i in profile:
								count += profile[i]
							if i+1 in profile:
								count += profile[i+1]
							if i+2 in profile:
								count += profile[i+2]
							if codon not in codon_dict_cond[condition]:
								codon_dict_cond[condition][codon] = {"ribo_count":0,"codon_count":0.0}
							#codon_dict_cond[codon]["ribo_count"] += count
							#codon_dict_cond[codon]["codon_count"] += 1

							if codon not in codon_usage_dict:
								codon_usage_dict[codon] = {"ribo_count":0,"codon_count":0.0}
							codon_usage_dict[codon]["ribo_count"] += count
							codon_usage_dict[codon]["codon_count"] += 1
					for codon in codon_usage_dict:
						codon_dict_cond[condition][codon]["ribo_count"] += codon_usage_dict[codon]["ribo_count"]
						condition_totals[condition] += codon_usage_dict[codon]["ribo_count"]
						codon_dict_cond[condition][codon]["codon_count"] = codon_usage_dict[codon]["codon_count"]
					sqlite_db["codon_usage_dict"] = codon_usage_dict
					sqlite_db.commit()
				sqlite_db.close()
		factor_diff = float(condition_totals["condition1"])/float(condition_totals["condition2"])
		for codon in codon_dict_cond["condition1"]:
			count1 = codon_dict_cond["condition1"][codon]["ribo_count"]
			count2 = codon_dict_cond["condition2"][codon]["ribo_count"]*factor_diff
			diff = count1-count2
			diff_codon_dict[codon] = {"ribo_count":diff, "codon_count":codon_dict_cond["condition1"][codon]["codon_count"]}
		return metainfo_plots.codon_usage(diff_codon_dict,short_code,str(title_size)+"pt", str(axis_label_size)+"pt", str(marker_size)+"pt")
	elif plottype == "tran_corr":
		master_list = []
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				tran1_count = 1.001
				tran2_count = 1.001
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				study = filepath.split("/")[-2]
				if filename not in master_dict:
					master_dict[filename] = {"tran1_count":0,
											"tran2_count":0}
				if os.path.isfile(filepath):
					#Add the counts to the profile
					sqlite_db = SqliteDict(filepath, autocommit=False)
					if "offsets" in sqlite_db:
						offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
					else:
						offsets = {}
					profile = {}
					#TRAN1
					if tran_corr_transcript1 in sqlite_db:
						sqlite_db_tran = sqlite_db[tran_corr_transcript1]
						for readlen in sqlite_db_tran["unambig"]:
							if readlen in offsets:
								offset = offsets[readlen]
							else:
								offset = 15
							for pos in sqlite_db_tran["unambig"][readlen]:
								count = sqlite_db_tran["unambig"][readlen][pos]
								offset_pos = offset+pos
								if offset_pos not in profile:
									profile[offset_pos] = 0
								profile[offset_pos] += count
					for pos in profile:
						tran1_count += profile[pos]
					#TRAN2
					profile = {}
					if tran_corr_transcript2 in sqlite_db:
						sqlite_db_tran = sqlite_db[tran_corr_transcript2]
						for readlen in sqlite_db_tran["unambig"]:
							if readlen in offsets:
								offset = offsets[readlen]
							else:
								offset = 15
							for pos in sqlite_db_tran["unambig"][readlen]:
								count = sqlite_db_tran["unambig"][readlen][pos]
								offset_pos = offset+pos
								if offset_pos not in profile:
									profile[offset_pos] = 0
								profile[offset_pos] += count
					for pos in profile:
						tran2_count += profile[pos]
				master_dict[filename]["tran1_count"] += tran1_count
				master_dict[filename]["tran2_count"] += tran2_count
				master_list.append((file_id, filename, log(tran1_count,2), log(tran2_count,2), study))
		sorted_master_list = sorted(master_list, key=lambda x:x[2])
		return metainfo_plots.tran_corr(tran_corr_transcript1, tran_corr_transcript2,sorted_master_list,organism, transcriptome)
	elif plottype == "mismatches":
		positive_hits = 0
		result_list = []
		file_string = ""
		total_trans = 0
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				file_string += "{},".format(file_id)
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(config.ANNOTATION_DIR,organism), autocommit=False)
		else:
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome), autocommit=False)
		if organism == "homo_sapiens" or organism == "homo_sapiens_polio":
			longest_tran_db = SqliteDict("{0}homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite".format(config.ANNOTATION_DIR), autocommit=False)
			longest_tran_list = longest_tran_db["transcripts"]
			longest_tran_db.close()
		else:
			longest_tran_list = traninfo_dict.keys()

		if mismatch_agg == True:
			for transcript in longest_tran_list:
				total_trans += 1
				if total_trans%100 == 0:
					print "total_trans", total_trans
				cds_start = traninfo_dict[transcript]["cds_start"]
				cds_stop = traninfo_dict[transcript]["cds_stop"]
				tranlen = traninfo_dict[transcript]["length"]
				if mismatch_region == "all":
					minpos = 0
					maxpos = tranlen
				elif mismatch_region == "fiveleader":
					minpos = 0
					maxpos = cds_start
				elif mismatch_region == "cds":
					minpos = cds_start
					maxpos = cds_stop
				elif mismatch_region == "threetrailer":
					minpos = cds_stop
					maxpos = tranlen
				elif mismatch_region == "cds_start":
					minpos = cds_start-0
					maxpos = cds_start+3
				elif mismatch_region == "cds_stop":
					minpos = cds_stop-1
					maxpos = cds_stop+2
				if positive_hits > mismatch_maxhit:
					break
				sequence = traninfo_dict[transcript]["seq"]
				profile = {}
				mismatch_profile = {"A":{},"T":{},"G":{},"C":{}}
				for filetype in file_paths_dict:
					for file_id in file_paths_dict[filetype]:

						filepath = file_paths_dict[filetype][file_id]
						if os.path.isfile(filepath):

							#Add the counts to the profile
							sqlite_db = SqliteDict(filepath, autocommit=False)
							if transcript in sqlite_db:
								sqlite_db_tran = sqlite_db[transcript]
								for readlen in sqlite_db_tran["unambig"]:
									for pos in sqlite_db_tran["unambig"][readlen]:
										count = sqlite_db_tran["unambig"][readlen][pos]
										for x in range(pos, pos+readlen):
											if x not in profile:
												profile[x] = 0
											profile[x] += count

								#Add the mismatch to the profile
								sqlite_db_seqvar = dict(sqlite_db[transcript]["seq"])
								for pos in sqlite_db_seqvar:
									#convert to one based
									fixed_pos = pos+1
									for char in sqlite_db_seqvar[pos]:
										if char != "N":
											if fixed_pos not in mismatch_profile[char]:
												mismatch_profile[char][fixed_pos] = 0
											count = sqlite_db_seqvar[pos][char]
											mismatch_profile[char][fixed_pos] += count
						else:
							return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
						sqlite_db.close()

				for pos in profile:
					if pos < minpos or pos > maxpos:
						continue
					if profile[pos] > mismatch_minreadcount:
						for char in mismatch_profile:
							if pos in mismatch_profile[char]:
								mismatch_count = mismatch_profile[char][pos]
								per = (float(mismatch_count)/float(profile[pos]))*100
								if per > 100:
									per = 100
								transition = "{}->{}".format(sequence[pos-1], char)
								if per > mismatch_minper and per < mismatch_maxper:
									#return "Mismatch in transcript {} at position {}, read count:{}, mismatch count {}, transition {}".format(transcript, pos, profile[pos], mismatch_count, transition)
									trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?&hili='+str(pos-10)+'_'+str(pos+10)+'&tran='+transcript+'&cov=T&nuc=T&files='+file_string+'" target="_blank_" >View on trips-viz</a>'
									result_list.append([transcript, pos, profile[pos], mismatch_count, transition, int(per),trips_link, "Aggregate"])
									positive_hits += 1
		else:
			for transcript in longest_tran_list:
				total_trans += 1
				if total_trans%100 == 0:
					print "total_trans", total_trans
				cds_start = traninfo_dict[transcript]["cds_start"]
				cds_stop = traninfo_dict[transcript]["cds_stop"]
				tranlen = traninfo_dict[transcript]["length"]
				if mismatch_region == "all":
					minpos = 0
					maxpos = tranlen
				elif mismatch_region == "fiveleader":
					minpos = 0
					maxpos = cds_start
				elif mismatch_region == "cds":
					minpos = cds_start
					maxpos = cds_stop
				elif mismatch_region == "threetrailer":
					minpos = cds_stop
					maxpos = tranlen
				elif mismatch_region == "cds_start":
					minpos = cds_start-0
					maxpos = cds_start+2
				elif mismatch_region == "cds_stop":
					minpos = cds_stop-0
					maxpos = cds_stop+2


				if positive_hits > mismatch_maxhit:
					break
				sequence = traninfo_dict[transcript]["seq"]
				file_string = ""
				for filetype in file_paths_dict:
					for file_id in file_paths_dict[filetype]:
						file_string += "{},".format(file_id)
						filepath = file_paths_dict[filetype][file_id]
						if os.path.isfile(filepath):

							#Add the counts to the profile
							sqlite_db = SqliteDict(filepath, autocommit=False)
							profile = {}
							mismatch_profile = {"A":{},"T":{},"G":{},"C":{}}
							if transcript in sqlite_db:
								sqlite_db_tran = sqlite_db[transcript]
								for readlen in sqlite_db_tran["unambig"]:
									for pos in sqlite_db_tran["unambig"][readlen]:
										count = sqlite_db_tran["unambig"][readlen][pos]
										for x in range(pos, pos+readlen):
											if x not in profile:
												profile[x] = 0
											profile[x] += count

								#Add the mismatch to the profile
								sqlite_db_seqvar = dict(sqlite_db[transcript]["seq"])
								for pos in sqlite_db_seqvar:
									#convert to one based
									fixed_pos = pos+1
									for char in sqlite_db_seqvar[pos]:
										if char != "N":
											if fixed_pos not in mismatch_profile[char]:
												mismatch_profile[char][fixed_pos] = 0
											count = sqlite_db_seqvar[pos][char]
											mismatch_profile[char][fixed_pos] += count
						else:
							return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page. ".format(filepath)
						sqlite_db.close()

						for pos in profile:
							if pos < minpos or pos > maxpos:
								continue
							if profile[pos] > mismatch_minreadcount:
								for char in mismatch_profile:
									if pos in mismatch_profile[char]:
										mismatch_count = mismatch_profile[char][pos]
										per = (float(mismatch_count)/float(profile[pos]))*100
										if per > 100:
											per = 100
										transition = "{}->{}".format(sequence[pos-1], char)
										if per > mismatch_minper and per < mismatch_maxper:
											#return "Mismatch in transcript {} at position {}, read count:{}, mismatch count {}, transition {}".format(transcript, pos, profile[pos], mismatch_count, transition)
											trips_link = '<a href="https://trips.ucc.ie/'+organism+'/'+transcriptome+'/interactive_plot/?&hili='+str(pos-10)+'_'+str(pos+10)+'&tran='+transcript+'&cov=T&nuc=T&files='+file_string+'" target="_blank_" >View on trips-viz</a>'
											result_list.append([transcript, pos, profile[pos], mismatch_count, transition, int(per),trips_link,(filepath.split("/"))[-1]  ])
											positive_hits += 1

		table_str = "<table class='prediction_table hover'><tr><th>Transcript</th><th>Position</th><th>Filename</th><th>Read Count</th><th>Mismatch Count</th><th>Percentage</th><th>Transition</th><th>View</th></td>"
		for item in result_list:
			table_str += "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(item[0],item[1],item[7],item[2],item[3],item[5], item[4],item[6])
		table_str += "</table>"
		return table_str
	elif plottype == "te":
		#traninfo_dict = shelve.open("{0}{1}/{1}.shelf".format(config.ANNOTATION_DIR,organism))
		region = data["region"]
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			#transhelve = SqliteDict("{0}{1}/{1}.sqlite".format(config.ANNOTATION_DIR,organism), autocommit=False)
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			#print "getting tran data", "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome)
			#transhelve = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome), autocommit=False)
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
			print "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome)
		cursor = transhelve.cursor()
		#print " this is the te tranlist :{}:".format(te_tranlist)
		
		if te_tranlist == "" or te_tranlist == ['']:
			te_tranlist = None
		if count_type != "tpm":
			if te_tranlist == None:
				if region == "all":
					cursor.execute("SELECT * from transcripts WHERE principal = 1")
				else:
					cursor.execute("SELECT * from transcripts WHERE principal = 1 and tran_type = 1")
			else:
				#print "SELECT * from transcripts WHERE transcript IN ('{}')".format(str(te_tranlist).strip("[]"))
				cursor.execute("SELECT * from transcripts WHERE transcript IN ({})".format(str(te_tranlist).strip("[]")))
		else:
			print "SELECTING EVERY TRASNSCIPRT"
			cursor.execute("SELECT * from transcripts")

		if count_type == "tpm":
			if total_files > 200:
				return "Error: A maximum of 200 files can be used when TPM is selected."

				
		if te_tranlist == None:
			if total_files >= 31 and count_agg == False:
				return "Error: A maximum of 30 files can be selected if not aggregating counts and using all transcripts. Reduce number of selected files or click the 'Aggregate counts' checkbox at the top right of the page. Alternatively input a list of transcripts in the 'Transcript list' box."

		te_minimum_reads = int(te_minimum_reads)
		transcript_list =  cursor.fetchall()
		#print transcript_list
		# User may have passed list of genes instead of transcripts
		if transcript_list == [] and te_tranlist != None:
			#print "SELECT * from transcripts WHERE gene IN ({}) AND principal = 1".format(str(te_tranlist).strip("[]"))
			cursor.execute("SELECT * from transcripts WHERE gene IN ({}) AND principal = 1".format(str(te_tranlist).strip("[]")))
			transcript_list =  cursor.fetchall()
		
		traninfo_dict = {}
		for result in transcript_list:
			if result[0] != "":
				traninfo_dict[result[0]] = {"transcript":result[0] , "gene":result[1].replace(",","_").replace(";","_"), "length":result[2] , "cds_start":result[3] , "cds_stop":result[4] , "seq":result[5].upper() ,
					"strand":result[6], "stop_list":result[7].split(","),"start_list":result[8].split(","), "exon_junctions":result[9].split(","),
					"tran_type":result[10], "principal":result[11]}

		transhelve.close()
		longest_tran_list = traninfo_dict.keys()
		print "Transcripts NUMBER", len(longest_tran_list)
		if count_agg == True:
			table_str = aggregate_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism,all_seq_types, te_minimum_reads, html_args,count_type,te_tranlist)
		elif count_agg == False:
			table_str = sample_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism,all_seq_types, te_minimum_reads, html_args,count_type,te_tranlist)
		return table_str

	elif plottype == "mrna_dist":
		longest_tran_list = []
		cds_dict = {}
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		cursor.execute("SELECT transcript,cds_start,cds_stop from transcripts where principal = 1 and tran_type = 1;")
		result = cursor.fetchall()
		for row in result:
			longest_tran_list.append(str(row[0]))
			cds_dict[str(row[0])] = {"cds_start":int(row[1]),"cds_stop":int(row[2])}
		mrna_dist_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				filename = (filepath.split("/")[-1]).replace(".sqlite","")
				mrna_dist_dict[filename] = {"5_leader":0,
											"start_codon":0,
											"cds":0,
											"stop_codon":0,
											"3_trailer":0,
											"total":0}
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "mrna_dist_dict" in sqlite_db:
					print "loading old results"
					mrna_dist_dict[filename]["5_leader"] = sqlite_db["mrna_dist_dict"]["5_leader"]
					mrna_dist_dict[filename]["start_codon"] = sqlite_db["mrna_dist_dict"]["start_codon"]
					mrna_dist_dict[filename]["cds"] = sqlite_db["mrna_dist_dict"]["cds"]
					mrna_dist_dict[filename]["stop_codon"] = sqlite_db["mrna_dist_dict"]["stop_codon"]
					mrna_dist_dict[filename]["3_trailer"] = sqlite_db["mrna_dist_dict"]["3_trailer"]
					mrna_dist_dict[filename]["total"] = float(sqlite_db["mrna_dist_dict"]["5_leader"]+sqlite_db["mrna_dist_dict"]["start_codon"]+sqlite_db["mrna_dist_dict"]["cds"]+sqlite_db["mrna_dist_dict"]["stop_codon"]+sqlite_db["mrna_dist_dict"]["3_trailer"])
				if mrna_dist_dict[filename]["total"] == 0:
					print "calculating new results"
					for transcript in longest_tran_list:
						try:
							transcript_dict = sqlite_db[transcript]["unambig"]
						except:
							continue
						try:
							cds_start = cds_dict[transcript]["cds_start"]
							cds_stop = cds_dict[transcript]["cds_stop"]
						except:
							continue
						for readlen in transcript_dict:
							for five_pos in transcript_dict[readlen]:
								three_pos = five_pos+readlen
								if three_pos <= cds_start+3:
									mrna_dist_dict[filename]["5_leader"] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
									mrna_dist_dict[filename]["start_codon"] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
									mrna_dist_dict[filename]["cds"] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
									mrna_dist_dict[filename]["stop_codon"] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_stop-9:
									mrna_dist_dict[filename]["3_trailer"] += transcript_dict[readlen][five_pos]
					sqlite_db["mrna_dist_dict"] = {"5_leader":mrna_dist_dict[filename]["5_leader"],
													"start_codon":mrna_dist_dict[filename]["start_codon"],
													"cds":mrna_dist_dict[filename]["cds"],
													"stop_codon":mrna_dist_dict[filename]["stop_codon"],
													"3_trailer":mrna_dist_dict[filename]["3_trailer"],
													"total":(mrna_dist_dict[filename]["5_leader"]+mrna_dist_dict[filename]["start_codon"]+mrna_dist_dict[filename]["cds"]+mrna_dist_dict[filename]["stop_codon"]+mrna_dist_dict[filename]["3_trailer"])}
					sqlite_db.commit()
				sqlite_db.close()
		connection.close()
		return metainfo_plots.mrna_dist(mrna_dist_dict,short_code, background_col,title_size, axis_label_size, subheading_size,marker_size,mrna_dist_per,md_start,md_stop)
	elif plottype == "mrna_dist_readlen":
		minreadlen = 15
		maxreadlen = 100
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(config.ANNOTATION_DIR,organism), autocommit=False)
		else:
			traninfo_dict = SqliteDict("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome), autocommit=False)
		longest_tran_list = []
		cds_dict = {}
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		cursor.execute("SELECT transcript,cds_start,cds_stop from transcripts where principal = 1 and tran_type = 1;")
		result = cursor.fetchall()
		for row in result:
			longest_tran_list.append(str(row[0]))
			cds_dict[str(row[0])] = {"cds_start":int(row[1]),"cds_stop":int(row[2])}
		mrna_dist_dict = {"5_leader":collections.OrderedDict(),
						  "start_codon":collections.OrderedDict(),
						  "cds":collections.OrderedDict(),
						  "stop_codon":collections.OrderedDict(),
						  "3_trailer":collections.OrderedDict(),
						  "total":collections.OrderedDict()}

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				try:
					sqlite_db = SqliteDict(filepath, autocommit=False)
					#opendict = dict(sqlite_db)
					#sqlite_db.close()
				except:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)

				if "mrna_dist_readlen_dict" in sqlite_db:
					for readlen in range(minreadlen,maxreadlen+1):
						if readlen not in mrna_dist_dict["5_leader"]:
							mrna_dist_dict['5_leader'][readlen]=0
							mrna_dist_dict['start_codon'][readlen]=0
							mrna_dist_dict['cds'][readlen]=0
							mrna_dist_dict['stop_codon'][readlen]=0
							mrna_dist_dict['3_trailer'][readlen]=0
							mrna_dist_dict['total'][readlen]=0
						mrna_dist_dict['5_leader'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["5_leader"][readlen]
						mrna_dist_dict['start_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["start_codon"][readlen]
						mrna_dist_dict['cds'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["cds"][readlen]
						mrna_dist_dict['stop_codon'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["stop_codon"][readlen]
						mrna_dist_dict['3_trailer'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["3_trailer"][readlen]
						mrna_dist_dict['total'][readlen] += sqlite_db["mrna_dist_readlen_dict"]["total"][readlen]
				else:
					file_specific_mrna_dist_dict = {"5_leader":{},
													"start_codon":{},
													"cds":{},
													"stop_codon":{},
													"3_trailer":{},
													"total":{}}
					for readlen in range(minreadlen,maxreadlen+1):
							file_specific_mrna_dist_dict['5_leader'][readlen]=0
							file_specific_mrna_dist_dict['start_codon'][readlen]=0
							file_specific_mrna_dist_dict['cds'][readlen]=0
							file_specific_mrna_dist_dict['stop_codon'][readlen]=0
							file_specific_mrna_dist_dict['3_trailer'][readlen]=0
							file_specific_mrna_dist_dict['total'][readlen]=0
					for transcript in longest_tran_list:
						try:
							transcript_dict = sqlite_db[transcript]["unambig"]
						except:
							continue
						try:
							cds_start = cds_dict[transcript]["cds_start"]
							cds_stop = cds_dict[transcript]["cds_stop"]
						except:
							continue
						for readlen in range(minreadlen,maxreadlen+1):
							if readlen not in mrna_dist_dict["5_leader"]:
								mrna_dist_dict['5_leader'][readlen]=0
								mrna_dist_dict['start_codon'][readlen]=0
								mrna_dist_dict['cds'][readlen]=0
								mrna_dist_dict['stop_codon'][readlen]=0
								mrna_dist_dict['3_trailer'][readlen]=0
								mrna_dist_dict['total'][readlen]=0
							if readlen not in transcript_dict:
								continue
							for five_pos in transcript_dict[readlen]:
								three_pos = five_pos+readlen
								if three_pos <= cds_start+3:
									mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["5_leader"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_start-4 and three_pos >= cds_start+4:
									mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["start_codon"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_start-3 and three_pos <= cds_stop-2:
									mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["cds"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos <= cds_stop-9 and three_pos >= cds_stop-1:
									mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["stop_codon"][readlen] += transcript_dict[readlen][five_pos]
								elif five_pos >= cds_stop-9:
									mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
									file_specific_mrna_dist_dict["3_trailer"][readlen] += transcript_dict[readlen][five_pos]
					sqlite_db["mrna_dist_readlen_dict"] = file_specific_mrna_dist_dict
					sqlite_db.commit()
				sqlite_db.close()


		# If mrna_readlen_per is true normalize everything over the max in it's category
		if mrna_readlen_per == True:
			# For each category find the max value, max value will default to 1 if there is no values in that category.
			max_five = max(1,float(max(mrna_dist_dict["5_leader"].values())))
			max_start = max(1,float(max(mrna_dist_dict["start_codon"].values())))
			max_cds = max(1,float(max(mrna_dist_dict["cds"].values())))
			max_stop = max(1,float(max(mrna_dist_dict["stop_codon"].values())))
			max_three = max(1,float(max(mrna_dist_dict["3_trailer"].values())))
			for readlen in mrna_dist_dict["5_leader"]:
				mrna_dist_dict["5_leader"][readlen] = (float(mrna_dist_dict["5_leader"][readlen])/max_five)*100
			for readlen in mrna_dist_dict["start_codon"]:
				mrna_dist_dict["start_codon"][readlen] = (float(mrna_dist_dict["start_codon"][readlen])/max_start)*100
			for readlen in mrna_dist_dict["cds"]:
				mrna_dist_dict["cds"][readlen] = (float(mrna_dist_dict["cds"][readlen])/max_cds)*100
			for readlen in mrna_dist_dict["stop_codon"]:
				mrna_dist_dict["stop_codon"][readlen] = (float(mrna_dist_dict["stop_codon"][readlen])/max_stop)*100
			for readlen in mrna_dist_dict["3_trailer"]:
				mrna_dist_dict["3_trailer"][readlen] = (float(mrna_dist_dict["3_trailer"][readlen])/max_three)*100
		connection.close()
		return metainfo_plots.mrna_dist_readlen(mrna_dist_dict, mrna_readlen_per,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "rust_dwell":
		traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(config.ANNOTATION_DIR,organism), autocommit=False)
		longest_tran_db = SqliteDict("/home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/principal_isoforms_5ldr3tlr_rnaseq.sqlite",autocommit=True)
		longest_tran_list = longest_tran_db["transcripts"]
		codon_count_dict = {"TTT":0, "TTC":0, "TTA":0, "TTG":0,
		"TCT":0, "TCC":0, "TCA":0, "TCG":0,
		"TAT":0, "TAC":0, "TAA":0, "TAG":0,
		"TGT":0, "TGC":0, "TGA":0, "TGG":0,
		"CTT":0, "CTC":0, "CTA":0, "CTG":0,
		"CCT":0, "CCC":0, "CCA":0, "CCG":0,
		"CAT":0, "CAC":0, "CAA":0, "CAG":0,
		"CGT":0, "CGC":0, "CGA":0, "CGG":0,
		"ATT":0, "ATC":0, "ATA":0, "ATG":0,
		"ACT":0, "ACC":0, "ACA":0, "ACG":0,
		"AAT":0, "AAC":0, "AAA":0, "AAG":0,
		"AGT":0, "AGC":0, "AGA":0, "AGG":0,
		"GTT":0, "GTC":0, "GTA":0, "GTG":0,
		"GCT":0, "GCC":0, "GCA":0, "GCG":0,
		"GAT":0, "GAC":0, "GAA":0, "GAG":0,
		"GGT":0, "GGC":0, "GGA":0, "GGG":0}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)

				#TODO CHANGE THIS SO THAT THE CODON COUNT DICT IS OUTSIDE THE FILEPATH FOR LOOP
				offsets = opendict["offsets"]["fiveprime"]["offsets"]
				position_dict = {}
				for transcript in longest_tran_list:
					if transcript in opendict:
						cds_start = traninfo_dict[transcript]["cds_start"]
						cds_stop = traninfo_dict[transcript]["cds_stop"]
						transeq = traninfo_dict[transcript]["seq"]
						for readlen in opendict[transcript]["unambig"]:
							offset = 15
							for pos in opendict[transcript]["unambig"][readlen]:
								a_site = (pos+offset)+1
								if a_site > cds_start+120 and a_site < cds_stop-60:
									codon = transeq[a_site:a_site+3]
									codon_count_dict[codon] += opendict[transcript]["unambig"][readlen][pos]
		return metainfo_plots.rust_dwell(codon_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "unmapped":
		master_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				if "frequent_unmapped_reads" not in sqlite_db:
					return "No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
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
		connection.close()
		return html_table
	elif plottype == "contamination":
		count_dict = {}
		master_sequence = ""
		contaminant_file = open("/home/DATA/www/tripsviz/tripsviz/static/contaminants/mycoplasma.fa")
		contaminant_lines = contaminant_file.read()
		contaminant_split = contaminant_lines.split(">")
		for entry in contaminant_split[1:]:
			header = entry.split("\n")[0]
			sequence = "".join(entry.split("\n")[1:])
			master_sequence += sequence
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				filename = filepath.split("/")[-1]
				if filename not in count_dict:
					count_dict[filename] = {"count":0,"coverage":[],"unique_reads":0}
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				if "frequent_unmapped_reads" not in sqlite_db:
					return "No unmapped reads data for {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath.split("/")[-1])
				#unmapped reads list is a list of tuples of length 100, first item in tuple is a sequence second is a count
				unmapped_reads_list = sqlite_db["frequent_unmapped_reads"]
				sqlite_db.close()
				for tup in unmapped_reads_list:
					read = tup[0]
					count = tup[1]
					readlen = len(read)
					for x in range(0,len(master_sequence)-readlen):
						mismatches = 0
						for y in range(0,readlen):
							if master_sequence[x+y] != read[y]:
								mismatches += 1
							if mismatches > 2:
								break
						if mismatches <= 2:
							count_dict[filename]["unique_reads"] += 1
							count_dict[filename]["count"] += count
							for i in range(x,x+readlen):
								if i not in count_dict[filename]["coverage"]:
									count_dict[filename]["coverage"].append(i)
		master_seq_len = len(master_sequence)
		for filename in count_dict:
			coverage = float(len(count_dict[filename]["coverage"]))/float(master_seq_len)
			coverage = round((coverage*100),2)
			count_dict[filename]["coverage"] = coverage
		title = "Contamination counts ({})".format(short_code)
		top_reads = (sorted(count_dict.items(), key=operator.itemgetter(1)))
		html_table = "<h1><center>{}</center></h1>".format(title)
		html_table += """<table class="unmapped_table">
		<thead><tr><th>Filename</th><th>Counts</th><th>Unique reads</th><th>Percentage coverage</th></tr></thead>"""
		for tup in top_reads[::-1]:
			html_table += ("<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>".format(tup[0], tup[1]["count"],tup[1]["unique_reads"],tup[1]["coverage"]))
		html_table += ("</table>")
		connection.close()
		return html_table
	elif plottype == "replicate_comp":
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
				transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
			else:
				return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		trancursor = transhelve.cursor()
		trancursor.execute("SELECT transcript,principal from transcripts")
		prin_tran_list = []
		result = trancursor.fetchall()
		for row in result:
			if row[1] == 1:
				prin_tran_list.append(row[0])
		minimum_reads = int(minimum_reads)
		mapped_reads_dict = {}
		factor_dict = {}
		if normalise == True:
			for filetype in file_paths_dict:
				for file_id in file_paths_dict[filetype]:
					filepath = file_paths_dict[filetype][file_id]
					sqlite_db = SqliteDict(filepath, autocommit=False)
					mapped_reads = sqlite_db["coding_counts"]
					if mapped_reads == None or mapped_reads == 0:
						return "Error: Mapped reads info missing for one or more files, cannot normalise"
					mapped_reads_dict[file_id] = float(mapped_reads)
			minval = min(mapped_reads_dict.values())
			for file_id in mapped_reads_dict:
				factor = minval/mapped_reads_dict[file_id]
				factor_dict[file_id] = factor
		if minimum_reads >0:
			min_log_val = log(minimum_reads,2)
		else:
			min_log_val = 0
		labels = []
		transcript_dict = {}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				cursor.execute("SELECT file_description from files where file_id = '{}';".format(file_id))
				result = cursor.fetchone();
				label = result[0]
				lbl_tag = 0
				while label in labels:
					lbl_tag += 1
					label = result[0] + str(lbl_tag)
				labels.append(label)
				filepath = file_paths_dict[filetype][file_id]

				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict =sqlite_db["unambiguous_all_totals"]
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				for transcript in prin_tran_list:
					if transcript not in transcript_dict:
						transcript_dict[transcript] = {}
					if transcript in opendict:
						try:
							count = opendict[transcript]
							if normalise == True:
								count = float(count)*factor_dict[file_id]
							if count >= min_log_val:
								transcript_dict[transcript][label] = log(count,2)
						except:
							pass
		connection.close()
		del_list = []
		for transcript in transcript_dict:
			if len(transcript_dict[transcript]) != len(labels):
				del_list.append(transcript)
		for tran in del_list:
			del transcript_dict[tran]
		return metainfo_plots.replicate_comp(labels, transcript_dict, min_log_val,short_code,background_col,str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt",str(marker_size)+"pt",corr_type)
	elif plottype == "nuc_comp":
		master_count_dict = {"A":collections.OrderedDict(),"T":collections.OrderedDict(),"G":collections.OrderedDict(),"C":collections.OrderedDict()}
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "nuc_counts" not in sqlite_db:
					return "No nucleotide counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				if nuc_comp_direction == "nuc_comp_five":
					if "nuc_counts" in sqlite_db:
						if nuccomp_reads in sqlite_db["nuc_counts"]:
							nuc_counts = sqlite_db["nuc_counts"][nuccomp_reads]
						else:
							nuc_counts = get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome)
					else:
						nuc_counts = get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome)
				elif nuc_comp_direction == "nuc_comp_three":
					if "threeprime_nuc_counts" in sqlite_db:
						if nuccomp_reads in sqlite_db["threeprime_nuc_counts"]:
							nuc_counts = sqlite_db["threeprime_nuc_counts"][nuccomp_reads]
						else:
							nuc_counts = get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome)
					else:
						nuc_counts = get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome)
				#return str(nuccomp_reads)
				if nuc_comp_direction == "nuc_comp_five":
					for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
						for i in range(0,readlen+1):
							for nuc in ["A","T","G","C"]:
								if i not in master_count_dict[nuc]:
									master_count_dict[nuc][i] = 0
								if readlen in nuc_counts:
									if i in nuc_counts[readlen]:
										master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
				elif nuc_comp_direction == "nuc_comp_three":
					for readlen in range(nuc_minreadlen, nuc_maxreadlen+1):
						for i in range(-1,-(nuc_maxreadlen),-1):
							for nuc in ["A","T","G","C"]:
								if i not in master_count_dict[nuc]:
									master_count_dict[nuc][i] = 0
								if readlen in nuc_counts:
									if i in nuc_counts[readlen]:
										master_count_dict[nuc][i] += nuc_counts[readlen][i][nuc]
				sqlite_db.commit()
				sqlite_db.close()
		master_dict = {"A":collections.OrderedDict(),
					   "T":collections.OrderedDict(),
					   "G":collections.OrderedDict(),
					   "C":collections.OrderedDict()}
		master_dict["A"][0] = 0
		master_dict["T"][0] = 0
		master_dict["G"][0] = 0
		master_dict["C"][0] = 0
		if nuc_comp_direction == "nuc_comp_five":
			for nuc in ["A","T","G","C"]:
				for i in range(0,nuc_maxreadlen):
					if i in master_count_dict[nuc]:
						thiscount = master_count_dict[nuc][i]
						othercount = 0.01
						for subnuc in ["A","T","G","C"]:
							othercount +=  master_count_dict[subnuc][i]
						if nuc_comp_type == "nuc_comp_per":
							master_dict[nuc][i] = float(thiscount)/float(othercount)
						elif nuc_comp_type == "nuc_comp_count":
							master_dict[nuc][i] = float(thiscount)
		elif nuc_comp_direction == "nuc_comp_three":
			for nuc in ["A","T","G","C"]:
				for i in range(-1,-(nuc_maxreadlen),-1):
					if i in master_count_dict[nuc]:
						thiscount = master_count_dict[nuc][i]
						othercount = 0.01
						for subnuc in ["A","T","G","C"]:
							othercount +=  master_count_dict[subnuc][i]
						if nuc_comp_type == "nuc_comp_per":
							master_dict[nuc][i] = float(thiscount)/float(othercount)
						elif nuc_comp_type == "nuc_comp_count":
							master_dict[nuc][i] = float(thiscount)
		title = "Nucleotide composition"
		connection.close()
		return metainfo_plots.nuc_comp(master_dict, nuc_maxreadlen,title, nuc_comp_type,nuc_comp_direction,short_code,background_col,a_col,t_col,g_col,c_col,title_size, axis_label_size, subheading_size,marker_size)
	elif plottype == "dinuc_bias":
		master_count_dict = collections.OrderedDict([("AA",0), ("AT",0), ("AG",0), ("AC",0),
													 ("TA",0), ("TT",0), ("TG",0), ("TC",0),
													 ("GA",0), ("GT",0), ("GG",0), ("GC",0),
													 ("CA",0), ("CT",0), ("CG",0), ("CC",0)])
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				dinuc_counts = sqlite_db["dinuc_counts"]
				for readlen in dinuc_counts:
					for dinuc in dinuc_counts[readlen]:
						master_count_dict[dinuc] += dinuc_counts[readlen][dinuc]

		connection.close()
		return metainfo_plots.dinuc_bias(master_count_dict,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)


	elif plottype == "fastq_screen":
		html_filepath = ""
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if html_filepath == "":
					html_filepath = filepath.replace(".sqlite","_lessrRNA_screen.html")
				else:
					return "Error: Only one dataset at a time can be selected for fastq screen"
		if os.path.isfile(html_filepath):
			openfile = open(html_filepath,"r")
			fastq_lines = openfile.readlines()
			#The base64 encoded png string in the header is too long for firefox, will work for one plot and then crash firefox, this is a hack to prevent that
			fixed_html = ""
			for line in fastq_lines:
				if "iVBORw0KGgoAAAANSUhEUgAAA4wAAAGVCAYAAAHC" in line:
					fixed_html += '<a style="float:left;" href="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen" target="_blank"><img width="50%" height="50%" alt="FastQ Screen"src="/static/fastq_screen.png"</a>'
					continue
				else:
					fixed_html += (line)
			# Replace serves two functions here, first remove the padding line in the body tag as this affects the header bar and everything else on trips,
			# but removal means fastq screen logo is slightly off screen
			# Second remove the max-width line in the .container class, replace it with the padding line removed from the body tag, as this will now be specific to the container
			# and fix the fastq screen logo.
			fixed_html = str(fixed_html.replace("padding:0 20px 20px","").replace("max-width:1200px;","padding:0 20px 20px").replace("<html>","").replace("</html>","").replace("<body>","").replace("</body>","")).replace("<!DOCTYPE html>","").replace("<head>","").replace("</head>","").replace("container","container2")
			return fixed_html
		else:
			return "No fastq_screen file available for this dataset"
	elif plottype == "explore_offsets":
		readlen_dict = {}
		traninfo_dict = SqliteDict("{0}{1}/{1}.sqlite".format(config.ANNOTATION_DIR,organism), autocommit=False)
		tranlist = traninfo_dict.keys()[:10000]
		labels = []
		f0_counts = []
		f1_counts = []
		f2_counts = []

		#For the first file in selected files
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					opendict = dict(sqlite_db)
					sqlite_db.close()
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page".format(filepath)
				tranlist = opendict.keys()

				#For each readlength we display on the final graph
				for readlen in range(25,36):
					try:
						chosen_offset = opendict["offsets"]["fiveprime"]["offsets"][readlen]
					except:
						chosen_offset = 15
					labels.append("")
					labels.append("{}_{}".format(readlen, chosen_offset))
					labels.append("")
					for offset in [chosen_offset-1,chosen_offset,chosen_offset+1]:
						trancount = 0
						inframe_counts = 0
						outframe_counts = 0
						if readlen not in readlen_dict:
							readlen_dict[readlen] = {chosen_offset-1:[0,0,0],chosen_offset:[0,0,0],chosen_offset+1:[0,0,0]}
						#for each transcript get the frame counts breakdown from the cds given this particular offset
						for tran in tranlist:
							if tran not in traninfo_dict:
								continue
							tempdict = dict(opendict[tran])
							trancount += 1
							if trancount > 5000:
								break

							if "cds_start" not in traninfo_dict[tran]:
								continue
							cds_start = traninfo_dict[tran]["cds_start"]
							cds_stop = traninfo_dict[tran]["cds_stop"]

							if cds_start == "NULL" or cds_stop == "NULL":
								continue
							if cds_start <= 1 or cds_stop <= 1:
								continue
							#to account for 0-based counts ,without this line the frame will be wrong
							cds_start += 1
							cds_frame = cds_start%3
							#first walk through this entry in the opendict for only the readlength in question applying the relevant offset
							count_dict = {}

							if readlen in tempdict["unambig"]:
								for fiveprime_pos in tempdict["unambig"][readlen]:
									count = tempdict["unambig"][readlen][fiveprime_pos]
									new_pos = fiveprime_pos + offset
									count_dict[new_pos] = count
								for i in range(cds_start,cds_stop):
									frame = i%3
									if i in count_dict:
										if frame == cds_frame:
											inframe_counts += count_dict[i]
										else:
											outframe_counts += count_dict[i]

						f0_counts.append(inframe_counts)
						f1_counts.append(outframe_counts)
						f2_counts.append(0)
		connection.close()
		return metainfo_plots.explore_offsets(f0_counts, f1_counts, f2_counts, labels,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)


	elif plottype == "metagene_plot":
		minpos = -300
		maxpos = 300
		pos_list = []
		mapped_reads_dict = {}
		for i in range(minpos, maxpos+1):
			pos_list.append(i)
		count_dict = {"fiveprime":{}, "threeprime":{}}
		for primetype in ["fiveprime","threeprime"]:
			for i in range(minpos,maxpos+1):
				count_dict[primetype][i] = 0

		if metagene_aggregate == True:
			fiveprime_counts = []
			threeprime_counts = []
		else:
			if total_files >= 7:
				return "Can only choose a maximum of 6 files if not using aggregate option"
			fiveprime_counts = {}
			threeprime_counts = {}

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				cursor.execute("SELECT file_description from files where file_id = '{}';".format(file_id))
				file_desc = cursor.fetchone()[0]
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
					mapped_reads = sqlite_db["coding_counts"]+sqlite_db["noncoding_counts"]
					mapped_reads_dict[file_desc] = float(mapped_reads)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "metagene_counts" not in sqlite_db:
					return "No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				if metagene_type == "metagene_start":
					mgc = sqlite_db["metagene_counts"]
				elif metagene_type == "metagene_stop":
					mgc = sqlite_db["stop_metagene_counts"]
				elif metagene_type == "metagene_second_aug":
					mgc = sqlite_db["secondary_metagene_counts"]
				elif metagene_type == "metagene_custom":
					mgc = create_custom_metagene(custom_seq_list,exclude_first_val,exclude_last_val,include_first_val,include_last_val,custom_search_region,exclude_first, exclude_last, include_first, include_last,sqlite_db,organism,metagene_tranlist,metagene_frame)

				if metagene_offsets == True:
					new_mgc = {"fiveprime":{},"threeprime":{}}
					offsets = sqlite_db["offsets"]
					for readlen in mgc["fiveprime"]:
						new_mgc["fiveprime"][readlen] = {}
						if readlen in offsets["fiveprime"]["offsets"]:
							offset = offsets["fiveprime"]["offsets"][readlen]
						else:
							offset = 15
						for pos in mgc["fiveprime"][readlen]:
							count = mgc["fiveprime"][readlen][pos]
							new_pos = pos+offset
							if new_pos not in new_mgc["fiveprime"][readlen]:
								new_mgc["fiveprime"][readlen][new_pos] = 0
							new_mgc["fiveprime"][readlen][new_pos] += count
					for readlen in mgc["threeprime"]:
						new_mgc["threeprime"][readlen] = {}
						if readlen in offsets["threeprime"]["offsets"]:
							offset = offsets["threeprime"]["offsets"][readlen]
						else:
							offset = -12
						for pos in mgc["threeprime"][readlen]:
							count = mgc["threeprime"][readlen][pos]
							new_pos = pos+offset
							if new_pos not in new_mgc["threeprime"][readlen]:
								new_mgc["threeprime"][readlen][new_pos] = 0
							new_mgc["threeprime"][readlen][new_pos] += count
					mgc = new_mgc
				sqlite_db.close()

				for primetype in ["fiveprime", "threeprime"]:
					for i in pos_list:
						for readlen in range(minreadlen, maxreadlen+1):
							if readlen in mgc[primetype]:
								if i in mgc[primetype][readlen]:
									count_dict[primetype][i] += (mgc[primetype][readlen][i])

				if metagene_aggregate == False:

					fiveprime_counts[file_desc] = []
					threeprime_counts[file_desc] = []
					for i in pos_list:
						fiveprime_counts[file_desc].append(count_dict["fiveprime"][i])
						threeprime_counts[file_desc].append(count_dict["threeprime"][i])
					#reset the count dict so that counts are not aggregated
					for primetype in ["fiveprime","threeprime"]:
						for i in range(minpos,maxpos+1):
							count_dict[primetype][i] = 0

		if metagene_normalise == True and metagene_aggregate == False:
			min_mapped_reads = 1000000000
			for file_desc in mapped_reads_dict:
				if mapped_reads_dict[file_desc] < min_mapped_reads:
					min_mapped_reads = mapped_reads_dict[file_desc]
			for file_desc in mapped_reads_dict:
				try:
					factor = float(min_mapped_reads/mapped_reads_dict[file_desc])
				except:
					return "Error, missing mapped reads value for one of the files so cannot normalize"
				mapped_reads_dict[file_desc] = factor
			for file_desc in fiveprime_counts:
				norm_counts = []
				for count in fiveprime_counts[file_desc]:
					factor = mapped_reads_dict[file_desc]
					normalised_count = count*factor
					norm_counts.append(normalised_count)
				fiveprime_counts[file_desc] = norm_counts
			for file_desc in threeprime_counts:
				norm_counts = []
				for count in threeprime_counts[file_desc]:
					norm_counts.append(count*mapped_reads_dict[file_desc])
				threeprime_counts[file_desc] = norm_counts

		if metagene_aggregate == True:
			for i in pos_list:
				fiveprime_counts.append(count_dict["fiveprime"][i])
				threeprime_counts.append(count_dict["threeprime"][i])
		title = "Metagene profile"
		connection.close()



		return metainfo_plots.metagene_plot(pos_list,fiveprime_counts,threeprime_counts,metagene_type,title,minpos, maxpos,short_code,background_col,metagene_fiveprime_col,metagene_threeprime_col,title_size, axis_label_size, subheading_size,marker_size,metagene_end,metagene_aggregate)

	elif plottype == "trip_periodicity":
		read_dict = {"readlengths":[],
					 "frame1":[],
					 "frame2":[],
					 "frame3":[]}
		if trip_maxreadlen < trip_minreadlen:
			return "Error: max read length less than min read length, increase max read length using the input at the top of the page."
		for i in range(trip_minreadlen, trip_maxreadlen+1):
			read_dict["readlengths"].append(i)
			read_dict["frame1"].append(0)
			read_dict["frame2"].append(0)
			read_dict["frame3"].append(0)
		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "trip_periodicity" not in sqlite_db:
					return "No triplet periodicity data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				trip_periodicity_dict = sqlite_db["trip_periodicity"]
				sqlite_db.close()
				readlen_index = 0
				for readlength in read_dict["readlengths"]:
					if readlength in trip_periodicity_dict["fiveprime"]:
						read_dict["frame1"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["0"]
						read_dict["frame2"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["1"]
						read_dict["frame3"][readlen_index] += trip_periodicity_dict["fiveprime"][readlength]["2"]
					readlen_index += 1
		title = "Triplet periodicity"
		connection.close()
		return metainfo_plots.trip_periodicity_plot(read_dict,title,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size)

	elif plottype == "mapped_reads_plot":
		labels = [""]
		unmapped = [0]
		mapped_coding = [0]
		mapped_noncoding = [0]
		ambiguous = [0]
		cutadapt_removed = [0]
		rrna_removed = [0]
		pcr_duplicates = [0]

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				filepath = file_paths_dict[filetype][file_id]

				cursor.execute("SELECT file_name,file_description from files where file_id = '{}';".format(file_id))
				result = cursor.fetchone();
				file_name = (result[0]).replace(".shelf","")
				labels.append(result[1])

				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)


				if "unmapped_reads" in sqlite_db:
					unmapped.append(sqlite_db["unmapped_reads"])
				else:
					unmapped.append(0)

				if "coding_counts" in sqlite_db:
					mapped_coding.append(sqlite_db["coding_counts"])
				else:
					mapped_coding.append(0)

				if "noncoding_counts" in sqlite_db:
					mapped_noncoding.append(sqlite_db["noncoding_counts"])
				else:
					mapped_noncoding.append(0)

				if "ambiguous_counts" in sqlite_db:
					ambiguous.append(sqlite_db["ambiguous_counts"])
				else:
					ambiguous.append(0)

				if "cutadapt_removed" in sqlite_db:
					if sqlite_db["cutadapt_removed"] != "NULL":
						cutadapt_removed.append(sqlite_db["cutadapt_removed"])
					else:
						cutadapt_removed.append(0)
				else:
					cutadapt_removed.append(0)

				if "rrna_removed" in sqlite_db:
					rrna_removed.append(sqlite_db["rrna_removed"])
				else:
					rrna_removed.append(0)

				if "pcr_duplicates" in sqlite_db:
					pcr_duplicates.append(sqlite_db["pcr_duplicates"])
				else:
					pcr_duplicates.append(0)
				sqlite_db.close()

		connection.close()
		labels.append("")

		#Append a 0 to the end of every list so that there will be an empty space on the plot at the right hand side
		for listname in [unmapped, mapped_coding, mapped_noncoding, ambiguous, cutadapt_removed, rrna_removed,pcr_duplicates]:
			listname.append(0)
		return metainfo_plots.mapped_reads_plot(unmapped, mapped_coding, mapped_noncoding, labels,ambiguous,cutadapt_removed,rrna_removed,short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,breakdown_per,pcr_duplicates)

	elif plottype == "heatmap":
		min_readlen = heatmap_minreadlen
		max_readlen = heatmap_maxreadlen
		min_pos = heatmap_startpos
		max_pos = heatmap_endpos
		master_count_list = []
		

		for filetype in file_paths_dict:
			for file_id in file_paths_dict[filetype]:
				positions = []
				readlengths = []
				count_list = []
				filepath = file_paths_dict[filetype][file_id]
				if os.path.isfile(filepath):
					sqlite_db = SqliteDict(filepath, autocommit=False)
				else:
					return "File not found: {}, please report this to tripsvizsite@gmail.com or via the contact page.".format(filepath)
				if "metagene_counts" not in sqlite_db:
					return "No metagene counts data for this file, please report this to tripsvizsite@gmail.com or via the contact page."
				if heatmap_metagene_type == "metagene_start":
					if metagene_tranlist == "":
						mgc = sqlite_db["metagene_counts"]
					else:
						mgc = create_custom_metagene("AUG", 0, 0, 3, 0, "cds", False, False, True, False,sqlite_db, organism,metagene_tranlist,"All",transcriptome)
				elif heatmap_metagene_type == "metagene_stop":
					if metagene_tranlist == "":
						mgc = sqlite_db["stop_metagene_counts"]
					else:
						mgc = create_custom_metagene("UAG,UAA,UGA", 0, 0, 0, 3, "cds", False, False, False, True,sqlite_db, organism,metagene_tranlist,"All",transcriptome)
				elif heatmap_metagene_type == "metagene_second_aug":
					mgc = sqlite_db["secondary_metagene_counts"]
				sqlite_db.close()

				for primetype in [heatmap_direction]:
					for readlen in range(max_readlen, min_readlen-1, -1):
						for i in range(min_pos, max_pos+1):
							if readlen in mgc[primetype]:
								if i in mgc[primetype][readlen]:
									if mgc[primetype][readlen][i] != 0:
										count_list.append(mgc[primetype][readlen][i])
									else:
										count_list.append(0)
								else:
									count_list.append(0)
								readlengths.append(readlen)
								positions.append(i)
							else:
								readlengths.append(readlen)
								positions.append(i)
								count_list.append(0)
				if master_count_list == []:
					master_count_list = count_list
				else:
					for i in range(0,len(count_list)):
						master_count_list[i] = master_count_list[i] + count_list[i]
								
								
								
								
								
		title = "Heatmap"
		connection.close()
		#Replace 0 values with None in master_count_list and apply log transformation if applicable
		fixed_master_count_list = []
		for i in range(0,len(master_count_list)):
			if master_count_list[i] == 0:
				fixed_master_count_list.append(None)
			else:
				if log_scale == True:
					fixed_master_count_list.append(log(master_count_list[i],2))
				else:
					fixed_master_count_list.append(master_count_list[i])
			
		return metainfo_plots.heatplot(min_readlen, max_readlen, min_pos, max_pos, positions, readlengths,fixed_master_count_list,heatmap_metagene_type,title,reverse_scale,color_palette,short_code,background_col,maxscaleval,str(title_size)+"pt", str(axis_label_size)+"pt", str(subheading_size)+"pt",str(marker_size)+"pt")

	else:
		if plottype not in ["replicate_comp"]:
			print "ERROR2 plottype is not in list",plottype
		if (plottype.strip(" ").replace("\n","")) not in ["replicate_comp"]:
			print "ERROR 3",plottype
	return "Error, unknown plot type selected: {}".format(plottype)


def get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome):
	if os.path.isfile("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome)):
		transhelve = sqlite3.connect("{0}{1}/{1}.{2}.sqlite".format(config.ANNOTATION_DIR,organism,transcriptome))
	else:
		return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	cursor = transhelve.cursor()
	cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1")
	result = cursor.fetchall()
	master_dict = {}
	offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]

	for row in result:
		tran = row[0]
		cds_start = int(row[1])
		cds_stop = int(row[2])
		seq = row[3].replace("T","U")
		if tran in sqlite_db:
			counts = sqlite_db[tran]["unambig"]
			for readlen in counts:
				if readlen not in master_dict:
					master_dict[readlen] = {}
					for i in range(0,int(readlen)):
						master_dict[readlen][i] = {"A":0,"T":0,"G":0,"C":0,"N":0}
				for pos in counts[readlen]:

					count = counts[readlen][pos]
					if readlen in offsets:
						offset_pos = pos+offsets[readlen]
					else:
						offset_pos = pos+15

					if offset_pos >= cds_start and offset_pos <= cds_stop:

						readframe = offset_pos%3
						cds_frame = (cds_start+2)%3


						if readframe == cds_frame:
							inframe = True
						else:
							inframe = False

						if nuccomp_reads == "inframe" and inframe == False:
							continue
						if nuccomp_reads == "offrame" and inframe == True:
							continue
						readseq = seq[pos:pos+readlen]
						for i in range(0,len(readseq)):
							char = readseq[i].replace("U","T")
							master_dict[readlen][i][char]+= count
	#save results so they won't have to be computed again later
	if "nuc_counts" in sqlite_db:
		new_nuc_counts = sqlite_db["nuc_counts"]
		new_nuc_counts[nuccomp_reads] = master_dict
		sqlite_db["nuc_counts"] = new_nuc_counts
	else:
		sqlite_db["nuc_counts"] = {nuccomp_reads:master_dict}


	sqlite_db.commit()
	return master_dict

# Groups together counts from different filepaths for the metainformation counts table
def aggregate_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism, all_seq_types, te_minimum_reads, html_args,count_type,te_tranlist):
	print "aggregate counts called"
	file_list = ""
	transcript_dict = {}
	table_str = ""
	filename = organism+"_translation_efficiencies_"+str(time.time())+".csv"
	table_str += filename+"?~"
	mapped_reads = 0.0
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			file_list += "{},".format(file_id)
			filepath = file_paths_dict[seq_type][file_id]
			if os.path.isfile(filepath):
				sqlite_db = SqliteDict(filepath, autocommit=False)
			else:
				return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
			try:
				mapped_reads += sqlite_db["coding_counts"]
				mapped_reads += sqlite_db["noncoding_counts"]
			except:
				pass
			if region == "all":
				opendict = sqlite_db["unambiguous_all_totals"]
			elif region == "cds":
				opendict = sqlite_db["unambiguous_cds_totals"]
			elif region == "fiveprime":
				opendict = sqlite_db["unambiguous_fiveprime_totals"]
			elif region == "threeprime":
				opendict = sqlite_db["unambiguous_threeprime_totals"]
			sqlite_db.close()
			#print"opendict", opendict
			for transcript in longest_tran_list:
				if transcript not in transcript_dict:
					transcript_dict[transcript] = {}
				if seq_type not in transcript_dict[transcript]:
					transcript_dict[transcript][seq_type] = 0
				if transcript in opendict:
					transcript_dict[transcript][seq_type] += opendict[transcript]

	#If count_type is rpkm or tpm, transform the raw counts
	if count_type == "rpkm":
		if mapped_reads == 0.0:
			return "Cannot calculate RPKM. No mapped reads data available"
		#per million scaling factor
		pmsf = mapped_reads/1000000
		for transcript in longest_tran_list:
			tranlength = float(traninfo_dict[transcript]["length"])
			if traninfo_dict[transcript]["cds_start"] != None:
				if region == "fiveprime":
					tranlength = float(traninfo_dict[transcript]["cds_start"]-1)
				elif region == "cds":
					tranlength = float((traninfo_dict[transcript]["cds_stop"]+3)-(traninfo_dict[transcript]["cds_start"]-1))
				elif region == "threeprime":
					tranlength = float(traninfo_dict[transcript]["length"]-(traninfo_dict[transcript]["cds_stop"]+3))
			# get tranlength in kilobases
			tranlength = tranlength/1000
			#print"AGGREGATE,pmsf, tranlength",pmsf, tranlength
			for seq_type in transcript_dict[transcript]:
				#print"COUNT ", transcript_dict[transcript][seq_type]
				transcript_dict[transcript][seq_type] = round((transcript_dict[transcript][seq_type]/pmsf)/tranlength,2)
	#If count_type is rpkm or tpm, transform the raw counts
	if count_type == "tpm":
		total_rpk = 0.0

		for transcript in longest_tran_list:
			tranlength = float(traninfo_dict[transcript]["length"])
			if traninfo_dict[transcript]["cds_start"] != None:
				if region == "fiveprime":
					tranlength = float(traninfo_dict[transcript]["cds_start"]-1)
				elif region == "cds":
					tranlength = float((traninfo_dict[transcript]["cds_stop"]+3)-(traninfo_dict[transcript]["cds_start"]-1))
				elif region == "threeprime":
					tranlength = float(traninfo_dict[transcript]["length"]-traninfo_dict[transcript]["cds_stop"]+3)
			# get tranlength in kilobases
			tranlength = tranlength/1000
			for seq_type in transcript_dict[transcript]:
				rpk = transcript_dict[transcript][seq_type]/tranlength
				transcript_dict[transcript][seq_type] = rpk 
				total_rpk += rpk 
		pmsf = total_rpk/1000000
		for transcript in transcript_dict:
			for seq_type in transcript_dict[transcript]:
				transcript_dict[transcript][seq_type] = round(transcript_dict[transcript][seq_type]/pmsf,2)
		print "total rpk", total_rpk


	total_rows = 0
	tmp_te_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
	tmp_te_file.write("Filename, Gene,Transcript,Region,Riboseq count, Rnaseq count, Translation efficiency")
	for seq_type in all_seq_types:
		if seq_type != "riboseq" and seq_type != "rnaseq":
			tmp_te_file.write(",{}".format(seq_type))
	tmp_te_file.write("\n")
	all_rows = []
	if te_tranlist != None:
		final_tranlist = te_tranlist
		print "'''final tranlist is te tranlist"
	else:
		final_tranlist = transcript_dict.keys()
		print "final tranlist is transcript dict.keys)("
	for transcript in final_tranlist:
		#print "FINAL TRANSCRIPT", final_tranlist
		seq_count_dict = {}
		for seq_type in transcript_dict[transcript]:
			try:
				gene = (traninfo_dict[transcript]["gene"]).replace(",","_").replace(";","_")
			except:
				gene = "Unknown"
			if seq_type not in seq_count_dict:
				seq_count_dict[seq_type] = transcript_dict[transcript][seq_type]
		if "riboseq" in seq_count_dict:
			riboseq_count = seq_count_dict["riboseq"]
		else:
			riboseq_count = 0
		if "rnaseq" in seq_count_dict:
			rnaseq_count = seq_count_dict["rnaseq"]
		else:
			rnaseq_count = 0
		if rnaseq_count < te_minimum_reads or riboseq_count < te_minimum_reads:
			continue
		if rnaseq_count == 0 or riboseq_count == 0:
			te = 0
		else:
			te = float(transcript_dict[transcript]["riboseq"])/float(transcript_dict[transcript]["rnaseq"])
			#te = round(te,2)
		tmp_te_file.write("Aggregate,{},{},{},{},{},{}".format(gene,transcript,region, riboseq_count, rnaseq_count, te))
		input_list = ["Aggregate", gene,transcript,region, riboseq_count, rnaseq_count, te]
		for seq_type in all_seq_types:
			if seq_type != "riboseq" and seq_type != "rnaseq":
				if seq_type in seq_count_dict:
					input_list.append(seq_count_dict[seq_type])
					tmp_te_file.write(",{}".format(seq_count_dict[seq_type]))
				else:
					if seq_type in file_paths_dict:
						input_list.append(0)
						tmp_te_file.write(",0")
		input_list.append("<a href='http://trips.ucc.ie/"+organism+"/"+html_args["transcriptome"]+"/interactive_plot/?tran="+transcript+"&files="+file_list+"' target='_blank_' >View plot</a>")
		all_rows.append(input_list)
		tmp_te_file.write("\n")
	tmp_te_file.close()
	os.chmod("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename), 0777)
	#if both rnaseq and riboseq files, sort by te, else sort by the relevant count
	anyfile = False
	if len(file_paths_dict["riboseq"]) != 0 and len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
	elif len(file_paths_dict["riboseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[4],reverse=True)
	elif len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[5],reverse=True)
	else:
		for seq_type in all_seq_types:
			if seq_type not in file_paths_dict:
				continue
			if len(file_paths_dict[seq_type]) != 0:
				anyfile = True
				all_sorted_rows = sorted(all_rows, key=lambda x: x[7],reverse=True)
		if anyfile == False:
			return "No files selected. Select a file by clicking on a study name in the studies section. Then select one of the files that appear in the files section."
	for row in all_sorted_rows:
		total_rows += 1
		if total_rows <= 1000:
			input_str = ""
			for item in row:
				input_str += "{}.;".format(item)
			input_str += "?~"
			table_str += input_str
	table_str = "TE?~"+str(total_rows)+"?~"+table_str
	return table_str



def sample_counts(file_paths_dict, traninfo_dict, longest_tran_list, region, organism, all_seq_types, te_minimum_reads, html_args,count_type,te_tranlist):
	file_list = ""
	transcript_dict = {}
	table_str = ""
	filename = organism+"_translation_efficiencies_"+str(time.time())+".csv"
	table_str += filename+"?~"
	mapped_reads_dict = {}
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			file_list += "{},".format(file_id)
			filepath = file_paths_dict[seq_type][file_id]
			inputfilename = (filepath.split("/")[-1]).replace(".sqlite","")
			if inputfilename not in mapped_reads_dict:
				mapped_reads_dict[inputfilename] = 0
			
			if os.path.isfile(filepath):
				sqlite_db = SqliteDict(filepath, autocommit=False)
			else:
				return "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
			try:
				mapped_reads_dict[inputfilename] += sqlite_db["coding_counts"]
				mapped_reads_dict[inputfilename] += sqlite_db["noncoding_counts"]
			except:
				pass
			if region == "all":
				try:
					opendict = sqlite_db["unambiguous_all_totals"]
				except:
					return "No unambigous totals availabel for file: {}".format("/".join(filepath.split("/")[-2:]))
			elif region == "cds":
				try:
					opendict = sqlite_db["unambiguous_cds_totals"]
				except:
					return "No unambigous totals availabel for file: {}".format("/".join(filepath.split("/")[-2:]))
			elif region == "fiveprime":
				try:
					opendict = sqlite_db["unambiguous_fiveprime_totals"]
				except:
					return "No unambigous totals availabel for file: {}".format("/".join(filepath.split("/")[-2:]))
			elif region == "threeprime":
				try:
					opendict = sqlite_db["unambiguous_threeprime_totals"]
				except:
					return "No unambigous totals availabel for file: {}".format("/".join(filepath.split("/")[-2:]))
			sqlite_db.close()
			for transcript in longest_tran_list:
				if transcript not in transcript_dict:
					transcript_dict[transcript] = {}
				if seq_type not in transcript_dict[transcript]:
					transcript_dict[transcript][seq_type] = {}
				if transcript in opendict:
					transcript_dict[transcript][seq_type][inputfilename] = {"count":opendict[transcript],"file_id":str(file_id)}
				else:
					transcript_dict[transcript][seq_type][inputfilename] = {"count":0,"file_id":str(file_id)}
					
					
					
	#If count_type is rpkm or tpm, transform the raw counts
	if count_type == "rpkm":
		for transcript in longest_tran_list:
			tranlength = float(traninfo_dict[transcript]["length"])
			if traninfo_dict[transcript]["cds_start"] != None:
				if region == "fiveprime":
					tranlength = float(traninfo_dict[transcript]["cds_start"]-1)
				elif region == "cds":
					tranlength = float((traninfo_dict[transcript]["cds_stop"]+3)-(traninfo_dict[transcript]["cds_start"]-1))
				elif region == "threeprime":
					tranlength = float(traninfo_dict[transcript]["length"]-(traninfo_dict[transcript]["cds_stop"]+3))
			# get tranlength in kilobases
			tranlength = tranlength/1000
			for seq_type in transcript_dict[transcript]:
				for inputfilename in transcript_dict[transcript][seq_type]:
					mapped_reads = float(mapped_reads_dict[inputfilename])
					if mapped_reads == 0.0:
						return "No mapped reads data for file {}".format(inputfilename)
					#per million scaling factor
					pmsf = mapped_reads/1000000
					transcript_dict[transcript][seq_type][inputfilename]["count"] = round((transcript_dict[transcript][seq_type][inputfilename]["count"]/pmsf)/tranlength,2)
	#skipped_files = []
	#If count_type is rpkm or tpm, transform the raw counts
	if count_type == "tpm":
		total_rpk_dict = {}
		for transcript in longest_tran_list:
			tranlength = float(traninfo_dict[transcript]["length"])
			if traninfo_dict[transcript]["cds_start"] != None:
				if region == "fiveprime":
					tranlength = float(traninfo_dict[transcript]["cds_start"]-1)
				elif region == "cds":
					tranlength = float((traninfo_dict[transcript]["cds_stop"]+3)-(traninfo_dict[transcript]["cds_start"]-1))
				elif region == "threeprime":
					tranlength = float(traninfo_dict[transcript]["length"]-(traninfo_dict[transcript]["cds_stop"]+3))
			# get tranlength in kilobases
			tranlength = tranlength/1000
			for seq_type in transcript_dict[transcript]:
				for inputfilename in transcript_dict[transcript][seq_type]:
					if inputfilename not in total_rpk_dict:
						total_rpk_dict[inputfilename] = 0
					rpk = transcript_dict[transcript][seq_type][inputfilename]["count"]/tranlength
					transcript_dict[transcript][seq_type][inputfilename]["count"] = rpk 
					total_rpk_dict[inputfilename] += rpk 

		for transcript in transcript_dict:
			for seq_type in transcript_dict[transcript]:
				for inputfilename in transcript_dict[transcript][seq_type]:
					pmsf = total_rpk_dict[inputfilename]/1000000
					transcript_dict[transcript][seq_type][inputfilename]["count"] = round(transcript_dict[transcript][seq_type][inputfilename]["count"]/pmsf,2)
						#except Exception as e:
						#	print "Error:",e,pmsf, total_rpk_dict[inputfilename],inputfilename, transcript,transcript_dict[transcript][seq_type][inputfilename]["count"]
		print "total rpk", total_rpk_dict
	total_rows = 0
	tmp_te_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
	tmp_te_file.write("Gene,Transcript,Region,")
	#write out headers to tmp_te_file
	transcript = transcript_dict.keys()[0]
	for seq_type in transcript_dict[transcript]:
		for inputfilename in transcript_dict[transcript][seq_type]:
			#if inputfilename not in skipped_files:
			tmp_te_file.write("{}_{},".format(inputfilename,seq_type))
	tmp_te_file.write("\n")
				
	all_rows = []

	if te_tranlist != None:
		final_tranlist = te_tranlist
	else:
		final_tranlist = transcript_dict.keys()
	for transcript in final_tranlist:
		try:
			gene = traninfo_dict[transcript]["gene"]
		except:
			gene = "Unknown"
		tmp_te_file.write("{},{},{},".format(gene,transcript,region))
		for seq_type in transcript_dict[transcript]:
			for inputfilename in transcript_dict[transcript][seq_type]:
				#if inputfilename not in skipped_files:
				count = transcript_dict[transcript][seq_type][inputfilename]["count"]
				tmp_te_file.write("{},".format(count))
		tmp_te_file.write("\n")



	for transcript in final_tranlist:
		for seq_type in transcript_dict[transcript]:
			if seq_type not in file_paths_dict:
				continue
			for inputfilename in transcript_dict[transcript][seq_type]:
				#if inputfilename not in skipped_files:
				count = transcript_dict[transcript][seq_type][inputfilename]["count"]
				file_id = transcript_dict[transcript][seq_type][inputfilename]["file_id"]
				if count < te_minimum_reads:
					continue
				try:
					gene = traninfo_dict[transcript]["gene"]
				except:
					gene = "Unknown"
				riboseq_count = 0.001
				rnaseq_count = 0.001
				if seq_type == "riboseq":
					riboseq_count = count
				elif seq_type == "rnaseq":
					rnaseq_count = count
				try:
					te = riboseq_count/rnaseq_count
				except:
					te = 0
				#tmp_te_file.write("{},{},{},{},{},{},{}".format(inputfilename,gene,transcript,region, riboseq_count, rnaseq_count, te))
				input_list = [inputfilename, gene,transcript,region, riboseq_count, rnaseq_count, te]
				if seq_type != "riboseq" and seq_type != "rnaseq":
					input_list.append(count)
				input_list.append("<a href='http://trips.ucc.ie/"+organism+"/"+html_args["transcriptome"]+"/interactive_plot/?tran="+transcript+"&files="+file_id+"' target='_blank_' >View plot</a>")
				all_rows.append(input_list)
				#tmp_te_file.write("\n")
	#tmp_te_file.close()
	os.chmod("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename), 0777)
	#if both rnaseq and riboseq files, sort by te, else sort by the relevant count
	anyfile = False
	if len(file_paths_dict["riboseq"]) != 0 and len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[6],reverse=True)
	elif len(file_paths_dict["riboseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[4],reverse=True)
	elif len(file_paths_dict["rnaseq"]) != 0:
		all_sorted_rows = sorted(all_rows, key=lambda x: x[5],reverse=True)
	else:
		for seq_type in all_seq_types:
			if seq_type not in file_paths_dict:
				continue
			if len(file_paths_dict[seq_type]) != 0:
				anyfile = True
				all_sorted_rows = sorted(all_rows, key=lambda x: x[7],reverse=True)
		if anyfile == False:
			return "No files selected, please report this to tripsvizsite@gmail.com or via the contact page."
	for row in all_sorted_rows:
		total_rows += 1
		if total_rows <= 1000:
			input_str = ""
			for item in row:
				input_str += "{}.;".format(item)
			input_str += "?~"
			table_str += input_str
	table_str = "TE?~"+str(total_rows)+"?~"+table_str
	return table_str




def calc_gc(seq):
	count_dict = {"A":0.0,"T":0.0,"G":0.0,"C":0.0}
	for char in seq:
		count_dict[char] += 1
	return count_dict
