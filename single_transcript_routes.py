from flask import Blueprint, render_template, abort, request, make_response, url_for, jsonify
from flask import current_app as app
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,fetch_user
import riboflask
import collections
from flask_login import current_user
import logging
import json
try:
	from orfQuant import incl_OPM_run_orfQuant
	from tripsTPM import TPM
except:
	pass

#This is the single transcript plot page, user chooses gene, files and other settings
single_transcript_plotpage_blueprint = Blueprint("interactiveplotpage", __name__, template_folder="templates")
@single_transcript_plotpage_blueprint.route('/<organism>/<transcriptome>/interactive_plot/')
def interactiveplotpage(organism,transcriptome):
	#global user_short_passed
	user_short_passed = True
	global local
	try:
		print (local)
	except:
		local = False

	organism = str(organism)

	
	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	
	user,logged_in = fetch_user()
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)

	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism,transcriptome))
	result = (cursor.fetchone())
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	studyinfo_dict = fetch_study_info(organism)
	user_transcript = request.args.get('tran')
	user_readscore = request.args.get('rs')
	user_hili = request.args.get('hili')
	user_generate_shorturl = request.args.get('genshort')
	user_files = request.args.get('files')
	user_minread = request.args.get('minread')
	user_maxread = request.args.get('maxread')
	user_dir = request.args.get('dir')
	user_line_graph = request.args.get('lg')
	user_ambig = request.args.get('ambig')
	user_cov = request.args.get('cov')
	user_nuc = request.args.get('nuc')
	user_short = request.args.get('short')
	user_crd = request.args.get('crd')

	if user_files != None:
		user_files = user_files.split(",")
		user_files = [str(x) for x in user_files]
	else:
		user_files = []

	user_ribo_studies = request.args.get('ribo_studies')
	if user_ribo_studies != None:
		user_ribo_studies = user_ribo_studies.split(",")
		user_ribo_studies = [str(x) for x in user_ribo_studies]
	else:
		user_ribo_studies = []
	user_proteomics_studies = request.args.get('proteomics_studies')
	if user_proteomics_studies != None:
		user_proteomics_studies = user_proteomics_studies.split(",")
		user_proteomics_studies = [str(x) for x in user_proteomics_studies]
	else:
		user_proteomics_studies = []
		
	user_rna_studies = request.args.get('rna_studies')
	if user_rna_studies != None:
		user_rna_studies = user_rna_studies.split(",")
		user_rna_studies = [str(x) for x in user_rna_studies]
	else:
		user_rna_studies = []

	if user_generate_shorturl == "F":
		user_generate_shorturl = False
	else:
		user_generate_shorturl = True

	user_hili_starts = []
	user_hili_stops = []
	try:
		for item in user_hili.split(","):
			user_hili_starts.append(int(item.split("_")[0]))
			user_hili_stops.append(int(item.split("_")[1]))
	except:
		user_hili_start = None
		user_hili_stop = None

	try:
		user_minread = int(user_minread)
		user_maxread = int(user_maxread)
	except:
		user_minread = None
		user_maxread = None
	advanced='True'
	connection.close()
	consent = request.cookies.get("cookieconsent_status")
	if consent == "deny":
		resp = make_response(render_template('index.html', gwips_info=gwips_info, gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,organism=organism,transcriptome=transcriptome,default_tran=default_tran,
						   user_transcript=user_transcript, user_readscore=user_readscore, user_hili_starts=user_hili_starts, user_hili_stops=user_hili_stops,local=local,studies_dict=accepted_studies,
						   accepted_files=accepted_files,user_files=user_files,user_ribo_studies=user_ribo_studies,user_proteomics_studies=user_proteomics_studies, user_rna_studies=user_rna_studies,user_minread=user_minread,user_maxread=user_maxread,
						   user_dir=user_dir,user_line_graph=user_line_graph,user_ambig=user_ambig,user_cov=user_cov,user_nuc=user_nuc,user_short=user_short, user_crd=user_crd,studyinfo_dict=studyinfo_dict,
						   advanced=advanced,seq_types=seq_types))
		for cookie_name in request.cookies:
			if cookie_name != "cookieconsent_status":
				resp.delete_cookie(cookie_name)
		return resp
	return render_template('index.html', gwips_info=gwips_info, gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,organism=organism,transcriptome=transcriptome,default_tran=default_tran,
						   user_transcript=user_transcript, user_readscore=user_readscore, user_hili_starts=user_hili_starts, user_hili_stops=user_hili_stops,local=local,studies_dict=accepted_studies,
						   accepted_files=accepted_files,user_files=user_files,user_ribo_studies=user_ribo_studies,user_proteomics_studies=user_proteomics_studies, user_rna_studies=user_rna_studies,user_minread=user_minread,user_maxread=user_maxread,
						   user_dir=user_dir,user_line_graph=user_line_graph,user_ambig=user_ambig,user_cov=user_cov,user_nuc=user_nuc,user_short=user_short, user_crd=user_crd,studyinfo_dict=studyinfo_dict,
						   advanced=advanced,seq_types=seq_types)



# Creates and serves the plots for the single transcript plot page
single_transcript_query_blueprint = Blueprint("query", __name__, template_folder="templates")
@single_transcript_query_blueprint.route('/query', methods=['POST'])
def query():
	#global user_short_passed
	try:
		user = current_user.name
	except:
		user = None
	#print "user", user
	tran_dict = {}
	gene_dict = {}
	ribo_user_files = {}
	print ("raw data", request.data)
	data = json.loads(request.data)
	print ("PROCESSED DATA",data)
	print (data.keys())

	tran = data['transcript'].upper().strip()
	readscore = data['readscore']
	secondary_readscore = data['secondary_readscore']
	minread = int(data['minread'])
	maxread = int(data['maxread'])
	minfiles = int(data['minfiles'])
	organism = data['organism']
	seqhili = data['seqhili'].split(",")
	hili_start = int(data['hili_start'])
	hili_stop = int(data['hili_stop'])
	transcriptome = data['transcriptome']
	advanced =  data["advanced"]
	logging.debug("FILE LIST")
	logging.debug(str(data["file_list"]))
	#logging.warn(len(data["file_list"]))
	#logging.debug("Length of alt file list is"+ len(data["alt_file_list"]))
	# Send file_list (a list of integers intentionally encoded as strings due to javascript), to be converted to a dictionary with riboseq/rnaseq lists of file paths.
	file_paths_dict = fetch_file_paths(data["file_list"],organism)

	primetype = data["primetype"]
	user_hili_starts = data["user_hili_starts"]
	user_hili_stops = data["user_hili_stops"]
	user_short = data["user_short"]

	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
	
	
	user,logged_in = fetch_user()

	if owner == 1:
		if os.path.isfile("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome)):
			sqlite_path_organism = "{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome)
			transhelve = sqlite3.connect(sqlite_path_organism)
		else:
			return_str =  "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
			if app.debug == True:
				return return_str, "NO_CELERY", {'Location': None}
			else:
				return jsonify({'current': 400, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""} 
	else:
		sqlite_path_organism = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome)
		transhelve = sqlite3.connect(sqlite_path_organism)
	cursor = transhelve.cursor()
	cursor.execute("SELECT * from transcripts WHERE transcript = '{}'".format(tran))
	
	result = cursor.fetchone()
	inputtran = True

	if result != None:
		newtran = result[0]
	else:
		inputtran = False
	if inputtran == False:
		cursor.execute("SELECT * from transcripts WHERE gene = '{}'".format(tran))
		result = cursor.fetchall()

		if result != []:
			if len(result) == 1:
				tran = str(result[0][0])
			else:
				return_str = "TRANSCRIPTS"
				if user == "test":
					return_str = "QUANT_TRANSCRIPTS"
					if len(file_paths_dict["riboseq"].values()) > 0:
						pre_orfQuant_res = incl_OPM_run_orfQuant(tran, sqlite_path_organism, file_paths_dict["riboseq"].values())
						pre_TPM_Ribo = TPM(tran, sqlite_path_organism, file_paths_dict["riboseq"].values(), "ribo")

						max_TPM_Ribo = max(pre_TPM_Ribo.values())
						TPM_Ribo = {transcript:round((pre_TPM_Ribo[transcript] / max_TPM_Ribo)*100, 2) for transcript in pre_TPM_Ribo}

						max_orf = max(pre_orfQuant_res.values())
						orfQuant_res = {transcript:round((pre_orfQuant_res[transcript] / max_orf)*100, 2) for transcript in pre_orfQuant_res}

					else:
						orfQuant_res = {transcript[0]:"Null" for transcript in result}
						TPM_Ribo = {transcript[0]:"Null" for transcript in result}
					
					if len(file_paths_dict["rnaseq"].values()) > 0:
						pre_TPM_RNA = TPM(tran, sqlite_path_organism, file_paths_dict["rnaseq"].values(), "rna")
						max_TPM_RNA = max(pre_TPM_RNA.values())
						TPM_RNA = {transcript:round((pre_TPM_RNA[transcript] / max_TPM_RNA)*100, 2) for transcript in pre_TPM_RNA}

					else:
						TPM_RNA = {transcript[0]:"Null" for transcript in result}


				for transcript in result:
					cursor.execute("SELECT length,cds_start,cds_stop,principal,version from transcripts WHERE transcript = '{}'".format(transcript[0]))
					tran_result = cursor.fetchone()
					tranlen = tran_result[0]
					cds_start = tran_result[1]
					cds_stop = tran_result[2]
					if str(tran_result[3]) == "1":
						principal = "principal"
					else:
						principal = ""
					version = tran_result[4]
					if cds_start == "NULL" or cds_start == None:
						cdslen = "NULL"
						threeutrlen = "NULL"
					else:
						cdslen = cds_stop-cds_start
						threeutrlen = tranlen - cds_stop
					if user == "test":
						if transcript[0] in orfQuant_res:
							OPM_coverage = orfQuant_res[transcript[0]]
						else:
							OPM_coverage = "NULL"

						if transcript[0] in TPM_RNA:
							RNA_coverage = TPM_RNA[transcript[0]]
						else:
							RNA_coverage = "NULL"

						if transcript[0] in TPM_Ribo:
							ribo_coverage = TPM_Ribo[transcript[0]]
						else:
							ribo_coverage = "NULL"

						return_str += (":{},{},{},{},{},{},{},{},{}".format(transcript[0],version, tranlen, cds_start, cdslen, threeutrlen, OPM_coverage, ribo_coverage, RNA_coverage))
					else:
						return_str += (":{},{},{},{},{},{},{}".format(transcript[0],version, tranlen, cds_start, cdslen, threeutrlen,principal))
				if app.debug == True:
					return return_str, "NO_CELERY", {'Location': None}
				else:
					if user == "test":
						return jsonify({'current': 400, 'total': 100, 'status': 'quant_tran_list','result': return_str}), 200, {'Location': ""} 
					else:
						return jsonify({'current': 400, 'total': 100, 'status': 'tran_list','result': return_str}), 200, {'Location': ""} 
				
		else:
			return_str =  "ERROR! Could not find any gene or transcript corresponding to {}".format(tran)
			logging.debug(return_str)
			if app.debug == True:
				return return_str, "NO_CELERY", {'Location': None}
			else:
				return jsonify({'current': 400, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""} 
	transhelve.close()
	if 'varlite' in data:
		lite = "y"
	else:
		lite="n"
	if 'preprocess' in data:
		preprocess = True
	else:
		preprocess = False
	if 'uga_diff' in data:
		uga_diff = True
	else:
		uga_diff = False
	if 'color_readlen_dist' in data:
		color_readlen_dist = True
	else:
		color_readlen_dist = False
	if 'ribocoverage' in data:
		ribocoverage = True
	else:
		ribocoverage = False
	if "nucseq" in data:
		nucseq = True
	else:
		nucseq = False
	if "mismatches" in data:
		mismatches = True
	else:
		mismatches = False
	if "ambiguous" in data:
		ambiguous = "ambig"
	else:
		ambiguous = "unambig"
	if "pcr" in data:
		pcr = True
	else:
		pcr = False
	if "noisered" in data:
		noisered = True
	else:
		noisered = False

	if "mismatch" in data:
		mismatch = True
	else:
		mismatch = False
	user_short_passed = False
	if data["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,data["transcriptome"],"interactive_plot")
	else:
		short_code = data["user_short"]
		user_short_passed = True

	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	background_col = config.BACKGROUND_COL
	uga_col = config.UGA_COL
	uag_col = config.UAG_COL
	uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	cds_marker_size = config.CDS_MARKER_SIZE
	cds_marker_colour = config.CDS_MARKER_COLOUR
	legend_size = config.LEGEND_SIZE
	ribo_linewidth = config.RIBO_LINEWIDTH
	#Put any publicly available seq types (apart from riboseq and rnaseq) here
	seq_rules = {"proteomics":{"frame_breakdown":1},"conservation":{"frame_breakdown":1},"tcpseq":{"frame_breakdown":0}}

	#get user_id
	if current_user.is_authenticated:
		user_name = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user_name))
		result = (cursor.fetchone())
		user_id = result[0]
		print ("current user id is", user_id)
		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,uga_col,uag_col,uaa_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size,ribo_linewidth from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		uga_col = result[1]
		uag_col = result[2]
		uaa_col = result[3]
		title_size = result[4]
		subheading_size = result[5]
		axis_label_size = result[6]
		marker_size = result[7]
		cds_marker_size = result[8]
		cds_marker_colour = result[9]
		legend_size = result[10]
		ribo_linewidth = result[11]
		#get rules for all custom seq types
		cursor.execute("SELECT * from seq_rules WHERE user_id = {};".format(user_id))
		result = (cursor.fetchall())
		for row in result:
			seq_name = row[1]
			frame_breakdown = row[2]
			seq_rules[seq_name] = {"frame_breakdown":frame_breakdown}
		connection.close()
	if tran != "":
		return riboflask.generate_plot(tran, ambiguous, minread, maxread, lite , ribocoverage, organism, readscore, noisered,primetype,
								minfiles,nucseq, user_hili_starts, user_hili_stops,uga_diff,file_paths_dict,short_code, color_readlen_dist,
								background_col,uga_col, uag_col, uaa_col,advanced,seqhili,seq_rules,title_size,
								subheading_size,axis_label_size,marker_size,transcriptome,config.UPLOADS_DIR,cds_marker_size,cds_marker_colour,
								legend_size,ribo_linewidth,secondary_readscore,pcr,mismatches,hili_start, hili_stop)


	else:
		return "ERROR! Could not find any transcript or gene corresponding to {}".format(tran)

