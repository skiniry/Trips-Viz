import os
import time
import sys
import sqlite3
from flask import Flask,get_flashed_messages, render_template, request, send_from_directory, flash, redirect, url_for
from flask_recaptcha import ReCaptcha
from flask_login import LoginManager, UserMixin, login_required, login_user, logout_user, current_user
from flask_security import Security
import ast
import stats_plots
from werkzeug.security import generate_password_hash, check_password_hash
from werkzeug.exceptions import BadRequest
import config
from werkzeug import secure_filename
from sqlitedict import SqliteDict
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import smtplib
from core_functions import fetch_file_paths,generate_short_code,base62_to_integer,build_profile,User
from metainfo_routes import metainfo_plotpage_blueprint, metainfoquery_blueprint
from comparison_routes import comparison_plotpage_blueprint, comparisonquery_blueprint
from single_transcript_routes import single_transcript_plotpage_blueprint, single_transcript_query_blueprint
from diff_exp_routes import diff_plotpage_blueprint,diffquery_blueprint

try:
	from orfquery_routes import translated_orf_blueprint,orfquery_blueprint
	import riboflask_datasets
	from gene_reg_routes import gene_regulation_page, gene_regulation_query
	from traninfo_routes import traninfo_plotpage_blueprint, traninfoquery_blueprint	
	app.register_blueprint(gene_regulation_page)
	app.register_blueprint(gene_regulation_query)
	app.register_blueprint(translated_orf_blueprint)
	app.register_blueprint(orfquery_blueprint)
	app.register_blueprint(traninfo_plotpage_blueprint)
	app.register_blueprint(traninfoquery_blueprint)
except:
	pass

from threading import Lock
lock = Lock()

user_short_passed = False
app = Flask(__name__, static_folder='static')
app.register_blueprint(metainfo_plotpage_blueprint)
app.register_blueprint(metainfoquery_blueprint)
app.register_blueprint(single_transcript_plotpage_blueprint)
app.register_blueprint(single_transcript_query_blueprint)
app.register_blueprint(comparison_plotpage_blueprint)
app.register_blueprint(comparisonquery_blueprint)
app.register_blueprint(diff_plotpage_blueprint)
app.register_blueprint(diffquery_blueprint)
app.config.from_pyfile('config.py')
recaptcha = ReCaptcha(app=app)
app.config['UPLOAD_FOLDER'] = '/static/tmp'
app.config['SECURITY_RECOVERABLE'] = True
app.config['SECURITY_REGISTER_URL'] = "/recoverpassword"
#change cookie name and path, this avoids cookies clashing with other flask apps on the same server 
try:
	if sys.argv[1] == "true":
		pass
except:
	app.config.update(
			SESSION_COOKIE_NAME = 'session_tripsviz',
			SESSION_COOKIE_PATH = '/'
			)

# flask-login
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"
login_manager.login_message = None



# Provides statistics on trips such as number of organisms, number of files, number of studies etc and lists updates.
@app.route('/stats/')
def statisticspage():
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()

	cursor.execute("SELECT organism_id from organisms WHERE private = 0;")
	result = (cursor.fetchall())
	no_organisms = len(result)

	public_studies = []
	cursor.execute("SELECT study_id from studies WHERE private = 0;")
	result = (cursor.fetchall())
	for row in result:
		study_id = row[0]
		public_studies.append(study_id)
	no_studies = len(public_studies)

	riboseq_files = 0
	cursor.execute("SELECT study_id,file_id from files WHERE file_type = 'riboseq';")
	result = (cursor.fetchall())
	for row in result:
		if row[0] in public_studies:
			riboseq_files += 1

	rnaseq_files = 0
	cursor.execute("SELECT study_id,file_id from files WHERE file_type = 'rnaseq';")
	result = (cursor.fetchall())
	for row in result:
		if row[0] in public_studies:
			rnaseq_files += 1

	# Create the graph which breaksdown the number of studies per organism
	org_dict = {}
	cursor.execute("SELECT organism_id,organism_name from organisms WHERE private = 0;")
	result = (cursor.fetchall())
	for row in result:
		if "_" in row[1]:
			org_name = (row[1].split("_")[0][0])+"."+(row[1].split("_")[1])
		else:
			org_name = row[1]
		org_dict[row[0]] = {"organism_name":org_name,
							"riboseq_files":0,
							"rnaseq_files":0}
	cursor.execute("SELECT file_type,organism_id from files;")
	result = (cursor.fetchall())
	for row in result:
		if row[1] in org_dict:
			if row[0] == "riboseq":
				org_dict[row[1]]["riboseq_files"] += 1
			elif row[0] == "rnaseq":
				org_dict[row[1]]["rnaseq_files"] += 1
	read_dict = {"organisms":[],
					 "riboseq_files":[],
					 "rnaseq_files":[]}
	for org in org_dict:
		read_dict["organisms"].append(org_dict[org]["organism_name"])
		read_dict["riboseq_files"].append(org_dict[org]["riboseq_files"])
		read_dict["rnaseq_files"].append(org_dict[org]["rnaseq_files"])
	org_breakdown_graph = stats_plots.org_breakdown_plot(read_dict)

	# Create the graph which breaks down studies published per year
	year_dict = {}
	cursor.execute("SELECT paper_year from studies WHERE private = 0;")
	result = cursor.fetchall()
	for row in result:
		if row[0] == "None" or row[0] == "NULL" or row[0] == "":
			continue
		if int(row[0]) not in year_dict:
			year_dict[int(row[0])] = 0
		year_dict[int(row[0])] += 1
	final_year_dict = {"year":[],
					   "no_studies":[]}
	min_year = min(year_dict.keys())
	max_year = max(year_dict.keys())
	for i in range(min_year,max_year+1):
		final_year_dict["year"].append(str(i))
		if i in year_dict:
			final_year_dict["no_studies"].append(year_dict[i])
		else:
			final_year_dict["no_studies"].append(0)
	year_plot = stats_plots.year_dist(final_year_dict)
	news_string = ""
	cursor.execute("SELECT * from updates ORDER BY date DESC;")
	result = cursor.fetchall()
	for row in result:
		news_string += "<tr><td>{}</td><td>{}</td></tr>".format(row[0], row[1])
	connection.close()
	return render_template('statistics.html', no_organisms=no_organisms, no_studies=no_studies, riboseq_files=riboseq_files, rnaseq_files=rnaseq_files,org_breakdown_graph=org_breakdown_graph,year_plot=year_plot,news_string=news_string)


# Contact page
@app.route('/contactus/',methods=["GET", "POST"])
def contactus():
	if request.method == "POST":
		name = str(request.form['name'])
		email = str(request.form['email'])
		fromaddr = "ribopipe@gmail.com"
		toaddr = "tripsvizsite@gmail.com"
		msg = MIMEMultipart()
		msg['From'] = fromaddr
		msg['To'] = toaddr
		msg['Subject'] = str(request.form['subject'])
		msg.attach(MIMEText("Name: {}\nEmail: {}\nMessage: {}".format(str(request.form['name']),str(request.form['email']),str(request.form['message']))))
		server = smtplib.SMTP('smtp.gmail.com', 587)
		server.starttls()
		server.login(fromaddr, config.EMAIL_PASS)
		text = msg.as_string()
		server.sendmail(fromaddr, toaddr, text)
		server.quit()
		flash("Message sent successfully")
	return render_template('contact.html')


#Allows users to change some global settings such as plot background colour, title size, tick label size, etc. 
@app.route('/settings/')
@login_required
def settingspage():
	global local
	try:
		print local
	except:
		local = False
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	cursor.execute("SELECT background_col,readlength_col,metagene_fiveprime_col,metagene_threeprime_col,nuc_comp_a_col,nuc_comp_t_col,nuc_comp_g_col,nuc_comp_c_col,uag_col,uaa_col,uga_col,comp_uag_col,comp_uaa_col,comp_uga_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size,ribo_linewidth from user_settings WHERE user_id = {};".format(user_id))
	result = (cursor.fetchone())
	background_colour = result[0]
	readlength_colour = result[1]
	metagene_fiveprime_colour = result[2]
	metagene_threeprime_colour = result[3]
	nuc_comp_a_col = result[4]
	nuc_comp_t_col = result[5]
	nuc_comp_g_col = result[6]
	nuc_comp_c_col = result[7]
	uag_col = result[8]
	uaa_col = result[9]
	uga_col = result[10]
	comp_uag_col = result[11]
	comp_uaa_col = result[12]
	comp_uga_col = result[13]
	title_size = result[14]
	subheading_size = result[15]
	axis_label_size = result[16]
	marker_size = result[17]
	cds_marker_size = result[18]
	cds_marker_colour = result[19]
	legend_size = result[20]
	ribo_linewidth = result[21]
	connection.close()
	return render_template('settings.html',
						   local=local,
						   background_colour=background_colour,
						   readlength_colour=readlength_colour,
						   metagene_fiveprime_colour=metagene_fiveprime_colour,
						   metagene_threeprime_colour=metagene_threeprime_colour,
						   nuc_comp_a_col=nuc_comp_a_col,
						   nuc_comp_t_col=nuc_comp_t_col,
						   nuc_comp_g_col=nuc_comp_g_col,
						   nuc_comp_c_col=nuc_comp_c_col,
						   uag_col=uag_col,
						   uaa_col=uaa_col,
						   uga_col=uga_col,
						   comp_uag_col=comp_uag_col,
						   comp_uaa_col=comp_uaa_col,
						   comp_uga_col=comp_uga_col,
						   title_size=title_size,
						   subheading_size=subheading_size,
						   axis_label_size=axis_label_size,
						   marker_size=marker_size,
						   cds_marker_width=cds_marker_size,
						   cds_marker_colour=cds_marker_colour,
						   legend_size=legend_size,
						   ribo_linewidth=ribo_linewidth)




# Allows users to download fasta files, as well as scripts needed to produce their own sqlite files
@app.route('/downloads/')
def downloadspage():
	global local
	try:
		print local
	except:
		local = False
	organism_dict = {"Scripts":["bam_to_sqlite.py","tsv_to_sqlite.py","create_annotation_sqlite.py","create_transcriptomic_to_genomic_sqlite.py"]}
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	try:
		user = current_user.name
	except:
		user = None

	cursor.execute("SELECT organism_name from organisms where private = 0")
	result = cursor.fetchall()
	for row in result:
		organism = row[0]
		organism_dict[organism] = []
	trips_annotation_dir = "{}/trips_annotations/".format(config.SCRIPT_LOC)
	for org in os.listdir(trips_annotation_dir):
		if org not in organism_dict:
			continue
		for filename in os.listdir(trips_annotation_dir+"/"+org):
			if "." in filename:
				ext = filename.split(".")[-1]
				if ext == "fa" or ext == "gtf":
					organism_dict[org].append(filename)
				elif ext == "sqlite":
					if "transcriptomic" in filename or org in filename:
						organism_dict[org].append(filename)

	connection.close()
	return render_template('downloads.html',
						   local=local,
						   user=user,
						   organism_dict=organism_dict
						   )

# Called when user downloads something from the downloads page
@app.route('/downloadquery', methods = ['GET', 'POST'])
def download_file():
	if request.method == 'POST':
		organism = request.form["organism"]
		assembly = request.form["assembly"]
		return send_from_directory("{}/trips_annotations/{}".format(config.SCRIPT_LOC, organism), assembly, as_attachment=True)



# Allows users to upload their own sqlite files and transcriptomes. 
@app.route('/uploads/')
@login_required
def uploadspage():
	global local
	try:
		print local
	except:
		local = False
	organism_dict = {}
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	user = current_user.name
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	org_id_dict = {}
	cursor.execute("SELECT organism_name,transcriptome_list,organism_id from organisms where private = 0 OR owner = {};".format(user_id))
	result = cursor.fetchall()
	for row in result:
		organism = row[0]
		transcriptome = row[1]
		print ("organism, transcriptome", organism, transcriptome)
		if organism not in organism_dict:
			organism_dict[organism] = [transcriptome]
		else:
			organism_dict[organism].append(transcriptome)
		org_id_dict[row[2]] = [organism,transcriptome]
	cursor.execute("SELECT organism_id from organism_access WHERE user_id = '{}';".format(user_id))
	result = (cursor.fetchall())
	for row in result:
		organism_id = row[0]
		cursor.execute("SELECT organism_name,transcriptome_list from organisms where organism_id = {};".format(organism_id))
		org_result = cursor.fetchone()
		organism_name = org_result[0]
		transcriptome = org_result[1]
		if organism_name not in organism_dict:
			organism_dict[organism_name] = [transcriptome]
		else:
			organism_dict[organism_name].append(transcriptome)
		org_id_dict[organism_id] = [organism_name,transcriptome]	
	#print ("organism dict0", organism_dict)
	study_dict = {}
	cursor.execute("SELECT study_id,study_name,organism_id from studies where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		#add to organism dict if not caught earlier (this only happens when database is modified manually)
		if row[2] not in org_id_dict:
			cursor.execute("SELECT organism_name,transcriptome_list,organism_id from organisms where organism_id =  {};".format(row[2]))
			#organism_dict[organism] = [transcriptome]
			org_id_dict[row[2]] = [row[0],transcriptome]
		study_dict[int(row[0])] = [row[1].replace("_{}".format(user_id),"",1),org_id_dict[row[2]][0],org_id_dict[row[2]][1],[]]
	#print ("organism dict1", organism_dict)
	transcriptome_dict = {}
	cursor.execute("SELECT organism_id,organism_name,transcriptome_list from organisms where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		transcriptome_dict[int(row[0])] = [row[1],row[2]]
	for study_id in study_dict:
		cursor.execute("SELECT user_id FROM study_access WHERE study_id = {}".format(study_id))
		result = cursor.fetchall()
		for row in result:
			cursor.execute("SELECT username FROM users WHERE user_id = {}".format(row[0]))
			result2 = cursor.fetchone()
			shared_user = result2[0]
			if shared_user not in study_dict[study_id][3]:
				study_dict[study_id][3].append(result2[0])
	file_dict = {}
	cursor.execute("SELECT file_name,study_id,file_id,file_description from files where owner = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		print "row", row, row[0], row[1], row[2], row[3]
		cursor.execute("SELECT study_name from studies where study_id = {}".format(row[1]))
		study_name = (cursor.fetchone())
		if study_name != None:
			file_dict[row[0]] = [study_name[0].replace("_{}".format(user_id),"",1),row[2],row[3]]
	seq_dict = {}
	cursor.execute("SELECT seq_name,frame_breakdown from seq_rules where user_id = {}".format(user_id))
	result = cursor.fetchall()
	for row in result:
		seq_dict[row[0]] = [row[1]]
	connection.close()
	#print ("passing file dict to uploads.html", file_dict)
	print "organism_dict", organism_dict
	print "transcriptome dict", transcriptome_dict
	return render_template('uploads.html',
						   local=local,
						   user=user,
						   organism_dict=organism_dict,
						   study_dict=study_dict,
						   transcriptome_dict=transcriptome_dict,
						   file_dict=file_dict,
						   seq_dict=seq_dict
						   )


# Called when user uploads something on the uploads page
@app.route('/uploadquery', methods = ['GET', 'POST'])
@login_required
def upload_file():
	if request.method == 'POST':
		uploaded_files = request.files.getlist("file")
		for f in uploaded_files:
			connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
			connection.text_factory = str
			cursor = connection.cursor()
			user = current_user.name
			cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
			result = (cursor.fetchone())
			user_id = result[0]
			filename = secure_filename(f.filename)
			file_ext = (filename.split("."))[-1]
			if file_ext != "sqlite":
				flash("Error: File extension should be sqlite not {}".format(file_ext))
				return redirect("https://trips.ucc.ie/uploads")
			if user != "public":
				foldername = "{}_{}".format(request.form["foldername"].replace(" ","_"),user_id)
			else:
				foldername = "{}".format(request.form["foldername"].replace(" ","_"))
			organism = request.form["organism"]
			assembly = request.form["assembly"]
			filetype_radio = request.form["filetype"]
			if filetype_radio == "riboseq":
				filetype = "riboseq"
			elif filetype_radio == "rnaseq":
				filetype = "rnaseq"
			elif filetype_radio == "other":
				filetype = (request.form["seq_type"]).lower().strip()
			
			#if this filetype is new for this user insert a new entry into seq_rules table
			if filetype != "riboseq" and filetype != "rnaseq":
				cursor.execute("SELECT * from seq_rules where user_id = {} and seq_name = '{}'".format(user_id,filetype))
				result = cursor.fetchone()
				if result == None:
					cursor.execute("INSERT INTO seq_rules VALUES ({},'{}',0)".format(user_id,filetype))

			if not os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC, foldername)):
				os.makedirs("{}/uploads/{}".format(config.SCRIPT_LOC, foldername))
			upload_file_path = "{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename)
			f.save("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename))
			sqlite_db = SqliteDict("{}/uploads/{}/{}".format(config.SCRIPT_LOC, foldername,filename))
			try:
				file_description = sqlite_db["description"]
			except:
				file_description = "NULL"
			sqlite_db.close()
			#get file id
			cursor.execute("SELECT MAX(file_id) FROM files;")
			result = cursor.fetchone();
			new_file_id = int(result[0])+1
			cursor.execute("SELECT organism_id FROM organisms WHERE organism_name = '{}' AND transcriptome_list = '{}';".format(organism, assembly))
			result = cursor.fetchone();
			organism_id = int(result[0])
			#get study_id
			cursor.execute("SELECT study_id FROM studies WHERE study_name = '{}' and organism_id = {} and owner = {}".format(foldername, organism_id, user_id))
			result = cursor.fetchone();
			if result != None:
				study_id = int(result[0])
			else:
				cursor.execute("SELECT MAX(study_id) FROM studies;")
				result = cursor.fetchone();
				study_id = int(result[0])+1
				if user != "public":
					cursor.execute("INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})".format(study_id,organism_id,foldername,'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',1,user_id))
				else:
					cursor.execute("INSERT INTO studies VALUES({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}',{},{})".format(study_id,organism_id,foldername,'NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL',0,user_id))
				#cursor.execute("SELECT study_access from users WHERE user_id = {}".format(user_id))
				#result = cursor.fetchone()
				#study_access_list = result[0]
				#if str(study_id) not in study_access_list.split(","):
					#study_access_list += ",{}".format(study_id)
					#cursor.execute("UPDATE users SET study_access = '{}' WHERE user_id = {}".format(study_access_list, user_id))
				cursor.execute("INSERT INTO study_access VALUES({},{});".format(study_id,user_id))
			cursor.execute("INSERT INTO files VALUES({},{},{},'{}','{}','{}',{},{},{},'{}')".format(new_file_id,organism_id,study_id,filename,file_description,filetype,user_id,0,0,""))
			connection.commit()
			connection.close()
			flash("File uploaded successfully")
		return redirect("https://trips.ucc.ie/uploads")



# Called when a user uploads a custom transcriptome
@app.route('/uploadtranscriptome', methods = ['GET', 'POST'])
@login_required
def upload_transcriptome():
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	user = current_user.name
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	if request.method == 'POST':
		organism = (request.form["organism"]).lower().strip().replace(" ","_")
		assembly = (request.form["assembly"]).lower().strip().replace(" ","_")
		default_tran = (request.form["default_tran"]).lower().strip().replace(" ","_")
		uploaded_annotation = request.files.getlist("anno_file")
		if not os.path.isdir("{}/uploads/transcriptomes/{}".format(config.SCRIPT_LOC, user_id)):
			os.makedirs("{}/uploads/transcriptomes/{}".format(config.SCRIPT_LOC, user_id))
		if not os.path.isdir("{}/uploads/transcriptomes/{}/{}".format(config.SCRIPT_LOC, user_id, organism)):
			os.makedirs("{}/uploads/transcriptomes/{}/{}".format(config.SCRIPT_LOC, user_id, organism))
		if not os.path.isdir("{}/uploads/transcriptomes/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly)):
			os.makedirs("{}/uploads/transcriptomes/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly))
		for f in uploaded_annotation:
			filename = secure_filename(f.filename)
			ext = filename.split(".")[-1]
			if ext != "sqlite":
				return """Error: Expecting extension sqlite but got extension {}. The file generated by the create_annotation_sqlite.py script should be uploaded here.
						This script can be gotten on the downloads page, by selecting the Scripts group.""".format(ext)
			#Instead of using filename of the uploaded file we rename it to organism_assembly.sqlite, to keep things consistent
			filename = "{}_{}.sqlite".format(organism, assembly)
			f.save("{}/uploads/transcriptomes/{}/{}/{}/{}".format(config.SCRIPT_LOC, user_id, organism, assembly,filename))
			cursor.execute("SELECT MAX(organism_id) from organisms;")
			max_org_id = (cursor.fetchone()[0])+1
			cursor.execute("INSERT INTO organisms VALUES({},'{}','{}','NULL','NULL','NULL','NULL','{}',1,{})".format(max_org_id, organism, assembly,default_tran,user_id))
			connection.commit()
			connection.close()
		flash("File uploaded successfully")
		return redirect("https://trips.ucc.ie/uploads")


# Called by flask in case of an error in the code, returns the exception so it can be displayed to user
@app.errorhandler(500)
def handle_bad_request(e):
	return 'ERROR: '+str(e)+" please report this to tripsvizsite@gmail.com or via the contact page. "

# This is the page where users login. 
@app.route("/user/login", methods=["GET", "POST"])
def login():
	global local
	#if user is already logged in then redirect to homepage
	if current_user.is_authenticated:
		return redirect("/")
	error=None
	if request.method == 'POST':
		username = str(request.form['username']).strip()
		password = str(request.form['password']).strip()
		if recaptcha.verify() or local == True or username == "developer":
			username_dict = {}
			connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT username,password from users;")
			result = (cursor.fetchall())
			connection.close()
			for row in result:
				username_dict[row[0]] = row[1]
			if username in username_dict:
				if check_password_hash(username_dict[username],password) == True or local == True:
					id = username
					user = User(id)
					login_user(user)
					nxt = request.args.get('next')
					if nxt != None:
						if "<function login" in nxt:
							nxt = "/"
					else:
						nxt = "/"
					return redirect(nxt)
				else:
					error = 'Either username or password incorrect. Please try again.'
					return render_template('login.html',error=error)
			else:
				error = 'Either username or password incorrect. Please try again.'
				return render_template('login.html',error=error)
		else:
			error = 'Invalid Captcha. Please try again.'
			return render_template('login.html',error=error)
	else:
		return render_template('login.html',error=error)

# This is the page where users create a new login. 
@app.route("/create", methods=["GET", "POST"])
def create():
	#if user is already logged in then redirect to homepage
	if current_user.is_authenticated:
		return redirect("/")
	error=None
	if request.method == 'POST':
		username = str(request.form['username'])
		password = str(request.form['password'])
		password2 = str(request.form['password2'])
		if recaptcha.verify() or local == True:
			username_dict = {}
			connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT username,password from users;")
			result = (cursor.fetchall())
			connection.close()
			for row in result:
				username_dict[row[0]] = row[1]
			if username in username_dict:
				error = "Error: {} is already registered".format(username)
				return render_template('create.html',error=error)
			if password == "":
				error = "Password cannot be empty"
				return render_template('create.html',error=error)
			if password != password2:
				error = "Passwords do not match"
				return render_template('create.html',error=error)
			hashed_pass = generate_password_hash(password)
			connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT MAX(user_id) from users;")
			result = cursor.fetchone()
			max_user_id = int(result[0])
			user_id = max_user_id+1
			#Add -1 to study access list, causes problems when adding study id's later if we don't
			cursor.execute("INSERT INTO users VALUES ({},'{}','{}','-1','',0);".format(user_id, username, hashed_pass))
			cursor.execute("INSERT INTO user_settings VALUES ('20','20','32','18',{},'#F2F2F7','#FF5F5B','#FF5F5B','#9ACAFF','#FF5F5B','#90E090','#9ACAFF','#FFFF91','gray','gray','gray','gray','gray','gray',2,'gray',17,1);".format(user_id))
			connection.commit()
			connection.close()
			return redirect("/")
		else:
			error = 'Invalid Captcha. Please try again.'
			return render_template('create.html',error=error)
	else:
		return render_template('create.html',error=error)

# Called when user presses the save button on the orf_translation page. 
@app.route('/anno_query', methods=['POST'])
def anno_query():
	data = ast.literal_eval(request.data)
	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	try:
		user = current_user.name
	except:
		user = None
		return "Error user not signed in"
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	cursor.execute("INSERT INTO users_saved_cases VALUES('{}','{}',{},{},{},{},'{}','START_CODON','CDS_OVERLAP','START_SCORE','STOP_SCORE','ENTROPY','TE','COVERAGE','CDS_RATIO','{}','FILE_LIST',{},'{}','{}','{}');".format(
																						data["gene"].strip(),
																						data["transcript"].strip(),
																						data["start"].strip(),
																						data["stop"].strip(),
																						data["length"].strip(),
																						data["score"].strip(),
																						data["type"].strip(),
																						data["trips_link"].strip(),
																						user_id,
																						data["label"].strip(),
																						data["organism"].strip(),
																						data["transcriptome"].strip()
																						))
	connection.commit()
	connection.close()
	return ""



# This page shows the saved ORFs specific to the signed in user
@app.route('/saved/')
def saved():
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except:
		user = None

	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	advanced = True
	if user != None:
		cursor.execute("SELECT advanced from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		if result[0] == 1:
			advanced = True
		else:
			advanced = False

	if user != None:
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
	cursor.execute("SELECT organism_id from organism_access WHERE user_id = '{}';".format(user_id))
	result = (cursor.fetchall())
	organism_access_list =[]
	for row in result:
		organism_access_list.append(row[0])
	cursor.execute("SELECT organism_id,organism_name,private from organisms;")
	organism_list = []
	# List of all rows returned
	result = (cursor.fetchall())
	for row in result:
		if row[2] == 0:
			organism_list.append(str(row[1]))
		elif row[2] == 1:
			if row[0] in organism_access_list:
				organism_list.append(str(row[1]))
	connection.close()
	return render_template('user_saved_cases.html',local=local,advanced=advanced,organism_list=organism_list)


# Retrieves saved ORFs
@app.route('/savedquery', methods=['POST'])
def savedquery():
	start_time = time.time()
	acceptable = 0
	data = ast.literal_eval(request.data)
	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	organism = data["organism"]
	label = data["label"]
	returnstr = ""
	try:
		user = current_user.name
	except:
		user = None
		return "Error user not signed in"
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	start_codons = []
	if "sc_aug" in data:
		start_codons.append("ATG")
	if "sc_cug" in data:
		start_codons.append("CTG")
	if "sc_gug" in data:
		start_codons.append("GTG")
	if "sc_none" in data:
		start_codons.append("None")

	# structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
	accepted_orf_dict = {}
	if organism != 'Select an Organism':
		if label == "":
			cursor.execute("SELECT * FROM users_saved_cases WHERE user_id = {} and organism = '{}';".format(user_id,organism))
		else:
			cursor.execute("SELECT * FROM users_saved_cases WHERE user_id = {} AND label = '{}' and organism = '{}';".format(user_id,label,organism))
	else:
		if label == "":
			cursor.execute("SELECT * FROM users_saved_cases WHERE user_id = {};".format(user_id))
		else:
			cursor.execute("SELECT * FROM users_saved_cases WHERE user_id = {} AND label = '{}';".format(user_id,label))
	result = cursor.fetchall()
	total_rows = 0
	for row in result:
		gene = str(row[0])
		tran = str(row[1])
		start = str(row[2])
		stop = str(row[3])
		length = str(row[4])
		score = str(row[5])
		label = str(row[18])
		start_codon = str(row[7])
		cds_overlap = str(row[8])
		start_score = str(row[9])
		stop_score = str(row[10])
		entropy = str(row[11])
		te = str(row[12])
		coverage = str(row[13])
		cds_ratio = str(row[14])
		trips_link = str(row[15])

		# Limit the number of returned cases to 1000, as datatable is memory intensive 
		if total_rows <1000:
			returnstr += "{},{},{},{},{},{},{},{},{}.,/".format(gene, tran, start, stop, length,score, label, start_codon, trips_link)
			total_rows += 1
		else:
			break
	return returnstr



# Allows users to delete previously saved cases
@app.route('/del_query', methods=['POST'])
def del_query():
	data = ast.literal_eval(request.data)
	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	try:
		user = current_user.name
	except:
		user = None
		return "Error user not signed in"
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	cursor.execute("DELETE from users_saved_cases WHERE tran = '{}' and start = '{}' and stop = '{}' and trips_link = '{}' and label = '{}' and user_id = {};".format(data["transcript"],data["start"],data["stop"],data["trips_link"],data["label"],user_id))
	connection.commit()
	connection.close()
	return ""


# Allows users to logout
@app.route("/user/logout")
@login_required
def logout():
	logout_user()
	return redirect(login)


# callback to reload the user object
@login_manager.user_loader
def load_user(userid):
	return User(userid)

# Points to robots.txt in static folder 
@app.route('/robots.txt')
def static_from_root():
	return send_from_directory(app.static_folder, request.path[1:])

#This is the help page, linked from various other pages to explain terms on that page.
@app.route('/help/')
def helppage():
	parent_acc = request.args.get('parent_acc')
	child_acc = request.args.get('child_acc')
	return render_template('help.html', parent_acc=parent_acc, child_acc=child_acc)




#This is the short url page, user supplies a short code which will be converted to a full url which user will then be redirected to
@app.route('/short/<short_code>/')
def short(short_code):
	try:
		user = current_user.name
	except:
		user = None
	#First convert short code to an integer
	integer = base62_to_integer(short_code)
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT url from urls WHERE url_id = '{}';".format(integer))
	result = cursor.fetchone()
	if result == None:
		return "Short code not recognized."
	url = ""
	url += result[0]
	#add a keyword to the url to prevent generating another shortcode
	url += "&short={}".format(short_code)
	connection.close()
	return redirect(url)

#This is the home page it show a list of organisms as defined by trips_dict
@app.route('/')
def homepage():
	try:
		user = current_user.name
	except:
		user = None
	organism_access_list = []
	organism_list = []
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	user_id = -1
	if user != None:
		flash("You are logged in as {}".format(user))
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		#get a list of organism id's this user can access
		cursor.execute("SELECT organism_id from organism_access WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchall())
		for row in result:
			organism_access_list.append(int(row[0]))
				
	#returns a tuple with each field as a seperate string
	cursor.execute("SELECT organism_id,organism_name,private,owner from organisms;")
	# List of all rows returned
	result = (cursor.fetchall())
	for row in result:
		if row[2] == 0:
			if row[1] not in organism_list:
				organism_list.append(row[1])
		elif row[2] == 1:
			if row[0] in organism_access_list or row[3] == user_id:
				if row[1] not in organism_list:
					organism_list.append(row[1])
	organism_list.sort()
	connection.close()
	return render_template('landing.html',organisms=organism_list)



#This is the home page it show a list of organisms as defined by trips_dict
@app.route('/predictions/')
def predictionspage():
	studies = ["park","xu","combo","battle"]
	return render_template('index_predictions.html',studies=studies)

#show a list of transcriptomes
@app.route('/<organism>/')
def transcriptomepage(organism):
	transcriptomes = []

	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()

	cursor.execute("SELECT transcriptome_list from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchall())
	if result:
		for row in result:
			transcriptomes.append(row[0])
	transcriptomes = str(transcriptomes).strip("[]").replace("'","")
	return render_template('transcriptomes.html', transcriptomes=transcriptomes)

#show a list of plot types
@app.route('/<organism>/<transcriptome>/')
def plogpage(organism,transcriptome):
	try:
		user = current_user.name
	except:
		user = None
	return render_template('plot_types.html',current_username=user)



# Updates the settings for a specific user
@app.route('/settingsquery', methods=['GET','POST'])
@login_required
def settingsquery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	new_password = data['new_password']
	new_password2 = data['new_password2']
	curr_password = data['curr_password']
	if new_password != "" or new_password2 != "":
		if new_password != new_password2:
			return "ERROR: New passwords do not match"
		curr_password_hash = generate_password_hash(curr_password)
		cursor.execute("SELECT password FROM users WHERE username = '{}'".format(user))
		result = cursor.fetchone()
		old_password_hash = result[0]
		if check_password_hash(old_password_hash,curr_password) == True:
			new_password_hash = generate_password_hash(new_password)
			cursor.execute("UPDATE users SET password = '{}' WHERE username = '{}'".format(new_password_hash, user))
			connection.commit()
		else:
			return "ERROR: Current password is not correct"
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]

	if "advanced" in data:
		cursor.execute("UPDATE users SET advanced = 1 WHERE user_id = '{}';".format(user_id))
		connection.commit()
	else:
		cursor.execute("UPDATE users SET advanced = 0 WHERE user_id = '{}';".format(user_id))
		connection.commit()

	#get a list of organism id's this user can access
	cursor.execute("UPDATE user_settings SET background_col = '{}' WHERE user_id = '{}';".format(data["background_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET title_size = '{}' WHERE user_id = '{}';".format(data["title_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET subheading_size = '{}' WHERE user_id = '{}';".format(data["subheading_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET axis_label_size = '{}' WHERE user_id = '{}';".format(data["axis_label_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET marker_size = '{}' WHERE user_id = '{}';".format(data["marker_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET legend_size = '{}' WHERE user_id = '{}';".format(data["legend_size"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET ribo_linewidth = '{}' WHERE user_id = '{}';".format(data["ribo_linewidth"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET cds_marker_width = '{}' WHERE user_id = '{}';".format(data["cds_marker_width"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET cds_marker_colour = '{}' WHERE user_id = '{}';".format(data["cds_marker_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET readlength_col = '{}' WHERE user_id = '{}';".format(data["readlength_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET metagene_fiveprime_col = '{}' WHERE user_id = '{}';".format(data["metagene_fiveprime_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET metagene_threeprime_col = '{}' WHERE user_id = '{}';".format(data["metagene_threeprime_colour"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_a_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_a_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_t_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_t_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_g_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_g_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET nuc_comp_c_col = '{}' WHERE user_id = '{}';".format(data["nuc_comp_c_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET uag_col = '{}' WHERE user_id = '{}';".format(data["uag_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET uga_col = '{}' WHERE user_id = '{}';".format(data["uga_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET uaa_col = '{}' WHERE user_id = '{}';".format(data["uaa_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET comp_uag_col = '{}' WHERE user_id = '{}';".format(data["comp_uag_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET comp_uga_col = '{}' WHERE user_id = '{}';".format(data["comp_uga_col"].strip(),user_id))
	connection.commit()
	cursor.execute("UPDATE user_settings SET comp_uaa_col = '{}' WHERE user_id = '{}';".format(data["comp_uaa_col"].strip(),user_id))
	connection.commit()
	connection.close()
	flash("Settings have been updated")
	return "Settings have been updated"

# Allows users to delete files
@app.route('/deletequery', methods=['GET','POST'])
@login_required
def deletequery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	#print ("DATA", data)
	for key in data:
		file_id = data[key]["file_id"]
		if "filecheck" in data[key]:
			cursor.execute("SELECT * FROM files WHERE file_id = {}".format(file_id))
			result = cursor.fetchone()
			study_id = result[2]
			filename = result[3]
			cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
			result = cursor.fetchone()
			study_name = result[2]
			full_path = "{}{}/{}".format(config.UPLOADS_DIR,study_name, filename)
			#os.remove(full_path)
			cursor.execute("DELETE FROM files WHERE file_id = {}".format(file_id))
		cursor.execute("UPDATE files SET file_description = '{}' WHERE file_id = {}".format(data[key]["file_desc"] ,file_id))
		if data[key]["cutadapt_removed"] != '0':
			cursor.execute("SELECT organism_id FROM files WHERE file_id = {}".format(file_id))
			result = cursor.fetchone()
			cursor.execute("SELECT organism_name FROM organisms WHERE organism_id = {}".format(result[0]))
			organism = cursor.fetchone()[0]
			filepath_dict = fetch_file_paths([file_id],organism)
			for seq_type in filepath_dict:
				if file_id in filepath_dict[seq_type]:
					filepath = filepath_dict[seq_type][file_id]
					print "filepath", filepath
					opendict = SqliteDict(filepath,autocommit=True)
					print "BEFORE", opendict["cutadapt_removed"]
					opendict["cutadapt_removed"] = int(data[key]["cutadapt_removed"])
					print "AFTER", opendict["cutadapt_removed"]
					print "cutadapt removed", int(data[key]["cutadapt_removed"])
					opendict.close()
			
		if data[key]["rrna_removed"] != '0':
			cursor.execute("SELECT organism_id FROM files WHERE file_id = {}".format(file_id))
			result = cursor.fetchone()
			cursor.execute("SELECT organism_name FROM organisms WHERE organism_id = {}".format(result[0]))
			organism = cursor.fetchone()[0]
			filepath_dict = fetch_file_paths([file_id],organism)
			for seq_type in filepath_dict:
				if file_id in filepath_dict[seq_type]:
					filepath = filepath_dict[seq_type][file_id]
					opendict = SqliteDict(filepath,autocommit=True)
					opendict["rrna_removed"] = int(data[key]["rrna_removed"])
					opendict.close()
		if data[key]["unmapped"] != '0':
			cursor.execute("SELECT organism_id FROM files WHERE file_id = {}".format(file_id))
			result = cursor.fetchone()
			cursor.execute("SELECT organism_name FROM organisms WHERE organism_id = {}".format(result[0]))
			organism = cursor.fetchone()[0]
			filepath_dict = fetch_file_paths([file_id],organism)
			for seq_type in filepath_dict:
				if file_id in filepath_dict[seq_type]:
					filepath = filepath_dict[seq_type][file_id]
					opendict = SqliteDict(filepath,autocommit=True)
					opendict["unmapped_reads"] = int(data[key]["unmapped"])
					opendict.close()
		
		
	connection.commit()
	connection.close()
	flash("Files have been deleted")
	return redirect("https://trips.ucc.ie/uploads")


# Allows users to delete studies,modify access, modify the organism/transcriptome assembly or study name
@app.route('/deletestudyquery', methods=['GET','POST'])
@login_required
def deletestudyquery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	for study_id in data:
		studycheck = data[study_id][0]
		#Delete studies where the "delete" checkbox is checked
		if studycheck.split("_")[-1] != "undefined":
			study_id = studycheck.split("_")[-1]
			#First delete all files on the server associated with this study, if there are any
			cursor.execute("SELECT * FROM files WHERE study_id = {}".format(study_id))
			result = cursor.fetchall()
			if result != None:
				for row in result:
					file_id = row[0]
					filename = row[3]
					cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
					result = cursor.fetchone()
					study_name = row[2]
					full_path = "{}{}/{}".format(config.UPLOADS_DIR,study_name, filename)
					if os.path.isfile(full_path):
						os.remove(full_path)
			#Now remove the study and the files associated with it from the db
			cursor.execute("DELETE FROM studies WHERE study_id = {}".format(study_id))
			cursor.execute("DELETE FROM files WHERE study_id = {}".format(study_id))
			connection.commit()
		#Modify access list next to studies 
		study_access = data[study_id][1].split(",")
		#check study_access against a list of all users
		all_users = {}
		cursor.execute("SELECT username,user_id FROM users;")
		result = cursor.fetchall()
		for row in result:
			all_users[row[0]] = row[1]
		# Check that all users exist
		for username in study_access:
			if username:
				if username not in all_users.keys():
					flash("Error: User {} is not registered on Trips-Viz".format(username))
					return str(get_flashed_messages())
				else:
					cursor.execute("SELECT * FROM study_access WHERE user_id = {} and study_id = {};".format(all_users[username],study_id))
					result = cursor.fetchone()
					if result == None:
						cursor.execute("\n\n\nINSERT INTO study_access VALUES({},{});".format(study_id,all_users[username]))
		
		#Modify study names if they have changed
		new_study_name = "{}_{}".format(data[study_id][2],user_id)
		cursor.execute("SELECT study_name FROM studies WHERE study_id = {}".format(study_id))
		old_study_name = cursor.fetchone()[0]
		if old_study_name != new_study_name:
			#Update study name in the sqlite
			cursor.execute("UPDATE studies SET study_name = '{}' WHERE study_id = {}".format(new_study_name,study_id))
			#If the new_study_name folder does not exist, rename the old study to the new study, else move all files from old folder to new folder
			if not os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC,new_study_name)):
				os.rename("{0}/uploads/{1}".format(config.SCRIPT_LOC,old_study_name),"{0}/uploads/{1}".format(config.SCRIPT_LOC,new_study_name))
			else:
				if os.path.isdir("{}/uploads/{}".format(config.SCRIPT_LOC,old_study_name)):
					for filename in os.listdir("{}/uploads/{}".format(config.SCRIPT_LOC,new_study_name)):
						if os.path.isfile("{}/uploads/{}/{}".format(config.SCRIPT_LOC, old_study_name,filename)):
							os.rename("{}/uploads/{}/{}".format(config.SCRIPT_LOC, old_study_name,filename),"{}/uploads/{}/{}".format(config.SCRIPT_LOC, new_study_name,filename))
		# Change organism/transcriptome assembly if applicable
		organism_name = data[study_id][3]
		assembly_name = data[study_id][4]
		#print "organism name, assembly_name", organism_name, assembly_name
		cursor.execute("SELECT organism_id FROM studies WHERE study_id = {}".format(study_id))
		org_id = cursor.fetchone()[0]
		cursor.execute("SELECT organism_name,transcriptome_list FROM organisms WHERE organism_id = {}".format(org_id))
		result = cursor.fetchone()
		if result != None:
			old_organism = result[0]
			old_assembly = result[1]
			if old_organism != organism_name or old_assembly != assembly_name:
				print "{}: Changing from {} to {} or {} to {}".format(old_study_name,old_organism,organism_name, old_assembly, assembly_name)
				#Check if the new orgnaism and new assembly are a valid combination
				cursor.execute("SELECT organism_id  FROM organisms WHERE organism_name = '{}' AND transcriptome_list = '{}'".format(organism_name, assembly_name))
				result = cursor.fetchone()
				if result == None:
					return "Invalid organism/transcriptome combo for study {}".format(new_study_name)
				else:
					cursor.execute("UPDATE studies SET organism_id = {} WHERE study_id = {}".format(result[0],study_id))
					print "Org_id is {}".format(result[0])
					pass
					#update study_id with new org_id
				

	connection.commit()
	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")


# Allows users to delete transcriptomes
@app.route('/deletetranscriptomequery', methods=['GET','POST'])
@login_required
def deletetranscriptomequery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	for organism_id in data:
		transcriptomecheck = data[organism_id][0]
		if transcriptomecheck.split("_")[-1] != "undefined":
			organism_id = transcriptomecheck.split("_")[-1]
			#Delete the annotation sqlite file 
			cursor.execute("SELECT organism_name, transcriptome_list FROM organisms WHERE organism_id = {}".format(organism_id))
			result = cursor.fetchone()
			sqlite_path = "{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR, user_id, result[0], result[1] )
			if os.path.isfile(sqlite_path):
				os.remove(sqlite_path)
			sqlite_dir = "{0}transcriptomes/{1}/{2}/{3}".format(config.UPLOADS_DIR, user_id, result[0],result[1])
			if os.path.isdir(sqlite_dir):
				os.rmdir(sqlite_dir)
			cursor.execute("DELETE FROM organisms WHERE organism_id = {}".format(organism_id))
			#delete all files on the server associated with this organism, if there are any
			cursor.execute("SELECT * FROM files WHERE organism_id = {}".format(organism_id))
			result = cursor.fetchall()
			study_ids = []
			study_names = []
			if result != None:
				for row in result:
					file_id = row[0]
					filename = row[3]
					study_id = row[2]
					study_ids.append(study_id)
					cursor.execute("SELECT * FROM studies WHERE study_id = {}".format(study_id))
					result = cursor.fetchone()
					study_name = result[2]
					study_names.append(study_name)
					full_path = "{}{}/{}".format(config.UPLOADS_DIR,study_name, filename)
					if os.path.isfile(full_path):
						os.remove(full_path)
			
			for study_name in study_names:
					full_path = "{}{}".format(config.UPLOADS_DIR,study_name)
					if os.path.isdir(full_path):
						os.rmdir(full_path)
				
			#Now remove the study and the files associated with it from the db
			for study_id in study_ids:
				cursor.execute("DELETE FROM studies WHERE study_id = {}".format(study_id))
				cursor.execute("DELETE FROM files WHERE study_id = {}".format(study_id))
			connection.commit()
	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")

# Updates the "sequence rules", for custom sequence types
@app.route('/seqrulesquery', methods=['GET','POST'])
@login_required
def seqrulesquery():
	data = ast.literal_eval(request.data)
	user = current_user.name
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	#get user_id
	cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
	result = (cursor.fetchone())
	user_id = result[0]
	for seq_type in data:
		if data[seq_type][0] == 'False':
			cursor.execute("UPDATE seq_rules SET frame_breakdown = 0 WHERE seq_name = '{}' and user_id = {};".format(seq_type, user_id))
		elif data[seq_type][0] == 'True':
			cursor.execute("UPDATE seq_rules SET frame_breakdown = 1 WHERE seq_name = '{}' and user_id = {};".format(seq_type, user_id))
	connection.commit()
	connection.close()
	flash("Update successful")
	return redirect("https://trips.ucc.ie/uploads")



# breaks down the counts from each file for a specific ORF
@app.route('/<organism>/<transcriptome>/dataset_breakdown/')
def dataset_breakdown(organism,transcriptome):
	#ip = request.environ['REMOTE_ADDR']
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	
	accepted_orftype = request.args.get("region")
	transcript = request.args.get("transcript")
	start = int(request.args.get("start"))
	stop = int(request.args.get("stop"))
	files = request.args.get("files").strip(",")
	file_list = []
	for file_id in files.split("_"):
		file_list.append(int(file_id))
	file_paths_dict = fetch_file_paths(file_list,organism)
	xlist = []
	ylist = []
	filenames = []
	file_descs = []
	studies = []
	raw_reads = []
	controls = []
	cell_lines = []
	control_colors = []
	study_colors = []
	cell_line_colors = []
	i = 0
	cell_color_dict = {"HEK293":"#e00025","HeLa":"#150ed1","BJ fibroblast":"#0cfffa","fibroblast":"#00c6d1","MCF10A-ER-Src":"#00d100","HCT116":"#ffe900",
						"U2OS":"#ffc042"}
	orfquery_connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	orfquery_cursor = orfquery_connection.cursor()
	for file_id in file_paths_dict["riboseq"]:
		i += 1
		table_name = accepted_orftype
		xlist.append(i)
		sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
		try:
			raw_count = sqlite_db[table_name][transcript][stop][start]["profile_count"]
			raw_reads.append(raw_count)
		except:
			raw_count = 0
			raw_reads.append(raw_count)
		orfquery_cursor.execute("SELECT study_id,file_name,file_description from files WHERE file_id = {};".format(file_id))
		result = orfquery_cursor.fetchone()
		study_id = int(result[0])
		filenames.append(result[1].replace(".sqlite",""))
		file_descs.append(result[2])
		if True:
			controls.append("True")
			control_colors.append("#4ef46f")
		else:
			controls.append("False")
			control_colors.append("#000000")
		cell_lines.append("HEK293")
		if True:
			cell_line_colors.append(cell_color_dict["HEK293"])
		else:
			cell_line_colors.append("#bababa")
			
		ylist.append(raw_count/(float(10000)/100000))
		orfquery_cursor.execute("SELECT * from studies WHERE study_id = {};".format(study_id))
		result = orfquery_cursor.fetchone()
		studies.append(result[0])
		study_colors.append('#BABABA')
	orfquery_cursor.close()
	return riboflask_datasets.generate_plot(xlist,ylist,filenames,file_descs,studies,raw_reads,controls, cell_lines,control_colors, study_colors, cell_line_colors,transcript,start,stop)



# Estimates the time taken to complete the orfquery search
@app.route('/estimate_orfquery', methods=['POST'])
def estimate_orfquery():
	
	return "Estimated time: < 15 minutes"
	#return "Estimated time: 10 minutes"
	connection = sqlite3.connect("{}/trips.sqlite".format(config.SCRIPT_LOC))
	cursor = connection.cursor()
	acceptable = 0
	data = ast.literal_eval(request.data)
	html_args = data["html_args"]
	organism = data["organism"]
	file_paths_dict = fetch_file_paths(data["file_list"],organism)
	tran_list = data["tran_list"]
	accepted_orftype = data["region"]
	try:
		user = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
	except:
		user_id = None
	if html_args["user_short"] == "None":
		short_code = generate_short_code(data,organism,data["transcriptome"],"orf_translation")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True	
	transcriptome = data["transcriptome"]
	start_codons = []

	if "sc_aug" in data:
		start_codons.append("ATG")
	if "sc_cug" in data:
		start_codons.append("CTG")
	if "sc_gug" in data:
		start_codons.append("GTG")
	if "sc_none" in data:
		start_codons.append("None")
		
		
	min_start_increase = float(data["min_start_increase"])
	max_start_increase = float(data["max_start_increase"])

	min_stop_decrease = float(data["min_stop_decrease"])
	max_stop_decrease = float(data["max_stop_decrease"])

	min_cds_ratio = float(data["min_cds_ratio"])
	max_cds_ratio = float(data["max_cds_ratio"])
	
	min_coverage = float(data["min_coverage"])
	max_coverage = float(data["max_coverage"])

	min_lowest_frame_diff = float(data["min_lowest_frame_diff"])
	max_lowest_frame_diff = float(data["max_lowest_frame_diff"])

	min_highest_frame_diff = float(data["min_highest_frame_diff"])
	max_highest_frame_diff = float(data["max_highest_frame_diff"])

	min_cds = float(data["min_cds"])
	max_cds = float(data["max_cds"])

	min_len = float(data["min_len"])
	max_len = float(data["max_len"])

	min_avg = float(data["min_avg"])
	max_avg = float(data["max_avg"])

	filename = short_code+".csv"
	if os.path.isfile("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)):
		return "Loading previous result."

	if tran_list != "":
		tran_list  = tran_list.replace(","," ")
		for item in tran_list.split(" "):
			user_defined_transcripts.append(item)

	if accepted_orftype == "":
		return "Error no accepted_orftypes selected:"
	tran_gene_dict = {}
	# structure of orf dict is transcript[stop][start] = {"length":x,"score":0,"cds_cov":0} each stop can have multiple starts
	accepted_orf_dict = {}
	user_defined_transcripts = []
	tran_list = data["tran_list"]
	ambig = False
	if "ambig_check" in data:
		ambig = True

	filtered_transcripts = {}
	if "saved_check" in data:
		#if filter previously saved cases is turned on, then we query the sqlite database here and remove hits from transcript_list
		cursor.execute("SELECT tran,stop FROM users_saved_cases WHERE user_id = '{}' and organism = '{}';".format(user_id,organism))
		result = cursor.fetchall()
		
		for tran in result:
			if str(tran[0]) not in filtered_transcripts:
				filtered_transcripts[str(tran[0])] = []
			filtered_transcripts[str(tran[0])].append(int(tran[1]))

	if start_codons == []:
		return "Error no start codon types selected:"

	traninfo_connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips_annotations/{0}/{0}.v2.sqlite".format(organism))
	traninfo_cursor = traninfo_connection.cursor()
	
	principal_transcripts = []
	traninfo_cursor.execute("SELECT transcript,gene FROM transcripts WHERE principal = 1;")
	result = traninfo_cursor.fetchall()
	for row in result:
		principal_transcripts.append(str(row[0]))
		tran_gene_dict[row[0]] = row[1]
	
	if accepted_orftype != "drops":
		table_name = "{}".format(accepted_orftype)
	else:
		table_name = "cds"
	
	traninfo_cursor.execute("SELECT * FROM {} WHERE start_codon IN ({}) AND cds_coverage >= {} AND cds_coverage <= {} AND length >= {} AND length <= {} AND transcript IN ({});".format(table_name, str(start_codons).strip("[]"),min_cds,max_cds,min_len, max_len, str(principal_transcripts).strip("[]")))
	result = traninfo_cursor.fetchall()
	for row in result:
		if user_defined_transcripts != []:
			if row[0] not in user_defined_transcripts:
				continue
		if str(row[0]) not in accepted_orf_dict:
			accepted_orf_dict[str(row[0])] = {}
		if row[4] not in accepted_orf_dict[str(row[0])]:
			accepted_orf_dict[str(row[0])][row[4]] = {}
		accepted_orf_dict[str(row[0])][row[4]][row[3]] = {"length":row[2],
													  "score":0,
													  "cds_cov":row[6],
													  "start_codon":str(row[1])}

	# This will be populated with the users chosen file_ids and passed to the table, so that the trips link can use these files aswell.
	file_string = ""

	#Now build a profile for every transcript in accepted_transcripts
	master_profile_dict = {}

	# string based list of all acceptable transcripts, square brackets removed so that can be passed to sql
	transcript_list = str(accepted_orf_dict.keys()).strip("[]")

	#Mysql call will fail if we pass an empty list, so append NULL if thats the case
	if transcript_list == []:
		transcript_list.append("NULL")

	all_scores = []
	all_te = []
	all_start_increases = []
	all_stop_decreases = []
	all_cds_ratios = []
	all_results = []

	cds_average_dict = {}
	score_dict = {}

	# keeps track of the number of hits per gene
	gene_count_dict = {}
	missing_file_ids = []

	if file_paths_dict["rnaseq"] == {}:
		if "te_check" in data:
			del data["te_check"]
	
	if file_paths_dict["riboseq"] == {} and file_paths_dict["rnaseq"] == {}:
		flash("Error no files selected")
		return "Error no files selected"
	returnstr = ""
	
	all_values = []
	offset_dict = {}
	for file_id in file_paths_dict["riboseq"]:
		file_string += "{};".format(file_id)
		sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
		try:
			offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
			offset_dict[file_id] = offsets
		except:
			offset_dict[file_id] = {}
		sqlite_db.close()
	tran_count = 0
	best_high_frame = 0
	best_low_frame = 0
	best_start_score = 0
	best_stop_score = 0
	best_inframe_cov = 0
	tran_samples = 0
	start_time = time.time()
	
	for transcript in accepted_orf_dict:
		tran_samples += 1
		if tran_samples == 50:
			break
			
		gene = tran_gene_dict[transcript]

		tran_count += 1
		if tran_count%100 == 0:
			print "tran_count", tran_count
		profile = {}
		for file_id in file_paths_dict["riboseq"]:
			sqlite_db = SqliteDict(file_paths_dict["riboseq"][file_id])
			if transcript not in sqlite_db:
				continue
			offsets = offset_dict[file_id]
		
			subprofile = build_profile(sqlite_db[transcript], offsets,ambig)
			for pos in subprofile:
				if pos not in profile:
					profile[pos] = 0
				profile[pos] += subprofile[pos]
			
		for stop in accepted_orf_dict[transcript]:
			best_values = {"start":-1,"high_frame_count":-1,"low_frame_count":-1,"start_score":-1,"stop_score":-1,"final_score":-1,"coverage":0}
			for start in accepted_orf_dict[transcript][stop]:
				length = stop-start
				inframe_count = 0
				minframe_count = 0
				highframe_count = 0

				if_cl = []
				mo_count = 0
				po_count = 0
				if_cov = 0.0
				if_len = 0.0
				
				for x in range(start+8, stop-9,3):
					curr_mo = 0
					curr_po = 0
					if_len += 1
					if x-1 in profile:
						mo_count += profile[x-1]
						curr_mo = profile[x-1]
					if x+1 in profile:
						po_count += profile[x+1]
						curr_po = profile[x+1]
					if x in  profile:
						if_cl.append(profile[x])
						if "highest_frame_diff_check" in data:
							if profile[x] > max(curr_mo, curr_po):
								if_cov += 1
						else:
							if profile[x] > min(curr_mo, curr_po):
								if_cov += 1
				if_cov = if_cov/if_len
				
				#In frame count discards the highest peak
				inframe_count = sum(sorted(if_cl)[:-1])
					
				lowest_frame_count = inframe_count - min(mo_count,po_count)
				high_frame_count = inframe_count - max(mo_count, po_count)
				before_start = 0
				after_start = 0
				for y in range(start-7,start+7,3):
					if y in profile:
						if y < start:
							before_start += profile[y]
						else:
							after_start += profile[y]
				start_score = after_start-before_start
				
				before_stop = 0
				after_stop = 0
				for z in range(stop-9,stop+7,3):
					if z in profile:
						if z < stop:
							before_stop += profile[z]
						else:
							after_stop += profile[z]
				stop_score = before_stop-after_stop
				final_score_values = []
				if "start_increase_check" in data:
					final_score_values.append(start_score)
				if "stop_decrease_check" in data:
					final_score_values.append(stop_score)
				if "lowest_frame_diff_check" in data:
					final_score_values.append(lowest_frame_count)
				if "highest_frame_diff_check" in data:
					final_score_values.append(high_frame_count)
				if "coverage_check" in data:
					final_score_values.append(if_cov)
				final_score = sum(final_score_values)
				if final_score > best_values["final_score"] or (final_score == best_values["final_score"] and start > best_values["start"]):
					best_values["start"] = start
					best_values["high_frame_count"] = high_frame_count
					best_values["low_frame_count"] = lowest_frame_count
					best_values["start_score"] = start_score
					best_values["stop_score"] = stop_score
					best_values["final_score"] = final_score
					best_values["coverage"] = if_cov
					if high_frame_count > best_high_frame:
						best_high_frame = high_frame_count
					if lowest_frame_count > best_low_frame:
						best_low_frame = lowest_frame_count
					if stop_score > best_stop_score:
						best_stop_score = stop_score
					if start_score > best_start_score:
						best_start_score = start_score
					if if_cov > best_inframe_cov:
						best_inframe_cov = if_cov
			if best_values["final_score"] > 0:
				all_values.append([gene,
								transcript,
								best_values["start"],
								stop,
								length,
								best_values["high_frame_count"],
								best_values["low_frame_count"],
								best_values["stop_score"],
								best_values["start_score"],
								best_values["coverage"]])
			
	total_time = time.time()-start_time
	total_trans = len(accepted_orf_dict.keys())/50
	estimated_seconds = total_trans*total_time
	
	if estimated_seconds > 3600:
		return "Estimated time: {} hours".format(round((estimated_seconds/60)/60,1))
	elif estimated_seconds > 60:
		return "Estimated time: {} minutes".format(round(estimated_seconds/60,0))
	else:
		return "Estimated time: < 1 minute"







if __name__ == '__main__':
	local=False
	try:
		if sys.argv[1] == "true":
			local = True
	except:
		pass
	try:
		port_no = int(sys.argv[2])
	except:
		port_no = 5000
	if local == False:
		app.run(host='0.0.0.0',debug=False)
	else:
		app.run(host='0.0.0.0', port=port_no, debug=True)
