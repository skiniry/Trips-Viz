function changebuttonlink(gwips_info){
		gwips_clade = gwips_info["database"];
		gwips_org = gwips_info["organism"];
		gwips_db = gwips_info["clade"];
		var userInput = document.getElementById('tran').value;
		var lnk = document.getElementById('gwips_link');
		console.log("GWIPS CLADE IS"+gwips_clade);
		console.log("GWIPS org IS"+gwips_org);
		console.log("GWIPS db IS"+gwips_db);
		lnk.href = "https://gwips.ucc.ie/cgi-bin/hgTracks?clade="+gwips_clade+"&org="+gwips_org+"&db="+gwips_db+"&position=" + userInput + "&Submit=submit&hgsid=39696_i99V6EU3BDqA1tdBqCHlfiiNDpJN&pix=1835";
};

function addRow(rownum, innerHTML) {
		console.log("Row num is"+rownum);
		var row = document.getElementById("tablerow"+rownum);
		var x = row.insertCell(0);
		x.innerHTML = innerHTML;
}

function addFileRow(rownum, innerHTML, studyname) {
		console.log(studyname+"tablerow"+rownum);
		var row = document.getElementById(studyname+"tablerow"+rownum);
		var x = row.insertCell(0);
		x.innerHTML = innerHTML;
}

function populate_jquery_dialog(splitdata) {
		var dialog_div = document.getElementById('dialog-2');
		dialog_div.innerHTML = "";

		tablehtml = "<table id='myTable' class='tablesorter tran_select_box hover-highlight'>  <thead> <tr><th><b>Transcript</b></th> <th><b>Version</b></th>   <th><b>Length</b></th>    <th><b>5' Length</b></th>    <th><b>Cds length</b></th>    <th><b>3' Length</b></th>  <th>Type</th>  </tr></thead>"
		for(var i = 1; i < splitdata.length; i++) {
				var sublist = splitdata[i].split(",");
				tablehtml += "<tr>"
				// for item in sublist add an element to the table
				for(var x = 0; x < sublist.length; x++) {
						str_input = sublist[x];
						//console.log("STR INPUT IS "+str_input);
						if (x == 0) {
								tablehtml += "<td><input type='radio' id='transcript_sel' value='"+str_input+"' name='transcript_sel' checked></input>"+str_input+"</td>"
						}
						else {
								tablehtml += "<td>"+str_input+"</td>"
						}
				}
				//finish this row
				tablehtml += "</tr>";
		}
		//finish the table
		tablehtml += "</table>";
		dialog_div.innerHTML += tablehtml;
}





//When the view study button is clicked get the currently selected study and display info in a popup
function view_study_info(){
		var study_dialog = document.getElementById('study_dialog');
		study_dialog.innerHTML = "";
		var selected_study = document.querySelector('input[name=content_type]:checked').value;
		console.log("This is what a reals tudy looks like"+selected_study);
		tablehtml = "<table class='study_dialog_table'><thead><tr><th>Study name</th><th>GSE</th><th>Title</th><th>Year</th><th>Paper</th><th>SRP</th><th>PMID</th><th>Adapters</th></tr></thead>";

		for (var key in studies_dict)  {
			if (($('.rnaseq_study'+key+':checked').length == 0) && ($('.riboseq_study'+key+':checked').length == 0)){
				//document.getElementById('label'+key+'').style.borderColor='#e2e2e2';
			}
			else {

				var study = "study"+key;
				var study_name = studyinfo_dict[study]["study_name"];
				var gse = studyinfo_dict[study]["gse_nos"];
				var title = studyinfo_dict[study]["paper_title"];
				var year = studyinfo_dict[study]["paper_year"];
				var paper_link = studyinfo_dict[study]["paper_link"];
				var srp = studyinfo_dict[study]["srp_nos"];
				var pmid = studyinfo_dict[study]["paper_pmid"];
				var adapters = studyinfo_dict[study]["adapters"];

				// Highlight the currently selected study in bold
				if (selected_study == study) {
					tablehtml += "<tr style='font-weight:bold'>";
				}
				else {
						tablehtml += "<tr>";
				}
				tablehtml += "<td>"+study_name+"</td>";
				tablehtml += "<td>  <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+gse+"'   target='_blank'>"+gse+"</a> </td>";
				tablehtml += "<td>"+title+"</td>";
				tablehtml += "<td>"+year+"</td>";
				tablehtml += "<td><a href='"+paper_link+"'   target='_blank'>"+paper_link+"</a></td>";
				tablehtml += "<td><a href='https://www.ncbi.nlm.nih.gov/Traces/study/?acc="+srp+"&go=go'    target='_blank'>"+srp+"</a></td>";
				tablehtml += "<td>"+pmid+"</td>";
				tablehtml += "<td>"+adapters+"</td></tr>";
				if (selected_study == study) {
					tablehtml += "</b>"
				}

			}
		};

		tablehtml += "</table>";
		//var study = document.querySelector('input[name=content_type]:checked').value;

		study_dialog.innerHTML += tablehtml;

		$( "#study_dialog" ).dialog( "open" );


		};




/*
		// First find what study is currently checked
		var study = document.querySelector('input[name=content_type]:checked').value;
		//console.log("Study is "+study);
		//console.log("studyinfo_dict "+ studyinfo_dict);

		// get the relevant info from the study info dict
		var year = studyinfo_dict[study]["paper_year"]
		var gse = studyinfo_dict[study]["gse_nos"]
		console.log("GSE "+gse)
		var title = studyinfo_dict[study]["paper_title"]
		var paper_link = studyinfo_dict[study]["paper_link"]
		var srp = studyinfo_dict[study]["srp_nos"]
		var desc = studyinfo_dict[study]["description"]
		var adapters = studyinfo_dict[study]["adapters"]
		//console.log("year "+year);

		// reset the html in study dialog and populate with the above info
		var study_dialog = document.getElementById('study_dialog');
		study_dialog.innerHTML = "";
		tablehtml = "<table class='study_dialog_table'>"
		tablehtml += "<tr><td><b>GSE</b></td><td>  <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+gse+"'   target='_blank'>"+gse+"</a> </td></tr>"
		tablehtml += "<tr><td><b>Title</b></td><td>"+title+"</td></tr>"
		tablehtml += "<tr><td><b>Year</b></td><td>"+year+"</td></tr>"
		tablehtml += "<tr><td><b>Paper</b></td><td><a href='"+paper_link+"'   target='_blank'>"+paper_link+"</a></td></tr>"
		tablehtml += "<tr><td><b>Srp</b></td><td><a href='https://www.ncbi.nlm.nih.gov/Traces/study/?acc="+srp+"&go=go'    target='_blank'>"+srp+"</a></td></tr>"
		tablehtml += "<tr><td><b>Description</b></td><td>"+desc+"</td></tr>"
		tablehtml += "<tr><td><b>Adapters</b></td><td>"+adapters+"</td></tr>"
		tablehtml += "</table>"
		study_dialog.innerHTML += tablehtml;

		//$( "#dialog-2" ).dialog( "open" );
		$( "#study_dialog" ).dialog( "open" );
}

*/


$("#UncheckAll").click(function(inclass){
		$("input.inclass[type='checkbox']").prop('checked',false);
});


$("#CheckAll").click(function(inclass){
		$("input.inclass[type='checkbox']").prop('checked',true);
})

// Unchecks every checkbox that has the class <seq_type>_file
function uncheck_all(){
		var seq_type = document.querySelector('input[name="seq_type"]:checked').value;
		$('input.'+seq_type+'_file').prop('checked',false);
		update_colours();
};

// Checks every checkbox that has the class <seq_type>_file
function check_all(){
		var seq_type = document.querySelector('input[name="seq_type"]:checked').value;
		$('input.'+seq_type+'_file').prop('checked',true);
		update_colours();
};


// checks all checkboxes with the <seq_type>_study<study_id> class, the study_id is gotten from the value of the currently selected study
function check_study(){
		var curstudy = $('input[name="content_type"]:checked').val();
		var seq_type = document.querySelector('input[name="seq_type"]:checked').value;
		$('input.'+seq_type+'_'+curstudy+'').prop('checked',true);
		update_colours();
};

function uncheck_study(){
		var curstudy = $('input[name="content_type"]:checked').val();
		var seq_type = document.querySelector('input[name="seq_type"]:checked').value;
		$('input.'+seq_type+'_'+curstudy+'').prop('checked',false);
		update_colours();
};




$(function() {
		$( "#study_dialog" ).dialog({
				autoOpen: false,
				modal: true,
				minWidth: 1700,
				title: "Study info",
				draggable: false,
				position: {
						my: "center",
						at: "center"
				}
		});
});



// This is needed to show/hide checkboxes when each studies radio button is selected/deselected, any files with a class of <seq_type><study_id>
$(document).ready(function () {
		$(".show").hide();
		var show_sel = $(".show")
		$('input[name=content_type]').on('change', function() {
				var n = $(this).val();
				var seq_type = document.querySelector('input[name="seq_type"]:checked').value;
				var selector = "#"+seq_type+ n; //Create a selector based on the id
				$(".show").hide(); //Hide all the checkbox containers
				$(selector).show(); //Show only the related container
		});
});


// This is needed to show/hide checkboxes when the seq_type button is selected/deselected, any files with a class of <seq_type><study_id>
$(document).ready(function () {
		$(".show").hide();
		$('input[name=seq_type]').on('change', function() {
				console.log("Seq TypE WaS CHAngEd");
				var n = document.querySelector('input[name="content_type"]:checked').value;
				var seq_type = document.querySelector('input[name="seq_type"]:checked').value;

				for (var study_id in studies_dict) {
					$("#label"+study_id+"").css({ opacity: 0.2 });
				}
				for (file_type in seq_study_dict) {
					if (file_type == seq_type) {
						console.log("File type is seq_type "+file_type+"  "+seq_type);
						for (var sub_study_id in seq_study_dict[file_type]) {
							$("#label"+seq_study_dict[file_type][sub_study_id]+"").css({ opacity: 1 });
						}
					}
				}

				var selector = "#"+seq_type+ n; //Create a selector based on the id
				$(".show").hide(); //Hide all the checkbox containers
				$(selector).show(); //Show only the related container
		});
});
