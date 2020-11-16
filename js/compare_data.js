//setInterval(refresh_annonces, 2000); // intervalle de temps pour rafraichir le contenu de la page
// fonction permettant d'ajouter des données et de vider le formulaire.

var folder_name = '' ;

var accept_analysis = true;

document.addEventListener('DOMContentLoaded', function() {
	if (document.getElementById('compare-data-form') != null) {
		// identification du formulaire (récupération sous forme d'une variable)
		var form = document.getElementById('compare-data-form');
		//fonction appelée lorsqu'on click sur le bouton de type submit
		form.addEventListener("submit", function(event) {

			event.preventDefault();
			// création du webservice
			var request = new XMLHttpRequest();

			request.addEventListener('load', function(data){
				folder_name = data.target.responseText;
				document.getElementById('comparison_type').value ='cvc';
				document.getElementById('dataset_1').value ='';
				document.getElementById('dataset_1_file').value ='';
				document.getElementById('dataset_2').value ='';
				document.getElementById('dataset_2_file').value ='';
				document.getElementById('genetic_code').value ='1';

				document.getElementById('content').innerHTML= ' <ul class="nav nav-pills nav-stacked left-menu" id="stacked-menu" id="test"> \
        <li> \
	          <a href="#" class = "btn btn-default" data-target="#cousin_log_system_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">See COUSIN system log<span class="caret arrow"></span></a> \
						<ul class="nav nav-pills nav-stacked collapse" id="cousin_log_system_container" style="overflow-x:auto;"> \
						</ul> \
	        </li> \
        <li><a class = "btn btn-success" href="./php/download_results.php?folder='+folder_name+'" target="_blank" style="font-size : large;">Download your results</a></li> \
      </ul>';
      			display_cousin_system_log(folder_name);
				//display_results(folder_name);
				//display_graphical_results(folder_name);
				//delete_folder(folder_name);
			});

			request.addEventListener('error', function(data) {
				console.log('error', data);

			});

			accept_analysis = true;

			var dataset_1_file_check = document.getElementById("dataset_1_file").value;
			var dataset_1_text_check = document.getElementById("dataset_1").value;


			if (dataset_1_file_check.match(/\S/)) {
				var dataset_1_file_check_size = document.getElementById("dataset_1_file").files[0].size;
			}

			else {
				var dataset_1_file_check_size = 0;
			}

			if ((dataset_1_file_check_size >= 1073741824)) {
				accept_analysis = false
				alert('Error, the size of your first dataset file is too big ! You should consider using the local version, reduce the size of your file, or reconsider your file.');
			}

			else if ((dataset_1_file_check.match(/\S/)) && (dataset_1_text_check.match(/\S/))) {
				alert('Warning, we found two kind of input for the first dataset : text and file. Only the data contained in the file is kept.');
			}

			else if (!(dataset_1_file_check.match(/\S/)) && !(dataset_1_text_check.match(/\S/))) {
				accept_analysis = false
				alert('There is no data for the first dataset. Please enter some data');
			}


			var dataset_2_file_check = document.getElementById("dataset_2_file").value;
			var dataset_2_text_check = document.getElementById("dataset_2").value;

			if (dataset_2_file_check.match(/\S/)) {
				var dataset_2_file_check_size = document.getElementById("dataset_2_file").files[0].size;
			}

			else {
				var dataset_2_file_check_size = 0;
			}

			if ((dataset_2_file_check_size >= 1073741824)) {
				accept_analysis = false
				alert('Error, the size of your second dataset file is too big ! You should consider using the local version, reduce the size of your file, or reconsider your file.');
			}


			else if ((dataset_2_file_check.match(/\S/)) && (dataset_2_text_check.match(/\S/))) {
				alert('Warning, we found two kind of input for the second dataset : text and file. Only the data contained in the file is kept.');
			}

			else if (!(dataset_2_file_check.match(/\S/)) && !(dataset_2_text_check.match(/\S/))) {
				accept_analysis = false
				alert('There is no data for the second dataset. Please enter some data');
			}


			if (accept_analysis == true) {
				document.getElementById('content').innerHTML='<center> \
				<i class="fa fa-cog fa-spin fa-5x fa-fw"></i> \
				<span class="sr-only">Loading...</span> \
				<h2> Loading </h2> \
				</center>';
				request.open("POST", "./php/compare_data_step.php", true); // ouverture du code ajout de données php
				request.send(new FormData(form)); // envoi des données au code ajout php

			}
			// display_results();
			// display_graphical_results();
		});
	}
});

// fonction permettant de récuperer le message, de créer le contenu de la variable ret
function display_results(folder) {
	var request = new XMLHttpRequest();

	request.addEventListener('load', function(data) {
		//console.log(recherche);
		var table_line = data.target.responseText;
		table_line = table_line.split('\n');
		var new_html = ' <font size="2"><table class="table table-bordered table-inverse table-responsive">';

		if (table_line.length < 150) {
			for (var i=0; i < table_line.length; i++) {
				new_html += build_results_html(table_line[i]);
			}
			new_html += '</tbody></table></font>'
		}
		else {
			for (var i=0; i < 150; i++) {
				new_html += build_results_html(table_line[i]);
			}
			new_html += '</tbody></table></font>'
			new_html += '<center> <i class="fa fa-ellipsis-v fa-4x"></i> </center>'
		}
		

		document.querySelector('#results_container').innerHTML = new_html;


	});

	request.open("POST", "./php/get_results.php",true );
	request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	request.send("folder=" + folder);

}

function build_results_html(line) {
	var table_element = line.split('\t');
	if (table_element[0] == 'Header') {
		html_table_line = '<thead><tr>';
		for (var i=0; i < table_element.length; i++) {
			html_table_line += '<th class="text-center">' + table_element[i].substring(0,50) + '</th>';


		}
		html_table_line += '</tr></thead>';
	}
	else {
		html_table_line = '<tbody><tr>';
		for (var i=0; i < table_element.length; i++) {
			if (i == 0) {
				html_table_line += '<td class="text-center" style ="font-weight: bold";>' + table_element[i].substring(0,50) + '</td>';
			}
			else {
				html_table_line += '<td class="text-center">' + table_element[i].substring(0,50) + '</td>';
			}
		}
		html_table_line += '</tr>';
	}

	return html_table_line;
}


function delete_folder(folder) {

	request = new XMLHttpRequest();

	request.addEventListener('error', function(data) {
		console.log('error', data);

	});

	request.open("POST", "./php/delete_dir.php",true );
	request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	request.send("folder=" + folder);

}


function display_basic_simulation_graphical_results(folder) {
	var request = new XMLHttpRequest();
	var new_html = "";

	request.addEventListener('load', function(data) {
		var array_files = JSON.parse(data.target.responseText);
		if (array_files.length < 20) {
			for(var i=0;i<array_files.length;i++){
				new_html += 	'<center><embed src="./php/' + array_files[i] + '" width="700px" height="550px"> </embed></center> <br>'
			}
		
		}
		else {
			for(var i=0;i<20;i++){
				new_html += 	'<center><embed src="./php/' + array_files[i] + '" width="700px" height="550px"> </embed></center> <br>'
			}
			new_html += '<center> <i class="fa fa-ellipsis-v fa-4x"></i> </center>'
		} 
			

		document.querySelector('#basic_simulation_graphical_results_container').innerHTML = new_html;

	});

	request.open("POST", "./php/get_basic_simulation_graphs.php",true );
	request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	request.send("folder=" + folder);

}

function display_cousin_system_log(folder) {
var request = new XMLHttpRequest();
	request.addEventListener('load', function(data) {
		//console.log(recherche);
		var table_line = data.target.responseText;
		table_line = table_line.split("\n");
		var new_html = ' <center> <p> <b> COUSIN LOG SYSTEM </b> </p> </center>';
		for (var i=0; i < table_line.length; i++) {
			if (i%4 == 0 && i != 0) {
				new_html += "<br> <span style='white-space: pre-wrap; font-family: monospace; line-height: 0.01;'>" + table_line[i] + "</span> <br>";
			}
			else {
				new_html +="<span style='white-space: pre-wrap; font-family: monospace; line-height: 0.01;'>" + table_line[i] + "</span> <br>";
			}
		}
		document.querySelector('#cousin_log_system_container').innerHTML = new_html;
	});

	request.open("POST", "./php/get_cousin_system_log.php",true );
	request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	request.send("folder=" + folder);
}
