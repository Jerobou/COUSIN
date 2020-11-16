//setInterval(refresh_annonces, 2000); // intervalle de temps pour rafraichir le contenu de la page
// fonction permettant d'ajouter des données et de vider le formulaire.

var folder_name = '' ;

var accept_analysis = true;

document.addEventListener('DOMContentLoaded', function() {

	var example_action = document.getElementById("example_button");

	example_action.addEventListener("click", function(event) {

		document.getElementById('fasta_sequences').value ='';
		document.getElementById('fasta_file').value ='';
		document.getElementById('genetic_code').value='1';
		document.getElementById('optimal_codons_text').value='';
		document.getElementById('optimal_codons_file').value='';
		document.getElementById('simulation_routine').checked = false;
		document.getElementById('specific_COUSIN').checked = false;
		document.getElementById('vector_analysis').checked = false;
		event.preventDefault();
		var request = new XMLHttpRequest();

		request.addEventListener('load', function(data){
			document.getElementById('fasta_sequences').value= data.target.responseText;
		});

		request.open("POST", "./php/give_example_query.php", true);
		request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
		request.send();


		document.getElementById('genetic_code').value='1';


		window.alert("An example has been provided with ten Homo sapiens CDSs (query). Press the 'Do the analysis' button to launch your example study ! ")

	});
	if (document.getElementById('create-table-form') != null) {
		// identification du formulaire (récupération sous forme d'une variable)
		var form = document.getElementById('create-table-form');
		//fonction appelée lorsqu'on click sur le bouton de type submit
		form.addEventListener("submit", function(event) {

			event.preventDefault();
			// création du webservice
			var request = new XMLHttpRequest();

			request.addEventListener('load', function(data){

				simulation_routine_checked = document.getElementById('simulation_routine').checked;

				folder_name = data.target.responseText;
				//window.alert(data.target.responseText);
				document.getElementById('fasta_sequences').value ='';
				document.getElementById('fasta_file').value ='';
				document.getElementById('genetic_code').value='1';
				document.getElementById('optimal_codons_text').value='';
				document.getElementById('optimal_codons_file').value='';
				document.getElementById('simulation_routine').checked = false;
				document.getElementById('specific_COUSIN').checked = false;
				document.getElementById('vector_analysis').checked = false;


				if (simulation_routine_checked == true) {
					document.getElementById('content').innerHTML= ' <ul class="nav nav-pills nav-stacked left-menu" id="stacked-menu" id="test"> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#cu_table_results_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">Display your Codon Usage Table<span class="caret arrow"></span></a> \
						<ul class="nav nav-pills nav-stacked collapse" id="cu_table_results_container" style="overflow-x:auto;"> \
						</ul> \
	        </li> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#results_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">Display your results<span class="caret arrow"></span></a> \
	          <ul class="nav nav-pills nav-stacked collapse" id="results_container" style="overflow-x:auto;"> \
	          </ul> \
	        </li> \
					<li> \
	          <a href="#" class = "btn btn-default" data-target="#basic_simulation_graphical_results_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">Display your graphical results (basic simulation)<span class="caret arrow"></span></a> \
	          <ul class="nav nav-pills nav-stacked collapse" id="basic_simulation_graphical_results_container" style="overflow-x:auto;"> \
	          </ul> \
	        </li> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#cousin_log_system_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">See COUSIN system log<span class="caret arrow"></span></a> \
						<ul class="nav nav-pills nav-stacked collapse" id="cousin_log_system_container" style="overflow-x:auto;"> \
						</ul> \
	        </li> \
	        <li><a href="./php/download_results.php?folder='+folder_name+'" class = "btn btn-success" target="_blank" style="font-size : large;">Download your results</a></li> \
	      </ul>';

				display_cu_table(folder_name);
				display_results(folder_name);
				display_basic_simulation_graphical_results(folder_name);
				display_cousin_system_log(folder_name);
					//delete_folder(folder_name);

				}

				else {

					document.getElementById('content').innerHTML= ' <ul class="nav nav-pills nav-stacked left-menu" id="stacked-menu" id="test"> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#cu_table_results_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">Display your Codon Usage Table<span class="caret arrow"></span></a> \
						<ul class="nav nav-pills nav-stacked collapse" id="cu_table_results_container" style="overflow-x:auto;"> \
						</ul> \
	        </li> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#results_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">Display your results<span class="caret arrow"></span></a> \
	          <ul class="nav nav-pills nav-stacked collapse" id="results_container" style="overflow-x:auto;"> \
	          </ul> \
	        </li> \
	        <li> \
	          <a href="#" class = "btn btn-default" data-target="#cousin_log_system_container" data-toggle="collapse" data-parent="#stacked-menu" style="font-size : large;">See COUSIN system log<span class="caret arrow"></span></a> \
						<ul class="nav nav-pills nav-stacked collapse" id="cousin_log_system_container" style="overflow-x:auto;"> \
						</ul> \
	        </li> \
	        <li><a href="./php/download_results.php?folder='+folder_name+'" class = "btn btn-success" target="_blank" style="font-size : large;">Download your results</a></li> \
	      </ul>';

				display_cu_table(folder_name);
				display_results(folder_name);
				display_basic_simulation_graphical_results(folder_name);
				display_cousin_system_log(folder_name);

				}
			});

			request.addEventListener('error', function(data) {
				console.log('error', data);

			});

			accept_analysis = true;

			var fasta_file_check = document.getElementById("fasta_file").value;	
			var fasta_sequences_text_check = document.getElementById("fasta_sequences").value;

			if (fasta_file_check.match(/\S/)) {
				var fasta_file_check_size = document.getElementById("fasta_file").files[0].size;
			}

			else {
				var fasta_file_check_size = 0;
			}

			if ((fasta_file_check_size >= 1073741824)) {
				accept_analysis = false
				alert('Error, the size of your fasta file is too big ! You should consider using the local version, or reduce the size of your file !');
			}

			else if ((fasta_file_check.match(/\S/)) && (fasta_sequences_text_check.match(/\S/))) {
				alert('Warning, we found two kind of input for the fasta sequences (query) : text and file. Only the data contained in the file is kept.');
			}

			else if (!(fasta_file_check.match(/\S/)) && !(fasta_sequences_text_check.match(/\S/))) {
				accept_analysis = false
				alert('There is no data for the fasta sequences (query). Please enter some data');
			}

			var optimal_codons_file_check = document.getElementById("optimal_codons_file").value;
			var optimal_codons_text_check = document.getElementById("optimal_codons_text").value;

			if (optimal_codons_file_check.match(/\S/)) {
				var optimal_codons_file_check_size = document.getElementById("optimal_codons_file").files[0].size;
			}

			else {
				var optimal_codons_file_check_size = 0;
			}

			if ((optimal_codons_file_check_size >= 1073741824)) {
				accept_analysis = false
				alert('Error, the size of your "optimal codons" file is, strikingly, too big ! You should consider using the local version, reduce the size of your file, or reconsider your optimal codons file.');
			}

			else if ((optimal_codons_file_check.match(/\S/)) && (optimal_codons_text_check.match(/\S/))) {
				alert('Warning, we found two kind of input for the optimal codons data : text and file. Only the data contained in the file is kept.');
			}

			else if (!(optimal_codons_file_check.match(/\S/)) && !(optimal_codons_text_check.match(/\S/))) {
				alert('There is no data for the optimal codons. A dataset of these codons will be constructed from the reference data.');
			}

			if (accept_analysis == true) {
				document.getElementById('content').innerHTML='<center> \
				<i class="fa fa-cog fa-spin fa-5x fa-fw"></i> \
				<span class="sr-only">Loading...</span> \
				<h2> Loading </h2> \
				</center>';
				request.open("POST", "./php/create_table.php", true); // ouverture du code ajout de données php
				request.send(new FormData(form)); // envoi des données au code ajout php

			}

			// display_results();
			// display_basic_simulation_graphical_results();
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
		var new_html = ' <br> <font size="1.5"><table class="table table-bordered table-inverse table-responsive">';

		if (table_line.length < 150) {
			for (var i=0; i < table_line.length; i++) {
				if (table_line[i] != "") {
					new_html += build_results_html(table_line[i]);	
				}
			}
			new_html += '</tbody></table></font>'
		}
		else {
			for (var i=0; i < 150; i++) {
				new_html += build_results_html(table_line[i]);
			}
			new_html += '</tbody></table></font> <br>'
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
			if (i == 0) {
				html_table_line += '<th class="text-center" style="width:150px;">' + table_element[i] + '</th>';
			}
			else {
				html_table_line += '<th class="text-center" style="width:75px;">' + table_element[i] + '</th>';
			}
		}
		html_table_line += '</tr></thead>';
	}
	else {
		html_table_line = '<tbody><tr>';
		for (var i=0; i < table_element.length; i++) {
			if (i == 0) {
				html_table_line += '<td class="text-center" style ="font-weight: bold;">' + table_element[i] + '</td>';
			}
			else {
				html_table_line += '<td class="text-center">' + table_element[i] + '</td>';
			}
		}
		html_table_line += '</tr>';
	}

	return html_table_line;
}


function display_cu_table(folder) {
var request = new XMLHttpRequest();

	request.addEventListener('load', function(data) {
		//console.log(recherche);
		var table_line = data.target.responseText;
		table_line = table_line.split("\n");
		var new_html = ' <center> <p> <b> CODON USAGE TABLE [triplet] [frequency: per thousand] ([number]) </b> </p>';
		for (var i=0; i < table_line.length; i++) {
			if (i%4 == 0 && i != 0) {
				new_html += "<br> <span style='white-space: pre-wrap; font-family: monospace; line-height: 0.01;'>" + table_line[i] + "</span> <br>";
			}

			else {
				new_html +="<span style='white-space: pre-wrap; font-family: monospace; line-height: 0.01;'>" + table_line[i] + "</span> <br>";
			}
		}
		new_html += '</center>'
		document.querySelector('#cu_table_results_container').innerHTML = new_html;
	});

	request.open("POST", "./php/get_cu_table.php",true );
	request.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	request.send("folder=" + folder);

}

function build_cu_table_html(line) {
	console.log(line);
	var regexp_cu_table = /([A-z]{3}[^\(]+\({1}\s*[0-9]*\){1})/g;
	var match = regexp_cu_table.exec(line);
	console.log(match);
	// var table_element = line.split(regexp_cu_table[1]);
	// if (table_element[0] == 'Header') {
	// 	html_table_line = '<thead><tr>';
	// 	for (var i=0; i < table_element.length; i++) {
	// 		html_table_line += '<th>' + table_element[i].substring(0,50) + '</th>';
	//
	//
	// 	}
	// 	html_table_line += '</tr></thead>';
	// }
	// else {
	// 	html_table_line = '<tbody><tr>';
	// 	for (var i=0; i < table_element.length; i++) {
	// 		html_table_line += '<td>' + table_element[i].substring(0,50) + '</td>';
	// 	}
	// 	html_table_line += '</tr>';
	// }
	//
	// return html_table_line;
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
