//setInterval(refresh_annonces, 2000); // intervalle de temps pour rafraichir le contenu de la page
// fonction permettant d'ajouter des données et de vider le formulaire.

var accept_download = 0 ;

document.addEventListener('DOMContentLoaded', function() {

	if (document.getElementById('local-version-form') != null) {
		// identification du formulaire (récupération sous forme d'une variable)
		var form = document.getElementById('local-version-form');
		//fonction appelée lorsqu'on click sur le bouton de type submit
		form.addEventListener("submit", function(event) {

			event.preventDefault();
			// création du webservice
			var request = new XMLHttpRequest();

			request.addEventListener('load', function(data){

				document.getElementById('first_name').value ='';
				document.getElementById('last_name').value ='';
				document.getElementById('status').value='Master';
				document.getElementById('email_adress').value ='';

				accept_download = data.target.responseText;

				if (accept_download == 1) {

					location.href = "./php/download_version.php";


				}
				else {

					window.alert("Error, something went wrong and you can't download COUSIN right now. Try again or contat us (jerome.bourret@ird.fr)");
				}

			});	

			request.addEventListener('error', function(data) {
				console.log('error', data);

			});

			request.open("POST", "./php/save_downloaders_infos.php", true); // ouverture du code ajout de données php
			request.send(new FormData(form)); // envoi des données au code ajout php

		});
	}
});
