<?php
session_start();

?>

<!doctype html>

<html>

  <!-- Ici est contenu l'ensemble de la partie "haute" du site. Soit la barre de navigation
        Latérale ainsi que celle horizontale-->
  <head>
    <link rel="icon" type="image/png" href="./images/COUSIN_FAVICON.png" />
    <meta charset="utf-8">
    <title>COUSIN analysis tool</title>

    <!-- external libs -->
	  <link href="css/bootstrap.min.css" rel="stylesheet" crossorigin="anonymous">
    <link rel="stylesheet" href="css/font-awesome/css/font-awesome.min.css">

    <!-- my libs -->
    <link rel="stylesheet" href="css/cousin.css">
    <script type="text/javascript" src="js/jquery.js"> </script>

  </head>
  <body>
	  <div class="container-title">

    <!-- Ici, création d'une navbar horizontale spécifique :
      Premier cas : si l'utilisateur est connecté, affiche son login et la possibilité de se déconnecter
      Second cas : si l'utilisateur n'est pas connecté, il peut se connecter (si il est enregistré au préalable)-->

      <div class="navbar navbar-inverse">
        <div class="container-fluid" >
          <div class="navbar-header">
            <a class="navbar-brand" id="menu-toggle" href=#><i id="menu_toggle_icon" class="fa fa-minus-square fa-2x"></i></a>
            <script>
            $('#menu-toggle').click(function(){
              $(this).find('i').toggleClass('fa-plus-square fa-2x fa-minus-square fa-2x')
            });
            </script>
            <a class="navbar-brand" href="./index.php"><img id="logo_image" src="images/COUSIN_LOGO.png" alt="Dispute Bills"> <esp_hor> COUSIN (COdon Usage Similarity INdex) </esp_hor></a>
          </div>
          <div id="navbar-titre" class="navbar-collapse collapse">
            <div class="nav navbar-nav navbar-right">
              <a class="navbar-brand" href="https://github.com/Jerobou/cousin" target="_blank"><img id="github_image" src="images/GITHUB_LOGO.png" alt="Dispute Bills"> <esp_hor> GitHub </esp_hor></a>
            </div>
  		</div>
		</div>
    </div>
  </div>

    <!--Barre latérale : -->

    <div id="wrapper">
        <div id="sidebar-wrapper">
            <ul class="sidebar-nav" id="test">
              <li class="sidebar-brand"><a href="index.php"><i class="fa fa-home"></i></i> Homepage</a></li>
              <br>
              <li id="bullet_triangle"><a href="./calculation.php"> <i class="fa fa-caret-right"></i> Basic calculation step</a></li>
              <li id="bullet_triangle"></i><a href="./pattern_analysis.php"> <i class="fa fa-caret-right"></i> Pattern Analysis</a></li>
              <li id="bullet_triangle"><a href="./optimization.php"> <i class="fa fa-caret-right"></i> Optimization Analysis</a></li>
              <li id="bullet_triangle"><a href="./simulation_analysis.php"> <i class="fa fa-caret-right"></i> Simulation Analysis</a></li>
              <li id="bullet_triangle"><a href="./clustering_analysis.php"> <i class="fa fa-caret-right"></i> Clustering Analysis</a></li>
              <li id="bullet_triangle"><a href="./create_table.php"> <i class="fa fa-caret-right"></i> Create a Codon Usage Table</a></li>
              <li id="bullet_triangle"><a href="./compare_data.php"> <i class="fa fa-caret-right"></i> Compare two datasets</a></li>
              <br>
              <li id="bullet_triangle"><a href="./tutorial.php"> <i class="fa fa-info-circle"></i> Tutorial</a></li>
              <li id="bullet_triangle"><a href="./local_version.php"> <i class="fa fa-download"></i> Local version</a></li>
              <li id="bullet_triangle"><a href="./faq.php"> <i class="fa fa-question-circle"></i> FAQ </a></li>
              <li id="bullet_triangle"><a href="./development.php"> <i class="fa fa-flask"></i> Development </a></li>
              <!--<li id="bullet_triangle"><a href="./beta_user_guide.php"> <i class="fa fa-user"></i> Beta-users guide </a></li> -->

              <center>
              <a href="https://cordis.europa.eu/project/rcn/199662/factsheet/en" target="blank"> <img class = "img-logo-sidebar" src="./images/logo_codovirevol.png"> </img> </a>
              </center>

            </ul>
        </div>
