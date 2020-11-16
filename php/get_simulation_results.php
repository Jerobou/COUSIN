<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$results = file_get_contents("./outputs/".$_REQUEST['folder']."/simulation_results/simulation_results.tsv");
echo ($results);



?>
