<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$results = file_get_contents("./outputs/".$_REQUEST['folder']."/results/results.tsv");
echo ($results);



?>
