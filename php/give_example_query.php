<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$example = file_get_contents("../H_sapiens_example/h_sapiens_query_example.fasta");
echo ($example);


?>
