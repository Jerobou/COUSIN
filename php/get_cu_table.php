<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$cu_table = file_get_contents("./outputs/".$_REQUEST['folder']."/codon_usage_tables/cu_table.txt");

echo ($cu_table);



?>
