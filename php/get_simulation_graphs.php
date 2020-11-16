<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$dir = "./outputs/".$_REQUEST['folder']."/simulation_graphs";

$graphs = preg_grep('~\.(pdf)$~', scandir($dir));
$graphs = array_values($graphs);

for ($i = 0 ; $i < sizeof($graphs) ; $i ++) {

  $graphs[$i] = $dir."/".$graphs[$i];

}

echo(json_encode($graphs));


?>
