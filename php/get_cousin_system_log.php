<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$cousin_logs_system = file_get_contents("./cousin_logs_system/".$_REQUEST['folder'].".txt");

echo ($cousin_logs_system);

?>
