<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);
// Définit le fuseau horaire par défaut à utiliser. Disponible depuis PHP 5.1
date_default_timezone_set('Europe/Paris');

putenv("PATH=/usr/local/bin/:" . exec('echo $PATH'));

if (isset($_REQUEST['first_name'])) {
  $first_name = $_REQUEST['first_name'];
}

if (isset($_REQUEST['last_name'])) {
  $last_name = $_REQUEST['last_name'];
}

if (isset($_REQUEST['status'])) {
  $status = $_REQUEST['status'];
}

if (isset($_REQUEST['email_adress'])) {
  $email_adress = $_REQUEST['email_adress'];
}

$infos = $first_name."\t".$last_name."\t".$status."\t".$email_adress.PHP_EOL;

$downloaders_list = fopen('./downloaders_list/downloaders_list.tsv', "a");   
fwrite($downloaders_list, $infos);
fclose($downloaders_list);

echo(TRUE);

?>
