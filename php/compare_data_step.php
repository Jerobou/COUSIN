<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);
// Définit le fuseau horaire par défaut à utiliser. Disponible depuis PHP 5.1
date_default_timezone_set('Europe/Paris');

putenv("PATH=/usr/local/bin/:" . exec('echo $PATH'));

function delete_dir($dir) {
   $files = array_diff(scandir($dir), array('.','..'));
    foreach ($files as $file) {
      (is_dir("$dir/$file")) ? delTree("$dir/$file") : unlink("$dir/$file");
    }
    return rmdir($dir);
  }

$folder_input_name = uniqid(date('y_m_d_H_i_s_'));

mkdir('./uploads/'.$folder_input_name."/");


$name       = $_FILES['dataset_1_file']['name'];
$temp_name  = $_FILES['dataset_1_file']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if(move_uploaded_file($temp_name, $location.$name)){
         $dataset_1_file = "./uploads/".$folder_input_name."/".$_FILES['dataset_1_file']['name'];
       }
   }
else {
  if ($_REQUEST['dataset_1'] != '') {
    $input_file = fopen('./uploads/'.$folder_input_name."/dataset_1.txt", "w");
    fwrite($input_file, $_REQUEST['dataset_1']);
    fclose($input_file);
    $dataset_1_file = "./uploads/".$folder_input_name."/dataset_1.txt";
  }
}

$name       = $_FILES['dataset_2_file']['name'];
$temp_name  = $_FILES['dataset_2_file']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if(move_uploaded_file($temp_name, $location.$name)){
         $dataset_2_file = "./uploads/".$folder_input_name."/".$_FILES['dataset_2_file']['name'];
       }
   }
else {
  if ($_REQUEST['dataset_2'] != '') {
    $input_file = fopen('./uploads/'.$folder_input_name."/dataset_2.txt", "w");
    fwrite($input_file, $_REQUEST['dataset_2']);
    fclose($input_file);
    $dataset_2_file = "./uploads/".$folder_input_name."/dataset_2.txt";
  }
}

$comparison_type = $_REQUEST['comparison_type'];

$genetic_code = intval($_REQUEST['genetic_code']);

$folder_output_name = uniqid(date('H_i_s_d_m_y_'));

$exec_code = "python3 ./cousin/cousin.py compare_data -c ".$dataset_1_file." -d ".$dataset_2_file." -C ".$comparison_type." -g ".$genetic_code." -o ./outputs/".$folder_output_name." 2>&1";

$cousin_system_log = shell_exec($exec_code);

$log_cousin_system_file = "./cousin_logs_system/".$folder_output_name.".txt";

file_put_contents($log_cousin_system_file, $cousin_system_log);

delete_dir('./uploads/'.$folder_input_name);

echo($folder_output_name);

?>
