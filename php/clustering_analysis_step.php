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

$name       = $_FILES['fasta_file']['name'];
$temp_name  = $_FILES['fasta_file']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if(move_uploaded_file($temp_name, $location.$name)){
         $fasta_file = "./uploads/".$folder_input_name."/".$_FILES['fasta_file']['name'];
       }
   }
else {
  if ($_REQUEST['fasta_sequences'] != '') {
    $input_file = fopen('./uploads/'.$folder_input_name.'/fasta_sequences.fasta', "w");
    fwrite($input_file, $_REQUEST['fasta_sequences']);
    fclose($input_file);
    $fasta_file = "./uploads/".$folder_input_name."/fasta_sequences.fasta";
  }
}

$name       = $_FILES['cu_file']['name'];
$temp_name  = $_FILES['cu_file']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if (move_uploaded_file($temp_name, $location.$name)){
         $cu_file = "./uploads/".$folder_input_name."/".$_FILES['cu_file']['name'];
       }
   }
else {
  if ($_REQUEST['cu_table'] != '') {
    $input_file = fopen('./uploads/'.$folder_input_name."/cu_table.txt", "w");
    fwrite($input_file, $_REQUEST['cu_table']);
    fclose($input_file);
    $cu_file = "./uploads/".$folder_input_name."/cu_table.txt";
  }
}

$name       = $_FILES['optimal_codons_file']['name'];
$temp_name  = $_FILES['optimal_codons_file']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if (move_uploaded_file($temp_name, $location.$name)){
         $optimal_codons_file = "./uploads/".$folder_input_name."/".$_FILES['optimal_codons_file']['name'];
       }
   }
else {
  if ($_REQUEST['optimal_codons_text'] != '') {
    $input_file = fopen('./uploads/'.$folder_input_name."/optimal_codons_text.csv", "w");
    fwrite($input_file, $_REQUEST['optimal_codons_text']);
    fclose($input_file);
    $optimal_codons_file = "./uploads/".$folder_input_name."/optimal_codons_text.csv";
  }
  else {
    $optimal_codons_file = "nothing";
  }
}

$genetic_code = intval($_REQUEST['genetic_code']);

$name       = $_FILES['pattern_file_clustering']['name'];
$temp_name  = $_FILES['pattern_file_clustering']['tmp_name'];
if ((isset($name)) && (!empty($name))){
       $location = './uploads/'.$folder_input_name."/";
       if(move_uploaded_file($temp_name, $location.$name)){
        $number_cluster = "./uploads/".$folder_input_name."/".$_FILES['pattern_file_clustering']['name'];
       }
   }
else {
  if ($_REQUEST['pattern_data_clustering'] != '') {
    fopen ('./uploads/'.$folder_input_name."/pattern_data_clustering.txt", "w");
    fwrite($input_file, $_REQUEST['/pattern_data_clustering']);
    fclose($input_file);
    $number_cluster = "./uploads/".$folder_input_name."/pattern_data_clustering.txt";
  }
  else {
    $number_cluster = intval($_REQUEST['number_cluster']);
  }
}

$folder_output_name = uniqid(date('H_i_s_d_m_y_'));

if (isset($_REQUEST['simulation_routine'])) {
  $simulation_routine = $_REQUEST['simulation_routine'];
}
else {
  $simulation_routine = "No";
}

if (isset($_REQUEST['type_of_var'])) {
  $type_of_var = $_REQUEST['type_of_var'];
}
else {
  $type_of_var = "variables";
}

if (isset($_REQUEST['specific_COUSIN'])) {
  $specific_COUSIN = $_REQUEST['specific_COUSIN'];
}
else {
  $specific_COUSIN = "Yes";
}

if (isset($_REQUEST['vector_analysis'])) {
  $vector_analysis = $_REQUEST['vector_analysis'];
}
else {
  $vector_analysis = "No";
}

$exec_code = "python3 ./cousin/cousin.py clustering_analysis -s ".$fasta_file." -c ".$cu_file." -g ".$genetic_code." -S ".$simulation_routine." -X ".$specific_COUSIN." -J ".$vector_analysis." -b ".$optimal_codons_file." -w ".$type_of_var." -n ".$number_cluster." -o ./outputs/".$folder_output_name." 2>&1";

$cousin_system_log = shell_exec($exec_code);

$log_cousin_system_file = "./cousin_logs_system/".$folder_output_name.".txt";

file_put_contents($log_cousin_system_file, $cousin_system_log);

delete_dir('./uploads/'.$folder_input_name);

echo($folder_output_name);

?>