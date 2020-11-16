<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$dir = "./cousin_tuto";

$tutorial = $dir."/cousin_tuto.pdf";

function getFileType($file){
		$path_chunks = explode("/", $file);
		$thefile = $path_chunks[count($path_chunks) - 1];
		$dotpos = strrpos($thefile, ".");
		return strtolower(substr($thefile, $dotpos + 1));
}

$type = getFileType($tutorial);

header("Content-type: application/".$type);
header("Content-Disposition: attachment;filename=cousin_tuto.".$type);
header("Content-Transfer-Encoding: binary");
header("Content-length: " . filesize($tutorial));
header('Pragma: no-cache');
header('Expires: 0');
readfile($tutorial);
?>

?>
