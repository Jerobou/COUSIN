<?php

ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

$zip_path = "./cousin.zip";
$zip = new ZipArchive();

if ($zip->open($zip_path, ZIPARCHIVE::CREATE | ZIPARCHIVE::OVERWRITE) !== TRUE) {
		die ("An error occurred creating your ZIP file.");
	}


$dir = "./cousin_for_download/";

$files_in_dir = new RecursiveIteratorIterator(
		new RecursiveDirectoryIterator($dir),
		RecursiveIteratorIterator::LEAVES_ONLY
);

foreach ($files_in_dir as $name => $file){
		// Skip directories (they would be added automatically)
	if (!$file->isDir())
	{
			// Add current file to archive
			$zip->addFile($file);
	}
	}

$zip->close();

function getFileType($file){
		$path_chunks = explode("/", $file);
		$thefile = $path_chunks[count($path_chunks) - 1];
		$dotpos = strrpos($thefile, ".");
		return strtolower(substr($thefile, $dotpos + 1));
}

$type = getFileType($zip_path);

header("Content-type: application/".$type);
header("Content-Disposition: attachment;filename=cousin.".$type);
header("Content-Transfer-Encoding: binary");
header("Content-length: " . filesize($zip_path));
header('Pragma: no-cache');
header('Expires: 0');
readfile($zip_path);
?>
