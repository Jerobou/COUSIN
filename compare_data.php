<?php
	include('header.php');
?>

<script src="./js/compare_data.js"></script>

<div id="page-content-wrapper">
	<div id="wrapper-options" class="collapse in">

	<!--Page contenant les listes d'offres de manière dynamique -->
	<div>
		<p style="right : 3px; top : 0; font-size : 90%; position : absolute"> * Theses parameters are mandatory </p>
	</div>
		<div id="options" class="panel panel-default">
			<center>
		<form id="compare-data-form">
			<div id = "option_1">
		<div class="container-fluid">
				<div class="form-comparison-type">
					<label for="comparison-type">Comparison type * </label>
					<select type="text" class="form-control" id="comparison_type" name="comparison_type" required>
						<option value="cvc">CVC - Codon Usage Table against another One</option>
						<option value="cvs">CVS - Codon Usage Table against a list of sequences</option>
						<option value="svc">SVC - List of sequences againt a Codon Usage Table</option>
						<option value="svs">SVS - Liste of sequences against another one</option>

					</select>
				</div>
						<div class="form-dataset-1">
							<label for="dataset-1">First dataset (CUT or FASTA file) * </label>
							<textarea class="form-control" rows="3" id="dataset_1" name="dataset_1"></textarea>
							<input type="file" id="dataset_1_file" name="dataset_1_file">
						</div>
						<div class="form-dataset-2">
							<label for="dataset-2">Second dataset (CUT or FASTA file) * </label>
							<textarea class="form-control" rows="3" id="dataset_2" name="dataset_2"></textarea>
							<input type="file" id="dataset_2_file" name="dataset_2_file">
						</div>
		</div>
	</div>

	<div id = "option_2">
<div class="container-fluid">
	<div class="form-genetic-code">
		<label for="genetic-code">Genetic Code * </label>
		<select type="text" class="form-control" id="genetic_code" name="genetic_code" required>
			<option value="1">1 - The Standard Code</option>
			<option value="2">2 - The Vertebrate Mitochondrial Code</option>
			<option value="3">3 - The Yeast Mitochondrial Code</option>
			<option value="4">4 - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code</option>
			<option value="5">5 - The Invertebrate Mitochondrial Code</option>
			<option value="6">6 - The Ciliate, Dasycladacean and Hexamita Nuclear Code</option>
			<option value="9">9 - The Echinoderm Mitochondrial Code</option>
			<option value="10">10 - The Euplotid Nuclear Code</option>
			<option value="11">11 - The Bacterial Code</option>
			<option value="12">12 - The Alternative Yeast Nuclear Code</option>
			<option value="13">13 - The Ascidian Mitochondrial Code</option>
			<option value="14">14 - The Flatworm Mitochondrial Code</option>
			<option value="15">15 - Blepharisma Nuclear Code</option>
		</select>
	</div>
</div>

</div>

</center>

</div>
</div>

<a id="collapse_options" class="btn btn-primary" data-toggle="collapse" href="#wrapper-options"><i class="fa fa-minus-square fa-lg"></i> <span>Hide options</span><span class="hidden">Show options</span></a>
<script>
$('#collapse_options').click(function(){
	$('.hidden').removeClass('hidden').hide();
	$(this).find('i').toggleClass('fa-plus-square fa-lg fa-minus-square fa-lg')
	$(this).find('span').each(function() { $(this).toggle(); });
});
</script>

<center>
<button type="submit" class="btn btn-primary">Do the analysis</button>
</center>
</form>

<div id=content>
</div>

</div>

<?php
	include('footer.php');
?>
