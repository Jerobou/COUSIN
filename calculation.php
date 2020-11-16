<?php
include('header.php');
?>

<script src="./js/calculation.js"></script>

<div id="page-content-wrapper">
	<div id="wrapper-options" class="collapse in">

		<!--Page contenant les listes d'offres de manière dynamique -->
		<div>
			<p style="right : 3px; top : 0; font-size : 90%; position : absolute"> * Theses parameters are mandatory </p>
		</div>
		<div id="options" class="panel panel-default">
			<center>
				<form id="calculation_form">
					<div id="option_1">
						<div class="container-fluid">
							<div class="form-fasta">
								<label for="fasta">FASTA sequences * </label>
								<textarea class="form-control" rows="3" id="fasta_sequences" name="fasta_sequences"></textarea>
								<input type="file" id="fasta_file" name="fasta_file">
							</div>
							<div class="form-cu-table">
								<label for="cu-table">Codon Usage Table (<a href="http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606" target="_blank">kazusa </a> format) * </label>
								<textarea class="form-control" rows="3" id="cu_table" name="cu_table"></textarea>
								<input type="file" id="cu_file" name="cu_file">
							</div>
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
									<option value="11">11 - The Bacterial, Archaeal and Plant Plastid Code</option>
									<option value="12">12 - The Alternative Yeast Nuclear Code</option>
									<option value="13">13 - The Ascidian Mitochondrial Code</option>
									<option value="14">14 - The Flatworm Mitochondrial Code</option>
									<option value="15">15 - Blepharisma Nuclear Code</option>
									<option value="16">16 - Chlorophycean Mitochondrial Code</option>
									<option value="21">21 - Trematode Mitochondrial Code</option>
									<option value="22">22 - Scenedesmus obliquus Mitochondrial Code</option>
									<option value="23">23 - Thraustochytrium Mitochondrial Code</option>
									<option value="24">24 - Pterobranchia Mitochondrial Code</option>
									<option value="25">25 - Candidate Division SR1 and Gracilibacteria Code</option>
									<option value="26">26 - Pachysolen tannophilus Nuclear Code</option>
									<option value="27">27 - Karyorelict Nuclear Code</option>
									<option value="28">28 - Condylostoma Nuclear Code</option>
									<option value="29">29 - Mesodinium Nuclear Code</option>
									<option value="30">30 - Peritrich Nuclear Code</option>
									<option value="31">31 - Blastocrithidia Nuclear Code</option>
									<option value="33">33 - Cephalodiscidae Mitochondrial UAA-Tyr Code</option>
								</select>
							</div>
						</div>
					</div>
					<div id="option_2">
						<div class="container-fluid">
							<div class="form-optimal-codons">
								<label for="optimal-codons"> Optimal codons (FOP, CBI methods) </label>
								<textarea class="form-control" rows="3" id="optimal_codons_text" name="optimal_codons_text"></textarea>
								<input type="file" id="optimal_codons_file" name="optimal_codons_file">

							</div>
							<div id="square_box_options" name="square_box_options" class="square_box_options">
								<div class="form-simulation-routine">
									<input type="checkbox" id="simulation_routine" name="simulation_routine" value="Yes">
									<label for="simulation_routine">Perform basic simulation ?</label>
								</div>
								<div class="form-specific-COUSIN">
									<input type="checkbox" id="specific_COUSIN" name="specific_COUSIN" value="No">
									<label for="specific_COUSIN">Remove COUSIN boundaries (default -3;4)</label>
								</div>
								<div class="form-specific-vector_analysis">
									<input type="checkbox" id="vector_analysis" name="vector_analysis" value="Yes">
									<label for="vector_analysis">Create vectors for each query ?</label>
								</div>
							</div>
							<br>
							<a id="example_button" class="btn btn-primary"><i class="fa fa-cogs fa-lg"></i> <span>Example</span></a>

						</div>
					</div>

			</center>
		</div>
	</div>

	<a id="collapse_options" class="btn btn-primary" data-toggle="collapse" href="#wrapper-options"><i class="fa fa-minus-square fa-lg"></i> <span>Hide options</span><span class="hidden">Show options</span></a>
	<script>
		$('#collapse_options').click(function() {
			$('.hidden').removeClass('hidden').hide();
			$(this).find('i').toggleClass('fa-plus-square fa-lg fa-minus-square fa-lg')
			$(this).find('span').each(function() {
				$(this).toggle();
			});
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