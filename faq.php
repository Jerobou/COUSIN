<?php
	include('header.php');
?>

<!-- Page Content -->
<div id="page-content-wrapper">
		<div class="container-fluid">
				<div class="row">
						<div class="col-lg-12">
								<center><h1>Frequently Asked Questions</h1>
								<h2>Got a question ? We'll gladly help !  </h2>
								</center>
								<br>	
								<br>							
								<h5><b><u> Q : Why is COUSIN running with given examples, but not with my dataset? </u></b></h5>
								<h5> A : Make sure that your query sequences contain: </h5>
								<ul>
									<li> <b>codons</b>: your sequence should only contains codons, meaning that it should be divided by 3 without any resulting modulo.</li>	
									<li> <b>Only ATGC alphabet</b>: CUPrefs analysis underly that your dataset should contain precise information on said codon usage. Thus, all codons with incomplete information (IUPAC nomenclature, gaps, unrecognized alphabet, ...) must be removed.</li>
									<li> <b>Only non-termination codon except at the end</b>: a STOP codon at another location than the end of the sequence suggest that all nucleotides after this codon are not related to the coding part of your sequence. To correct that, remove the STOP codon, and if necessary, the nucleotides following this codon. </li>
								</ul>

								<br><br><br>

								<h5><b><u> Q : It's still not working! What about the Codon Usage Table (CUT)? </u></b></h5>
								<h5> A : your CUT should be formatted in the <a href="http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606" target="_blank">kazusa-like format</a>. Here is an example of a CUT made thanks to the COUSIN <a href="./create_table.php">CUT creation step</a>: </h5><br>
									<p style="font-family: monospace;">
										UUU 17.6(714298)  UCU 15.2(618711)  UAU 12.2(495699)  UGU 10.6(430311)<br>
										UUC 20.3(824692)  UCC 17.7(718892)  UAC 15.3(622407)  UGC 12.6(513028)<br>
										UUA  7.7(311881)  UCA 12.2(496448)  UAA  1.0( 40285)  UGA  1.6( 63237)<br>
										UUG 12.9(525688)  UCG  4.4(179419)  UAG  0.8( 32109)  UGG 13.2(535595)<br>
<br>
										CUU 13.2(536515)  CCU 17.5(713233)  CAU 10.9(441711)  CGU  4.5(184609)<br>
										CUC 19.6(796638)  CCC 19.8(804620)  CAC 15.1(613713)  CGC 10.4(423516)<br>
										CUA  7.2(290751)  CCA 16.9(688038)  CAA 12.3(501911)  CGA  6.2(250760)<br>
										CUG 39.6(1611801)  CCG  6.9(281570)  CAG 34.2(1391973)  CGG 11.4(464485)<br>
<br>
										AUU 16.0(650473)  ACU 13.1(533609)  AAU 17.0(689701)  AGU 12.1(493429)<br>
										AUC 20.8(846466)  ACC 18.9(768147)  AAC 19.1(776603)  AGC 19.5(791383)<br>
										AUA  7.5(304565)  ACA 15.1(614523)  AAA 24.4(993621)  AGA 12.2(494682)<br>
										AUG 22.0(896005)  ACG  6.1(246105)  AAG 31.9(1295568)  AGG 12.0(486463)<br>
<br>
										GUU 11.0(448607)  GCU 18.4(750096)  GAU 21.8(885429)  GGU 10.8(437126)<br>
										GUC 14.5(588138)  GCC 27.7(1127679)  GAC 25.1(1020595)  GGC 22.2(903565)<br>
										GUA  7.1(287712)  GCA 15.8(643471)  GAA 29.0(1177632)  GGA 16.5(669873)<br>
										GUG 28.1(1143534)  GCG  7.4(299495)  GAG 39.6(1609975)  GGG 16.5(669768)<br>
									</p>
								<br><br><br>

								<h5><b><u> Q : I'm a little lost. What should I do to begin? </u></b></h5>
								<h5> A : A good COUSIN analysis should be performed as following: </h5>
								<ol>
									<li>
									You should prepare your dataset. First of all, check that all your sequences are coding sequences that can be interpreted by COUSIN (see first FAQ question). Second, check if your Codon Usage Table (CUT) is in kazusa-like format.
									</li>
									<li>
										If you don't have any CUT, you can construct one using a dataset of coding sequences. For example, you could create the CUT of a organism using:
										<ul>
											<li>The set of its highly expressed genes</li>
											<li>Its whole set of genes</li>
										</ul>
										Please note that these are only examples. You should create a reference CUT from a comprehensive and study-fitting dataset.
									</li>
									<li>
										Perform a first COUSIN calculation on your query sequences. This should give you some first answers on GC-content and CUPrefs of your sequences.
									</li>
									<li>
										A clustering analysis using the "vector of codons" clustering option. This would aggregate your sequences following their CUPrefs.
									</li>
									<li>
										Enjoy your results.
									</li>
								</ol>


								<br><br><br><br>

								<h5><b> Got a specific question ? Ask me ! jerome.bourret@ird.fr </b></h5>
								<br><br><br>

						</div>
				</div>
		</div>
</div>
<!-- /#page-content-wrapper -->


<?php
	include('footer.php');
?>
