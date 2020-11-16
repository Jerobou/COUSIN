<?php
	include('header.php');
?>

<!-- Page Content -->
<div id="page-content-wrapper">
		<div class="container-fluid">
				<div class="row">
						<div class="col-lg-12">
								<center><h1>User's guide </h1>
									<br><br>
								<h4>Dear user, this page is for you !  </h4>
								<h4>First of all, thanks a lot for testing COUSIN. </h4>
								</center>
								<br><br>
								<h5> COUSIN is still in beta-test, and can show non user-friendly features ! </h5>
								<h5> This page has been created to answer some of your questions. If you see any bug while using COUSIN, please check the following requirements list: </h5>
							<br>
								
								<ul>
									<li> <h5> Are your query sequences CDSs ? Can they be divided by 3 and if they got a STOP codon, is it at the end of the sequences ?  </h5> </li>
									<li> <h5> Is your reference dataset in a <a href="http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606"> kazusa-style </a> format ? Is it informative ? </h5> </li>
									<li> <h5> Are your sure you selected a genetic code related to your queries and your reference table ? Are all your queries sequences following the same genetic code ? </h5> </li>
								</ul>


							<br>

							<h5> If you're not sure about what's happening, go check the "COUSIN system log" given after an analysis (if you can access to it !). This will give you the terminal informations. Please note that you may see Python3 warning messages in this system log. This is perfectly normal. </h5>

							<p> If you see "strange" error message while using COUSIN, do not worry. It could be due to a bug listed above, or to something else (sorry...). Since the website is still under construction, no description of bugs are given right now. Try following the data formats given by the examples to perform your own analysis !</p>

						</div>
				</div>
		</div>
</div>
<!-- /#page-content-wrapper -->


<?php
	include('footer.php');
?>
