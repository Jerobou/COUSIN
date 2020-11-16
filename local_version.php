<?php
	include('header.php');
?>

<script src="./js/local_version.js"></script>


<div>
	<div id="page-content-wrapper">

<center>
<h2> COUSIN local version </h2>
</center>
<br>

	<h5> The COUSIN tool can be donwloaded and installed on a local machine for a personal use. Please be aware that a Python3 version is needed with the following modules: </h5>

	<ul>
		<h5> <li> pandas, numpy & scipy </li> </h5>
		<h5> <li> sklearn & pyclustering </li> </h5>
		<h5> <li> matplotlib & seaborn </li> </h5>
	</ul>

	<br>

	<h5> We strongly advise to check the <a href="./tutorial.php"> COUSIN tutorial </a> before using the local version. </h5>
	<br>

	<h5> To download the COUSIN tool, please fill this form : </h5>
	<br>
	<br>
		<div>
		<form id="local-version-form">
			<div class="form-group">
    			<label for="first_name">First name</label>
    			<input type="text" name="first_name" class="form-control" id="first_name" placeholder="First name" required>
  			</div>
  			<div class="form-group">
    			<label for="last_name">Last name</label>
    			<input type="text" name="last_name" class="form-control" id="last_name" placeholder="Last name" required>
  			</div>
  			  			<div class="form-group">
    			<label for="status">Status</label>
    			<select type="text" name="status" class="form-control" id="status" required>
    				<option value="Master">Md. student</option>
    				<option value="PhD-student">PhD. student</option>
    				<option value="technician">technician</option>
    				<option value="engineer">engineer</option>
    				<option value="PhD-researcher">PhD. researcher</option>
    				<option value="Prof">Professor</option>
    				<option value="other">Other</option>					
    			</select>
    			
  			<div class="form-group">
    			<label for="email_adress">Email address</label>
    			<input type="email" name="email_adress" class="form-control" id="email_adress" aria-describedby="emailHelp" placeholder="Enter email" required>
    			<small id="emailHelp" class="form-text text-muted">We'll never share your data with anyone else.</small>
  			</div>
  			<center>
  			<button type="submit" id="download_COUSIN" class="btn btn-success""><i class="fa fa-download"></i>   DOWNLOAD THE COUSIN TOOL (V 1.0)</button>
  			</center>
		</form>
	</div>

</div>
</div>
<?php
	include('footer.php');
?>
