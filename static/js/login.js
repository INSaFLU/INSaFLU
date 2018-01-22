


$('#id_login_anonymous').click( function(){
	if ($('#id_login_anonymous').is(":checked")){
		document.getElementById('id_username').value = "demo";
		document.getElementById('id_password').value = "demo_user";
	}
	else{
		document.getElementById('id_username').value = "";
		document.getElementById('id_password').value = "";
	}
	
});