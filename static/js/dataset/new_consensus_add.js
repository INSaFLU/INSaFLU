
/// 
$().ready(function(){
	var name = $("#id_user_name").text().replace(/^\s+|\s+$/g, '');
	var lst_split = name.split(" ");
	if (lst_split.length == 3 & lst_split[0] == 'demo' & lst_split[1] == '-' & lst_split[2] == 'Logout'){
		$("#id_consensus_fasta").attr("disabled", "disabled");
	}
	else{
		$("#id_consensus_fasta").removeAttr("disabled");
	}
	$('#div_id_name').append('<div id="id_project_name_error"></div>')
	
});


//
//This test on the name change
//
$("#id_name").on("change paste keyup", function () {
  $.ajax({
    url: $("#id_form_consensus").attr("data-validate-consensus-url"),
   // data: form.serialize(),
    // dataType: 'json',
    data : { 
    	consensus_name : $('#id_name').val(),
		csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
    success: function (data) {
      $(document).find('#error_1_id_name').remove();
      if (data.is_taken) {
         var error_message = "<strong>" + data.error_message + "</strong>";
         $("#id_project_name_error").append(
            $('<span/>', { 										// creates a dynamic div element on the fly
	        id: 'error_1_id_name',
	        class: 'fields_error',
	        html: error_message,
    	 }));
      }
    },
    
    // handle a non-successful response
    error : function(xhr,errmsg,err) {
        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
    }
});
});