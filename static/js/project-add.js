		//
//	This test on the submit
//
$('#form_table_project_reference').on('submit', function () {
    var total_checked = 0;
    
    var form = $(this).closest("form");
    $.ajax({
        url: form.attr("data-validate-project-reference-url"), // the endpoint
        data : { project_name : $('#id_project_name').val(),
        		user_name : $('#id_user_name').text() }, // data sent with the post request

		success: function (data) {
          if (data.is_taken) {
          		alert(data.error_message);
          		return False;
          }
        }
    });
    
    $('#form_table_project_reference input:checked').each(function() {
		total_checked += 1;
    });
    if (total_checked == 0){
     	alert('You need to select one reference!');
     	return false;
    }
    else if (total_checked > 1){
    	alert('There is more than one reference selected!');
     	return false;
    }
    return true;
});

//
//	This test on the name change
//
$("#id_project_name").on("change paste keyup", function () {
      var form = $(this).closest("form");
      $.ajax({
        url: form.attr("data-validate-project-reference-url"),
       // data: form.serialize(),
        // dataType: 'json',
        data : { project_name : $('#id_project_name').val(),
        		user_name : $('#id_user_name').text() }, // data sent with the post request
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
            alert(errmsg);
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
   });
});
