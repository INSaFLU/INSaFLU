//
//	This test on the name change
//
$("#id_project_name").on("change paste keyup", function () {
      $.ajax({
        url: $("#id_form_table_project_reference").attr("data-validate-project-reference-url"),
       // data: form.serialize(),
        // dataType: 'json',
        data : { 
        	project_name : $('#id_project_name').val(),
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


/// add function to toggle the checkbox in tables
/// set the Listener and get the checked in the server, and set the box check in the client
$(document).ready(function(){

	var elements = document.getElementsByName("select_ref");
	for(var i = 0, n = elements.length; i < n; i++){
		elements[i].addEventListener('click', toggle_check_box, false);
	}

	// get if there's any checked in the server
	$.ajax({
		url: $('#table_with_check_id').attr("set-check-box-values-url"),
		data : { 
			get_check_box_single : '1',
			csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
		success: function (data) {
			for (key in data){
				if (key === 'is_ok'){ continue; }
				var remember = document.getElementById(key);
				if (remember != null){
					remember.checked = data[key];
				}
			}
		},
	});
});


/// toggle checkbox
function toggle_check_box(source) {
	$(document).find('#id_reference_error').remove();
	$.ajax({
		url: $('#table_with_check_id').attr("set-check-box-values-url"),
		data : { get_change_check_box_single : source.target.checked,
				 value : source.target.value }, // data sent with the post request
		success: function (data) { 
			for (key in data){
				if (key === 'is_ok'){ continue; }
				var remember = document.getElementById(key);
				if (remember != null){
					remember.checked = data[key];
				}
				else{
					remember.checked = false;
				}
			}
		},
	});
};