////////////////////////////
///
///		Only the check box table
///


$('#id-set-turn-on-off-button').on('click', function(){

	var software_id = $('#id-set-turn-on-off-button').attr('software_id');
	var project_id = $('#id-set-turn-on-off-button').attr('project_id');
	var project_sample_id = $('#id-set-turn-on-off-button').attr('project_sample_id');
	var sample_id = $('#id-set-turn-on-off-button').attr('sample_id');
	var televir_project_id= $('#id-set-turn-on-off-button').attr('televir_project_id');
	var type_of_use_id = $('#id-set-turn-on-off-button').attr('type_of_use_id');
	//block all page
	wait_screen();
	
	$.ajax({
        url: '/settings/ajax/turn_on_off_software',
        data : {
        	software_id : software_id,
        	project_id : project_id,
        	project_sample_id : project_sample_id,
        	sample_id : sample_id,
			type_of_use_id: type_of_use_id,
			televir_project_id: televir_project_id,
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
	    	if (data['is_ok']) {
				var element = document.getElementById('check_box_' + software_id);
				element.checked = data['is_to_run'];
	        	  
	        	/// add message with informaton
	        	$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-success">' +
	        		data['message'] +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
	        }
	        else{
	        	/// add message with warning
				$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        			data['message'] +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
	        		$('#id-label-tur-on-off').text(data['message']);
					$('#id-set-turn-on-off-button').attr('disabled','disabled');
	         }
			 $.unblockUI();
	        },
	        
	        // handle a non-successful response
	        error : function(xhr,errmsg,err) {
	            alert(errmsg);
	            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
				$.unblockUI();
	        }
	});
});

///////////////////////
///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_show_turn_on_off_modal'){
		var software_id = $(this).attr('software_id');
		var type_of_use_id = $(this).attr('type_of_use_id');
        var technology_id = $(this).attr('technology_id');
		var televir_project_id= $(this).attr('televir_project_id');
		var project_id = $('#id_show_turn_on_off_modal').attr('project_id');
		var project_sample_id = $('#id_show_turn_on_off_modal').attr('project_sample_id');
		var sample_id = $('#id_show_turn_on_off_modal').attr('sample_id');
	
		//block all page
		wait_screen();
		
		$.ajax({
	        url: '/settings/ajax/get_software_name_to_turn_on_off',
	        data : {
	        	software_id : software_id,
				type_of_use_id: type_of_use_id,
				project_id : project_id,
        		project_sample_id : project_sample_id,
        		sample_id : sample_id,
				televir_project_id: televir_project_id,
                technology_id: technology_id,
	    		csrfmiddlewaretoken: '{{ csrf_token }}'
	        }, // data sent with the post request
	        		
	        success: function (data) {
	          	if (data['is_ok']) {
					$('#id-set-turn-on-off-button').attr('software_id', software_id);
					$('#id-set-turn-on-off-button').attr('project_id', project_id);
					$('#id-set-turn-on-off-button').attr('project_sample_id', project_sample_id);
					$('#id-set-turn-on-off-button').attr('sample_id', sample_id);
					$('#id-set-turn-on-off-button').attr('televir_project_id', televir_project_id);
					$('#id-set-turn-on-off-button').attr('type_of_use_id', type_of_use_id);
					$('#id-set-turn-on-off-button').attr('technology_id', technology_id);
					$('#id-label-turn-on-off').text(data['message']);
					$('#id-set-turn-on-off-button').removeAttr('disabled');
	         	}
				else{
					$('#id-label-turn-on-off').text(data['message']);
					$('#id-set-turn-on-off-button').attr('disabled','disabled');			
				}
				$.unblockUI();
	        },
	        
	        // handle a non-successful response
	        error : function(xhr,errmsg,err) {
	            alert(errmsg);
	            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
				$.unblockUI();
	        }
		});
	

	}
});


