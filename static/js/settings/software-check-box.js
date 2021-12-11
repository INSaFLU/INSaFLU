////////////////////////////
///
///		Only the check box table
///


$('#id-set-turn-on-off-button').on('click', function(){

	var software_id = $('#id-set-turn-on-off-button').attr('software_id');
	//block all page
	$.blockUI({ css: { 
	            border: 'none', 
	            padding: '15px', 
	            backgroundColor: '#000', 
	            '-webkit-border-radius': '10px', 
	            '-moz-border-radius': '10px', 
	            opacity: .5, 
	            color: '#fff' 
	        } }); 
	
	$.ajax({
        url: '/settings/ajax/turn_on_off_software',
        data : {
        	software_id : software_id,
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
	        	/// add message with informaton
	        	  $('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
	        		'Something went rong.' +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
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
	var software_id = $(this).attr('software_id');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_show_turn_on_off_modal'){
		//block all page
		$.blockUI({ css: { 
		            border: 'none', 
		            padding: '15px', 
		            backgroundColor: '#000', 
		            '-webkit-border-radius': '10px', 
		            '-moz-border-radius': '10px', 
		            opacity: .5, 
		            color: '#fff' 
		        } });
		
		$.ajax({
	        url: '/settings/ajax/get_software_name_to_turn_on_off',
	        data : {
	        	software_id : software_id,
	    		csrfmiddlewaretoken: '{{ csrf_token }}'
	        }, // data sent with the post request
	        		
	        success: function (data) {
	          	if (data['is_ok']) {
					$('#id-set-turn-on-off-button').attr('software_id', software_id);
					$('#id-label-tur-on-off').text(data['message']);
					$('#id-set-turn-on-off-button').removeAttr('disabled');
	         	}
				else{
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
	

	}
});


