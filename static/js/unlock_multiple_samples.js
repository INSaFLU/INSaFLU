

/// remove all files
$('#id-unlock-button').on('click', function(){

	$.ajax({
        url: $('#id-modal-body-unlock-file-sample').attr("remove-single-value-url"),
        data : { 
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
      	  
        	  /// refresh page
        	  location.reload();
        	  
        	  /// add message with information
        	//  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
        	//	data['message_number_of_changes'] +
			//	'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
			//	'</div>');
          }
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            alert(errmsg);
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
	});
});