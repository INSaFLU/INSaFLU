//
//	Update the Dataset Name
//
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	var ref_name = $(this).attr('ref_name');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_rebuild_dataset_modal'){
		$('#id-label-rebuild').text('Do you want to rebuild results for  \'' + ref_name + '\'?');
		$('#id-modal-body-dataset_rebuild').attr('pk', $(this).attr('pk'));
		$('#id-modal-body-dataset_rebuild').attr('ref_name', ref_name);
	}
});

// rebuild
$('#id-rebuild-button').on('click', function(){

	$.ajax({
        url: $('#id-modal-body-dataset_rebuild').attr("dataset_rebuild-url"),
        data : {
        	dataset_id : $('#id-modal-body-dataset_rebuild').attr('pk'),
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	          	  
        	  /// add message with informaton
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
        		'The Dataset \'' + $('#id-modal-body-dataset_rebuild').attr('ref_name') + '\' is being rebuilt.' +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
          }
          else{
        	/// add message with informaton
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
        		'The Dataset \'' + $('#id-modal-body-dataset_rebuild').attr('ref_name') + '\' could not be rebuilt.' +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
          }
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            alert(errmsg);
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
	});
});