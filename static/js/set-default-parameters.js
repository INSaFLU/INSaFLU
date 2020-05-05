///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	var ref_name = $(this).attr('ref_name');
	var td_to_update = e.target.parentNode.parentNode.parentNode.parentNode.id;
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_default_parameter'){
		$('#id-label-set-default').text('Do you want to set default values for \'' + ref_name + '\'?');
		$('#id-modal-body-set-default').attr('pk', $(this).attr('pk'));
		$('#id-modal-body-set-default').attr('ref_name', ref_name);
		$('#id-modal-body-set-default').attr('td_to_update', td_to_update);
	}
});


$('#id-set-default-button').on('click', function(){

	$.ajax({
        url: $('#id-modal-body-set-default').attr("remove-single-value-url"),
        data : { 
        	software_id : $('#id-modal-body-set-default').attr('pk'),
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  /// add message with informaton
        	  $('#id_messages_set_default').append('<div class="alert alert-dismissible alert-success">' +
        		'The software \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' was set to default values.' +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
        	  var element_to_update = document.getElementById($('#id-modal-body-set-default').attr('td_to_update'))
        	  element_to_update.getElementsByTagName("td")[2].textContent = data['default']
          }
          else{
        	/// add message with informaton
        	  $('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        		'The software \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' was not set to default values.' +
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