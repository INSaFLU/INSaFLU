

///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
  var ref_name = $(this).attr('ref_name');
  var sample_id = $(this).attr('sample_id');

	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_remove_reference_modal'){
		$('#id-label-remove').text('Do you want to remove reference \'' + ref_name + '\'?');
		$('#id-modal-body-remove-sample').attr('pk', $(this).attr('pk'));
    $('#id-modal-body-remove-sample').attr('ref_name', ref_name);
    $('#id-modal-body-remove-sample').attr('sample_id', sample_id);
	}
});


$('#id-remove-button').on('click', function () {
  
  url = $('#id-modal-body-remove-sample').attr("remove-single-value-url");
  var csrftoken = $('#id-modal-body-remove-sample').attr('data-csrf');

	$.ajax({
      url: $('#id-modal-body-remove-sample').attr("remove-single-value-url"),
      type: 'POST',
      data : { 
        reference_id: $('#id-modal-body-remove-sample').attr('pk'),
        sample_id: $('#id-modal-body-remove-sample').attr('sample_id'),
      csrfmiddlewaretoken: csrftoken,
      }, // data sent with the post request
          
      success: function (data) {
        if (data['is_ok']) {
          
          /// remove line
          ///document.getElementById($('#id-modal-body-remove-sample').attr('tr_to_remove')).remove();
          

          /// add message with informaton
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
          'The reference \'' + $('#id-modal-body-remove-sample').attr('ref_name') + '\' was successfully removed.' +
          '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
          '</div>');
          
          update_ref_table();

          
        }
        else{
          /// special message in case present_in_televir_project == True
          if (data['is_error']){
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                  'The reference \'' + $('#id-modal-body-remove-sample').attr('ref_name') + '\' was not removed.' +
                  '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                  '</div>');
          }
          }
      },
      
      // handle a non-successful response
      error : function(xhr,errmsg,err) {
          alert(errmsg);
          console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
      }
	});
});