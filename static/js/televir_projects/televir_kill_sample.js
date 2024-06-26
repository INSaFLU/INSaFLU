


///
$(document).on("click", "a", function (e) {
    if ($(this).attr("id") === "id_kill_reference_modal") {
    
        var attr = $(this).attr('id');
        var ref_name = $(this).attr('ref_name');
        var sample_pk = $(this).attr('pk');
        
        // For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
        if (attr === 'id_kill_reference_modal'){
            $('#id-label-kill').text('Terminate \'' + ref_name + '\' jobs? Running and Queued jobs will be terminated.');
            $('#id-modal-body-kill-sample').attr('pk', sample_pk);
            $('#id-modal-body-kill-sample').attr('ref_name', ref_name);
        }
        else if (attr === 'id_add_kill_message'){
            $('#id_messages_kill').append('<div class="alert alert-dismissible alert-warning">' +
                    'No jobs to deploy.' +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
        }
    }
});


$('#id-kill-button').on('click', function(){

    url= $('#id-modal-body-kill-sample').attr("kill-single-value-url");
    sample_id = $('#id-modal-body-kill-sample').attr('pk');
    token = $('#id-modal-body-kill-sample').attr('csrfmiddlewaretoken');

	$.ajax({
        url: url,
        type: 'POST',
        data : {
        	sample_id : sample_id,
    		csrfmiddlewaretoken: token,
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
            if (data['is_empty']) {
              $('#id_messages_kill').append('<div class="alert alert-dismissible alert-warning">' +
                'No jobs to kill.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            } else {
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
                'The sample \'' + $('#id-modal-body-kill-sample').attr('ref_name') + '\' was successfully terminated.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            }
          }
          else{
        	/// add message with informaton
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
        		'The sample \'' + $('#id-modal-body-kill-sample').attr('ref_name') + '\' was not terminated.' +
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


$('#id-kill-all-button').on('click', function(){

  url= $('#id-modal-body-kill-all-sample').attr("remove-all-value-url");
  project_id = $('#id-kill-all-button').attr('project_id');
  token = $('#id-modal-body-kill-sample').attr('csrfmiddlewaretoken');

  $.ajax({
        url: url,
        type: 'POST',
        data : {
          project_id : project_id,
        csrfmiddlewaretoken: token,
        }, // data sent with the post request
            
    success: function (data) {
          if (data['is_ok']) {
            if (data['is_empty']) {
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                'No jobs to kill.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            } else {
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
                'Project + \'' + $('#id-kill-all-button').attr('project-name') + '\' jobs were successfully terminated.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
              location.reload();
            }
          }
          else{
          /// add message with informaton
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The sample \'' + $('#id-modal-body-kill-sample').attr('ref_name') + '\' was not terminated.' +
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