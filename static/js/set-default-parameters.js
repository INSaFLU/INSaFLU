///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	var ref_name = $(this).attr('ref_name');
	var td_to_update = e.target.parentNode.parentNode.parentNode.parentNode.id;
	var proj_name = $(this).attr('proj_name');
	var pk_proj = $(this).attr('pk_proj');
	var pk_proj_sample = $(this).attr('pk_proj_sample');
	var pk_sample = $(this).attr('pk_sample');
	var type_software = $(this).attr('type_software');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_default_parameter'){
		$('#id-label-set-default').text('Do you want to set default values for \'' + ref_name + '\'?');
		$('#id-modal-body-set-default').attr('pk', $(this).attr('pk'));		/// software ID
		$('#id-modal-body-set-default').attr('ref_name', ref_name);
		$('#id-modal-body-set-default').attr('td_to_update', td_to_update);
		$('#id-modal-body-set-default').attr('type_software', type_software);
		
		// info about project or project_sample to set defaults in software
		if ( typeof proj_name !== typeof undefined && proj_name !== false ) {
			$('#id-modal-body-set-default').attr('proj_name', proj_name );
			
			if ( typeof pk_proj !== typeof undefined && pk_proj !== false ) {
				$('#id-modal-body-set-default').attr('pk_proj', pk_proj );
			}
			if ( typeof pk_proj_sample !== typeof undefined && pk_proj_sample !== false ) {
				$('#id-modal-body-set-default').attr('pk_proj_sample', pk_proj_sample );
			}
			if ( typeof pk_sample !== typeof undefined && pk_sample !== false ) {
				$('#id-modal-body-set-default').attr('pk_sample', pk_sample );
			}
		}
	}
});


$('#id-set-default-button').on('click', function(){

	$.ajax({
        url: $('#id-modal-body-set-default').attr("remove-single-value-url"),
        data : { 
        	software_id : $('#id-modal-body-set-default').attr('pk'),
        	project_name : $('#id-modal-body-set-default').attr('proj_name'),			/// can be Null
        	project_id : $('#id-modal-body-set-default').attr('pk_proj'),					/// can be Null
        	project_sample_id : $('#id-modal-body-set-default').attr('pk_proj_sample'),	/// can be Null
        	sample_id : $('#id-modal-body-set-default').attr('pk_sample'),	/// can be Null
    		csrfmiddlewaretoken : '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  /// add message with informaton
        	  $('#id_messages_set_default').append('<div class="alert alert-dismissible alert-success">' +
        		'The ' + $('#id-modal-body-set-default').attr('type_software') + 
        		' \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' was set to default values.' +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
        	  var element_to_update = document.getElementById($('#id-modal-body-set-default').attr('td_to_update'))
        	  if ($('#id-modal-body-set-default').attr('type_software') === 'software'){
        		  element_to_update.getElementsByTagName("td")[3].textContent = data['default']
        	  }
        	  else{
        		  element_to_update.getElementsByTagName("td")[2].textContent = data['default']
        	  }
          }
          else{
        	/// add message with informaton
        	  $('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        		'The ' + $('#id-modal-body-set-default').attr('type_software') + 
        		' \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' was not set to default values.' +
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