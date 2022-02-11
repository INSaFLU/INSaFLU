///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_default_parameter'){
		var ref_name = $(this).attr('ref_name');
		var td_to_update = e.target.parentNode.parentNode.parentNode.parentNode.id;
		var proj_name = $(this).attr('proj_name');
		var pk_proj = $(this).attr('pk_proj');
		var pk_proj_sample = $(this).attr('pk_proj_sample');
		var pk_sample = $(this).attr('pk_sample');
		var type_software = $(this).attr('type_software');
		var tooltip_title = $(this).attr('title');
		var id_image = $(this).attr('id_image');
	
		var message = 'Do you want to ' + tooltip_title.toLowerCase() + ' for \'' + ref_name + '\'?';
		/// Only show this message in projects
		if (ref_name === 'Mask consensus by sites' && typeof pk_proj_sample === typeof undefined) {
			message += "\nThis will clean all sites that are defined even in samples attached to this project.";
		}
		$('#id-label-set-default').text(message);
		$('#id-modal-title-remove').text(tooltip_title);
		$('#id-modal-body-set-default').attr('pk', $(this).attr('pk'));		/// software ID
		$('#id-modal-body-set-default').attr('ref_name', ref_name);
		$('#id-modal-body-set-default').attr('td_to_update', td_to_update);
		$('#id-modal-body-set-default').attr('type_software', type_software);
		$('#id-modal-body-set-default').attr('id_image', id_image);
		
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

// set default values
$('#id-set-default-button').on('click', function(){

	/// set wait screen 
	wait_screen();

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
					' \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' ' + data['message'] +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
				var element_to_update = document.getElementById($('#id-modal-body-set-default').attr('td_to_update'))
				if ($('#id-modal-body-set-default').attr('type_software') === 'software'){
					element_to_update.getElementsByTagName("td")[4].textContent = data['default']
				}
				else{
					element_to_update.getElementsByTagName("td")[3].textContent = data['default']
				}

				/// clean yellow from the icon in edit MaskConsensus,
				/// Only works for "Mask consensus by sites" software
				var software_name = $('#id-modal-body-set-default').attr('ref_name');
				if (software_name === 'Mask consensus by sites'){
					var for_image_id = $('#id-modal-body-set-default').attr('pk_proj');
					if (typeof for_image_id === typeof undefined || for_image_id === null){
						for_image_id = $('#id-modal-body-set-default').attr('pk_proj_sample');
					}
					
					if (typeof for_image_id !== typeof undefined && for_image_id !== null){
			 			var element_image = document.getElementById('icon_mask_consensus_' + for_image_id);
						if ( typeof element_image !== typeof undefined && element_image !== null ){
							var class_data = element_image.getAttribute('class');
							element_image.setAttribute('class', class_data.replace("warning_fa_icon", ""));
						}
					}
					
					// clean the tooltip
					var element = document.getElementById($('#id_set_positions_to_mask_regions').attr('i_to_change_id'));
					if (element != null && typeof element !== "undefined" && element.value == '') {
						element.setAttribute('title', 'Edit parameters');
					}					
				}
				$.unblockUI();
				/// END clean yellow from the icon in edit MaskConsensus
          	}
          	else{
        		/// add message with informaton
        	  	$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        			'The ' + $('#id-modal-body-set-default').attr('type_software') + 
        			' \'' + $('#id-modal-body-set-default').attr('ref_name') + '\' was NOT set to default values.' +
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