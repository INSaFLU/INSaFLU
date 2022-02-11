
// to keep data to mask
var dict_mask = {};
var mask_element_was_selected = "";

///
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'showMaskModal'){
		var i_to_change_id = $(this).attr('id_image');		// image in the table, to change color is has data...
		var reference_name = $(this).attr('reference_name');
		var project_name = $(this).attr('project_name');
		var td_to_update = e.target.parentNode.parentNode.parentNode.parentNode.id;
	
		//block all page
		wait_screen();

		$('#id_set_positions_to_mask_regions').attr('project_id', $(this).attr('project_id'));
		$('#id_set_positions_to_mask_regions').attr('project_sample_id', $(this).attr('project_sample_id'));
		$('#id_set_positions_to_mask_regions').attr('reference_name', reference_name);
		$('#id_set_positions_to_mask_regions').attr('project_name', project_name);
		$('#id_set_positions_to_mask_regions').attr('td_to_update', td_to_update);
		$('#id_set_positions_to_mask_regions').attr('i_to_change_id', i_to_change_id);	// image in the table, to change color is has data...
		$('#id-label-name').text('Reference: ' + reference_name);
		
		$.ajax({
	        url: $('#id_set_positions_to_mask_regions').attr("mask-consensus-actual-values-url"),
	        data : {
				project_id : $('#id_set_positions_to_mask_regions').attr('project_id'),
				project_sample_id : $('#id_set_positions_to_mask_regions').attr('project_sample_id'),
				csrfmiddlewaretoken: '{{ csrf_token }}',
	        }, // data sent with the post request
			success: function (data) {
          		if (data['is_ok']) {
					const combo_node = document.getElementById("combo_select_elements_mask_id");
					combo_node.textContent = '';
					mask_element_was_selected = '';

					dict_mask = JSON.parse(data['all_data']);
					for (element_name in dict_mask['__GeneticElement__']['dt_elements_mask']) {
						var element_to_add = '<option value="' + element_name + '">' + element_name + '</option>';
						$("#combo_select_elements_mask_id").append(element_to_add);
						if (mask_element_was_selected.length == 0){
							mask_element_was_selected = element_name;
							$("#id_select_length").text(" " + dict_mask['__GeneticElement__']['dt_elements_size'][mask_element_was_selected] + " bases");
						}
					}
					// warning project, it clean all data in project samples
					const parent = document.getElementById("id_message_warning_project")
					while (parent.firstChild) {
    					parent.firstChild.remove()
					}
					$('#id_message_warning_project')
					if (data['warning_project']){
						$('#id_message_warning_project').append('<div class="alert alert-dismissible alert-warning">' +
        					data['warning_project'] +
						'</div>');
					}
					
					$("#id_select_length").text(" " + dict_mask['__GeneticElement__']['dt_elements_size'][mask_element_was_selected] + " bases");
					$("#id_tile_mask_consensus").text("Mask consensus: " + project_name);
					document.getElementById("id_mask_sites").value = dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_sites'];
					document.getElementById("id_mask_from_beginning").value = dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_beginning'];
					document.getElementById("id_mask_till_ends").value = dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_ends'];
					document.getElementById("id_mask_regions").value = dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_regions'];
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
	else if (attr === 'id_add_sample_message'){
		$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        		'It is not possible for now to add samples to the projects because there is a big backlog. We are solving this problem right now.' +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
	}
});


//catch the element from tree nucleotides combo box
$("#combo_select_elements_mask_id").change(function () {
	change_mask_values();
});

// keep actual values and show the others
function change_mask_values() {
	var element_selected = $('#combo_select_elements_mask_id option:selected').val();

	// keep the data thar was setted
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_sites'] = $('#id_mask_sites').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_beginning'] = $('#id_mask_from_beginning').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_ends'] = $('#id_mask_till_ends').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_regions'] = $('#id_mask_regions').val();

	//// new data	
	mask_element_was_selected = element_selected; 

	/// set new data
	$("#id_select_length").text(" " + dict_mask['__GeneticElement__']['dt_elements_size'][element_selected] + " bases");
	document.getElementById("id_mask_sites").value = dict_mask['__GeneticElement__']['dt_elements_mask'][element_selected]['__MaskingConsensus__']['mask_sites'];
	document.getElementById("id_mask_from_beginning").value = dict_mask['__GeneticElement__']['dt_elements_mask'][element_selected]['__MaskingConsensus__']['mask_from_beginning'];
	document.getElementById("id_mask_till_ends").value = dict_mask['__GeneticElement__']['dt_elements_mask'][element_selected]['__MaskingConsensus__']['mask_from_ends'];
	document.getElementById("id_mask_regions").value = dict_mask['__GeneticElement__']['dt_elements_mask'][element_selected]['__MaskingConsensus__']['mask_regions'];
}


$('#id-mask-consensus').on('click', function(){

	/// set wait screen 
	wait_screen();

	/// get data that was in the modal
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_sites'] = $('#id_mask_sites').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_beginning'] = $('#id_mask_from_beginning').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_from_ends'] = $('#id_mask_till_ends').val();
	dict_mask['__GeneticElement__']['dt_elements_mask'][mask_element_was_selected]['__MaskingConsensus__']['mask_regions'] = $('#id_mask_regions').val();
	
	/// new to be post because of the size
	$.ajax({
		type: "POST",
        url: $('#id_set_positions_to_mask_regions').attr("mask-consensus-url"),
        data : {
			all_data : JSON.stringify(dict_mask, null, 4),
			project_id : $('#id_set_positions_to_mask_regions').attr('project_id'),
			project_sample_id : $('#id_set_positions_to_mask_regions').attr('project_sample_id'),
			csrfmiddlewaretoken: '{{ csrf_token }}',
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  
        	  /// add new class
			  var element = document.getElementById($('#id_set_positions_to_mask_regions').attr('i_to_change_id'));
			  if (element != null && typeof element !== "undefined" && element.value == '') {
        	  	element.setAttribute('title', data['new_title_i']);
			  	element.setAttribute('class', data['new_class_i']);
		      }
        	  
        	  /// add message with informaton
			  if (data['is_going_to_mask']){
        	  	$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-success">' +
        			data['message'] +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
			  }
			  else{
				$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
        			data['message'] +
					'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
					'</div>');
				}
				var element_to_update = document.getElementById($('#id_set_positions_to_mask_regions').attr('td_to_update'))
				element_to_update.getElementsByTagName("td")[4].textContent = data['default']
          	}
          	else{
				/// add message with informaton
				var project_sample_id = $('#id_set_positions_to_mask_regions').attr('project_sample_id');
				if ( project_sample_id !== null){
	        	  	$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
	        			'A project sample inside project \'' + $('#id_set_positions_to_mask_regions').attr('project_name') + '\' was not going to masked consensus.' +
						'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
						'</div>');
				}
				else{
					$('#id_messages_set_default').append('<div class="alert alert-dismissible alert-warning">' +
	        			'The project \'' + $('#id_set_positions_to_mask_regions').attr('project_name') + '\' was not going to masked consensus.' +
						'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
						'</div>');
				}
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