

$("#id_path_name_1").click( function(){
	$('#error_1_id_path_name_1').empty();
});
$("#id_path_name_2").click( function(){
	$('#error_1_id_path_name_2').empty();
});

// for data picker
$('#id_date_of_onset').datepicker({
    uiLibrary: 'bootstrap4',
    iconsLibrary: 'fontawesome',
    format: 'dd/mm/yyyy',
});

$('#id_date_of_collection').datepicker({
    uiLibrary: 'bootstrap4',
    iconsLibrary: 'fontawesome',
    format: 'dd/mm/yyyy',
});

$('#id_date_of_receipt_lab').datepicker({
    uiLibrary: 'bootstrap4',
    iconsLibrary: 'fontawesome',
    format: 'dd/mm/yyyy',
    placeholder: 'First name',
});


/// 
$(document).on("click", "a", function(){
	var attr = $(this).attr('id');
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_data_set_add_modal'){
		$('h4.modal-title').text('Add a data set');
		$('h4.modal-title').attr('id', 'id_data_set_add_modal');
		$('#id_name_to_insert').attr('placeholder', 'Data set name');
		$('#id-label-replace').text('Set a new Data set:');
	}
	else if (attr === 'id_vaccine_add_modal'){
		$('h4.modal-title').text('Set a new vaccine status');
		$('h4.modal-title').attr('id', 'id_vaccine_add_modal');
		$('#id_name_to_insert').attr('placeholder', 'Vaccine status');
		$('#id-label-replace').text('Vaccine status:');
	}
	$('#id_name_to_insert').val('');
});


$('#id-save-button').on('click', function(){
	//id_name_to_insert
	$.ajax({
        url: $('#modal-body-add-sample').attr("add-single-value-url"),
        /// need to add crfs
        data : { type_data : $('h4.modal-title').attr('id'),
        		value : $('#id_name_to_insert').val(),
    			csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  /// add the value to the combo box
        	  if ($('h4.modal-title').attr('id') == 'id_data_set_add_modal'){
        		  // add to dataset
        		  $('#id_data_set').append($('<option>', {
        			    value: data['value'],
        			    text: data['text'],
        			}));
        	  }
        	  else if ($('h4.modal-title').attr('id') == 'id_vaccine_add_modal'){
        		  // add to vaccine
        		  $('#id_vaccine_status').append($('<option>', {
        			    value: data['value'],
        			    text: data['text'],
        		  }));
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

// to remove
$(document).on("click", "a", function(){
	var attr = $(this).attr('id');
	var value_to_remove = "";
	// For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
	if (attr === 'id_data_set_remove_modal'){
		$('h4.modal-title-remove').text('Remove a data set');
		$('h4.modal-title-remove').attr('id', 'id_data_set_remove_modal');
		value_to_remove = $( "#id_data_set option:selected" ).text();
		$('#id-label-remove').text('Do you want to remove ' + $( "#id_data_set option:selected" ).text() + '?');
	}
	else if (attr === 'id_vaccine_remove_modal'){
		$('h4.modal-title-remove').text('Remove a vaccine status');
		$('h4.modal-title-remove').attr('id', 'id_vaccine_remove_modal');
		value_to_remove = $( "#id_vaccine_status option:selected" ).text();
		$('#id-label-remove').text('Do you want to remove ' + $( "#id_vaccine_status option:selected" ).text() + '?');
	}
	// ajax to remove
	if (attr === 'id_data_set_remove_modal' || attr === 'id_vaccine_remove_modal'){
		$.ajax({
	        url: $('#modal-body-remove-sample').attr("remove-single-value-url"),
	        /// need to add crfs
	        data : { type_data : attr,
	        		is_to_test : true,
	        		value : value_to_remove,
	    			csrfmiddlewaretoken: '{{ csrf_token }}' 
	        }, // data sent with the post request
	        		
	        success: function (data) {
	          if (data['is_ok']) {
	        	  /// add the value to the combo box
    			  if (data['is_can_remove']){
    				  $('#id-remove-button').remove()
    				  $('#id-modal-footer-remove').append('<button type="button" id="id-remove-button" class="btn btn-primary" data-dismiss="modal">Remove</button>');
    			  }
    			  else{
    				  $('#id-label-remove').text(data['message']);
    				  /// remove the OK button
    				  $('#id-remove-button').remove()
    			  }
	          }
	        },
	        
	        // handle a non-successful response
	        error : function(xhr,errmsg,err) {
	            alert(errmsg);
	            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	        }
		});
	}
});

$('#id-modal-footer-remove').on('click', '#id-remove-button', function(){
	
	if ($('h4.modal-title-remove').attr('id') == 'id_data_set_remove_modal'){
		value_to_remove = $( "#id_data_set option:selected" ).text();
	}
	else if ($('h4.modal-title-remove').attr('id') == 'id_vaccine_remove_modal'){
		value_to_remove = $( "#id_data_set option:selected" ).text();
	}
	
	$.ajax({
        url: $('#modal-body-remove-sample').attr("remove-single-value-url"),
        /// need to add crfs
        data : { type_data : $('h4.modal-title-remove').attr('id'),
        		is_to_test : false,
        		value : value_to_remove,
    			csrfmiddlewaretoken: '{{ csrf_token }}' 
        }, // data sent with the guest request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  /// add the value to the combo box
        	  if ($('h4.modal-title-remove').attr('id') == 'id_data_set_remove_modal'){
        		  // remove to dataset
        		  $("#id_data_set option[value='" + data['value_to_remove'] + "']").remove();
        	  }
        	  else if ($('h4.modal-title-remove').attr('id') == 'id_vaccine_remove_modal'){
        		  // remove to vaccine
        		  $("#id_vaccine_status option[value='" + data['value_to_remove'] + "']").remove();
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




